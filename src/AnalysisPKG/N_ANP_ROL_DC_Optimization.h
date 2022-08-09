//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Purpose       : 
// Special Notes :
// Creator       : 
// Creation Date : 
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_ROL_DC_Optimizationh
#define Xyce_N_ANP_ROL_DC_Optimizationh

#include <Xyce_config.h>

#ifndef FD_HESSIAN
#define FD_HESSIAN 0 // 0 to set all hessvec to zero (Gauss-Newton Hessian); else FD Hessian
#endif

#ifndef IDENTITY_PRECONDITIONER
#define IDENTITY_PRECONDITIONER 1
#endif

#ifdef Xyce_ROL

namespace Xyce {
namespace Nonlinear {
template <typename ScalarT>
class objectiveFunctionData;
}
}

class ROL_Objective;

// ---------- ROL Includes ------------//

#include "ROL_StdVector.hpp"
#include "ROL_Vector.hpp"
#include "ROL_XyceVector.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <map>

// ------------ Xyce Includes -----------//

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_ROL.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_OutputMgrAdapter.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Solver.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Builder.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_NLS_Sensitivity.h>
#include <N_NLS_fwd.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace Analysis {

enum sensDiffMode
{
  SENS_FWD,
  SENS_REV,
  SENS_CNT,
  NUM_DIFF_MODES
};

//-------------------------------------------------------------------------
// Class         : ROL_DC
// Purpose       : Thin inheritance layer of DCSweep class for ROL
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 07/29/22
//-------------------------------------------------------------------------
class ROL_DC: public DCSweep
{
public:

  ROL_DC(
      AnalysisManager &analysis_manager, 
      Nonlinear::Manager &nonlinear_manager,
      Loader::Loader &loader, 
      Linear::System & linear_system,
      Topo::Topology & topology,
      IO::InitialConditionsManager & initial_conditions_manager)
  : DCSweep( analysis_manager, &linear_system, nonlinear_manager, loader, topology, initial_conditions_manager ),
    analysisManager_( analysis_manager ),
    nonlinearManager_( nonlinear_manager ),
    loader_( loader ),
    topology_( topology ),
    initialConditionsManager_( initial_conditions_manager ),
    linearSystem_( linear_system ),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    stepLoopSize_(0),
    numParams_(0)
  {}

  virtual ~ROL_DC() { doFree(); }

  void setSweepValue(int step);
  bool doAllocations(int nc, int nz);
  int  getLoopSize() { return dcLoopSize_; }
  std::map< std::string, std::vector<std::string> >& getDataNamesMap() { return dataNamesMap_; }
  std::map< std::string, std::vector< std::vector<double> > >& getDataTablesMap() { return dataTablesMap_; }

  bool setAnalysisParams(const std::vector<Util::OptionBlock>& paramsBlock);
  bool createObjectives(const std::vector<ROL_Objective>& objVec);

  using DCSweep::setTimeIntegratorOptions;

  using DCSweep::doFinish;
  using DCSweep::doLoopProcess;
  using DCSweep::doProcessFailedStep;
  using DCSweep::doHandlePredictor;

  bool doInit();
  bool doProcessSuccessfulStep();

  Teuchos::RCP<::ROL::Objective_SimOpt<RealT> > obj_;

  std::vector<Linear::Vector *>         solutionPtrVector_;
  std::vector<Linear::Vector *>         statePtrVector_;
  std::vector<Linear::Vector *>         constraintPtrVector_;
  std::vector<Linear::Vector *>         mydfdpPtrVector_;
  std::vector<Linear::Vector *>         mydqdpPtrVector_;
  std::vector<Linear::Vector *>         mydbdpPtrVector_;
  std::vector<Linear::Vector *>         mysensRHSPtrVector_;

protected:
  using DCSweep::doRun;

private:
  bool doFree();

  std::vector<ROL_Objective>            objVec_;

  AnalysisManager &                     analysisManager_;
  Nonlinear::Manager &                  nonlinearManager_; // TT
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  Linear::System &                      linearSystem_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  int                                   stepLoopSize_;
  int                                   numParams_;
};


// **************** EqualityConstraint class *************************************//



template<class Real>
class EqualityConstraint_ROL_DC : public ::ROL::Constraint_SimOpt<Real> 
{

 private:
  Real                                  nu_, nc_, nz_;
  Nonlinear::NonLinearSolver &          nls_;
  Linear::System &                      linearSystem_;
  Loader::NonlinearEquationLoader &     nEqLoader_;
  bool                                  forceFD_;
  std::vector<std::string>              paramNameVec_;
  Linear::Vector *                      origFVectorPtr_;
  Linear::Vector *                      pertFVectorPtr_;
  Linear::Vector *                      origQVectorPtr_;
  Linear::Vector *                      pertQVectorPtr_;
  Linear::Vector *                      origBVectorPtr_;
  Linear::Vector *                      pertBVectorPtr_;
  std::vector<Real>                     paramValue_;
  // TODO (asjavee): areParamsChanged saves us from unnecessary solves. It 
  //   should be obsolete once we remove "false" from the initialization of our
  //   reduced objective (currently line 396 of N_ANP_ROL.C).
  bool                                  areParamsChanged_;

  // HKT: Limit the parameter changes to when the parameters are actually changing.
  void update_2( const ::ROL::Vector<Real> &z, ::ROL::UpdateType type, int iter = -1 ) 
  {
    Teuchos::RCP<const std::vector<Real>> zp = (Teuchos::dyn_cast<::ROL::StdVector<Real>>(const_cast<::ROL::Vector<Real>&>(z))).getVector();

    std::string parameterName; 
    for (int i = 0; i < nz_; i++)
    {
      if ( paramValue_[i] != (*zp)[i] )
      {
        parameterName = paramNameVec_[i];
        paramValue_[i] = (*zp)[i];
        nEqLoader_.setParam(parameterName, paramValue_[i]);
        //std::cout << "setControlParams: " << parameterName << " " << paramValue_[i] << std::endl;
        areParamsChanged_ = true;
      }
    }
  }

  bool rhs_()
  {
    Stats::StatTop _residualStat("Residual");
    Stats::TimeBlock _residualTimer(_residualStat);

    nEqLoader_.loadRHS();
    ++numResidualLoads_;
    totalResidualLoadTime_ += nEqLoader_.getResidualTime();

    return true;
  }

  bool jacobian_()
  {
    Stats::StatTop _jacobianStat("Jacobian");
    Stats::TimeBlock _jacobianTimer(_jacobianStat);

    nEqLoader_.loadJacobian();
    ++numJacobianLoads_;
    totalJacobianLoadTime_ += nEqLoader_.getJacobianTime();
 
    return true;
  }

 protected:
  ROL_DC &                              rolSweep_;
  AnalysisManager &                     analysisManager_;
      
 public:
  int numResidualLoads_, numJacobianLoads_;
  double totalResidualLoadTime_, totalJacobianLoadTime_;

  EqualityConstraint_ROL_DC(
    Real                                nu, // # of solution variables
    Real                                nc, // # of constraint equations
    Real                                nz, // # of optimization variables
    AnalysisManager &                   analysisManager,
    Nonlinear::NonLinearSolver &        nls,
    Linear::System &                    linearSystem,
    std::vector<std::string>            paramNameVec,
    ROL_DC &                            rolSweep)
    : nu_(nu),
      nc_(nc),
      nz_(nz),
      forceFD_(false),
      analysisManager_(analysisManager),
      nEqLoader_(analysisManager.getNonlinearEquationLoader()),
      rolSweep_(rolSweep),
      nls_(nls),
      linearSystem_(linearSystem),
      paramNameVec_(paramNameVec),
      origFVectorPtr_(0),
      pertFVectorPtr_(0),
      origQVectorPtr_(0),
      pertQVectorPtr_(0),
      origBVectorPtr_(0),
      pertBVectorPtr_(0),
      paramValue_(nz, 0.0),
      numResidualLoads_(0),
      numJacobianLoads_(0),
      totalResidualLoadTime_(0.0),
      totalJacobianLoadTime_(0.0),
      areParamsChanged_(true)
  {}

  ~EqualityConstraint_ROL_DC() {}

  void value(::ROL::Vector<Real> &c, 
             const ::ROL::Vector<Real> &u,
             const ::ROL::Vector<Real> &z, 
             Real &tol)
  {
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > cp = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(c)).getVector());
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    
    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);  // Turns off voltage limiting
    
    bool success;
    for (int i = 0; i < nc_; i++) 
    {
      rolSweep_.setSweepValue(i);  // TT: Updates source term. Note: A more direct approach is to use setParam.
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[i]);
      rhs_();
      *((*cp)[i]) = *(linearSystem_.getRHSVector());
      // (*cp)[i]->print(Xyce::dout());
    }
  }

  void solve(::ROL::Vector<Real> &c, 
             ::ROL::Vector<Real> &u, 
             const ::ROL::Vector<Real> &z, 
             Real &tol)
  { 
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(u)).getVector();

    if (areParamsChanged_)
    {
      // Sweep to obtain the vector of solutions.
      bool success = rolSweep_.doInit() && rolSweep_.doLoopProcess() && rolSweep_.doFinish(); 
   
      // Write the vector of solutions into u.
      for (int i = 0; i < nc_ ; i++) {
        *((*up)[i]) = *(rolSweep_.solutionPtrVector_[i]);
        // (*up)[i]->print(Xyce::dout());
      }
    }
    areParamsChanged_ = false;
    
    value(c, u, z, tol);
  }

  void applyJacobian_1(::ROL::Vector<Real> &jv,
                      const ::ROL::Vector<Real> &v,
                      const ::ROL::Vector<Real> &u,
                      const ::ROL::Vector<Real> &z, 
                      Real &tol)
  {
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > jvp = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vp = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    
    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);// turns off voltage limiting

    bool success;  
    for (int i = 0; i < nc_; i++)
    {
      // Set next solution vector pointer to the values given in vector u[i] 
      // Linear::Vector * temp = (*up)[i].getRawPtr();
      // success = analysisManager_.getDataStore()->setNextSolVectorPtr(temp); //segfault

      // Using newly added function in DataStore class
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[i]);

      //rolSweep_.setSweepValue(i);// not needed here
      // Load rhs: appears to be necessary
      success = rhs_(); 
      
      // Load Jacobian
      success = jacobian_();
     
      // Get Jacobian
      Linear::Matrix & Jac = *(linearSystem_.getJacobianMatrix());
      
      // Apply Jacobian
      Jac.matvec(false, const_cast< const Linear::Vector& >( *((*vp)[i]) ), *((*jvp)[i]) );
      (*jvp)[i]->scale( -1.0 );
     
      // // Test
      // std::cout << "Jac_1 Success = " << success << std::endl;
      // (*jvp)[i]->print(Xyce::dout());
    }
    nEqLoader_.setVoltageLimiterStatus(vstatus);
  }

  void applyAdjointJacobian_1(::ROL::Vector<Real> &jv,
                              const ::ROL::Vector<Real> &v,
                              const ::ROL::Vector<Real> &u,
                              const ::ROL::Vector<Real> &z, 
                              Real &tol)
  {
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > jvp = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vp = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();

    bool success;

    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);// turns off voltage limiting

    for (int i=0;i<nc_;i++)
    {
      // Set next solution vector pointer to the values given in vector u[i] 
      // Linear::Vector * temp = (*up)[i].getRawPtr();
      // success = analysisManager_.getDataStore()->setNextSolVectorPtr(temp); //segfault

      // Using newly added function in DataStore class
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[i]);

      //rolSweep_.setSweepValue(i);// not needed here
      // Load rhs: appears to be necessary
      success = rhs_(); 
      
      // Load Jacobian
      success = jacobian_();
      
      // Get Jacobian
      Linear::Matrix & Jac = *(linearSystem_.getJacobianMatrix());

      // Apply Jacobian transpose
      Jac.matvec(true, const_cast< const Linear::Vector& >( *((*vp)[i]) ), *((*jvp)[i]) );
      (*jvp)[i]->scale( -1.0 );
      
    }
    nEqLoader_.setVoltageLimiterStatus(vstatus);
  }

  void applyInverseJacobian_1(::ROL::Vector<Real> &jv,
                              const ::ROL::Vector<Real> &v,
                              const ::ROL::Vector<Real> &u,
                              const ::ROL::Vector<Real> &z, 
                              Real &tol)
  {
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > jvp = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vp = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();

    bool success;

    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);// turns off voltage limiting

    // save rhs and Newton vectors
    Linear::Vector * savedRHSVectorPtr_ = nls_.getRHSVector().cloneCopyVector();
    Linear::Vector * savedNewtonVectorPtr_ = nls_.getNewtonVector().cloneCopyVector();

    for (int i=0;i<nc_;i++)
    {
      // Set next solution vector pointer to the values given in vector u[i] 
      // Linear::Vector * temp = (*up)[i].getRawPtr();
      // success = analysisManager_.getDataStore()->setNextSolVectorPtr(temp); //segfault

      // Using newly added function in DataStore class
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[i]);

      //rolSweep_.setSweepValue(i);// not needed here
      // Load rhs: appears to be necessary
      success = rhs_();
      
      // Load Jacobian
      success = jacobian_();
      
      // save rhs and Newton vectors
      // NOTE (HKT): this doesn't save a COPY of the vector, just copies the pointer.
      // When the vector is overwritten by vp[i], the saved vectors will be overwritten as well.
      // Restoring them below will have no impact, the vectors were never saved to begin with.
      // NOTE 2 (HKT): this save and restore vector work has been moved outside the loop.
      //Linear::Vector * savedRHSVectorPtr_ = nls_.rhsVectorPtr_;
      //Linear::Vector * savedNewtonVectorPtr_ = nls_.NewtonVectorPtr_;
      //savedRHSVectorPtr_->update(1.0, *(nls_.rhsVectorPtr_),0.0);
      //savedNewtonVectorPtr_->update(1.0, *(nls_.NewtonVectorPtr_),0.0);

      // Linear::Vector * temp1 = (*vp)[i].getRawPtr();
      nls_.setRHSVector( *(*vp)[i], -1.0 );  // solver expects negative rhs

      int status = nls_.getLinearSolver()->solve(false);
      if (status!=0)
      {
        std::string msg("Sensitivity::solveAdjoint.  Solver failed\n");
        Report::DevelFatal() << msg;
      }

      *((*jvp)[i]) = nls_.getNewtonVector();

    }

    // Restore the RHS and Newton vectors.
    nls_.setRHSVector( *savedRHSVectorPtr_ );
    nls_.setNewtonVector( *savedNewtonVectorPtr_ );
   
    // Delete saved vectors.        
    delete savedRHSVectorPtr_;
    delete savedNewtonVectorPtr_;
 
    nEqLoader_.setVoltageLimiterStatus(vstatus);
  }

  void applyInverseAdjointJacobian_1(::ROL::Vector<Real> &jv, 
                                   const ::ROL::Vector<Real> &v,
                                   const ::ROL::Vector<Real> &u,
                                   const ::ROL::Vector<Real> &z,
                                   Real &tol)
  {
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > jvp = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vp = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(z))).getVector();

    bool success;
    
    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);// turns off voltage limiting

    // save rhs and Newton vectors
    // [see notes in applyInverseJacobian_1]    
    Linear::Vector * savedRHSVectorPtr_ = nls_.getRHSVector().cloneCopyVector();
    Linear::Vector * savedNewtonVectorPtr_ = nls_.getNewtonVector().cloneCopyVector();

    for (int i=0;i<nc_;i++)
    {
      // Set next solution vector pointer to the values given in vector u[i] 
      // Linear::Vector * temp = (*up)[i].getRawPtr();
      // success = analysisManager_.getDataStore()->setNextSolVectorPtr(temp); //segfault

      // Using newly added function in DataStore class
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[i]);

      //rolSweep_.setSweepValue(i);// not needed here
      // Load rhs: appears to be necessary
      success = rhs_();

      // Load Jacobian
      success = jacobian_();
      
      // save rhs and Newton vectors
      // [see notes in applyInverseJacobian_1]    
      //Linear::Vector * savedRHSVectorPtr_ = nls_.rhsVectorPtr_;
      //Linear::Vector * savedNewtonVectorPtr_ = nls_.NewtonVectorPtr_;
      //savedRHSVectorPtr_->update(1.0, *(nls_.rhsVectorPtr_),0.0);
      //savedNewtonVectorPtr_->update(1.0, *(nls_.NewtonVectorPtr_),0.0);

      // Linear::Vector * temp1 = (*vp)[i].getRawPtr();
      // nls_.rhsVectorPtr_->update(-1.0, *temp1, 0.0);
      nls_.setRHSVector( *(*vp)[i], -1.0 );

      int status = nls_.getLinearSolver()->solveTranspose(false);// some reusefactors
      if (status!=0)
      {
        std::string msg("Sensitivity::solveAdjoint.  Solver failed\n");
        Report::DevelFatal() << msg;
      }
      
      *((*jvp)[i]) = nls_.getNewtonVector();
    }

    // Restore the RHS and Newton vectors.
    nls_.setRHSVector( *savedRHSVectorPtr_ );
    nls_.setNewtonVector( *savedNewtonVectorPtr_ );
   
    // Delete saved vectors.        
    delete savedRHSVectorPtr_;
    delete savedNewtonVectorPtr_;
 
    nEqLoader_.setVoltageLimiterStatus(true);
  }

  void applyJacobian_2(::ROL::Vector<Real> &jv,
                       const ::ROL::Vector<Real> &v,
                       const ::ROL::Vector<Real> &u,
                       const ::ROL::Vector<Real> &z,
                       Real &tol)
  {
    Teuchos::RCP<std::vector<Teuchos::RCP<Linear::Vector> > > jvp = Teuchos::rcp_const_cast<std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(jv)).getVector());
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > vp = (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(v))).getVector();
    
    bool success;    
    jv.zero();
    
    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);// turns off voltage limiting
    
    // loop over constraint equations
    for (int k=0;k<nc_;k++)
    {
      rolSweep_.setSweepValue(k); // update source term; Important!
      
      // Note: perhaps all source terms need be updated
      
      // Set next solution vector pointer to the values given in vector u[i] 
      // Linear::Vector * temp = (*up)[k].getRawPtr();
      // success = analysisManager_.getDataStore()->setNextSolVectorPtr(temp); //segfault
      
      // Using newly added function in DataStore class
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[k]);
      
      // Load rhs
      rhs_();
      
      int iparam;
      std::string msg;
      
      analysisManager_.getDataStore()->paramOrigVals_.clear();
      
      // it is necessary to load the Jacobian here to make sure we have the most
      // up-to-date matrix.  The Jacobian is not loaded for the final 
      // evaluation of the residual in the Newton solve.
      success = jacobian_();
      
      // Loop over the vector of parameters.  For each parameter, find the
      // device entity (a model or an instance) which corresponds to it, and
      // perform the finite difference calculation.
      std::vector<std::string>::iterator firstParam = paramNameVec_.begin ();
      std::vector<std::string>::iterator lastParam  = paramNameVec_.end ();
      std::vector<std::string>::iterator iterParam;
      for ( iterParam=firstParam, iparam=0;
            iterParam!=lastParam; ++iterParam, ++iparam )
      {
        std::string paramName(*iterParam);
        //std::cout << "Computing derivative wrt " << paramNameVec_[iparam] << std::endl;  
        // get the original value of this parameter, to be used for scaling and/or 
        // numerical derivatives later.
        double paramOrig = 0.0;
        bool found = nEqLoader_.getParamAndReduce(analysisManager_.getPDSManager()->getPDSComm()->comm(),paramName, paramOrig);
        
        if (!found)
        {
          std::string msg("EqualityConstraint::applyJacobian_2: cannot find parameter ");
          msg += paramName;
          Report::DevelFatal() << msg;
        }
        
        analysisManager_.getDataStore()->paramOrigVals_.push_back(paramOrig);
        
        // check if derivative is available analytically.  If not, take FD.
        bool analyticAvailable = false;
        if (!forceFD_)
        {
          analyticAvailable = nEqLoader_.analyticSensitivitiesAvailable (paramName);
        }
        if (analyticAvailable)
        {
          //std::cout << "Analytic derivs available" << std::endl;
          std::vector<double> dfdpVec;
          std::vector<double> dqdpVec;
          std::vector<double> dbdpVec;
          
          std::vector<int> FindicesVec;
          std::vector<int> QindicesVec;
          std::vector<int> BindicesVec;
          
          nEqLoader_.getAnalyticSensitivities(paramName,dfdpVec,dqdpVec,dbdpVec,FindicesVec, QindicesVec, BindicesVec);
          
          rolSweep_.mydfdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydqdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydbdpPtrVector_[iparam]->putScalar(0.0);
          
          int Fsize=FindicesVec.size();
          for (int i=0;i<Fsize;++i)
          {
            Linear::Vector & dfdpRef = *(rolSweep_.mydfdpPtrVector_[iparam]);
            dfdpRef[FindicesVec[i]]  += dfdpVec[i];
          }
          
          int Qsize=QindicesVec.size();
          for (int i=0;i<Qsize;++i)
          {
            Linear::Vector & dqdpRef = *(rolSweep_.mydqdpPtrVector_[iparam]);
            dqdpRef[QindicesVec[i]]  += dqdpVec[i];
          }
          
          int Bsize=BindicesVec.size();
          for (int i=0;i<Bsize;++i)
          {
            Linear::Vector & dbdpRef = *(rolSweep_.mydbdpPtrVector_[iparam]);
            dbdpRef[BindicesVec[i]]  += dbdpVec[i];
          }
          
          rolSweep_.mydfdpPtrVector_[iparam]->fillComplete();
          rolSweep_.mydqdpPtrVector_[iparam]->fillComplete();
          rolSweep_.mydbdpPtrVector_[iparam]->fillComplete();
          
          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            Xyce::dout() << *iterParam << ": ";
            Xyce::dout().setf(std::ios::scientific);
            Xyce::dout() << std::endl;
            
            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {
              Xyce::dout() 
                <<"dfdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydfdpPtrVector_[iparam]))[k1]
                <<" dqdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydqdpPtrVector_[iparam]))[k1]
                <<" dbdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydbdpPtrVector_[iparam]))[k1]
                <<std::endl;
            }
          }
        }
        else
        {
          //Xyce::dout() << "Analytical derivatives NOT available" << std::endl;
          
          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            Xyce::dout() << std::endl << "  Calculating numerical df/dp, dq/dp and db/dp for: ";
            Xyce::dout() << *iterParam << std::endl;
          }

          // save a copy of the DAE vectors
          delete origFVectorPtr_;
          origFVectorPtr_ = analysisManager_.getDataStore()->daeFVectorPtr->cloneCopyVector();
          pertFVectorPtr_ = linearSystem_.builder().createVector();

          delete origQVectorPtr_;
          origQVectorPtr_ = analysisManager_.getDataStore()->daeQVectorPtr->cloneCopyVector();
          pertQVectorPtr_ = linearSystem_.builder().createVector();

          delete origBVectorPtr_;
          origBVectorPtr_ = analysisManager_.getDataStore()->daeBVectorPtr->cloneCopyVector();
          pertBVectorPtr_ = linearSystem_.builder().createVector();
          
          // now perturb the value of this parameter.
          double sqrtEta_ = pow(10,int(log10(fabs(paramOrig))));; // TT: this should be parameter specific (see Sensitivity class)
          double dp = sqrtEta_ * 1.e-5 * (1.0 + fabs(paramOrig));
          double paramPerturbed = paramOrig;
          
          int difference=SENS_FWD; // TT: this should be optional
          if (difference==SENS_FWD)
          {
            paramPerturbed += dp;
          }
          else if (difference==SENS_REV)
          {
            paramPerturbed -= dp;
          }
          else if (difference==SENS_CNT)
          {
            static std::string tmp = "difference=central not supported.\n";
            Report::DevelFatal() << tmp;
          }
          else
          {
            static std::string tmp = "difference not recognized!\n";
            Report::UserFatal0() << tmp;
          }
          
          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            int maxParamStringSize_ = 16; // TT: not sure what this should be
            Xyce::dout() << std::setw(maxParamStringSize_)<< *iterParam
                        << " dp = " << std::setw(11)<< std::scientific<< std::setprecision(4) << dp 
                        << " original value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramOrig 
                        << " modified value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramPerturbed 
                        <<std::endl;
          }
          
          nEqLoader_.setParam (paramName, paramPerturbed);

          // Now that the parameter has been perturbed,
          // calculate the numerical derivative.

          rolSweep_.setSweepValue(k);// needed?
          // Load F,Q and B.
          rhs_();

          // save the perturbed DAE vectors
          (*pertFVectorPtr_) = *(analysisManager_.getDataStore()->daeFVectorPtr);
          (*pertQVectorPtr_) = *(analysisManager_.getDataStore()->daeQVectorPtr);
          (*pertBVectorPtr_) = *(analysisManager_.getDataStore()->daeBVectorPtr);

          // calculate the df/dp vector.  
          double rdp=1/dp;
          rolSweep_.mydfdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydfdpPtrVector_[iparam]->update (+1.0, *(pertFVectorPtr_)); //+Fperturb
          rolSweep_.mydfdpPtrVector_[iparam]->update (-1.0, *(origFVectorPtr_)); //-Forig
          rolSweep_.mydfdpPtrVector_[iparam]->scale(rdp);

          // calculate the dq/dp vector.  
          rolSweep_.mydqdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydqdpPtrVector_[iparam]->update (+1.0, *(pertQVectorPtr_)); //+Fperturb
          rolSweep_.mydqdpPtrVector_[iparam]->update (-1.0, *(origQVectorPtr_)); //-Forig
          rolSweep_.mydqdpPtrVector_[iparam]->scale(rdp);

          // calculate the db/dp vector.  
          rolSweep_.mydbdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydbdpPtrVector_[iparam]->update (+1.0, *(pertBVectorPtr_)); //+Fperturb
          rolSweep_.mydbdpPtrVector_[iparam]->update (-1.0, *(origBVectorPtr_)); //-Forig
          rolSweep_.mydbdpPtrVector_[iparam]->scale(rdp);

          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            Xyce::dout() << *iterParam << ": ";
            Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
            Xyce::dout() << "deviceSens_dp = " << dp << std::endl;

            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {

              Xyce::dout() 
                <<"fpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertFVectorPtr_))[k1]
                <<" forig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origFVectorPtr_))[k1]
                <<" dfdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydfdpPtrVector_[iparam]))[k1]
                <<std::endl;
            }

            Xyce::dout() << std::endl;
            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {
              Xyce::dout() 
                <<"qpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertQVectorPtr_))[k1]
                <<" qorig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origQVectorPtr_))[k1]
                <<" dqdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydqdpPtrVector_[iparam]))[k1]
                <<std::endl;
            }

            Xyce::dout() << std::endl ;
            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {
              Xyce::dout() 
                <<"bpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertBVectorPtr_))[k1]
                <<" borig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origBVectorPtr_))[k1]
                <<" dbdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydbdpPtrVector_[iparam]))[k1]
                <<std::endl;

            }

            // std::ostringstream filename; 
            // filename << netlistFilename_ << "_dfdp";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // dfdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_fpert";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // pertFVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_dqdp";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // dqdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_qpert";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // pertQVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_dbdp";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // dbdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));
  
            // filename.str("");
            // filename << netlistFilename_ << "_bpert";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // pertBVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));
          }

          // now reset the parameter and rhs to previous values.
          nEqLoader_.setParam (paramName, paramOrig);

          analysisManager_.getDataStore()->daeFVectorPtr->update(1.0, *(origFVectorPtr_), 0.0);
          analysisManager_.getDataStore()->daeQVectorPtr->update(1.0, *(origQVectorPtr_), 0.0);
          analysisManager_.getDataStore()->daeBVectorPtr->update(1.0, *(origBVectorPtr_), 0.0);
        }

        // Now collect dfdp, dqdp, dbdp into dFdp (as done in obtainSensitivityResiduals in NoTimeIntegration, that is we ignore dqdp)

        Linear::Vector & sensRHS = *(rolSweep_.mysensRHSPtrVector_[iparam]);
        Linear::Vector & dfdpRef = *(rolSweep_.mydfdpPtrVector_[iparam]);
        Linear::Vector & dbdpRef = *(rolSweep_.mydbdpPtrVector_[iparam]);

        //sensRHS.linearCombo(0.0,sensRHS,+1.0,dfdpRef); saving old code
        sensRHS.update(+1.0,dfdpRef,0.0);
        //sensRHS.linearCombo(1.0,sensRHS,-1.0,dbdpRef); saving old code
        sensRHS.update(-1.0,dbdpRef);

        sensRHS.scale(-1.0);

        // multiply by v[iparam]
        sensRHS.scale( (*vp)[iparam] );

        // add to jvp[k]
        (*jvp)[k]->update(1.0,sensRHS);

      }// end of parameters loop
        
    }// end of constraints loop
    nEqLoader_.setVoltageLimiterStatus(vstatus);
  }

  void applyAdjointJacobian_2(::ROL::Vector<Real> &jv,
                              const ::ROL::Vector<Real> &v,
                              const ::ROL::Vector<Real> &u,
                              const ::ROL::Vector<Real> &z,
                              Real &tol)
  {
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vp = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<std::vector<Real> > jvp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(jv)).getVector());

    bool success;    
    jv.zero();
    
    bool vstatus = nEqLoader_.getVoltageLimiterStatus();
    nEqLoader_.setVoltageLimiterStatus(false);// turns off voltage limiting

    // loop over constraint equations
    for (int k=0;k<nc_;k++)
    {
      rolSweep_.setSweepValue(k); // update source term; Important!
      
      // Set next solution vector pointer to the values given in vector u[i] 
      // Linear::Vector * temp = (*up)[k].getRawPtr();
      // success = analysisManager_.getDataStore()->setNextSolVectorPtr(temp); //segfault

      // Using newly added function in DataStore class
      success = analysisManager_.getDataStore()->setNextSolVectorPtr(*(*up)[k]);

      // Load rhs
      success = rhs_();
      
      int iparam;
      std::string msg;

      analysisManager_.getDataStore()->paramOrigVals_.clear();

      // it is necessary to load the Jacobian here to make sure we have the most
      // up-to-date matrix.  The Jacobian is not loaded for the final 
      // evaluation of the residual in the Newton solve.
      jacobian_();

      // Loop over the vector of parameters.  For each parameter, find the
      // device entity (a model or an instance) which corresponds to it, and
      // perform the finite difference calculation.
      std::vector<std::string>::iterator firstParam = paramNameVec_.begin ();
      std::vector<std::string>::iterator lastParam  = paramNameVec_.end ();
      std::vector<std::string>::iterator iterParam;
      for ( iterParam=firstParam, iparam=0;
            iterParam!=lastParam; ++iterParam, ++iparam )
      {
        std::string paramName(*iterParam);
        // std::cout << "Computing derivative wrt " << paramNameVec_[iparam] << std::endl;  
        // get the original value of this parameter, to be used for scaling and/or 
        // numerical derivatives later.
        double paramOrig = 0.0;
        bool found = nEqLoader_.getParamAndReduce(analysisManager_.getPDSManager()->getPDSComm()->comm(),paramName, paramOrig);

        if (!found)
        {
          std::string msg("EqualityConstraint::applyJacobian_2: cannot find parameter ");
          msg += paramName;
          Report::DevelFatal() << msg;
        }

        analysisManager_.getDataStore()->paramOrigVals_.push_back(paramOrig);

        // check if derivative is available analytically.  If not, take FD.
        bool analyticAvailable = false;
        if (!forceFD_)
        {
          analyticAvailable = nEqLoader_.analyticSensitivitiesAvailable (paramName);
        }
        if (analyticAvailable)
        {
          std::vector<double> dfdpVec;
          std::vector<double> dqdpVec;
          std::vector<double> dbdpVec;

          std::vector<int> FindicesVec;
          std::vector<int> QindicesVec;
          std::vector<int> BindicesVec;

          nEqLoader_.getAnalyticSensitivities(paramName,dfdpVec,dqdpVec,dbdpVec,FindicesVec, QindicesVec, BindicesVec);

          rolSweep_.mydfdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydqdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydbdpPtrVector_[iparam]->putScalar(0.0);

          int Fsize=FindicesVec.size();
          for (int i=0;i<Fsize;++i)
          {
            Linear::Vector & dfdpRef = *(rolSweep_.mydfdpPtrVector_[iparam]);
            dfdpRef[FindicesVec[i]]  += dfdpVec[i];
          }

          int Qsize=QindicesVec.size();
          for (int i=0;i<Qsize;++i)
          {
            Linear::Vector & dqdpRef = *(rolSweep_.mydqdpPtrVector_[iparam]);
            dqdpRef[QindicesVec[i]]  += dqdpVec[i];
          }

          int Bsize=BindicesVec.size();
          for (int i=0;i<Bsize;++i)
          {
            Linear::Vector & dbdpRef = *(rolSweep_.mydbdpPtrVector_[iparam]);
            dbdpRef[BindicesVec[i]]  += dbdpVec[i];
          }

          rolSweep_.mydfdpPtrVector_[iparam]->fillComplete();
          rolSweep_.mydqdpPtrVector_[iparam]->fillComplete();
          rolSweep_.mydbdpPtrVector_[iparam]->fillComplete();

          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            Xyce::dout() << *iterParam << ": ";
            Xyce::dout().setf(std::ios::scientific);
            Xyce::dout() << std::endl;

            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {
              Xyce::dout() 
                <<"dfdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydfdpPtrVector_[iparam]))[k1]
                <<" dqdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydqdpPtrVector_[iparam]))[k1]
                <<" dbdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydbdpPtrVector_[iparam]))[k1]
                <<std::endl;
            }
          }
        }
        else
        {
          // std::cout << "Analytical derivatives NOT available" << std::endl;

          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            Xyce::dout() << std::endl << "  Calculating numerical df/dp, dq/dp and db/dp for: ";
            Xyce::dout() << *iterParam << std::endl;
          }

          // save a copy of the DAE vectors
          origFVectorPtr_->update(1.0, *(analysisManager_.getDataStore()->daeFVectorPtr), 0.0);
          origQVectorPtr_->update(1.0, *(analysisManager_.getDataStore()->daeQVectorPtr), 0.0);
          origBVectorPtr_->update(1.0, *(analysisManager_.getDataStore()->daeBVectorPtr), 0.0);

          // now perturb the value of this parameter.
          double sqrtEta_ = pow(10,int(log10(fabs(paramOrig))));; // TT: this should be parameter specific
          double dp = sqrtEta_ * 1.e-4 * (1.0 + fabs(paramOrig));
          double paramPerturbed = paramOrig;

          int difference=SENS_FWD; // TT: this should be optional
          if (difference==SENS_FWD)
          {
            paramPerturbed += dp;
          }
          else if (difference==SENS_REV)
          {
            paramPerturbed -= dp;
          }
          else if (difference==SENS_CNT)
          {
            static std::string tmp = "difference=central not supported.\n";
            Report::UserFatal0() << tmp;
          }
          else
          {
            static std::string tmp = "difference not recognized!\n";
            Report::UserFatal0() << tmp;
          }

          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            int maxParamStringSize_ = 16; // TT: not sure what this should be
            Xyce::dout() << std::setw(maxParamStringSize_)<< *iterParam
                         << " dp = " << std::setw(11)<< std::scientific<< std::setprecision(4) << dp 
                         << " original value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramOrig 
                         << " modified value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramPerturbed 
                         <<std::endl;
          }

          nEqLoader_.setParam (paramName, paramPerturbed);

          // Now that the parameter has been perturbed,
          // calculate the numerical derivative.

          rolSweep_.setSweepValue(k);// needed?
          // Load F,Q and B.
          rhs_();

          // save the perturbed DAE vectors
          pertFVectorPtr_->update(1.0, *(analysisManager_.getDataStore()->daeFVectorPtr), 0.0);
          pertQVectorPtr_->update(1.0, *(analysisManager_.getDataStore()->daeQVectorPtr), 0.0);
          pertBVectorPtr_->update(1.0, *(analysisManager_.getDataStore()->daeBVectorPtr), 0.0);

          // calculate the df/dp vector.  
          double rdp=1/dp;
          rolSweep_.mydfdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydfdpPtrVector_[iparam]->update (+1.0, *(pertFVectorPtr_)); //+Fperturb
          rolSweep_.mydfdpPtrVector_[iparam]->update (-1.0, *(origFVectorPtr_)); //-Forig
          rolSweep_.mydfdpPtrVector_[iparam]->scale(rdp);

          // calculate the dq/dp vector.  
          rolSweep_.mydqdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydqdpPtrVector_[iparam]->update (+1.0, *(pertQVectorPtr_)); //+Fperturb
          rolSweep_.mydqdpPtrVector_[iparam]->update (-1.0, *(origQVectorPtr_)); //-Forig
          rolSweep_.mydqdpPtrVector_[iparam]->scale(rdp);

          // calculate the db/dp vector.  
          rolSweep_.mydbdpPtrVector_[iparam]->putScalar(0.0);
          rolSweep_.mydbdpPtrVector_[iparam]->update (+1.0, *(pertBVectorPtr_)); //+Fperturb
          rolSweep_.mydbdpPtrVector_[iparam]->update (-1.0, *(origBVectorPtr_)); //-Forig
          rolSweep_.mydbdpPtrVector_[iparam]->scale(rdp);

          if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
          {
            Xyce::dout() << *iterParam << ": ";
            Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
            Xyce::dout() << "deviceSens_dp = " << dp << std::endl;

            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {

              Xyce::dout() 
                <<"fpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertFVectorPtr_))[k1]
                <<" forig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origFVectorPtr_))[k1]
                <<" dfdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydfdpPtrVector_[iparam]))[k1]
                <<std::endl;
            }

            Xyce::dout() << std::endl;
            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {
              Xyce::dout() 
                <<"qpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertQVectorPtr_))[k1]
                <<" qorig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origQVectorPtr_))[k1]
                <<" dqdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydqdpPtrVector_[iparam]))[k1]
                <<std::endl;
            }

            Xyce::dout() << std::endl ;
            for (int k1 = 0; k1 < analysisManager_.getDataStore()->solutionSize; ++k1)
            {
              Xyce::dout() 
                <<"bpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(pertBVectorPtr_))[k1]
                <<" borig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(origBVectorPtr_))[k1]
                <<" dbdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*(rolSweep_.mydbdpPtrVector_[iparam]))[k1]
                <<std::endl;

            }

            // std::ostringstream filename; 
            // filename << netlistFilename_ << "_dfdp";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // dfdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));
 
            // filename.str("");
            // filename << netlistFilename_ << "_fpert";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // pertFVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_dqdp";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // dqdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_qpert";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // pertQVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_dbdp";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // dbdpPtrVector_[iparam]->writeToFile(const_cast<char *>(filename.str().c_str()));

            // filename.str("");
            // filename << netlistFilename_ << "_bpert";
            // filename << std::setw(3) << std::setfill('0') << iparam;
            // filename << ".txt";
            // pertBVectorPtr_->writeToFile(const_cast<char *>(filename.str().c_str()));
          }

          // now reset the parameter and rhs to previous values.
          nEqLoader_.setParam (paramName, paramOrig);

          (*analysisManager_.getDataStore()->daeFVectorPtr) = *(origFVectorPtr_);
          (*analysisManager_.getDataStore()->daeQVectorPtr) = *(origQVectorPtr_);
          (*analysisManager_.getDataStore()->daeBVectorPtr) = *(origBVectorPtr_);
        }

        // Now collect dfdp, dqdp, dbdp into dFdp (as done in obtainSensitivityResiduals in NoTimeIntegration, that is we ignore dqdp)

        Linear::Vector & sensRHS = *(rolSweep_.mysensRHSPtrVector_[iparam]);
        Linear::Vector & dfdpRef = *(rolSweep_.mydfdpPtrVector_[iparam]);
        Linear::Vector & dbdpRef = *(rolSweep_.mydbdpPtrVector_[iparam]);

        // Compute sensRHS = dfdpRef - dbdpRef 
        sensRHS.update(+1.0, dfdpRef, -1.0, dbdpRef, 0.0);
        //sensRHS.linearCombo(0.0,sensRHS,+1.0,dfdpRef); saving old code
        //sensRHS.update(+1.0,dfdpRef,0.0);
        //sensRHS.linearCombo(1.0,sensRHS,-1.0,dbdpRef); saving old code
        //sensRHS.update(-1.0,dbdpRef);

        sensRHS.scale(-1.0);  // Why not integrate scaling in update above?

        // Take dot product with k-th vector in v
        Real dot = sensRHS.dotProduct(*((*vp)[k]));

        // Add to jvp
        (*jvp)[iparam] += dot;

      }// end of parameters loop
       
    }// end of constraints loop
    nEqLoader_.setVoltageLimiterStatus(vstatus);
  }

#if FD_HESSIAN==0  
  void applyAdjointHessian_11(::ROL::Vector<Real> &ahwv, const ::ROL::Vector<Real> &w, const ::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &u, const ::ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }
  void applyAdjointHessian_12(::ROL::Vector<Real> &ahwv, const ::ROL::Vector<Real> &w, const ::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &u, const ::ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }
  void applyAdjointHessian_21(::ROL::Vector<Real> &ahwv, const ::ROL::Vector<Real> &w, const ::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &u, const ::ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }
  void applyAdjointHessian_22(::ROL::Vector<Real> &ahwv, const ::ROL::Vector<Real> &w, const ::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &u, const ::ROL::Vector<Real> &z, Real &tol)
  {
    ahwv.zero();
  }
#endif

#if IDENTITY_PRECONDITIONER==1
  void applyPreconditioner(::ROL::Vector<Real> &pv, const ::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &x, const::ROL::Vector<Real> &g, Real &tol)
  {
    pv.set(v.dual());
  }
#endif

};


template<class Real>
class Objective_DC_L2Norm : public ::ROL::Objective_SimOpt<Real>
{
private:
  Real alpha_; // Penalty Parameter
  int ns_,nz_;
  Teuchos::RCP<std::vector<Real> > Imeas_;
  int solindex_;

public: 
  
  Objective_DC_L2Norm(Real alpha, int ns, int nz, std::vector<ROL_Objective_Arg<Real> >& objArgs)
  {
    ns_ = ns; // number of sources (same as number of constraints)
    nz_ = nz; // number of optimization variables
    alpha_ = alpha;
    Imeas_ = Teuchos::rcp(new std::vector<Real>(ns_,0.0));

    for ( int iObj = 0; iObj < objArgs.size(); ++iObj )
    {
      // Check if the argument points at data from .DATA or a solution index
      if (objArgs[iObj].objDataPtr != 0)
      {
        int index = objArgs[iObj].objIndex;
        std::vector< std::vector< double > > & dataPtr = *(objArgs[iObj].objDataPtr);
        if ( dataPtr.size() != ns_ )  
          Report::UserError0() << "ROL DC data objective is not the same dimension as the number of constraints: " << ns_;
        for (int i=0; i<ns_; ++i)
        {
          (*Imeas_)[i] = dataPtr[i][index];
        }
      }
      else
      {
        solindex_ = objArgs[iObj].objIndex;
      }
    }
  }
    
  Real value(const ::ROL::Vector<Real> &u, 
             const ::ROL::Vector<Real> &z,
             Real &tol)
  {
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    Real val = 0.0;
    Real current;
    // Here we need to know which index corresponds to the state variable we need (e.g., current);
    // solution index (solindex) is specified in constructor
    for (int i=0;i<ns_;i++)
    {
      current = (*(*up)[i])[solindex_]; // vector
      val += (current - (*Imeas_)[i])*(current - (*Imeas_)[i]);
    }
    return 0.5*val;
  }
  
  void gradient_1(::ROL::Vector<Real> &g,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    // Unwrap g
    Teuchos::RCP< std::vector<Teuchos::RCP<Linear::Vector> > > gup = Teuchos::rcp_const_cast< std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT U
    g.zero();
    Real current;
    for (int i=0; i<ns_; i++)
    {
      current = (*(*up)[i])[solindex_]; // vector
      (*(*gup)[i])[solindex_] = (current-(*Imeas_)[i]);
    }
    
  }

  void gradient_2(::ROL::Vector<Real> &g,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    // Unwrap g
    Teuchos::RCP<std::vector<Real> > gzp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<const ::ROL::StdVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT Z
    for (int i=0; i<nz_; i++)
    {
      (*gzp)[i] = 0.0;
    }
  }

  void hessVec_11(::ROL::Vector<Real> &hv, 
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    Teuchos::RCP< std::vector<Teuchos::RCP<Linear::Vector> > > hvup = Teuchos::rcp_const_cast< std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(hv)).getVector());
    // Unwrap v
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vup = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    
    hv.zero();
    // COMPUTE GRADIENT WRT U
    for(int i=0; i<ns_; i++)
    {
      (*(*hvup)[i])[solindex_] = (*(*vup)[i])[solindex_]; // picking the same component of vup[i] that corresponds to current in up[i]
    }
  }

  void hessVec_12(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }
  
  void hessVec_21(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }
  
  void hessVec_22(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }

};

template<class Real>
class Objective_DC_AMP : virtual public ::ROL::Objective_SimOpt<Real>
{
private:
  int ns_,nz_,first_,middle_,last_;
  int solindex_;
  
public: 
  
  Objective_DC_AMP(int ns, int nz)
  {
    ns_ = ns; // number of sources (same as number of constraints)
    nz_ = nz; // number of optimization variables
    first_ = 0;
    middle_ = ns_/2;
    last_ = ns_-1;
    solindex_ = 3; // V(3)
  }
    
  Real value(const ::ROL::Vector<Real> &u,
             const ::ROL::Vector<Real> &z,
             Real &tol)
  {
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    Real val = 0.0;
    Real V3,t;
    Real V3_0  = (*(*up)[first_ ])[solindex_];
    Real V3_10 = (*(*up)[middle_])[solindex_];
    Real V3_20 = (*(*up)[last_  ])[solindex_];
    Real A = (V3_0 - V3_20)/2.0;
    std::cout << "A = " << A << ' ' << (*zp)[0] << ' ' << (*zp)[1] << std::endl;
    for (int i=0;i<ns_;i++)
    {
      V3 = (*(*up)[i])[solindex_];
      t = -1.0 + i*(2.0/(ns_-1));
      val += ( A*t - V3 + V3_10 )*( A*t - V3 + V3_10 );
    }
    return 0.5*val;
  }
  
  void gradient_1(::ROL::Vector<Real> &g,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    // Unwrap g
    Teuchos::RCP< std::vector<Teuchos::RCP<Linear::Vector> > > gup = Teuchos::rcp_const_cast< std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT U
    g.zero();
    Real V3,t;
    Real V3_0  = (*(*up)[first_ ])[solindex_];
    Real V3_10 = (*(*up)[middle_])[solindex_];
    Real V3_20 = (*(*up)[last_  ])[solindex_];
    Real A = (V3_0 - V3_20)/2.0;
    (*(*gup)[first_ ])[solindex_] += ( -A - V3_0 + V3_10  )*(-3.0/2.0) + ( A - V3_20 + V3_10 )/2.0;
    (*(*gup)[last_  ])[solindex_] += (  A - V3_20 + V3_10 )*(-3.0/2.0) + ( -A - V3_0 + V3_10 )/2.0;
    (*(*gup)[middle_])[solindex_] += 2*V3_10 - V3_0 - V3_20;
    for (int i=1; i<ns_-1; i++)
    {
      V3 = (*(*up)[i])[solindex_];
      t = -1.0 + i*(2.0/(ns_-1));
      (*(*gup)[i ])[solindex_] -= ( A*t - V3 + V3_10 );
      (*(*gup)[first_ ])[solindex_] += ( A*t - V3 + V3_10 )*(t/2.0);
      (*(*gup)[last_  ])[solindex_] -= ( A*t - V3 + V3_10 )*(t/2.0);
      (*(*gup)[middle_])[solindex_] += ( A*t - V3 + V3_10 );
    }
  }

  void gradient_2(::ROL::Vector<Real> &g,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    // Unwrap g
    Teuchos::RCP<std::vector<Real> > gzp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<const ::ROL::StdVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT Z
    for (int i=0; i<nz_; i++)
    {
      (*gzp)[i] = 0.0;
    }
  }

  void hessVec_11(::ROL::Vector<Real>& hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    Teuchos::RCP< std::vector<Teuchos::RCP<Linear::Vector> > > hvup = Teuchos::rcp_const_cast< std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(hv)).getVector());
    // Unwrap v
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vup = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();    
    hv.zero();
    Real V3,t;
    Real V3_0  = (*(*vup)[first_ ])[solindex_];
    Real V3_10 = (*(*vup)[middle_])[solindex_];
    Real V3_20 = (*(*vup)[last_  ])[solindex_];
    Real sum1 = 0.0;
    for(int i=1;i<ns_-1;i++)
    {
      sum1 += std::pow((-1.0 + i*2.0/(ns_-1)),2);
    }
    sum1 /= 4.0;
    (*(*hvup)[first_ ])[solindex_] +=  (sum1+5.0/2.0)*V3_0 - V3_10 - (sum1+3.0/2.0)*V3_20;
    (*(*hvup)[middle_])[solindex_] = - V3_0 - V3_20;
    (*(*hvup)[last_  ])[solindex_] += -(sum1+3.0/2.0)*V3_0 - V3_10 + (sum1+5.0/2.0)*V3_20;
    for(int i=1; i<ns_-1; i++)
    {
      t = -1.0 + i*2.0/(ns_-1);
      V3 = (*(*vup)[i])[solindex_];
      (*(*hvup)[i])[solindex_] += V3
        + (-t/2.0)*V3_0
        -          V3_10
        + ( t/2.0)*V3_20;
      (*(*hvup)[first_ ])[solindex_] -= (t/2.0)*V3;
      (*(*hvup)[middle_])[solindex_] += (i==middle_) ? (ns_-1)*V3 : -V3;
      (*(*hvup)[last_  ])[solindex_] += (t/2.0)*V3;
    }
  }

  void hessVec_12(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }
  
  void hessVec_21(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }
  
  void hessVec_22(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }

};

template<class Real>
class Penalty_DC_AMP : virtual public ::ROL::Objective_SimOpt<Real>
{
private:
  int ptype_; // penalty type
  Real alpha_; // penalty parameter
  Real epsilon_; // penalty parameter
  Real ampl_; // amplification parameter
  int ns_,nz_,first_,middle_,last_;
  int solindex_;
  
public: 
  
  Penalty_DC_AMP(int ptype, Real alpha, Real ampl, int ns, int nz)
  {
    ns_ = ns; // number of sources (same as number of constraints)
    nz_ = nz; // number of optimization variables
    first_ = 0;
    middle_ = ns_/2;
    last_ = ns_-1;
    ptype_ = ptype;
    alpha_ = alpha;
    solindex_ = 3; // V(3)
    epsilon_ = 10.0;
    ampl_ = ampl;
  }

  void switch_ptype(int ptype)
  {
    ptype_ = ptype;
  }
    
  Real value(const ::ROL::Vector<Real> &u,
             const ::ROL::Vector<Real> &z,
             Real &tol)
  {
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    Real val = 0;
    Real A = ((*(*up)[first_])[solindex_] - (*(*up)[last_])[solindex_] )/2.0; // negative gain
    if(ptype_==1)
      val += alpha_*(0.5+0.5*tanh(epsilon_*(ampl_ - A*A)));
    if(ptype_==2)    
      val += alpha_*( (ampl_ - A*A > 0) ? ampl_-A*A : 0 );
    if(ptype_==3)
      val += alpha_*std::pow(( (ampl_ - A*A > 0) ? ampl_-A*A : 0 ),2);
    if(ptype_==4)
      val += alpha_*std::pow(( (A - ampl_ > 0) ? A - ampl_ : 0 ),2); // Moreau-Yosida penalty
    if(ptype_==5)
      val += alpha_*2.*A; // maximize gain
    return 0.5*val;
  }
  
  void gradient_1(::ROL::Vector<Real> &g,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    // Unwrap g
    Teuchos::RCP< std::vector<Teuchos::RCP<Linear::Vector> > > gup = Teuchos::rcp_const_cast< std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp =
      (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT U
    g.zero();
    Real A = ((*(*up)[first_])[solindex_] - (*(*up)[last_])[solindex_] )/2.0;
    if(ptype_==1)
    {
      (*(*gup)[first_])[solindex_] -= 0.5*(0.5*alpha_*epsilon_*A)/std::pow(cosh(epsilon_*(ampl_-A*A)),2);
      (*(*gup)[last_ ])[solindex_] += 0.5*(0.5*alpha_*epsilon_*A)/std::pow(cosh(epsilon_*(ampl_-A*A)),2);
    }

    if(ptype_==2)
    {    
      (*(*gup)[first_])[solindex_] -= 0.5*alpha_*( (ampl_ - A*A > 0) ? A : 0 );
      (*(*gup)[last_ ])[solindex_] += 0.5*alpha_*( (ampl_ - A*A > 0) ? A : 0 );
    }

    if(ptype_==3)
    {
      (*(*gup)[first_])[solindex_] -= 0.5*alpha_*( (ampl_ - A*A > 0) ? 8.0*A - 2.0*A*A*A : 0 );
      (*(*gup)[last_ ])[solindex_] += 0.5*alpha_*( (ampl_ - A*A > 0) ? 8.0*A - 2.0*A*A*A : 0 );
    }

    if(ptype_==4)
    {
      (*(*gup)[first_])[solindex_] += alpha_*( (A - ampl_ > 0) ? 0.5*(A - ampl_) : 0 );
      (*(*gup)[last_ ])[solindex_] -= alpha_*( (A - ampl_ > 0) ? 0.5*(A - ampl_) : 0 );
    }

    if(ptype_==5)
    {
      (*(*gup)[first_])[solindex_] = 0.5*alpha_;
      (*(*gup)[last_ ])[solindex_] = -0.5*alpha_;
    }
  }

  void gradient_2(::ROL::Vector<Real> &g,
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    // Unwrap g
    Teuchos::RCP<std::vector<Real> > gzp = Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<const ::ROL::StdVector<Real> >(g)).getVector());
    // Unwrap x
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Real> > zp = (Teuchos::dyn_cast< ::ROL::StdVector<Real> >(const_cast< ::ROL::Vector<Real> &>(z))).getVector();
    // COMPUTE GRADIENT WRT Z
    for (int i=0; i<nz_; i++)
    {
      (*gzp)[i] = 0.0;
    }
  }

  void hessVec_11(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    Teuchos::RCP< std::vector<Teuchos::RCP<Linear::Vector> > > hvup = Teuchos::rcp_const_cast< std::vector<Teuchos::RCP<Linear::Vector> > >((Teuchos::dyn_cast<Linear::ROL_XyceVector<Real> >(hv)).getVector());
    // Unwrap v
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > up = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(u))).getVector();
    Teuchos::RCP<const std::vector<Teuchos::RCP<Linear::Vector> > > vup = (Teuchos::dyn_cast< Linear::ROL_XyceVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v))).getVector();
    
    hv.zero();
    Real V3_0  = (*(*vup)[first_ ])[solindex_];
    Real V3_20 = (*(*vup)[last_  ])[solindex_];
    Real A = ((*(*up)[first_])[solindex_] - (*(*up)[last_])[solindex_] )/2.0;

    if(ptype_==1)
    {
      (*(*hvup)[first_])[solindex_] -= 0.5*alpha_*( (sinh(epsilon_*(ampl_-A*A))/std::pow(cosh(epsilon_*(ampl_-A*A)),3.0))*std::pow(epsilon_*A,2.0) + epsilon_/(std::pow(cosh(epsilon_*(ampl_-A*A)),2.0)*4.0) )*(V3_0 - V3_20);
      (*(*hvup)[last_ ])[solindex_] += 0.5*alpha_*( (sinh(epsilon_*(ampl_-A*A))/std::pow(cosh(epsilon_*(ampl_-A*A)),3.0))*std::pow(epsilon_*A,2.0) + epsilon_/(std::pow(cosh(epsilon_*(ampl_-A*A)),2.0)*4.0) )*(V3_0 - V3_20);
    }

    if(ptype_==2)
    {    
      (*(*hvup)[first_])[solindex_] += 0.5*alpha_*( (ampl_ - A*A > 0) ? -0.5 : 0 )*(V3_0 - V3_20);
      (*(*hvup)[last_ ])[solindex_] -= 0.5*alpha_*( (ampl_ - A*A > 0) ? -0.5 : 0 )*(V3_0 - V3_20);
    }

    if(ptype_==3)
    {
      (*(*hvup)[first_])[solindex_] += 0.5*alpha_*( (ampl_ - A*A > 0) ? 3.0*A*A - 4.0 : 0 )*(V3_0 - V3_20);
      (*(*hvup)[last_])[solindex_] -= 0.5*alpha_*( (ampl_ - A*A > 0) ? 3.0*A*A - 4.0 : 0 )*(V3_0 - V3_20);
    }

    if(ptype_==4)
    {
      (*(*hvup)[first_])[solindex_] += alpha_*( (A - ampl_ > 0) ? 0.25 : 0 )*(V3_0 - V3_20);
      (*(*hvup)[last_ ])[solindex_] += alpha_*( (A - ampl_ > 0) ? 0.25 : 0 )*(V3_20 - V3_0);
    }

  }

  void hessVec_12(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }
  
  void hessVec_21(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }
  
  void hessVec_22(::ROL::Vector<Real> &hv,
                  const ::ROL::Vector<Real> &v, 
                  const ::ROL::Vector<Real> &u,
                  const ::ROL::Vector<Real> &z,
                  Real &tol)
  {
    hv.zero();
  }

};

template<class Real>
class SumObjective : public virtual ::ROL::Objective<Real>
{

private:
  std::vector< Teuchos::RCP< ::ROL::Objective<Real> > > objVec_;
  std::vector<bool>                                     types_;

public:

  SumObjective(std::vector< Teuchos::RCP< ::ROL::Objective<Real> > > & objVec,
               std::vector<bool>                                     & types)
    : objVec_(objVec), types_(types)
  {}
  
  virtual void update( const ::ROL::Vector<Real> &x, bool flag = true, int iter = -1 )
  {
    for(int i=0;i<objVec_.size();i++)
    {
      if(types_[i]==0)
      {
        // ::ROL::CVaRVector<Real> &xc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(const_cast< ::ROL::Vector<Real> &>(x));
        // objVec_[i]->update(*(xc.getVector()),flag,iter); // CVaR vector
      }
      else if(types_[i]==1)
      {
        objVec_[i]->update(x,flag,iter);
      }
    }
  }

  virtual Real value(const ::ROL::Vector<Real> &x, 
                     Real &tol )
  {
    Real val = 0;
    for(int i=0;i<objVec_.size();i++)
    {
      if(types_[i]==0)
      {
        // ::ROL::CVaRVector<Real> &xc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(const_cast< ::ROL::Vector<Real> &>(x));
        // val += objVec_[i]->value(*xc.getVector(),tol); // CVaR vector
      }
      else if(types_[i]==1)
      {
        objVec_[i]->update(x,true,-1);
        val += objVec_[i]->value(x,tol);
      }
    }
    return val;
  }

  virtual void gradient(::ROL::Vector<Real> &g,
                        const ::ROL::Vector<Real> &x,
                        Real &tol)
  {
    g.zero();
    Teuchos::RCP< ::ROL::Vector<Real> > temp = g.clone();    
    for(int i=0;i<objVec_.size();i++)
    {
      if(types_[i]==0)
      {
        // ::ROL::CVaRVector<Real> &tempc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(*temp);
        // ::ROL::CVaRVector<Real> &xc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(const_cast< ::ROL::Vector<Real> &>(x));
        // objVec_[i]->gradient(const_cast< ::ROL::Vector<Real> &>(*tempc.getVector()),*xc.getVector(),tol); // CVaR vector
        // g.plus(tempc);
      }
      else if(types_[i]==1)
      {
        objVec_[i]->gradient(*temp,x,tol);
        g.plus(*temp);
      }
    }
  }
  
  virtual void hessVec(::ROL::Vector<Real> &hv,
                       const ::ROL::Vector<Real> &v,
                       const ::ROL::Vector<Real> &x,
                       Real &tol)
  {
    hv.zero();
    Teuchos::RCP< ::ROL::Vector<Real> > temp = hv.clone();    
    for(int i=0;i<objVec_.size();i++)
    {
      if(types_[i]==0)
      {
        // ::ROL::CVaRVector<Real> &tempc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(*temp);
        // ::ROL::CVaRVector<Real> &xc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(const_cast< ::ROL::Vector<Real> &>(x));
        // ::ROL::CVaRVector<Real> &vc = Teuchos::dyn_cast< ::ROL::CVaRVector<Real> >(const_cast< ::ROL::Vector<Real> &>(v));
        // objVec_[i]->hessVec(const_cast< ::ROL::Vector<Real> &>(*tempc.getVector()),*vc.getVector(),*xc.getVector(),tol); // CVaR vector
        // hv.plus(tempc);
      }
      else if(types_[i]==1)
      {
        objVec_[i]->update(x,true,-1);
        objVec_[i]->hessVec(*temp,v,x,tol);
        hv.plus(*temp);
      }
    }
  }
  
};

template<class Real>
class BoundConstraint_ROL_DC : public ::ROL::BoundConstraint<Real>
{

private:
  /// Vector of lower bounds
  std::vector<Real> x_lo_;
  /// Vector of upper bounds
  std::vector<Real> x_up_;
  /// Half of the minimum distance between upper and lower bounds
  Real min_diff_;
  /// Scaling for the epsilon margin
  Real scale_;
  int n;
  
public:
  BoundConstraint_ROL_DC( Real scale, std::vector<Real> lo_bounds, std::vector<Real> up_bounds )
  {
    //assume lo_bounds and up_bounds are same size            
    n = lo_bounds.size();
    for (int i=0;i<n;i++)
    {
      x_lo_.push_back(lo_bounds[i]);
    }
    for (int i=0;i<n;i++)
    {
      x_up_.push_back(up_bounds[i]);
    }

    scale_ = scale;
    min_diff_ = 1.e+200;
    for(int i=0;i<n;i++)
    {
      min_diff_ = std::min(x_up_[i]-x_lo_[i],min_diff_);
    }
    min_diff_ *= 0.5;
  }

  void project( ::ROL::Vector<Real> &x )
  {
    Teuchos::RCP<std::vector<Real> > ex =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(x)).getVector());
    for (int i=0;i<n;i++)
    {
      (*ex)[i] = std::max(x_lo_[i],std::min(x_up_[i],(*ex)[i]));
    }
  }
  bool isFeasible( const ::ROL::Vector<Real> &x )
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    bool bs = true;
    for (int i=0;i<n;i++)
    {
      bs = bs && (((*ex)[i] >= this->x_lo_[i])?true:false);
      bs = bs && (((*ex)[i] <= this->x_up_[i])?true:false);
    }
    return bs;
  }

  void pruneLowerActive(::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &x, Real eps)
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    for ( int i = 0; i < n; i++ )
    {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) )
      {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneUpperActive(::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &x, Real eps)
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    for ( int i = 0; i < n; i++ )
    {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn) )
      {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneActive(::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &x, Real eps)
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(this->scale_*eps,this->min_diff_);
    for ( int i = 0; i < n; i++ )
    {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ||
           ((*ex)[i] >= this->x_up_[i]-epsn) )
      {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneLowerActive(::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &g, const ::ROL::Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0))
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(this->scale_*xeps,this->min_diff_);
    for ( int i = 0; i < n; i++ )
    {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > geps) )
      {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneUpperActive(::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &g, const ::ROL::Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0))
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(this->scale_*xeps,this->min_diff_);
    for ( int i = 0; i < n; i++ )
    {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < -geps) )
      {
        (*ev)[i] = 0.0;
      }
    }
  }

  void pruneActive(::ROL::Vector<Real> &v, const ::ROL::Vector<Real> &g, const ::ROL::Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0))
  {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<::ROL::StdVector<Real> >(const_cast<::ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<::ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(this->scale_*xeps,this->min_diff_);
    for ( int i = 0; i < n; i++ )
    {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] >  geps) ||
           ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < -geps) )
      {
        (*ev)[i] = 0.0;
      }
    }
  }

  // void setVectorToUpperBound( ::ROL::Vector<Real> &u )
  // {
  //   Teuchos::RCP<std::vector<Real> > us = Teuchos::rcp( new std::vector<Real>(2,0.0) );
  //   us->assign(this->x_up_.begin(),this->x_up_.end());
  //   Teuchos::RCP<::ROL::Vector<Real> > up = Teuchos::rcp( new ::ROL::StdVector<Real>(us) );
  //   u.set(*up);
  // }

  // void setVectorToLowerBound( ::ROL::Vector<Real> &l )
  // {
  //   Teuchos::RCP<std::vector<Real> > ls = Teuchos::rcp( new std::vector<Real>(2,0.0) );
  //   ls->assign(this->x_lo_.begin(),this->x_lo_.end());
  //   Teuchos::RCP<::ROL::Vector<Real> > lp = Teuchos::rcp( new ::ROL::StdVector<Real>(ls) );
  //   l.set(*lp);
  // }

  const Teuchos::RCP<const ::ROL::Vector<Real>> getLowerBound(void) const
  {
    Teuchos::RCP<std::vector<Real>> xp = Teuchos::rcp_const_cast<std::vector<Real>>(Teuchos::rcp(&x_lo_));
    return ::ROL::makePtr<::ROL::StdVector<Real>>(xp);
  }

  const Teuchos::RCP<const ::ROL::Vector<Real>> getUpperBound(void) const
  {
    Teuchos::RCP<std::vector<Real>> xp = Teuchos::rcp_const_cast<std::vector<Real>>(Teuchos::rcp(&x_up_));
    return ::ROL::makePtr<::ROL::StdVector<Real>>(xp);
  }

};

} // namespace Analysis
} // namespace Xyce`

#endif // Xyce_ROL

#endif // Xyce_N_ANP_ROL_DC_Optimization

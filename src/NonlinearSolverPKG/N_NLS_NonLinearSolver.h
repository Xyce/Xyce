//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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

//-------------------------------------------------------------------------
//
// Purpose        : Specification file which declares an interface common to
//                  all supported nonlinear solver algorithms.  The Manager
//                  class uses this interface to call a concrete algorithm.
//
// Special Notes  : This is the "Strategy" class in the Strategy design
//                  pattern.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NonLinearSolver_h
#define Xyce_N_NLS_NonLinearSolver_h

// ---------- Standard Includes ----------

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_NLS_fwd.h>
#include <N_TIA_fwd.h>

#include <N_NLS_ReturnCodes.h>
#include <N_NLS_NonLinInfo.h>
#include <N_UTL_Stats.h>

// ---------- Using Declarations ----------
using Teuchos::RCP;

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : NonLinearSolver
// Purpose       : Nonlinear Solver Abstract Class
// Special Notes : Many of the virtual functions should not, in general,
//                 be redefined in derived classes. Check the Virtual Notes
//                 on each function.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
class NonLinearSolver
{
  friend class ConductanceExtractor;
  friend class Sensitivity;
  friend class TwoLevelNewton;

public:
  NonLinearSolver(const IO::CmdParse & cp);
  virtual ~NonLinearSolver();

protected:
  const Analysis::AnalysisManager &getAnalysisManager() const
  {
    return *analysisManager_;
  }

  Analysis::AnalysisManager &getAnalysisManager()
  {
    return *analysisManager_;
  }

public:
  bool getMatrixFreeFlag() const
  {
    return matrixFreeFlag_;
  }

  void setMatrixFreeFlag(bool matrixFreeFlag)
  {
    matrixFreeFlag_ = matrixFreeFlag;
  }

  virtual bool setOptions(const Util::OptionBlock& OB) = 0;
  virtual bool setTranOptions(const Util::OptionBlock& OB) = 0;
  virtual bool setHBOptions(const Util::OptionBlock& OB) = 0;
  virtual bool setNLPOptions(const Util::OptionBlock& OB) = 0;
  virtual bool setLocaOptions(const Util::OptionBlock& OB);
  virtual bool setTwoLevelLocaOptions(const Util::OptionBlock& OB);
  virtual bool setTwoLevelOptions    (const Util::OptionBlock& OB);
  virtual bool setTwoLevelTranOptions(const Util::OptionBlock& OB);
  virtual bool setLinsolOptions(const Util::OptionBlock& OB);
  virtual bool setICOptions(const Util::OptionBlock& OB);
  virtual bool setNodeSetOptions(const Util::OptionBlock& OB);

  virtual bool registerLinearSystem(Linear::System* ptr);
  virtual bool registerAnalysisManager(Analysis::AnalysisManager* tmp_anaIntPtr);
  virtual bool registerNonlinearEquationLoader(Loader::NonlinearEquationLoader* ptr);
  virtual bool registerTIADataStore(TimeIntg::DataStore * ptr);
  virtual bool registerParallelMgr(N_PDS_Manager * ptr);

  bool registerTwoLevelSolver (TwoLevelNewton * ptr);

  virtual bool registerSolverFactory(const Linear::SolverFactory *ptr);
  virtual bool registerPrecondFactory(const Linear::PrecondFactory *ptr);

  virtual bool registerOutputMgr (IO::OutputMgr * outPtr);
  virtual bool registerInitialConditionsManager(IO::InitialConditionsManager * ptr);

  virtual bool initializeAll();

  virtual int solve(NonLinearSolver * nlsTmpPtr = NULL) = 0;
  virtual inline int takeFirstSolveStep (NonLinearSolver * nlsTmpPtr = NULL);
  virtual inline int takeOneSolveStep ();

  virtual int getNumIterations() const = 0;

  virtual int getDebugLevel() const = 0;
  virtual bool getScreenOutputFlag () const = 0;
  virtual double getDebugMinTime() const = 0;
  virtual double getDebugMaxTime() const = 0;
  virtual int getDebugMinTimeStep() const = 0;
  virtual int getDebugMaxTimeStep() const = 0;
  virtual bool getMMFormat () const = 0;

  virtual bool isFirstContinuationParam() const = 0;
  virtual bool isFirstSolveComplete() const = 0;
  virtual int getContinuationStep() const = 0;
  virtual int getParameterNumber() const = 0;

  virtual bool getLocaFlag ();

  virtual inline int getNumResidualLoads();
  virtual inline int getNumJacobianLoads();
  virtual inline int getNumLinearSolves();
  virtual inline int getNumFailedLinearSolves();
  virtual inline int getNumJacobianFactorizations();
  virtual inline unsigned int getTotalNumLinearIters();
  virtual inline double getTotalLinearSolveTime();
  virtual inline double getTotalResidualLoadTime();
  virtual inline double getTotalJacobianLoadTime();

  virtual TwoLevelNewtonMode getCouplingMode ();

  virtual void setAnalysisMode(AnalysisMode mode) = 0;
  virtual void resetAll (AnalysisMode mode);
  virtual void setReturnCodes (const ReturnCodes & retCodesTmp);
  virtual bool enableSensitivity () {return true;}
  virtual double getMaxNormF() const = 0;
  virtual int getMaxNormFindex () const = 0;

  Teuchos::RCP<Linear::Solver> getLinearSolver();
  void setLinearSolver( const Teuchos::RCP<Linear::Solver>& lasSolver );

  // use for debugging:
  void debugOutput1(Linear::Matrix & jacobian, Linear::Vector & rhs);
  void debugOutput3(Linear::Vector & dxVector, Linear::Vector & xVector);

  void debugOutputDAE();

  void  setDebugFlags(int output_step_number, double        time);
  
  virtual bool applyJacobian(const Linear::Vector& input, Linear::Vector& result);

public:
  // for now.
#ifdef Xyce_ROL
  Linear::Vector* rhsVectorPtr_;
  Linear::Vector* NewtonVectorPtr_;
  Teuchos::RCP<Linear::Solver> lasSolverRCPtr_;
#endif
  
protected:

  virtual void resetCountersAndTimers_();
  virtual bool setX0_();
  virtual bool rhs_();
  virtual bool jacobian_();
  virtual bool newton_();
  virtual bool gradient_();

protected:
  const IO::CmdParse & commandLine_;
  std::string netlistFilename_;
#ifndef Xyce_ROL
  Linear::Vector* rhsVectorPtr_;
#endif

  Linear::Matrix* jacobianMatrixPtr_;
  Linear::Vector* gradVectorPtr_;
#ifndef Xyce_ROL
  Linear::Vector* NewtonVectorPtr_;
#endif
  Linear::Vector* solWtVectorPtr_;
  Linear::System* lasSysPtr_;
#ifndef Xyce_ROL
  Teuchos::RCP<Linear::Solver> lasSolverRCPtr_;
#endif
  Linear::Problem * lasProblemPtr_;
  const Linear::SolverFactory *        lasSolverFactoryPtr_;
  const Linear::PrecondFactory *        lasPrecFactoryPtr_;
  Util::OptionBlock* linsolOptionBlockPtr_;
  Loader::NonlinearEquationLoader* nonlinearEquationLoader_;
  TwoLevelNewton * tlnPtr_;
  ParamMgr * nonlinearParameterManager_;
  IO::OutputMgr * outMgrPtr_;
  IO::InitialConditionsManager * initialConditionsManager_;
  N_PDS_Manager * pdsMgrPtr_;
  TimeIntg::DataStore * dsPtr_;

  int numJacobianLoads_;
  int numJacobianFactorizations_;
  int numLinearSolves_;
  int numFailedLinearSolves_;
  int numResidualLoads_;
  unsigned int totalNumLinearIters_;
  double totalLinearSolveTime_;
  double totalResidualLoadTime_;
  double totalJacobianLoadTime_;
  ReturnCodes retCodes_;

  bool debugTimeFlag_;

  int contStep_;

private:
  Analysis::AnalysisManager *   analysisManager_;
  int                           outputStepNumber_;  // this is either the time step number or the dc sweep step number,
                                                    // depending on the mode.  It is only used in setting up output file
                                                    // names.
  bool matrixFreeFlag_;
};

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumResidualLoads
// Return Type   : Integer (Get the total number of residual loads)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumResidualLoads()
{
  return numResidualLoads_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumJacobianLoads
// Return Type   : Integer (Get the total number of Jacobian loads)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumJacobianLoads()
{
  return numJacobianLoads_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumLinearSolves
// Return Type   : Integer (total number of successful linear solves)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumLinearSolves()
{
  return numLinearSolves_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumFailedLinearSolves
// Return Type   : Integer (total number of failed linear solves)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumFailedLinearSolves()
{
  return numFailedLinearSolves_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getNumJacobianFactorizations
// Return Type   : Integer (total number of Jacobian factorizations)
//---------------------------------------------------------------------------
inline int NonLinearSolver::getNumJacobianFactorizations()
{
  return numJacobianFactorizations_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalNumLinearIters
// Return Type   : unsigned int (total number of iterative linear solver
//                 iterations)
//---------------------------------------------------------------------------
inline unsigned int NonLinearSolver::getTotalNumLinearIters()
{
  return totalNumLinearIters_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalLinearSolveTime
// Return Type   : double (total linear solve time in seconds)
//---------------------------------------------------------------------------
inline double NonLinearSolver::getTotalLinearSolveTime()
{
  return totalLinearSolveTime_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalResidualLoadTime
// Return Type   : double (total residual load (claculation) time in seconds)
//---------------------------------------------------------------------------
inline double NonLinearSolver::getTotalResidualLoadTime()
{
  return totalResidualLoadTime_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::getTotalJacobianLoadTime
// Return Type   : double (total Jacobian load (claculation) time in seconds)
//---------------------------------------------------------------------------
inline double NonLinearSolver::getTotalJacobianLoadTime()
{
  return totalJacobianLoadTime_;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::takeFirstSolveStep
// Return Type   : int
//---------------------------------------------------------------------------
inline int NonLinearSolver::takeFirstSolveStep(
  NonLinearSolver *     nlsTmpPtr)
{
  return -1;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::takeOneSolveStep
// Return Type   : int
//---------------------------------------------------------------------------
inline int NonLinearSolver::takeOneSolveStep ()
{
  return -1;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::resetAll
// Return Type   : void
//---------------------------------------------------------------------------
inline void NonLinearSolver::resetAll (AnalysisMode mode)
{
  setAnalysisMode(mode);
}
//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setReturnCodes
// Purpose       : Allows the user to set return codes.
//
// Special Notes : This was put in place mainly to help with continuation
//                 solves.  If running continuation, I don't want the
//                 solver to consider "nearConvergence" to be a success.
//                 There are other circumstances that I may wish to turn
//                 this off as well.  For example, if running with the
//                 predictor/corrector step error control turned off (this
//                 is an option), I also would not want "nearConverged" to
//                 be good enough.
//
//                 The rational for having a "nearConverged" option  is the
//                 idea that we can just let the time integrator handle
//                 things and let it use the predictor/corrector stuff to
//                 determine whether or not to reject the step.  If that
//                 stuff is turned off, you otherwise unavailable, you
//                 don't want to do this.
//
//                 This way, the calling code can, if it wants, decide for
//                 itself what senarios will be considered a solver success.
//
// Scope         : public
// Creator       : Eric R. Keiter, 9233
// Creation Date : 10/30/04
//-----------------------------------------------------------------------------
inline void
NonLinearSolver::setReturnCodes(
  const ReturnCodes &   ret_codes)
{
  retCodes_ = ret_codes;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::getLocaFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/09/05
//-----------------------------------------------------------------------------
inline bool NonLinearSolver::getLocaFlag ()
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::getLinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi K. Thornquist, SNL
// Creation Date : 10/4/16
//-----------------------------------------------------------------------------
inline Teuchos::RCP<Linear::Solver> NonLinearSolver::getLinearSolver()
{
  return lasSolverRCPtr_;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setLinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi K. Thornquist, SNL
// Creation Date : 10/4/16
//-----------------------------------------------------------------------------
inline void NonLinearSolver::setLinearSolver( const Teuchos::RCP<Linear::Solver>& lasSolver )
{
  lasSolverRCPtr_ = lasSolver;
}

} // namespace Nonlinear
} // namespace Xyce

#endif


//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Body for the two level Newton class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/20/02
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------
#include <cstdio>
#include <string>
#include <list>

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Solver.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_DampedNewton.h>
#include <N_NLS_Manager.h>
#include <N_NLS_NOX_Interface.h>
#include <N_NLS_ReturnCodes.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>

#include <N_PDS_Comm.h>

// ----------   Static Declarations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::TwoLevelNewton
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
TwoLevelNewton::TwoLevelNewton(
    bool                  noxFlag,
    bool                  noxFlagInner,
    const IO::CmdParse &  cp)
  : NonLinearSolver(cp),
  maxOuterSteps_(20),
  maxContSteps_(10),
  contStep_(0),
  increaseContScalar_(1.5),
  decreaseContScalar_(0.2), // was 0.125
  twoLevelAlgorithm_(3),
  twoLevelAlgorithmTran_(0),
  outerLoopActiveFlag_(true),
  externalAnalysisMode(DC_OP),
  setupOuterLoopParamsFlag_(false),
  setupTranParamsFlag_(false),
  noxFlag_(noxFlag),
  noxFlagInner_(noxFlagInner),//inner loop possibly different than outer loop.
  numInterfaceNodesSetup_(false),
  twoLevelCouplingMode_(FULL_PROBLEM),
  savedRHSPtr_(0),
  savedNextSolPtr_(0),
  jdxpVectorPtr_(0),
  numSubProblems_(0),
  continuationType_(1),
  innerLoopFailFatal_(true),
  totalSolveFailFatal_(false),
  doFullNewtonFinalEnforcement_(true),
  nlsPassingPtr_(0),
  continuationCalledBefore_(false),
  voltLimTol_(1.0e-6),
  reuseFactors_(true)
{
  // allocate the "outer loop" and "inner loop" solvers.
  if (noxFlag_)
  {
    nlsOuterPtr_ = new N_NLS_NOX::Interface(commandLine_);
  }
  else
  {
    nlsOuterPtr_ = new DampedNewton(commandLine_);
  }

  if (noxFlagInner_)
  {
    nlsInnerPtr_ = new N_NLS_NOX::Interface(commandLine_);
  }
  else
  {
    nlsInnerPtr_ = new DampedNewton(commandLine_);
  }

  nlsOuterPtr_->registerTwoLevelSolver(this);
  nlsInnerPtr_->registerTwoLevelSolver(this);

}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::~TwoLevelNewton
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
TwoLevelNewton::~TwoLevelNewton()
{
  delete nlsOuterPtr_;
  delete nlsInnerPtr_;

  if (savedRHSPtr_!=0) delete savedRHSPtr_;
  if (savedNextSolPtr_!=0) delete savedNextSolPtr_;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerAnalysisManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerAnalysisManager(Analysis::AnalysisManager * analysis_manager)
{
  bool bsuccess = true;
  bool tmpBool = true;

  analysisManager_ = analysis_manager;
  tmpBool = nlsOuterPtr_->registerAnalysisManager(analysisManager_);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerAnalysisManager(analysisManager_);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerTIADataStore
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerTIADataStore(TimeIntg::DataStore * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  dsPtr_ = ptr;

  tmpBool = nlsOuterPtr_->registerTIADataStore(ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerTIADataStore(ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerParallelMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerParallelMgr (Parallel::Manager * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  pdsMgrPtr_ = ptr;

  tmpBool = nlsOuterPtr_->registerParallelMgr(ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerParallelMgr(ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerLinearSystem (Linear::System * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  lasSysPtr_ = ptr;
  tmpBool = nlsOuterPtr_->registerLinearSystem (ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerLinearSystem (ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerNonlinearEquationLoader(Loader::NonlinearEquationLoader * loaderPtr_tmp)
{
  bool bsuccess = true;
  bool tmpBool = true;

  nonlinearEquationLoader_ = loaderPtr_tmp;
  tmpBool = nlsOuterPtr_->registerNonlinearEquationLoader(loaderPtr_tmp);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerNonlinearEquationLoader(loaderPtr_tmp);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/23/03
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerOutputMgr (IO::OutputMgr * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  outMgrPtr_ = ptr;
  tmpBool = nlsOuterPtr_->registerOutputMgr (ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerOutputMgr (ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerInitialConditionsManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/25/2019
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerInitialConditionsManager(IO::InitialConditionsManager * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  initialConditionsManager_ = ptr;
  tmpBool = nlsOuterPtr_->registerInitialConditionsManager (ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerInitialConditionsManager (ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerSolverFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/25/2019
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerSolverFactory(const Linear::SolverFactory * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  lasSolverFactoryPtr_ = ptr;
  tmpBool = nlsOuterPtr_->registerSolverFactory (ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerSolverFactory (ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::registerPrecondFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/25/2019
//-----------------------------------------------------------------------------
bool TwoLevelNewton::registerPrecondFactory(const Linear::PrecondFactory * ptr)
{
  bool bsuccess = true;
  bool tmpBool = true;

  lasPrecFactoryPtr_ = ptr;
  tmpBool = nlsOuterPtr_->registerPrecondFactory (ptr);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->registerPrecondFactory (ptr);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getCouplingMode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/05/02
//-----------------------------------------------------------------------------
TwoLevelNewtonMode TwoLevelNewton::getCouplingMode ()
{
  return twoLevelCouplingMode_;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getNumResidualLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::getNumResidualLoads()
{
  int numResLoads = 0;
  numResLoads += nlsOuterPtr_->getNumResidualLoads ();
  numResLoads += numResidualLoads_;
  return numResLoads;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getNumJacobianLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::getNumJacobianLoads()
{
  int numJacLoads = 0;
  numJacLoads += nlsOuterPtr_->getNumJacobianLoads ();
  numJacLoads += numJacobianLoads_;
  return numJacLoads;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getNumLinearSolves
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::getNumLinearSolves()
{
  int numLinSolves = 0;
  numLinSolves += nlsOuterPtr_->getNumLinearSolves ();
  numLinSolves += numLinearSolves_;
  return numLinSolves;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getNumFailedLinearSolves
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::getNumFailedLinearSolves()
{
  int numFLinSolves = 0;
  numFLinSolves += nlsOuterPtr_->getNumFailedLinearSolves ();
  numFLinSolves += numFailedLinearSolves_;
  return numFLinSolves;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getNumJacobianFactorizations
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::getNumJacobianFactorizations()
{
  int numJacFact = 0;
  numJacFact += nlsOuterPtr_->getNumJacobianFactorizations ();
  numJacFact += numJacobianFactorizations_;
  return numJacFact;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getTotalNumLinearIters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
unsigned int TwoLevelNewton::getTotalNumLinearIters()
{
  int numLinIters  = 0;
  numLinIters  += nlsOuterPtr_->getTotalNumLinearIters ();
  numLinIters  += totalNumLinearIters_;
  return numLinIters ;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getTotalLinearSolveTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
double TwoLevelNewton::getTotalLinearSolveTime()
{
  double totalLinSolveTime = 0.0;
  totalLinSolveTime += nlsOuterPtr_->getTotalLinearSolveTime();
  totalLinSolveTime += totalLinearSolveTime_;
  return totalLinSolveTime;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getTotalResidualLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
double TwoLevelNewton::getTotalResidualLoadTime()
{
  double totalResLoadTime = 0.0;
  totalResLoadTime += nlsOuterPtr_-> getTotalResidualLoadTime();
  totalResLoadTime += totalResidualLoadTime_;
  return totalResLoadTime;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getTotalJacobianLoadTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
double TwoLevelNewton::getTotalJacobianLoadTime()
{
  double totalJacLoadTime = 0.0;
  totalJacLoadTime += nlsOuterPtr_->getTotalJacobianLoadTime();
  totalJacLoadTime += totalJacobianLoadTime_;
  return totalJacLoadTime;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getNumIterations
// Purpose       : Returns the current number of nonlinear iterations, for
//                 the solver currently being used (for two-level more than
//                 one solver is usually being invoked).
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Compuational Sciences
// Creation Date : 10/24/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::getNumIterations() const
{
  int numIters = 0;
  if (outerLoopActiveFlag_)
  {
    numIters = nlsOuterPtr_->getNumIterations ();
  }
  else
  {
    numIters = nlsInnerPtr_->getNumIterations ();
  }

  return numIters;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getMaxNormF
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Memes Modeling
// Creation Date : 9/27/2009
//-----------------------------------------------------------------------------
double TwoLevelNewton::getMaxNormF() const
{
  double result = nlsInnerPtr_->getMaxNormF() + nlsOuterPtr_->getMaxNormF();
  return result;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::getMaxNormFindex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Memes Modeling
// Creation Date : 9/27/2009
//-----------------------------------------------------------------------------
int TwoLevelNewton::getMaxNormFindex() const
{
  int indexInner = nlsInnerPtr_->getMaxNormFindex();
  int indexOuter = nlsOuterPtr_->getMaxNormFindex();

  return indexInner; // usually the inner will be the one you want,
  // but there should be more of a detailed comparison
  // here just in case.  FIX THIS
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::initializeAll
//
// Purpose       : This serves the same purpose as the initializeAll
//                 function in the other solvers.
//
// Special Notes : Each of the solvers (inner, outer, and this one) will
//                 ultimately call the base class initialize all function,
//                 meaning that there are potentially 3 of each linear
//                 solver (3 iterative, 3 direct).  This sort of makes
//                 sense for the inner and outer, as they may want
//                 different linear solvers.  But, for the wrapper class,
//                 TwoLevelNewton, I'm only doing this so that it is
//                 possible to pass a solver pointer to Sensitivity,
//                 if needed.
//
//                 In the future, it may make sense to clean this up a
//                 litte.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::initializeAll ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  tmpBool = nlsInnerPtr_->initializeAll ();
  bsuccess = bsuccess && tmpBool;

  tmpBool= nlsOuterPtr_->initializeAll ();
  bsuccess = bsuccess && tmpBool;

  tmpBool = NonLinearSolver::initializeAll();
  bsuccess = bsuccess && tmpBool;

  savedRHSPtr_ = lasSysPtr_->builder().createVector ();
  savedNextSolPtr_ = lasSysPtr_->builder().createVector ();
  jdxpVectorPtr_ = lasSysPtr_->getJDXPVector ();

  // set up the return codes so that the "inner" solver is subject to
  // greater restrictions than the outter solver.
  ReturnCodes retCode;
  retCode.nearConvergence = -3;
  nlsInnerPtr_->setReturnCodes (retCode);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setOptions
//
// Purpose       : This function processes the options set in the
//                 ".options NONLIN" line of the netlist.
//
// Special Notes : Mostly, these options are just passed on through to the
//                 "outer" nonlinear solver.  However, since the outer
//                 solver is incremented one Newton step at a time, and
//                 the actual outer control loop sits here in the two level
//                 class, a few of the parameters are needed.  In
//                 particular, the maximum number of Newton steps.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setOptions(const Util::OptionBlock & OB)
{
  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for ( ; it_tpL != end_tpL; ++it_tpL)
  {
    if (it_tpL->uTag() == "MAXSTEP")
    {
      maxOuterSteps_ = static_cast<int>(it_tpL->getImmutableValue<int>());
    }
  }

  return nlsOuterPtr_->setOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setTranOptions(const Util::OptionBlock & OB)
{
  return nlsOuterPtr_->setTranOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 04/07/2015
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setHBOptions(const Util::OptionBlock & OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setNLPOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/07/2015
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setNLPOptions(const Util::OptionBlock & OB)
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/10/03
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setLocaOptions (const Util::OptionBlock & OB)
{
  outerLocaOptions_ = OB;
  return nlsOuterPtr_->setLocaOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setTwoLevelLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/10/03
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setTwoLevelLocaOptions (const Util::OptionBlock & OB)
{
  innerLocaOptions_ = OB;
  return nlsInnerPtr_->setLocaOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setTwoLevelOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setTwoLevelOptions(const Util::OptionBlock & OB)
{
  Util::OptionBlock OBtmp;

  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for ( ; it_tpL != end_tpL; ++it_tpL)
  {
    if (it_tpL->uTag() == "ALGORITHM")
    {
      twoLevelAlgorithm_ = static_cast<int>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "NOX")
    {
      noxFlagInner_ = it_tpL->getImmutableValue<int>();
    }
    else if (it_tpL->uTag() == "MAXCONTSTEPS")
    {
      maxContSteps_ = static_cast<int>(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "CONTINUATIONFLAG")
    {
      int tmp = static_cast<int>(it_tpL->getImmutableValue<int>());
      continuationType_ = tmp;
    }
    else if (it_tpL->uTag() == "INNERFAIL")
    {
      int tmp = static_cast<int>(it_tpL->getImmutableValue<int>());
      innerLoopFailFatal_ = (tmp!=0);
    }
    else if (it_tpL->uTag() == "EXITWITHFAILURE")
    {
      int tmp = static_cast<int>(it_tpL->getImmutableValue<int>());
      totalSolveFailFatal_ = (tmp!=0);
    }
    else if (it_tpL->uTag() == "FULLNEWTONENFORCE")
    {
      int tmp = static_cast<int>(it_tpL->getImmutableValue<int>());
      doFullNewtonFinalEnforcement_ = (tmp!=0);
    }
    else if (it_tpL->uTag() == "CONPARAM")
    {
      paramNameList.push_back(it_tpL->stringValue());
    }
    else if (it_tpL->uTag() == "VOLTLIMTOL")
    {
      voltLimTol_ = static_cast<double>(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "REUSEFACTORS")
    {
      reuseFactors_ = static_cast<bool>(it_tpL->getImmutableValue<bool>());
    }
    else // anything that is not a "special" two-level param, push
      // back on the the tmp param list.
    {
      OBtmp.addParam(*it_tpL);
    }
  }

  innerSolverOptions_ = OBtmp;  // keep a copy of these.

  nlsInnerPtr_->setOptions(OBtmp);

  if (twoLevelAlgorithm_ < 0 || twoLevelAlgorithm_ > 5)
  {
    Report::UserWarning0() << "Now the only algorithms 0 to 5 are two-level algorithm.  Resetting to 0";

    twoLevelAlgorithm_ = 0;
  }

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << "\n" << std::endl
      << Xyce::section_divider << std::endl
      << "\n***** 2-level Inner Loop Nonlinear solver options:\n" << std::endl
      << "\talgorithm:\t\t\t" <<  twoLevelAlgorithm_ << std::endl
      << "\toutersteps:\t\t\t" <<  maxOuterSteps_ << std::endl
      << "\tmaxContSteps:\t\t\t" <<  maxContSteps_ << std::endl;

    Xyce::dout() << "\n***** Done printing Inner Loop Params:\n" << std::endl
      << Xyce::section_divider << std::endl
      << "\n" << std::endl;
  }

  // Now that the loop is done allocate the param val array:
  paramFinalVal.resize(paramNameList.size(),0.0);
  paramCurrentVal.resize(paramNameList.size(),0.0);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setTwoLevelTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::setTwoLevelTranOptions(const Util::OptionBlock & OB)
{
  setupTranParamsFlag_ = true;
  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << "In TwoLevelNewton::setTwoLevelTranOptions" << std::endl;
  }

  Util::OptionBlock OBtmp;

  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for ( ; it_tpL != end_tpL; ++it_tpL)
  {
    if (it_tpL->uTag() == "ALGORITHM")
    {
      twoLevelAlgorithmTran_ = static_cast<int>(it_tpL->getImmutableValue<int>());
    }
    else if ( it_tpL->uTag() == "MAXCONTSTEPS" )
    {
      maxContStepsTran_ = static_cast<int>(it_tpL->getImmutableValue<int>());
    }
    else // anything that is not a "special" two-level param, push
      // back on the the tmp param list.
    {
      OBtmp.addParam(*it_tpL);
    }
  }

  nlsInnerPtr_->setTranOptions(OBtmp);

  if (twoLevelAlgorithmTran_ < 0 || twoLevelAlgorithmTran_ > 3)
  {
    Report::UserWarning0() << "Now the only algorithms 0 to 3 are two-level algorithm.  Resetting to 0";

    twoLevelAlgorithmTran_ = 0;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setAnalysisMode
// Purpose       : This function is slightly different than the
//                 function of the same name in NLParams.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/21/02
//-----------------------------------------------------------------------------
void TwoLevelNewton::setAnalysisMode (AnalysisMode mode)
{
  // this should be 0,1 or 2  (DC_OP, DC_SWEEP, or TRANSIENT)
  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl
      << "Setting the externalAnalysisMode = " << mode << std::endl;
  }

  externalAnalysisMode = mode;

  nlsOuterPtr_->setAnalysisMode(mode);
  nlsInnerPtr_->setAnalysisMode(mode);
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::setLinsolOptions
// Purpose       : see header file
// Special Notes : see header file
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/9/00
//-----------------------------------------------------------------------------

bool TwoLevelNewton::setLinsolOptions(const Util::OptionBlock & OB)
{
  bool bsuccess = true;
  bool tmpBool = true;
  tmpBool = nlsOuterPtr_->setLinsolOptions(OB);
  bsuccess = bsuccess && tmpBool;

  tmpBool = nlsInnerPtr_->setLinsolOptions(OB);
  bsuccess = bsuccess && tmpBool;

  tmpBool = NonLinearSolver::setLinsolOptions(OB);
  bsuccess = bsuccess && tmpBool;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::printStepInfo_
// Purpose       : Print out 2-level Newton step information.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/21/02
//-----------------------------------------------------------------------------
void TwoLevelNewton::printStepInfo_
  (int step, int success, TwoLevelNewtonMode solveType)
  {
    char tmpChar[128];
    if (solveType==FULL_PROBLEM)
    {
      Xyce::lout() << "---------- 2LNiter: " << step << "\t" << success << "\tFULL PROBLEM --------------------------------" << std::endl;
    }
    else if (solveType==INNER_PROBLEM)
    {
      Xyce::lout() << "---------- 2LNiter: " << step << "\t" << success << "\tINNER PROBLEM ----------------------------" << std::endl;
    }
    else
    {
      Xyce::lout() << "---------- 2LNiter: " << step << "\t" << success << "\tOUTER PROBLEM ----------------------------" << std::endl;
    }
  }

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::zeroInnerLoopStatistics_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/21/02
//-----------------------------------------------------------------------------
void TwoLevelNewton::zeroInnerLoopStatistics_ ()
{
  numResidualLoads_ = 0;
  numJacobianLoads_ = 0;
  numLinearSolves_ = 0;
  numFailedLinearSolves_ = 0;
  numJacobianFactorizations_ = 0;
  totalNumLinearIters_ = 0;
  totalLinearSolveTime_ = 0.0;
  totalResidualLoadTime_ = 0.0;
  totalJacobianLoadTime_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::calcInnerLoopStatistics_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/21/02
//-----------------------------------------------------------------------------
void TwoLevelNewton::calcInnerLoopStatistics_ ()
{
  numResidualLoads_          += nlsInnerPtr_->getNumResidualLoads ();
  numJacobianLoads_          += nlsInnerPtr_->getNumJacobianLoads ();
  numLinearSolves_           += nlsInnerPtr_->getNumLinearSolves ();
  numFailedLinearSolves_     += nlsInnerPtr_->getNumFailedLinearSolves ();
  numJacobianFactorizations_ += nlsInnerPtr_->getNumJacobianFactorizations();
  totalNumLinearIters_       += nlsInnerPtr_->getTotalNumLinearIters ();
  totalLinearSolveTime_      += nlsInnerPtr_->getTotalLinearSolveTime ();
  totalResidualLoadTime_     += nlsInnerPtr_->getTotalResidualLoadTime ();
  totalJacobianLoadTime_     += nlsInnerPtr_->getTotalJacobianLoadTime ();
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::calcOuterLoopStatistics_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/21/02
//-----------------------------------------------------------------------------
void TwoLevelNewton::calcOuterLoopStatistics_ ()
{
  numResidualLoads_          += nlsOuterPtr_->getNumResidualLoads ();
  numJacobianLoads_          += nlsOuterPtr_->getNumJacobianLoads ();
  numLinearSolves_           += nlsOuterPtr_->getNumLinearSolves ();
  numFailedLinearSolves_     += nlsOuterPtr_->getNumFailedLinearSolves ();
  numJacobianFactorizations_ += nlsOuterPtr_->getNumJacobianFactorizations();
  totalNumLinearIters_       += nlsOuterPtr_->getTotalNumLinearIters ();
  totalLinearSolveTime_      += nlsOuterPtr_->getTotalLinearSolveTime ();
  totalResidualLoadTime_     += nlsOuterPtr_->getTotalResidualLoadTime ();
  totalJacobianLoadTime_     += nlsOuterPtr_->getTotalJacobianLoadTime ();
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::algorithm0_
// Purpose       : This algorithm is the full newton algorithm.
//
// Special Notes : The main thing that is different about running this
//                 function, rather than just not using the 2-level class
//                 at all, is that here the conductance and capacitance
//                 calculations are performed at the end of the solve.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/07/03
//-----------------------------------------------------------------------------
int TwoLevelNewton::algorithm0_(bool nl_poisson_dcop)
{
  int status = -1;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl << "Running algorithm 0:" << std::endl;
  }

  // set up the local algorithm variable:
  int algorithm = twoLevelAlgorithm_;
  if ( externalAnalysisMode ==2) algorithm = twoLevelAlgorithmTran_;

  // This step is neccessary in case we've just switched algorithms.
  if ( algorithm == 0) twoLevelCouplingMode_=FULL_PROBLEM;

  status = nlsOuterPtr_->solve ();

  // Now do conductance/capacitance extractions, if it makes sense.
  if (!nl_poisson_dcop)
  {
    calcCouplingTerms_ ();
  }

  // Do this so that I can test the ConductanceExtractor class.
  twoLevelCouplingMode_ = FULL_PROBLEM;
  nonlinearEquationLoader_->loadJacobian ();

  return status;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::algorithm1_
// Purpose       : This is the first two-level algorithm implemented.
//                 The idea is to have nested Newton solves.
//
//                 Inner loop  = PDE-device only
//                 Outer loop  = Full Newton.
//
// Special Notes : Doesn't work with NOX yet...
//
//                 This will exit immediately if the inner loop fails.
//                 The reason for this is that the point of 2-level
//                 Newton is to solve problems that are difficult
//                 as a result of two very different problem types
//                 being coupled together.
//
//                 Ideally, for the two-level Newton, you have a
//                 circuit problem that can be easily solved on its
//                 own, and a PDE-device problem that can easily be
//                 solved on its own, and the only problem is that
//                 they can't easily be solved together.  Two level
//                 Newton separates them, and applies different
//                 (hopefully ideal) solver settings to each.
//
//                 If the PDE problem
//                 (the inner loop) can't be solved with its
//                 current settings, even in stand-alone mode, 2-level
//                 Newton won't help.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::algorithm1_()
{
  int status = -1;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl << "Running algorithm 1:" << std::endl;
  }

  bool firstOuterStepTaken = false;
  nlsPassingPtr_ = 0;

  if (status < 0)
  {
    int twoLevelStep;
    for (twoLevelStep=0;twoLevelStep<maxOuterSteps_; ++twoLevelStep)
    {
      // Device Only:
      twoLevelCouplingMode_       = INNER_PROBLEM;
      outerLoopActiveFlag_= false;
      int statInner = nlsInnerPtr_->solve (nlsPassingPtr_);
      nlsPassingPtr_ = 0;
      calcInnerLoopStatistics_ ();

      // Output the nonlinear solver step information:
      if (VERBOSE_NONLINEAR)
      {
        printStepInfo_(twoLevelStep+1,statInner, twoLevelCouplingMode_);
      }

      // If we can't even get the inner loop to converge,
      // just give up.
      if (statInner < 0)
      {
        break;
      }

      // Full Problem:
      twoLevelCouplingMode_ = FULL_PROBLEM;
      outerLoopActiveFlag_= true;
      if (firstOuterStepTaken)
      {
        status = nlsOuterPtr_->takeOneSolveStep ();
      }
      else
      {
        firstOuterStepTaken = true;
        status = nlsOuterPtr_->takeFirstSolveStep (nlsInnerPtr_);
      }

      nlsPassingPtr_ = nlsOuterPtr_;

      // Output the nonlinear solver step information:
      if (VERBOSE_NONLINEAR)
      {
        printStepInfo_(twoLevelStep+1,status, twoLevelCouplingMode_);
      }

      // exit?
      // Use the total "full Newton" error to evaluate this step.
      if (status >= 0)
      {
        break;
      }

    } // end of for loop
  } // end of status if statement

  if (VERBOSE_NONLINEAR && status >=0)
  {
    Xyce::dout() << "TWO LEVEL Newton succeeded!" << std::endl;
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::algorithm2_
// Purpose       : Similar to algorithm1, but this one allows for the
//                 inner PDE Newton solve to be gradually stepped up to
//                 the circuit imposed boundary conditions.  Essentially,
//                 the inner loop Newton solve is a continuation method.
//
//                 This could almost be called a "3-level" Newton, as there
//                 are now 3 nested loops.
//
// Special Notes : not finished yet...
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::algorithm2_()
{
  int status = -1;
  int statInner = 0;
  bool statusFull = false;
  bool firstOuterStepTaken = false;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl << "Running algorithm 2:" << std::endl;
  }

  // Start out with a single step of the full problem:
  // Full Problem:
  twoLevelCouplingMode_       = FULL_PROBLEM;
  firstOuterStepTaken = true;
  outerLoopActiveFlag_= true;
  status = nlsOuterPtr_->takeFirstSolveStep (nlsInnerPtr_);

  // Output the nonlinear solver step information:
  if (VERBOSE_NONLINEAR)
  {
    printStepInfo_(0,status,twoLevelCouplingMode_);
  }

  nlsPassingPtr_ = 0;
  if (status <= 0)
  {
    int twoLevelStep;
    for (twoLevelStep=0;twoLevelStep<maxOuterSteps_; ++twoLevelStep)
    {
      // Device Only:
      twoLevelCouplingMode_       = INNER_PROBLEM;
      outerLoopActiveFlag_= false;

      if (continuationType_ == 1)
      {
        statInner = continuationLoop_ ();
      }
      else if (continuationType_ == 2)
      {
        statInner = locaLoop_ ();
      }
      else
      {
        statInner = nlsInnerPtr_->solve (nlsPassingPtr_);
        nlsPassingPtr_ = 0;
        calcInnerLoopStatistics_ ();
      }

      // Output the nonlinear solver step information:
      if (VERBOSE_NONLINEAR)
      {
        printStepInfo_(twoLevelStep+1,statInner,twoLevelCouplingMode_);
      }

      if (innerLoopFailFatal_)
      {
        // If we can't even get the inner loop to converge,
        // just give up.
        if (statInner <= 0)
        {
          break;
        }
      }

      // Full Problem:
      twoLevelCouplingMode_       = FULL_PROBLEM;
      outerLoopActiveFlag_= true;

      if (firstOuterStepTaken)
      {
        status = nlsOuterPtr_->takeOneSolveStep ();
      }
      else
      {
        firstOuterStepTaken = true;
        status = nlsOuterPtr_->takeFirstSolveStep (nlsInnerPtr_);
      }
      nlsPassingPtr_ = nlsOuterPtr_;

      // Output the nonlinear solver step information:
      if (VERBOSE_NONLINEAR)
      {
        printStepInfo_(twoLevelStep+1,status,twoLevelCouplingMode_);
      }

      // exit?
      // Use the total "full Newton" error to evaluate this step.
      // If the ckt considers itself converged at this point, do a
      // final reality check by taking a full Newton step.
      statusFull = false;

      if (status > 0 && statInner>0) statusFull = true;

      if (statusFull) break;

    } // end of for loop
  } // end of status if statement

  if (VERBOSE_NONLINEAR && status >0 && statInner > 0)
  {
    Xyce::dout() << "TWO LEVEL Newton succeeded!" << std::endl;
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::algorithm3_
//
// Purpose       : This is another 2-level algorithm.  For this
//                 algorithm, the de-coupling is greater than for
//                 algorithms 1 and 2.
//
//                 The circuit Newton steps do *not* include any
//                 PDE device elements, and never do.  The PDE
//                 devices are replaced by approximated conductances.
//
// Special Notes : This is the best algorithm.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/25/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::algorithm3_()
{
  int status = -1;
  int statInner = 0;
  bool statusFull = false;
  bool firstOuterStepTaken = false;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl << "Running algorithm 3:" << std::endl;
  }

  nlsPassingPtr_ = 0;

  int twoLevelStep;
  if (status <= 0)
  {
    for (twoLevelStep=0;twoLevelStep<maxOuterSteps_; ++twoLevelStep)
    {
      // Device Only:
      // Solve device problem with continuation:
      twoLevelCouplingMode_       = INNER_PROBLEM;
      outerLoopActiveFlag_= false;
      statInner = 0;

      if (continuationType_ == 1)
      {
        statInner = continuationLoop_ ();
      }
      else if (continuationType_ == 2)
      {
        statInner = locaLoop_ ();
      }
      else
      {
        statInner = nlsInnerPtr_->solve (nlsPassingPtr_);
        nlsPassingPtr_ = 0;
        calcInnerLoopStatistics_ ();
      }

      // Output the nonlinear solver step information:
      if (VERBOSE_NONLINEAR)
      {
        printStepInfo_(twoLevelStep+1,statInner,twoLevelCouplingMode_);
      }

      if (innerLoopFailFatal_)
      {
        // If we can't even get the inner loop to converge,
        // just give up.
        if (statInner <= 0)
        {
          break;
        }
      }

      // Now do all the conductance extractions:
      calcCouplingTerms_ ();

      // ckt-only problem:
      twoLevelCouplingMode_       = OUTER_PROBLEM;
      outerLoopActiveFlag_= true;
      if (firstOuterStepTaken)
      {
        status = nlsOuterPtr_->takeOneSolveStep ();
      }
      else
      {
        firstOuterStepTaken = true;
        status = nlsOuterPtr_->takeFirstSolveStep (nlsInnerPtr_);
      }

      if (noxFlag_)
      {
        nlsPassingPtr_ = nlsOuterPtr_;
      }

      // Output the nonlinear solver step information:
      if (VERBOSE_NONLINEAR)
      {
        printStepInfo_(twoLevelStep+1,status,twoLevelCouplingMode_);
      }

      // check if voltage limiting is still active by getting the norm of the
      // jdxp vector:
      double twoNormJDXP_ = 0.0;
      jdxpVectorPtr_->lpNorm(2, &twoNormJDXP_);

      bool voltLimStat = (twoNormJDXP_ <= voltLimTol_);

      if (VERBOSE_NONLINEAR)
      {
        Xyce::dout() << std::endl
          << "   2-norm of voltage limiting vector: " << twoNormJDXP_ << std::endl;    
      }

      // exit?
      statusFull = false;
      if (status > 0 && statInner>0 && voltLimStat) statusFull = true;

      if (statusFull) break;

    } // end of outer steps "for" loop

  } // end of status if statement

  // Do a final few "full Newton" steps, to get a final consistency between
  // the two solvers.
  int statFinal = 2;
  if (doFullNewtonFinalEnforcement_ && statusFull)
  {
    twoLevelCouplingMode_=FULL_PROBLEM;
    statFinal = nlsOuterPtr_->solve ();

    // Output the nonlinear solver step information:
    if (VERBOSE_NONLINEAR)
    {
      ++twoLevelStep;
      printStepInfo_(twoLevelStep+1,status,twoLevelCouplingMode_);
    }

    // Do one final conductance/capacitance extraction.  Only necessary if 
    // we've done a final "full Newton" solve.  If we haven't done a full Newton 
    // solve, then the most recent solve was the "outer" solve, and attempting to 
    // compute the conductance/capacitance with the "outer" Jacobian will not work.
    calcCouplingTerms_ ();
  }

  if (VERBOSE_NONLINEAR && status >0 && statInner > 0 && statFinal > 0)
  {
    Xyce::dout() << "TWO LEVEL Newton succeeded!" << std::endl;
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::algorithm4_
//
// Purpose       : This is mostly a test algorithm.
//
//                 This is not a two-level class algorithm, but rather a
//                 continuation algorithm to test out
//                 some of the new continuation capabilities.
//
//                 Mostly, I'm putting this here to test out the new
//                 "setParam" function in the device package, which is
//                 needed by LOCA.
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/26/03
//-----------------------------------------------------------------------------
int TwoLevelNewton::algorithm4_()
{
  int status = -1;
  int statNL;
  nlsPassingPtr_ = 0;
  bool continuationLoop = false;
  bool successBool;
  int contMaxTmp = 100;  // this is just a guess for number of steps.
  int stepsLeft = contMaxTmp;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl << "Running algorithm 4:" << std::endl
      << std::endl << "Initial continuation steps: "
      << contMaxTmp << std::endl;
  }

  // The solves should be FULL newton solves.
  twoLevelCouplingMode_       = FULL_PROBLEM;

  // These values assume that any parameter will start at zero and ramp up
  // to a final value.  The final value is the value that is set for this
  // parameter in the netlist.  So, for example, if a voltage source is set
  // to 5.0 volts in the netlist:
  //
  //  vid 4 0 5.00
  //
  // and one of the continuation parameters is vid:v0, then the final value
  // for vid:v0 is 5.00.

  double initVal = 0.0;
  double currVal = initVal;
  double prevVal = initVal;
  double stepSizeEst;

  std::vector<std::string>::iterator iter;
  std::vector<std::string>::iterator begin = paramNameList.begin ();
  std::vector<std::string>::iterator end   = paramNameList.end ();

  std::vector<double>::iterator iterFinalVal;
  std::vector<double>::iterator beginFinalVal = paramFinalVal.begin ();
  std::vector<double>::iterator endFinalVal   = paramFinalVal.end ();

  std::vector<double>::iterator iterCurrentVal;
  std::vector<double>::iterator beginCurrentVal = paramCurrentVal.begin ();
  std::vector<double>::iterator endCurrentVal   = paramCurrentVal.end ();

  // Get the final values for each param., then set each to zero.
  for (iter=begin, iterFinalVal=beginFinalVal, iterCurrentVal=beginCurrentVal;
      iter!=end;
      ++iter, ++iterFinalVal, ++iterCurrentVal)
  {
    //*iterFinalVal = nonlinearEquationLoader_->getParam(*iter);
    *iterFinalVal = 1.0;
    nonlinearEquationLoader_->setParam(*iter, 0.0);
  }

  // now, loop over each parameter, and do a continuation loop for each.
  // For each param, the loop will adjust the param from 0.0 to the final
  // value.  Once each parameter has achieved a final value, the solve is
  // finished.
  for (iter=begin, iterFinalVal=beginFinalVal;
      iter!=end;
      ++iter, ++iterFinalVal)
  {
    // get the initial stepsize:
    stepSizeEst = (*iterFinalVal-0.0)/(static_cast<double>(contMaxTmp));
    currVal = 0.0;
    prevVal = 0.0;
    contStep_=1;

    if (VERBOSE_NONLINEAR)
    {
      Xyce::dout() << "Parameter = " << *iter
        << " finalVal = " << *iterFinalVal << std::endl;
    }

    // begin continuation Loop for the current parameter, *iter:
    bool continuationLoopFinished = false;
    int numTotalFailures = 0;
    while (!continuationLoopFinished)
    {
      bool stepFinished = false;
      int numFailures = 0;

      while(!stepFinished)
      {
        if (stepSizeEst != 0.0)
        {
          stepsLeft = static_cast<int>((*iterFinalVal-currVal)/stepSizeEst) + 1;
        }
        else
        {
          stepsLeft = 1;
        }

        if (VERBOSE_NONLINEAR)
        {
          Xyce::dout() << std::endl << "Continuation Step: " << contStep_
            << "   Estimated Remaining Steps: " << stepsLeft
            << "  " << *iter
            << std::endl
            << "currVal= " << currVal
            << "  prevVal= " << prevVal
            << "  step= " << stepSizeEst
            << std::endl;
        }

        if (stepsLeft < 0)
        {
          Report::UserFatal0() << "Continuation step estimate broken.  Exiting";
        }

        // save a copy of solution:
        (*savedNextSolPtr_) = (*dsPtr_->nextSolutionPtr);

        // set the continuation-parameter:
        nonlinearEquationLoader_->setParam (*iter, currVal);
        *iterCurrentVal = currVal;

        // perform the nonlinear solve:
        statNL = nlsOuterPtr_->solve (nlsPassingPtr_);

        int stepsTaken = nlsOuterPtr_->getNumIterations ();
        int stepsNotTaken = abs(maxOuterSteps_ - stepsTaken);

        nlsPassingPtr_ = 0;

        // Add to the stats:
        calcOuterLoopStatistics_ ();

        // Earlier (in initializeAll), return codes should have been set,
        // hopefully correctly...
        successBool = (statNL > 0);

        if (successBool) // success!
        {
          stepFinished = true;

          //if (numFailures <= 0 && stepsNotTaken > 7)
          if (numFailures <= 0)
          {
            stepSizeEst *= increaseContScalar_;
          }

          --numFailures;
          if (numFailures < 0) numFailures = 0;

          prevVal = currVal;
          currVal += stepSizeEst;

          if ( (*iterFinalVal >= 0 && currVal > *iterFinalVal) ||
              (*iterFinalVal <  0 && currVal < *iterFinalVal) )
          {
            currVal = *iterFinalVal;
            stepSizeEst = currVal - prevVal;
          }

          if (DEBUG_NONLINEAR)
          {
            Xyce::dout() << "\nRight before outputHOMOTOPY:" << std::endl;
            for (int ieric=0;ieric<paramNameList.size();++ieric)
            {
              Xyce::dout() << paramNameList[ieric] << "\t";
              Xyce::dout() << paramCurrentVal[ieric] << std::endl;
            }
          }


          outMgrPtr_->outputHomotopy(pdsMgrPtr_->getPDSComm()->comm(), paramNameList, paramCurrentVal, *dsPtr_->nextSolutionPtr);
        }
        else // failure!
        {
          stepSizeEst *= decreaseContScalar_;

          // restore the solution:
          (*dsPtr_->nextSolutionPtr) = (*savedNextSolPtr_);

          ++numFailures;
          ++numTotalFailures;

          currVal = prevVal + stepSizeEst;
        }
      } // end of stepFinished while loop

      if (!successBool) // failure...
      {
        break;
      }

      ++contStep_;

      continuationLoopFinished =
        ( (*iterFinalVal >= 0 && prevVal >= *iterFinalVal) ||
          (*iterFinalVal <  0 && prevVal <= *iterFinalVal) );

    } // end of continuation loop.

    if (VERBOSE_NONLINEAR)
    {
      Xyce::dout() << "currVal= " << currVal
        << "  prevVal= " << prevVal
        << std::endl
        << std::endl
        << "Total number of failures = " << numTotalFailures << std::endl
        << "Number of actual steps   = " << contStep_-1
        << std::endl;
    }

  } // end of parameter loop.

  status = statNL;

  return status;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::algorithm5_
//
// Purpose       : Same as algorithm0, but calls the inner solver rather
//                 than the outer.
//
// Special Notes : This is best used if you want to use LOCA most of the
//                 time, but not on the nonlinear poisson solve
//                 (firstDCOPstep).
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/07/03
//-----------------------------------------------------------------------------
int TwoLevelNewton::algorithm5_()
{
  int status = -1;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl << "Running algorithm 5:" << std::endl;
  }

  // set up the local algorithm variable:
  int algorithm = twoLevelAlgorithm_;
  if ( externalAnalysisMode ==2) algorithm = twoLevelAlgorithmTran_;

  twoLevelCouplingMode_=FULL_PROBLEM;

  status = nlsInnerPtr_->solve ();

  // Now do conductance/capacitance extractions, if it makes sense.
  // if (!nl_poisson_dcop)
  // {
  //   calcCouplingTerms_ ();
  // }

  return status;
}


//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::solve
// Purpose       :
// Special Notes : Doesn't work with NOX yet...
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
int
  TwoLevelNewton::solve(
      NonLinearSolver *     nlsTmpPtr)
  {
    int status = -1;

    bool nl_poisson_dcop = analysisManager_->getDCOPSolve() == Analysis::NL_POISSON;

    zeroInnerLoopStatistics_ ();

    // set up the local algorithm variable:
    int algorithm = twoLevelAlgorithm_;
    if ( externalAnalysisMode ==2) algorithm = twoLevelAlgorithmTran_;

    // Some intitial setup:
    if (!numInterfaceNodesSetup_)
    {
      numInterfaceNodesSetup_ = true;
      nonlinearEquationLoader_->getNumInterfaceNodes (numInterfaceNodes_);
      numSubProblems_ = numInterfaceNodes_.size ();

      if (VERBOSE_NONLINEAR)
      {
        Xyce::dout() << std::endl
          << "numSubProblems_ = " << numSubProblems_ << std::endl;
      }
    }

    // Algorithm 0:
    // Full Newton, same as if two level was never called.
    if (algorithm == 0 || nl_poisson_dcop)
    {
      status = algorithm0_(nl_poisson_dcop);
    }
    // else if algorithm 1:
    // outter loop is full Newton, inner loop PDE device only.
    else if (algorithm == 1)
    {
      status = algorithm1_ ();
    }
    // same as algorithm 1, but with continuation applied to the inner loop.
    else if (algorithm == 2)
    {
      status = algorithm2_ ();
    }
    // outter loop is circuit only, inner loop PDE device only, with
    // continuation.
    else if (algorithm == 3)
    {
      status = algorithm3_ ();
    }
    else if (algorithm == 4)
    {
      status = algorithm4_ ();
    }
    else if (algorithm == 5)
    {
      status = algorithm5_ ();
    }
    else
    {
      Report::UserFatal0()
        << "Two-Level Newton Algorithm set to invalid number.";
    }

    // Sometimes, when trying to debug, you want the code to exit the first
    // time the solver bombs.  The totalSolveFailFlag gives that option.
    // Right before exiting, however, it instructs the device package to
    // output all its stuff. (tecplot files, text files, etc.)
    if (totalSolveFailFatal_ && status<=0)
    {
#ifndef Xyce_PARALLEL_MPI
      nonlinearEquationLoader_->getLoader().outputPlotFiles();
#endif

      Report::UserFatal0()
        << "Two-Level Newton Algorithm failed to converge.  Exiting.";
    }

    return status;
  }

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::calcCouplingTerms_
//
// Purpose       : The purpose of this function is to manage the extraction
//                 of conductances (coupling terms) from each PDE device.
//                 These conductances may be neccessary if the outter Newton
//                 loop is to be a ckt-only loop.  In that case, the
//                 PDE devices are replaced by much smaller conductance
//                 based models.
//
// Special Notes : Performing a conductance calculation on a PDE device
//                 requires a chain rule calculation, which includes several
//                 linear solves of the Jacobian.  It is because of the
//                 required linear solve that some of the work is done
//                 up here in the nonlinear solver, instead of having all
//                 of it be done by the individual PDE devices.
//
//                 There may be a better solution later, but for the time
//                 being I didn't want linear solves to be performed down
//                 in the device package.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool TwoLevelNewton::calcCouplingTerms_ ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  char filename1[256];

  for (int ich = 0; ich < 256; ++ich)
  { filename1[ich] = 0; }

  // If the current mode is not "INNER_PROBLEM", then we need to re-load
  // the matrix, as the ckt part of the matrix needs to be kept out of
  // this calculation.  Also, re-load Jacobian to make sure it is
  // completely up-to-date with the current solution.
  TwoLevelNewtonMode savedMode = twoLevelCouplingMode_;
  if (twoLevelCouplingMode_ != INNER_PROBLEM)
  {
    twoLevelCouplingMode_ = INNER_PROBLEM;
  }
  nonlinearEquationLoader_->loadJacobian ();

  if (DEBUG_NONLINEAR)
  {
    sprintf(filename1,"%s","tmpJac.txt");
    Linear::Matrix *A = lasSysPtr_->getJacobianMatrix();
    A->writeToFile(filename1);
  }

  // save a copy of the RHS, as it is going to get mangled.
  Linear::Vector *rhsVecPtr = lasSysPtr_->getRHSVector();
  Linear::Vector *newtVecPtr = nlsOuterPtr_->NewtonVectorPtr_;
  *savedRHSPtr_ = *rhsVecPtr;

  bool firstLinearSolve=true;
  lasSolverRCPtr_ = nlsOuterPtr_->lasSolverRCPtr_;

  // Loop over each sub-problem (PDE device).
  // Within each sub-problem, loop over each coupling term. (electrode)
  for (int iSubProblem=0; iSubProblem<numSubProblems_; ++iSubProblem)
  {
    int iCouple;
    int numCoupleTerms = numInterfaceNodes_[iSubProblem];

    if (VERBOSE_NONLINEAR)
    {
      Xyce::dout() << "\n  numCoupleTerms = " << numCoupleTerms << std::endl;
    }

    for (iCouple=0;iCouple<numCoupleTerms;++iCouple)
    {
      // first zero out the RHS vector
      rhsVecPtr->putScalar(0.0);

      // load RHS vector with dFdV.
      tmpBool = nonlinearEquationLoader_->loadCouplingRHS (iSubProblem, iCouple, rhsVecPtr );
      bsuccess = bsuccess &&  tmpBool;

      if (DEBUG_NONLINEAR)
      {
        sprintf(filename1,"dfdv%02d.txt", iCouple);
        rhsVecPtr->writeToFile(filename1);
      }

      // solve linear system to get dXdV.
      int solutionStatus;
      if(firstLinearSolve)
      {
        solutionStatus = lasSolverRCPtr_->solve(false);
        firstLinearSolve=false;

        Util::Param param( "Refactored", 0 );
        lasSolverRCPtr_->getInfo( param );
        if( param.getImmutableValue<int>() ) ++numJacobianFactorizations_;
      }
      else
      {
        solutionStatus = lasSolverRCPtr_->solve(reuseFactors_);
        firstLinearSolve=false;
      }

      if( solutionStatus ) 
      {
        ++numFailedLinearSolves_;
        tmpBool = false;
      }
      else
      {
        tmpBool = true;
      }
      bsuccess = bsuccess && tmpBool;

      if (!tmpBool) 
      {
        Report::UserFatal0()
          << "Linear solve for conductance extraction failed.  iSubProblem = " << iSubProblem
          << " solutionStatus = " << solutionStatus
          << " Exiting.";
      }
      else
      {
        if (VERBOSE_NONLINEAR)
        {
          Xyce::dout() 
            << "Linear solve for conductance extraction succeeded.  iSubProblem = " << iSubProblem
            << " solutionStatus = " << solutionStatus<<std::endl;
        }
      }

      if (!tmpBool) ++numFailedLinearSolves_;
      totalLinearSolveTime_ += lasSolverRCPtr_->solutionTime();
      numLinearSolves_ += 1;

      if (DEBUG_NONLINEAR)
      {
        sprintf(filename1,"dxdv%02d.txt", iCouple);
        newtVecPtr->writeToFile(filename1);
      }

      // copy the newton vector (result of linear solve) into the
      // RHSvector.  Doing this b/c device package knows about rhs, but
      // doesn't know about newton.
      *rhsVecPtr = *newtVecPtr;

      // instruct each PDE device to finish the conductance calc.
      tmpBool = nonlinearEquationLoader_->calcCouplingTerms (iSubProblem, iCouple, rhsVecPtr);
      bsuccess = bsuccess && tmpBool;
    }
  }

  // now restore everything, just in case:
  *rhsVecPtr = *savedRHSPtr_;
  twoLevelCouplingMode_ = savedMode;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::continuationLoop_
//
// Purpose       : Boundary condition continuation loop:
//
// Special Notes : This continuation loop allows for variable step size.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/06/02
//-----------------------------------------------------------------------------
int TwoLevelNewton::continuationLoop_ ()
{
  int statInner;
  bool successBool;
  int contMaxTmp;

  // Instruct the PDE devices to set up the incremental
  // boundary conditions
  int suggestedSteps = 10; // this is a guess to the number of steps.
  suggestedSteps = nonlinearEquationLoader_->enablePDEContinuation ();

  if (suggestedSteps < 1) suggestedSteps = 1;
  contMaxTmp = suggestedSteps;

  double stepSizeEst = 1.0/(static_cast<double>(contMaxTmp));
  double currentAlpha = 0.0;
  double previousAlpha = 0.0;
  int stepsLeft = contMaxTmp;

  // If the continuation loop has never been called before, then leave the
  // initial alpha at zero, because we have no solution to start from,
  // probably.  If it has been called before, go ahead and take the first
  // step, as it should just be one step beyond an already obtained
  // solution.
  if (continuationCalledBefore_)
  {
    currentAlpha = stepSizeEst;
  }
  continuationCalledBefore_ = true;

  contStep_=1;
  bool continuationLoopFinished = false;

  int numTotalFailures = 0;

  while (!continuationLoopFinished)
  {
    bool stepFinished = false;
    int numFailures = 0;

    while(!stepFinished)
    {

      stepsLeft = static_cast<int>((1.0-currentAlpha)/stepSizeEst) + 1;

      if (VERBOSE_NONLINEAR)
      {
        Xyce::dout() << std::endl << "Continuation Step: " << contStep_
          << "   Estimated Remaining Steps: " << stepsLeft
          << std::endl
          << "current alpha = " << currentAlpha
          << "  prev. alpha = " << previousAlpha
          << "  step = " << stepSizeEst
          << std::endl;
      }

      if (stepsLeft < 0)
      {
        Report::UserFatal0()
          << "Continuation step estimate broken.  Exiting.";
      }

      // save a copy of solution:
      (*savedNextSolPtr_) = (*dsPtr_->nextSolutionPtr);

      std::string paramName = "pdealpha";
      nonlinearEquationLoader_->setParam (paramName, currentAlpha);

      // perform the nonlinear solve:
      statInner = nlsInnerPtr_->solve (nlsPassingPtr_);
      nlsPassingPtr_ = 0;
      calcInnerLoopStatistics_ ();

      // Earlier (in initializeAll), return codes for the inner solver
      // should have been set so that "nearConvergence" is not considered
      // an adequate success.  In other words, it returns a negative
      // number, not positive.

      if (DEBUG_NONLINEAR)
      {
        Xyce::dout() << "Status of inner loop solve:  " << statInner << std::endl;
      }

      successBool = (statInner > 0);

      if (successBool) // success!
      {
        stepFinished = true;

        if (numFailures <= 0)
        {
          stepSizeEst *= increaseContScalar_;
        }

        --numFailures;
        if (numFailures < 0) numFailures = 0;

        previousAlpha = currentAlpha;
        currentAlpha += stepSizeEst;

        if (currentAlpha > 1.0)
        {
          currentAlpha = 1.0;
          stepSizeEst = currentAlpha - previousAlpha;
        }
      }
      else // failure!
      {
        stepSizeEst *= decreaseContScalar_;

        // restore the solution:
        (*dsPtr_->nextSolutionPtr) = (*savedNextSolPtr_);

        ++numFailures;
        ++numTotalFailures;

        currentAlpha = previousAlpha + stepSizeEst;
      }
    } // end of stepFinished while loop

    if (!successBool) // failure...
    {
      break;
    }

    ++contStep_;
    continuationLoopFinished = (previousAlpha >= 1.0);

  } // end of continuation loop.

  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << "current alpha = " << currentAlpha
      << "  previous alpha = " << previousAlpha
      << std::endl
      << std::endl
      << "Total number of failures = " << numTotalFailures << std::endl
      << "Number of actual steps   = " << contStep_-1
      << std::endl;
  }

  nonlinearEquationLoader_->disablePDEContinuation ();
  contStep_=0;

  return statInner;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::locaLoop_
//
// Purpose       : Boundary condition continuation (using LOCA) loop:
//
// Special Notes : This continuation loop allows for variable step size.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/26/04
//-----------------------------------------------------------------------------
int TwoLevelNewton::locaLoop_ ()
{
  int statInner;

  int suggestedSteps = nonlinearEquationLoader_->enablePDEContinuation ();

  statInner = nlsInnerPtr_->solve (nlsPassingPtr_);

  nlsPassingPtr_ = 0;
  calcInnerLoopStatistics_ ();

  nonlinearEquationLoader_->disablePDEContinuation ();

  return statInner;
}

//-----------------------------------------------------------------------------
// Function      : TwoLevelNewton::enableSensitivity
//
// Purpose       : This re-sets the code for a sensitivity calculation.
//                 Mainly, it loads the jacobian and rhs vectors in FULL
//                 newton mode, if neccessary.
//
// Special Notes : Since the implementation of the "fullNewtonEnforce"
//                 parameter, this has become less crucial.  When this option
//                 is enabled (the default), then the last loads prior to this
//                 function being called were full newton loads, so there is
//                 no need to re-do them.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/03
//-----------------------------------------------------------------------------
bool TwoLevelNewton::enableSensitivity ()
{
  if (VERBOSE_NONLINEAR)
  {
    Xyce::dout() << std::endl
      << Xyce::section_divider << std::endl
      << "TwoLevelNewton::enableSensitivity " << std::endl;
  }

  bool bsuccess = true;
  bool tmpBool = true;

  twoLevelCouplingMode_=FULL_PROBLEM;
  tmpBool = NonLinearSolver::rhs_ (); bsuccess = bsuccess && tmpBool;
  tmpBool = NonLinearSolver::jacobian_ (); bsuccess = bsuccess && tmpBool;

  if (VERBOSE_NONLINEAR) 
  {
    // print out the norm info for this "full newton" residual:
    double maxNormRHS_=0, twoNormRHS_ = 0.0;
    rhsVectorPtr_->infNorm(&maxNormRHS_);
    rhsVectorPtr_->lpNorm(2, &twoNormRHS_);
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << std::endl;
    Xyce::dout() << "Max. norm of full Newton RHS: " << maxNormRHS_ << std::endl;
    Xyce::dout() << "   2-norm of full Newton RHS: " << twoNormRHS_ << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  if (DEBUG_NONLINEAR)
  {
    static int callsSens = 0;
    char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
    char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;

    sprintf(filename1, "matrixTmp%d.txt",callsSens);
    Linear::Matrix *A = lasSysPtr_->getJacobianMatrix();
    A->writeToFile(filename1);


    Linear::Vector *b = lasSysPtr_->getRHSVector();
    sprintf(filename2, "rhsTmp%d.txt", callsSens);
    int size, i;
    FILE *fp1;

#ifndef Xyce_PARALLEL_MPI
    fp1 = fopen(filename2,"w");
    size = b->globalLength();
    for (i=0;i<size;++i)
    {
      double output = b->getElementByGlobalIndex(i);
      fprintf(fp1,"%25.18e\n",output);
    }
    fclose(fp1);
#endif

    Linear::Vector *x = dsPtr_->nextSolutionPtr;
    sprintf(filename2, "solTmp%d.txt", callsSens);

#ifndef Xyce_PARALLEL_MPI
    fp1 = fopen(filename2,"w");
    size = x->globalLength();
    for (i=0;i<size;++i)
    {
      double output = x->getElementByGlobalIndex(i);
      fprintf(fp1,"%25.18e\n",output);
    }
    fclose(fp1);
#endif

    ++callsSens;
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with 2-level. (for the outer loop)
//-----------------------------------------------------------------------------
bool TwoLevelNewton::isFirstContinuationParam() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Dummy function ... for now.  This probably needs to get revamped.
//-----------------------------------------------------------------------------
int TwoLevelNewton::getContinuationStep() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function ... for now.  This probably needs to get revamped.
//-----------------------------------------------------------------------------
int TwoLevelNewton::getParameterNumber () const
{
  return 0;
}

} // namespace Nonlinear
} // namespace Xyce

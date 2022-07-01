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

//-------------------------------------------------------------------------
//
// Purpose        : Body file for the implemenation of the Newton trust-region
//                  related methods.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/28/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Standard Includes   ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Solver.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_DampedNewton.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_ParamMgr.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_PDS_Comm.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>
#include <N_UTL_Math.h>

// ---------- Static Initializations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : DampedNewton::DampedNewton
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
DampedNewton::DampedNewton(
  const IO::CmdParse &  cp)
  : NonLinearSolver (cp),
    nlParams(DC_OP,cp),
    basicNewton_(true),
    isNormRHS_NaN_(false),
    normRHS_(0.0),
    maxNormRHS_(0.0),
    maxNormRHSindex_(-1),
    normRHS_init_(0.0),
    normRHS_old_(0.0),
    wtNormDX_(0.0),
    normSoln_(0.0),
    stepLength_(1.0),
    nlStep_(0),
    newtonStep_(0),
    searchStep_(0),
    searchDirectionPtr_(0),
    iNumCalls_(0),
    firstTime(true),
    initialDeltaXTol(0.0),
    etaOld(0.1),
    nlResNormOld(0.0),
    tmpConvRate(0.0),
    linearStatus_(true),
    count(0)
{
  nonlinearParameterManager_ = new ParamMgr (commandLine_);

  resetCountersAndTimers_();
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::~DampedNewton
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
DampedNewton::~DampedNewton()
{
  if (!basicNewton_)
  {
    delete searchDirectionPtr_;
  }
  delete nonlinearParameterManager_;

}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::setOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/29/00
//-----------------------------------------------------------------------------
bool DampedNewton::setOptions(const Util::OptionBlock & OB)
{
  bool bsuccess = nlParams.setOptions(OB);

  nonlinearParameterManager_->addParameterSet(DC_OP, nlParams);
  nonlinearParameterManager_->addParameterSet(DC_SWEEP, nlParams);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::setTranOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/05/01
//-----------------------------------------------------------------------------
bool DampedNewton::setTranOptions(const Util::OptionBlock & OB)
{
  NLParams nlTranParams(TRANSIENT, commandLine_);
  bool bsuccess = nlTranParams.setOptions(OB);

  nonlinearParameterManager_->addParameterSet(TRANSIENT, nlTranParams);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::setHBOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/03/2009
//-----------------------------------------------------------------------------
bool DampedNewton::setHBOptions(const Util::OptionBlock & OB)
{
  NLParams nlHBParams(HB_MODE, commandLine_);
  bool bsuccess = nlHBParams.setOptions(OB);

  nonlinearParameterManager_->addParameterSet(HB_MODE, nlHBParams);
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DampedNewton::setNLPOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/07/2015
//-----------------------------------------------------------------------------
bool DampedNewton::setNLPOptions(const Util::OptionBlock & OB)
{
  NLParams nlNLPParams(DC_NLPOISSON, commandLine_);
  bool bsuccess = nlNLPParams.setOptions(OB);

  nonlinearParameterManager_->addParameterSet(DC_NLPOISSON, nlNLPParams);
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::initializeAll
// Purpose       : Called after all register and set functions.
//                 Once the various registrations have taken place,
//                 this function sets the remaining pointers.
// Special Notes :
// Scope         : public
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
// Creation Date : 1/31/02
//-----------------------------------------------------------------------------
bool DampedNewton::initializeAll()
{
  bool bsuccess = NonLinearSolver::initializeAll();

  // make sure the current nlParams is correct.
  nonlinearParameterManager_->getCurrentParams(nlParams);

  if ( nlParams.getNLStrategy() != NEWTON )
  {
    basicNewton_ = false;

    searchDirectionPtr_ = lasSysPtr_->builder().createVector();
  }
  else
  {
    if (!basicNewton_)
    {
      delete searchDirectionPtr_;
    }
    basicNewton_ = true;
    searchDirectionPtr_ = NewtonVectorPtr_;
  }

  if (!getMatrixFreeFlag())
  {
    bsuccess = bsuccess && (jacobianMatrixPtr_ != 0);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::printHeader_
// Purpose       : Print out header for step information.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
void DampedNewton::printHeader_(std::ostream &os)
{
  os << std::endl
     << "  Iter           Step         Wt DX        Inf-Norm      2-Norm (rel)\n"
     << "  -------------------------------------------------------------------\n";
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::printFooter_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/11
//-----------------------------------------------------------------------------
void DampedNewton::printFooter_(std::ostream &os)
{
  os << Xyce::section_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::printStepInfo_
// Purpose       : Print out Newton step information
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
void DampedNewton::printStepInfo_(std::ostream &os, int step)
{
  os << "Niter: " << std::fixed << step 
    << "     " << std::setprecision(4) << std::scientific << stepLength_ 
    << "     " << std::setprecision(4) << std::scientific << wtNormDX_ 
    << "     " << std::setprecision(4) << std::scientific << maxNormRHS_ 
    << "     " << std::setprecision(4) << std::scientific << normRHS_rel_ 
    << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::updateWeights_
// Purpose       : Updates the error weight vector using the previous solution
//                 values as "typical" and the absolute and relative
//                 tolerances.
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
void DampedNewton::updateWeights_()
{
  // ***** Weighting Vector for Solution *****

  // On the first call to the nonlinear solver,
  // weigh based on the minimums.
  double solnNorm = 0.0;
  dsPtr_->nextSolutionPtr->infNorm(&solnNorm);

  if ((iNumCalls_ == 0) && (solnNorm <= Util::MachineDependentParams::DoubleMin()))
  {
    solWtVectorPtr_->putScalar(nlParams.getRelTol() + nlParams.getAbsTol());
  }
  else
  {
    double newSoln, oldSoln;

    int length = dsPtr_->nextSolutionPtr->localLength();
    for (int i = 0; i < length; ++i)
    {
      newSoln = (*dsPtr_->nextSolutionPtr)[i];
      oldSoln = (*dsPtr_->currSolutionPtr)[i];

      (*(solWtVectorPtr_))[i] =
        nlParams.getRelTol() * std::max(fabs(oldSoln), fabs(newSoln)) +
        nlParams.getAbsTol();
    }
  }

  if (nlParams.getMaskingFlag())
  {
    int length = dsPtr_->nextSolutionPtr->localLength();
    Linear::Vector& mask = *(lasSysPtr_->getDeviceMaskVector());

    for (int i = 0; i < length; ++i)
    {
      if (mask[i] == 0.0)
        (*(solWtVectorPtr_))[i] = Util::MachineDependentParams::MachineBig();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::solve
// Purpose       : solve implements a damped Newton method nonlinear equations
//                 solver.
// Special Notes : The implementation is based upon the original damped Newton
//                 method implemented in Xyce(TM) by Eric Keiter.
//
// Selected variables and methods:
//
// nextSolVector      the approximate solution to the nonlinear equation
//                    The rhsVector and jacobianMatrix load nextSolVector.
//
// equateTmpVectors_  tmpSolVector := nextSolVector
// updateTmpSol_      tmpSolVector := nextVector + searchDirection*stepLength_
//
// switchTmpVectors   swap &tmpSolVector and &nextSolVector, makes it possible
//                    for the nonlinear solver to preserve the entire
//                    circuit state associated with nextSolVector.
//
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
int DampedNewton::solve(NonLinearSolver * nlsTmpPtr)
{
  static const char *trace = "DampedNewton::solve: ";

  if (DEBUG_NONLINEAR)
  {
    setDebugFlags(getAnalysisManager().getStepNumber() + 1, getAnalysisManager().getTime());

    if (debugTimeFlag_ && nlParams.getDebugLevel() > 0 )
      Xyce::dout() << retCodes_;
  }

  resetCountersAndTimers_();

  // Get the current analysis mode
  AnalysisMode mode1 = nonlinearParameterManager_->getAnalysisMode();

  // Change the nlparams so that they are appropriate for the current
  // time integration mode, if neccessary.
  nonlinearParameterManager_->getCurrentParams(nlParams);

  // Output the nonlinear solver information header:
  if (VERBOSE_NONLINEAR)
  {
    printHeader_(Xyce::lout());
  }

  // For the initial RHS load, the step number needs to be zero.  The Xyce
  // device package needs to know this - either in the event that it might want
  // to set initial conditions.
  nlStep_ = newtonStep_ = 0;

  // Recall that prior to this solver being called, the new solution has been
  // predicted and resides in nextSolVectorPtr.
  rhs_();

  // Relative RHS norm (normRHS_rel_).  On the first evaluation this is one
  // since normRHS_init = normRHS.
  normRHS_rel_ = 1.0;

  // Weighted norm of the change (wtNormDX_) is set to zero initially.
  wtNormDX_    = 0.0;

  if (VERBOSE_NONLINEAR)
  {
    // Max RHS norm (maxNormRHS_).
    if (nlParams.getMaskingFlag())
    {
      rhsVectorPtr_->wMaxNorm(*(getPNormWeights()), &maxNormRHS_, &maxNormRHSindex_);
    }
    else
    { 
      rhsVectorPtr_->infNorm(&maxNormRHS_, &maxNormRHSindex_);
    }

    // Print out the starting point information.
    printStepInfo_(Xyce::lout(), nlStep_);
  }
  
  // Used to calculate resConvRate which is used by converged_(). This changes
  // each nonlinear iteration to be the norm of the RHS from the previous
  // iteration.
  normRHS_old_ = normRHS_;

  // Used to calculate normRHS_rel_, which is output by printStepInfo_().
  normRHS_init_ = normRHS_;

  // Update the error weighting vector for use with the weighted norms.
  if (mode1 == TRANSIENT)
    updateWeights_();

  // Nonlinear loop. Loop continues until convergedStatus is nonzero. A
  // positive value indicates that the method has converged. A negative value
  // that the method has failed. This is updated using the subroutine
  // converged_() at the end of each nonlinear solver iteration.
  int convergedStatus = 0;
  while (convergedStatus == 0)
  {
    // Increment step counters.  This needs to be updated *before* the rhs_
    // call, which is inside of computeStepLength_.
    nlStep_++;

    // Calculate the Jacobian for the current iterate.
    jacobian_();

    if (DEBUG_NONLINEAR && !getMatrixFreeFlag())
    {
      debugOutput1( *(lasSysPtr_->getJacobianMatrix()), *(lasSysPtr_->getRHSVector()));
    }

    // Calculate a direction.
    direction_();

    if (!basicNewton_)
    {
      setX0_();
    }

    // Calculate a step length (damping or backtracking for Newton, line search
    // for others...) and take that step and recalulate the RHS.
    computeStepLength_();

    if (DEBUG_NONLINEAR)
    {
      debugOutput3 (*dsPtr_->nextSolutionPtr, *searchDirectionPtr_ );
    }

    // Test for convergence based on the weighted norms of the RHS and DX.
    // Also update convergence rate and counters.

    // Update the error weighting vector for use with the weighted norms.
    if (mode1 != TRANSIENT)
      updateWeights_();

    // Scale the search direction so that we can calculate some norms on it and
    // use them in the convergence tests and output.
    if (!basicNewton_)
    {
      searchDirectionPtr_->scale(stepLength_);
    }

    // Check our iteration status. Returns positive if done, negative if error,
    // zero otherwise.
    convergedStatus = converged_();

    if (VERBOSE_NONLINEAR)
    {
      printStepInfo_(Xyce::lout(), nlStep_);
    }

    // Increment diagnostic step counters.  These need to be updated after
    // direction_, in which the current direction.
    newtonStep_++;

  } // while (convergedStatus == 0)

  // Increment the number of calls to the nonlinear solver.
  iNumCalls_++;

  if (VERBOSE_NONLINEAR)
  {
    printFooter_(Xyce::lout());
  }

  // A positive value of convergedStatus indicates a successful nonlinear
  // solutions. Otherwise, there was an error.
  return convergedStatus;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::takeFirstSolveStep
//
// Purpose       : This is the same as the function "solve", except that
//                 it only does the various initializations at the top
//                 of the "::solve" function, and only takes the first
//                 NL step.
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
int DampedNewton::takeFirstSolveStep(NonLinearSolver * nlsTmpPtr)
{
  static const char *trace = "DampedNewton::takeFirstSolveStep: ";

  if (DEBUG_NONLINEAR)
  {
    setDebugFlags(getAnalysisManager().getStepNumber() + 1, getAnalysisManager().getTime());
  }

  // resolve the counters issue later.  It may be neccessary to save the
  // counters for later.  Or not.
  resetCountersAndTimers_();

  // Change the nlparams so that they are appropriate for the current
  // time integration mode, if neccessary.
  nonlinearParameterManager_->getCurrentParams(nlParams);

  // Output the nonlinear solver information header:
  if (VERBOSE_NONLINEAR)
  {
    printHeader_(Xyce::lout());
  }

  // For the initial RHS load, the step number needs to be zero.  The Xyce
  // device package needs to know this - either in the event that it might want
  // to set initial conditions.
  nlStep_ = newtonStep_ = 0;

  // Recall that prior to this solver being called, the new solution has been
  // predicted and resides in nextSolVectorPtr.
  rhs_();

  // Relative RHS norm (normRHS_rel_).  On the first evaluation this is one
  // since normRHS_init = normRHS.
  normRHS_rel_ = 1.0;

  // Weighted norm of the change (wtNormDX_) is set to zero initially.
  wtNormDX_    = 0.0;

  if (VERBOSE_NONLINEAR)
  {
    // Max RHS norm (maxNormRHS_).
    rhsVectorPtr_->infNorm(&maxNormRHS_);

    // Print out the starting point information.
    printStepInfo_(Xyce::lout(), nlStep_);
  }

  // Check to ensure that we haven't been passed a converged solution already
  // (i.e., normRHS < absTol)
  if (normRHS_ < Util::MachineDependentParams::MachineEpsilon()) return  1;

  // Used to calculate resConvRate which is used by converged_(). This changes
  // each nonlinear iteration to be the norm of the RHS from the previous
  // iteration.
  normRHS_old_ = normRHS_;

  // Used to calculate normRHS_rel_, which is output by printStepInfo_().
  normRHS_init_ = normRHS_;

  // The variable convergedStatus: positive value indicates that
  // the method has converged. A negative value
  // that the method has failed. This is updated using the subroutine
  // converged_() at the end of each nonlinear solver iteration.
  int convergedStatus = 0;

  // Increment step counters.  This needs to be updated *before* the rhs_
  // call, which is inside of computeStepLength_.
  nlStep_++;

  // Calculate the Jacobian for the current iterate.
  jacobian_();

  if (DEBUG_NONLINEAR && !getMatrixFreeFlag())
  {
    debugOutput1( *(lasSysPtr_->getJacobianMatrix()), *(lasSysPtr_->getRHSVector()));
  }

  // Calculate a direction.
  direction_();

  if (!basicNewton_)
  {
    setX0_();
  }

  // Calculate a step length (damping or backtracking for Newton, line search
  // for others...) and take that step and recalulate the RHS.
  computeStepLength_();

  if (DEBUG_NONLINEAR)
  {
    debugOutput3 (*dsPtr_->nextSolutionPtr, *searchDirectionPtr_ );
  }

  // Test for convergence based on the weighted norms of the RHS and DX.
  // Also update convergence rate and counters.  1. Update the error
  // weighting vector for use with the weighted norms.
  updateWeights_();

  // Scale the search direction so that we can calculate some norms on it and
  // use them in the convergence tests and output.
  if (!basicNewton_)
  {
    searchDirectionPtr_->scale(stepLength_);
  }

  // Check our iteration status. Returns positive if done, negative if error,
  // zero otherwise.
  convergedStatus = converged_();

  if (VERBOSE_NONLINEAR)
  {
    printStepInfo_(Xyce::lout(), nlStep_);
  }

  // Increment diagnostic step counters.  These need to be updated after
  // direction_, in which the current direction.
  newtonStep_++;

  // Increment the number of calls to the nonlinear solver.
  iNumCalls_++;

  // A positive value of convergedStatus indicates a successful nonlinear
  // solutions. Otherwise, there was an error.
  return convergedStatus;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::takeOneSolveStep
//
// Purpose       : Just like "::initializeSolve" is essentially the top
//                 of the "::solve" function, this function is similar
//                 to the bottom of the "::solve" function.  This
//                 function also only takes a single NL step.
//
//                The main differences between this function and solve are:
//                  1) There is no while loop. - this is only a single
//                      step.
//
//                  2) Most tasks performed prior to the while loop are not
//                     performed here, as they are initialization tasks.
//                     These include:

//                      - step counters, such as "nlStep_" and "newtonStep_"
//                         are not reset to zero.
//
//                      - The function "resetCountersAndTimers" is not called.
//                        This may be fixed up later.
//
//                      -  The resetting of the deltaXtol, which is performed
//                         on the "first call" to ::solve, is not done.
//                         I'll have to sort this out later as well.
//
//                      - the loadJacobianFlag is not updated - assumed to
//                         have not changed.
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 01/09/01
//-----------------------------------------------------------------------------
int DampedNewton::takeOneSolveStep()
{
  static const char *trace = "DampedNewton::takeOneSolveStep";

  if (DEBUG_NONLINEAR)
  {
    setDebugFlags(getAnalysisManager().getStepNumber() + 1, getAnalysisManager().getTime());
  }

  // Change the nlparams so that they are appropriate for the current
  // time integration mode, if neccessary.
  nonlinearParameterManager_->getCurrentParams(nlParams);

  // Recall that prior to this solver step being called, a new solution has been
  // predicted and resides in nextSolVectorPtr.
  rhs_();

  // Used to calculate resConvRate which is used by converged_(). This changes
  // each nonlinear iteration to be the norm of the RHS from the previous
  // iteration.
  normRHS_old_ = normRHS_;

  // Nonlinear loop. Loop continues until convergedStatus is nonzero. A
  // positive value indicates that the method has converged. A negative value
  // that the method has failed. This is updated using the subroutine
  // converged_() at the end of each nonlinear solver iteration.
  int convergedStatus = 0;

  // Increment step counters.  This needs to be updated *before* the rhs_
  // call, which is inside of computeStepLength_.
  nlStep_++;

  // Calculate the Jacobian for the current iterate.
  jacobian_();

  if (DEBUG_NONLINEAR && !getMatrixFreeFlag())
  {
    debugOutput1( *(lasSysPtr_->getJacobianMatrix()), *(lasSysPtr_->getRHSVector()));
  }

  // Calculate a direction.
  direction_();

  if (!basicNewton_)
  { 
    setX0_();
  }

  // Calculate a step length (damping or backtracking for Newton, line search
  // for others...) and take that step and recalulate the RHS.
  computeStepLength_();

  if (DEBUG_NONLINEAR)
  {
    debugOutput3 (*dsPtr_->nextSolutionPtr, *searchDirectionPtr_ );
  }

  // Test for convergence based on the weighted norms of the RHS and DX.
  // Also update convergence rate and counters.  1. Update the error
  // weighting vector for use with the weighted norms.
  updateWeights_();

  // Scale the search direction so that we can calculate some norms on it and
  // use them in the convergence tests and output.
  if (!basicNewton_)
  {
    searchDirectionPtr_->scale(stepLength_);
  }

  // Check our iteration status. Returns positive if done, negative if error,
  // zero otherwise.
  convergedStatus = converged_();

  if (VERBOSE_NONLINEAR)
  {
    printStepInfo_(Xyce::lout(), nlStep_);
  }

  // Increment diagnostic step counters.  These need to be updated after
  // direction_, in which the current direction.
  newtonStep_++;

  if (DEBUG_NONLINEAR && !getMatrixFreeFlag() && convergedStatus > 0)
  {
    // ERK Note: this one needs the nl solve step incremented by 1!
//    nlStep_++;
    debugOutput1(*lasSysPtr_->getJacobianMatrix(), *lasSysPtr_->getRHSVector());
//    nlStep_--; // restore original value.
  }

  // Increment the number of calls to the nonlinear solver.
  iNumCalls_++;

  // A positive value of convergedStatus indicates a successful nonlinear
  // solutions. Otherwise, there was an error.
  return convergedStatus;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::updateX_()
//
// Purpose       : Update the value of nextSolVectorPtr_ using
//                 steplength and the searchDirectionPtr_.  On input, it is
//                 assumed that tmpSolVectorPtr_ contains the initial guess for
//                 the nonlinear iteration, i.e., setX0_ was called earlier.
//
// Special Notes : Also modifies errorEstVectorPtr_. Only this and
//                 updateX_() should modify the solution directly.
//
//                 ERK: This is new version of this function, hopefully less
//                 confusing.  The roles of "tmp" and "next" have been
//                 reversed, which removes the need for some of the switch
//                 functions.
//
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
//
//                 Eric Keiter, SNL, Parallel Computational Sciences (9233)
//
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
void DampedNewton::updateX_()
{
  if (!basicNewton_)
    dsPtr_->nextSolutionPtr->update(1.0, *dsPtr_->tmpSolVectorPtr,
                                    stepLength_,
                                    *searchDirectionPtr_, 0.0);
  else
    dsPtr_->nextSolutionPtr->update(1.0, *searchDirectionPtr_);
  
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::rhs_()
// Purpose       : Updates the RHS based on nextSolVectorPtr_ and
//                 calculates normRHS_ based on 2-norm.
// Special Notes : The rhsVectorPtr_ is really the NEGATIVE of F(x).
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool DampedNewton::rhs_()
{
  Teuchos::RCP<Linear::Vector> rhsCopy = Teuchos::rcp( rhsVectorPtr_->cloneCopyVector() );
  bool status = NonLinearSolver::rhs_();

  if (DEBUG_NONLINEAR)
  {
    debugOutput3 (*dsPtr_->nextSolutionPtr, *searchDirectionPtr_ );
  }

  isNormRHS_NaN_ = false;

  if (nlParams.getMaskingFlag())
  {
    int length = dsPtr_->nextSolutionPtr->globalLength();
    rhsVectorPtr_->wRMSNorm(*(getPNormWeights()), &normRHS_);
    normRHS_ *= std::sqrt( length );  // Undo scaling of RMS norm
  }
  else
  {
    rhsVectorPtr_->lpNorm(2, &normRHS_);
  }
 
  if (std::isnan(normRHS_) || std::isinf(normRHS_))
  {
    if ( VERBOSE_NONLINEAR )
    {
      std::vector<int> nanEntries;
      Linear::checkVectorForNaNs( *rhsVectorPtr_, nanEntries );
      if ( nanEntries.size() )
      {
        Xyce::lout() << "NaN/Inf check found " << nanEntries.size() << " entries in the residual vector!" << std::endl;
        for (int i=0; i<nanEntries.size(); ++i)
        {
          Xyce::lout() << "Residual entry: [" << nanEntries[i] << "]" << std::endl;
        }
      }
    }

    isNormRHS_NaN_ = true;
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::newtonDirection_
// Purpose       : Computes the Newton Direction inexactly using
//                 jacobianMatrixPtr_ and the rhsVectorPtr_.
//                 On output, NewtonVectorPtr_ contains the Newton
//                 Direction.
// Special Notes : This is *not* meant to override newton_ in
//                 NonLinearSolver.
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool DampedNewton::newton_()
{
  // Set the solver tolerance based on the residual norm.
  if (nlParams.getForcingFlag()) setForcing_(normRHS_);

  return NonLinearSolver::newton_();
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::direction_
// Purpose       : The function calculates the direction vector used in
//                 the nonlinear solver (e.g., Newton direction).
// Special Notes : If the search direction is a Newton direction, then a linear
//                 system is solved inexactly.  To ensure that the approximate
//                 solution is a descent direction, the projection of the
//                 solution along the gradient is computed to working
//                 precision.  On output, searchDirectionPtr_ points to the
//                 search direction vector.  The vector allocated to store the
//                 NewtonVector is always used to store the
//                 SearchDirectionVector, regardless of the actual algorithm
//                 used to determine the search direction.
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences Department
// Creation Date : 02/21/01
//-----------------------------------------------------------------------------
void DampedNewton::direction_()
{
  static const char *trace = "DampedNewton::direction_: ";

  // --- Calculate SearchDirection ---
  // Compute the Newton direction
  linearStatus_ = newton_();

  // Copy the Newton direction into the search direction
  if (!basicNewton_)
  {
    *searchDirectionPtr_ = *NewtonVectorPtr_;
  }
} // end direction_

//-----------------------------------------------------------------------------
// Function      : DampedNewton::computeStepLength_
// Purpose       : This function calculates the step length for a given search
//                 direction. If the search direction is Newton, it may use one
//                 of a variety of backtracking methods including on suggested
//                 by
// Special Notes : !!!!!NOTE:  In this method, we use the term FULLSTEP to mean
//                 the maximum allowed Newton step (0..1).  fullStep = 1.
//                 corresponds to a Newton step, but when computeStepLength_
//                 is coupled with a constraint back tracking, 0<fullStep<=1.
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences Department
// Creation Date : 02/21/01
//-----------------------------------------------------------------------------
bool DampedNewton::computeStepLength_()

{
  static const char *trace = "DampedNewton::computeStepLength_";

  searchStep_ = 0;  // number of search steps

  // Choose step length approach based on the current solution method.
  switch (nlParams.getSearchMethod())
  {
    case DIVIDE:

      return divide_();
      break;

    case BACKTRACK:

      return backtrack_();
      break;

    default:

      return fullNewton_();
      break;
  }

  // If the code gets to this point in the function there has been an
  // error...
  return false;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::divide_()
// Purpose       : Eric's "Divide and Conquer" method. Subsequently halves
//                 the steplength until there is improvement in the RHS or
//                 the steplength becomes smaller than Util::MachineDependentParams::MachineEpsilon().
// Special Notes : Does not verify that the search direction is a descent
//                 direction.
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/29/01
//-----------------------------------------------------------------------------
bool DampedNewton::divide_()
{
  static const char *trace = "DampedNewton::divide__";

  // Want to get reduction in the norm, so we first copy the norm from
  // the previous iterate
  const double oldNormRHS = normRHS_;

  // The upper bound on the backtracking (damping) is determined by
  // the applied constraints.
  const double fullStep = 1.0;

  // Starting value
  stepLength_ = fullStep;

  // Update the solution.
  updateX_();

  // Evaluate the new residual and take its norm:
  rhs_();

  if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
  {
    Xyce::dout() << "\n\tSearch Step: " << searchStep_ << ", Step Size: " << stepLength_ << std::endl
                 << "\toldNormRHS: " << oldNormRHS << ", normRHS: " << normRHS_ << std::endl;
  }

  // Linesearch...
  bool searchDone = (normRHS_ < oldNormRHS);
  while (!searchDone)
  {
    // Cut the step length
    stepLength_ *= 0.5;

    // If the step length gets too small, try full step and finish.
    if (stepLength_ < Util::MachineDependentParams::MachineEpsilon())
    {
      if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
      {
        Xyce::Report::UserWarning0() << "\tStep size too small: " << stepLength_ << "\n\tTaking a full step.\n";
      }

      stepLength_ = fullStep;
      searchDone = true;
    }

    updateX_();

    // Evaluate the new RHS vector and take norm:
    rhs_();

    searchStep_++;

    // Check to see if we're done:
    searchDone = ((normRHS_ < oldNormRHS) || (searchDone) ||
                  (searchStep_ >= nlParams.getMaxSearchStep()));

    if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
    {
      Xyce::dout() << "\n\tSearch Step: " << searchStep_ << ", Step Size: " << stepLength_ << std::endl
                   << "\toldNormRHS: " << oldNormRHS << ", normRHS: " << normRHS_ << std::endl;
    }
  }

  return (normRHS_ < oldNormRHS);
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::backtrack_()
// Purpose       : Backtracking method based on Dennis & Shnabel.
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences Department
// Creation Date : 07/28/01
//-----------------------------------------------------------------------------
bool DampedNewton::backtrack_()
{
  static const char *trace = "DampedNewton::computeStepLength_";

  // Want to get reduction in the norm, so we first copy the norm from
  // the previous iterate
  const double oldNormRHS = normRHS_;

  // The upper bound on the backtracking (damping) is determined by
  // the applied constraints.
  const double fullStep = 1.0;

  // Starting value
  stepLength_ = fullStep;

  // Update the solution.
  updateX_();

  // Evaluate the new residual and take its norm:
  rhs_();

  // Linesearch...
  double rho, delta, theta, lambda;
  const double t = 1.0e-04, thetaMin = 0.1, thetaMax = 0.5;
  const double minStep = pow(Util::MachineDependentParams::MachineEpsilon(), 0.33);

  // If we're using the inexact-Newton forcing, we use the forcing value
  // (linear-solver convergence tolerance) to set the initial "lambda".
  if (nlParams.getForcingFlag())
    lambda = 1.0 - nlParams.getForcingTerm();
  else
    lambda = 1.0;

  rho = normRHS_ / oldNormRHS;

  while ((rho > 1.0 - t * lambda) &&
         (searchStep_ < nlParams.getMaxSearchStep()))
  {
    delta = rho * rho - 1.0 + 2.0 * lambda;
    if (delta <= 0.0)
      theta = thetaMax;
    else
    {
      theta = lambda / delta;

      if (theta > thetaMax)
        theta = thetaMax;
      else if (theta < thetaMin)
        theta = thetaMin;
    }

    lambda      *= theta;
    stepLength_ *= theta;

    // Jump out if the step-size gets too small...
    if (stepLength_ < minStep)
    {
      stepLength_ = minStep;
      searchStep_ = nlParams.getMaxSearchStep();
    }

    // Update solution & evaluate the new RHS vector and take norm:
    updateX_();
    rhs_();

    rho = normRHS_ / oldNormRHS;
    searchStep_++;

    if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
    {
      Xyce::dout() << "\n\tSearch Step: " << searchStep_ << ", Step Size: " << stepLength_ << std::endl
                   << "\toldNormRHS: " << oldNormRHS << ", normRHS: " << normRHS_ << std::endl;
    }
  }

  return (rho <= 1.0 - t * lambda);

}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::fullNewton_()
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/29/01
//-----------------------------------------------------------------------------
bool DampedNewton::fullNewton_()
{
  // ***** Full Newton method *****
  stepLength_ = 1.0;

  updateX_();

  // Evaluate the new residual and take its norm:
  rhs_();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::converged_
// Purpose       : This method checks for convergence of the Newton solver.
//                 Much of this is currently based on the methods in both
//                 P.N. Brown, A.C. Hindmarsh and L.R. Petzold, "Consistent
//                 initial condition calculations for differential-algebraic
//                 systems", SIAM J. Sci. Comput., Vol. 19, No. 5,
//                 pp. 1495-1512, Sep. 1998 and J.E. Dennis, Jr. and
//                 R.B. Schnabel, "Numerical Methods for Unconstrained
//                 Optimization and Nonlinear Equations", Prentice-Hall, 1983.
//
// Special Notes : Returns 0 if not converged, positive if converged, negative
//                 if not converged but some other error.
//
//                 Eric Keiter, 9233, 03/27/03: I modified the code so that
//                 it would use values in the ReturnCodes class as
//                 the integer values for success/failure, rather than
//                 hardwired numbers.  The defaults of this class are the
//                 same numbers as the old hardwired ones.
//
//                 The reason for doing this is to allow the user, or
//                 calling code to set which circumstances should be
//                 considered adequate convergence.  If
//                 running a two-level continuation loop, within the transient mode,
//                 there needed to be a way to easily disable the
//                 "nearConverged" senario.
//
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences Department
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
int DampedNewton::converged_()

{
  static const char *trace = "DampedNewton::converged_";

  const double machineEpsilon = Util::MachineDependentParams::MachineEpsilon();
  const double convRateMax  = 0.5 * Util::MachineDependentParams::DoubleMax();
  const double minResReduct = 9.0e-1;
  const double stagnationTol = 1.0e-3;

  AnalysisMode mode1 = nonlinearParameterManager_->getAnalysisMode();

  // check if the linear solver failed (-9 return code)
  if (!linearStatus_)
  {
    // Check if NaNs were detected
    std::vector< std::pair<int, int> >& nanEntries = lasSolverRCPtr_->getNaNEntries();
    if ( VERBOSE_NONLINEAR )
    {
      if ( nanEntries.size() )
      {
        Xyce::lout() << "NaN check found " << nanEntries.size() << " entries in the residual or Jacobian!" << std::endl;
        for (int i=0; i<nanEntries.size(); ++i)
        {
          if (nanEntries[i].second >= 0)
            Xyce::lout() << "Jacobian entry: [" << nanEntries[i].first << ", " << nanEntries[i].second << "]" << std::endl;
          else
            Xyce::lout() << "Residual entry: [" << nanEntries[i].first << "]" << std::endl;
        }
      }
    }

    // If NaNs were found in the residual or Jacobian, it is not a linear solver failure.
    if ( nanEntries.size() )
      return retCodes_.nanFail; // default = -6
    else
      return retCodes_.linearSolverFailed; // default: -9
  }

  // Devices need to satisfy their own convergence criteria
  if (nlParams.getEnforceDeviceConvFlag ()) 
  {
    bool allDevicesConverged_ = nonlinearEquationLoader_->allDevicesConverged(pdsMgrPtr_->getPDSComm()->comm());
    if (!allDevicesConverged_ && (nlStep_ < nlParams.getMaxNewtonStep() ))
    {
      return 0;
    }
    
    if (!allDevicesConverged_ && (nlStep_ >= nlParams.getMaxNewtonStep() ))
      return retCodes_.tooManySteps;
  }

  // This test is for 2-level solves ONLY.
  bool innerDevicesConverged_ = nonlinearEquationLoader_->innerDevicesConverged(pdsMgrPtr_->getPDSComm()->comm());
  if (!innerDevicesConverged_ )
  {
    return 0;
  }

  // 2-norm of RHS is used in converged_(), check if it's a NaN
  if (isNormRHS_NaN_)
  {
    return retCodes_.nanFail; // default = -6
  }

  // Max RHS norm (maxNormRHS_) is used in converged_() and output by
  // printStepInfo_().
  if (nlParams.getMaskingFlag())
  {
    rhsVectorPtr_->wMaxNorm(*(getPNormWeights()), &maxNormRHS_, &maxNormRHSindex_);
  }
  else
  { 
    rhsVectorPtr_->infNorm(&maxNormRHS_, &maxNormRHSindex_);
  }

  // This parameter "takes-out" any damping induced in the size of the norm by
  // a line-search or other globalization method.
  // Weighted norm of the change (wtNormDX_) is used in converged_() and
  // output by printStepInfo_().

  searchDirectionPtr_->wMaxNorm(*solWtVectorPtr_, &wtNormDX_);

  // If the RHS norm is so small already, just say we're converged,
  // and move on.  If we don't do this, we may wind up dividing by a
  // zero later. 
  if (normRHS_ < machineEpsilon)
  {
    if (normRHS_old_ < machineEpsilon)
    {
      resConvRate_ = 1.0;
    }
    else
    {
      resConvRate_ = 0.0;
    }
    if (normRHS_init_ < machineEpsilon)
    {
      normRHS_rel_ = 1.0;
    }
    else
    {
      normRHS_rel_ = 0.0;
    }
    normRHS_old_ = normRHS_;
    
    return retCodes_.normTooSmall;
  }
  
  // Relative RHS norm (normRHS_rel_) is output by printStepInfo_().
  normRHS_rel_ = normRHS_ / normRHS_init_;

  // Residual convergence rate (resConvRate) used by converged_(),
  // also printed out below.
  resConvRate_ = normRHS_ / normRHS_old_;

  // Reset old norm for next nonlinear iteration.
  normRHS_old_ = normRHS_;

  double updateSize = wtNormDX_ / stepLength_;

  // Check for "normal" convergence
  if (maxNormRHS_ <= nlParams.getRHSTol() &&
      updateSize <= nlParams.getDeltaXTol())
    return retCodes_.normalConvergence; // 2;

  // Transient "Near Converged"...
  if (nlStep_ >= nlParams.getMaxNewtonStep() &&
      mode1 == TRANSIENT)
  {
    // Check for "near" convergence and let the time integrator handle the
    // error analysis and decide on whether or not to accept its step
    if ((normRHS_rel_ <= minResReduct) && (resConvRate_ <= 1.0))
      return retCodes_.nearConvergence; // 3;
  }

  // Check to see if update is really small
  if (updateSize <= nlParams.getSmallUpdateTol ())
    return retCodes_.smallUpdate; // 4;

  // Next check the number of steps
  if (nlStep_ >= nlParams.getMaxNewtonStep())
    return retCodes_.tooManySteps; // -1;

  // Make sure that we haven't had too big an update
  else if (resConvRate_ > convRateMax)
    return retCodes_.updateTooBig; // -2;

  // Check for a stall in the convergence rate - if it is near one for five
  // consecutive steps, we've stalled so jump out.
  if (mode1 == TRANSIENT && fabs(resConvRate_ - 1.0) <= stagnationTol)
  {
    if (count == 0 || resConvRate_ < tmpConvRate)
      tmpConvRate = resConvRate_;

    count++;
    if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
    {
      Xyce::dout() << "\tcount:\t" <<  count << std::endl
                   << "\tresConvRate_:\t" <<  resConvRate_ <<  "\n" << std::endl;
    }
  }
  else
  {
    count = 0;
  }

  if (mode1 == TRANSIENT && count == 5)
  {
    count = 0;
    if ((normRHS_rel_ < minResReduct) && (tmpConvRate <= 1.0))
    {
      if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
      {
        Xyce::dout() << "\ttmpConvRate: " <<  tmpConvRate << std::endl
                           << "\tnormRHS_rel_:\t" <<  normRHS_rel_ << std::endl
                           << "\tReturning 3\n" << std::endl;
      }

      return retCodes_.nearConvergence; // 3;
    }
    else
    {
      if (DEBUG_NONLINEAR && debugTimeFlag_ && isActive(Diag::NONLINEAR_PARAMETERS) )
      {
        Xyce::dout() << "\ttmpConvRate: " <<  tmpConvRate << std::endl
                     << "\tnormRHS_rel_:\t" <<  normRHS_rel_ << std::endl
                     << "\tReturning -3\n" << std::endl;
      }

      return retCodes_.stalled;  // -3;
    }
  }

  return 0;

}

//-----------------------------------------------------------------------------
// Function      : DampedNewton::setForcing_
// Purpose       : This method calculates the forcing term (i.e., linear
//                 residual tolerance) for iterative solvers based on the
//                 method of Walker and Pernice (RefXXX)
// Special Notes :
// Scope         : private
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences Department
// Creation Date : 05/03/01
//-----------------------------------------------------------------------------
void DampedNewton::setForcing_(const double nlResNorm)
{
  const double etaExp  = 0.5 * (1.0 + sqrt(5.0));
  const double etaMax  = 0.1;
  const double etaMin  = 1.0e-12;
  const double etaInit = 1.0e-01;
  double eta;
  double etaSafe;

  if (newtonStep_ == 0)
  {
    eta          = etaInit;
    nlResNormOld = nlResNorm;
  }

  else if (nlResNormOld > Util::MachineDependentParams::DoubleMin())
  {
    Util::Param linRes( "RESIDUAL", 0.0 );
    lasSolverRCPtr_->getInfo( linRes );
    eta = fabs(nlResNorm - linRes.getImmutableValue<double>()) / nlResNormOld;
    eta *= eta;

    // First safeguard...
    etaSafe = pow(etaOld, etaExp);
    if (etaSafe > 0.1)
      eta = std::max(eta, etaSafe);

    // Second safeguard...
    eta = std::min(etaMax, eta);

    // Not too small...
    eta = std::max(etaMin, eta);

    if (DEBUG_NONLINEAR)
    {
      Xyce::dout() << "\tnlResNorm: " << nlResNorm
                   << " linRes: " << linRes.getImmutableValue<double>()
                   << " nlResNormOld: " << nlResNormOld
                   << " calculated eta: " << fabs(nlResNorm - linRes.getImmutableValue<double>())/nlResNormOld << std::endl;
    }
  }

  else
    eta = etaMax;

  etaOld       = eta;
  nlResNormOld = nlResNorm;

  if (VERBOSE_NONLINEAR)
  {
    Xyce::lout() << "\t\teta:\t" <<  eta <<  "\n" << std::endl;
  }

  // Set the parameter
  nlParams.setForcingTerm(eta);

}


//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
bool DampedNewton::isFirstContinuationParam() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
int DampedNewton::getContinuationStep() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
int DampedNewton::getParameterNumber() const
{
  return 0;
}

} // namespace Nonlinear
} // namespace Xyce

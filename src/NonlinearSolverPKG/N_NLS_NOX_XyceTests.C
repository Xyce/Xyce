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
// Purpose        : Status test.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 04/15/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include "N_NLS_NOX_Group.h"
#include "N_NLS_NOX_XyceTests.h"
#include "N_NLS_NOX_Vector.h"
#include "N_LAS_Vector.h"
#include "N_LAS_Solver.h"
#include "N_LAS_SystemHelpers.h"
#include "N_LOA_NonlinearEquationLoader.h"
#include "NOX.H"
#include "NOX_Solver_LineSearchBased.H"

#include <N_UTL_MachDepParams.h>

// ----------   Namespaces   ----------

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

// ----------   Code   ----------

//-----------------------------------------------------------------------------
// Function      : XyceTests::XyceTests
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
XyceTests::XyceTests(
  Parallel::Machine                     comm,
  bool                                  isTransient,
  double                                normF,
  double                                machPrec,
  TimeIntg::DataStore*                  data_store,
  double                                epsilon_a, 
  double                                epsilon_r, 
  double                                tol,
  int                                   maxIters,
  double                                convRate,
  double                                relConvRate,
  double                                maxConvRate,
  double                                stagnationTol,
  int                                   maxBadSteps,
  int                                   checkDeviceConvergence,
  double                                smallUpdateTol,
  Loader::NonlinearEquationLoader *     loader, 
  Linear::Solver *                      lsolver,
  bool                                  maskingFlag,
  Linear::Vector *                      maskVectorPtr)
  : comm_(comm),
    status_(NOX::StatusTest::Unconverged),
    returnTest_(0),
    isTransient_(isTransient),
    niters_(-1),
    maxNormFindex_(-1),
    maxNormF_(0.0),
    normF_curr_(0.0),
    requestedMaxNormF_(normF),
    requestedMachPrecTol_(machPrec),
    pWeightsVectorPtr_(0),
    dsPtr_(data_store),
    weightsVectorPtr_(0),
    updateVectorPtr_(0),
    epsilon_a_(epsilon_a),
    epsilon_r_(epsilon_r),
    tol_(tol),
    weightedUpdate_(0.0),
    maxIters_(maxIters),
    requestedConvRate_(convRate),
    currentConvRate_(1.0),
    requestedRelativeConvRate_(relConvRate),
    currentRelativeConvRate_(1.0),
    normResidualInit_(1.0),
    smallUpdateTol_(smallUpdateTol),
    maxConvRate_(maxConvRate),
    lastIteration_(-1),
    badStepCount_(0),
    maxBadSteps_(maxBadSteps),
    minConvRate_(1.0),
    stagnationTol_(stagnationTol),
    xyceReturnCode_(0),
    checkDeviceConvergence_(checkDeviceConvergence),
    loaderPtr_(loader),
    lasSolverPtr_(lsolver),
    maskingFlag_( maskingFlag ),
    weightMaskVectorPtr_( maskVectorPtr),
    allDevicesConverged_(false),
    innerDevicesConverged_(false)
{
}

//-----------------------------------------------------------------------------
// Function      : XyceTests::XyceTests
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
XyceTests::~XyceTests()
{
  if ( pWeightsVectorPtr_ != 0 )
  {
    delete pWeightsVectorPtr_;
  }

  if( weightsVectorPtr_ != 0 )
  {
    delete weightsVectorPtr_;
    delete updateVectorPtr_;
  }
}


//-----------------------------------------------------------------------------
// Function      : XyceTests::checkStatus
// Purpose       : main testing function for this class
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
NOX::StatusTest::StatusType
XyceTests::checkStatus(
  const NOX::Solver::Generic&   problem,
  NOX::StatusTest::CheckType    checkType)
{
  status_ = NOX::StatusTest::Unconverged;
  xyceReturnCode_ = 0;
  niters_ = problem.getNumIterations();

  returnTest_ = 0;

  const double machineEpsilon = Util::MachineDependentParams::MachineEpsilon();

  // Get the current and previous solutions
  const Linear::Vector& x = (dynamic_cast<const Vector&>
    (problem.getSolutionGroup().getX())).getNativeVectorRef();
  const Linear::Vector& oldX = (dynamic_cast<const Vector&>
    (problem.getPreviousSolutionGroup().getX())).getNativeVectorRef();

  // Test #9 (-9):  if linear solver failed, reject.
  // Need to check if the reason for the failure is NaNs in the linear system.
  const N_NLS_NOX::Group & grp = 
      (dynamic_cast<const N_NLS_NOX::Group &>(problem.getSolutionGroup()));

  const Linear::Vector& F = (dynamic_cast<const Vector&>
    (problem.getSolutionGroup().getF())).getNativeVectorRef();

  bool linearSolverStatus = grp.linearSolverStatus();
  if (!linearSolverStatus)
  {
    status_ = NOX::StatusTest::Failed;

    // Check if NaNs were detected
    std::vector< std::pair<int, int> >& nanEntries = lasSolverPtr_->getNaNEntries();
    if ( VERBOSE_NONLINEAR )
    {
      if ( nanEntries.size() )
      {
        Xyce::lout() << "NaN/Inf check found " << nanEntries.size() << " entries in the residual or Jacobian!" << std::endl;
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
    {
      xyceReturnCode_ = retCodes_.nanFail; // default: -6
      returnTest_ = 0;
    }
    else
    {
      xyceReturnCode_ = retCodes_.linearSolverFailed; // default: -9
      returnTest_ = 9;
    }
    return status_;
  }

  // Test0 - NaN/Inf checker
  NOX::StatusTest::StatusType check = finiteTest_.checkStatus(problem, checkType);
  if (check == NOX::StatusTest::Failed) 
  {
    status_ = check;
    returnTest_ = 0;
    xyceReturnCode_ = retCodes_.nanFail; // default: -6

    if ( VERBOSE_NONLINEAR )
    {
      std::vector<int> nanEntries;
      Linear::checkVectorForNaNs( F, nanEntries );
      if ( nanEntries.size() )
      {
        Xyce::lout() << "NaN/Inf check found " << nanEntries.size() << " entries in the residual vector!" << std::endl;
        for (int i=0; i<nanEntries.size(); ++i)
        {
          Xyce::lout() << "Residual entry: [" << nanEntries[i] << "]" << std::endl;
        }
      }
    }

    return status_;
  }

  // Test 8 - Devices need to satisfy their own convergence criteria
  if (checkDeviceConvergence_) 
  {
    allDevicesConverged_ = loaderPtr_->allDevicesConverged(comm_);
    if (!allDevicesConverged_)
    {
      if (niters_ < maxIters_) 
      {
        status_ = NOX::StatusTest::Unconverged;
        returnTest_ = 8;
        xyceReturnCode_ = 0;
        return status_;
      }
      else
      {
        status_ = NOX::StatusTest::Failed;
        returnTest_ = 8;
        xyceReturnCode_ = retCodes_.tooManySteps;
        return status_;
      }
    }
  }

  // This test is for 2-level solves ONLY.
  innerDevicesConverged_ = loaderPtr_->innerDevicesConverged(comm_);
  if (!innerDevicesConverged_ )
  {
    status_ = NOX::StatusTest::Unconverged;
    returnTest_ = 8;
    xyceReturnCode_ = 0;
    return status_;
  }

  // Compute a few norms needed by various tests.
  // Max F norm (maxNormF_).
  if (maskingFlag_)
  {
    if (pWeightsVectorPtr_ == 0)
    {
      pWeightsVectorPtr_ = x.cloneVector();
      pWeightsVectorPtr_->putScalar( 1.0 );

      for (int i = 0; i < x.localLength(); ++i ) 
      {
        if ( (*weightMaskVectorPtr_)[i] == 0.0 )
        {
          (*pWeightsVectorPtr_)[i] = Util::MachineDependentParams::MachineBig();
        }
      }
    }
    F.wMaxNorm( *pWeightsVectorPtr_, &maxNormF_, &maxNormFindex_ );
  }
  else
  {
    F.infNorm( &maxNormF_, &maxNormFindex_ );
  }

  // Compute 2-norms of residual
  double normF_prev = 0.0;

  if (maskingFlag_)
  {
    int length = F.globalLength();
    F.wRMSNorm(*pWeightsVectorPtr_, &normF_curr_);
    normF_curr_ *= std::sqrt( length );  // Undo scaling of RMS norm

    const Linear::Vector& prevF = (dynamic_cast<const Vector&>
    (problem.getPreviousSolutionGroup().getF())).getNativeVectorRef();
    prevF.wRMSNorm(*pWeightsVectorPtr_, &normF_prev);
    normF_prev *= std::sqrt( length );  // Undo scaling of RMS norm
  }
  else
  {
    normF_curr_ = (problem.getSolutionGroup().getNormF());
    normF_prev = (problem.getPreviousSolutionGroup().getNormF());
  }

  // Test 1 - Residual (2-)norm too small 
  // If the RHS norm is so small already, just say we're converged,
  // and move on.  If we don't do this, we may wind up dividing by a
  // zero later.  Only perform this test in transient, not valid for DC.
  if (isTransient_ && (normF_curr_ < machineEpsilon))
  { 
    currentConvRate_ = normF_curr_ / normF_prev;

    if (normF_prev < machineEpsilon)
    { 
      currentConvRate_ = 1.0;
    }
    else
    { 
      currentConvRate_ = 0.0;
    }
    if (normResidualInit_ < machineEpsilon)
    { 
      currentRelativeConvRate_ = 1.0;
    }
    else
    { 
      currentRelativeConvRate_ = 0.0;
    }
   
    status_ = NOX::StatusTest::Converged;
    returnTest_ = 1;
    xyceReturnCode_ = retCodes_.normTooSmall; // default: 1
    return status_;
  }

  // Test 2 - Normal convergence based on rhs residual (2a) and 
  // update norm (2b).

  // Allocate space if necessary
  if (weightsVectorPtr_ == 0) 
  {
    weightsVectorPtr_ = x.cloneCopyVector();
    // when we create weightsVectorPtr_ from the latest solution, there
    // is a chance that one of the values will be zero.  If this isn't 
    // the DC op step or the first iteration of a time step, then
    // we'll end up dividing by zero in the wMaxNorm function below.
    // So to be safe we'll just add epsilon_a_ on the when we create 
    // this vector.
    for (int i=0; i< x.localLength() ; ++i ) 
    {
      (*weightsVectorPtr_)[i] += epsilon_a_;
    }    
    updateVectorPtr_ = x.cloneCopyVector();
  }

  // Local references
  Linear::Vector& weights = *weightsVectorPtr_;

  // Compute local portion of weights vector
  // Weights are recomputed at each nonlinear iteration of a DC Op calc
  // but only at the beginning of a transient nonlinear solve.
  if ((!isTransient_) || (niters_ == 0)) 
  {
    int length = x.localLength();
    for (int i = 0; i < length; ++i ) 
    {
      //update[i] = x[i] - oldX[i];
      weights[i] =
        epsilon_r_ * std::max(fabs(x[i]), fabs((*dsPtr_->currSolutionPtr)[i])) + epsilon_a_;
      if (maskingFlag_ && ((*weightMaskVectorPtr_)[i] == 0.0) )
      {
        weights[i] = Util::MachineDependentParams::MachineBig();
      }
    }
  }

  if (niters_ < 1) 
  {
    weightedUpdate_ = 1.0;
  }
  else 
  {
    // Next compute the update
    updateVectorPtr_->update(1.0, x, -1.0, oldX, 0.0);

    // Compute final result
    updateVectorPtr_->wMaxNorm(weights,&weightedUpdate_);

    // RPP: If a line search is being used, we must account for any 
    // damping of the step length.  Otherwise delta X could be small due 
    // the line search and not due to being close to a solution.
    const NOX::Solver::LineSearchBased* test = 0;
    test = dynamic_cast<const NOX::Solver::LineSearchBased*>(&problem);
    if (test != 0) 
    {
      weightedUpdate_ = weightedUpdate_/(test->getStepSize());
    }

    if ((weightedUpdate_ < tol_) &&(maxNormF_ < requestedMaxNormF_)) 
    {
      status_ = NOX::StatusTest::Converged;
      returnTest_ = 2;
      xyceReturnCode_ = retCodes_.normalConvergence; // default: 2
      return status_;
    }
  }
 
  // Test 3 - Near Convergence - Hit max iterations but residual 
  // and convergence rate indicate we may be near a converged solution.  
  // Therefore, let the time stepper decide whether or  not the step is ok.
  // Transient mode ONLY!
  // NOTE: Convergence rates are based on the 2-Norm, not max norm!
  if (niters_ > 0) 
  {
    // ||F(x_current)|| / ||F(x_previous)||
    currentConvRate_ = normF_curr_ / normF_prev;

    //  ||F(x)|| / ||F(x_init)||
    currentRelativeConvRate_ = normF_curr_ / normResidualInit_;
  }    
  else 
  {
    currentConvRate_ = 1.0;
    currentRelativeConvRate_ = 1.0;
  }

  if (isTransient_) 
  {
    if (niters_ == 0) 
    {
      if (maskingFlag_)
      {
        int length = F.globalLength();
        F.wRMSNorm(*pWeightsVectorPtr_, &normResidualInit_);
        normResidualInit_ *= std::sqrt( length );  // Undo scaling of RMS norm
      }
      else
        normResidualInit_ = problem.getSolutionGroup().getNormF();
    }
 
    // Test only if we hit the max number of iterations
    if (niters_ >= maxIters_) 
    {
      if ((currentConvRate_ <= requestedConvRate_) && 
          (currentRelativeConvRate_ <= requestedRelativeConvRate_)) 
      {
        status_ = NOX::StatusTest::Converged;
        xyceReturnCode_ = retCodes_.nearConvergence; // default: 3

        if (xyceReturnCode_ < 0)
          status_ = NOX::StatusTest::Failed;
        else
          status_ = NOX::StatusTest::Converged;
      }
      else
      {
        status_ = NOX::StatusTest::Failed;
        xyceReturnCode_ = retCodes_.tooManySteps; // default: -1
      }
      
      returnTest_ = 3;
      return status_;
    }
  } // end test 3

  // Test 4 - Update is too small
  if ((niters_ > 0) && (weightedUpdate_ < smallUpdateTol_) && (niters_ < maxIters_)) 
  {
    if (isTransient_) // Let the time integrator determine convergence (+4)
    {
      xyceReturnCode_ = retCodes_.smallUpdate; // default: 4
      status_ = NOX::StatusTest::Failed;
    }
    else
    {
      xyceReturnCode_ = 0; // neither pass nor fail.
      status_ = NOX::StatusTest::Unconverged;
    }

    returnTest_ = 4;

    return status_;
  }

  // Test 5 - Max nonlinear iterations (if transient, this will be checked 
  // in the NearConvergence test (#3)
  if ((!isTransient_) && (niters_ >= maxIters_)) 
  {
    status_ = NOX::StatusTest::Failed;
    returnTest_ = 5;
    xyceReturnCode_ = retCodes_.tooManySteps; // default: -1
    return status_;
  }

  // Test 6 - update is too big
  if (currentConvRate_ > maxConvRate_) 
  {
    status_ = NOX::StatusTest::Failed;
    returnTest_ = 6;
    xyceReturnCode_ = retCodes_.updateTooBig; // default: -2
    return status_;
  }

  // Test 7 - Stall in the convergence rate. Transient mode ONLY! 
  if (isTransient_) 
  {  
    // First time through we don't do anything but reset the counters
    if (niters_ == 0) 
    {
      badStepCount_ = 0;
      lastIteration_ = 0;
      //minConvRate = 1.0;  // Don't reset this.  Xyce solver never does.
    } 

    // Make sure we have not already counted the last nonlinear iteration.
    // This protects against multiple calls to checkStatus() in between 
    // nonlinear iterations.
    bool isCounted = false;
    if (niters_ == lastIteration_) 
    {
      isCounted = true;
    }
    else
    {
      lastIteration_ = niters_;
    }
    
    // Set counter appropriately
    if (!isCounted) 
    {
      if (fabs(currentConvRate_ - 1.0) <= stagnationTol_) 
      {
        if ((badStepCount_ == 0) || (currentConvRate_ < minConvRate_)) 
        {
          minConvRate_ = currentConvRate_;
        }
        ++badStepCount_ ;
      }
      else
      {
        badStepCount_ = 0;
      }
    }

    if (badStepCount_ >= maxBadSteps_) 
    {
      if ((currentRelativeConvRate_ <= 0.9) && (minConvRate_ <= 1.0)) 
      {
        status_ = NOX::StatusTest::Converged;
        returnTest_ = 7;
        xyceReturnCode_ = retCodes_.nearConvergence;   // default: 3
                // note - I'm not sure if this is 
               // really a near convergece test - but 3 is the code for it...

        if (xyceReturnCode_ < 0)
          status_ = NOX::StatusTest::Failed;
        else
          status_ = NOX::StatusTest::Converged;
      }
      else 
      {
        status_ = NOX::StatusTest::Failed;
        returnTest_ = 7;
        xyceReturnCode_ = retCodes_.stalled; // default: -3
      }
    }
  }

  return status_;
}

//-----------------------------------------------------------------------------
// Function      : XyceTests::print
// Purpose       : verbose output
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
std::ostream& XyceTests::print(std::ostream& stream, int indent) const
{
  // precision
  int p = 5;

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << status_ << "by Test #" << returnTest_ << "\n";

  indent += 4;

  //for (int j = 0; j < indent; ++j )
  // stream << ' ';
  finiteTest_.print(stream, indent);

  if (checkDeviceConvergence_) {
    for (int j = 0; j < indent; ++j )
      stream << ' ';
    stream << "8. Devices are Converged: ";
    if (allDevicesConverged_)
      stream << "true" << "\n";
    else
      stream << "false" << "\n";
  }

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "1. Two-Norm F too small" << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Machine Precision: " << NOX::Utils::sciformat(normF_curr_, p)
	 << " < " << NOX::Utils::sciformat(requestedMachPrecTol_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "2. Normal Convergence" << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Inf-Norm F: " << NOX::Utils::sciformat(maxNormF_, p)
	 << " < " << NOX::Utils::sciformat(requestedMaxNormF_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Weighted Update: " << NOX::Utils::sciformat(weightedUpdate_, p)
	 << " < " << NOX::Utils::sciformat(tol_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "3. Near Convergence" << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Max Iters: " << niters_
	 << " < " << maxIters_ << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Convergence Rate: " 
	 << NOX::Utils::sciformat(currentConvRate_, p)
	 << " < " << NOX::Utils::sciformat(requestedConvRate_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Relative Convergence Rate: " 
	 << NOX::Utils::sciformat(currentRelativeConvRate_, p)
	 << " < " << NOX::Utils::sciformat(requestedRelativeConvRate_, p) 
	 << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "4. Small Weighted Update: " 
	 << NOX::Utils::sciformat(weightedUpdate_, p)
	 << " < " << NOX::Utils::sciformat(smallUpdateTol_, p) << "\n";

  if (!isTransient_) {

    for (int j = 0; j < indent; ++j )
      stream << ' ';
    stream << "5. Maximum Iterations: " 
	   << niters_
	   << " < " << maxIters_ << "\n";
    
    for (int j = 0; j < indent; ++j )
      stream << ' ';
    stream << "6. Large Conv Rate: " 
	   << NOX::Utils::sciformat(currentConvRate_, p)
	   << " < " << NOX::Utils::sciformat(maxConvRate_, p) << "\n";
  }

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "7. Stagnation " << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Bad Step Count: " 
	 << badStepCount_ << " < " << maxBadSteps_ << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "   Stagnation Tolerance: " 
	 << NOX::Utils::sciformat(fabs(currentConvRate_ - 1.0), p)
	 << " < " << NOX::Utils::sciformat(stagnationTol_, p) << "\n";

  for (int j = 0; j < indent; ++j )
    stream << ' ';
  stream << "9. Linear solver failed." << "\n";

  stream << std::endl;
  return stream;
}

//-----------------------------------------------------------------------------
// Function      : XyceTests::getXyceReturnCode
// Purpose       : verbose output
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
int XyceTests::getXyceReturnCode() const
{
  return xyceReturnCode_;
}

//-----------------------------------------------------------------------------
// Function      : XyceTests::getMaxNormF
// Purpose       : verbose output
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
double XyceTests::getMaxNormF() const
{
  return maxNormF_;
}

//-----------------------------------------------------------------------------
// Function      : XyceTests::getMaxNormFindex
// Purpose       : verbose output
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
int XyceTests::getMaxNormFindex() const
{
  return maxNormFindex_;
}

}}}

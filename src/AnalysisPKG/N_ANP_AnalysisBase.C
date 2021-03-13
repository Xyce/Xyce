//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose       : Base class analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_AnalysisManager.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ERH_Message.h>
#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h> 
#include <N_PDS_Comm.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Analysis {

StatCounts::StatCounts()
  : successfulStepsTaken_(0),
    successStepsThisParameter_(0),
    failedStepsAttempted_(0),
    jacobiansEvaluated_(0),
    iterationMatrixFactorizations_(0),
    linearSolves_(0),
    failedLinearSolves_(0),
    linearIters_(0),
    residualEvaluations_(0),
    nonlinearConvergenceFailures_(0),
    linearSolutionTime_(0.0),
    residualLoadTime_(0.0),
    jacobianLoadTime_(0.0)
{}

StatCounts &
StatCounts::operator+=(
  const StatCounts &   stats)
{
  successfulStepsTaken_ += stats.successfulStepsTaken_;
  successStepsThisParameter_ += stats.successStepsThisParameter_;
  failedStepsAttempted_ += stats.failedStepsAttempted_;
  jacobiansEvaluated_ += stats.jacobiansEvaluated_;
  iterationMatrixFactorizations_ += stats.iterationMatrixFactorizations_;
  linearSolves_ += stats.linearSolves_;
  failedLinearSolves_ += stats.failedLinearSolves_;
  linearIters_ += stats.linearIters_;
  residualEvaluations_ += stats.residualEvaluations_;
  nonlinearConvergenceFailures_ += stats.nonlinearConvergenceFailures_;
  linearSolutionTime_ += stats.linearSolutionTime_;
  residualLoadTime_ += stats.residualLoadTime_;
  jacobianLoadTime_ += stats.jacobianLoadTime_;

  return *this;
}


Util::JSON &operator<<(Util::JSON &json, const StatCounts &s)
{
  json << Util::JSON::open
       << Util::nameValuePair("successfulStepsTaken", s.successfulStepsTaken_) << Util::JSON::sep
       << Util::nameValuePair("failedStepsAttempted", s.failedStepsAttempted_) << Util::JSON::sep
       << Util::nameValuePair("jacobiansEvaluated", s.jacobiansEvaluated_) << Util::JSON::sep
       << Util::nameValuePair("iterationMatrixFactorizations", s.iterationMatrixFactorizations_) << Util::JSON::sep
       << Util::nameValuePair("linearSolves", s.linearSolves_) << Util::JSON::sep
       << Util::nameValuePair("failedLinearSolves", s.failedLinearSolves_) << Util::JSON::sep
       << Util::nameValuePair("linearIters", s.linearIters_) << Util::JSON::sep
       << Util::nameValuePair("residualEvaluations", s.residualEvaluations_) << Util::JSON::sep
       << Util::nameValuePair("nonlinearConvergenceFailures", s.nonlinearConvergenceFailures_) << Util::JSON::sep
       << Util::nameValuePair("residualLoadTime", s.residualLoadTime_) << Util::JSON::sep
       << Util::nameValuePair("jacobianLoadTime", s.jacobianLoadTime_) << Util::JSON::sep
       << Util::nameValuePair("linearSolutionTime", s.linearSolutionTime_)
       << Util::JSON::close;

  return json;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::AnalysisBase
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, SNL.
// Creation Date : 01/29/08
//-----------------------------------------------------------------------------
AnalysisBase::AnalysisBase(AnalysisManager &analysis_manager, const char *name )
  : name_(name),
    NOOP_(false), 
    inputOPFlag_(false),
    doubleDCOPEnabled_(false),
    doubleDCOPStep_(0),
    firstDCOPStep_(-1),
    lastDCOPStep_(-1),
    beginningIntegration(true),
    baseIntegrationMethod_(TimeIntg::NoTimeIntegration::type),
    stepNumber(0),
    tranStepNumber(0),
    dataSpecification_(false),
    saveStatCountsVector_(),
    stats_()
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::~AnalysisBase()
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Todd S. Coffey, SNL.
// Creation Date : 01/29/08
//-----------------------------------------------------------------------------
AnalysisBase::~AnalysisBase()
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::run
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool AnalysisBase::run()
{
  Stats::StatTop _analysis(name_);
  Stats::TimeBlock _analysisTimer(_analysis);

  return doRun();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::init
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool AnalysisBase::init()
{
  Stats::StatTop _analysis(name_);
  Stats::TimeBlock _analysisTimer(_analysis);

  return doInit();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::processSuccessfulStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool AnalysisBase::processSuccessfulStep()
{
  return doProcessSuccessfulStep();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::processFailedStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool AnalysisBase::processFailedStep()
{
  return doProcessFailedStep();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::finish
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool AnalysisBase::finish()
{
  return doFinish();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::handlePredictor
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool AnalysisBase::handlePredictor()
{
  return doHandlePredictor();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void AnalysisBase::registerLinearSystem(Linear::System * linear_system_ptr)
{
  static std::string tmp = "This analysis type doesn't support embedded sampling\n";
  Report::UserFatal0() << tmp;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::resetForStepAnalysis()
// Purpose       : When doing a .STEP sweep, some data must be reset to its
//                 initial state.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool AnalysisBase::resetForStepAnalysis()
{
  stats_.successStepsThisParameter_ = 0;
  stepNumber = 0;
  beginningIntegration = true;

  return true;
}



//-----------------------------------------------------------------------------
// Function      : AnalysisBase::resetAll
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/16/10
//-----------------------------------------------------------------------------
void AnalysisBase::resetAll()
{
  stepNumber = 0;
  tranStepNumber = 0;

  stats_ = StatCounts();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::saveLoopInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 12/16/10
//-----------------------------------------------------------------------------
int AnalysisBase::saveLoopInfo ()
{
  if (saveStatCountsVector_.empty()) // push back empty stats as sentinel
    saveStatCountsVector_.push_back(StatCounts());

  saveStatCountsVector_.push_back(stats_);

  return saveStatCountsVector_.size() - 1;
}

StatCounts
operator-(
  const StatCounts &   s0,
  const StatCounts &   s1)
{
  StatCounts s;

  s.successfulStepsTaken_ = s0.successfulStepsTaken_ - s1.successfulStepsTaken_;
  s.failedStepsAttempted_ = s0.failedStepsAttempted_ - s1.failedStepsAttempted_;
  s.jacobiansEvaluated_ = s0.jacobiansEvaluated_ - s1.jacobiansEvaluated_;
  s.iterationMatrixFactorizations_ = s0.iterationMatrixFactorizations_ - s1.iterationMatrixFactorizations_;
  s.linearSolves_ = s0.linearSolves_ - s1.linearSolves_;
  s.failedLinearSolves_ = s0.failedLinearSolves_ - s1.failedLinearSolves_;
  s.linearIters_ = s0.linearIters_ - s1.linearIters_;
  s.residualEvaluations_ = s0.residualEvaluations_ - s1.residualEvaluations_;
  s.nonlinearConvergenceFailures_ = s0.nonlinearConvergenceFailures_ - s1.nonlinearConvergenceFailures_;
  s.residualLoadTime_ = s0.residualLoadTime_ - s1.residualLoadTime_;
  s.jacobianLoadTime_ = s0.jacobianLoadTime_ - s1.jacobianLoadTime_;
  s.linearSolutionTime_ = s0.linearSolutionTime_ - s1.linearSolutionTime_;

  return s;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisBase::printLoopInfo
// Purpose       : Prints out time loop information.
// Special Notes : Prints stats from save point start to save point finish.
//                 Special case 0,0 is entire run to this point
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/26/00
//-----------------------------------------------------------------------------
bool AnalysisBase::printLoopInfo(int t1, int t2)
{
  bool bsuccess = true;

  if (t1 == 0 && t2 == 0)
  {
    t2 = saveLoopInfo();
  }

  lout() << "\tNumber Successful Steps Taken:\t\t" << saveStatCountsVector_[t2].successfulStepsTaken_ - saveStatCountsVector_[t1].successfulStepsTaken_ << std::endl
         << "\tNumber Failed Steps Attempted:\t\t" << saveStatCountsVector_[t2].failedStepsAttempted_ - saveStatCountsVector_[t1].failedStepsAttempted_ << std::endl
         << "\tNumber Jacobians Evaluated:\t\t" << saveStatCountsVector_[t2].jacobiansEvaluated_ - saveStatCountsVector_[t1].jacobiansEvaluated_ << std::endl;

  // This statistic is only collected when two-level Newton is being used.
  if (saveStatCountsVector_[t2].iterationMatrixFactorizations_ > saveStatCountsVector_[t1].iterationMatrixFactorizations_)
  {
    lout() << "\tNumber Iteration Matrix Factorizations:\t" << saveStatCountsVector_[t2].iterationMatrixFactorizations_ - saveStatCountsVector_[t1].iterationMatrixFactorizations_ << std::endl;
  }

  lout() << "\tNumber Linear Solves:\t\t\t" << saveStatCountsVector_[t2].linearSolves_ - saveStatCountsVector_[t1].linearSolves_ << std::endl
         << "\tNumber Failed Linear Solves:\t\t" << saveStatCountsVector_[t2].failedLinearSolves_ - saveStatCountsVector_[t1].failedLinearSolves_ << std::endl;

  // This statistic is only collected when iterative linear solvers are being used.
  if (saveStatCountsVector_[t2].linearIters_ > saveStatCountsVector_[t1].linearIters_)
  {
    lout() << "\tNumber Linear Solver Iterations:\t" << saveStatCountsVector_[t2].linearIters_ - saveStatCountsVector_[t1].linearIters_ << std::endl;
  }
  lout() << "\tNumber Residual Evaluations:\t\t" << saveStatCountsVector_[t2].residualEvaluations_ - saveStatCountsVector_[t1].residualEvaluations_ << std::endl
         << "\tNumber Nonlinear Convergence Failures:\t" << saveStatCountsVector_[t2].nonlinearConvergenceFailures_ - saveStatCountsVector_[t1].nonlinearConvergenceFailures_ << std::endl
         << "\tTotal Residual Load Time:\t\t" << saveStatCountsVector_[t2].residualLoadTime_ - saveStatCountsVector_[t1].residualLoadTime_ << " seconds" << std::endl
         << "\tTotal Jacobian Load Time:\t\t" << saveStatCountsVector_[t2].jacobianLoadTime_ - saveStatCountsVector_[t1].jacobianLoadTime_ << " seconds" << std::endl
         << "\tTotal Linear Solution Time:\t\t" << saveStatCountsVector_[t2].linearSolutionTime_ - saveStatCountsVector_[t1].linearSolutionTime_ << " seconds" << std::endl << std::endl;

  return bsuccess;
}

void gatherStepStatistics(StatCounts &stats, Nonlinear::NonLinearSolver &nonlinear_solver, int newton_convergence_status) 
{
  if (newton_convergence_status <= 0) {
    ++stats.nonlinearConvergenceFailures_;
  }

  stats.jacobiansEvaluated_      += nonlinear_solver.getNumJacobianLoads();
  stats.linearSolves_            += nonlinear_solver.getNumLinearSolves();
  stats.failedLinearSolves_      += nonlinear_solver.getNumFailedLinearSolves();
  stats.linearIters_             += nonlinear_solver.getTotalNumLinearIters();
  stats.residualEvaluations_     += nonlinear_solver.getNumResidualLoads();
  stats.iterationMatrixFactorizations_ += nonlinear_solver.getNumJacobianFactorizations();
  stats.linearSolutionTime_            += nonlinear_solver.getTotalLinearSolveTime();
  stats.residualLoadTime_              += nonlinear_solver.getTotalResidualLoadTime();
  stats.jacobianLoadTime_              += nonlinear_solver.getTotalJacobianLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : setTimeIntegratorOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool
AnalysisBase::setDCOPOption(
  const Util::Param &           param)
{
  if (param.uTag() == "DOUBLEDCOP") {
    const std::string &double_dcop_options = param.stringValue();
    if (equal_nocase(double_dcop_options, "NL_POISSON,DRIFT_DIFFUSION")) {
      firstDCOPStep_ = 0;
      lastDCOPStep_ = 1;
      doubleDCOPStep_ = firstDCOPStep_;
      return true;
    }
    else if (equal_nocase(double_dcop_options, "NL_POISSON")) {
      firstDCOPStep_ = 0;
      lastDCOPStep_ = 0;
      doubleDCOPStep_ = firstDCOPStep_;
      return true;
    }
    else if (equal_nocase(double_dcop_options, "DRIFT_DIFFUSION")) {
      firstDCOPStep_ = 1;
      lastDCOPStep_ = 1;
      doubleDCOPStep_ = firstDCOPStep_;
      return true;
    }
    else
      Report::UserError0() << "Unknown DOUBLEDCOP option " << double_dcop_options;
  }

  return false;
}

} // namespace Analysis
} // namespace Xyce

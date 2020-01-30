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

//-----------------------------------------------------------------------------
//
// Purpose        : Base class for Analysis types
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AnalysisBase_h
#define Xyce_N_ANP_AnalysisBase_h

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_JSON.h>
#include <N_UTL_Stats.h>

namespace Xyce {
namespace Analysis {

struct StatCounts
{
  StatCounts();
  StatCounts &operator+=(const StatCounts &stats);

  unsigned int        successfulStepsTaken_;        ///< Number of consecutive successful time-integration steps.
  unsigned int        successStepsThisParameter_;
  unsigned int        failedStepsAttempted_;        ///< Total number of failed time-integration steps.
  unsigned int        jacobiansEvaluated_;
  unsigned int        iterationMatrixFactorizations_;
  unsigned int        linearSolves_;
  unsigned int        failedLinearSolves_;
  unsigned int        linearIters_;
  unsigned int        residualEvaluations_;
  unsigned int        nonlinearConvergenceFailures_;
  double              linearSolutionTime_;
  double              residualLoadTime_;
  double              jacobianLoadTime_;
};

Util::JSON &operator<<(Util::JSON &json, const StatCounts &s);

class ProcessorBase
{
public:
  ProcessorBase()
  {}

  virtual ~ProcessorBase()
  {}
};


//-------------------------------------------------------------------------
// Class         : AnalysisBase
// Purpose       : Base class for common analysis functions
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class AnalysisBase : public ProcessorBase
{
public:
  AnalysisBase(AnalysisManager &analysis_manager, const char *name);
  virtual ~AnalysisBase();

private:
  AnalysisBase(const AnalysisBase &);                   ///< No copying
  AnalysisBase &operator=(const AnalysisBase &);        ///< No assignment

public:
  virtual void setTranStepNumber(int step)
  {
    tranStepNumber = step;
  }

  virtual int getTranStepNumber()
  {
    return tranStepNumber;
  }

  virtual int getStepNumber()
  {
    return stepNumber;
  }

  void setStepNumber(int step)
  {
    stepNumber=step;
  }

  virtual double getCurrentFreq()
  {
    return 0.0;
  }

  virtual bool getSparcalc()
  {
    return 0;
  }

  virtual void setRFParamsRequested(const std::string & type)
  {
    // This is a no op for every analysis mode, except for AC
  }

  virtual const TimeIntg::TIAParams &getTIAParams() const = 0;
  virtual TimeIntg::TIAParams &getTIAParams() = 0;

  virtual bool outputFailureStats(std::ostream &os) {return true;}

  bool run();
  bool init();
  bool processSuccessfulStep();
  bool processFailedStep();
  bool finish();
  bool handlePredictor();

  virtual void registerLinearSystem(Linear::System * linear_system_ptr);

  virtual void stepCallBack() {};
  virtual void registerParentAnalysis(AnalysisBase * parentPtr) {};

protected:
  virtual bool doRun() = 0;
  virtual bool doInit() = 0;
  virtual bool doProcessSuccessfulStep() = 0;
  virtual bool doProcessFailedStep() = 0;
  virtual bool doFinish() = 0;
  virtual bool doHandlePredictor() = 0;

public:
  //-----------------------------------------------------------------------------
  // Function      : AnalysisBase::printStepHeader()
  // Purpose       : Prints out step information.
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 6/26/00
  //-----------------------------------------------------------------------------
  virtual void printStepHeader(std::ostream &os)
  {}

  //-----------------------------------------------------------------------------
  // Function      : AnalysisBase::printProgress()
  // Purpose       : Outputs run completion percentage and estimated
  //                 time-to-completion.
  // Special Notes :
  // Scope         : public
  // Creator       : Scott A. Hutchinson, SNL, Computational Sciences
  // Creation Date : 06/07/2002
  //-----------------------------------------------------------------------------
  virtual void printProgress(std::ostream &os)
  {}

  // Two Level specific
  virtual bool twoLevelStep() { return true; }

  virtual bool isAnalysis(int analysis_type) const { return false; }

  // functions related to PDE problems requiring two DCOP solves
  bool getDoubleDCOPEnabled() const
  {
    return doubleDCOPEnabled_;
  }

  void setDoubleDCOPEnabled(bool enable)
  {
    doubleDCOPEnabled_ = enable;
  }

  virtual bool getDCOPFlag() const = 0;

  virtual int getDoubleDCOPStep() const
  {
    return doubleDCOPStep_;
  }

  //-----------------------------------------------------------------------------
  // Function      : AnalysisBase::firstDoubleDCOPStep
  // Purpose       : If the current step is the first step of
  //                  a "doubleDCOP", then return "true".
  //
  //  Explanation:
  //
  //  If there are PDE semiconductor devices as part of this problem,
  //  there may need to be a "double-pass"
  //
  //  first pass  = nonlinear poisson solution
  //  second pass = drift diffusion solution
  //
  // Special Notes : Only PDE problems can ever return true.
  // Scope         :
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 11/21/04
  //-----------------------------------------------------------------------------
  bool firstDoubleDCOPStep() const
  {
    return getDoubleDCOPEnabled() && getDoubleDCOPStep() != lastDCOPStep_;
  }

  void nextDCOPStep()
  {
    doubleDCOPStep_ = lastDCOPStep_;
  }

  bool setDCOPOption(const Util::Param &param);

  bool getNOOP() const
  {
    return NOOP_;
  }

  void setNOOP(bool noop)
  {
    NOOP_ = noop;
  }

  virtual bool printLoopInfo(int start, int finish);

  virtual void setBeginningIntegrationFlag(bool bif)
  {
    beginningIntegration = bif;
  }

  virtual bool getBeginningIntegrationFlag() const
  {
    return beginningIntegration;
  }

  virtual void setIntegrationMethod(int im)
  {
    baseIntegrationMethod_ = im;
  }

  virtual int getIntegrationMethod() const
  {
    return baseIntegrationMethod_;
  }

  void setInputOPFlag(bool initial_conditions_loaded)
  {
    inputOPFlag_ = initial_conditions_loaded;
  }
  
  bool getInputOPFlag() const
  {
    return inputOPFlag_;
  }

  bool getDataSpecification() const
  {
    return dataSpecification_;
  }

  bool resetForStepAnalysis();
  void resetAll();
  int saveLoopInfo ();

  // step statistic functions
  double getTotalLinearSolutionTime() const;
  double getTotalResidualLoadTime() const;
  double getTotalJacobianLoadTime() const;

  const StatCounts &getStatCounts(int index = -1) const
  {
    if (index == -1)
      return saveStatCountsVector_.back();
    else
      return saveStatCountsVector_[index];
  }

  const char *getName() const
  {
    return name_;
  }

private:
  const char *          name_;
  bool                  NOOP_;                          ///< if true, disable DCOP
  bool                  inputOPFlag_;                   ///< true if starting from an initial condition.
  bool                  doubleDCOPEnabled_;             ///< true if doing a double-DCOP is possible.
  int                   doubleDCOPStep_;                ///< current step in the DCOP loop.
  int                   firstDCOPStep_;
  int                   lastDCOPStep_;

protected:
  bool                  beginningIntegration;
  unsigned int          baseIntegrationMethod_;         ///< Current time-integration method flag.

    // NOTE: For now tranStepNumber is the same as stepNumber, but later I will
    //       change it.  stepNumber will later include both dcop and tran steps.
    //       I haven't changed it yet b/c I need to check what devices call
    //       getStepNumber, and what they expect to get.
  unsigned int          stepNumber;                     ///< Time-integration step number counter.

protected:
  unsigned int          tranStepNumber;
  bool                  dataSpecification_;

  std::vector<StatCounts>     saveStatCountsVector_;

public:
  StatCounts                    stats_;
};

StatCounts operator-(const StatCounts &s0, const StatCounts &s1);

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalLinearSolutionTime() const
{
  return stats_.linearSolutionTime_;
}

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalResidualLoadTime() const
{
  return stats_.residualLoadTime_;
}

//-----------------------------------------------------------------------------
inline double AnalysisBase::getTotalJacobianLoadTime() const
{
  return stats_.jacobianLoadTime_;
}

void gatherStepStatistics(StatCounts &stats, Nonlinear::NonLinearSolver &nonlinear_solver, int newton_convergence_status);

} // namespace Analysis
} // namespace Xyce


#endif // Xyce_N_ANP_AnalysisBase_h

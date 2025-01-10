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


//-----------------------------------------------------------------------------
//
// Purpose       : This file defines the class for the time integration
//                 stepsize control algorithms.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_STEP_ERROR_CONTROL_H
#define Xyce_N_TIA_STEP_ERROR_CONTROL_H

// ---------- Standard Declarations ----------
#include <iosfwd>
#include <set>
#include <vector>
#include <N_UTL_Math.h>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>

#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace TimeIntg {

//-----------------------------------------------------------------------------
// Class         : StepErrorControl
// Purpose       : This is the time integration step & error control class.
//
// Special Notes :  ERK. 6/19/2010  This class was originally designed to
//                  handle much of the time step selection logic.  This
//                  design intent was associated with the original
//                  "old-DAE(ODE)" version of the integrator, which
//                  is no longer in the code.  However numerous artifacts
//                  of this design are still here.  For example, the step size
//                  variables (such as currentTimeStep, currentTime)
//                  are still here.  So at this point, this is mostly a data
//                  storage class for step variables that have a global
//                  scope.
//
//                  The implementation of the "new-DAE" form
//                  moved away from this design idea, and placed most of the
//                  step-size selection directly inside of specific algorithms.
//                  So, OneStep and Gear handle stepsize selection specific to it.
//
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class StepErrorControl
{
  friend class Gear12;
  friend class OneStep;

  public:
  typedef std::vector<Util::BreakPoint> BreakPointVector;

  StepErrorControl(
      const std::string &         netlist_filename,
      Analysis::AnalysisManager & analysis_manager,
      WorkingIntegrationMethod &  working_integration_method,
      const TIAParams &           tia_params);

  virtual ~StepErrorControl();

  private:
  StepErrorControl(const StepErrorControl &);
  StepErrorControl &operator=(const StepErrorControl &);

  public:
  const Util::BreakPointLess &getBreakPointLess() const {
    return breakPointLess_;
  }

  // The "stop time" is either the next discontinuity point, or the final time,
  // whichever comes first.
  void updateStopTime(
      Parallel::Machine   comm, 
      bool                breakpoints_enabled,
      double              initial_time,
      bool                min_time_steps_breakpoint_given,
      double              min_time_steps_breakpoint);

  // Requests dynamic time step information from the loader.
  bool updateMaxTimeStep(Parallel::Machine comm, Loader::Loader &loader, const TIAParams &tia_params, double suggestedMaxTimeStep=0.0);

  void updateTwoLevelTimeInfo(
      Parallel::Machine   comm,
      double              nextTimeStep,
      double              nextTime,
      int                 currentOrder,
      bool                breakpoints_enabled,
      double              initial_time,
      bool                min_time_steps_breakpoint_given,
      double              min_time_steps_breakpoint);

  // Print out the time-step information.
  void outputTimeInfo(std::ostream &os);

  int getNumberOfSteps() const { return numberOfSteps_; }

  // Requests dynamic breakpoint information from the loader.  Adds, subtracts
  // from the breakpoints array.
  bool updateBreakPoints(const Loader::Loader &loader, double initial_time);

  // add individual breakpoints
  void setBreakPoint(const Util::BreakPoint &breakpoint, double initial_time);
  void setBreakPoint(double bp);
  void doubleCheckEndBreakPoint();

  // add a set of breakpoints
  void setBreakPoints(const std::vector<double>& times);

  // signal that pause breakpoint has been reached
  void simulationPaused(double initial_time);

  double getBreakPointEqualTolerance() const { return breakPointEqual_.tolerance_;}

  // is the current time a pause time?
  bool isPauseTime();

  // Sets the minimum time step based on machine precision.
  bool updateMinTimeStep();

  bool isFinished()
  {
    double delta = currentTime - finalTime;
    delta = std::fabs(delta);
    return (delta < 1.0e-10 * (finalTime - initialTime));
  }

  void evaluateStepError(const Loader::Loader &loader, const TIAParams &tia_params);

  // Gets the size of the restart data.
  int getRestartDataSize( bool pack );

  // Output restart data.
  bool dumpRestartData(char * buf, int bsize, int & pos, Parallel::Communicator * comm, bool pack );

  // Load restart data.
  bool restoreRestartData(char * buf, int bsize, int & pos, Parallel::Communicator * comm, bool pack, double &initial_time);

  // called if you want to start the time integration over from scratch.
  bool resetAll(const TIAParams &tia_params);

  // Method which sets the TIA parameters object pointer.
  bool setFromTIAParams(const TIAParams &tia_params);

  void printBreakPoints (std::ostream & os) const;

  double getEstOverTol() const;

  void setTimeStep(double newTimeStep);

  int getCurrentOrder() { return currentOrder_; }

  private:
  bool initializeBreakPoints(double start_time, double initial_time, double final_time);

  void updatePauseTime(Util::BreakPoint breakpoint, double initial_time);

  void integrationStepReport(std::ostream &os, bool step_attempt_status, bool sAStatus, bool testedError, const TIAParams &tia_params);

  void terseIntegrationStepReport(std::ostream &os, bool step_attempt_status, bool sAStatus, bool testedError, const TIAParams &tia_params);


  // member data:
  private:
  Analysis::AnalysisManager &   analysisManager_;       /// Analysis manager
  WorkingIntegrationMethod &    wimPtr_;                /// Working integration method.
  const std::string             netlistFilename_;

  public:
  // StepSize Variables
  double startingTimeStep;
  double currentTimeStep;
  double lastAttemptedTimeStep;
  double lastTimeStep;
  double oldeTimeStep;
  double minTimeStep;
  double maxTimeStep;
  double maxTimeStepUser;
  double maxTimeStepBP;

  double savedTimeStep;

  double lastTime;
  double currentTime;
  double nextTime;
  double stopTime;
  double initialTime;
  double finalTime;

  double currentTimeStepRatio;
  double currentTimeStepSum;

  double lastTimeStepRatio;
  double lastTimeStepSum;

  int newtonConvergenceStatus;

  int         nIterations;                    ///< Number of newton iterations

  int numberSuccessiveFailures;
  bool stepAttemptStatus;

  // Needed for 2-level:
  bool previousCallStepSuccessful;
  double estOverTol_;
  bool TimeStepLimitedbyBP;

  double pauseTime;                             ///< Time step value at which to "pause" the simulation.
  bool pauseSetAtZero;                          ///< Flag used to indicate that a pause was specifically set at time zero and thus should not be ignored.

  private:
  double minStepPrecisionFac_;
  double newtonStepReduction_;
  double restartTimeStepScale_;
  double tolAimFac_;

  Util::BreakPointLess breakPointLess_;
  Util::BreakPointEqual breakPointEqual_;

  BreakPointVector breakPoints_;
  BreakPointVector::iterator currentPauseBP;

  private :
  // 03/08/04 tscoffe:  Local data for BDF 1-5 method
  int currentOrder_;      // Current order of integration
  int minOrder_;          // minimum order = max(1,user option minord)
  int maxOrder_;          // maximum order = min(5,user option maxord)
  int usedOrder_;         // order used in current step (used after currentOrder is updated)
  double alphas_;         // $\alpha_s$ fixed-leading coefficient of this BDF method
  std::vector<double> alpha_;  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
  // note:   $h_n$ = current step size, n = current time step
  double alpha0_;         // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  double cj_;             // $-\alpha_s/h_n$ coefficient used in local error test
  double ck_;             // local error coefficient
  std::vector<double> sigma_;  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
  std::vector<double> gamma_;  // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  // calculate time derivative of history array for predictor
  std::vector<double> beta_;   // coefficients used to evaluate predictor from history array
  std::vector<double> psi_;    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to
  // compute $\beta_j(n)$
  int numberOfSteps_;     // number of total time integration steps taken
  int nef_;               // number of successive error failures (yes I know this is duplicated above)
  double  usedStep_;      // step-size used in current step
  int nscsco_;            // Number of Steps taken with Constant Step-size and Constant Order
  double Ek_,Ekm1_,Ekm2_,Ekp1_; // Estimated LTE at currentOrder, currentOrder-1, and currentOrder-2
  double Est_;            // Represents which of Ek,Ekm1,Ekm2 will be used for step-size selection
  double Tk_,Tkm1_,Tkm2_,Tkp1_; // Estimated $\|h_n^{k+1} y_n^{(k+1)}\|$,  basically, order * E
  int  newOrder_;         // order of next step
  bool initialPhase_;     // determines if we're in the initial phase of integration where we
  // double the step-size & increase the order.
  double h0_safety_;      // safety factor in picking initial step-size
  double h0_max_factor_;  // h0 <= h0_max_factor * length of integration time
  double h_phase0_incr_;  // amount to increase step-sizes during initial phase
  double h_max_inv_;      // inverse of maximum step size
  double Tkm1_Tk_safety_; // magic number for determining order reduction
  double Tkp1_Tk_safety_; // magic number for determining order reduction
  double r_factor_;       // basic reduction factor for rr
  double r_safety_;       // safety factor in computing rr
  double r_fudge_;        // fudge factor in computing rr
  double r_min_;          // minimum reduction factor
  double r_max_;          // maximum reduction factor
  double r_hincr_test_;   // threshold for increasing the step-size
  double r_hincr_;        // factor used for increasing the step-size
  int max_LET_fail_;      // max number of LTE test failures before quitting.

  int maxNumfail_;        // max number of error test failure
  bool reportedPauseBP;   // true if any devices have produced pause breakpoints.
};

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for step error control class.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const StepErrorControl & sec);

  } // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_STEP_ERROR_CONTROL_H

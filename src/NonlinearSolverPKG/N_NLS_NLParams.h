//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : Basic user-specifid parameters data structure.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NLParams_h
#define Xyce_N_NLS_NLParams_h

// ----------   Standard Includes   ----------
#include <N_UTL_Math.h>
#include <climits>

#include <N_IO_fwd.h>
#include <N_NLS_Manager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace Nonlinear {

// ---------- Enumerated Types ----------

// Support line-search methods.
enum LineSearchMethod {
  FULL,            // Newton's Method, undamped.
  DIVIDE,          // Reduce step-size by successively halving the previous
                   // step.
  BACKTRACK,       // Backtrack using interpolation (ala Dennis & Schnabel).
  BANK_ROSE,       // Bank and Rose algorithm - see method notes.
  DESCENT,         // More'-Thuente line search
  SIMPLE_BACKTRACK // Simple backtracking method
};

// Nonlinear solution "strategies".
enum NLStrategy {
  NEWTON,             // Pure Newton's method
  GRADIENT,           // Pure gradient method
  NEWTON_GRADIENT,    // Combined Newton/gradient method
  MOD_NEWTON,         // Pure modified Newton's method
  MOD_NEWTON_GRADIENT // Combined modified-Newton/gradient method
};

//-----------------------------------------------------------------------------
// Class         : NLParams
// Purpose       : Stores nonlinear solver user settings.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
class NLParams
{
public:
  static void populateMetadata(IO::PkgOptionsMgr &options_manager);

  //NLParams(AnalysisMode mode=DC_OP);
  NLParams(AnalysisMode mode, const IO::CmdParse & cp);
  NLParams(const NLParams & right);

  ~NLParams();

  NLParams & operator=(const NLParams & right);

  // ***** Accessor functions *****

  void setPrintParamsFlag();
  void clearPrintParamsFlag();
  bool getPrintParamsFlag() const;

  bool setOptions(const Util::OptionBlock & OB);
  bool setCmdLineOptions ();

  inline void setNLStrategy(NLStrategy strategy);
  inline void setNLStrategy(int strategy);
  inline void resetNLStrategy();
  inline NLStrategy getNLStrategy() const;

  inline void setSearchMethod(LineSearchMethod method);
  inline void setSearchMethod(int method);
  inline void resetSearchMethod();
  inline LineSearchMethod getSearchMethod() const;

  inline void   setDeltaXTol(double Tolerance);
  inline void   resetDeltaXTol();
  inline double getDeltaXTol() const;

  inline void   setSmallUpdateTol (double Tolerance);
  inline void   resetSmallUpdateTol ();
  inline double getSmallUpdateTol () const;

  inline void   setEnforceDeviceConvFlag (bool flag);
  inline void   resetEnforceDeviceConvFlag ();
  inline bool   getEnforceDeviceConvFlag () const;

  inline void   setRHSTol(double Tolerance);
  inline void   resetRHSTol();
  inline double getRHSTol() const;

  inline void   setAbsTol(double Tolerance);
  inline void   resetAbsTol();
  inline double getAbsTol() const;

  inline void   setRelTol(double Tolerance);
  inline void   resetRelTol();
  inline double getRelTol() const;

  inline void     setMaxNewtonStep(unsigned maxNewtonStep);
  inline void     resetMaxNewtonStep();
  inline unsigned getMaxNewtonStep() const;

  inline void     setMaxSearchStep(unsigned maxSearchStep);
  inline void     resetMaxSearchStep();
  inline unsigned getMaxSearchStep() const;

  inline void     setForcingFlag(bool flag);
  inline void     resetForcingFlag();
  inline bool     getForcingFlag() const;

  inline void     setForcingTerm(double value);
  inline void     resetForcingTerm();
  inline double   getForcingTerm() const;

  void printParams(std::ostream &os);

  inline void setDebugLevel(int value);
  inline void resetDebugLevel();
  inline int  getDebugLevel() const;

  inline void setDebugMinTimeStep(int value);
  inline void resetDebugMinTimeStep();
  inline int getDebugMinTimeStep() const;

  inline void setDebugMaxTimeStep(int value);
  inline void resetDebugMaxTimeStep();
  inline int getDebugMaxTimeStep() const;

  inline void setDebugMinTime(double value);
  inline void resetDebugMinTime();
  inline double getDebugMinTime() const;

  inline void setDebugMaxTime(double value);
  inline void resetDebugMaxTime();
  inline double getDebugMaxTime() const;

  inline void setScreenOutputFlag (bool bval);
  inline void resetScreenOutputFlag ();
  inline bool getScreenOutputFlag () const;

  inline void setMaskingFlag (bool bval);
  inline void resetMaskingFlag ();
  inline bool getMaskingFlag () const;

  inline void setMMFormat(bool value);
  inline void resetMMFormat();
  inline bool getMMFormat() const;

protected:
private:

public:

protected:
  // Print control flag
  bool printParamsFlag_;

  // command line parser
  const IO::CmdParse * commandLine_;

  // Calling solution method (e.g., DC Operation Point, Transient, etc.)
  AnalysisMode analysisMode_;
  bool         modeToggled_;

  // Nonlinear solution strategy
  NLStrategy nlStrategy_;

  // Damping method (e.g., Bank and Rose)
  LineSearchMethod searchMethod_;

  // Absolute convergence tolerance for the norm of the residual.
  double absTol_;

  // Relative convergence tolerance.
  double relTol_;

  // Weighted deltaX (update) norm tolerance (Petzold)
  double deltaXTol_;

  // Special deltaX (update) norm tolerance, used for the 
  // "small update" test.   
  double smallUpdateTol_;

  // Weighted RHS (residual) norm tolerance
  double RHSTol_;

  // Maximum number of solution steps to attempt.
  unsigned maxNewtonStep_;

  // maximum number of damping steps:
  unsigned maxSearchStep_;

  // inexact-Newton forcing flag:
  bool INForcingFlag_;

  // check device convergence flag:
  bool enforceDeviceConvFlag_;

  // inexact-Newton forcing term:
  double eta_;

  // norm level:
  int normLevel_;

  // linear optimization flag
  bool linearOptimization_;

  // Debug output options:
  int debugLevel_;
  int debugMinTimeStep_;
  int debugMaxTimeStep_;
  double debugMinTime_;
  double debugMaxTime_;
  bool screenOutputFlag_;
  bool matrixMarketFormat_;
  bool maskingFlag_;
};

//-----------------------------------------------------------------------------
// Function      : NLParams::setPrintParamsFlag
// Purpose       : Sets the flag which controls the printing of the nonlinear
//                 solver parameters list.  This flag should be set anytime one
//                 of the parameters is modified as the code runs.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
inline void NLParams::setPrintParamsFlag()
{
  printParamsFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::clearPrintParamsFlag
// Purpose       : Clears the flag which controls the printing of the nonlinear
//                 solver parameters list.  This flag is cleared once the
//                 output is printed.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
inline void NLParams::clearPrintParamsFlag()
{
  printParamsFlag_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getPrintParamsFlag
// Purpose       : Returns the flag value
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
inline bool NLParams::getPrintParamsFlag() const
{
  return printParamsFlag_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setNLStrategy
// Purpose       : Accessor method to set the nonlinear solver strategy.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline void NLParams::setNLStrategy(NLStrategy strategy)
{
  nlStrategy_ = strategy;

  // Make sure the method selected is supported, otherwise revert to the
  // default.
  if (nlStrategy_ != NEWTON &&
      nlStrategy_ != GRADIENT &&
      nlStrategy_ != NEWTON_GRADIENT &&
      nlStrategy_ != MOD_NEWTON &&
      nlStrategy_ != MOD_NEWTON_GRADIENT)
    resetNLStrategy();
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setNLStrategy
// Purpose       : Accessor method to set the nonlinear solver strategy.
// Special Notes : Takes an integer argument
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline void NLParams::setNLStrategy(int strategy)
{
  nlStrategy_ = static_cast<NLStrategy> (strategy);

  // Make sure the method selected is supported, otherwise revert to the
  // default.
  if (nlStrategy_ != NEWTON &&
      nlStrategy_ != GRADIENT &&
      nlStrategy_ != NEWTON_GRADIENT &&
      nlStrategy_ != MOD_NEWTON &&
      nlStrategy_ != MOD_NEWTON_GRADIENT)
    resetNLStrategy();
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetNLStrategy
// Purpose       : Accessor method to reset the nonlinear solver strategy to
//                 its default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline void NLParams::resetNLStrategy()
{
  nlStrategy_ = NEWTON;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getNLStrategy
// Purpose       : Accessor method to return the nonlinear solver strategy.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
inline NLStrategy NLParams::getNLStrategy() const
{
  return nlStrategy_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setSearchMethod
// Purpose       : Accessor method to set the line-search method
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
inline void NLParams::setSearchMethod(LineSearchMethod method)
{
  searchMethod_ = method;

  // Revert to the default if the method selected is not supported
  if (searchMethod_ != FULL &&
      searchMethod_ != DIVIDE &&
      searchMethod_ != BACKTRACK &&
      searchMethod_ != SIMPLE_BACKTRACK &&
      searchMethod_ != BANK_ROSE &&
      searchMethod_ != DESCENT)
    searchMethod_ = FULL;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setSearchMethod
// Purpose       : Accessor method to set the line-search method
// Special Notes : Takes and integer argument
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
inline void NLParams::setSearchMethod(int method)
{
  searchMethod_ = static_cast<LineSearchMethod> (method);

  // Revert to the default if the method selected is not supported
  if (searchMethod_ != FULL &&
      searchMethod_ != DIVIDE &&
      searchMethod_ != BACKTRACK &&
      searchMethod_ != SIMPLE_BACKTRACK &&
      searchMethod_ != BANK_ROSE &&
      searchMethod_ != DESCENT)
    searchMethod_ = FULL;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetSearchMethod
// Purpose       : Accessor method to reset the default line-search method
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::resetSearchMethod()
{
  searchMethod_ = FULL;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getSearchMethod
// Purpose       : Accessor method to return the line-search method
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/13/01
//-----------------------------------------------------------------------------
inline LineSearchMethod NLParams::getSearchMethod() const
{
  return searchMethod_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setDeltaXTol
// Purpose       : Accessor method to set the deltaX (update) tolerance.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/20/01
//-----------------------------------------------------------------------------
inline void NLParams::setDeltaXTol(double value)
{
  deltaXTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetDeltaXTol
// Purpose       : Accessor method to reset the deltaX (update) tolerance to
//                 the default value.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/20/01
//-----------------------------------------------------------------------------
inline void NLParams::resetDeltaXTol()
{
  deltaXTol_ = 1.0;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getDeltaXTol
// Purpose       : Accessor method to return the deltaX (update) tolerance
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/20/01
//-----------------------------------------------------------------------------
inline double NLParams::getDeltaXTol() const
{
  return deltaXTol_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setSmallUpdateTol
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/02/03
//-----------------------------------------------------------------------------
inline void NLParams::setSmallUpdateTol (double value)
{
  smallUpdateTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetSmallUpdateTol
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/02/03
//-----------------------------------------------------------------------------
inline void NLParams::resetSmallUpdateTol ()
{
  smallUpdateTol_ = 1.0e-6;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getSmallUpdateTol
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 07/02/03
//-----------------------------------------------------------------------------
inline double NLParams::getSmallUpdateTol () const
{
  return smallUpdateTol_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setRHSTol
// Purpose       : Accessor method to set the RHS (residual) tolerance.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::setRHSTol(double value)
{
  RHSTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetRHSTol
// Purpose       : Accessor method to reset the RHS (residual) tolerance to
//                 the default value.
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::resetRHSTol()
{
  RHSTol_ = 1.0e-06;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getRHSTol
// Purpose       : Accessor method to return the RHS (residual) tolerance
// Special Notes : This is used in combination with a WRMS norm.
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline double NLParams::getRHSTol() const
{
  return RHSTol_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setAbsTol
// Purpose       : Accessor method to set the absTol tolerance.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/22/01
//-----------------------------------------------------------------------------
inline void NLParams::setAbsTol(double value)
{
  absTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetAbsTol
// Purpose       : Accessor method to reset the absTol tolerance to
//                 the default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline void NLParams::resetAbsTol()
{
  absTol_ = 1.0e-12;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getAbsTol
// Purpose       : Accessor method to return the absTol tolerance
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline double NLParams::getAbsTol() const
{
  return absTol_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setRelTol
// Purpose       : Accessor method to set the relTol tolerance.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/22/01
//-----------------------------------------------------------------------------
inline void NLParams::setRelTol(double value)
{
  relTol_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetRelTol
// Purpose       : Accessor method to reset the relTol tolerance to
//                 the default value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline void NLParams::resetRelTol()
{
  relTol_ = 1.0e-03;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getRelTol
// Purpose       : Accessor method to return the relTol tolerance
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/21/01
//-----------------------------------------------------------------------------
inline double NLParams::getRelTol() const
{
  return relTol_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setMaxNewtonStep
// Purpose       : Accessor method to set the maximum number of Newton steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::setMaxNewtonStep(unsigned value)
{
  maxNewtonStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetMaxNewtonStep
// Purpose       : Accessor method to reset the maximum number of Newton steps
//                 to the default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::resetMaxNewtonStep()
{
  maxNewtonStep_ = 200;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getMaxNewtonStep
// Purpose       : Accessor method to get the maximum number of Newton steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline unsigned NLParams::getMaxNewtonStep() const
{
  return maxNewtonStep_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setMaxSearchStep
// Purpose       : Accessor method to set the maximum number of line-search
//                 steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::setMaxSearchStep(unsigned value)
{
  maxSearchStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetMaxSearchStep
// Purpose       : Accessor method to reset the maximum number of line-search
//                 steps to the default.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline void NLParams::resetMaxSearchStep()
{
  maxSearchStep_ = 9;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getMaxSearchStep
// Purpose       : Accessor method to get the maximum number of line-search
//                 steps.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/24/01
//-----------------------------------------------------------------------------
inline unsigned NLParams::getMaxSearchStep() const
{
  return maxSearchStep_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setForcingFlag
// Purpose       : Accessor method to set the inexact-Newton forcing flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/04/01
//-----------------------------------------------------------------------------
inline void NLParams::setForcingFlag(bool value)
{
  INForcingFlag_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetForcingFlag
// Purpose       : Accessor method to reset the default inexact-Newton forcing
//                 flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/04/01
//-----------------------------------------------------------------------------
inline void NLParams::resetForcingFlag()
{
    INForcingFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getForcingFlag
// Purpose       : Accessor method to get the inexact-Newton forcing flag.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/04/01
//-----------------------------------------------------------------------------
inline bool NLParams::getForcingFlag() const
{
  return INForcingFlag_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setEnforceDeviceConvFlag
// Purpose       : Accessor method to set the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 06/30/05
//-----------------------------------------------------------------------------
inline void NLParams::setEnforceDeviceConvFlag (bool flag)
{
  enforceDeviceConvFlag_ = flag;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetEnforceDeviceConvFlag
// Purpose       : Accessor method to reset the default inexact-Newton forcing
//                 value.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 06/30/05
//-----------------------------------------------------------------------------
inline void NLParams::resetEnforceDeviceConvFlag ()
{
  enforceDeviceConvFlag_ = true;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getEnforceDeviceConvFlag
// Purpose       : Accessor method to get the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 06/30/05
//-----------------------------------------------------------------------------
inline bool NLParams::getEnforceDeviceConvFlag () const
{
  return enforceDeviceConvFlag_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setForcingTerm
// Purpose       : Accessor method to set the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/29/01
//-----------------------------------------------------------------------------
inline void NLParams::setForcingTerm(double value)
{
  eta_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetForcingTerm
// Purpose       : Accessor method to reset the default inexact-Newton forcing
//                 value.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/29/01
//-----------------------------------------------------------------------------
inline void NLParams::resetForcingTerm()
{
  eta_ = 1.0e-01;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getForcingTerm
// Purpose       : Accessor method to get the inexact-Newton forcing term.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 06/29/01
//-----------------------------------------------------------------------------
inline double NLParams::getForcingTerm() const
{
  return eta_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::setDebugLevel(int value)
{
  debugLevel_ = value;
  setNonlinearDebugLevel(value);
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::resetDebugLevel()
{
  debugLevel_ = 1;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline int  NLParams::getDebugLevel() const
{
  return debugLevel_;
}


//-----------------------------------------------------------------------------
// Function      : NLParams::setDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::setDebugMinTimeStep(int value)
{
  debugMinTimeStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::resetDebugMinTimeStep()
{
  debugMinTimeStep_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline int NLParams::getDebugMinTimeStep() const
{
  return debugMinTimeStep_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::setDebugMaxTimeStep(int value)
{
  debugMaxTimeStep_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::resetDebugMaxTimeStep()
{
  debugMaxTimeStep_ = Util::MachineDependentParams::IntMax(); 
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline int NLParams::getDebugMaxTimeStep() const
{
  return debugMaxTimeStep_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::setDebugMinTime(double value)
{
  debugMinTime_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::resetDebugMinTime()
{
  debugMinTime_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline double NLParams::getDebugMinTime() const
{
  return debugMinTime_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::setDebugMaxTime(double value)
{
  debugMaxTime_ = value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline void NLParams::resetDebugMaxTime()
{
  debugMaxTime_ = Util::MachineDependentParams::DoubleMax();
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
inline double NLParams::getDebugMaxTime() const
{
  return debugMaxTime_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 04/25/04
//-----------------------------------------------------------------------------
inline void NLParams::setScreenOutputFlag (bool bval)
{
  screenOutputFlag_ = bval;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 04/25/04
//-----------------------------------------------------------------------------
inline void NLParams::resetScreenOutputFlag ()
{
  screenOutputFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 04/25/04
//-----------------------------------------------------------------------------
inline bool NLParams::getScreenOutputFlag () const
{
  return screenOutputFlag_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setMaskingFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 10/31/2014
//-----------------------------------------------------------------------------
inline void NLParams::setMaskingFlag (bool bval)
{
  maskingFlag_ = bval;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetMaskingFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 10/31/2014
//-----------------------------------------------------------------------------
inline void NLParams::resetMaskingFlag ()
{
  maskingFlag_ = false;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getMaskingFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 10/31/2014
//-----------------------------------------------------------------------------
inline bool NLParams::getMaskingFlag () const
{
  return maskingFlag_;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
inline void NLParams::setMMFormat(bool value)
{
  matrixMarketFormat_=value;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::resetMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
inline void NLParams::resetMMFormat()
{
  matrixMarketFormat_=false;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::getMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
inline bool NLParams::getMMFormat() const
{
  return matrixMarketFormat_;
}

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_NLParams_h

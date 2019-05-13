//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose       : This file implements the class associated with all user
//                 specified parameters which relate to the time integration
//                 algorithms and problem definition.
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

#include <Xyce_config.h>

#include <iostream>

#include <N_NLS_fwd.h>

#include <N_TIA_TIAParams.h>
#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace TimeIntg {

int maxOrder(
  const IO::CmdParse &  command_line)
{
  return command_line.getArgumentIntValue("-maxord", 2);
}


//-----------------------------------------------------------------------------
// Function      : TIAParams::TIAParams
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------

TIAParams::TIAParams()
  : initialTime(0.0),
    finalTime(0.0),
    initialTimeStep(1.0e-10),
    minTimeStep(0.0),
    minTimeStepGiven(false),
    minTimeStepsBP(10),
    minTimeStepsBPGiven(false),
    maxTimeStep(1.0e+99),
    maxTimeStepGiven(false),
    constantTimeStepFlag(false),
    restartTimeStepScale(0.005),
    initialOutputTime(0.0),
    initialOutputTimeGiven(false),
    useDeviceTimeStepMaxFlag(true),
    errorAnalysisOption(LOCAL_TRUNCATED_ESTIMATES),
    errorAnalysisOptionResetCount(0),
    bpEnable(true),
    NLmin(3),
    NLmax(8),
    delmax(1.0e+99),
    delmaxGiven(false),
    timestepsReversal(false),
    testFirstStep(false),
    newBPStepping(true),
    maskIVars(false),
    newLte(1),
    relErrorTol(1.0e-3),
    relErrorTolGiven(false),
    absErrorTol(1.0e-6),
    errTolAcceptance(1.0),
    jacLimitFlag(false),
    jacLimit(1.0e+17),
    maxOrder(2),
    minOrder(1),
    interpOutputFlag(true),
    minTimeStepRecoveryCounter(0),
    bpPrune(true)
{}

//-----------------------------------------------------------------------------
// Function      : TIAParams::TIAParams
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
TIAParams::~TIAParams()
{}

//-----------------------------------------------------------------------------
// Function      : TIAParams::printParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233.
// Creation Date : 7/12/01
//-----------------------------------------------------------------------------
void TIAParams::printParams(std::ostream &os, int analysis) const
{
  os << "\n" << std::endl;
  os << Xyce::section_divider << std::endl;
  os << "\n***** Time Integration solver options:\n" << std::endl;

  if (analysis == Xyce::Nonlinear::TRANSIENT)
  {
    os << "\tAnalysis:\t\t\tTRANSIENT" << std::endl
       << "\tInitial Time (sec):\t\t" <<  initialTime << std::endl
       << "\tFinal Time (sec):\t\t" <<  finalTime << std::endl
       << "\tStarting Time Step(sec):\t" <<  initialTimeStep << std::endl
       << "\tRestart Time Step Scale:\t" <<  restartTimeStepScale << std::endl
       << "\tError Analysis option:\t" <<  errorAnalysisOption << std::endl
       << "\tInitial Output Time:\t" <<  initialOutputTime << std::endl
      // << "\tTime Integration method:\t" << integrationMethod << std::endl

       << (constantTimeStepFlag ? "\tUsing Constant Step Size" : "\tUsing Variable Step Size") << std::endl
       << (useDeviceTimeStepMaxFlag ? "\tUsing Device specified maximum stepsize" : "\tNOT using Device specified maximum stepsize") << std::endl;
      // << (nlNearConvFlag ? "\tNL Near Convergence Flag is ON" : "\tNL Near Convergence Flag is OFF") << std::endl
      // << (passNLStall ? "\tNL Pass Non-linear Stalls is ON" : "\tNL Pass Non-linear Stalls is OFF") << std::endl;
  }
  else
  {
    os << "\tAnalysis:\t\t\tDC SWEEP" << std::endl;
  }

  os << "\tabsErrorTol:\t\t\t" <<  absErrorTol  << std::endl
     << "\trelErrorTol:\t\t\t" <<  relErrorTol  << std::endl
     << "\tMaximum Order:\t\t\t" <<  maxOrder << std::endl
     << "\tMinimum Order:\t\t\t" <<  minOrder << std::endl
     << "\tInterpolated Output Flag:\t\t " << (interpOutputFlag ? "true": "false") << std::endl;
  
    // << "\tConductance Test Flag:\t\t" << (condTestFlag  ? "true": "false") << std::endl;

  // if (condTestDeviceNames.empty())
  // {
  //   os << "\tConductance Test Device Name List is:\t\tEMPTY" << std::endl;
  // }
  // else
  // {
  //   os << "\tConductance Test Device Name List contains:  " << std::endl;
  //   for (std::list< std::string >::iterator it = condTestDeviceNames.begin(); it != condTestDeviceNames.end(); ++it)
  //   {
  //     os << " \"" << *it << "\"";
  //   }
  // }
  os << Xyce::section_divider << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : TIAParams::operator=
// Purpose       : "=" operator.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/17/05
//-----------------------------------------------------------------------------
TIAParams::TIAParams(const TIAParams &right)
  : initialTime(right.initialTime),
    finalTime(right.finalTime),
    initialTimeStep(right.initialTimeStep),
    minTimeStep(right.minTimeStep),
    minTimeStepGiven(right.minTimeStepGiven),
    maxTimeStep(right.maxTimeStep),
    maxTimeStepGiven(right.maxTimeStepGiven),
    constantTimeStepFlag(right.constantTimeStepFlag),
    restartTimeStepScale(right.restartTimeStepScale),
    initialOutputTime(right.initialOutputTime),
    initialOutputTimeGiven(right.initialOutputTimeGiven),
    useDeviceTimeStepMaxFlag(right.useDeviceTimeStepMaxFlag),
    errorAnalysisOption(right.errorAnalysisOption),
    errorAnalysisOptionResetCount(right.errorAnalysisOptionResetCount),
    bpEnable(right.bpEnable),
    relErrorTol(right.relErrorTol),
    absErrorTol(right.absErrorTol),
    errTolAcceptance(right.errTolAcceptance),
    jacLimitFlag(right.jacLimitFlag),
    jacLimit(right.jacLimit),
    maxOrder(right.maxOrder),
    minOrder(right.minOrder),
    interpOutputFlag(right.interpOutputFlag)
{}

//-----------------------------------------------------------------------------
// Function      : TIAParams::operator=
// Purpose       : "=" operator.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/17/05
//-----------------------------------------------------------------------------
TIAParams & TIAParams::operator=(const TIAParams &right)
{
  if (this != &right)
  {
    initialOutputTime = right.initialOutputTime;
    initialOutputTimeGiven = right.initialOutputTimeGiven;
    initialTime = right.initialTime;
    finalTime = right.finalTime;
    initialTimeStep = right.initialTimeStep;
    maxTimeStep = right.maxTimeStep;
    maxTimeStepGiven = right.maxTimeStepGiven;
    minTimeStep = right.minTimeStep;
    minTimeStepGiven = right.minTimeStepGiven;
    constantTimeStepFlag = right.constantTimeStepFlag;
    useDeviceTimeStepMaxFlag = right.useDeviceTimeStepMaxFlag;
    errorAnalysisOption = right.errorAnalysisOption;
    errorAnalysisOptionResetCount = right.errorAnalysisOptionResetCount;
    bpEnable = right.bpEnable;
    restartTimeStepScale = right.restartTimeStepScale;
    relErrorTol = right.relErrorTol;
    absErrorTol = right.absErrorTol;
    errTolAcceptance = right.errTolAcceptance;
    jacLimitFlag = right.jacLimitFlag;
    jacLimit = right.jacLimit;
    maxOrder = right.maxOrder;
    minOrder = right.minOrder;
    interpOutputFlag = right.interpOutputFlag;
  }

  return *this;
}

void
TIAParams::setMaxOrder(
  int           max_order)
{
  if (max_order < 1)
    max_order = 1;
  else if (max_order > 2)
    max_order = 2;

  maxOrder = max_order;
}

//-----------------------------------------------------------------------------
// Function      : setTimeIntegratorOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool TIAParams::setTimeIntegratorOption(
  const Util::Param &           param)
{
  bool value_set =
    Util::setValue(param, "USEDEVICEMAX", useDeviceTimeStepMaxFlag)
    || setValue(param, "RELTOL", relErrorTol, relErrorTolGiven)
    || setValue(param, "ABSTOL", absErrorTol)
    || setValue(param, "BPENABLE", bpEnable)
    || setValue(param, "RESTARTSTEPSCALE", restartTimeStepScale)
    || setValue(param, "MINTIMESTEPSBP", minTimeStepsBP, minTimeStepsBPGiven)
    || setValue(param, "JACLIMITFLAG", jacLimitFlag)
    || setValue(param, "JACLIMIT", jacLimit)
    || setValue(param, "MAXORD", maxOrder)
    || setValue(param, "MINORD", minOrder)
    || setValue(param, "TIMESTEPSREVERSAL", timestepsReversal)
    || setValue(param, "TESTFIRSTSTEP", testFirstStep)
    || setValue(param, "DELMAX", delmax, delmaxGiven)
    || setValue(param, "NLMIN", NLmin)
    || setValue(param, "NLMAX", NLmax)
    || setValue(param, "NEWLTE", newLte)
    || setValue(param, "NEWBPSTEPPING", newBPStepping) 
    || setValue(param, "MASKIVARS", maskIVars)
    || setValue(param, "INTERPOUTPUT", interpOutputFlag)
    || setValue(param, "DTMIN", minTimeStep, minTimeStepGiven)
    || setValue(param, "MINTIMESTEPRECOVERY", minTimeStepRecoveryCounter)
    || setValue(param, "CONSTSTEP", constantTimeStepFlag)
    || setValue(param, "ERROPTION", errorAnalysisOption)
    || setValue(param, "BPPRUNE", bpPrune);

  if (value_set)
    ;
  else
  {
    return false;
  }

  if (NLmin > NLmax)
  {
    Report::UserError() << ".options timeint NLMIN = " << NLmin << " > " << NLmax << " = NLMAX!";
  }

  if (NLmin > NLmax)
  {
    Report::UserError() << ".options timeint NLMIN = " << NLmin << " > " << NLmax << " = NLMAX!";
  }

  if  (newLte < 0 || newLte > 3)
    Report::UserError() <<  "Unsupported NEWLTE type";

  // tscoffe/tmei 12/7/07:  If error option = 1 (NO LTE) then make sure minTimeStepBP is enabled.
  if (errorAnalysisOption == TimeIntg::NO_LOCAL_TRUNCATED_ESTIMATES)
  {
    minTimeStepsBPGiven = true;
  }

  if (minTimeStepRecoveryCounter > 0)
  {
    Xyce::lout() << "The time step recovery algorithm, set by .OPTIONS TIMEINT MINTIMESTEPRECOVERY>0, is considered deprecated.  It will be removed from a future version of Xyce." <<std::endl;
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : setAnalysisOptions
// Purpose       :
// Special Notes : These are from '.options timeint'
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool TIAParams::setAnalysisOption(
  const Util::Param &           param)
{
  bool value_set =
    Util::setValue(param,"TSTART", initialOutputTime, initialOutputTimeGiven)
    || Util::setValue(param, "TSTOP", finalTime)
    || Util::setValue(param, "TSTEP", initialTimeStep)
    || Util::setValue(param, "DTMAX", maxTimeStep, maxTimeStepGiven);

  if (value_set)
    ;
  else
  {
    return false;
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << "setting maxTimeStep = " << maxTimeStep << std::endl;
  }

  return true;
}

void
TIAParams::populateMetadata(
  IO::PkgOptionsMgr &options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("TIMEINT");

    parameters.insert(Util::ParamMap::value_type("METHOD", Util::Param("METHOD", 1)));
    if (DEBUG_ANALYSIS)
      parameters.insert(Util::ParamMap::value_type("CONSTSTEP", Util::Param("CONSTSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("USEDEVICEMAX", Util::Param("USEDEVICEMAX", 1)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("RESTARTSTEPSCALE", Util::Param("RESTARTSTEPSCALE", .005)));
//  parameters.insert(Util::ParamMap::value_type("NLNEARCONV", Util::Param("NLNEARCONV", 1)));
    parameters.insert(Util::ParamMap::value_type("NLNEARCONV", Util::Param("NLNEARCONV", 0)));
    parameters.insert(Util::ParamMap::value_type("NLSMALLUPDATE", Util::Param("NLSMALLUPDATE", 1)));
    parameters.insert(Util::ParamMap::value_type("DOUBLEDCOP", Util::Param("DOUBLEDCOP", "")));
    parameters.insert(Util::ParamMap::value_type("RESETTRANNLS", Util::Param("RESETTRANNLS", 1)));
    parameters.insert(Util::ParamMap::value_type("BPENABLE", Util::Param("BPENABLE", 1)));
    parameters.insert(Util::ParamMap::value_type("EXITTIME", Util::Param("EXITTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("EXITSTEP", Util::Param("EXITSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("ERROPTION", Util::Param("ERROPTION", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("JACLIMITFLAG", Util::Param("JACLIMITFLAG", 0)));
    parameters.insert(Util::ParamMap::value_type("JACLIMIT", Util::Param("JACLIMIT", 1.0e17)));
    parameters.insert(Util::ParamMap::value_type("DAESTATEDERIV", Util::Param("DAESTATEDERIV", 0)));
    parameters.insert(Util::ParamMap::value_type("TESTFIRSTSTEP", Util::Param("TESTFIRSTSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DTMIN", Util::Param("DTMIN", 0.0)));
    parameters.insert(Util::ParamMap::value_type("NEWBPSTEPPING", Util::Param("NEWBPSTEPPING", 0)));
    parameters.insert(Util::ParamMap::value_type("MASKIVARS", Util::Param("MASKIVARS",  0)));  
    parameters.insert(Util::ParamMap::value_type("MINTIMESTEPSBP", Util::Param("MINTIMESTEPSBP", 10)));
    parameters.insert(Util::ParamMap::value_type("NEWLTE", Util::Param("NEWLTE", 1)));
    parameters.insert(Util::ParamMap::value_type("MAXORD", Util::Param("MAXORD", 2)));
    parameters.insert(Util::ParamMap::value_type("MINORD", Util::Param("MINORD", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTINTERPMPDE", Util::Param("OUTPUTINTERPMPDE", 1)));
    parameters.insert(Util::ParamMap::value_type("INTERPOUTPUT", Util::Param("INTERPOUTPUT", 1)));
    parameters.insert(Util::ParamMap::value_type("CONDTEST", Util::Param("CONDTEST", 0)));
    parameters.insert(Util::ParamMap::value_type("CONDTESTDEVICENAME", Util::Param("CONDTESTDEVICENAME", "dev_name")));
    parameters.insert(Util::ParamMap::value_type("ISOCONDTEST", Util::Param("ISOCONDTEST", 0)));
    parameters.insert(Util::ParamMap::value_type("ISOCONDTESTDEVICENAME", Util::Param("ISOCONDTESTDEVICENAME", "dev_name")));
    parameters.insert(Util::ParamMap::value_type("PASSNLSTALL", Util::Param("PASSNLSTALL", false)));
    parameters.insert(Util::ParamMap::value_type("NLMIN", Util::Param("NLMIN", 3)));
    parameters.insert(Util::ParamMap::value_type("NLMAX", Util::Param("NLMAX", 8)));
    parameters.insert(Util::ParamMap::value_type("DELMAX", Util::Param("DELMAX", 1.0e+99)));
    parameters.insert(Util::ParamMap::value_type("TIMESTEPSREVERSAL", Util::Param("TIMESTEPSREVERSAL", false)));
    parameters.insert(Util::ParamMap::value_type("MINTIMESTEPRECOVERY", Util::Param("MINTIMESTEPRECOVERY", 0)));
    parameters.insert(Util::ParamMap::value_type("HISTORYTRACKINGDEPTH", Util::Param("HISTORYTRACKINGDEPTH", 50)));

    parameters.insert(Util::ParamMap::value_type("BREAKPOINTS", Util::Param("BREAKPOINTS", "VECTOR")));

    parameters.insert(Util::ParamMap::value_type("BPPRUNE", Util::Param("BPPRUNE", false)));
  }
}

} // namespace TimeIntg
} // namespace Xyce

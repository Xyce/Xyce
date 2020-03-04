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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/13/01
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_UTL_Param.h>
#include <N_NLS_DampedNewton.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : NLParams::NLParams
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/09/00
//-----------------------------------------------------------------------------
NLParams::NLParams(AnalysisMode mode, const IO::CmdParse & cp)
  : printParamsFlag_(true),
    commandLine_(&cp),
    analysisMode_(mode),
    modeToggled_(true),
    debugLevel_(1),
    debugMinTimeStep_(0),
    debugMaxTimeStep_(Util::MachineDependentParams::IntMax()),
    debugMinTime_(0.0),
    debugMaxTime_(Util::MachineDependentParams::DoubleMax()),
    screenOutputFlag_(false),
    matrixMarketFormat_(false),
    maskingFlag_(false)
{
  // Set default update norm tolerance
  resetDeltaXTol();

  // Set default small update norm tolerance
  resetSmallUpdateTol ();

  // Set default residual norm tolerance
  resetRHSTol();

  // Set default absolute tolerance value for use in weighted norm
  resetAbsTol();

  // Set default relative tolerance value for use in weighted norm
  resetRelTol();

  resetEnforceDeviceConvFlag ();

  resetSearchMethod();
  resetDirection();
  resetNLStrategy();
  resetMaxNewtonStep();
  resetMaxSearchStep();
  resetForcingFlag();
  resetForcingTerm();
  resetNormLevel();
  resetConstraintBT();
  resetGlobalBTMax();
  resetGlobalBTMin();
  resetGlobalBTChange();

  // Set the default parameters for transient, if the specified mode
  // is TRANSIENT.
  // The defaults set in the NLParams constructor are for DC_OP,
  // so the transient ones should be reset here.
  if  (mode==TRANSIENT)
  {
    setNLStrategy(NEWTON);
    setSearchMethod(FULL);
    setMaxSearchStep(2);
    setMaxNewtonStep(20);
    setDeltaXTol(0.33);
    setForcingFlag(false);
    setAbsTol(1.0e-06);
    setRelTol(1.0e-02);
    setRHSTol(1.0e-02);
  }

  if  (mode == HB_MODE)
  {
    setNLStrategy(NEWTON);
    setSearchMethod(FULL);
    setAbsTol(1.0e-9);
    setRHSTol(1.0e-4);
  }


  setCmdLineOptions ();
}

//-----------------------------------------------------------------------------
// Function      : NLParams::NLParams
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/09/00
//-----------------------------------------------------------------------------
NLParams::NLParams(const NLParams & right)
  : printParamsFlag_(right.printParamsFlag_),
    commandLine_(right.commandLine_),
    analysisMode_(right.analysisMode_),
    modeToggled_(right.modeToggled_),
    nlStrategy_(right.nlStrategy_),
    searchMethod_(right.searchMethod_),
    direction_(right.direction_),
    absTol_(right.absTol_),
    relTol_(right.relTol_),
    deltaXTol_(right.deltaXTol_),
    RHSTol_(right.RHSTol_),
    maxNewtonStep_(right.maxNewtonStep_),
    maxSearchStep_(right.maxSearchStep_),
    INForcingFlag_(right.INForcingFlag_),
    eta_(right.eta_),
    normLevel_(right.normLevel_),
    linearOptimization_(right.linearOptimization_),
    constraintBT_(right.constraintBT_),
    globalBTMax_(right.globalBTMax_),
    globalBTMin_(right.globalBTMin_),
    globalBTChange_(right.globalBTChange_),
    debugLevel_(right.debugLevel_),
    debugMinTimeStep_(right.debugMinTimeStep_),
    debugMaxTimeStep_(right.debugMaxTimeStep_),
    debugMinTime_(right.debugMinTime_),
    debugMaxTime_(right.debugMaxTime_)
{
}

//-----------------------------------------------------------------------------
// Function      : NLParams::~NLParams
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/09/00
//-----------------------------------------------------------------------------
NLParams::~NLParams()
{
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setOptions
// Purpose       : This function takes an .options statement option block,
//                 and uses it to set the various parameters of
//                 the NLParams class.
//
// Special Notes : This was originally in DampedNewton, but it makes
//                 more sense to have it here.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool NLParams::setOptions (const Util::OptionBlock & OB)
{
  for (Util::ParamList::const_iterator it_tpL = OB.begin();
       it_tpL != OB.end(); ++ it_tpL)
  {
    if (it_tpL->uTag() == "ABSTOL")
    {
      setAbsTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "RELTOL")
    {
      setRelTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "DELTAXTOL")
    {
      setDeltaXTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "SMALLUPDATETOL")
    {
      setSmallUpdateTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "ENFORCEDEVICECONV")
    {
      setEnforceDeviceConvFlag(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "RHSTOL")
    {
      setRHSTol(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "MAXSTEP")
    {
      setMaxNewtonStep(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTBT")
    {
      setConstraintBT(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTMAX")
    {
      setGlobalBTMax(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTMIN")
    {
      setGlobalBTMin(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "CONSTRAINTCHANGE")
    {
      setGlobalBTChange(it_tpL->getImmutableValue<double>());
    }
    else if (it_tpL->uTag() == "NLSTRATEGY")
    {
      setNLStrategy(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "SEARCHMETHOD")
    {
      setSearchMethod(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "MAXSEARCHSTEP")
    {
      setMaxSearchStep(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "IN_FORCING")
    {
      setForcingFlag(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "NORMLVL")
    {
      setNormLevel(it_tpL->getImmutableValue<int>());
    }
    else if (it_tpL->uTag() == "NOX")
    {
      // do nothing.
    }
    else if (it_tpL->uTag() == "MATRIXMARKET")
    {
      if (DEBUG_NONLINEAR)
      {
        setMMFormat (static_cast<bool>(it_tpL->getImmutableValue<double>()));
      }
    }
    else if (it_tpL->uTag() == "DEBUGLEVEL")
    {
      if (DEBUG_NONLINEAR)
      {
        setDebugLevel(it_tpL->getImmutableValue<int>());
      }
    }
    else if (it_tpL->uTag() == "DEBUGMINTIMESTEP")
    {
      if (DEBUG_NONLINEAR)
      {
        setDebugMinTimeStep(it_tpL->getImmutableValue<int>());
      }
    }
    else if (it_tpL->uTag() == "DEBUGMAXTIMESTEP")
    {
      if (DEBUG_NONLINEAR)
      {
        setDebugMaxTimeStep(it_tpL->getImmutableValue<int>());
      }
    }
    else if (it_tpL->uTag() == "DEBUGMINTIME")
    {
      if (DEBUG_NONLINEAR)
      {
        setDebugMinTime(it_tpL->getImmutableValue<double>());
      }
    }
    else if (it_tpL->uTag() == "DEBUGMAXTIME")
    {
      if (DEBUG_NONLINEAR)
      {
        setDebugMaxTime(it_tpL->getImmutableValue<double>());
      }
    }
    else if (it_tpL->uTag() == "SCREENOUTPUT")
    {
      if (DEBUG_NONLINEAR)
      {
        setScreenOutputFlag (static_cast<bool>(it_tpL->getImmutableValue<double>()));
      }
    }
    else if (it_tpL->uTag() == "USEMASKING")
    {
      setMaskingFlag(static_cast<bool>(it_tpL->getImmutableValue<double>()));
    }
    else
    {
      Xyce::Report::UserFatal0() <<  it_tpL->uTag()
                                 << " is not a recognized nonlinear solver option.";
    }
  }

  setCmdLineOptions ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::printParams
// Purpose       : Print out the nonlinear solver parameter values.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 01/16/01
//-----------------------------------------------------------------------------
void NLParams::printParams(std::ostream &os)
{
  os << "\n" << std::endl
     << Xyce::section_divider << std::endl;
  os << "\n***** Nonlinear solver options:\n" << std::endl
     << "\tabsTol:\t\t\t" << getAbsTol() << std::endl
     << "\trelTol:\t\t\t" << getRelTol() << std::endl
     << "\tdeltaXTol (weighted):\t" << getDeltaXTol() << std::endl
     << "\tRHSTol:\t\t\t" << getRHSTol() << std::endl
     << "\tSmall Update Tol:\t" << getSmallUpdateTol() << std::endl
     << "\tmax NL Steps:\t\t" << getMaxNewtonStep() << std::endl;

  if (analysisMode_ == DC_OP)
    os << "\tAnalysis Mode:\t\t" << analysisMode_ << "\t(DC Op)" << std::endl;
  else if (analysisMode_ == DC_SWEEP)
    os << "\tAnalysis Mode:\t\t" << analysisMode_ <<  "\t(DC Sweep)" << std::endl;
  else if (analysisMode_ == TRANSIENT)
    os << "\tAnalysis Mode:\t\t" << analysisMode_ << "\t(Transient)" << std::endl;

  NLStrategy strategy = getNLStrategy();
  if (strategy == NEWTON)
    os << "\tNL Strategy:\t\t" << strategy << "\t(None => Newton)" << std::endl;
  else if (strategy == GRADIENT)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Gradient)" << std::endl;
  else if (strategy == NEWTON_GRADIENT)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Newton/Gradient)" << std::endl;
  else if (strategy == MOD_NEWTON)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Modified-Newton)" << std::endl;
  else if (strategy == MOD_NEWTON_GRADIENT)
    os << "\tNL Strategy:\t\t" << strategy << "\t(Modified-Newton/Gradient)" << std::endl;

  LineSearchMethod searchMethod = getSearchMethod();
  if (searchMethod == FULL)
    os << "\tsearch method:\t\t" << searchMethod << "\t(None => Full Newton Steps)" << std::endl;

  else if (searchMethod == DIVIDE)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Divide)" << std::endl;

  else if (searchMethod == BACKTRACK)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Backtrack)" << std::endl;

  else if (searchMethod == SIMPLE_BACKTRACK)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Simple Backtrack)" << std::endl;

  else if (searchMethod == BANK_ROSE)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Bank and Rose Algorithm)" << std::endl;

  else if (searchMethod == DESCENT)
    os << "\tsearch method:\t\t" << searchMethod << "\t(Line Search)" << std::endl;

  os << "\tmax search steps:\t" << getMaxSearchStep() << std::endl
     << "\tinexact-Newton forcing:\t" << getForcingFlag() << std::endl
     << "\tnorm level:\t\t" << getNormLevel() << std::endl
     << "\tconstraint backtrack:\t" << getConstraintBT() << std::endl;

  if (DEBUG_NONLINEAR)
  {
    os << "\tdebugLevel:\t\t" << getDebugLevel () <<std::endl
       << "\tdebugMinTimeStep:\t" << getDebugMinTimeStep () <<std::endl
       << "\tdebugMaxTimeStep:\t" << getDebugMaxTimeStep () <<std::endl;
  }

  os << Xyce::section_divider << "\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::operator=
// Purpose       : "=" operator.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 09/05/01
//-----------------------------------------------------------------------------
NLParams & NLParams::operator=(const NLParams & right)
{
  if (this != &right)
  {
    commandLine_   = right.commandLine_;
    nlStrategy_    = right.nlStrategy_;
    searchMethod_  = right.searchMethod_;
    direction_     = right.direction_;
    deltaXTol_     = right.deltaXTol_;
    RHSTol_        = right.RHSTol_;
    absTol_        = right.absTol_;
    relTol_        = right.relTol_;
    maxNewtonStep_ = right.maxNewtonStep_;
    maxSearchStep_ = right.maxSearchStep_;
    INForcingFlag_ = right.INForcingFlag_;
    eta_           = right.eta_;
    normLevel_     = right.normLevel_;

    analysisMode_  = right.analysisMode_;

    linearOptimization_ = right.linearOptimization_;

    constraintBT_   = right.constraintBT_;
    globalBTMax_    = right.globalBTMax_;
    globalBTMin_    = right.globalBTMin_;
    globalBTChange_ = right.globalBTChange_;

    // Debug output options:
    debugLevel_       = right.debugLevel_;
    debugMinTimeStep_ = right.debugMinTimeStep_;
    debugMaxTimeStep_ = right.debugMaxTimeStep_;
    debugMinTime_     = right.debugMinTime_;
    debugMaxTime_     = right.debugMaxTime_;
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : NLParams::setCmdLineOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 05/24/05
//-----------------------------------------------------------------------------

bool NLParams::setCmdLineOptions ()
{
  // set (or override) debug levels based on command line options
  if (DEBUG_NONLINEAR && commandLine_->argExists( "-ndl" ))
    setDebugLevel(commandLine_->getArgumentIntValue("-ndl", getDebugLevel()));

  return true;
}

void NLParams::populateMetadata(IO::PkgOptionsMgr& options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("NONLIN");

  parameters.insert(Util::ParamMap::value_type("NLSTRATEGY", Util::Param("NLSTRATEGY", 0)));
  parameters.insert(Util::ParamMap::value_type("SEARCHMETHOD", Util::Param("SEARCHMETHOD", 0)));
  parameters.insert(Util::ParamMap::value_type("NOX", Util::Param("NOX", 1)));
  parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-12)));
  parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-3)));
  parameters.insert(Util::ParamMap::value_type("DELTAXTOL", Util::Param("DELTAXTOL", 1.0)));
  parameters.insert(Util::ParamMap::value_type("SMALLUPDATETOL", Util::Param("SMALLUPDATETOL", 1.0e-6)));
  parameters.insert(Util::ParamMap::value_type("RHSTOL", Util::Param("RHSTOL", 1.0E-6)));
  parameters.insert(Util::ParamMap::value_type("MAXSTEP", Util::Param("MAXSTEP", 200)));
  parameters.insert(Util::ParamMap::value_type("MAXSEARCHSTEP", Util::Param("MAXSEARCHSTEP", 0)));
  parameters.insert(Util::ParamMap::value_type("NORMLVL", Util::Param("NORMLVL", 2)));
  parameters.insert(Util::ParamMap::value_type("LINOPT", Util::Param("LINOPT", 0)));
  parameters.insert(Util::ParamMap::value_type("CONSTRAINTBT", Util::Param("CONSTRAINTBT", 0)));
  parameters.insert(Util::ParamMap::value_type("CONSTRAINTMAX", Util::Param("CONSTRAINTMAX", Util::MachineDependentParams::DoubleMax())));
  parameters.insert(Util::ParamMap::value_type("CONSTRAINTMIN", Util::Param("CONSTRAINTMIN", -Util::MachineDependentParams::DoubleMax())));
  parameters.insert(Util::ParamMap::value_type("CONSTRAINTCHANGE", Util::Param("CONSTRAINTCHANGE", 0.0)));
  parameters.insert(Util::ParamMap::value_type("IN_FORCING", Util::Param("IN_FORCING", 0)));
  parameters.insert(Util::ParamMap::value_type("AZ_TOL", Util::Param("AZ_TOL", 1.0E-12)));
  parameters.insert(Util::ParamMap::value_type("DLSDEBUG", Util::Param("DLSDEBUG", 0)));
  parameters.insert(Util::ParamMap::value_type("MATRIXMARKET", Util::Param("MATRIXMARKET", 0)));
  parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 99999999)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 1.0E99)));
  parameters.insert(Util::ParamMap::value_type("SCREENOUTPUT", Util::Param("SCREENOUTPUT", 0)));
  parameters.insert(Util::ParamMap::value_type("USEMASKING", Util::Param("USEMASKING", 0)));
  parameters.insert(Util::ParamMap::value_type("RECOVERYSTEPTYPE", Util::Param("RECOVERYSTEPTYPE", 0)));
  parameters.insert(Util::ParamMap::value_type("RECOVERYSTEP", Util::Param("RECOVERYSTEP", 1.0)));
  parameters.insert(Util::ParamMap::value_type("MEMORY", Util::Param("MEMORY", 400)));
  parameters.insert(Util::ParamMap::value_type("CONTINUATION", Util::Param("CONTINUATION", 0)));
  parameters.insert(Util::ParamMap::value_type("ENFORCEDEVICECONV", Util::Param("ENFORCEDEVICECONV", 1)));
}

} // namespace Nonlinear
} // namespace Xyce


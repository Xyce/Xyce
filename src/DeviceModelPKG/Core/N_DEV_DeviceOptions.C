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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/06/01
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <climits>

#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_ExtendedString.h>

#ifndef Xyce_NEW_EXCESS_PHASE
static const bool EXCESS_PHASE = false;
#else
static const bool EXCESS_PHASE = true;
#endif

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::DeviceOptions
// Purpose       : constructor
//
// Special Notes : The are initialized to default values, but are reset
//                 in the setupDefaultOptions function.  Confusing I know,
//                 but I had a reason for doing this, I think.
//
//                 Consider  setupDefaultOptions function to be the
//                 ultimate authority.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/01/02
//-----------------------------------------------------------------------------
DeviceOptions::DeviceOptions()
  : defad (0.0e+0),  // MOS drain diffusion area.
    defas (0.0e+0),  // MOS source diffusion area.
    defl  (1.0e-4),  // MOS channel length.
    defw  (1.0e-4),  // MOS channel width.
    modelBinningFlag(false),
    lengthScale(1.0),
    lengthScaleGiven(false),
    abstol(1.0e-12), // absolute current error tol.
    reltol(1.0e-4),  // relative current error tol.
    chgtol(1.0e-12), // absolute charge error tol.
    gmin  (1.0e-12), // minimum allowed conductance.
    gmin_orig (1.0e-12), // minimum allowed conductance, final
    gmin_init (1.0e-02), // minimum allowed conductance, initial
    gmin_scalar(1.0e10),
    gmax  (1.0e20),   // maximum allowed conductance.
    testJac_relTol(0.01),   // reltol for num. jacobian diagnostic
    testJac_absTol(1.0e-8), // abstol for num. jacobian diagnostic.
    testJac_SqrtEta(1.0e-8), // dx = numJacSqrtEta * (1.0 + fabs(soln[i]));
    deviceSens_dp(1.0e-8),  // similar to eta, but for numerical device sensitivities
    tnom  (CONSTREFTEMP),
    temp(Util::Param("TEMP", CONSTREFTEMP)),
    matrixSensitivityFlag(false),
    testJacobianFlag (false),
    testJacStartStep(0),
    testJacStopStep(Util::MachineDependentParams::IntMax()),
    testJacWarn (false),
    testJacDeviceNameGiven( false ),
    testJacDeviceName(""),
    voltageLimiterFlag (true),
    b3soiVoltageLimiterFlag (true),
    lambertWFlag (0),
    newMeyerFlag(false),
    icMultiplier (10000.0),
    defaultMaxTimeStep (1.0e99),
    vdsScaleMin(0.3),
    vgstConst(4.5),
    length0(5.0e-6), // used in mosfet "size" homotopy
    width0(200.0e-6), // used in mosfet "size" homotopy
    tox0(6.0e-8), // used in mosfet "size" homotopy
    minRes(0.0),
    minCap(0.0),
    exp_order(100.0),
    zeroResistanceTol(1.0e-100),
    checkForZeroResistance(true),
    debugMinTimestep (0),
    debugMaxTimestep (Util::MachineDependentParams::IntMax()), // later, this should be MAX_INT
    debugMinTime (0),
    debugMaxTime (Util::MachineDependentParams::DoubleMax()),
    verboseLevel (0),
    rc_const(1e-9),
    newABM(0),
    newExcessPhase(EXCESS_PHASE),
    defaultNewExcessPhase(EXCESS_PHASE),
    excessPhaseScalar1 (1.0),
    excessPhaseScalar2 (1.0),
    photocurrentFormulation(Xyce::Device::photoForm::ORIGINAL),
    //photocurrentFormulation(Xyce::Device::photoForm::REDUCED),
    photocurrent_dx_reltol(1.0e-4),
    photocurrent_reltol(1.0e-8),
    photocurrent_abstol(1.0e-8),
    photocurrent_fixed_tau(false),
    photocurrent_FE_predictor(true),
    maskPhotocurrentDelayVars(false),
    disableInitJctFlag(false),
    randomSeed (0),
    tryToCompact (false),
    calculateAllLeadCurrents (false),
    digInitState(3),
    separateLoad(true),
    pwl_BP_off(false)
{
  setSensitivityDebugLevel(0);
  setDeviceDebugLevel(1);
}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::~DeviceOptions
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/01/02
//-----------------------------------------------------------------------------
DeviceOptions::~DeviceOptions ()
{}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::registerOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
bool DeviceOptions::setOptions(const Util::OptionBlock & option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const std::string &tag = (*it).uTag();

    if (tag == "DEFAD")
      defad     = (*it).getImmutableValue<double>();
    else if (tag == "DEFAS")
      defas     = (*it).getImmutableValue<double>();
    else if (tag == "DEFL")
      defl      = (*it).getImmutableValue<double>();
    else if (tag == "DEFW")
      defw      = (*it).getImmutableValue<double>();
    else if (tag == "ABSTOL")
      abstol    = (*it).getImmutableValue<double>();
    else if (tag == "RELTOL")
      reltol    = (*it).getImmutableValue<double>();
    else if (tag == "CHGTOL")
      chgtol    = (*it).getImmutableValue<double>();
    else if (tag == "GMIN")
      gmin      = (*it).getImmutableValue<double>();
    else if (tag == "GMINSCALAR")
      gmin_scalar = (*it).getImmutableValue<double>();
    else if (tag == "GMAX")
      gmax      = (*it).getImmutableValue<double>();
    else if (tag == "TJRELTOL")
      testJac_relTol = (*it).getImmutableValue<double>();
    else if (tag == "TJABSTOL")
      testJac_absTol = (*it).getImmutableValue<double>();
    else if (tag == "TJSQRTETA")
      testJac_SqrtEta = (*it).getImmutableValue<double>();
    else if (tag == "SENSDP")
      deviceSens_dp = (*it).getImmutableValue<double>();
    else if (tag == "TNOM")
      tnom      = (*it).getImmutableValue<double>() + CONSTCtoK;
    else if (tag == "TEMP")
    {
      temp = *it;
      if ( !(*it).isTimeDependent() )
        temp.setVal((*it).getImmutableValue<double>() + CONSTCtoK);
      else
      {
        temp.setVal( (*it).stringValue() );
        temp.setTimeDependent( true );
      }
    }
    else if (tag=="MATRIXSENS")
      matrixSensitivityFlag = static_cast<bool>((*it).getImmutableValue<int>());
    else if (tag == "TESTJAC")
      testJacobianFlag = static_cast<bool>((*it).getImmutableValue<int>());
    else if (tag == "TESTJACSTARTSTEP")
      testJacStartStep = (*it).getImmutableValue<int>();
    else if (tag == "TESTJACSTOPSTEP")
      testJacStopStep = (*it).getImmutableValue<int>();
    else if (tag == "TESTJACWARN")
      testJacWarn = static_cast<bool> ((*it).getImmutableValue<int>());
    else if (tag == "TESTJACDEVICENAME")
    {
      testJacDeviceName = (*it).stringValue();
      testJacDeviceNameGiven = true;
    }
    else if (tag == "VOLTLIM")
      voltageLimiterFlag    = static_cast<bool>((*it).getImmutableValue<int>());
    else if (tag == "B3SOIVOLTLIM")
      b3soiVoltageLimiterFlag = static_cast<bool>((*it).getImmutableValue<int>());
    else if (tag == "LAMBERTW")
      lambertWFlag = static_cast<int>((*it).getImmutableValue<int>());
    else if (tag == "ICFAC" )
      icMultiplier = (*it).getImmutableValue<double>();
    else if (tag == "MAXTIMESTEP" )
      defaultMaxTimeStep = (*it).getImmutableValue<double>();

    else if (tag == "VDSSCALEMIN" )
      vdsScaleMin = (*it).getImmutableValue<double>();
    else if (tag == "VGSTCONST" )
      vgstConst  = (*it).getImmutableValue<double>();
    else if (tag == "LENGTH0" )
      length0 = (*it).getImmutableValue<double>();
    else if (tag == "WIDTH0" )
      width0  = (*it).getImmutableValue<double>();
    else if (tag == "TOX0" )
      tox0    = (*it).getImmutableValue<double>();
    else if (tag == "MINRES" )
      minRes  = (*it).getImmutableValue<double>();
    else if (tag == "MINCAP" )
      minCap  = (*it).getImmutableValue<double>();
    else if (tag == "NEWMEYER" )
      newMeyerFlag = (*it).getImmutableValue<bool>();
    else if (tag == "SENSDEBUGLEVEL")
    {
      setDeviceSensitivityDebugLevel((*it).getImmutableValue<int>());
    }
    else if (tag == "DEBUGLEVEL")
    {
      setDeviceDebugLevel((*it).getImmutableValue<int>());
    }
    else if (tag == "VERBOSELEVEL")
    {
      verboseLevel     = ((*it).getImmutableValue<int>());
    }
    else if (tag == "DEBUGMINTIMESTEP")
    {
      debugMinTimestep = ((*it).getImmutableValue<int>());
    }
    else if (tag == "DEBUGMAXTIMESTEP")
    {
      debugMaxTimestep = ((*it).getImmutableValue<int>());
    }
    else if (tag == "DEBUGMINTIME")
    {
      debugMinTime     = ((*it).getImmutableValue<double>());
    }
    else if (tag == "DEBUGMAXTIME")
    {
      debugMaxTime     = ((*it).getImmutableValue<double>());
    }
    else if (tag == "NEWEXCESSPHASE")
    {
      newExcessPhase = static_cast<bool> ((*it).getImmutableValue<int>());
    }
    else if (tag == "EXCESSPHASESCALAR1")
    {
      excessPhaseScalar1 = ((*it).getImmutableValue<double>());
    }
    else if (tag == "EXCESSPHASESCALAR2")
    {
      excessPhaseScalar2 = ((*it).getImmutableValue<double>());
    }
    else if (tag == "ZERORESISTANCETOL")
    {
      zeroResistanceTol = ((*it).getImmutableValue<double>());
    }
    else if (tag == "CHECKFORZERORESISTANCE")
    {
      checkForZeroResistance = ((*it).getImmutableValue<bool>());
    }
    else if (tag == "RANDOMSEED" )
    {
      randomSeed = ((*it).getImmutableValue<long>());
    }
    else if (tag == "TRYTOCOMPACT")
    {
      tryToCompact = static_cast<bool> ((*it).getImmutableValue<int>());
    }
    else if (tag == "CALCULATEALLLEADCURRENTS")
    {
      calculateAllLeadCurrents = ((*it).getImmutableValue<bool>());
    }
    else if (tag == "RCCONST")
    {
      rc_const = ((*it).getImmutableValue<double>());
    }
    else if (tag ==   "SMOOTHBSRC" )
    {
      newABM = static_cast<bool> ((*it).getImmutableValue<int>());
    }
    else if (tag == "DIGINITSTATE")
    {
      digInitState = ((*it).getImmutableValue<int>());
    }
    else if (tag == "SEPARATELOAD")
    {
      separateLoad = static_cast<bool> ((*it).getImmutableValue<int>());
    }
    else if (tag == "PWLBPOFF")
    {
      pwl_BP_off = static_cast<bool> ((*it).getImmutableValue<int>());
    }
#ifdef Xyce_RAD_MODELS
    else if (tag == "PHOTOCURRENT_FORMULATION")
    {
      if ((*it).isNumeric())
      {
        photocurrentFormulation = ((*it).getImmutableValue<int>());
      }
      else
      {
        ExtendedString p((*it).stringValue()); p.toUpper();
        if (p == "ORIGINAL")
        {
          photocurrentFormulation = Xyce::Device::photoForm::ORIGINAL;
        }
        else if (p == "REDUCED")
        {
          photocurrentFormulation = Xyce::Device::photoForm::REDUCED;
        }
        else if (p == "TVR")
        {
          photocurrentFormulation = Xyce::Device::photoForm::TVR;
        }
        else
        {
           Report::UserError() << "Unrecognized photocurrent formulation " << p;
        }
      }
    }
    else if (tag == "PHOTOCURRENT_DX_RELTOL")
    {
      photocurrent_dx_reltol = ((*it).getImmutableValue<double>());
    }
    else if (tag == "PHOTOCURRENT_RELTOL")
    {
      photocurrent_reltol = ((*it).getImmutableValue<double>());
    }
    else if (tag == "PHOTOCURRENT_ABSTOL")
    {
      photocurrent_abstol = ((*it).getImmutableValue<double>());
    }
    else if (tag == "PHOTOCURRENT_FIXED_TAU")
    {
      photocurrent_fixed_tau = ((*it).getImmutableValue<bool>());
    }
    else if (tag == "PHOTOCURRENT_FE_PREDICTOR")
    {
      photocurrent_FE_predictor = ((*it).getImmutableValue<bool>());
    }
    else if (tag == "PHOTOCURRENT_MASKING")
    {
      maskPhotocurrentDelayVars = ((*it).getImmutableValue<bool>());
    }
#endif
    else if (tag == "ALL_OFF")
    {
      disableInitJctFlag = ((*it).getImmutableValue<bool>());
    }
    else
    {
      Report::UserError0() << tag << " is not a recognized device package option.";
    }
  }

  gmin_orig = gmin;
  gmin_init = gmin*gmin_scalar;  // by default, 10 orders of magnitude larger.

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    dout() << *this << std::endl;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::registerParserOptions
// Purpose       :
// Special Notes : this is mostly to support the HSPICE/ngspice option "scale"
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/29/20
//-----------------------------------------------------------------------------
bool DeviceOptions::setParserOptions(const Util::OptionBlock & option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    const std::string &tag = (*it).uTag();

    if (tag == "MODEL_BINNING")
    {
      modelBinningFlag = static_cast<bool> ((*it).getImmutableValue<int>());
    }
    else if (tag == "SCALE")
    {
      lengthScale = (*it).getImmutableValue<double>();
      lengthScaleGiven = true;
    }
  }
  return true;
}

void
DeviceOptions::populateMetadata(
  IO::PkgOptionsMgr &           options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("DEVICE");

  parameters.insert(Util::ParamMap::value_type("DEFAS", Util::Param("DEFAS", 0.0)));
  parameters.insert(Util::ParamMap::value_type("DEFAD", Util::Param("DEFAD", 0.0)));
  parameters.insert(Util::ParamMap::value_type("DEFL", Util::Param("DEFL", 1.0e-4)));
  parameters.insert(Util::ParamMap::value_type("DEFW", Util::Param("DEFW", 1.0e-4)));
  parameters.insert(Util::ParamMap::value_type("GMIN", Util::Param("GMIN", 1.0e-12)));
  parameters.insert(Util::ParamMap::value_type("GMINSCALAR", Util::Param("GMINSCALAR", 1.0e+10)));
  parameters.insert(Util::ParamMap::value_type("GMAX", Util::Param("GMAX", 1.0e20)));
  parameters.insert(Util::ParamMap::value_type("TEMP", Util::Param("TEMP", 27.0)));
  parameters.insert(Util::ParamMap::value_type("TNOM", Util::Param("TNOM", 27.0)));
  parameters.insert(Util::ParamMap::value_type("SCALESRC", Util::Param("SCALESRC", 0.0)));
  parameters.insert(Util::ParamMap::value_type("MATRIXSENS", Util::Param("MATRIXSENS", 0)));
  parameters.insert(Util::ParamMap::value_type("TESTJAC", Util::Param("TESTJAC", 0)));
  parameters.insert(Util::ParamMap::value_type("TESTJACSTARTSTEP", Util::Param("TESTJACSTARTSTEP", 0)));
  parameters.insert(Util::ParamMap::value_type("TESTJACSTOPSTEP", Util::Param("TESTJACSTOPSTEP", 0)));
  parameters.insert(Util::ParamMap::value_type("TJRELTOL", Util::Param("TJRELTOL", 0.01)));
  parameters.insert(Util::ParamMap::value_type("TJABSTOL", Util::Param("TJABSTOL", 1.0e-8)));
  parameters.insert(Util::ParamMap::value_type("TJSQRTETA", Util::Param("TJSQRTETA", 1.0e-8)));
  parameters.insert(Util::ParamMap::value_type("SENSDP", Util::Param("SENSDP", 1.0e-8)));
  parameters.insert(Util::ParamMap::value_type("TESTJACWARN", Util::Param("TESTJACWARN", 0)));
  parameters.insert(Util::ParamMap::value_type("TESTJACDEVICENAME", Util::Param("TESTJACDEVICENAME", "")));
  parameters.insert(Util::ParamMap::value_type("VOLTLIM", Util::Param("VOLTLIM", 1)));
  parameters.insert(Util::ParamMap::value_type("B3SOIVOLTLIM", Util::Param("B3SOIVOLTLIM", 1)));
  parameters.insert(Util::ParamMap::value_type("ICFAC", Util::Param("ICFAC", 10000.0)));
  parameters.insert(Util::ParamMap::value_type("LAMBERTW", Util::Param("LAMBERTW", 0)));
  parameters.insert(Util::ParamMap::value_type("MAXTIMESTEP", Util::Param("MAXTIMESTEP", 1.0e99)));
  parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 65536)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
  parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 100.0)));
  parameters.insert(Util::ParamMap::value_type("VERBOSELEVEL", Util::Param("VERBOSELEVEL", 0)));
  parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0e-12)));
  parameters.insert(Util::ParamMap::value_type("CHGTOL", Util::Param("CHGTOL", 1.0e-12)));
  parameters.insert(Util::ParamMap::value_type("VDSSCALEMIN", Util::Param("VDSSCALEMIN", 0.3)));
  parameters.insert(Util::ParamMap::value_type("VGSTCONST", Util::Param("VGSTCONST", 4.5)));
  parameters.insert(Util::ParamMap::value_type("LENGTH0", Util::Param("LENGTH0", 5.0e-6)));
  parameters.insert(Util::ParamMap::value_type("WIDTH0", Util::Param("WIDTH0", 200.0e-6)));
  parameters.insert(Util::ParamMap::value_type("TOX0", Util::Param("TOX0", 6.0e-8)));
  parameters.insert(Util::ParamMap::value_type("MINRES", Util::Param("MINRES", 0.0)));
  parameters.insert(Util::ParamMap::value_type("MINCAP", Util::Param("MINCAP", 0.0)));
  parameters.insert(Util::ParamMap::value_type("SENSDEBUGLEVEL", Util::Param("SENSDEBUGLEVEL", 0)));
  parameters.insert(Util::ParamMap::value_type("NEWEXCESSPHASE", Util::Param("NEWEXCESSPHASE", 1)));
  parameters.insert(Util::ParamMap::value_type("EXCESSPHASESCALAR1", Util::Param("EXCESSPHASESCALAR1", 1.0)));
  parameters.insert(Util::ParamMap::value_type("EXCESSPHASESCALAR2", Util::Param("EXCESSPHASESCALAR2", 1.0)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_FORMULATION", Util::Param("PHOTOCURRENT_FORMULATION", 0)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_DX_RELTOL", Util::Param("PHOTOCURRENT_DX_RELTOL", 1.0e-4)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_RELTOL", Util::Param("PHOTOCURRENT_RELTOL", 0)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_ABSTOL", Util::Param("PHOTOCURRENT_ABSTOL", 0)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_FIXED_TAU", Util::Param("PHOTOCURRENT_FIXED_TAU", 0)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_FE_PREDICTOR", Util::Param("PHOTOCURRENT_FE_PREDICTOR", 1)));
  parameters.insert(Util::ParamMap::value_type("PHOTOCURRENT_MASKING", Util::Param("PHOTOCURRENT_MASKING", 0)));
  parameters.insert(Util::ParamMap::value_type("ALL_OFF", Util::Param("ALL_OFF", 0)));
  parameters.insert(Util::ParamMap::value_type("RANDOMSEED", Util::Param("RANDOMSEED", 0)));
  parameters.insert(Util::ParamMap::value_type("TRYTOCOMPACT", Util::Param("TRYTOCOMPACT", false)));
  parameters.insert(Util::ParamMap::value_type("CALCULATEALLLEADCURRENTS", Util::Param("CALCULATEALLLEADCURRENTS", false)));
  parameters.insert(Util::ParamMap::value_type("NEWMEYER", Util::Param("NEWMEYER", false )));
  parameters.insert(Util::ParamMap::value_type("ZERORESISTANCETOL", Util::Param("ZERORESISTANCETOL", 1.0e-100 )));
  parameters.insert(Util::ParamMap::value_type("CHECKFORZERORESISTANCE", Util::Param("CHECKFORZERORESISTANCE", true )));
  parameters.insert(Util::ParamMap::value_type("SMOOTHBSRC", Util::Param("SMOOTHBSRC", 0)));
  parameters.insert(Util::ParamMap::value_type("DIGINITSTATE", Util::Param("DIGINITSTATE", 3 )));
   
  parameters.insert(Util::ParamMap::value_type("RCCONST", Util::Param("RCCONST", 1e-9 )));
  parameters.insert(Util::ParamMap::value_type("SEPARATELOAD", Util::Param("SEPARATELOAD", 1)));
  parameters.insert(Util::ParamMap::value_type("PWLBPOFF", Util::Param("PWLBPOFF", 0)));
}

//-----------------------------------------------------------------------------
// Function      : DeviceOptions::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/23/05
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const DeviceOptions & devOp)
{
  os << "\n\n-----------------------------------------" << std::endl
     << "\tDevice Options:\n"
     << "\t\tdefad                 = " << devOp.defad << "\n"
     << "\t\tdefas                 = " << devOp.defas << "\n"
     << "\t\tdefl                  = " << devOp.defl << "\n"
     << "\t\tdefw                  = " << devOp.defw << "\n"
     << "\t\tabstol                = " << devOp.abstol << "\n"
     << "\t\treltol                = " << devOp.reltol << "\n"
     << "\t\tchgtol                = " << devOp.chgtol << "\n"
     << "\t\tgmin                  = " << devOp.gmin << "\n"
     << "\t\tgmin_orig             = " << devOp.gmin_orig << "\n"
     << "\t\tgmin_init             = " << devOp.gmin_init << "\n"
     << "\t\tgmin_scalar           = " << devOp.gmin_scalar << "\n"
     << "\t\tgmax                  = " << devOp.gmax << "\n"
     << "\t\ttnom                  = " << devOp.tnom << "\n"
     << "\t\ttestJacobianFlag      = " << devOp.testJacobianFlag << "\n"
     << "\t\ttestJacStartStep      = " << devOp.testJacStartStep << "\n"
     << "\t\ttestJacStopStep       = " << devOp.testJacStopStep << "\n"
     << "\t\ttestJacWarn           = " << devOp.testJacWarn     << "\n"
     << "\t\ttestJacDeviceName     = " << devOp.testJacDeviceName << "\n"

     << "\t\ttestJac_relTol        = " << devOp.testJac_relTol << "\n"
     << "\t\ttestJac_absTol        = " << devOp.testJac_absTol << "\n"
     << "\t\ttestJac_SqrtEta       = " << devOp.testJac_SqrtEta << "\n"
     << "\t\tdeviceSens_dp         = " << devOp.deviceSens_dp << "\n"

     << "\t\tvoltageLimiterFlag    = " << devOp.voltageLimiterFlag << "\n"
     << "\t\tb3soiVoltageLimiterFlag    = " << devOp.b3soiVoltageLimiterFlag << "\n"
     << "\t\tlambertWFlag          = " << devOp.lambertWFlag << "\n"
     << "\t\ticMultiplier          = " << devOp.icMultiplier << "\n"
     << "\t\tdefaultMaxTimeStep    = " << devOp.defaultMaxTimeStep << "\n"
     << "\t\tvdsScaleMin           = " << devOp.vdsScaleMin << "\n"
     << "\t\tvgstConst             = " << devOp.vgstConst << "\n"
     << "\t\tlength0               = " << devOp.length0 << "\n"
     << "\t\twidth0                = " << devOp.width0 << "\n"
     << "\t\ttox0                  = " << devOp.tox0 << "\n"
     << "\t\tminres                = " << devOp.minRes << "\n"
     << "\t\tmincap                = " << devOp.minCap << "\n"
     << "\t\tdebugMinTimestep      = " << devOp.debugMinTimestep << "\n"
     << "\t\tdebugMaxTimestep      = " << devOp.debugMaxTimestep << "\n"
     << "\t\tdebugMinTime          = " << devOp.debugMinTime << "\n"
     << "\t\tdebugMaxTime          = " << devOp.debugMaxTime << "\n"
     << "\t\tverboseLevel          = " << devOp.verboseLevel << "\n"
     << "\t\tnewExcessPhase        = " << devOp.newExcessPhase << "\n"
     << "\t\tdefaultNewExcessPhase = " << devOp.defaultNewExcessPhase << "\n"
     << "\t\texcessPhaseScalar1    = " << devOp.excessPhaseScalar1 << "\n"
     << "\t\texcessPhaseScalar2    = " << devOp.excessPhaseScalar2 << "\n"
     << "\t\tnewMeyerFlag    = " << devOp.newMeyerFlag << "\n"
     << "\t\tdigInitState    = " << devOp.digInitState << "\n"
     << "\t\tseparateLoad    = " << devOp.separateLoad << "\n"
     << Xyce::section_divider
     << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

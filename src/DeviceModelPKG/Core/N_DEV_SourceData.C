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
// Purpose        : This file  contains the member functions of the
//                  N_DEV_SourceData class, which is used by the Vsrc
//                  and ISRC devices.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <fstream>

#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>
#include <N_IO_DeviceBlock.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

void
sourceFunctionMetadata(DeviceParamMap &sourceFcnMap)
{
  std::vector<Param> sourceFcnParamList;

  // Source function metadata.
  sourceFcnParamList.resize(7);
  sourceFcnParamList[0].set("V1", 0.0);
  sourceFcnParamList[1].set("V2", 0.0);
  sourceFcnParamList[2].set("TD", 0.0);
  sourceFcnParamList[3].set("TR", 0.0);
  sourceFcnParamList[4].set("TF", 0.0);
  sourceFcnParamList[5].set("PW", 0.0);
  sourceFcnParamList[6].set("PER", 0.0);
  sourceFcnMap[ "PULSE" ] = sourceFcnParamList;

  sourceFcnParamList.resize(7);
  sourceFcnParamList[0].set("VHI", 0.0);
  sourceFcnParamList[1].set("VLO", 0.0);
  sourceFcnParamList[2].set("TD", 0.0);
  sourceFcnParamList[3].set("TR", 0.0);
  sourceFcnParamList[4].set("TF", 0.0);
  sourceFcnParamList[5].set("TSAMPLE", 0.0);
  sourceFcnParamList[6].set("DATA", "");
  sourceFcnMap[ "PAT" ] = sourceFcnParamList;

  sourceFcnParamList.resize(6);
  sourceFcnParamList[0].set("V0", 0.0);
  sourceFcnParamList[1].set("VA", 0.0);
  sourceFcnParamList[2].set("FREQ", 0.0);
  sourceFcnParamList[3].set("TD", 0.0);
  sourceFcnParamList[4].set("THETA", 0.0);
  sourceFcnParamList[5].set("PHASE", 0.0);
  sourceFcnMap[ "SIN" ] = sourceFcnParamList;

  sourceFcnParamList.resize(6);
  sourceFcnParamList[0].set("V1", 0.0);
  sourceFcnParamList[1].set("V2", 0.0);
  sourceFcnParamList[2].set("TD1", 0.0);
  sourceFcnParamList[3].set("TAU1", 0.0);
  sourceFcnParamList[4].set("TD2", 0.0);
  sourceFcnParamList[5].set("TAU2", 0.0);
  sourceFcnMap[ "EXP" ] = sourceFcnParamList;

  sourceFcnParamList.resize(5);
  sourceFcnParamList[0].set("V0", 0.0);
  sourceFcnParamList[1].set("VA", 0.0);
  sourceFcnParamList[2].set("FC", 0.0);
  sourceFcnParamList[3].set("MDI", 0.0);
  sourceFcnParamList[4].set("FS", 0.0);
  sourceFcnMap[ "SFFM" ] = sourceFcnParamList;

  // sourceFcnParamList is empty for source function "PWL" indicating
  // that it needs special handling in NetlistHandler.C
  sourceFcnParamList.clear();
  sourceFcnMap[ "PWL" ] = sourceFcnParamList;

}

//----------------------------------------------------------------------------
// Function       : CircuitMetadata::getSourceFunctionParameters
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : 
// Creation Date  : 
//----------------------------------------------------------------------------
const std::vector<Param> & getSourceFunctionParameters(
    const std::string &sourceFcn, 
    const IO::DeviceBlock & device_block,
    const IO::TokenVector & parsedInputLine )
{
  static DeviceParamMap sourceFcnMap;

  if (sourceFcnMap.empty())
    sourceFunctionMetadata(sourceFcnMap);

  DeviceParamMap::const_iterator it = sourceFcnMap.find(sourceFcn);
  if (it == sourceFcnMap.end())
  {
    Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_)
      << "No such source function " << sourceFcn << " in " << device_block.getInstanceName();
  }

  return (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : getSourceFunctionID
// Purpose       : Return the integer value corresponding to the
//                 given source function. Return NUM_SRC_DATA if the given string
//                 is unknown.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 09/30/2003
//-----------------------------------------------------------------------------
int getSourceFunctionID(const std::string & sourceFcn)
{
  if (sourceFcn == "PULSE")     return _PULSE_DATA;
  else if (sourceFcn == "SIN")  return _SIN_DATA;
  else if (sourceFcn == "EXP")  return _EXP_DATA;
  else if (sourceFcn == "SFFM") return _SFFM_DATA;
  else if (sourceFcn == "PWL")  return _PWL_DATA;
  else if (sourceFcn == "PAT")  return _PAT_DATA;
  else if (sourceFcn == "DC")   return _DC_DATA;
  else if (sourceFcn == "AC")   return _AC_DATA;  //tmei: 05/02
  //else if (sourceFcn == "DISTOF1")   return _DISTOF1_DATA;
  //else if (sourceFcn == "DISTOF2")   return _DISTOF2_DATA;
  else if (sourceFcn == "PORT" || sourceFcn == "Z0"  )  return _PORT_DATA;
  else return _NUM_SRC_DATA;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::SourceData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
SourceData::SourceData(const SolverState & ss1,
                       const DeviceOptions & do1)
  : sourceName_(""),
  typeName_(""),
  defaultParamName_(""),
  useLocalTime_(false),
  localTime_(0.0),
  time(0.0),
  SourceValue(0.0),
  initializeFlag_(false),
  solState_(ss1),
  devOptions_(do1),
  fastTimeScaleFlag_(false),
  realFlag_(true)
{}

//-----------------------------------------------------------------------------
// Function      : SourceData::~SourceData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
SourceData::~SourceData()
{}

//-----------------------------------------------------------------------------
// Function      : SourceData::initializeSource
// Purpose       : Base class initialization function.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool SourceData::initializeSource ()
{
  initializeFlag_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::returnSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
double SourceData::returnSource()
{
  return SourceValue;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::returnSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
std::string SourceData::getSourceTypeName()
{
  return typeName_;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void SourceData::printOutParams()
{}

//-----------------------------------------------------------------------------
// Function      : SourceData::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/01/02
//-----------------------------------------------------------------------------
double SourceData::getMaxTimeStepSize ()
{
  return devOptions_.defaultMaxTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : SourceData::getTime_
// Purpose       :
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/25/04
//-----------------------------------------------------------------------------
double   SourceData::getTime_()
{
  double tmpTime = 0.0;

  if (useLocalTime_ )
    tmpTime = localTime_;
  else
  {
    if (fastTimeScaleFlag_)
      tmpTime = solState_.currFastTime_;
    else
      tmpTime = solState_.currTime_;
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "SourceData::getTime.  time = " << tmpTime << std::endl
      << "SourceData::getTime.  currFastTime = " << solState_.currFastTime_ << std::endl
      << "SourceData::getTime.  currTime = " << solState_.currTime_ << std::endl;
  }

  return tmpTime;
}

// Class SinData

//-----------------------------------------------------------------------------
// Function      : SinData::SinData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//
// From 3f5:  remember to incorporate this later:
//
//-----------------------------------------------------------------------------
SinData::SinData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState   &         ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  V0(0.0),
  VA(0.0),
  FREQ(0.0),
  TD(0.0),
  THETA(0.0),
  PHASE(0.0),
  V0given(false),
  VAgiven(false),
  FREQgiven(false),
  TDgiven(false),
  THETAgiven(false),
  PHASEgiven(false)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V0")    { V0    = iter->getMutableValue<double>(); V0given = iter->given();}
    if (tmpname == "VA")    { VA    = iter->getMutableValue<double>(); VAgiven = iter->given();}
    if (tmpname == "FREQ")  { FREQ  = iter->getMutableValue<double>(); FREQgiven = iter->given();}
    if (tmpname == "TD")    { TD    = iter->getMutableValue<double>(); TDgiven = iter->given();}
    if (tmpname == "THETA") { THETA = iter->getMutableValue<double>(); THETAgiven = iter->given();}
    if (tmpname == "PHASE") { PHASE = iter->getMutableValue<double>(); PHASEgiven = iter->given();}
  }

  if (!(V0given && VAgiven && FREQgiven))
    UserError(device) << "V0, VA and FREQ are required for the SIN source function";

  typeName_ = "SIN";
  defaultParamName_ = "V0";
}

//-----------------------------------------------------------------------------
// Function      : SinData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool SinData::initializeSource ()
{
  // If neccessary, set defaults:
  double tstop = solState_.finalTime_;

  if (!FREQgiven)  FREQ  = 1.0/tstop;
  if (!TDgiven)    TD    = 0.0;
  if (!THETAgiven) THETA = 0.0;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SinData::~SinData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SinData::~SinData()
{}

//-----------------------------------------------------------------------------
// Function      : SinData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void SinData::printOutParams()
{
  Xyce::dout() << "SinData:\n";
  Xyce::dout() << "V0 = "    << V0 << std::endl;
  Xyce::dout() << "VA = "    << VA << std::endl;
  Xyce::dout() << "FREQ = "  << FREQ << std::endl;
  Xyce::dout() << "TD = "    << TD << std::endl;
  Xyce::dout() << "THETA = " << THETA << std::endl;
  Xyce::dout() << "PHASE = " << PHASE << std::endl;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : SinData::updateSource
// Purpose       : Update the sinwave source.
// Special Notes :
//
//   V0    - offset  (V or A)
//   VA    - Amplitude  (V or A)
//   FREQ  - frequency in Hz
//   TD    - delay in seconds
//   THETA - damping factor (Hz).
//   PHASE - phase (degrees)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool SinData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  time = getTime_();

  time -= TD;
  double mpi = M_PI;
  if (time <= 0)
  {
    //SourceValue = V0;
    SourceValue = V0 + VA * sin (2.0*mpi*(PHASE/360)) ;
  }
  else
  {
    // 2PI to convert from hz to radians/sec
    SourceValue = V0 + VA * sin (2.0*mpi*(FREQ*time + PHASE/360)) *
      exp( -(time*THETA));
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : SinData::getParams
// Purpose       : Pass back the sine source params.
// Special Notes : TD and FREQ are interchanged from their normal order
//                 to accomodate the storage layout in the device classes
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SinData::getParams(double *params)
{
  params[0] = V0;
  params[1] = VA;
  params[2] = TD;
  params[3] = FREQ;
  params[4] = THETA;
  params[5] = PHASE;
}

//-----------------------------------------------------------------------------
// Function      : SinData::setParams
// Purpose       : Update the sine source params.
// Special Notes : TD and FREQ are interchanged from their normal order
//                 to accomodate the storage layout in the device classes
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SinData::setParams(double *params)
{
  bool reset=false;
  if (V0 != params[0])
  {
    V0 = params[0];
    reset = true;
  }
  if (VA != params[1])
  {
    VA = params[1];
    reset = true;
  }
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (FREQ != params[3])
  {
    FREQ = params[3];
    reset = true;
  }
  if (THETA != params[4])
  {
    THETA = params[4];
    reset = true;
  }
  if (PHASE != params[5])
  {
    PHASE = params[5];
    reset = true;
  }
  if (reset)
  {
    updateSource();
  }
}

// Class ExpData

//-----------------------------------------------------------------------------
// Function      : ExpData::ExpData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ExpData::ExpData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState   &         ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  V1   (0.0),
  V2   (0.0),
  TD1  (0.0),
  TAU1 (0.0),
  TD2  (0.0),
  TAU2 (0.0),
  V1given (false),
  V2given (false),
  TD1given (false),
  TAU1given (false),
  TD2given (false),
  TAU2given (false)
{
  // Set the user-defined params:
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V1")    { V1    = iter->getMutableValue<double>(); V1given   = iter->given();}
    if (tmpname == "V2")    { V2    = iter->getMutableValue<double>(); V2given   = iter->given();}
    if (tmpname == "TD1")   { TD1   = iter->getMutableValue<double>(); TD1given  = iter->given();}
    if (tmpname == "TAU1")  { TAU1  = iter->getMutableValue<double>(); TAU1given = iter->given();}
    if (tmpname == "TD2")   { TD2   = iter->getMutableValue<double>(); TD2given  = iter->given();}
    if (tmpname == "TAU2")  { TAU2  = iter->getMutableValue<double>(); TAU2given = iter->given();}
  }

  if (!(V1given && V2given))
    UserError(device) << "V1 and V2 are required for the EXP source function";

  typeName_ = "EXP";
  defaultParamName_ = "V1";
}

//-----------------------------------------------------------------------------
// Function      : ExpData::~ExpData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ExpData::~ExpData()
{}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : ExpData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void ExpData::printOutParams()
{
  Xyce::dout() << "ExpData:\n";

  Xyce::dout() << "V1 = " << V1 << std::endl;
  Xyce::dout() << "V2 = " << V2 << std::endl;

  Xyce::dout() << "TD1 = " << TD1 << std::endl;
  Xyce::dout() << "TAU1 = " << TAU1 << std::endl;

  Xyce::dout() << "TD2 = " << TD2 << std::endl;
  Xyce::dout() << "TAU2 = " << TAU2 << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : ExpData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool ExpData::initializeSource ()
{
  // If neccessary, set defaults:
  double tstep = solState_.startingTimeStep_;

  if (!TD1given)  TD1 = 0.0;
  if (!TAU1given) TAU1 = tstep;
  if (!TD2given)  TD2 = TD1 + tstep;
  if (!TAU2given) TAU2 = tstep;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ExpData::updateSource
// Purpose       : Updates an exponential source:
// Special Notes :
//
//    V1    -  Initial value (V or A)
//    V2    -  Pulsed value (V or A).
//    TD1   -  Rise delay time (seconds).
//    TAU1  -  Rise time constant (seconds)
//    TD2   -  Fall delay time (seconds).
//    TAU2  -  Fall time constant (seconds)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool ExpData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  time = getTime_();

  if (time <= TD1)
  {
    SourceValue = V1;
  }
  else if (time <= TD2)
  {
    SourceValue = V1 + (V2-V1)*(1-exp(-(time-TD1)/TAU1));
  }
  else
  {
    SourceValue = V1 + (V2-V1)*(1-exp(-(time-TD1)/TAU1)) +
      (V1-V2)*(1-exp(-(time-TD2)/TAU2)) ;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ExpData::getParams
// Purpose       : Pass back the exponential source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void ExpData::getParams(double *params)
{
  params[0] = V1;
  params[1] = V2;
  params[2] = TD1;
  params[3] = TAU1;
  params[4] = TD2;
  params[5] = TAU2;
}

//-----------------------------------------------------------------------------
// Function      : ExpData::setParams
// Purpose       : Update the exponential source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void ExpData::setParams(double *params)
{
  bool reset=false;
  if (V1 != params[0])
  {
    V1 = params[0];
    reset = true;
  }
  if (V2 != params[1])
  {
    V2 = params[1];
    reset = true;
  }
  if (TD1 != params[2])
  {
    TD1 = params[2];
    reset = true;
  }
  if (TAU1 != params[3])
  {
    TAU1 = params[3];
    reset = true;
  }
  if (TD2 != params[4])
  {
    TD2 = params[4];
    reset = true;
  }
  if (TAU2 != params[5])
  {
    TAU2 = params[5];
    reset = true;
  }
  if (reset)
  {
    updateSource();
  }
}

// Class PulseData

//-----------------------------------------------------------------------------
// Function      : PulseData::PulseData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PulseData::PulseData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState   &         ss1,
  const DeviceOptions &         do1)
: SourceData (ss1,do1),
  V1  (0.0),
  V2  (0.0),
  TD  (0.0),
  TR  (0.0),
  TF  (0.0),
  PW  (0.0),
  PER (0.0),
  V1given (false),
  V2given (false),
  TDgiven (false),
  TRgiven (false),
  TFgiven (false),
  PWgiven (false),
  PERgiven (false)
{
  // Get the user-defined values:
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V1")  { V1    = iter->getMutableValue<double>(); V1given = iter->given();}
    if (tmpname == "V2")  { V2    = iter->getMutableValue<double>(); V2given = iter->given();}
    if (tmpname == "TD")  { TD    = iter->getMutableValue<double>(); TDgiven = iter->given();}
    if (tmpname == "TR")  { TR    = iter->getMutableValue<double>(); TRgiven = iter->given();}
    if (tmpname == "TF")  { TF    = iter->getMutableValue<double>(); TFgiven = iter->given();}
    if (tmpname == "PW")  { PW    = iter->getMutableValue<double>(); PWgiven = iter->given();}
    if (tmpname == "PER") { PER   = iter->getMutableValue<double>(); PERgiven = iter->given();}
  }

  // For HSpice compatibility, at least V1 must be given.  Note: that in
  // other Spicen (e.g, PSpice and ngspice), both V1 and V2 must be given.
  if (!V1given)
    UserError(device) << "V1 is required for the PULSE source function";

  typeName_ = "PULSE";
  defaultParamName_ = "V2";
}

//-----------------------------------------------------------------------------
// Function      : PulseData::~PulseData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PulseData::~PulseData()
{}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : PulseData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void PulseData::printOutParams()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "  PulseData::printOutParams\n";
  Xyce::dout() << "  V1  = "    << V1 << std::endl;
  Xyce::dout() << "  V2  = "    << V2 << std::endl;

  Xyce::dout() << "  TD  = "    << TD << std::endl;
  Xyce::dout() << "  TR  = "    << TR << std::endl;
  Xyce::dout() << "  TF  = "    << TF << std::endl;
  Xyce::dout() << "  PW  = "    << PW << std::endl;
  Xyce::dout() << "  PER = "    << PER << std::endl;
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
bool PulseData::initializeSource ()
{
  // If neccessary, set the defaults:

  double tstep = solState_.startingTimeStep_;
  double tstop = solState_.finalTime_;

  if (!TDgiven)  TD  = 0.0;
  if (!TRgiven)  TR  = tstep;
  if (!TFgiven)  TF  = tstep;
  if (!PWgiven)  PW  = tstop;
  if (!PERgiven) PER = tstop;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool PulseData::updateSource()
{
  //notused:  double tstep = solState_.startingTimeStep;
  //notused:  double tstop = solState_.finalTime_;
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  double basetime = 0;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  PulseData::updateSources\n";
    printOutParams();
  }

  time = getTime_();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  Time = " << time << std::endl;
  }

  time -= TD;

  if (time > PER && PER != 0.0)
  {
    // repeating signal - figure out where we are in period
    basetime = PER * floor(time/PER);
    time -= basetime;
  }

  // This section got ugly because of a nasty roundoff bug.
  // Instead of doing "time > X" you need also check that time
  // is not within bptol of X.
  // So the following translation is used:
  // Instead of:                           we do:
  //  time > X                            time>X && fabs(time-x)>bptol
  //  time <= X                           time<X || fabs(time-x)<bptol

  if (time <= 0 || (time > (TR + PW + TF) &&
        (fabs (time - (TR+PW+TF)) > solState_.bpTol_) ) )
  {
    SourceValue = V1;
  }
  else if ((time > TR && fabs(time-TR) > solState_.bpTol_)
      && (time < (TR + PW) || fabs (time-(TR+PW))<solState_.bpTol_) )
  {
    SourceValue = V2;
  }
  else if (time > 0 &&
      (time < TR || fabs(time-TR) < solState_.bpTol_))
  {
    if (TR != 0.0)
    {
      SourceValue = V1 + (V2 - V1) * (time) / TR;
    }
    else
    {
      SourceValue = V1;
    }
  }
  else
  { // time > (TR + PW) && <= (TR + PW + TF)
    if (TF != 0.0)
    {
      SourceValue = V2 + (V1 - V2) * (time - (TR + PW)) / TF;
    }
    else
    {
      SourceValue = V2; //      SourceValue = 0.5 * (V1 + V2);
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  SourceValue = " << SourceValue << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::getParams
// Purpose       : Pass back the pulse source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void PulseData::getParams(double *params)
{
  params[0] = V1;
  params[1] = V2;
  params[2] = TD;
  params[3] = TR;
  params[4] = TF;
  params[5] = PW;
  params[6] = PER;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::setParams
// Purpose       : Update the pulse source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void PulseData::setParams(double *params)
{
  bool reset=false;
  if (V1 != params[0])
  {
    V1 = params[0];
    reset = true;
  }
  if (V2 != params[1])
  {
    V2 = params[1];
    reset = true;
  }
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (TR != params[3])
  {
    TR = params[3];
    reset = true;
  }
  if (TF != params[4])
  {
    TF = params[4];
    reset = true;
  }
  if (PW != params[5])
  {
    PW = params[5];
    reset = true;
  }
  if (PER != params[6])
  {
    PER = params[6];
    reset = true;
  }
  if (reset)
  {
    updateSource();
  }
}

//-----------------------------------------------------------------------------
// Function      : PulseData::getBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : Like much of this file, this is adapted from spice 3f5.
//                 Some of the stuff in it is a little hokey, and may get
//                 removed or modified later.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/08/01
//-----------------------------------------------------------------------------
bool PulseData::getBreakPoints
(std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  // currPeriodIndex is the integer index representing the period of
  // the current circuit time.
  int currPeriodIndex = 0;

  // subtract out the delay.
  double basetime = 0.0;

  time = getTime_() - TD;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  In PulseData::getBreakPoints\n";
    Xyce::dout() << "  time = " << time << std::endl;
    Xyce::dout() << "  TD   = " << TD  << std::endl;
    Xyce::dout() << "  PER  = " << PER << std::endl;
  }

  // repeating signal - figure out where we are in period
  if(time >= PER)
  {
    if (PER != 0.0)
    {
      currPeriodIndex = (static_cast<int> (floor(time/PER)));
      basetime = PER * (static_cast<double> (currPeriodIndex));
      time -= basetime;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  time            = " << time << std::endl;
    Xyce::dout() << " basetime         = " << basetime << std::endl;
    Xyce::dout() << "  currPeriodIndex = " << currPeriodIndex << std::endl;
    Xyce::dout() << std::endl;
  }

  // now that we know which period this is, push_back all breakpoints
  // in this period and the next.  If we are still in the delay, then
  // just use first two periods.

  // current period:
  breakPointTimes.push_back(basetime+TD);
  breakPointTimes.push_back(basetime+TD+TR);
  breakPointTimes.push_back(basetime+TD+TR+PW);
  breakPointTimes.push_back(basetime+TD+TR+PW+TF);

  if (PER != 0.0)
  {
    breakPointTimes.push_back(basetime+TD+PER);

    // next period:
    breakPointTimes.push_back(basetime+TD+PER+TR);
    breakPointTimes.push_back(basetime+TD+PER+TR+PW);
    breakPointTimes.push_back(basetime+TD+PER+TR+PW+TF);
    breakPointTimes.push_back(basetime+TD+PER+PER);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PulseData::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/25/01
//-----------------------------------------------------------------------------
double PulseData::getMaxTimeStepSize ()
{
  double maxTimeStep = devOptions_.defaultMaxTimeStep;

  // check if we are still in the delay or not.
  time = getTime_();

  if (time < TD) maxTimeStep = (0.1*TD );
  else           maxTimeStep = (0.1*PER);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "\nIn PulseData::getMaxTimeStepSize.  ";
    Xyce::dout() << " maxTimeStep = "<< maxTimeStep;
    Xyce::dout() << "  TD = " << TD << "  PER = " <<PER;
    Xyce::dout() << "  time = "<< time << std::endl;
  }

  return maxTimeStep;
}

// Class PWLinData
//-----------------------------------------------------------------------------
// Function      : debugOutput1
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/22/2023
//-----------------------------------------------------------------------------
void debugOutput1( const std::string funcName, const Param & param)
{
    const std::string & tmpname = param.tag();
    const bool & tmpgiven = param.given();

    if (tmpgiven)
    {
      std::cout << funcName;

      switch (param.getType()) 
      {
        case Xyce::Util::STR:
          std::cout<<" "<<param.uTag()<<" is STR type; value =  "<< param.stringValue()<<std::endl;
          break;
        case Xyce::Util::DBLE:
          std::cout<<" "<<param.uTag()<<" is DBLE type; value =  "<< param.getMutableValue<double>()<<std::endl;
          break;
        case Xyce::Util::EXPR:
          std::cout<<" "<<param.uTag()<<" is EXPR type; value =  "<<param.getValue<Util::Expression>().get_expression()<<std::endl;
          break;
        case Xyce::Util::BOOL:
          std::cout<<" "<<param.uTag()<<" is BOOL type; value =  "<<param.stringValue()<<std::endl;
          break;
        case Xyce::Util::STR_VEC:
          std::cout<<" "<<param.uTag()<<" is STR_VEC type; value =  "<<param.stringValue()<<std::endl;
          break;
        case Xyce::Util::INT_VEC:
          std::cout<<" "<<param.uTag()<<" is INT_VEC type; value =  "<<param.stringValue()<<std::endl;
          break;
        case Xyce::Util::DBLE_VEC:
          std::cout<<" "<<param.uTag()<<" is DBLE_VEC type; value =  "<<param.stringValue()<<std::endl;
          break;
        case Xyce::Util::DBLE_VEC_IND:
          std::cout<<" "<<param.uTag()<<" is DBLE_VEC_IND type; value =  "<<param.stringValue()<<std::endl;
          break;
        case Xyce::Util::COMPOSITE:
          std::cout<<" "<<param.uTag()<<" is COMPOSITE type; value =  "<<param.stringValue()<<std::endl;
          break;
        default:
          std::cout<<" "<<param.uTag()<<" is default type (whatever that is) value = : "<<param.stringValue()<<std::endl;
      }

      std::cout << std::endl;
    }
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::PWLinData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PWLinData::PWLinData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState &           ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  NUM(0),
  REPEAT(false),
  REPEATTIME(0.0),
  TD(0.0),
  loc_(0),
  preComputedBreakpointsDone(false)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();
    const bool & tmpgiven = iter->given();

    if (tmpname == "NUM")        NUM        = iter->getMutableValue<int>();
    if (tmpname == "R" && tmpgiven == true)
    {
      REPEAT = true;
      REPEATTIME = iter->getMutableValue<double>();
    }
    if (tmpname == "TD")         TD         = iter->getMutableValue<double>();

    if ( tmpname == "T" && iter->given() )
    {
      time = iter->getMutableValue<double>();

      if (iter->getType() == Xyce::Util::EXPR)
      {
        Util::Expression expression = iter->getValue<Util::Expression>();
        int index = TVVEC.size();
        timeExprList.push_back(std::pair<int,Util::Expression> (index,expression));
      }

      ++iter;

      // The if statement handles the case where the first time-value pair is not (0,0).
      // It inserts a time-value pair at time=0, with the voltage level of the 
      // first user specified time-voltage pair.  This is compatible with what
      // PSpice and HSpice do.  The number of time-value pairs (NUM) is also 
      // incremented by 1, if that time-voltage pair at time=0 needs to be inserted.
      if ( TVVEC.empty() && time > 0.0 )
      {
        if (iter->getType() == Xyce::Util::EXPR)
        {
          Util::Expression expression = iter->getValue<Util::Expression>();
          int index = TVVEC.size();
          valExprList.push_back(std::pair<int,Util::Expression> (index,expression));
        }

        TVVEC.push_back(std::pair<double,double>(0.0,iter->getMutableValue<double>()));
        NUM++;
      }

      if (iter->getType() == Xyce::Util::EXPR)
      {
        Util::Expression expression = iter->getValue<Util::Expression>();
        int index = TVVEC.size();
        valExprList.push_back(std::pair<int,Util::Expression> (index,expression));
      }
      TVVEC.push_back(std::pair<double,double>(time, iter->getMutableValue<double>()));
    }
  }

  // an empty list of time-voltage pairs will likely fail parsing in 
  // extractSourceFields().  So, this error test may never be used.
  if (NUM == 0)
    UserError(device) << "At least one time-voltage pair must be specified for the PWL source function";

  // only do this test if the REPEAT parameter (R) is true
  if ( REPEAT && (REPEATTIME < 0 || REPEATTIME >=TVVEC[NUM-1].first) )
    UserError(device) << "PWL source repeat value (R) must be >= 0 and < last value in time-voltage list";

  typeName_ = "PWL";
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::PWLinData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
PWLinData::~PWLinData()
{
  timeExprList.clear();
  valExprList.clear();
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void PWLinData::printOutParams()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "  NUM  = "    << NUM << std::endl;
  Xyce::dout() << "  REPEAT  = "    << REPEAT << std::endl;
  Xyce::dout() << "  REPEATTIME  = "    << REPEATTIME << std::endl;
  Xyce::dout() << "  TD  = "    << TD << std::endl;
  Xyce::dout() << "  loc_  = "    << loc_ << std::endl;

  Xyce::dout() << "  Time    Voltage" << std::endl;
  for( int i = 0; i < NUM; ++i )
    Xyce::dout() << " " << TVVEC[i].first << "  " << TVVEC[i].second << std::endl;

  Xyce::dout() << std::endl;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : PWLinData::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
bool PWLinData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  PWLinData::updateSource\n";
    printOutParams();
  }

  time = getTime_();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  Time = " << time << std::endl;
  }

  double time1, time2, voltage1, voltage2;

  if( time >= TD )
  {
    time -= TD;

    if( time <= TVVEC[NUM-1].first )
    {
      for( int i = 0; i < NUM; ++i )
        if( time < TVVEC[i].first ) {loc_ = i;break;}

      if( loc_ == 0 )
      {
        time1 = 0.0;
        voltage1 = 0.0;
      }
      else
      {
        time1 = TVVEC[loc_-1].first;
        voltage1 = TVVEC[loc_-1].second;
      }
      time2 = TVVEC[loc_].first;
      voltage2 = TVVEC[loc_].second;

    }
    else if( !REPEAT )
    {
      // this has the net-effect of setting SourceValue = voltage2
      // when combined with the next if-else block
      time1 = 0.0;
      time2 = 1.0;
      voltage1 = voltage2 = TVVEC[NUM-1].second;
    }
    else
    {
      double looptime = TVVEC[NUM-1].first - REPEATTIME;
      time -= TVVEC[NUM-1].first;
      time -= looptime * floor(time / looptime);
      time += REPEATTIME;

      for( int i = 0; i < NUM; ++i )
      {
        if( time < TVVEC[i].first ) {loc_ = i;break;}
      }
      if (time == REPEATTIME)
      {
        time1 = 0.0;
        time2 = 1.0;
        voltage1 = voltage2 = TVVEC[NUM-1].second;
      }
      else
      {
        if( loc_ == 0 )
        {
          time1 = REPEATTIME;
          voltage1 = TVVEC[NUM-1].second;
        }
        else
        {
          time1 = TVVEC[loc_-1].first;
          voltage1 = TVVEC[loc_-1].second;
        }
        time2 = TVVEC[loc_].first;
        voltage2 = TVVEC[loc_].second;
      }

    }

    if( time1 == time2 )
      SourceValue = voltage2;
    else
    {
      double length = time2 - time1;
      SourceValue = ( time2 - time ) * voltage1 / length;
      SourceValue += ( -time1 + time ) * voltage2 / length;
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
    {
      Xyce::dout() << "time: " << time << std::endl;
      Xyce::dout() << "time1: " << time1 << std::endl;
      Xyce::dout() << "time2: " << time2 << std::endl;
      Xyce::dout() << "voltage1: " << voltage1 << std::endl;
      Xyce::dout() << "voltage2: " << voltage2 << std::endl;
      Xyce::dout() << "Src: " << SourceValue << std::endl;
      Xyce::dout() << Xyce::section_divider << std::endl;
    }
  }
  else
  {
    SourceValue = 0.0;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/08/01
//-----------------------------------------------------------------------------
bool PWLinData::getBreakPoints
(std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  In PWLinData::getBreakPoints\n";
  }

  if (!initializeFlag_) bsuccess = initializeSource ();

  time = solState_.currTime_;

  time -= TD;

  if (devOptions_.pwl_BP_off) { return true; }

  // if it's a repeating signal, figure out what period we are in
  if (REPEAT && time >= TVVEC[NUM - 1].first)
  {
    double loopBaseTime = 0.0;
    double bp_time;

    double looptime = TVVEC[NUM-1].first - REPEATTIME;
    loopBaseTime = looptime * (1.0 + floor((time - TVVEC[NUM - 1].first) / looptime));
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
    {
      Xyce::dout() << "loopBaseTime: " << loopBaseTime << std::endl;
      Xyce::dout() << "floor function: " << floor((time - TVVEC[NUM - 1].first) / (TVVEC[NUM - 1].first - REPEATTIME)) << std::endl;
    }

    // now that we know which period this is, push_back all breakpoints
    // in this period and the next.
    for (int i = 0; i < NUM; ++i)
    {
      bp_time = TVVEC[i].first;
      if (bp_time >= REPEATTIME)
      {
        breakPointTimes.push_back(bp_time + loopBaseTime + TD);
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
        {
          Xyce::dout() << "bp_time + loopBaseTime + TD: " << bp_time + loopBaseTime + TD << std::endl;
        }
      }
    }
  }
  else
  {
    // if this is not periodic, then only need to send the full set of BP once.
    if (!preComputedBreakpointsDone)
    {
      int origSize = breakPointTimes.size();
      breakPointTimes.reserve(origSize+NUM);
      for (int i = 0; i < NUM; ++i)
      {
        double bp_time = TVVEC[i].first;
        breakPointTimes.push_back(bp_time+TD);
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
        {
          Xyce::dout() << "bp_time + TD: " << bp_time + TD << std::endl;
        }
      }
      preComputedBreakpointsDone = true;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::getParams
// Purpose       : Pass back the PWL source params.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/22/2023
//-----------------------------------------------------------------------------
void PWLinData::getParams (double *params)
{
// TD is par2
// R or REPEATTIME is not part of the "par" array.  However, it is right after par6.
// T and V are weird, so not part of double *
  params[2] = TD;
  params[7] = REPEATTIME;
}

//-----------------------------------------------------------------------------
// Function      : PWLinData::getParams
// Purpose       : Set the PWL source params.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/22/2023
//-----------------------------------------------------------------------------
void PWLinData::setParams (double *params)
{
  bool reset=false;
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (REPEATTIME != params[7])
  {
    REPEATTIME = params[7];
    reset = true;
  }

  // now handle the vector of T,V pairs.  
  //
  // ERK. 1/23/2023.  The T,V pairs cannot be updated via the params argument to this 
  // function. The reason is that the params argument, while it is a double* and thus 
  // can represent an array is not conveniently set up to handle arrays of arbitrary 
  // length. It is assuming that we have a fixed set of variables to potentially updated,
  // and they are in a certain order in the N_DEV_Vsrc and/or N_DEV_ISRC classes.  The
  // relevant variables in the N_DEV_Vsrc are double T and double V, which are not arrays 
  // or vectors.    They have historically just served as placeholders.  When the user
  // specifies a PWL list of numbers, you wind up with an arbibrary number of 
  // parameters, all named either "T" or "V".  Unlike some vector variables, these are not
  // given unique names.  Thus, even if there are a lot of them, they only map to a 
  // single double-precision variable.
  //
  // With enough hacking it could be made to work, by change T and V to be true vectors, and 
  // giving them unique names, but I've chosen a different, less cumbersome approach.
  //
  // The approach is this.  If any of the entries in the PWL source are EXPR type, 
  // then store a copy of their expression, as well as their index into the T,V array.  
  // Since these expressions are now stored, they can be re-evaluated whenever this 
  // function is called.  The external dependencies for expressions are handled automatically, 
  // so this update should work.  If the re-evaluated expression produces a different
  // number than the one already stored in the T,V vector, then the T,V vector is updated.
  {
    std::vector< std::pair<int,Util::Expression> >::iterator iter = timeExprList.begin();
    std::vector< std::pair<int,Util::Expression> >::iterator end = timeExprList.end();
    for (;iter!=end;iter++)
    {
      Util::Expression & expr = iter->second;
      double time;
      expr.evaluateFunction(time);
      if( time != TVVEC[iter->first].first )
      {
        TVVEC[iter->first].first = time;
        reset = true;
      }
    }
  }

  {
    std::vector< std::pair<int,Util::Expression> >::iterator iter = valExprList.begin();
    std::vector< std::pair<int,Util::Expression> >::iterator end = valExprList.end();
    for (;iter!=end;iter++)
    {
      Util::Expression & expr = iter->second;
      double val;
      expr.evaluateFunction(val);
      if( val != TVVEC[iter->first].second)
      {
        TVVEC[iter->first].second = val;
      }
    }
  }

  // reset is only true when the timing of the PWL source has changed.
  // changes to the "V" part of the TV pairs doesn' require anything to be
  // re-set.
  //
  // need to reset some things like breakpoints if times have changed.
  if (reset)
  {
    setupBreakPoints() ;
  }
}


// Class PatData

//-----------------------------------------------------------------------------
// Function      : PatData::PatData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 3/19/2019
//-----------------------------------------------------------------------------
PatData::PatData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState &           ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  VHI(0),
  VLO(0),
  TD(0),
  TR(0),
  TF(0),
  TSAMPLE(0),
  DATA(""),
  RB(0),
  R(0),
  PATREPEATTIME(0),
  NUMDATA(0),
  NUMBP(0),
  loc_(0), 
  starttime_(0.0)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    // required parameters
    if (tmpname == "VHI")     { VHI = iter->getMutableValue<double>(); VHIgiven = iter->given(); }
    if (tmpname == "VLO")     { VLO = iter->getMutableValue<double>(); VLOgiven = iter->given(); }
    if (tmpname == "TD")      { TD = iter->getMutableValue<double>(); TDgiven = iter->given(); }
    if (tmpname == "TR")      { TR = iter->getMutableValue<double>(); TRgiven = iter->given(); }
    if (tmpname == "TF")      { TF = iter->getMutableValue<double>(); TFgiven = iter->given(); }
    if (tmpname == "TSAMPLE") { TSAMPLE = iter->getMutableValue<double>(); TSAMPLEgiven = iter->given(); }

    if (tmpname == "DATA")
    {
      DATA = iter->getMutableValue<std::string>();
      DATAgiven = iter->given();
    }

    // optional parameters
    if (tmpname == "RB")      { RB = iter->getMutableValue<int>(); }
    if (tmpname == "R")       { R = iter->getMutableValue<int>(); }
  }

  // check for required parameters
  if (!VHIgiven || !VLOgiven || !TDgiven || !TRgiven || !TFgiven || !TSAMPLEgiven || !DATAgiven)
    UserError(device) << "VHI, VLO, TD, TR, TF, TSAMPLE and DATA parameters are all required for PAT source";

  // negative sample, rise and fall times are not allowed.
  if ( (TSAMPLE <= 0) || (TR <= 0) || (TF <= 0) )
    UserError(device) << "TR, TF and TSAMPLE must all be non-negative for the PAT source function";

  // HSPICE sets R values less than -1 to 0.  R value of -1 means repeat the pattern "forever".
  if (R < -1) { R=0; }

  // HSPICE sets RB values less than 1 to 1.
  if (RB < 1) { RB=1; }
  if ( RB != 1 )
    UserError(device) << "Only RB=1 is supported for the PAT source function";

  // There must be at least one bit
  if ( ((DATA[0]!='B') && (DATA[0]!='b')) || DATA.size() < 2)
  {
    UserError(device) << "Invalid DATA field for the PAT source function";
  }
  else
  {
    NUMDATA = DATA.size()-1;
    size_t zeroCount = std::count(DATA.begin()+1, DATA.end(), '0');
    size_t oneCount  = std::count(DATA.begin()+1, DATA.end(), '1');
    if (zeroCount + oneCount != NUMDATA)
    {
      UserError(device) << "Invalid bit symbol in DATA field for the PAT source function";
    }
    else
    {
      updateTVVEC();
    }
  }

  NUMBP = TVVEC.size();

  typeName_ = "PAT";
}

//-----------------------------------------------------------------------------
// Function      : PatData::PatData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 3/19/19
//-----------------------------------------------------------------------------
PatData::~PatData()
{}

//-----------------------------------------------------------------------------
// Function      : PatData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 3/19/2019
//-----------------------------------------------------------------------------
void PatData::printOutParams()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "  VHI = "     << VHI << std::endl;
  Xyce::dout() << "  VLO = "     << VLO << std::endl;
  Xyce::dout() << "  TD = "      << TD << std::endl;
  Xyce::dout() << "  TR = "      << TR << std::endl;
  Xyce::dout() << "  TF = "      << TF << std::endl;
  Xyce::dout() << "  TSAMPLE = " << TSAMPLE << std::endl;
  Xyce::dout() << "  DATA = "    << DATA << std::endl;
  Xyce::dout() << "  R = "       << R << std::endl;
  Xyce::dout() << "  RB = "      << RB << std::endl;

  Xyce::dout() << "  Time    Voltage" << std::endl;
  for( int i = 0; i < NUMBP; ++i )
    Xyce::dout() << " " << TVVEC[i].first << "  " << TVVEC[i].second << std::endl;

  Xyce::dout() << std::endl;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : PatData::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 3/19/2019
//-----------------------------------------------------------------------------
bool PatData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  PatData::updateSource\n";
    printOutParams();
  }

  time = getTime_();
  time -= TD;  // account for the time delay

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  Time = " << time << std::endl;
  }

  double time1, time2, voltage1, voltage2;

  if( time <= TVVEC[1].first )
  {
    // Use the value in the second breakpoint for the initial value of the
    // first instance of the bit pattern.  The value in the first breakpoint
    // was set to 0.5*(VHI+VLO) in the constructor, if the first and last bits
    // in the the pattern are different.
    SourceValue = TVVEC[1].second;
  }
  else if ( (R >= 0) && (time >= R*NUMDATA*TSAMPLE + TVVEC[NUMBP-2].first) )
  {
    // After the waveform stops repeating, the value in the second-to-last breakpoint
    // is "held" for the rest for the simulation.  (Note: the value in the last
    // breakpoint was set to 0.5*(VHI+VLO) in the constructor, if the first and
    // last bits in the the pattern are different. So, the second-to-last value
    // is used instead of the last value.)
    SourceValue = TVVEC[NUMBP-2].second;
  }
  else
  {
    if (time <= TVVEC[NUMBP-1].first)
    {
      for( int i = 0; i < NUMBP; ++i )
        if( time < TVVEC[i].first ) {loc_ = i;break;}

      time1 = TVVEC[loc_-1].first;
      voltage1 = TVVEC[loc_-1].second;

      time2 = TVVEC[loc_].first;
      voltage2 = TVVEC[loc_].second;
    }
    else
    {
      double looptime = TVVEC[NUMBP-1].first - PATREPEATTIME;
      time -= TVVEC[NUMBP-1].first;
      time -= looptime * floor(time / looptime);
      time += PATREPEATTIME;

      for( int i = 0; i < NUMBP; ++i )
      {
        if( time < TVVEC[i].first ) {loc_ = i;break;}
      }

      if (time == PATREPEATTIME)
      {
        time1 = 0.0;
        time2 = 1.0;
        voltage1 = voltage2 = TVVEC[NUMBP-1].second;
      }
      else
      {
        if( loc_ == 0 )
        {
          time1 = PATREPEATTIME;
          voltage1 = TVVEC[NUMBP-1].second;
        }
        else
        {
          time1 = TVVEC[loc_-1].first;
          voltage1 = TVVEC[loc_-1].second;
        }
        time2 = TVVEC[loc_].first;
        voltage2 = TVVEC[loc_].second;
      }
    }

    if( time1 == time2 )
    {
      SourceValue = voltage2;
    }
    else
    {
      double length = time2 - time1;
      SourceValue = ( time2 - time ) * voltage1 / length;
      SourceValue += ( -time1 + time ) * voltage2 / length;
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
    {
      Xyce::dout() << "loc_: " << loc_ << std::endl;
      Xyce::dout() << "time: " << time << std::endl;
      Xyce::dout() << "time1: " << time1 << std::endl;
      Xyce::dout() << "time2: " << time2 << std::endl;
      Xyce::dout() << "voltage1: " << voltage1 << std::endl;
      Xyce::dout() << "voltage2: " << voltage2 << std::endl;
      Xyce::dout() << "Src: " << SourceValue << std::endl;
      Xyce::dout() << Xyce::section_divider << std::endl;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PatData::updateTVVEC
// Purpose       : Update the TVVEC, based on string in the DATA parameter
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/4/2019
//-----------------------------------------------------------------------------
void PatData::updateTVVEC()
{
  // Make TVVEC from the DATA field.  The first and last bits are special cases, as they
  // are made based on the assumption that the waveforms repeats.  The first/last
  // breakpoints will then be ignored (in updateSource()) for the first/last
  // repeat of the waveform.  (Note: R=0 is essentially repeating the waveform "once".)
  // Also: Note that TD will be added in both updateSource() and getBreakPoints(),
  // so it is not used here. Finally, this approach does not support RB yet.  So,
  // it assumes that RB=1.

  // clear() is needed for the .STEP case
  TVVEC.clear();

  if (DATA[1] == '1')
  {
    if (DATA[NUMDATA] == '0')
    {
      // Add the first breakpoint at the intermediate value if there will be
      // a bit flip at the start of a repeated cycle.
      TVVEC.push_back(std::pair<double,double>(0,0.5*(VHI-VLO)));
    }
    else
    {
      // This introduces extra breakpoints, but makes the "accounting"
      // for repeated waveforms easier in updateSource()
      TVVEC.push_back(std::pair<double,double>(0,VHI));
    }
    TVVEC.push_back(std::pair<double,double>(0.5*TR,VHI));
  }
  else if (DATA[1] == '0')
  {
    if (DATA[NUMDATA] == '0')
    {
      TVVEC.push_back(std::pair<double,double>(0,VLO));
    }
    else
    {
      // Add the first breakpoint at the intermediate value if there will be
      // a bit flip at the start of a repeated cycle.
      TVVEC.push_back(std::pair<double,double>(0,0.5*(VHI-VLO)));
    }
    TVVEC.push_back(std::pair<double,double>(0.5*TF,VLO));
  }

  for (int i=2;i<=NUMDATA;i++)
  {
    // Add breakpoints when data values change
    if (DATA[i] != DATA[i-1])
    {
      if (DATA[i] == '0')
      {
        // falling edge
        TVVEC.push_back(std::pair<double,double>(((i-1)*TSAMPLE) - 0.5*TF,VHI));
        TVVEC.push_back(std::pair<double,double>(((i-1)*TSAMPLE) + 0.5*TF,VLO));
      }
      else if (DATA[i] == '1')
      {
        // rising edge
        TVVEC.push_back(std::pair<double,double>(((i-1)*TSAMPLE) - 0.5*TR,VLO));
        TVVEC.push_back(std::pair<double,double>(((i-1)*TSAMPLE) + 0.5*TR,VHI));
      }
    }
  }

  // Handle the end of the last bit in the pattern, in a way that supports the
  // REPEAT feature.  This approach does not support RB yet.  So, it assumes that RB=1.
  if (DATA[NUMDATA] == '1')
  {
    TVVEC.push_back(std::pair<double,double>((NUMDATA*TSAMPLE)-0.5*TF,VHI));
    if ( DATA[1] != DATA[NUMDATA] )
    {
      TVVEC.push_back(std::pair<double,double>(NUMDATA*TSAMPLE,0.5*(VHI-VLO)));
    }
    else
    {
      TVVEC.push_back(std::pair<double,double>(NUMDATA*TSAMPLE,VHI));
    }
  }
  else if (DATA[NUMDATA] == '0')
  {
    TVVEC.push_back(std::pair<double,double>((NUMDATA*TSAMPLE)-0.5*TR,VLO));
    if ( DATA[1] != DATA[NUMDATA] )
    {
      TVVEC.push_back(std::pair<double,double>(NUMDATA*TSAMPLE,0.5*(VHI-VLO)));
    }
    else
    {
      TVVEC.push_back(std::pair<double,double>(NUMDATA*TSAMPLE,VLO));
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PatData::getParams
// Purpose       : Pass back the pat source params.
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/4/2019
//-----------------------------------------------------------------------------
void PatData::getParams(double *params)
{
  params[0] = VHI;
  params[1] = VLO;
  params[2] = TD;
  params[3] = TR;
  params[4] = TF;
  params[5] = TSAMPLE;
}

//-----------------------------------------------------------------------------
// Function      : PatData::setParams
// Purpose       : Update the pat source params.
// Special Notes :
//
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/4/2019
//-----------------------------------------------------------------------------
void PatData::setParams(double *params)
{
  bool reset=false;
  if (VHI != params[0])
  {
    VHI = params[0];
    reset = true;
  }
  if (VLO != params[1])
  {
    VLO = params[1];
    reset = true;
  }
  if (TD != params[2])
  {
    TD = params[2];
    reset = true;
  }
  if (TR != params[3])
  {
    TR = params[3];
    reset = true;
  }
  if (TF != params[4])
  {
    TF = params[4];
    reset = true;
  }
  if (TSAMPLE != params[5])
  {
    TSAMPLE = params[5];
    reset = true;
  }

  if (reset)
  {
    updateTVVEC();
    updateSource();
  }
}

//-----------------------------------------------------------------------------
// Function      : PatData::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 03/19/2019
//-----------------------------------------------------------------------------
bool PatData::getBreakPoints
(std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
  {
    Xyce::dout() << "  In PatData::getBreakPoints\n";
  }

  if (!initializeFlag_) bsuccess = initializeSource ();

  time = solState_.currTime_;
  time -= TD;

  // if it's a repeating signal, figure out what period we are in
  if ( (R != 0) && (time >= TVVEC[NUMBP - 1].first))
  {
    double loopBaseTime = 0.0;
    double bp_time;

    double looptime = TVVEC[NUMBP-1].first - PATREPEATTIME;
    loopBaseTime = looptime * (1.0 + floor((time - TVVEC[NUMBP - 1].first) / looptime));
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
    {
      Xyce::dout() << "loopBaseTime: " << loopBaseTime << std::endl;
      Xyce::dout() << "floor function: " << floor((time - TVVEC[NUMBP - 1].first) / (TVVEC[NUMBP - 1].first - PATREPEATTIME)) << std::endl;
    }

    for (int i = 0; i < NUMBP; ++i)
    {
      bp_time = TVVEC[i].first;
      if (bp_time >= PATREPEATTIME)
      {
        breakPointTimes.push_back(bp_time + loopBaseTime + TD);
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
        {
          Xyce::dout() << "bp_time + loopBaseTime + TD: " << bp_time + loopBaseTime + TD << std::endl;
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < NUMBP; ++i)
    {
      double bp_time = TVVEC[i].first;
      breakPointTimes.push_back(bp_time+TD);
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && solState_.debugTimeFlag)
      {
        Xyce::dout() << "bp_time + TD: " << bp_time << std::endl;
      }
    }
  }

  return bsuccess;
}

// Class SFFMData

//-----------------------------------------------------------------------------
// Function      : SFFMData::SFFMData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

SFFMData::SFFMData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState   &         ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  V0  (0.0),
  VA  (0.0),
  FC  (0.0),
  MDI (0.0),
  FS  (0.0),
  V0given  (false),
  VAgiven  (false),
  FCgiven  (false),
  MDIgiven (false),
  FSgiven  (false)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "V0")  { V0    = iter->getMutableValue<double>(); V0given = iter->given(); }
    if (tmpname == "VA")  { VA    = iter->getMutableValue<double>(); VAgiven = iter->given(); }
    if (tmpname == "FC")  { FC    = iter->getMutableValue<double>(); FCgiven = iter->given(); }
    if (tmpname == "MDI") { MDI   = iter->getMutableValue<double>(); MDIgiven = iter->given(); }
    if (tmpname == "FS")  { FS    = iter->getMutableValue<double>(); FSgiven = iter->given(); }
  }

  if (!(V0given && VAgiven))
    UserError(device) << "V0 and VA are required for the SFFM source function";

  typeName_ = "SFFM";
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::~SFFMData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SFFMData::~SFFMData()
{}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : SFFMData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

void SFFMData::printOutParams()
{
  Xyce::dout() << "SFFMData:\n";
  Xyce::dout() << "V0 = " << V0 << std::endl;
  Xyce::dout() << "VA = " << VA << std::endl;
  Xyce::dout() << "FC = " << FC << std::endl;
  Xyce::dout() << "MDI = " << MDI << std::endl;
  Xyce::dout() << "FS = " << FS << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::initializeSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------

bool SFFMData::initializeSource ()
{
  // If neccessary, set the defaults:
  double tstop = solState_.finalTime_;

  if (!FCgiven) FC = 1.0/tstop;
  if (!FSgiven) FS = 1.0/tstop;

  initializeFlag_ = true;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::updateSource
// Purpose       :
// Special Notes :
//
//   V0    - offset  (V or A)
//   VA    - Amplitude  (V or A)
//   FC    - carrier frequency
//   MDI   - modulation index
//   FS    - signal frequency
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------

bool SFFMData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  time = getTime_();

  double mpi = M_PI;
  SourceValue = V0 + VA * sin((2 * mpi * FC * time) +
      MDI * sin (2 * mpi * FS * time));

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::getParams
// Purpose       : Pass back the SFFM source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SFFMData::getParams(double *params)
{
  params[0] = V0;
  params[1] = VA;
  params[2] = FC;
  params[3] = MDI;
  params[4] = FS;
}

//-----------------------------------------------------------------------------
// Function      : SFFMData::setParams
// Purpose       : Update the SFFM source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/14/05
//-----------------------------------------------------------------------------
void SFFMData::setParams(double *params)
{
  bool reset=false;
  if (V0 != params[0])
  {
    V0 = params[0];
    reset = true;
  }
  if (VA != params[1])
  {
    VA = params[1];
    reset = true;
  }
  if (FC != params[2])
  {
    FC = params[2];
    reset = true;
  }
  if (MDI != params[3])
  {
    MDI = params[3];
    reset = true;
  }
  if (FS != params[4])
  {
    FS = params[4];
    reset = true;
  }
  if (reset)
    updateSource();
}

// Class ACData

//-----------------------------------------------------------------------------
// Function      : ACData::ACData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
ACData::ACData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState   &         ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  ACMAG(1.0),
  ACPHASE(0.0),
  ACMAGgiven(false),
  ACPHASEgiven(false)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();

    if (tmpname == "ACMAG")    { ACMAG    = iter->getMutableValue<double>(); ACMAGgiven = iter->given();}
    if (tmpname == "ACPHASE") { ACPHASE = iter->getMutableValue<double>(); ACPHASEgiven = iter->given();}
  }

  typeName_ = "AC";
  defaultParamName_ = "ACMAG";

  if (ACMAG == 0.0)
  {
    UserWarning(device) << "AC magnitude is set to 0.0";
  }
}

//-----------------------------------------------------------------------------
// Function      : ACData::~ACData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ACData::~ACData()
{}

//-----------------------------------------------------------------------------
// Function      : ACData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
void ACData::printOutParams()
{
  Xyce::dout() << "ACData:\n";
  Xyce::dout() << "ACMAG = " << ACMAG << std::endl;
  Xyce::dout() << "ACPHASE = " << ACPHASE << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : ACData::updateSource
// Purpose       : Update the sinwave source.
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool ACData::updateSource()
{
  bool bsuccess = true;

  double mpi = M_PI;

  if (!initializeFlag_) bsuccess = initializeSource ();

  if (realFlag_)
  {
    SourceValue = ACMAG * cos(2.0*mpi*ACPHASE/360);
  }
  else
  {
    SourceValue = ACMAG * sin(2.0*mpi*ACPHASE/360);
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  SourceValue = " << SourceValue << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ACData::getParams
// Purpose       : Pass back the AC source params.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date :
//-----------------------------------------------------------------------------
void ACData::getParams(double *params)
{
  params[0] = ACMAG;
  params[1] = ACPHASE;
}

//-----------------------------------------------------------------------------
// Function      : ACData::setParams
// Purpose       : Update the AC source params.
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei
// Creation Date :
//-----------------------------------------------------------------------------
void ACData::setParams(double *params)
{
  bool reset=false;
  if (ACMAG!=  params[0])
  {
    ACMAG = params[0];
    reset = true;

    if (ACMAG == 0.0)
    {
      Report::UserWarning() << "AC magnitude is set to 0.0";
    }
  }
  if ( ACPHASE != params[1])
  {
    ACPHASE = params[1];
    reset = true;
  }
  if (reset)
    updateSource();
}


// Class ConstData

//-----------------------------------------------------------------------------
// Function      : ConstData::ConstData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------

ConstData::ConstData(
  const DeviceEntity &          device,
  const std::vector<Param> &    paramRef,
  const SolverState   &         ss1,
  const DeviceOptions &         do1)
: SourceData(ss1,do1),
  V0(0.0)
{
  std::vector<Param>::const_iterator iter = paramRef.begin();
  std::vector<Param>::const_iterator last = paramRef.end();

  for ( ; iter != last; ++iter)
  {
    const std::string & tmpname = iter->tag();
    if (tmpname == "DCV0")
    {
      V0 = iter->getMutableValue<double>();
    }
  }

  typeName_ = "CONST";
  defaultParamName_ = "DCV0";
  //  SourceValue = V0; // updateSource function is a no-op, essentially.
}

//-----------------------------------------------------------------------------
// Function      : ConstData::~ConstData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------
ConstData::~ConstData()
{}

//-----------------------------------------------------------------------------
// Function      : ConstData::printOutParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------
void ConstData::printOutParams()
{
  Xyce::dout() << "ConstData:\n";
  Xyce::dout() << "V0: " << V0 << std::endl;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : ConstData::updateSource
// Purpose       : Update the const source.
// Special Notes : ERK: this is now a no-op, as the source value is set in
//                 the constructor.  The value for this source never changes,
//                 so there isn't any point is re-setting the same value
//                 over and over.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/5/00
//-----------------------------------------------------------------------------
bool ConstData::updateSource()
{
  bool bsuccess = true;

  if (!initializeFlag_) bsuccess = initializeSource ();

  SourceValue = V0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : ConstData::getParams
// Purpose       : Pass back the const source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/13/05
//-----------------------------------------------------------------------------
void ConstData::getParams(double *params)
{
  params[0] = V0;
}

//-----------------------------------------------------------------------------
// Function      : ConstData::setParams
// Purpose       : Update the const source params.
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/13/05
//-----------------------------------------------------------------------------
void ConstData::setParams(double *params)
{
  if (V0 != params[0])
  {
    V0 = params[0];
    updateSource();
  }
}

//-----------------------------------------------------------------------------
// Function      : findSourceFieldPosition
// Purpose       : Find the position of a given field of an independent source
//                 in parsedLine_. Return the position if the field is found,
//                 otherwise return 0.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/02/2001
//-----------------------------------------------------------------------------
int
findSourceFieldPosition(
    const IO::TokenVector &       parsed_line,
    const std::string &           fieldToFind,
    int                           startPosition )
{
  size_t numFields = parsed_line.size();
  int fieldPosition = 0;

  if ( fieldToFind == "SOURCEFCN" ) // any transient source type.
  {
    ExtendedString field("");
    for ( size_t i = startPosition; i < numFields; ++i )
    {
      // The test here checks the field at position i to see if it
      // corresponds to a source function type registered in metadata.
      field = parsed_line[i].string_;
      field.toUpper();
      if ( getSourceFunctionID(field) != _NUM_SRC_DATA && // ie, not found
          getSourceFunctionID(field) != _AC_DATA  &&   // AC
          getSourceFunctionID(field) != _PORT_DATA ) //   port 
      {
        fieldPosition = i;
        break;
      }
    }
  }
  else
  {
    ExtendedString field("");
    for ( size_t i = startPosition; i < numFields; ++i )
    {
      field = parsed_line[i].string_;
      field.toUpper();
      if ( field == fieldToFind )
      {
        fieldPosition = i;
        break;
      }
    }
  }

  return fieldPosition;
}

//-----------------------------------------------------------------------------
// Function      : extractSourceData
// Purpose       : Extract the device data from parsedLine for independent
//                 source devices.
//
// Special Notes : Independent sources have the general form:
//
//   VXXXXXXX N+ N- <DC<> DC/TRAN VALUE> <AC <ACMAG <ACPHASE>>>
// + <DISTOF1 <F1MAG <F1PHASE>>> <DISTOF2 <F2MAG <F2PHASE>>>
//
// N+ and N- are the positive and negative nodes, respectively.
//
// DC/TRAN is the dc and transient analysis value of the source.
// The DC keyword is optional.
//
// ACMAG is the AC magnitude and ACPHASE is the AC phase. AC keyword and
// values only present if this is an AC source.
//
// DISTOF1 and DISTOF2 are the keywords that specify that the independent
// source has distortion inputs at the frequencies F1 and F2 respectively.
// These are not supported yet.
//
// Any independent source can be assigned a time-dependent value for
// transient analysis. If a source is assigned a time-dependent value,
// the time-zero value is used for dc analysis. There are five independent
// source functions: PULSE, EXP, SIN, PWL, and single-frequency FM (SFFM).
//
// Examples:
// VCC 10 0 DC 6
// VIN 13 2 0.001 AC 1 SIN(0 1 1MEG)
// ISRC 23 21 AC 0.333 45.0 SFFM(0 1 10K 5 1K)
// VMEAS 12 9
// VCARRIER 1 0 DISTOF1 0.1 -90.0
// VMODULATOR 2 0 DISTOF2 0.01
// IIN1 1 5 AC 1 DISTOF1 DISTOF2 0.001
//
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/25/2001
//-----------------------------------------------------------------------------
bool
extractSourceData(
    const IO::TokenVector &   parsedInputLine,
    IO::DeviceBlock &         device_block,
    const std::string &       primaryDeviceParameter)
{
  int numFields = parsedInputLine.size();
  bool tranGiven=false;
  bool acGiven=false;
//  bool dcGiven=true;

  bool dcGiven=false;

  bool portGiven=false;           // check to see if needed      ?? 
  bool z0Given=false;

  // Extract the source device nodes.
  device_block.extractNodes(parsedInputLine, -1, 0);

  // Set up the "primary" parameter.  (for vsrc's this is V0 by default, and is DC)
  Util::Param primaryParameter( "", "" );
  primaryParameter.setTag(primaryDeviceParameter);
  primaryParameter.setVal("0");
  int tranSourceType = _DC_DATA;

  // Set position on line to first field after the nodes.
  int linePosition = device_block.getNumberOfNodes() + 1;

  // Check for DC field.  If present, then reset DC parameter (default parameter)
  // to the value from the netlist.  The DC value is "special" because
  // it doesn't require a DC keyword in front of it.  Also, it is always the
  // first type specified on the instance line.
  //
  // All other src types require a keyword.  So, the handling is different.
  if ( linePosition < numFields )
  {
    ExtendedString field(parsedInputLine[linePosition].string_);
    field.toUpper();
    int sourceID = getSourceFunctionID(field);
    // The field is the DC field, advance past "DC" keyword.
    if ( sourceID == _DC_DATA )
    {
      ++linePosition ;
      if (linePosition < numFields)
      {
        field = parsedInputLine[linePosition].string_;
        field.toUpper();
        sourceID = getSourceFunctionID(field);
        dcGiven=true;

      }
      else
      {
        Report::UserError().at(device_block.getNetlistFilename(), device_block.getLineNumber())
          << "Line for device " << device_block.getInstanceName() << " is missing a value or expression after DC field";
        return false;
      }
    }

    // If the field is NOT related to any recognized non-DC
    // source, then proceed.
    if ( sourceID == _NUM_SRC_DATA )  // ie, type not found
    {
      // Get the DC value.
      if ( linePosition < numFields )
      {
        primaryParameter.setVal( parsedInputLine[linePosition].string_ );
        ++linePosition;

        dcGiven = true;

        // generate error message if DC value is invalid
        if( !( primaryParameter.isNumeric() || primaryParameter.hasExpressionValue() ) )
        {
          Report::UserError().at(device_block.getNetlistFilename(), device_block.getLineNumber())
            << "Invalid DC value \"" << primaryParameter.stringValue() << "\" for device " << device_block.getInstanceName();
          return false;
        }
      }
    }
    else
    {
      // no-op.  primary param was already set to zero.
    }
  }

  // Find the positions of the "AC", "DISTOF1", "DISTOF2", and source
  // function fields in parsedInputLine. Save names and positions of the
  // fields that appear in parsedLine in fieldNames and fieldPositions.
  std::vector<std::string> fieldNames;
  std::vector<int> fieldPositions;

  int fieldPosition = findSourceFieldPosition(parsedInputLine, "AC", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "AC" );
    fieldPositions.push_back( fieldPosition );
    acGiven=true;
    //    tranSourceType = getSourceFunctionID("AC"); //tmei: 05/02
    //    ++linePosition;
  }

  fieldPosition = findSourceFieldPosition(parsedInputLine, "DISTOF1", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "DISTOF1" );
    fieldPositions.push_back( fieldPosition );
  }

  fieldPosition = findSourceFieldPosition(parsedInputLine, "DISTOF2", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "DISTOF2" );
    fieldPositions.push_back( fieldPosition );
  }

  fieldPosition = findSourceFieldPosition(parsedInputLine, "SOURCEFCN", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "SOURCEFCN" );
    fieldPositions.push_back( fieldPosition );
    ExtendedString field( parsedInputLine[fieldPosition].string_ );
    field.toUpper();
    tranSourceType = getSourceFunctionID(field);
    tranGiven=true;
  }
  else
  {
    tranGiven=false;
  }

  // Set the type(s) of the source.
  fieldPosition = findSourceFieldPosition(parsedInputLine, "PORT", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "PORT" );
    fieldPositions.push_back( fieldPosition );
    portGiven=true; 
  }   


  fieldPosition = findSourceFieldPosition(parsedInputLine, "Z0", linePosition );
  if ( fieldPosition > 0 )
  {
    fieldNames.push_back( "Z0" );
    fieldPositions.push_back( fieldPosition );
    z0Given=true;  
  }

  // Set the type(s) of the source.

  Util::Param tranSourceTypeParameter( "TRANSIENTSOURCETYPE", tranSourceType );
  Util::Param acSourceTypeParameter( "ACSOURCETYPE", (int) _AC_DATA);
  Util::Param dcSourceTypeParameter( "DCSOURCETYPE", (int) _DC_DATA);

  // Add the device instance parameters and their default values
  // to instanceParameters and extract instance parameters from
  // parsedInputLine.

  int parameterStartPosition;
  int parameterEndPosition;
  device_block.extractInstanceParameters( parsedInputLine,
      linePosition,
      parameterStartPosition,
      parameterEndPosition,
      "continue");
  fieldNames.push_back( "PARAMS" );
  fieldPositions.push_back( parameterStartPosition );

  // ERK.  Before extractSourceFields can be called, the field vectors
  // need to be in order of ascending field position.
  //
  // This loop below does a quick insertion sort.
  // (see the insort routine from: http://www.yendor.com/programming/sort/)
  //
  // It is quick // as long as the field arrays are less than 15 elements.  I did
  // not use STL sort b/c I am sorting using the contents of fieldPositions,
  // but fieldNames must also be sorted the same way.
  int len = fieldPositions.size();
  for (int i = 1; i < len; i++)
  {
    int j = i;
    int temp = fieldPositions[j];
    std::string tmpName =  fieldNames[j];
    while (j > 0 && (fieldPositions[j-1] > temp))
    {
      fieldPositions[j] = fieldPositions[j-1];
      fieldNames[j] = fieldNames[j-1];
      j--;
    }
    fieldPositions[j] = temp;
    fieldNames[j] = tmpName;
  }

  // Extract the "AC", "DISTOF1", "DISTOF2", and transient
  // source function fields if they appear in parsedInputLine.
  bool result = extractSourceFields(parsedInputLine, device_block, fieldNames, fieldPositions );
  if ( !result )
  {
    return false;
  }

  // Set the DC value.
  Param * parameterPtr;
  if ( primaryParameter.tag() != "" )
  {
    parameterPtr = device_block.findInstanceParameter( primaryParameter );
    parameterPtr->setVal( primaryParameter );
    parameterPtr->setGiven( true );
  }

  // Set the transient source type.
  parameterPtr = device_block.findInstanceParameter( tranSourceTypeParameter );
  parameterPtr->setVal( tranSourceTypeParameter );
  parameterPtr->setGiven( tranGiven );

  // Set the AC source type, if it was given
  parameterPtr = device_block.findInstanceParameter( acSourceTypeParameter );
  parameterPtr->setVal( acSourceTypeParameter );
  parameterPtr->setGiven( acGiven );

  // Set the DC source type, if it was given
  parameterPtr = device_block.findInstanceParameter( dcSourceTypeParameter );
  parameterPtr->setVal( dcSourceTypeParameter );
  parameterPtr->setGiven( dcGiven );

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : extractSourceFields
//
// Purpose       : This function loops through the fields on a source instance
//                 line, processes them, adding defaults when necessary, and
//                 then adds them to the instance parameter database for
//                 this device.
//
// Special Notes : ERK:  this function does at least one really strange thing.
//                 Instead of plugging values into existing metadata, the
//                 parameters for the given source type are added to the
//                 metadata container.  So, for example, if the user has
//                 set V1, and V1 is already in the metadata, this function
//                 does not set the existing V1.  It adds a new one so that
//                 two instances of V1 are present.
//
//                 In practice, this still works, as the devices (vsrc and
//                 isrc) will use the last value set for any parameter,
//                 and the "added" version is the last value.  However,
//                 for a develoer it is confusing.
//
//                 I think the reason for the "adding" is the PWL source,
//                 which can have an arbitrary number of entries.
//
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/03/2001
//-----------------------------------------------------------------------------
bool
extractSourceFields(
    const IO::TokenVector &               parsedInputLine,
    IO::DeviceBlock &                     device_block,
    const std::vector<std::string> &      fieldNames,
    const std::vector<int> &              fieldPositions )
{
  if ( fieldNames[0] == "PARAMS" )
  {
    // No fields to extract.
    return true;
  }

  // Iterate through the fieldNames and extract the corresponding data
  // for the field. The fieldPositions vector gives the position in
  // parsedLine at which the fieldNames appear.
  int i;
  int numFields = fieldNames.size();
  for ( i = 0; i < numFields; ++i )
  {
    //    if ( fieldNames[i] == "AC" || fieldNames[i] == "DISTOF1" ||
    //         fieldNames[i] == "DISTOF2" )
    if ( fieldNames[i] == "AC")
    {
      // The code here inserts instance parameters for magnitude and phase
      // of AC stimulus.
      Param mag("", "");
      Param phase("", "");

      int numTerms = fieldPositions[i+1] - fieldPositions[i] - 1;

      std::string magName(fieldNames[i] + "MAG");
      std::string phaseName(fieldNames[i] + "PHASE");
      mag.set( magName, "1.0" );
      phase.set( phaseName, "0.0" );

      if ( numTerms > 2 )
      {
        device_block.issueUnrecognizedParameterError(fieldNames[i]+" Too Many Terms");
      }

      if ( numTerms >= 1 )
      {
        mag.setVal( parsedInputLine[fieldPositions[i]+1].string_ );
        mag.setGiven( true );
      }

      if ( numTerms == 2 )
      {
        phase.setVal( parsedInputLine[fieldPositions[i]+2].string_ );
        phase.setGiven( true );
      }

      // Add magnitude and phase parameters to instanceParameters.
      device_block.addInstanceParameter( mag );
      device_block.addInstanceParameter( phase );

    }
    else if ( fieldNames[i] == "SOURCEFCN" )
    {
      // Determine the source function.
      const std::string &sourceFunction = parsedInputLine[ fieldPositions[i] ].string_;

      // Check for optional parentheses around source function parameters,
      // and set start and end position of source function parameters on
      // the input line (skip the initial parenthese).
      int sourceFunctionParamStart;
      int sourceFunctionParamEnd;

      // skip forward over nested parentheses, and omit trailing parentheses
      int frontOffset=1;
      int backOffset=1;
      while (fieldPositions[i] + frontOffset < parsedInputLine.size() &&  parsedInputLine[ fieldPositions[i] + frontOffset ].string_ == "(" )
      {
        ++frontOffset;
      }
      while ( parsedInputLine[ fieldPositions[i+1] - backOffset ].string_ == ")" )
      {
        ++backOffset;
      }

      sourceFunctionParamStart = fieldPositions[i] + frontOffset;
      sourceFunctionParamEnd = fieldPositions[i + 1] - backOffset;

      if (!equal_nocase(sourceFunction, "PWL"))
      {
        // Get the source function parameters from metadata, set their
        // values to the given input value and add to the device instance
        // parameter list.

        // Set the number of parameters for this source function on
        // the input line.
        size_t numInputSourceFunctionParams = sourceFunctionParamEnd - sourceFunctionParamStart + 1;

        // const std::vector<Param> &sourceFunctionParameters = metadata_.getSourceFunctionParameters(sourceFunction);
        const std::vector<Param> &sourceFunctionParameters = getSourceFunctionParameters(sourceFunction,
            device_block, parsedInputLine);

        for ( size_t k = 0; k < sourceFunctionParameters.size(); ++k )
        {
          if ( k < numInputSourceFunctionParams )
          {
            Param parameter(sourceFunctionParameters[k].uTag(), 0.0);
            parameter.setVal( parsedInputLine[sourceFunctionParamStart + k].string_ );
            parameter.setGiven( true );
            device_block.addInstanceParameter( parameter );
          }
        }
      }
      else // PWL source
      {
        if ( parsedInputLine.size() <= sourceFunctionParamStart ) 
        {
          // this test handle instance lines of the form V1 1 0 PWL, without any other source data.
          // That instance line would cause a segfault later in this IF-ELSE block.
          Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_)
            << "PWL device missing source parameters: " << device_block.getInstanceName();
          return false;         
        }
        else
        { 
          ExtendedString field(parsedInputLine[sourceFunctionParamStart].string_);
          if (field.toUpper() != "FILE")
          {
            // Xyce expects the "R" parameter to be tagged "REPEATTIME". Find
            // that parameter and rename it.
            //Param* parameterPtr = findInstanceParameter( Param("R", "") );
            //parameterPtr->setTag( "REPEATTIME" );

            // Locate the start of the (time, value) pairs in the PWL source.
            int timeValPairStart = sourceFunctionParamStart;
            while ( parsedInputLine[timeValPairStart].string_ == "R" || parsedInputLine[timeValPairStart].string_ == "TD" )
            {
              timeValPairStart += 3;
              // next if statement accounts for a missing list of time-voltage pairs.  
              // Otherwise the next iteration of the while loop will cause a seg-fault.
              if (timeValPairStart >=  parsedInputLine.size()) {break;}
            }

            // In the PWL source case, sourceFunctionParamEnd could have been
            // miscaluated if the repeat time or time delay appeared with no
            // other tagged parameters (which are required to be at the end of
            // the line. In this case recalculate teh value of
            // sourceFunctionParamEnd.
            if ( sourceFunctionParamEnd < timeValPairStart )
            {
              //sourceFunctionParamEnd = (int) parsedInputLine.size() - 1;
              sourceFunctionParamEnd = static_cast<int> (parsedInputLine.size()) - 1;
              if ( parsedInputLine[sourceFunctionParamEnd].string_ == ")" )
              {
                --sourceFunctionParamEnd;
              }
            }

            // Add (time, value) pairs to the instance parameters. Also,
            // reset the parameter ("NUM") that indicates the number of
            // such pairs and "REPEAT" which is a boolean that indicates
            // that PWL function is periodic.
            Param time( "T", "" );
            Param value( "V", "" );
            int numTimeValuePairs = 0;
            int i = timeValPairStart;
            while ( i < sourceFunctionParamEnd )
            {
              ++numTimeValuePairs;

              time.setVal( parsedInputLine[i].string_ );
              time.setGiven( true );
              device_block.addInstanceParameter( time );
              ++i;

              value.setVal( parsedInputLine[i].string_ );
              value.setGiven( true );
              device_block.addInstanceParameter( value );
              ++i;

              // HSpice compatibility, where time-value pairs are
              // separated with a comma
              if ( i < sourceFunctionParamEnd )
              {
                if ( parsedInputLine[i].string_ == "," ) ++i;
              }
            }

            if ( numTimeValuePairs > 0 )
            {
              Param * parameterPtr = device_block.findInstanceParameter( Param("NUM", "") );
              parameterPtr->setVal( numTimeValuePairs );
              parameterPtr->setGiven( true );
            }
            else
            {
              // no FILE parameter given (that's the next else block)
              // so there should have been some time, value pairs given
              // however, numTimeValueParis == 0.
              Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_)
                << "Could not parse time/value pairs for PWL function in device: " << device_block.getInstanceName();
              return false;
            }
          }
          else
          {
            std::string tvFileNameIN;
            // PWL time value pairs are given in specified file. Get the
            // file name, open and read (time, value) pairs.
            if (parsedInputLine.size() >= sourceFunctionParamStart+2)
            {
              tvFileNameIN = parsedInputLine[sourceFunctionParamStart+1].string_;
            }
            else
            {
              Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_)
                << "Could not parse FILE specification for PWL function in device: " << device_block.getInstanceName();
              return false;
            }

            // If the file name is enclosed in double quotes, strip them.
            std::string tvFileName;
            if (tvFileNameIN[0] == '"' &&
                tvFileNameIN[tvFileNameIN.length()-1] =='"')
            {
              tvFileName = tvFileNameIN.substr(1, tvFileNameIN.length()-2);
            }
            else
            {
              tvFileName = tvFileNameIN;
            }

            // Error out if the user-specified PWL specification file does not exist, cannot
            // be opened, or is a directory name rather than a file name.  See SON Bugs 730 
            // and 785 for more details.
            if ( !(Util::checkIfValidFile(tvFileName)) )
            {
              Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_) 
                << "Could not find file " << tvFileName << " for PWL function in device: " << device_block.getInstanceName();
              return false;
            }

            std::ifstream tvDataIn;
            tvDataIn.open(tvFileName.c_str(), std::ios::in);
            if ( !tvDataIn.is_open() )
            {
              Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_) 
                << "Could not open file " << tvFileName << " for PWL function in device: " << device_block.getInstanceName();
              return false;
            }

            Param time( "T", "" );
            Param value( "V", "" );

            int numTimeValuePairs = 0;
            double timeIn;
            double valueIn;
            while (tvDataIn >> timeIn)
            {
              char ch;
              tvDataIn.get(ch);

              if (tvDataIn >> valueIn)
              {
                ++numTimeValuePairs;

                time.setVal(timeIn);
                time.setGiven( true );
                device_block.addInstanceParameter(time);

                value.setVal(valueIn);
                value.setGiven( true );
                device_block.addInstanceParameter(value);
              }
              else
              {
                Report::UserError().at(device_block.getNetlistFilename(), parsedInputLine[0].lineNumber_)
                  << "Problem reading " << tvFileName << std::endl
                  << "File format must be comma, tab or space separated value. There should be no extra spaces or tabs "
                  << "around the comma if it is used as the separator.";
                return false;
              }
            }

            tvDataIn.close();

            if ( numTimeValuePairs > 0 )
            {
              Param* parameterPtr = device_block.findInstanceParameter( Param("NUM", "") );
              parameterPtr->setVal( numTimeValuePairs );
              parameterPtr->setGiven( true );
            }

            else
            {
              Report::UserError() << "Failed to successfully read " << tvFileName;
              return false;
            }
          }
        }
      }
    }
    else if ( fieldNames[i] == "PORT")
    {
      // The code here inserts instance parameters for port 
      Param port(fieldNames[i], "1");

      int numTerms = fieldPositions[i+1] - fieldPositions[i];

      if ( numTerms > 3 )
      {
        device_block.issueUnrecognizedParameterError(fieldNames[i]+" Too Many Terms");
      }


      int linePosition = fieldPositions[i] + 1;

      if ( parsedInputLine[linePosition].string_ == "=" )
      {
        ++linePosition; // if field is "=", advance to next field
      }

      port.setVal( parsedInputLine[linePosition].string_ );
      port.setGiven( true );

      // Add port parameters to instanceParameters.
      device_block.addInstanceParameter( port );
    }
    else if ( fieldNames[i] == "Z0")
    {
      // The code here inserts instance parameters for Z0 
      Param Z0(fieldNames[i], "50");

      int numTerms = fieldPositions[i+1] - fieldPositions[i];

      if ( numTerms > 3 )
      {
        device_block.issueUnrecognizedParameterError(fieldNames[i]+" Too Many Terms");
      }
      int linePosition = fieldPositions[i] + 1;

      if ( parsedInputLine[linePosition].string_ == "=" )
      {
        ++linePosition; // if field is "=", advance to next field
      }

      Z0.setVal( parsedInputLine[linePosition].string_ );
      Z0.setGiven( true );

      // Add Z0 parameters to instanceParameters.
      device_block.addInstanceParameter( Z0 );
    }  
  }

  return true; // Only get here on success.
}

} // namespace Device
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Source data containers.  Used by the vsrc and isrc
//                  devices.
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
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_SourceData_h
#define Xyce_N_DEV_SourceData_h

// ---------- Standard Includes ----------
#include <vector>
#include <list>

// ----------   Xyce Includes   ----------
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>

#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_Param.h>

enum Src_index {
  _SIN_DATA,
  _EXP_DATA,
  _PULSE_DATA,
  _PWL_DATA,
  _PAT_DATA,
  _SFFM_DATA,
  _DC_DATA,
  _AC_DATA,
  _PORT_DATA,
  _NUM_SRC_DATA
};


namespace Xyce {
namespace Device {

typedef std::map<std::string, std::vector<Param>, LessNoCase> DeviceParamMap;

void sourceFunctionMetadata(DeviceParamMap &map);
int getSourceFunctionID(const std::string & sourceFcn);

//-----------------------------------------------------------------------------
// Class         : SourceData
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 4/24/00
//-----------------------------------------------------------------------------
class SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
  friend class SourceInstance;

public:
  SourceData(const SolverState & ss1, const DeviceOptions & do1);

private:
  SourceData(const SourceData &right);
  SourceData &operator=(const SourceData &right);

public:
  virtual ~SourceData();

  virtual void getSensitivityParams (
      std::vector<std::string> & sensParams,
      std::vector<double> & origVals) {};

  virtual bool getAnalyticSensitivityDevice ( int iparam, double & deriv) { return true; };

  virtual bool initializeSource ();

  virtual bool updateSource() = 0;
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv) = 0;

  virtual void setupBreakPoints() {return;}
  virtual bool getBreakPoints (std::vector<Util::BreakPoint> & breakPointTimes)
  { return true; }

  virtual double getMaxTimeStepSize ();

  virtual void setRealFlag(bool flag) { realFlag_ = true;}

  virtual double period() { return 0.0; }

  double returnSource ();

  std::string getSourceTypeName ();

  virtual void getParams (double *) {}
  virtual void setParams (double *) {}
  virtual void printOutParams();

  bool getFastTimeScaleFlag() const
  {
    return fastTimeScaleFlag_;
  }

  void setFastTimeScaleFlag(bool fastTimeScaleFlag)
  {
    fastTimeScaleFlag_ = fastTimeScaleFlag;
  }


  void setUseLocalTimeFlag(bool useLocalTime )
  {
    useLocalTime_ = useLocalTime;
  }


  void setTime(double time)
  {
    localTime_ = time;
  }

protected:
  double getTime_();

private:
  SourceData ();

protected:
  std::string sourceName_;
  std::string typeName_;
  std::string defaultParamName_;


  bool useLocalTime_;
  double localTime_;

  double time;
  double SourceValue;

  bool initializeFlag_;

  const SolverState & solState_;
  const DeviceOptions & devOptions_;

  bool fastTimeScaleFlag_;

  bool realFlag_;
};

//-----------------------------------------------------------------------------
// Class         : SinData
// Purpose       : This class contains data and functions associated with
//                 sinusoidal independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class SinData : public SourceData
{

public:
  SinData(
    const DeviceEntity &        device,
    const std::vector<Param> &  paramRef,
    const SolverState   &       ss1,
    const DeviceOptions &       do1);

  ~SinData();

  void getSensitivityParams (
      std::vector<std::string> & sensParams,
      std::vector<double> & origVals);

  bool getAnalyticSensitivityDevice (int iparam, double & deriv);

private:
  SinData(const SinData &right);
  SinData &operator=(const SinData &right);

public:
  bool initializeSource ();
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);

  void printOutParams();

  double getMaxTimeStepSize () { return (0.1/FREQ); }

  double period() { return (1.0/FREQ); }

private:
  // Data Members for Class Attributes

  double V0;    // Offset (V or A)
  double VA;    // Amplitude (V or A)
  double FREQ;  // Frequency (Hz)
  double TD;    // Delay (seconds)
  double THETA; // Damping factor (1/seconds)
  double PHASE; // Phase (degrees)

  bool   V0given;
  bool   VAgiven;
  bool   FREQgiven;
  bool   TDgiven;
  bool   THETAgiven;
  bool   PHASEgiven;

  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;
};

//-----------------------------------------------------------------------------
// Class         : ExpData
// Purpose       : This class contains data and functions associated with
//                 exponential independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class ExpData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  ExpData(const DeviceEntity & device, const std::vector<Param> & paramRef,
          const SolverState & ss1,
          const DeviceOptions & do1);

  ~ExpData();

  void getSensitivityParams (
      std::vector<std::string> & sensParams,
      std::vector<double> & origVals);

  bool getAnalyticSensitivityDevice (int iparam, double & deriv);

private:
  ExpData(const ExpData & right);
  ExpData &operator=(const ExpData & right);

public:
  bool initializeSource ();
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);

  void printOutParams ();

private:
  double V1;   // Initial value (V or A)
  double V2;   // Pulsed value (V or A).
  double TD1;  // Rise delay time (seconds).
  double TAU1; // Rise time constant (seconds)
  double TD2;  // Fall delay time (seconds).
  double TAU2; // Fall time constant (seconds)

  bool  V1given;
  bool  V2given;
  bool  TD1given;
  bool  TAU1given;
  bool  TD2given;
  bool  TAU2given;
};


//-----------------------------------------------------------------------------
// Class         : ACData
// Purpose       : This class contains data and functions associated with
//                 AC independent sources.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :
//-----------------------------------------------------------------------------
class ACData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  ACData(const DeviceEntity & device, const std::vector<Param> & paramRef,
         const SolverState   & ss1,
         const DeviceOptions & do1);

  ~ACData();

private:
  ACData(const ACData &right);
  ACData &operator=(const ACData &right);

public:
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);

  void printOutParams();

  void setRealFlag(bool flag) { realFlag_ = flag; }

private:
  double ACMAG;    // Amplitude (V or A)
  double ACPHASE; // Phase (degrees)

  bool ACMAGgiven;
  bool ACPHASEgiven;
};

//-----------------------------------------------------------------------------
// Class         : PulseData
// Purpose       : This class contains data and functions associated with
//                 pulsed independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class PulseData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  PulseData(const DeviceEntity & device, const std::vector<Param> & paramRef,
            const SolverState   & ss1,
            const DeviceOptions & do1);

  ~PulseData();

  void getSensitivityParams (
      std::vector<std::string> & sensParams,
      std::vector<double> & origVals);

  bool getAnalyticSensitivityDevice (int iparam, double & deriv);

private:
  PulseData(const PulseData  & right);
  PulseData &operator=(const PulseData  & right);

public:
  bool initializeSource();
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);
  bool getBreakPoints(std::vector<Util::BreakPoint> & breakPointTimes);

  void printOutParams();

  double getMaxTimeStepSize ();

  double period() { return PER; }

public:
  double V1;  // Initial value  (for a voltage source, units are Volts,
  // For a current source, units are Amps)
  double V2;  // Pulsed value. (Volts or Amps)
  double TD;  // Delay time  (seconds)
  double TR;  // Rise time (seconds)
  double TF;  // Fall Time  (seconds)
  double PW;  // Pulse Width (seconds)
  double PER; // Period (seconds)

  bool V1given;
  bool V2given;
  bool TDgiven;
  bool TRgiven;
  bool TFgiven;
  bool PWgiven;
  bool PERgiven;
};


//-----------------------------------------------------------------------------
// Class         : PWLinData
// Purpose       : This class contains the data and functions associated
//                 with piece-wise linear independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class PWLinData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  PWLinData(const DeviceEntity & device, const std::vector<Param> & paramRef,
            const SolverState   & ss1,
            const DeviceOptions & do1);

  ~PWLinData();

  void getSensitivityParams (
      std::vector<std::string> & sensParams,
      std::vector<double> & origVals);

  bool getAnalyticSensitivityDevice ( int iparam, double & deriv);

private:
  PWLinData(const PWLinData &right);
  PWLinData &operator=(const PWLinData &right);

public:
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  virtual void setupBreakPoints() { preComputedBreakpointsDone = false; }
  bool getBreakPoints( std::vector<Util::BreakPoint> & breakPointTimes);
  void getParams (double *);
  void setParams (double *);
  void printOutParams ();

private:
  // Data Members for Class Attributes
  int NUM; //number of time,voltage pairs
  bool REPEAT; //repeat cycle?
  double REPEATTIME; //start time in cycle for repeat
  double TD; //time delay
  std::vector< std::pair<double,double> > TVVEC;  // Array (time,voltage)

  std::vector< std::pair<int,Util::Expression> >  valExprList;
  std::vector< std::pair<int,Util::Expression> >  timeExprList;

  int loc_; //current location in time vector
  bool preComputedBreakpointsDone;
};

//-----------------------------------------------------------------------------
// Class         : PatData
// Purpose       : This class contains the data and functions associated
//                 with the Pattern independent sources.
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 3/18/19
//-----------------------------------------------------------------------------
class PatData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  PatData(const DeviceEntity & device, const std::vector<Param> & paramRef,
            const SolverState   & ss1,
            const DeviceOptions & do1);

  ~PatData();

private:
  PatData(const PatData &right);
  PatData &operator=(const PatData &right);

public:
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);
  void updateTVVEC();

  bool getBreakPoints( std::vector<Util::BreakPoint> & breakPointTimes);

  void printOutParams ();

private:
  // Data Members for Class Attributes
  double VHI;        // high voltage (or current) value
  double VLO;        // low voltage (or current) value
  double TD;         // time delay
  double TR;         // rise time
  double TF;         // fall time
  double TSAMPLE;    // bit period
  std::string DATA;  // bit pattern (as a string)
  int  RB;           // starting bit (in pattern) when repeating
  int  R;            // number of times pattern repeats
  double PATREPEATTIME; // Time at which the pattern repeats
  std::vector< std::pair<double,double> > TVVEC;  // Array (time,voltage)

  bool VHIgiven;
  bool VLOgiven;
  bool TDgiven;
  bool TRgiven;
  bool TFgiven;
  bool TSAMPLEgiven;
  bool DATAgiven;  

  int NUMDATA;           // number of values in string DATA
  int NUMBP;             // number of breakpoints generated from DATA string
  int loc_; //current location in time vector
  double starttime_; //absolute start time of current cycle
};

//-----------------------------------------------------------------------------
// Class         : SFFMData
// Purpose       : This class contains data and functions associated with
//                 Single-frequency FM independent sources.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class SFFMData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  SFFMData(const DeviceEntity & device, const std::vector<Param> & paramRef,
           const SolverState   & ss1,
           const DeviceOptions & do1);

  ~SFFMData();

  void getSensitivityParams (
      std::vector<std::string> & sensParams,
      std::vector<double> & origVals);

  bool getAnalyticSensitivityDevice (int iparam, double & deriv);

private:
  SFFMData(const SFFMData   & right);
  SFFMData &operator=(const SFFMData   & right);

public:
  bool initializeSource ();
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);

  void printOutParams ();

private:
  double V0;  // Offset. (V or A)
  double VA;  // Amplitude (V or A)
  double FC;  // Carrier frequency (Hz)
  double MDI; // Modulation index
  double FS;  // Signal frequency (Hz)

  bool V0given;
  bool VAgiven;
  bool FCgiven;
  bool MDIgiven;
  bool FSgiven;
};

//-----------------------------------------------------------------------------
// Class         : ConstData
// Purpose       : This class contains data and functions associated with
//                 const DC sources.  It is not yet implemented.
// Special Notes :
// Creator       : Robert Hoekstra
// Creation Date : 10/4/00
//-----------------------------------------------------------------------------
class ConstData : public SourceData
{
  friend class VsrcModel;
  friend class VsrcInstance;
  friend class ISRCModel;
  friend class ISRCInstance;

public:
  ConstData(const DeviceEntity & device, const std::vector<Param> & paramRef,
            const SolverState   & ss1,
            const DeviceOptions & do1);

  ~ConstData();

private:
  ConstData(const ConstData   & right);
  ConstData &operator=(const ConstData   & right);

public:
  virtual bool updateSource();
  virtual bool updateSourceDeriv(const std::string & paramName, double & deriv);
  void getParams (double *);
  void setParams (double *);

  void printOutParams();

private:
  double V0;
};

bool extractSourceData(
  const IO::TokenVector &               parsedInputLine,
  IO::DeviceBlock &                     device_block,
  const std::string &                   primaryDeviceParameter);

bool extractSourceFields(
  const IO::TokenVector &               parsedInputLine,
  IO::DeviceBlock &                     device_block,
  const std::vector<std::string> &      fieldNames,
  const std::vector<int> &              fieldPositions);

} // namespace Device
} // namespace Xyce

#endif

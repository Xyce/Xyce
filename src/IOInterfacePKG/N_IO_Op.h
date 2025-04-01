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
// Purpose        : Provide tools for calculating values for operators like
//                  V(), I(), N(), P() and W().
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 11/15/2013
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Op_h
#define Xyce_N_IO_Op_h

#include <iterator>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_DEV_Op.h>
#include <N_IO_Measure_fwd.h>
#include <N_UTL_fwd.h>

#include <N_IO_OutputMgr.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_UTL_Op.h>
#include <N_UTL_ExpressionData.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : CurrentIndexOp
// Purpose       : Operator for getting the "index" from the currently 
//                 active outputter
// Special Notes : The index is basically the line number of output, starting
//                 at zero for the first line.   Note: This class might be 
//                 removed when the lead currents are removed from the store
//                 vector.  It is still used though (as of Xyce 6.7).  
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class CurrentIndexOp : public Util::Op::Op<CurrentIndexOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  CurrentIndexOp(const std::string &name)
    : Base(name)
  {}

  virtual ~CurrentIndexOp()
  {}

  static complex get(const CurrentIndexOp &op, const Util::Op::OpData &op_data);
};

//-----------------------------------------------------------------------------
// Class         : LeadCurrentIndexOp
// Purpose       : Operator for getting the "index" from the currently active
//                 outputter
// Special Notes : The index is basically the line number of output, starting
//                 at zero for the first line
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class LeadCurrentIndexOp : public Util::Op::Op<LeadCurrentIndexOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  LeadCurrentIndexOp(const std::string &name)
    : Base(name)
  {}

  virtual ~LeadCurrentIndexOp()
  {}

  static complex get(const LeadCurrentIndexOp &op, const Util::Op::OpData &op_data);
};

//-----------------------------------------------------------------------------
// Class         : PowerIndexOp
// Purpose       : Operator for getting the "index" from the currently active
//                 outputter
// Special Notes : The index is basically the line number of output, starting
//                 at zero for the first line
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class PowerIndexOp : public Util::Op::Op<PowerIndexOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  PowerIndexOp(const std::string &name)
    : Base(name)
  {}

  virtual ~PowerIndexOp()
  {}

  static complex get(const PowerIndexOp &op, const Util::Op::OpData &op_data);
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrTimeOp
// Purpose       : Operator for getting the current simulation time being output
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class OutputMgrTimeOp : public Util::Op::Op<OutputMgrTimeOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrTimeOp(const std::string &name, const OutputMgr &output_manager, const double time_scale_factor)
    : Base(name),
      outputMgr_(output_manager),
      timeScaleFactor_(time_scale_factor)
  {}

  virtual ~OutputMgrTimeOp()
  {}

  static complex get(const OutputMgrTimeOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &     outputMgr_;
  const double          timeScaleFactor_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrFrequencyOp
// Purpose       : Operator for getting the current frequency being output
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class OutputMgrFrequencyOp : public Util::Op::Op<OutputMgrFrequencyOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrFrequencyOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrFrequencyOp()
  {}

  static complex get(const OutputMgrFrequencyOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrOutputNoiseOp
// Purpose       : Operator for output noise
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/29/2015
//-----------------------------------------------------------------------------
class OutputMgrOutputNoiseOp : public Util::Op::Op<OutputMgrOutputNoiseOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrOutputNoiseOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrOutputNoiseOp()
  {}

  static complex get(const OutputMgrOutputNoiseOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrInputNoiseOp
// Purpose       : Operator for input noise
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 1/29/2015
//-----------------------------------------------------------------------------
class OutputMgrInputNoiseOp : public Util::Op::Op<OutputMgrInputNoiseOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrInputNoiseOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrInputNoiseOp()
  {}

  static complex get(const OutputMgrInputNoiseOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrOutputNoiseContOp
// Purpose       : Operator for individual output noise contributions from
//               : specific devices
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 11/20/2017
//-----------------------------------------------------------------------------
class OutputMgrOutputNoiseContOp : public Util::Op::Op<OutputMgrOutputNoiseContOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
 OutputMgrOutputNoiseContOp(const std::string &name, int devIndex, const std::vector<int> typeIndex, const OutputMgr &output_manager)
    : Base(name),
      devIndex_(devIndex),
      typeIndex_(typeIndex),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrOutputNoiseContOp()
  {}

  static complex get(const OutputMgrOutputNoiseContOp &op, const Util::Op::OpData &op_data);

  const int               devIndex_;  // index of the requested device in the noiseDataVec_ of the NOISE object
  const std::vector<int>  typeIndex_; // vector of indices of the noise type within a given device's vector in the noiseDataVec_
  const OutputMgr &       outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrInputNoiseContOp
// Purpose       : Operator for individual input noise contributions from
//               : specific devices
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 11/20/2017
//-----------------------------------------------------------------------------
class OutputMgrInputNoiseContOp : public Util::Op::Op<OutputMgrInputNoiseContOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
 OutputMgrInputNoiseContOp(const std::string &name, int devIndex, const std::vector<int> typeIndex, const OutputMgr &output_manager)
    : Base(name),
      devIndex_(devIndex),
      typeIndex_(typeIndex),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrInputNoiseContOp()
  {}

  static complex get(const OutputMgrInputNoiseContOp &op, const Util::Op::OpData &op_data);

  int                     devIndex_;  // index of the noise source in the noiseDataVec_ of the NOISE object
  const std::vector<int>  typeIndex_; // vector of the indices of the noise type within a given device's vector in the noiseDataVec_
  const OutputMgr &       outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrTemperatureOp
// Purpose       : Operator for getting the current simulation temperature
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class OutputMgrTemperatureOp : public Util::Op::Op<OutputMgrTemperatureOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrTemperatureOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrTemperatureOp()
  {}

  static complex get(const OutputMgrTemperatureOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrStepSweepOp
// Purpose       : Operator for getting the current value of the step 
//                 parameter being swept
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class OutputMgrStepSweepOp : public Util::Op::Op<OutputMgrStepSweepOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrStepSweepOp(const std::string &name, const OutputMgr &output_manager, int index)
    :Base(name),
     index_(index),
     outputMgr_(output_manager)
  {}

  virtual ~OutputMgrStepSweepOp()
  {}

  static complex get(const OutputMgrStepSweepOp &op, const Util::Op::OpData &op_data);

  const int           index_;
  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrDCSweepOp
// Purpose       : Operator for getting the current value of the DC 
//                 voltage being swept
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class OutputMgrDCSweepOp : public Util::Op::Op<OutputMgrDCSweepOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrDCSweepOp(const std::string &name, const OutputMgr &output_manager, int index)
    : Base(name),
      index_(index),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrDCSweepOp()
  {}

  static complex get(const OutputMgrDCSweepOp &op, const Util::Op::OpData &op_data);

  const int           index_;
  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgrDCSweepOpCurrentValueOp
// Purpose       : Operator for getting the current value of the DC 
//                 voltage being output
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class OutputMgrDCSweepCurrentValueOp : public Util::Op::Op<OutputMgrDCSweepCurrentValueOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  OutputMgrDCSweepCurrentValueOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~OutputMgrDCSweepCurrentValueOp()
  {}

  static complex get(const OutputMgrDCSweepCurrentValueOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &   outputMgr_;
};

//-----------------------------------------------------------------------------
// Class         : StepNumOp
// Purpose       : Operator for getting the current value of the Step Number
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 8/19/2019
//-----------------------------------------------------------------------------
class StepNumOp : public Util::Op::Op<StepNumOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  StepNumOp(const std::string &name, const OutputMgr &output_manager)
    : Base(name),
      outputMgr_(output_manager)
  {}

  virtual ~StepNumOp()
  {}

  static complex get(const StepNumOp &op, const Util::Op::OpData &op_data);

  const OutputMgr &  outputMgr_;
};


//-----------------------------------------------------------------------------
// Class         : SolutionOp
// Purpose       : Operator for getting the current value of a solution 
//                 vector element
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SolutionOp : public Util::Op::Op<SolutionOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  SolutionOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionOp()
  {}

  static complex get(const SolutionOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : SolutionRealOp
// Purpose       : Operator for getting a solution variable and then 
//                 computing the real part
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SolutionRealOp : public Util::Op::Op<SolutionRealOp, Util::Op::ReduceSum>
{
public:
  SolutionRealOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}


  virtual ~SolutionRealOp()
  {}

  static complex get(const SolutionRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : SolutionImaginaryOp
// Purpose       : Operator for getting a solution variable and then 
//                 computing the imaginary part
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SolutionImaginaryOp : public Util::Op::Op<SolutionImaginaryOp, Util::Op::ReduceSum>
{
public:
  SolutionImaginaryOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionImaginaryOp()
  {}

  static complex get(const SolutionImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : SolutionMagnitudeOp
// Purpose       : Operator for getting a solution variable and then computing
//                 the magnitude
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SolutionMagnitudeOp : public Util::Op::Op<SolutionMagnitudeOp, Util::Op::ReduceSum>
{
public:
  SolutionMagnitudeOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionMagnitudeOp()
  {}

  static complex get(const SolutionMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : SolutionPhaseDegOp
// Purpose       : Operator for getting a solution vector element, and then
//                 computing the phase (in degrees)
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
class SolutionPhaseDegOp : public Util::Op::Op<SolutionPhaseDegOp, Util::Op::ReduceSum>
{
public:
  SolutionPhaseDegOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionPhaseDegOp()
  {}

  static complex get(const SolutionPhaseDegOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : SolutionPhaseRadOp
// Purpose       : Operator for getting a solution vector element, and then
//                 computing the phase (in radians)
// Special Notes :
// Creator       : Dave Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SolutionPhaseRadOp : public Util::Op::Op<SolutionPhaseRadOp, Util::Op::ReduceSum>
{
public:
  SolutionPhaseRadOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionPhaseRadOp()
  {}

  static complex get(const SolutionPhaseRadOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : SolutionDecibelOp
// Purpose       : Operator for getting a solution vector element, and then
//                 computing the magnitude (in dB)
// Special Notes : 
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SolutionDecibelsOp : public Util::Op::Op<SolutionDecibelsOp, Util::Op::ReduceSum>
{
public:
  SolutionDecibelsOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SolutionDecibelsOp()
  {}

  static complex get(const SolutionDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferenceOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the difference between them.  This is needed by 
//                 constructs like V(A,B).
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class VoltageDifferenceOp : public Util::Op::Op<VoltageDifferenceOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  VoltageDifferenceOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceOp()
  {}

  static complex get(const VoltageDifferenceOp &op, const Util::Op::OpData &op_data);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferenceRealOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the difference between their real parts.  This 
//                 is needed by constructs like VR(A,B)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class VoltageDifferenceRealOp : public Util::Op::Op<VoltageDifferenceRealOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceRealOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceRealOp()
  {}

  static complex get(const VoltageDifferenceRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferenceImaginaryOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the difference between their imaginary parts. 
//                 This is needed by constructs like VI(A,B)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class VoltageDifferenceImaginaryOp : public Util::Op::Op<VoltageDifferenceImaginaryOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceImaginaryOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceImaginaryOp()
  {}

  static complex get(const VoltageDifferenceImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferenceMagnitudeOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the magnitude of their difference. 
//                 This is needed by constructs like VM(A,B).
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class VoltageDifferenceMagnitudeOp : public Util::Op::Op<VoltageDifferenceMagnitudeOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceMagnitudeOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceMagnitudeOp()
  {}

  static complex get(const VoltageDifferenceMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferencePhaseDegOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the phase of their difference, in degrees.
//                 This is needed by constructs like VP(A,B).
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
class VoltageDifferencePhaseDegOp : public Util::Op::Op<VoltageDifferencePhaseDegOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferencePhaseDegOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferencePhaseDegOp()
  {}

  static complex get(const VoltageDifferencePhaseDegOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferencePhaseRadOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the phase of their difference, in radians.
//                 This is needed by constructs like VP(A,B).
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class VoltageDifferencePhaseRadOp : public Util::Op::Op<VoltageDifferencePhaseRadOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferencePhaseRadOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferencePhaseRadOp()
  {}

  static complex get(const VoltageDifferencePhaseRadOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : VoltageDifferenceDecibelsOp
// Purpose       : Operator for getting two solution variables, and then
//                 computing the magnitude (in decibels) of their difference. 
//                 This is needed by constructs like VDB(A,B).
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class VoltageDifferenceDecibelsOp : public Util::Op::Op<VoltageDifferenceDecibelsOp, Util::Op::ReduceSum>
{
public:
  VoltageDifferenceDecibelsOp(const std::string &name, int index1, int index2)
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~VoltageDifferenceDecibelsOp()
  {}

  static complex get(const VoltageDifferenceDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsOp
// Purpose       : Operator for getting S-parameter, Y-parameter and
//                 Z-parameter data.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
class RFparamsOp : public Util::Op::Op<RFparamsOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  RFparamsOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsOp()
  {}

  static complex get(const RFparamsOp &op, const Util::Op::OpData &op_data);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsRealOp
// Purpose       : Operator for getting the real part of S-parameter,
//                 Y-parameter and Z-parameter data.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
class RFparamsRealOp : public Util::Op::Op<RFparamsRealOp, Util::Op::ReduceNone>
{
public:
  RFparamsRealOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsRealOp()
  {}

  static complex get(const RFparamsRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsImaginaryOp
// Purpose       : Operator for getting the imaginary part of S-parameter,
//                 Y-parameter and Z-parameter data.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
class RFparamsImaginaryOp : public Util::Op::Op<RFparamsImaginaryOp, Util::Op::ReduceNone>
{
public:
  RFparamsImaginaryOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsImaginaryOp()
  {}

  static complex get(const RFparamsImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsMagnitudeOp
// Purpose       : Operator for getting the magnitude of S-parameter,
//                 Y-parameter and Z-parameter data.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
class RFparamsMagnitudeOp : public Util::Op::Op<RFparamsMagnitudeOp, Util::Op::ReduceNone>
{
public:
  RFparamsMagnitudeOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsMagnitudeOp()
  {}

  static complex get(const RFparamsMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsPhaseDegOp
// Purpose       : Operator for getting the phase of S-parameter,
//                 Y-parameter and Z-parameter data (in degrees).
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
class RFparamsPhaseDegOp : public Util::Op::Op<RFparamsPhaseDegOp, Util::Op::ReduceNone>
{
public:
  RFparamsPhaseDegOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsPhaseDegOp()
  {}

  static complex get(const RFparamsPhaseDegOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsPhaseRadOp
// Purpose       : Operator for getting a solution vector element, and then
//                 computing the phase (in radians)
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/2/2019
//-----------------------------------------------------------------------------
class RFparamsPhaseRadOp : public Util::Op::Op<RFparamsPhaseRadOp, Util::Op::ReduceNone>
{
public:
  RFparamsPhaseRadOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsPhaseRadOp()
  {}

  static complex get(const RFparamsPhaseRadOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : RFparamsDecibelsOp
// Purpose       : Operator for getting the magnitude (in dB) of S-parameter,
//                 Y-parameter and Z-parameter data.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
class RFparamsDecibelsOp : public Util::Op::Op<RFparamsDecibelsOp, Util::Op::ReduceNone>
{
public:
  RFparamsDecibelsOp(const std::string &name, const std::string &type, int index1, int index2)
    : Base(name),
      type_(type),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~RFparamsDecibelsOp()
  {}

  static complex get(const RFparamsDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const std::string   type_;
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : StateOp
// Purpose       : Operator for getting a value out of the state vector.
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StateOp : public Util::Op::Op<StateOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  StateOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StateOp()
  {}

  static complex get(const StateOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : StoreOp
// Purpose       : Operator for getting a value out of the Store vector.
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StoreOp : public Util::Op::Op<StoreOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  StoreOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreOp()
  {}

  static complex get(const StoreOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : StoreRealOp
// Purpose       : Operator for getting a store variable, and then computing 
//                 the real part of that variable
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StoreRealOp : public Util::Op::Op<StoreRealOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  StoreRealOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}


  virtual ~StoreRealOp()
  {}

  static complex get(const StoreRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : StoreImaginaryOp
// Purpose       : Operator for getting a store variable, and then computing 
//                 the imaginary part of that variable
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StoreImaginaryOp : public Util::Op::Op<StoreImaginaryOp, Util::Op::ReduceSum>
{
public:
  StoreImaginaryOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreImaginaryOp()
  {}

  static complex get(const StoreImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : StoreMagnitudeOp
// Purpose       : Operator for getting a store variable, and then computing 
//                 the magnitude of that variable
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StoreMagnitudeOp : public Util::Op::Op<StoreMagnitudeOp, Util::Op::ReduceSum>
{
public:
  StoreMagnitudeOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreMagnitudeOp()
  {}

  static complex get(const StoreMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : StorePhaseDegOp
// Purpose       : Operator for getting a store variable, and then computing
//                 the phase of that variable in degrees
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
class StorePhaseDegOp : public Util::Op::Op<StorePhaseDegOp, Util::Op::ReduceSum>
{
public:
  StorePhaseDegOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StorePhaseDegOp()
  {}

  static complex get(const StorePhaseDegOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};
//-----------------------------------------------------------------------------
// Class         : StorePhaseRadOp
// Purpose       : Operator for getting a store variable, and then computing
//                 the phase of that variable in radians
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StorePhaseRadOp : public Util::Op::Op<StorePhaseRadOp, Util::Op::ReduceSum>
{
public:
  StorePhaseRadOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StorePhaseRadOp()
  {}

  static complex get(const StorePhaseRadOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : StoreDecibelsOp
// Purpose       : Operator for getting a store variable, and then computing 
//                 the magnitude in dB of that variable
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class StoreDecibelsOp : public Util::Op::Op<StoreDecibelsOp, Util::Op::ReduceSum>
{
public:
  StoreDecibelsOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~StoreDecibelsOp()
  {}

  static complex get(const StoreDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataCurrentOp
// Purpose       : Operator for getting a value out of the lead current vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class BranchDataCurrentOp : public Util::Op::Op<BranchDataCurrentOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataCurrentOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataCurrentOp()
  {}

  static complex get(const BranchDataCurrentOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataCurrentRealOp
// Purpose       : Operator for getting a value out of the lead current vector
//                 and then computing its real part
// Special Notes :
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
class BranchDataCurrentRealOp : public Util::Op::Op<BranchDataCurrentRealOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataCurrentRealOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}


  virtual ~BranchDataCurrentRealOp()
  {}

  static complex get(const BranchDataCurrentRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataCurrentImaginaryOp
// Purpose       : Operator for getting a value out of the lead current vector
//                 and then computing its imaginary part
// Special Notes :
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
class BranchDataCurrentImaginaryOp : public Util::Op::Op<BranchDataCurrentImaginaryOp, Util::Op::ReduceSum>
{
public:
  BranchDataCurrentImaginaryOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataCurrentImaginaryOp()
  {}

  static complex get(const BranchDataCurrentImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataCurrentMagnitudeOp
// Purpose       : Operator for getting a value out of the lead current vector
//                 and then computing its magnitude
// Special Notes :
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
class BranchDataCurrentMagnitudeOp : public Util::Op::Op<BranchDataCurrentMagnitudeOp, Util::Op::ReduceSum>
{
public:
  BranchDataCurrentMagnitudeOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataCurrentMagnitudeOp()
  {}

  static complex get(const BranchDataCurrentMagnitudeOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentPhaseDegOp
// Purpose       : Operator for getting a value out of the lead current vector
//                 and then computing its phase, in degrees
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
class BranchDataCurrentPhaseDegOp : public Util::Op::Op<BranchDataCurrentPhaseDegOp, Util::Op::ReduceSum>
{
public:
  BranchDataCurrentPhaseDegOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataCurrentPhaseDegOp()
  {}

  static complex get(const BranchDataCurrentPhaseDegOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentPhaseRadOp
// Purpose       : Operator for getting a value out of the lead current vector
//                 and then computing its phase in radians
// Special Notes :
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
class BranchDataCurrentPhaseRadOp : public Util::Op::Op<BranchDataCurrentPhaseRadOp, Util::Op::ReduceSum>
{
public:
  BranchDataCurrentPhaseRadOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataCurrentPhaseRadOp()
  {}

  static complex get(const BranchDataCurrentPhaseRadOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

class BranchDataCurrentDecibelsOp : public Util::Op::Op<BranchDataCurrentDecibelsOp, Util::Op::ReduceSum>
{
public:
  BranchDataCurrentDecibelsOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataCurrentDecibelsOp()
  {}

  static complex get(const BranchDataCurrentDecibelsOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Function      : BranchDataVoltageOp
// Purpose       : Operator for getting a value out of the branch data 
//                 junction voltage vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class BranchDataVoltageOp : public Util::Op::Op<BranchDataVoltageOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataVoltageOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataVoltageOp()
  {}

  static complex get(const BranchDataVoltageOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);
  
  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataPosNegPowerOp
// Purpose       : Operator for power calculation using I*V
// Special Notes :
// Creator       : Richard Schiek, Sandia
// Creation Date : 05/05/2016
//-----------------------------------------------------------------------------
class BranchDataPosNegPowerOp : public Util::Op::Op<BranchDataPosNegPowerOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataPosNegPowerOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~BranchDataPosNegPowerOp()
  {}

  static complex get(const BranchDataPosNegPowerOp &op, const Util::Op::OpData &op_data);
  
  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataBJTPowerOp
// Purpose       : Operator for power computation for BJTs
// Special Notes :
// Creator       : Richard Schiek, Sandia
// Creation Date : 05/05/2016
//-----------------------------------------------------------------------------
class BranchDataBJTPowerOp : public Util::Op::Op<BranchDataBJTPowerOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataBJTPowerOp(const std::string &name, int indexB, int indexC, int indexE, int indexS )
    : Base(name),
      indexB_(indexB),
      indexC_(indexC),
      indexE_(indexE),
      indexS_(indexS)
  {}

  virtual ~BranchDataBJTPowerOp()
  {}

  static complex get(const BranchDataBJTPowerOp &op, const Util::Op::OpData &op_data);

  const int           indexB_;
  const int           indexC_;
  const int           indexE_;
  const int           indexS_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataMOSFETPowerOp
// Purpose       : Operator for power computation for MOSFETs
// Special Notes :
// Creator       : Pete Sholander, Sandia
// Creation Date : 8/18/2016
//-----------------------------------------------------------------------------
class BranchDataMOSFETPowerOp : public Util::Op::Op<BranchDataMOSFETPowerOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataMOSFETPowerOp(const std::string &name, int indexD, int indexG, int indexS, int indexB )
    : Base(name),
      indexD_(indexD),
      indexG_(indexG),
      indexS_(indexS),
      indexB_(indexB)
  {}

  virtual ~BranchDataMOSFETPowerOp()
  {}

  static complex get(const BranchDataMOSFETPowerOp &op, const Util::Op::OpData &op_data);

  const int           indexD_;
  const int           indexG_;
  const int           indexS_;
  const int           indexB_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataMESFETPowerOp
// Purpose       : Operator for power calculation for MESFETs
// Special Notes :
// Creator       : Pete Sholander, Sandia
// Creation Date : 08/24/2016
//-----------------------------------------------------------------------------
class BranchDataMESFETPowerOp : public Util::Op::Op<BranchDataMESFETPowerOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataMESFETPowerOp(const std::string &name, int indexD, int indexG, int indexS )
    : Base(name),
      indexD_(indexD),
      indexG_(indexG),
      indexS_(indexS)
  {}

  virtual ~BranchDataMESFETPowerOp()
  {}

  static complex get(const BranchDataMESFETPowerOp &op, const Util::Op::OpData &op_data);

  const int           indexD_;
  const int           indexG_;
  const int           indexS_;
};

//-----------------------------------------------------------------------------
// Function      : BranchDataJFETPowerOp
// Purpose       : Operator for power calculation for JFETs
// Special Notes : 
// Creator       : Pete Sholander, Sandia
// Creation Date : 08/27/2016
//-----------------------------------------------------------------------------
class BranchDataJFETPowerOp : public Util::Op::Op<BranchDataJFETPowerOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataJFETPowerOp(const std::string &name, int indexD, int indexG, int indexS )
    : Base(name),
      indexD_(indexD),
      indexG_(indexG),
      indexS_(indexS)
  {}

  virtual ~BranchDataJFETPowerOp()
  {}

  static complex get(const BranchDataJFETPowerOp &op, const Util::Op::OpData &op_data);

  const int           indexD_;
  const int           indexG_;
  const int           indexS_;
};

//-----------------------------------------------------------------------------
// Class         : BranchDataTRAPowerOp
// Purpose       : Operator for power calculation (using I1*V1+I2*V2) for 
//                 the Lossless Transmission Line, which is the T device or TRA
// Special Notes :
// Creator       : Pete Sholander, Sandia
// Creation Date : 9/15/2016
//-----------------------------------------------------------------------------
class BranchDataTRAPowerOp : public Util::Op::Op<BranchDataTRAPowerOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  BranchDataTRAPowerOp(const std::string &name, int index1, int index2 )
    : Base(name),
      index1_(index1),
      index2_(index2)
  {}

  virtual ~BranchDataTRAPowerOp()
  {}

  static complex get(const BranchDataTRAPowerOp &op, const Util::Op::OpData &op_data);
  
  const int           index1_;
  const int           index2_;
};

//-----------------------------------------------------------------------------
// Class         : SensitivityObjFunctionOp
// Purpose       : Operator for getting a value out of the 
//                 SensitivityObjFunction vector.
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class SensitivityObjFunctionOp : public Util::Op::Op<SensitivityObjFunctionOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivityObjFunctionOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivityObjFunctionOp()
  {}

  static complex get(const SensitivityObjFunctionOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpDirectOp : public Util::Op::Op<SensitivitydOdpDirectOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpDirectOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpDirectOp()
  {}

  static complex get(const SensitivitydOdpDirectOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpDirectScaledOp : public Util::Op::Op<SensitivitydOdpDirectScaledOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpDirectScaledOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpDirectScaledOp()
  {}

  static complex get(const SensitivitydOdpDirectScaledOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpAdjointOp : public Util::Op::Op<SensitivitydOdpAdjointOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpAdjointOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpAdjointOp()
  {}

  static complex get(const SensitivitydOdpAdjointOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};


class SensitivitydOdpAdjointScaledOp : public Util::Op::Op<SensitivitydOdpAdjointScaledOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  SensitivitydOdpAdjointScaledOp(const std::string &name, int index)
    : Base(name),
      index_(index)
  {}

  virtual ~SensitivitydOdpAdjointScaledOp()
  {}

  static complex get(const SensitivitydOdpAdjointScaledOp &op, const Util::Op::OpData &op_data);

  const int           index_;
};

//-----------------------------------------------------------------------------
// Class         : MeasureOp
// Purpose       : Operator for getting a .measure result
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class MeasureOp : public Util::Op::Op<MeasureOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  MeasureOp(const std::string &name, const Measure::Base &measure)
    : Base(name),
      measure_(measure)
  {}

  virtual ~MeasureOp()
  {}

  static complex get(const MeasureOp &op, const Util::Op::OpData &op_data);

  const Measure::Base &       measure_;
};

//-----------------------------------------------------------------------------
// Function      : ExpressionOp
// Purpose       : Class for constructing and initalizing an expression
//               : operator, and then getting values from that operator.
// Special Notes : 
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
class ExpressionOp : public Util::Op::Op<ExpressionOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  ExpressionOp(
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, 
    const std::string &name, 
    const std::string &expression, 
    Parallel::Machine comm, 
    const OutputMgr &output_manager);

  virtual ~ExpressionOp()
  {}

  void init(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager, const IO::OutputMgr &output_manager);

  static complex get(const ExpressionOp &op, const Util::Op::OpData &op_data);

  bool getIsComplex () const;

  mutable Util::ExpressionData          expressionData_;
  Parallel::Machine                     comm_;
  const OutputMgr &                     outputMgr_;

  const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp_; 
};

//-----------------------------------------------------------------------------
// Function      : ExpressionRealOp
// Purpose       : Class for constructing and initalizing an expression
//               : operator, and then getting values from that operator.
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 6/25/2020
//-----------------------------------------------------------------------------
class ExpressionRealOp : public Util::Op::Op<ExpressionRealOp, Util::Op::ReduceNone>
{
public:
  ExpressionRealOp(
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, 
    const std::string &name, 
    const std::string &expression, 
    Parallel::Machine comm, 
    const OutputMgr &output_manager);

  ExpressionRealOp( const ExpressionOp & op );

  virtual ~ExpressionRealOp()
  {}

  void init(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager, const IO::OutputMgr &output_manager);

  static complex get(const ExpressionRealOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  mutable Util::ExpressionData          expressionData_;
  Parallel::Machine                     comm_;
  const OutputMgr &                     outputMgr_;
  const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp_; 
};

//-----------------------------------------------------------------------------
// Function      : ExpressionImaginaryOp
// Purpose       : Class for constructing and initalizing an expression
//               : operator, and then getting values from that operator.
// Special Notes : 
// Creator       : Eric Keiter, SNL
// Creation Date : 6/25/2020
//-----------------------------------------------------------------------------
class ExpressionImaginaryOp : public Util::Op::Op<ExpressionImaginaryOp, Util::Op::ReduceNone>
{
public:
  ExpressionImaginaryOp(
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, 
    const std::string &name, 
    const std::string &expression, 
    Parallel::Machine comm, 
    const OutputMgr &output_manager);

  ExpressionImaginaryOp( const ExpressionOp & op);

  virtual ~ExpressionImaginaryOp()
  {}

  void init(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager, const IO::OutputMgr &output_manager);

  static complex get(const ExpressionImaginaryOp &op, const Util::Op::OpData &op_data);
  static complex eval(complex result);

  mutable Util::ExpressionData          expressionData_;
  Parallel::Machine                     comm_;
  const OutputMgr &                     outputMgr_;
  const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp_; 
};



} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Op_h

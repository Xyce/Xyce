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

#include <Xyce_config.h>

#include <N_IO_fwd.h>

#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_DEV_Op.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_MeasureBase.h>
#include <N_IO_Op.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : CurrentIndexOp::get
// Purpose       : get the "index" from the currently active outputter
// Special Notes : the index is basically the line number of output, starting
//                 at zero for the first line
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex CurrentIndexOp::get(const CurrentIndexOp &op, const Util::Op::OpData &op_data)
{
  return op_data.currentIndex_;
}

//-----------------------------------------------------------------------------
// Function      : LeadCurrentIndexOp::get
// Purpose       : get the "index" from the currently active outputter
// Special Notes : the index is basically the line number of output, starting
//                 at zero for the first line
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex LeadCurrentIndexOp::get(const LeadCurrentIndexOp &op, const Util::Op::OpData &op_data)
{
  return op_data.currentIndex_;
}

//-----------------------------------------------------------------------------
// Function      : PowerIndexOp::get
// Purpose       : get the "index" from the currently active outputter
// Special Notes : the index is basically the line number of output, starting
//                 at zero for the first line
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex PowerIndexOp::get(const PowerIndexOp &op, const Util::Op::OpData &op_data)
{
  return op_data.currentIndex_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrTimeOp::get
// Purpose       : get the current simulation time being output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex OutputMgrTimeOp::get(const OutputMgrTimeOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getTime()*op.timeScaleFactor_;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrOutputNoiseOp::get
// Purpose       : get value for total output noise spectral density
//               : at a given frequency
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 1/29/2015
//-----------------------------------------------------------------------------
complex OutputMgrOutputNoiseOp::get(const OutputMgrOutputNoiseOp &op, 
    const Util::Op::OpData &op_data)
{
  return op_data.onoise_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrInputNoiseOp::get
// Purpose       : get value for total input noise spectral density
//               : at a given frequency
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 1/29/2015
//-----------------------------------------------------------------------------
complex OutputMgrInputNoiseOp::get(const OutputMgrInputNoiseOp &op, 
    const Util::Op::OpData &op_data)
{
  return op_data.inoise_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrOutputNoiseContOp::get
// Purpose       : Get either the total noise output contribution for a device,
//               : or the noise output contribution from a specified noise-type 
//               : in that device, at a given frequency
// Special Notes : The ADMS devices may have duplicate entries for a given 
//                 noise type.  So, the typeIndex is a vector-of-ints rather
//                 than a single index.  Only one index should have a non-zero
//                 contribution though.
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 11/20/2017
//-----------------------------------------------------------------------------
complex OutputMgrOutputNoiseContOp::get(const OutputMgrOutputNoiseContOp &op, 
    const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.devIndex_ != -1 && op_data.noiseDataVec_ != 0) 
  {
    if (!op.typeIndex_.empty())
    {
      // noise output contribution from a specified noise-type in the device
      for (std::vector<int>::const_iterator it=op.typeIndex_.begin(); it!=op.typeIndex_.end() ;it++)
      { 
        result += (*op_data.noiseDataVec_)[op.devIndex_]->outputNoiseDens[*it];
      }
    }
    else
    {
      // total noise output contribution for the device
      result = (*op_data.noiseDataVec_)[op.devIndex_]->totalOutputNoise;
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrInputNoiseContOp::get
// Purpose       : Get either the total noise input contribution for a device,
//               : or the noise input contribution from a specified noise-type 
//               : in the device, at a given frequency
// Special Notes : The ADMS devices may have duplicate entries for a given 
//                 noise type.  So, the typeIndex is a vector-of-ints rather
//                 than a single index.  Only one index should have a non-zero
//                 contribution though.
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 11/20/2017
//-----------------------------------------------------------------------------
complex OutputMgrInputNoiseContOp::get(const OutputMgrInputNoiseContOp &op, 
    const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.devIndex_ != -1 && op_data.noiseDataVec_ != 0) 
  {
    if (!op.typeIndex_.empty() )
    {
      // noise input contribution from a specified noise-type in the device
      for (std::vector<int>::const_iterator it=op.typeIndex_.begin(); it!=op.typeIndex_.end() ;it++)
      { 
        result += (*op_data.noiseDataVec_)[op.devIndex_]->inputNoiseDens[*it];
      }
    }
    else
    {
      // total input output contribution for the device
      result = (*op_data.noiseDataVec_)[op.devIndex_]->totalInputNoise;
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrFrequencyOp::get
// Purpose       : get the current frequency being output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex OutputMgrFrequencyOp::get(const OutputMgrFrequencyOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getFrequency();
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrTemperatureOp::get
// Purpose       : get the current simulation temperature
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex OutputMgrTemperatureOp::get(const OutputMgrTemperatureOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getTemperature();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrStepSweepOp::get
// Purpose       : get the current value of the step parameter being swept.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex OutputMgrStepSweepOp::get(const OutputMgrStepSweepOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getStepSweep(op.index_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrDCSweepOp::get
// Purpose       : get the current value of the DC voltage being swept.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex OutputMgrDCSweepOp::get(const OutputMgrDCSweepOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getDCSweep(op.index_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrDCSweepCurrentValueOp::get
// Purpose       : get the current value of the DC voltage being output.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex OutputMgrDCSweepCurrentValueOp::get(const OutputMgrDCSweepCurrentValueOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getPRINTDCvalue();
}

//-----------------------------------------------------------------------------
// Function      : StepNumOp::get
// Purpose       : get the current step number
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/19/2019
//-----------------------------------------------------------------------------
complex
StepNumOp::get(const StepNumOp &op, const Util::Op::OpData &op_data)
{
  return op.outputMgr_.getStepNumber();
}

//-----------------------------------------------------------------------------
// Function      : SolutionOp::get
// Purpose       : get the current value of a solution vector element
// Special Notes : The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionOp::get(const SolutionOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0) 
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionRealOp::get
// Purpose       : get a solution variable in preparation for computing the
//                 real part
// Special Notes : Actually just gets the solution variable.  It does not
//                 take the real part.  The eval function does that.
//                 The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionRealOp::get(const SolutionRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0) 
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionRealOp::eval
// Purpose       : take the real part of a solution vector element
// Special Notes : Actually just returns the real part of a given complex
//                 value.  It does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionRealOp::eval(complex result)
{
  return result.real();
}

//-----------------------------------------------------------------------------
// Function      : SolutionImaginaryOp::get
// Purpose       : get a solution variable in preparation for computing
//                 the imaginary part
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the imaginary part.  The eval function does that.
//                 The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionImaginaryOp::get(const SolutionImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0) 
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionImaginaryOp::eval
// Purpose       : take the imaginary part of a complex number
// Special Notes : Actually just returns the imaginary part of a given complex
//                 value.  It does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : SolutionMagnitudeOp::get
// Purpose       : get the magnitude of a solution vector element
// Special Notes : Actually just gets the solution variable.  It does not
//                 take the magnitude. The eval function does that.
//                 The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionMagnitudeOp::get(const SolutionMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionMagnitudeOp::eval
// Purpose       : take the magnitude of a solution vector element
// Special Notes : Actually just takes the magnitude of a given complex
//                 value.  It does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}

//-----------------------------------------------------------------------------
// Function      : SolutionPhaseDegOp::get
// Purpose       : get the phase of a solution vector element in degrees
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the phase.  The eval function does that.
//                 The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
complex SolutionPhaseDegOp::get(const SolutionPhaseDegOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionPhaseDegOp::eval
// Purpose       : compute the phase of a solution vector element in degrees
// Special Notes : Actually just computes the phase of a given complex
//                 value.  It does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionPhaseDegOp::eval(complex result)
{
  return std::arg(result)*180.0/M_PI;
}


//-----------------------------------------------------------------------------
// Function      : SolutionPhaseRadOp::get
// Purpose       : get the phase of a solution vector element in radians
// Special Notes : Actually just gets the solution variable.  Does not
//                 take the phase.  The eval function does that.
//                 The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionPhaseRadOp::get(const SolutionPhaseRadOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionPhaseRadOp::eval
// Purpose       : compute the phase of a solution vector element in radians
// Special Notes : Actually just computes the phase of a given complex
//                 value.  It does NOT access the solution vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionPhaseRadOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : SolutionDecibelsOp::get
// Purpose       : get the magnitude (in dB) of a solution vector element
// Special Notes : Actually just gets the solution variable.  Does not
//                 find the magnitude.  The eval function does that.
//                 The op.index_ will be:
//                    a) -1 if it is a Ground node.
//                    b) the node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionDecibelsOp::get(const SolutionDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1 && op_data.realSolutionVector_ != 0)
  {
    result = complex((*op_data.realSolutionVector_)[op.index_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SolutionDecibelsOp::eval
// Purpose       : compute the magnitude (in dB) of a solution vector element
// Special Notes : Actually just computes the magnitude (in dB) of a given 
//                 complex value.  It does NOT access the solution vector 
//                 itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SolutionDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(abs(result));
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like V(A,B)
// Special Notes : op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceOp::get(const VoltageDifferenceOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);
  
  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceRealOp::get
// Purpose       : get the difference between two solution variables, as 
//                 needed by constructs like VR(A,B)
// Special Notes : Computes the difference, but does not take the real part.
//                 The real part is taken by the eval function.
//                 op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceRealOp::get(const VoltageDifferenceRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceRealOp::eval
// Purpose       : take the real part of a voltage difference. This is needed
//                 by constructs like VR(A,B)
// Special Notes : Actually just takes the real part of a complex number.
//                 It does NOT access the solution vector itself.  Must 
//                 "get" the complex voltage-difference first, with the get 
//                 function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceImagOp::get
// Purpose       : get the difference between two solution variables, as 
//                 needed by constructs like VI(A,B)
// Special Notes : Computes the difference, but does not take the imaginary
//                 part.  The eval function takes the imaginary part.
//                 op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceImaginaryOp::get(const VoltageDifferenceImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceImagOp::eval
// Purpose       : take the imaginary part of a voltage difference
// Special Notes : Actually just takes the imaginary part of a complex number.
//                 It does NOT access the solution vector itself.  Must 
//                 "get" the complex voltage-difference first, with the get 
//                 function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like VM(A,B)
// Special Notes : Computes the difference, but does not compute the magnitude.
//                 The eval function does that.
//                 op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceMagnitudeOp::get(const VoltageDifferenceMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceMagnitudeOp::eval
// Purpose       : Compute the magnitude of a voltage difference
// Special Notes : Actually just computes the magnitude of a complex number.
//                 It does NOT access the solution vector itself.  Must 
//                 "get" the complex voltage-difference first, with the get 
//                 function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferencePhaseDegOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like VP(A,B)
// Special Notes : Computes the difference, but does not take the phase.
//                 The eval function takes the phase in degrees.
//                 op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
complex VoltageDifferencePhaseDegOp::get(const VoltageDifferencePhaseDegOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferencePhaseDegOp::eval
// Purpose       : Compute the phase of a voltage difference
// Special Notes : Actually just takes the phase of a complex number in degrees.
//                 It does NOT access the solution vector itself.  Must
//                 "get" the complex voltage-difference first, with the get
//                 function.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
complex VoltageDifferencePhaseDegOp::eval(complex result)
{
  return std::arg(result)*180.0/M_PI;
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferencePhaseRadOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like VP(A,B)
// Special Notes : Computes the difference, but does not take the phase.
//                 The eval function takes the phase in radians.
//                 op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferencePhaseRadOp::get(const VoltageDifferencePhaseRadOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferencePhaseRadOp::eval
// Purpose       : Compute the phase of a voltage difference
// Special Notes : Actually just takes the phase of a complex number in radians.
//                 It does NOT access the solution vector itself.  Must
//                 "get" the complex voltage-difference first, with the get
//                 function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferencePhaseRadOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceDecibelsOp::get
// Purpose       : get the difference between two solution variables, as
//                 needed by constructs like VDB(A,B)
// Special Notes : Computes the difference, but does not compute the magnitude
//                 in dB.  The eval function does that.
//                 op.index1_ (or op.index2_) will be:
//                   a) -2 if its node is not found on this processor.
//                   b) -1 if its node is a Ground node.
//                   c) its node index (0 ...N), otherwise.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceDecibelsOp::get(const VoltageDifferenceDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.realSolutionVector_ != 0)
  {
    if (op.index1_ > -1)
    {
      result = complex((*op_data.realSolutionVector_)[op.index1_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index1_]);
    }
    if (op.index2_ > -1)
    {
      result -= complex((*op_data.realSolutionVector_)[op.index2_], op_data.imaginarySolutionVector_ == 0 ? 0.0 : (*op_data.imaginarySolutionVector_)[op.index2_]);
    }
  }
  return result;
}

//-----------------------------------------------------------------------------
// Function      : VoltageDifferenceDecibelsOp::eval
// Purpose       : Compute the magnitude of a voltage difference (in dB)
// Special Notes : Actually just computes the magnitude (in db) of a complex
//                 number. It does NOT access the solution vector itself. Must 
//                 "get" the complex voltage-difference first, with the get 
//                 function.
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex VoltageDifferenceDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(std::abs(result));
}

//-----------------------------------------------------------------------------
// Function      : RFparamsOp::get
// Purpose       : get RF parameter values, such as S(1,2), Y(1,2) or Z(1,2)
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsOp::get(const RFparamsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

   if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

  return result;
}


//-----------------------------------------------------------------------------
// Function      : RFparamsRealOp::get
// Purpose       : get a variable out of the RFparams map, in preparation for
//                 computing its real part with the eval function.  This is
//                 used by operators such as SR(1,2), YR(1,2) or Z(1,2).
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsRealOp::get(const RFparamsRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : RFparamsRealOp::eval
// Purpose       : Take the real part of an RF parameter. This is used
//                 by constructs like SR(1,2), YR(1,2) and ZR(1,2).
// Special Notes : Actually just takes the real part of a complex number.
//                 It does NOT access the RFparams map itself.  The get
//                 function does that.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : RFparamsImaginaryOp::get
// Purpose       : get a variable out of the RFparams map, in preparation for
//                 computing its imaginary part with the eval function.  This
//                 is used by operators such as SI(1,2), YI(1,2) or ZI(1,2).
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsImaginaryOp::get(const RFparamsImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

    return result;
}

//-----------------------------------------------------------------------------
// Function      : RFparamsImaginaryOp::eval
// Purpose       : Take the imaginary part of an RF parameter. This is used
//                 by constructs like SI(1,2), YI(1,2) and ZI(1,2).
// Special Notes : Actually just takes the imaginary part of a complex number.
//                 It does NOT access the RFparams map itself.  The get
//                 function does that.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsImaginaryOp::eval(complex result)
{
  return result.imag();
}


//-----------------------------------------------------------------------------
// Function      : RFparamsMagnitudeOp::get
// Purpose       : get a variable out of the RFparams map, in preparation for
//                 computing its magnitude with the eval function.  This
//                 is used by operators such as SI(1,2), YI(1,2) or ZI(1,2).
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsMagnitudeOp::get(const RFparamsMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

   if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : RFparamsMagnitudeOp::eval
// Purpose       : Take the magnitude of an RF parameter. This is used
//                 by constructs like SM(1,2), YM(1,2) and ZM(1,2).
// Special Notes : Actually just takes the magnitude of a complex number.
//                 It does NOT access the RFparams map itself.  The get
//                 function does that.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : RFparamsPhaseDegOp::get
// Purpose       : get a variable out of the RFparams map, in preparation for
//                 computing its phase (in degrees) with the eval function.  This
//                 is used by operators such as SP(1,2), YP(1,2) or ZP(1,2).
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsPhaseDegOp::get(const RFparamsPhaseDegOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : RFparamsPhaseDegOp::eval
// Purpose       : Take the phase (in degrees) of an RF parameter. This is used
//                 by constructs like SP(1,2), YP(1,2) and ZP(1,2).
// Special Notes : Actually just takes the phase of a complex number.
//                 It does NOT access the RFparams map itself.  The get
//                 function does that.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsPhaseDegOp::eval(complex result)
{
  return std::arg(result)*180.0/M_PI;
}


//-----------------------------------------------------------------------------
// Function      : RFparamsPhaseRadOp::get
// Purpose       : get a variable out of the RFparams map, in preparation for
//                 computing its phase (in radians) with the eval function.  This
//                 is used by operators such as SP(1,2), YP(1,2) or ZP(1,2).
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsPhaseRadOp::get(const RFparamsPhaseRadOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : RFparamsPhaseRadOp::eval
// Purpose       : Take the phase (in radians) of an RF parameter. This is used
//                 by constructs like SP(1,2), YP(1,2) and ZP(1,2).
// Special Notes : Actually just takes the phase of a complex number.
//                 It does NOT access the RFparams map itself.  The get
//                 function does that.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsPhaseRadOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : RFparamsDecibelsOp::get
// Purpose       : get a variable out of the RFparams map, in preparation for
//                 computing its magnitude (in dB) with the eval function.  This
//                 is used by operators such as SDB(1,2), YDB(1,2) or ZDB(1,2).
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsDecibelsOp::get(const RFparamsDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op_data.RFparams_ != 0)
  {
    Util::Op::RFparamsData::const_iterator it;
    it = (*op_data.RFparams_).find(op.type_);
    const Teuchos::SerialDenseMatrix<int, std::complex<double> > & param = *it->second;
    if ( op.index1_ > 0  && op.index2_ > 0 && op.index1_ <= param.numRows() && op.index2_ <= param.numRows() )
    {
      // Teuchos matrices start at (0,0), so subtract 1 from index1_ and index2_
      result = param(op.index1_-1,op.index2_-1);
    }
    else
    {
      Report::UserError0() << "Indices for " << op.name_ << " operator must be <= number of ports";
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : RFparamsDecibelsOp::eval
// Purpose       : Take the magnitude (in dB) of an RF parameter. This is used
//                 by constructs like SDB(1,2), YDB(1,2) and ZDB(1,2).
// Special Notes : Actually just takes the magnitude (in dB) of a complex number.
//                 It does NOT access the RFparams map itself.  The get
//                 function does that.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/01/2019
//-----------------------------------------------------------------------------
complex RFparamsDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(std::abs(result));
}


//-----------------------------------------------------------------------------
// Function      : StateOp::get
// Purpose       : Get a value out of the state vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StateOp::get(const StateOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ > -1)
  {
    result = complex(op_data.stateVector_ == 0 ? 0.0 : (*op_data.stateVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreOp::get
// Purpose       : Get a value out of the Store vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreOp::get(const StoreOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realStoreVector_ == 0 ? 0.0 : (*op_data.realStoreVector_)[op.index_], 0.0);
  }

  return result;
}


//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentOp::get
// Purpose       : Get a value out of the lead current vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex BranchDataCurrentOp::get(const BranchDataCurrentOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataVoltageOp::get
// Purpose       : Get a value out of the branch data junction voltage vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex BranchDataVoltageOp::get(const BranchDataVoltageOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realLeadCurrentDeltaVVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentDeltaVVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataPosNegPowerOp::get
// Purpose       : Compute the power using I*V
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Sandia
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex BranchDataPosNegPowerOp::get(const BranchDataPosNegPowerOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    if( (op_data.realLeadCurrentVector_ == 0) || (op_data.realLeadCurrentDeltaVVector_ == 0) )
    {
      result = 0.0;
    }
    else
    {
      result = (*op_data.realLeadCurrentVector_)[op.index_]*(*op_data.realLeadCurrentDeltaVVector_)[op.index_];
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataBJTPowerOp::get
// Purpose       : Power computation for BJTs
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Sandia
// Creation Date : 05/05/2016
//-----------------------------------------------------------------------------
complex BranchDataBJTPowerOp::get(const BranchDataBJTPowerOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.indexB_ != -1 && op.indexC_ != -1 )
  {
    if( (op_data.realLeadCurrentVector_ == 0) || (op_data.realLeadCurrentDeltaVVector_ == 0) )
    {
      result = 0.0;
    }
    else
    {
      result = std::abs((*op_data.realLeadCurrentVector_)[op.indexB_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexB_]) +
               std::abs((*op_data.realLeadCurrentVector_)[op.indexC_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexC_]);
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataJFETPowerOp::get
// Purpose       : Power computation for JFETs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Sandia
// Creation Date : 8/27/2016
//-----------------------------------------------------------------------------
complex BranchDataJFETPowerOp::get(const BranchDataJFETPowerOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.indexD_ != -1 && op.indexG_ != -1)
  {
    if( (op_data.realLeadCurrentVector_ == 0) || (op_data.realLeadCurrentDeltaVVector_ == 0) )
    {
      result = 0.0;
    }
    else
    {
      result = (*op_data.realLeadCurrentVector_)[op.indexD_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexD_] +
               (*op_data.realLeadCurrentVector_)[op.indexG_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexG_];     
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataMESFETPowerOp::get
// Purpose       : Power computation for MESFETs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Sandia
// Creation Date : 8/24/2016
//-----------------------------------------------------------------------------
complex BranchDataMESFETPowerOp::get(const BranchDataMESFETPowerOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.indexD_ != -1 && op.indexG_ != -1)
  {
    if( (op_data.realLeadCurrentVector_ == 0) || (op_data.realLeadCurrentDeltaVVector_ == 0) )
    {
      result = 0.0;
    }
    else
    {
      result = (*op_data.realLeadCurrentVector_)[op.indexD_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexD_] +
               (*op_data.realLeadCurrentVector_)[op.indexG_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexG_];     
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataMOSFETPowerOp::get
// Purpose       : Power computation for MOSFETs
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Sandia
// Creation Date : 8/18/2016
//-----------------------------------------------------------------------------
complex BranchDataMOSFETPowerOp::get(const BranchDataMOSFETPowerOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.indexD_ != -1 && op.indexG_ != -1)
  {
    if( (op_data.realLeadCurrentVector_ == 0) || (op_data.realLeadCurrentDeltaVVector_ == 0) )
    {
      result = 0.0;
    }
    else
    {
      result = (*op_data.realLeadCurrentVector_)[op.indexD_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexD_] +
               (*op_data.realLeadCurrentVector_)[op.indexG_]*(*op_data.realLeadCurrentDeltaVVector_)[op.indexG_];     
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataTRAPowerOp::get
// Purpose       : Compute the power using I1*V1+I2*V2 for Lossless Transmission 
//                 Line, which is the T device or TRA
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Sandia
// Creation Date : 9/15/2016
//-----------------------------------------------------------------------------
complex BranchDataTRAPowerOp::get(const BranchDataTRAPowerOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index1_ != -1 && op.index1_ != -1)
  {
    if( (op_data.realLeadCurrentVector_ == 0) || (op_data.realLeadCurrentDeltaVVector_ == 0) )
    {
      result = 0.0;
    }
    else
    {
      result = (*op_data.realLeadCurrentVector_)[op.index1_]*(*op_data.realLeadCurrentDeltaVVector_)[op.index1_] +
               (*op_data.realLeadCurrentVector_)[op.index2_]*(*op_data.realLeadCurrentDeltaVVector_)[op.index2_];     
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreRealOp::get
// Purpose       : get a store variable in preparation for computing the real
//                 part with the eval function.
// Special Notes : This method only works if a valid real store vector pointer
//                 is always passed.  If a null pointer is passed,
//                 a nonsense value is stored instead.  
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreRealOp::get(const StoreRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1 && op_data.realStoreVector_ != 0) {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreRealOp::eval
// Purpose       : take the real part of a complex number.  
// Special Notes : It does NOT access access the store vector itself.  
//                 Must "get" the value first, with the get function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreRealOp::eval(complex result)
{
  return result.real();
}


//-----------------------------------------------------------------------------
// Function      : StoreImaginaryOp::get
// Purpose       : get a store variable in preparation for computing
//                 imaginary part with the eval function.
// Special Notes : This method only works if a valid real store vector pointer
//                 is always passed.  If a null pointer is passed,
//                 a nonsense value is stored instead.  
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreImaginaryOp::get(const StoreImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1 && op_data.realStoreVector_ != 0) 
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreImaginaryOp::eval
// Purpose       : take the imaginary part of a complex number
// Special Notes : It does NOT access access the store vector itself.  
//                 Must "get" the value first, with the get function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreImaginaryOp::eval(complex result)
{
  return result.imag();
}

//-----------------------------------------------------------------------------
// Function      : StoreMagnitudeOp::get
// Purpose       : get the magnitude of a store vector element
// Special Notes : Actually just gets the store variable.  Does not
//                 take the magnitude.  The eval function does that.
//                 This method only works if a valid real store vector pointer
//                 is always passed.  If a null pointer is passed,
//                 a nonsense value is stored instead.  
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreMagnitudeOp::get(const StoreMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1 && op_data.realStoreVector_ != 0)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreMagnitudeOp::eval
// Purpose       : take the magnitude of a store vector element
// Special Notes : Actually just takes the magnitude of a given complex
//                 value, does NOT access the store vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : StorePhaseDegOp::get
// Purpose       : get the phase of a store vector element in degrees
// Special Notes : Actually just gets the store variable.  Does not
//                 take the phase.  The eval function does that.
//                 This method only works if a valid real store vector pointer
//                 is always passed.  If a null pointer is passed,
//                 a nonsense value is stored instead.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
complex StorePhaseDegOp::get(const StorePhaseDegOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1 && op_data.realStoreVector_ != 0)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StorePhaseDegOp::eval
// Purpose       : compute the phase of a store vector element in degrees
// Special Notes : Actually just computes the phase of a given complex
//                 value, does NOT access the store vector itself.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 7/12/2019
//-----------------------------------------------------------------------------
complex StorePhaseDegOp::eval(complex result)
{
  return std::arg(result)*180.0/M_PI;
}


//-----------------------------------------------------------------------------
// Function      : StorePhaseRadOp::get
// Purpose       : get the phase of a store vector element in radians
// Special Notes : Actually just gets the store variable.  Does not
//                 take the phase.  The eval function does that.
//                 This method only works if a valid real store vector pointer
//                 is always passed.  If a null pointer is passed,
//                 a nonsense value is stored instead.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StorePhaseRadOp::get(const StorePhaseRadOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1 && op_data.realStoreVector_ != 0)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StorePhaseRadOp::eval
// Purpose       : compute the phase of a store vector element in radians
// Special Notes : Actually just computes the phase of a given complex
//                 value, does NOT access the store vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StorePhaseRadOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : StoreDecibelsOp::get
// Purpose       : get the magnitude (in dB) of a store vector element
// Special Notes : Actually just gets the store variable.  Does not
//                 find the magnitude.  The eval function does that.
//                 This method only works if a valid real store vector pointer
//                 is always passed.  If a null pointer is passed,
//                 a nonsense value is stored instead.  
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreDecibelsOp::get(const StoreDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1 && op_data.realStoreVector_ != 0)
  {
    result = complex((*op_data.realStoreVector_)[op.index_], op_data.imaginaryStoreVector_ == 0 ? 0.0 : (*op_data.imaginaryStoreVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : StoreDecibelsOp::eval
// Purpose       : compute the magnitude (in dB) of a store vector element
// Special Notes : Actually just computes the magnitude (in dB) of a given
//                 complex value, does NOT access the store vector itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex StoreDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(abs(result));
}

//-----------------------------------------------------------------------------
// Function      : SensitivityObjFunctionOp::get
// Purpose       : Get a value out of the SensitivityObjFunction vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SensitivityObjFunctionOp::get(const SensitivityObjFunctionOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.objectiveVector_ == 0 || op.index_ >= op_data.objectiveVector_->size() ? 0.0 : (*op_data.objectiveVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpDirectOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SensitivitydOdpDirectOp::get(const SensitivitydOdpDirectOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpDirectVector_ == 0 || op.index_ >= op_data.dOdpDirectVector_->size() ? 0.0 : (*op_data.dOdpDirectVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpDirectScaledOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SensitivitydOdpDirectScaledOp::get(const SensitivitydOdpDirectScaledOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpDirectScaledVector_ == 0 || op.index_ >= op_data.dOdpDirectScaledVector_->size() ? 0.0 : (*op_data.dOdpDirectScaledVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpAdjointOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SensitivitydOdpAdjointOp::get(const SensitivitydOdpAdjointOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpAdjointVector_ == 0 || op.index_ >= op_data.dOdpAdjointVector_->size() ? 0.0 : (*op_data.dOdpAdjointVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : SensitivitydOdpAdjointScaledOp::get
// Purpose       : Get a value out of the SensitivityParameter vector.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex SensitivitydOdpAdjointScaledOp::get(const SensitivitydOdpAdjointScaledOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.dOdpAdjointScaledVector_ == 0 || op.index_ >= op_data.dOdpAdjointScaledVector_->size() ? 0.0 : (*op_data.dOdpAdjointScaledVector_)[op.index_], 0.0);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : MeasureOp::get
// Purpose       : Get a .measure result
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex MeasureOp::get(const MeasureOp &op, const Util::Op::OpData &op_data)
{
  complex result(const_cast<Measure::Base &>(op.measure_).getMeasureResult(), 0.0);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::ExpressionOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes string as second argument.  expressionData_
//                 constructor will process the string into an expression
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
ExpressionOp::ExpressionOp(
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, 
    const std::string &name, 
    const std::string &expression, 
    Parallel::Machine comm, 
    const OutputMgr &output_manager)
  : Base(name),
    expressionData_(grp, expression),
    comm_(comm),
    outputMgr_(output_manager),
    grp_(grp)
{
  init(comm, output_manager.getOpBuilderManager(), output_manager);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::init
// Purpose       : initialize an expression
// Special Notes : runs the ExpressionData::setup method to resolve
//                 symbols
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void ExpressionOp::init(
  Parallel::Machine                     comm,
  const Util::Op::BuilderManager &      op_builder_manager,
  const IO::OutputMgr &                 output_manager)
{
  expressionData_.setup(comm,
                        op_builder_manager,
                        output_manager.getMainContextFunctionMap(),
                        output_manager.getMainContextParamMap(),
                        output_manager.getMainContextGlobalParamMap());
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::get
// Purpose       : evaluate an expression
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex ExpressionOp::get(const ExpressionOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0,0.0);
  op.expressionData_.evaluate(op.comm_, 
        op.outputMgr_.getCircuitTime(), 
        op.outputMgr_.getCircuitTimeStep(), 
        op_data,result);
  return result;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionOp::
// Purpose       : evaluate an expression
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/14/2023
//-----------------------------------------------------------------------------
bool ExpressionOp::getIsComplex() const
{
  return expressionData_.getIsComplex();
}

//-----------------------------------------------------------------------------
// Function      : ExpressionRealOp::ExpressionRealOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes string as second argument.  expressionData_
//                 constructor will process the string into an expression
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
ExpressionRealOp::ExpressionRealOp(
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, 
    const std::string &name, 
    const std::string &expression, 
    Parallel::Machine comm, 
    const OutputMgr &output_manager)
  : Base(name),
    expressionData_(grp, expression),
    comm_(comm),
    outputMgr_(output_manager),
    grp_(grp)
{
  init(comm, output_manager.getOpBuilderManager(), output_manager);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionRealOp::ExpressionRealOp
// Purpose       : Constructor for expression Op
// Special Notes : 
//                 
// Scope         : 
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
ExpressionRealOp::ExpressionRealOp( const ExpressionOp & op )
  : Base(std::string("Re(" + op.getName() + ")")),
    expressionData_(op.grp_, op.expressionData_.getExpression()),
    comm_(op.comm_),
    outputMgr_(op.outputMgr_),
    grp_(op.grp_)
{
  init(comm_, outputMgr_.getOpBuilderManager(), outputMgr_);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionRealOp::init
// Purpose       : initialize an expression
// Special Notes : runs the ExpressionData::setup method to resolve
//                 symbols
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void ExpressionRealOp::init(
  Parallel::Machine                     comm,
  const Util::Op::BuilderManager &      op_builder_manager,
  const IO::OutputMgr &                 output_manager)
{
  expressionData_.setup(comm,
                        op_builder_manager,
                        output_manager.getMainContextFunctionMap(),
                        output_manager.getMainContextParamMap(),
                        output_manager.getMainContextGlobalParamMap());
}

//-----------------------------------------------------------------------------
// Function      : ExpressionRealOp::get
// Purpose       : evaluate an expression
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex ExpressionRealOp::get(const ExpressionRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0,0.0);
  op.expressionData_.evaluate(op.comm_, 
        op.outputMgr_.getCircuitTime(), 
        op.outputMgr_.getCircuitTimeStep(), 
        op_data,result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionRealOp::eval
// Purpose       : take the real part of a solution vector element
// Special Notes : Actually just returns the real part of a given complex
//                 value.  It does NOT access the solution vector itself.
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
complex ExpressionRealOp::eval(complex result)
{
  return result.real();
}

//-----------------------------------------------------------------------------
// Function      : ExpressionImaginaryOp::ExpressionImaginaryOp
// Purpose       : Constructor for expression Op
// Special Notes : Takes string as second argument.  expressionData_
//                 constructor will process the string into an expression
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
ExpressionImaginaryOp::ExpressionImaginaryOp(
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & grp, 
    const std::string &name, 
    const std::string &expression, 
    Parallel::Machine comm, 
    const OutputMgr &output_manager)
  : Base(name),
    expressionData_(grp, expression),
    comm_(comm),
    outputMgr_(output_manager),
    grp_(grp)
{
  init(comm, output_manager.getOpBuilderManager(), output_manager);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionImaginaryOp::ExpressionImaginaryOp
// Purpose       : Constructor for expression Op
// Special Notes : 
//                 
// Scope         : 
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
ExpressionImaginaryOp::ExpressionImaginaryOp( const ExpressionOp & op )
  : Base(std::string("Im(" + op.getName() + ")")),
    expressionData_(op.grp_, op.expressionData_.getExpression()),
    comm_(op.comm_),
    outputMgr_(op.outputMgr_),
    grp_(op.grp_)
{
  init(comm_, outputMgr_.getOpBuilderManager(), outputMgr_);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionImaginaryOp::init
// Purpose       : initialize an expression
// Special Notes : runs the ExpressionData::setup method to resolve
//                 symbols
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void ExpressionImaginaryOp::init(
  Parallel::Machine                     comm,
  const Util::Op::BuilderManager &      op_builder_manager,
  const IO::OutputMgr &                 output_manager)
{
  expressionData_.setup(comm,
                        op_builder_manager,
                        output_manager.getMainContextFunctionMap(),
                        output_manager.getMainContextParamMap(),
                        output_manager.getMainContextGlobalParamMap());
}

//-----------------------------------------------------------------------------
// Function      : ExpressionImaginaryOp::get
// Purpose       : evaluate an expression
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
complex ExpressionImaginaryOp::get(const ExpressionImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0,0.0);
  op.expressionData_.evaluate(op.comm_, 
        op.outputMgr_.getCircuitTime(), 
        op.outputMgr_.getCircuitTimeStep(), 
        op_data,result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionImaginaryOp::eval
// Purpose       : take the imaginary part of a complex number
// Special Notes : It does NOT access access the store vector itself.  
//                 Must "get" the value first, with the get function.
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
complex ExpressionImaginaryOp::eval(complex result)
{
  return result.imag();
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentRealOp::get
// Purpose       : get a variable out of the lead current vector, in 
//                 preparation for computing its real part with the eval
//                 function.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentRealOp::get(const BranchDataCurrentRealOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) 
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], op_data.imaginaryLeadCurrentVector_ == 0 ? 0.0 : (*op_data.imaginaryLeadCurrentVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentRealOp::eval
// Purpose       : take the real part of a lead current vector
// Special Notes : Actually just returns the real part of a given complex
//                 value.  It does NOT access the lead current vector itself.
//                 The get function does that access.
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentRealOp::eval(complex result)
{
  return result.real();
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentImaginaryOp::get
// Purpose       : get a variable out of the lead current vector, in 
//                 preparation for computing its imaginary part with the
//                 eval function. 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentImaginaryOp::get(const BranchDataCurrentImaginaryOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) 
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], op_data.imaginaryLeadCurrentVector_ == 0 ? 0.0 : (*op_data.imaginaryLeadCurrentVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentImaginaryOp::eval
// Purpose       : take the imaginary part of a lead current vector
// Special Notes : Actually just returns the imaginary part of a given complex
//                 value.  It does NOT access the lead current vector itself.
//                 The get function does that access.
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentImaginaryOp::eval(complex result)
{
  return result.imag();
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentMagnitudeOp::get
// Purpose       : get a variable out of the lead current vector, in 
//                 preparation for computing its magnitude with the eval
//                 function. 
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentMagnitudeOp::get(const BranchDataCurrentMagnitudeOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1) 
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], op_data.imaginaryLeadCurrentVector_ == 0 ? 0.0 : (*op_data.imaginaryLeadCurrentVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentMagnitudeOp::eval
// Purpose       : take the magnitude of a lead current vector
// Special Notes : Actually just returns the magnitude of a given complex
//                 value.  It does NOT access the lead current vector itself.
//                 The get function does that access.
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentMagnitudeOp::eval(complex result)
{
  return std::abs(result);
}


//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentPhaseDegOp::get
// Purpose       : get a variable out of the lead current vector, in
//                 preparation for computing its phase (in degrees) with the eval
//                 function.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015 
//-----------------------------------------------------------------------------
complex BranchDataCurrentPhaseDegOp::get(const BranchDataCurrentPhaseDegOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], op_data.imaginaryLeadCurrentVector_ == 0 ? 0.0 : (*op_data.imaginaryLeadCurrentVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentPhaseDegOp::eval
// Purpose       : compute the phase (in degrees) of a lead current vector element
// Special Notes : Actually just returns the magnitude of a given complex
//                 value.  It does NOT access the lead current vector itself.
//                 The get function does that access.
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015 
//-----------------------------------------------------------------------------
complex BranchDataCurrentPhaseDegOp::eval(complex result)
{
  return std::arg(result)*180.0/M_PI;
}


//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentPhaseRadOp::get
// Purpose       : get a variable out of the lead current vector, in
//                 preparation for computing its phase (in radians) with the eval
//                 function.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015 
//-----------------------------------------------------------------------------
complex BranchDataCurrentPhaseRadOp::get(const BranchDataCurrentPhaseRadOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], op_data.imaginaryLeadCurrentVector_ == 0 ? 0.0 : (*op_data.imaginaryLeadCurrentVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentPhaseRadOp::eval
// Purpose       : compute the phase (in radians) of a lead current vector element
// Special Notes : Actually just returns the magnitude of a given complex
//                 value.  It does NOT access the lead current vector itself.
//                 The get function does that access.
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 02/24/2015
//-----------------------------------------------------------------------------
complex BranchDataCurrentPhaseRadOp::eval(complex result)
{
  return std::arg(result);
}


//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentDecibelsOp::get
// Purpose       : get a variable out of the lead current vector, in 
//                 preparation for computing its magnitude (in dB) 
//                 with the eval function. 
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015 
//-----------------------------------------------------------------------------
complex BranchDataCurrentDecibelsOp::get(const BranchDataCurrentDecibelsOp &op, const Util::Op::OpData &op_data)
{
  complex result(0.0, 0.0);

  if (op.index_ != -1)
  {
    result = complex(op_data.realLeadCurrentVector_ == 0 ? 0.0 : (*op_data.realLeadCurrentVector_)[op.index_], op_data.imaginaryLeadCurrentVector_ == 0 ? 0.0 : (*op_data.imaginaryLeadCurrentVector_)[op.index_]);
  }

  return result;
}

//-----------------------------------------------------------------------------
// Function      : BranchDataCurrentDecibelsOp::eval
// Purpose       : compute the magnitude (in dB) of a lead current vector element
// Special Notes : Actually just returns the magnitude (in dB) of a given complex
//                 value.  It does NOT access the lead current vector itself.
//                 The get function does that access.  
// Scope         : public
// Creator       : Richard Schiek, SNL 
// Creation Date : 02/24/2015 
//-----------------------------------------------------------------------------
complex BranchDataCurrentDecibelsOp::eval(complex result)
{
  return 20.0*std::log10(abs(result));
}

} // namespace IO
} // namespace Xyce

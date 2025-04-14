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

//-------------------------------------------------------------------------
//
// Purpose        : This file contains specializations for the AppleFFT interface
//                  for various vector types.
//
// Special Notes  : 
//
// Creator        : Richard Schiek
//
// Creation Date  : 3/16/23
//
//
//
//
//-------------------------------------------------------------------------
// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <N_UTL_AppleFFT_Interface.hpp>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>
#include <cmath>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Apple FFT Interface specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_AppleFFT_Interface<std::vector<double> >::calculateDFT()
  {
    if(nextLargestPowerOf2_ == signalLength_)
    {
      // copy dftInData_ to yTdReal
      for( int i=0;i<nextLargestPowerOf2_;i++)
      {
        if( i<dftInData_->size())
        {
          yTdReal[i] = (*dftInData_)[i];
        }
        else
        {
          yTdReal[i] = 0.0;
        }
        yTdImag[i] = 0.0;
      }
      
      vDSP_DFT_ExecuteD(forwardSetup_, yTdReal.data(), yTdImag.data(), yFdReal.data(), yFdImag.data());
    
      // now interleave the output vectors for return.
      for( int i=0;i<nextLargestPowerOf2_;i++)
      {
        if( (2*i+1)<dftOutData_->size())
        {
          (*dftOutData_)[2*i] = yFdReal[i];
          (*dftOutData_)[2*i+1] = yFdImag[i];
        }
      }
    }
    else
    {
      for( auto k=0; k<(1+signalLength_/2); k++)
      {
        (*dftOutData_)[2*k] = 0.0;
        (*dftOutData_)[2*k+1] = 0.0;
        for( auto j=0; j<signalLength_; j++)
        {
          (*dftOutData_)[2*k] += std::cos(-2*M_PI*j*k/signalLength_)*(*dftInData_)[j];
          (*dftOutData_)[2*k+1] += std::sin(-2*M_PI*j*k/signalLength_)*(*dftInData_)[j];
        }
      }
    }
  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_AppleFFT_Interface<std::vector<double> >::calculateIFT()
  {
    if(nextLargestPowerOf2_ == signalLength_)
    {
      // copy iftInData_ to yFdReal and yFdImag
      for( int i=0;i<(nextLargestPowerOf2_/2);i++)
      {
        if( (2*i)<iftInData_->size())
        {
          yFdReal[i] = (*iftInData_)[2*i];
          yFdImag[i] = (*iftInData_)[2*i+1];
        }
        else
        {
          yFdReal[i] = 0.0;
          yFdImag[i] = 0.0;
        }
      }
      vDSP_DFT_ExecuteD(inverseSetup_, yFdReal.data(), yFdImag.data(), yTdReal.data(), yTdImag.data() );
      // now combine the output arrays into an interleaved result.
      for( int i=0;i<nextLargestPowerOf2_;i++)
      {
        if( i<iftOutData_->size())
        {
          (*iftOutData_)[i] = scaleFactor_ * yTdReal[i];
        }
      }
    }
    else
    {
      for( auto k=0; k<signalLength_; k++)
      {
        (*iftOutData_)[k] = 0.0;
        for( auto j=0; j<((signalLength_/2)+1); j++)
        {
          (*iftOutData_)[k] += (std::cos(2*M_PI*j*k/signalLength_))*(*iftInData_)[2*j] - (std::sin(2*M_PI*j*k/signalLength_))*(*iftInData_)[2*j+1];
        }
        for( auto j=(signalLength_/2), jc=(signalLength_/2)+1; j>0; j--, jc++)
        {
          (*iftOutData_)[k] += (std::cos(2*M_PI*jc*k/signalLength_))*(*iftInData_)[2*j] + (std::sin(2*M_PI*jc*k/signalLength_))*(*iftInData_)[2*j+1];
        }
        (*iftOutData_)[k] = (*iftOutData_)[k] / signalLength_;
      }
    }
  }

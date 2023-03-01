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

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Apple FFT Interface specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_AppleFFT_Interface<std::vector<double> >::calculateDFT()
  {
    // Although we used a const on input to show that we aren't changing dftInData_
    // we need to cast that away as the FFT library takes non-const pointers.
    //std::vector<double>::const_iterator inDataItr = (this->dftInData_)->begin();
    //double * inDataPtr = const_cast< double * >( &(*inDataItr) );
    //std::vector<double>::iterator outResultItr = (this->dftOutData_)->begin();
    //double * outResultPtr = &(*outResultItr);

    //long status = DftiComputeForward( fftDescriptor, inDataPtr, outResultPtr);
    //checkAndTrapErrors( status );
    
    // need to split interleaved input vector into two input vectors
    // or if the input vector is 1/2 the length of the output vector 
    // then leave the input as it is and make a zero vector of the same length 
    //  Also need to two output vectors for the FFT results.
    
    
    if( interleavedRoutines_ )
    {
      // copy dftInData_ to yTdInterleaved
      for( int i=0;i<nextLargestPowerOf2_;i++)
      {
        if( i<dftInData_->size())
        {
          yTdInterleaved[2*i] = (*dftInData_)[i];
          yTdInterleaved[2*i+1] = 0.0;
        }
        else
        {
          yTdInterleaved[2*i] = 0.0;
          yTdInterleaved[2*i+1] = 0.0;
        }
      }
      vDSP_DFT_Interleaved_ExecuteD(forwardInterleavedSetup_, (DSPDoubleComplex*)(yTdInterleaved.data()), (DSPDoubleComplex*)(yFdInterleaved.data()));
      
      for( int i=0;i<nextLargestPowerOf2_;i++)
      {
        if( i<dftOutData_->size())
        {
          // magnitude of frequency domain values here are a factor of 2 larger
          // than the results from FFTW.  Factor out that 2 here for consistency.
          (*dftOutData_)[i] = yFdInterleaved[i]/2.0;
        }
      }
      // yFdInterleaved[0] is the DC component and yFdInterleaved[1] 
      // is the Nyquist component.  Move the Nyquist component to the 
      // last real position.  Zero the original spot where the nyquist value was stored.
      //(*dftOutData_)[nextLargestPowerOf2_] = (*dftOutData_)[1];
      (*dftOutData_)[1] = 0.0;
    }
    else
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
  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_AppleFFT_Interface<std::vector<double> >::calculateIFT()
  {
    // Although we used a const on input to show that we aren't changing iftInData_
    // we need to cast that away as the FFT library takes non-const pointers.
    //std::vector<double>::const_iterator inDataItr = (this->iftInData_)->begin();
    //double * inDataPtr = const_cast< double * >( &(*inDataItr) );
    //std::vector<double>::iterator outResultItr = (this->iftOutData_)->begin();
    //double * outResultPtr = &(*outResultItr);

    //long status = DftiComputeBackward( fftDescriptor, inDataPtr, outResultPtr);
    //checkAndTrapErrors( status );
    
    // need to break the interleaved input vector into two vectors 
    // make destination arrays for the real and imaginary results 
    
    if( interleavedRoutines_ )
    {
      for( int i=0;i<(2*nextLargestPowerOf2_);i++)
      {
        if( i<iftInData_->size())
        {
          yFdInterleaved[i] = (*iftInData_)[i];
        }
        else
        {
          yFdInterleaved[i] = 0.0;
        }
      }
      
      vDSP_DFT_Interleaved_ExecuteD(inverseInterleavedSetup_, (DSPDoubleComplex*)(yFdInterleaved.data()), (DSPDoubleComplex*)(yTdInterleaved.data()));
      
      for( int i=0;i<nextLargestPowerOf2_;i++)
      {
        if( i<iftOutData_->size())
        {
          (*iftOutData_)[i] = scaleFactor_ * yTdInterleaved[2*i];
        }
      } 
      
    }
    else
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
        // this allows for a complex result in the time domain.  
        // not sure if it is needed.
        /*
        if( (2*i+1)<iftOutData_->size())
        {
          (*iftOutData_)[2*i] = scaleFactor_ * yTdReal[i];
          (*iftOutData_)[2*i+1] = scaleFactor_ * yTdImag[i];
        }
        */
      }
    }
  }
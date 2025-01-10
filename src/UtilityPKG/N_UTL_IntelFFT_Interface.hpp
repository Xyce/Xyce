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

//-------------------------------------------------------------------------
//
// Purpose        : This class is an implementation of the abstract FFT interface
//                  for the Intel FFT library.
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 5/27/08
//
//-------------------------------------------------------------------------
#ifndef N_UTL_INTEL_FFT_INTERFACE_HPP
#define N_UTL_INTEL_FFT_INTERFACE_HPP


// ---------- Standard Includes ----------

#include <N_UTL_FFTInterfaceDecl.hpp>
#include <mkl_dfti.h>
#include <iostream>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface for Intel FFT 
// Purpose       : This class acts as an abstract interface to an FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  
// Special Notes :
// Creator       : Heidi Thornquist
// Creation Date : 11/11/10
// Last Modified : 11/11/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_IntelFFT_Interface: public N_UTL_FFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.
    N_UTL_IntelFFT_Interface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      : N_UTL_FFTInterfaceDecl<VectorType>(length, numSignals, reqStride, overwrite)
    {
      // create the fft descriptor structor
      int fftDimension = 1;  // 1D fft's
      long status = DftiCreateDescriptor( &fftDescriptor, DFTI_DOUBLE, DFTI_REAL, fftDimension, this->signalLength_ );
      checkAndTrapErrors( status );

      // configure the fft library to do numberSignals of 1D FFT's at the same time  
      status = DftiSetValue( fftDescriptor, DFTI_NUMBER_OF_TRANSFORMS, this->numberSignals_);
      checkAndTrapErrors( status );

      // The real input is signalLength long and the output is twice as long, complex
      status = DftiSetValue( fftDescriptor, DFTI_INPUT_DISTANCE, this->signalLength_);
      checkAndTrapErrors( status );
      status = DftiSetValue( fftDescriptor, DFTI_OUTPUT_DISTANCE, 2*this->signalLength_);
      checkAndTrapErrors( status );

      if( !overwrite  )
      {
        // Don't overwrite the input with the results
        status = DftiSetValue( fftDescriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        checkAndTrapErrors( status );
      }

      // I had set this up to use a slightly tighter packing scheme, but FFTW
      // does not support such schemes.  So, we'll generalize and stick to 
      // one that is supported by both to make it easier for Xyce to use 
      // this inteface regardless of the underlying FFT library.
      //
      // Packing scheme is CSS which for real input of:
      // in0, in1, ... , inN-1
      //
      // produces
      // r0, 0.0, r1, c1, ... , rN, cN  (two extra element if N is even)
      // r0, 0.0, r1, c1, ... , rN, cN  (one extra element if N is odd)
      //
      // see
      //
      // http://www.intel.com/software/products/mkl/docs/WebHelp/mklrefman.htm
      //
      // for more details
    
      // The forward and backward transform must have a consistent scale factor
      // so that the inverse of a forward transform is the same signal.  By default
      // we'll use 1.0 for the forward scale factor and then 1/n for the inverse transform.
      double scaleFactor = 1.0 / this->signalLength_;
      status = DftiSetValue( fftDescriptor, DFTI_BACKWARD_SCALE, scaleFactor);
      checkAndTrapErrors( status );
   
      if( this->stride_ != 0 )
      {
        // try to overwrite the default stride
        int strideArray[2];
        strideArray[0] = 0;        // this is the offset from the start.  We'll fix at zero
        strideArray[1] = this->stride_;   // this is the offset to the next value.
        status = DftiSetValue(fftDescriptor, DFTI_INPUT_STRIDES, strideArray);
        checkAndTrapErrors( status );
        status = DftiSetValue(fftDescriptor, DFTI_OUTPUT_STRIDES, strideArray);
        checkAndTrapErrors( status );
      }
   
      // commit it so that the library can do any needed allocations
      status = DftiCommitDescriptor( fftDescriptor );
      checkAndTrapErrors( status );
    }

    // Basic destructor 
    virtual ~N_UTL_IntelFFT_Interface() 
    {
      // free the descriptor
      long status = DftiFreeDescriptor( &fftDescriptor );
      checkAndTrapErrors( status );
    }

    // Calculate FFT with the vectors that have been registered.
    void calculateDFT();

    // Calculate IFT with the vectors that have been registered.
    void calculateIFT();

    void calculateDFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      N_UTL_DFTInterfaceDecl<VectorType>::calculateDFT( inData, outData );
    }

    // Calculate IFT with new vectors, not the ones that have been registered or used in the constructor.
    void calculateIFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      N_UTL_DFTInterfaceDecl<VectorType>::calculateIFT( inData, outData );
    }

  private:
    // Check and trap errors.
    // NOTE: The Intel Math Library returns a status code after most FFT operations. 
    //       This method checks the status code for an error signal and then prints out 
    //       the text error message if there is one.
    void checkAndTrapErrors( long fftStatus )
    {
      long classError = DftiErrorClass(fftStatus, DFTI_NO_ERROR);
      if (! classError)
      {
        std::cout << "Error in FFT operation \"";
        char* errorMessage = DftiErrorMessage(fftStatus);
        std::cout << errorMessage << "\"" << std::endl
          << "Exiting." << std::endl;
        exit(-1);
      }
    }

    // Data structor which holds info about the fft (size, dimension, etc)
    DFTI_DESCRIPTOR * fftDescriptor;
};

#endif

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
// Purpose        : This class is an implementation of the abstract FFT interface
//                  for the Apple FFT library.
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 5/27/08
//
//-------------------------------------------------------------------------
#ifndef N_UTL_APPLE_FFT_INTERFACE_HPP
#define N_UTL_APPLE_FFT_INTERFACE_HPP


// ---------- Standard Includes ----------

#include <N_UTL_FFTInterfaceDecl.hpp>
#include <Accelerate/Accelerate.h>
#include <iostream>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface for Apple FFT 
// Purpose       : This class acts as an abstract interface to an FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  
// Special Notes : See 
//                 https://developer.apple.com/documentation/accelerate/discrete_fourier_transforms?language=objc
//                 for details on Apple's implementation.  
// 
// Creator       : Rich Schiek
// Creation Date : 02/16/23
// Last Modified : 02/16/23
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_AppleFFT_Interface: public N_UTL_FFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.
    N_UTL_AppleFFT_Interface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      : N_UTL_FFTInterfaceDecl<VectorType>(length, numSignals, reqStride, overwrite)
    {
      // Note, Apple's implementation doesn't support numSignals != 1 and a reqStride !=0
      // so error out in that case.  Other Apple DFT routines can handled numSignals > 1
      // but they are limited to vector lengths of 2^12 or less as the number of signals goes up.
      // see: 
      //  https://developer.apple.com/documentation/accelerate/3746721-vdsp_dft_interleaved_createsetup?language=objc
      //
      if( (numSignals != 1) || (reqStride != 0))
      {
        Xyce::dout() << "Apple's implementation of FFT/IFT does not allow for numSignals != 1 or reqStrid != 0";
        exit(-1);
      }
      
      // for the Apple routines, vector length must be a power of two.  If it is not then
      // the routine may segfault.  So the interface class must handled vector padding 
      // when the input vector is not a power of 2 in length 
      nextLargestPowerOf2_ = length;
      double log2Length = std::log2(vecLengthProvided);
      if( std::floor( log2Length)  < log2Length )
      {
        nextLargestPowerOf2_ = (int) std::pow(2,(floor(log2Length)+1));
        paddingRequired_ = true;
      }
    
      // create forward and inverse setup objects
      vDSP_DFT_SetupD forwardSetup = vDSP_DFT_zop_CreateSetupD(NULL, nextLargestPowerOf2_, vDSP_DFT_FORWARD);
      vDSP_DFT_SetupD inverseSetup = vDSP_DFT_zop_CreateSetupD(NULL, nextLargestPowerOf2_, vDSP_DFT_INVERSE);
      
      // The forward and backward transform must have a consistent scale factor
      // so that the inverse of a forward transform is the same signal.  By default
      // we'll use 1.0 for the forward scale factor and then 1/n for the inverse transform.
      scaleFactor_ = 1.0 / this->nextLargestPowerOf2;
    }

    // Basic destructor 
    virtual ~N_UTL_AppleFFT_Interface() 
    {
      // free the setup data 
      vDSP_DFT_DestroySetupD(forwardSetup);
      vDSP_DFT_DestroySetupD(inverseSetup);
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
    // NOTE: The Apple Math Library returns a status code after most FFT operations. 
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

    // vector size and next largest power of two size if the input vector is not
    // a power of 2 in length.  In that case the calculation vectors will need 
    // to be padded.
    unsigned int nextLargestPowerOf2_;
    bool paddingRequired_;
    // Data structure which holds info about the fft (size, direction)
    vDSP_DFT_SetupD forwardSetup_;
    vDSP_DFT_SetupD inverseSetup_;
    
    // scale factor needed on back transform 
    double scaleFactor_;
};

#endif

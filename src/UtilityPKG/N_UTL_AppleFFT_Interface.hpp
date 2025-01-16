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
#define __CLAPACK_H
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
        std::cout << "Apple's implementation of FFT/IFT does not allow for numSignals != 1 or reqStrid != 0";
        exit(-1);
      }
      
      // for the Apple routines, vector length must be a power of two.  If it is not then
      // the routine may segfault.  So the interface class must handled vector padding 
      // when the input vector is not a power of 2 in length 
      nextLargestPowerOf2_ = length;
      double log2Length = std::log2(length);
      if( std::floor( log2Length)  < log2Length )
      {
        nextLargestPowerOf2_ = (int) std::pow(2,(floor(log2Length)+1));
        paddingRequired_ = true;
      }
      //std::cout << "Apple FFT length = " << length << " nextLargestPowerOf2_ = " << nextLargestPowerOf2_;
      
      interleavedRoutines_ = false;
      // force non-interleaved routines for now
      /*
      if( nextLargestPowerOf2_ <= std::pow(2,14))
      {
        interleavedRoutines_ = true;
        forwardInterleavedSetup_ = vDSP_DFT_Interleaved_CreateSetupD(NULL, nextLargestPowerOf2_, vDSP_DFT_FORWARD, vDSP_DFT_Interleaved_RealtoComplex);
        inverseInterleavedSetup_ = vDSP_DFT_Interleaved_CreateSetupD(NULL, nextLargestPowerOf2_, vDSP_DFT_INVERSE, vDSP_DFT_Interleaved_ComplextoComplex);
        //std::cout << " interleavedRoutines_ = true" << std::endl;
      }
      else
      */
      {
        // create forward and inverse setup objects
        forwardSetup_ = vDSP_DFT_zop_CreateSetupD(NULL, nextLargestPowerOf2_, vDSP_DFT_FORWARD);
        inverseSetup_ = vDSP_DFT_zop_CreateSetupD(NULL, nextLargestPowerOf2_, vDSP_DFT_INVERSE);
        //std::cout << " interleavedRoutines_ = false" << std::endl;
      }
      
      // The forward and backward transform must have a consistent scale factor
      // so that the inverse of a forward transform is the same signal.  By default
      // we'll use 1.0 for the forward scale factor and then 1/n for the inverse transform.
      scaleFactor_ = 1.0 / this->nextLargestPowerOf2_;
      
      // allocate the input/output vectors
      yTdReal.resize(nextLargestPowerOf2_);
      yTdImag.resize(nextLargestPowerOf2_);
      yFdReal.resize(nextLargestPowerOf2_);
      yFdImag.resize(nextLargestPowerOf2_);
      yTdInterleaved.resize(2*nextLargestPowerOf2_);
      yFdInterleaved.resize(2*nextLargestPowerOf2_);
      
    }

    // Basic destructor 
    virtual ~N_UTL_AppleFFT_Interface() 
    {
      if( interleavedRoutines_ )
      {
        vDSP_DFT_Interleaved_DestroySetupD(forwardInterleavedSetup_);
        vDSP_DFT_Interleaved_DestroySetupD(inverseInterleavedSetup_);
      }
      else
      {
        // free the setup data 
        vDSP_DFT_DestroySetupD(forwardSetup_);
        vDSP_DFT_DestroySetupD(inverseSetup_);
      }
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
    // vector size and next largest power of two size if the input vector is not
    // a power of 2 in length.  In that case the calculation vectors will need 
    // to be padded.
    unsigned int nextLargestPowerOf2_;
    bool paddingRequired_;
    bool interleavedRoutines_;
    // Data structure which holds info about the fft (size, direction)
    vDSP_DFT_SetupD forwardSetup_;
    vDSP_DFT_SetupD inverseSetup_;
    vDSP_DFT_Interleaved_SetupD forwardInterleavedSetup_;
    vDSP_DFT_Interleaved_SetupD inverseInterleavedSetup_; 
    
    VectorType yTdReal, yTdImag, yFdReal, yFdImag;  
    VectorType yTdInterleaved, yFdInterleaved;
    
    // scale factor needed on back transform 
    double scaleFactor_;
};

#endif

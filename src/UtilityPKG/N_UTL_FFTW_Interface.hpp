//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
//                  for the FFTW library.
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 5/27/08
//
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTW_INTERFACE_HPP
#define N_UTL_FFTW_INTERFACE_HPP


// ---------- Standard Includes ----------

#include <N_UTL_FFTInterfaceDecl.hpp>

// ----------   Other Includes   ----------

#include <fftw3.h>
#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface for FFTW
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
class N_UTL_FFTW_Interface: public N_UTL_FFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[DFT/IFT] methods.
    N_UTL_FFTW_Interface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
      : N_UTL_FFTInterfaceDecl<VectorType>(length, numSignals, reqStride, overwrite),
        firstForwardFFT_(true),
        firstInverseFFT_(true)
    {}

    // Basic destructor 
    virtual ~N_UTL_FFTW_Interface() 
    {
      if (!firstForwardFFT_)
        fftw_destroy_plan(forwardPlan_);
      if (!firstInverseFFT_)
        fftw_destroy_plan(inversePlan_);
      // Only clean up when the object is being destroyed, save any wisdom FFTW may have obtained.
      fftw_cleanup();
    }

    // Register new vectors for the FFT/IFT interface to use.
    void registerVectors( const Teuchos::RCP<VectorType>& dftInData, const Teuchos::RCP<VectorType>& dftOutData,
                          const Teuchos::RCP<VectorType>& iftInData, const Teuchos::RCP<VectorType>& iftOutData )
    { 
      this->dftInData_ = dftInData; 
      this->dftOutData_ = dftOutData; 
      this->iftInData_ = iftInData; 
      this->iftOutData_ = iftOutData; 
      if (!firstForwardFFT_)
      {
        fftw_destroy_plan(forwardPlan_);
        firstForwardFFT_ = true;
      }
      if (!firstInverseFFT_)
      {
        fftw_destroy_plan(inversePlan_);
        firstInverseFFT_ = true;
      }
    }

    // Calculate FFT with new vectors, not the ones that have been registered or used in the constructor.
    void calculateDFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      // Check if the vectors are new, if so destroy the current plan an rebuild it in the next call to calculateDFT()
      if ((inData != this->dftInData_) || (outData != this->dftOutData_))
      {
        this->dftInData_ = inData;
        this->dftOutData_ = outData;
        if (!firstForwardFFT_)
        {
          fftw_destroy_plan(forwardPlan_);
          firstForwardFFT_ = true;
        }
      }
      calculateDFT();
    }

    // Calculate IFT with new vectors, not the ones that have been registered or used in the constructor.
    void calculateIFT( const Teuchos::RCP<VectorType>& inData, const Teuchos::RCP<VectorType>& outData )
    {
      // Check if the vectors are new, if so destroy the current plan an rebuild it in the next call to calculateIFT()
      if ((inData != this->iftInData_) || (outData != this->iftOutData_))
      {
        this->iftInData_ = inData;
        this->iftOutData_ = outData;
        if (!firstInverseFFT_)
        {
          fftw_destroy_plan(inversePlan_);
          firstInverseFFT_ = true;
        }
      }
      calculateIFT();
    }

    // Calculate DFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateDFT();
    // Calculate IFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateIFT();

  private:
    bool firstForwardFFT_, firstInverseFFT_;
    Teuchos::RCP<VectorType> inDataTmp_, outResultTmp_;
    fftw_plan forwardPlan_;
    fftw_plan inversePlan_;
};

#endif

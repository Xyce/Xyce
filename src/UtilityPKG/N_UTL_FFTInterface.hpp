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
// Purpose        : This class acts as an interface to an FFT library
//                  for FFT and IFT calculations.  This class should isolate
//                  Xyce from the specifics of a given FFT library so 
//                  that multiple libraries can be used.  
//
// Special Notes  : 
//
// Creator        : Richard Schiek 
//
// Creation Date  : 5/27/08
//
//-------------------------------------------------------------------------
#ifndef N_UTL_FFTINTERFACE_HPP
#define N_UTL_FFTINTERFACE_HPP 1


// ---------- Standard Includes ----------
#include <Xyce_config.h>

#ifdef Xyce_USE_INTEL_FFT
#include <N_UTL_IntelFFT_Interface.hpp>
#endif

#ifdef Xyce_USE_FFTW
#include <N_UTL_FFTW_Interface.hpp>
#endif

#include <N_UTL_FFTInterfaceDecl.hpp>

#include <N_ERH_ErrorMgr.h>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : FFTInterface
// Purpose       : This class acts as a templated interface to any FFT library
//                 for FFT and IFT calculations.  This class should isolate
//                 Xyce from the specifics of a given FFT library so 
//                 that multiple libraries can be used.  It is originally
//                 implemented for Intel's Math Library but may be extended
//                 to FFTW at some time in the future.
// Special Notes :
// Creator       : Richard Schiek (templated by Heidi Thornquist)
// Creation Date : 5/27/08
// Last Modified : 11/17/10
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_FFTInterface
{
  public:
    N_UTL_FFTInterface( int length, int numSignals=1, int reqStride=0, bool overwrite=false )
    {
#ifdef Xyce_USE_INTEL_FFT
      fftInterface_ = Teuchos::rcp( new N_UTL_IntelFFT_Interface<VectorType>( length, numSignals, reqStride, overwrite ) );
#elif defined(Xyce_USE_FFTW)
      fftInterface_ = Teuchos::rcp( new N_UTL_FFTW_Interface<VectorType>( length, numSignals, reqStride, overwrite ) );
#else
      Xyce::Report::DevelFatal0() 
        <<  "Xyce has not been configured with an FFT library! Please reconfigure with FFT enabled to perform any frequency analysis!";
#endif
    }
    
    virtual ~N_UTL_FFTInterface() {}
   
    // Register new vectors for the FFT/IFT interface to use.
    void registerVectors( VectorType& fftInData, VectorType* fftOutData,
                          VectorType& iftInData, VectorType* iftOutData )
    {
      fftInterface_->registerVectors( Teuchos::rcp( &fftInData, false ), Teuchos::rcp( fftOutData, false ),
                                      Teuchos::rcp( &iftInData, false ), Teuchos::rcp( iftOutData, false ) );
    }

    // Return the vectors that were registered for the FFT interface to use.
    void getFFTVectors( Teuchos::RCP<VectorType>& fftInData,  Teuchos::RCP<VectorType>& fftOutData )
    {
      fftInterface_->getDFTVectors( fftInData, fftOutData );
    }

    // Return the vectors that were registered for the IFT interface to use.
    // NOTE:  Consider using the already registered vectors to avoid recreation of the FFT object.
    void getIFTVectors( Teuchos::RCP<VectorType>& iftInData,  Teuchos::RCP<VectorType>& iftOutData )
    {
      fftInterface_->getIFTVectors( iftInData, iftOutData );
    }
 
    void calculateFFT( VectorType& inData, VectorType* outResult)
    {
      fftInterface_->calculateDFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
    }
    
    void calculateIFT( VectorType& inData, VectorType* outResult)
    {
      fftInterface_->calculateIFT( Teuchos::rcp( &inData, false ), Teuchos::rcp( outResult, false ) );
    }
   
    void calculateFFT()
    {
      fftInterface_->calculateDFT();
    }
    
    void calculateIFT()
    {
      fftInterface_->calculateIFT();
    }

    Teuchos::RCP< N_UTL_FFTInterfaceDecl<VectorType> > getFFTInterface()
    {
      return fftInterface_;
    }

  private:

  // A pointer to the FFT interface used by this class to compute the forward and backward transforms.
  Teuchos::RCP< N_UTL_FFTInterfaceDecl<VectorType> > fftInterface_;
};

#endif


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
// Purpose        : This file contains specializations for the FFTW interface
//                  for various vector types.
//
// Special Notes  : 
//
// Creator        : Heidi Thornquist
//
// Creation Date  : 5/27/08
//
//
//
//
//-------------------------------------------------------------------------
// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <N_UTL_FFTW_Interface.hpp>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// FFTW Interface specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_FFTW_Interface<std::vector<double> >::calculateDFT()
  {
    // If the plan needs to be constructed do that first, then execute the plan
    if (firstForwardFFT_)
    {
      // Although we used a const on input to show that we aren't changing dftInData_
      // we need to cast that away as the FFT library takes non-const pointers.
      std::vector<double>::const_iterator inDataItr = (this->dftInData_)->begin();
      double * inDataPtr = const_cast< double * >( &(*inDataItr) );

      // We need to create a temp vector for converting the storage format 
      outResultTmp_ = Teuchos::rcp( new std::vector<double>((this->dftOutData_)->size(),0.0) );

      forwardPlan_ = fftw_plan_r2r_1d(signalLength_, inDataPtr, &(*outResultTmp_)[0], 
                                     FFTW_R2HC, FFTW_ESTIMATE );
      firstForwardFFT_ = false;
    }

    // Execute the FFT.
    fftw_execute(forwardPlan_);
 
    // Copy the data into the appropriate place in dftOutData_. 
    int n2 = (int)(signalLength_/2);
    (*(this->dftOutData_))[0] = (*outResultTmp_)[0];
    (*(this->dftOutData_))[1] = 0.0;
    for(int i=1; i<=n2; ++i)
    { 
      (*(this->dftOutData_))[2*i] = (*outResultTmp_)[i];

      if ( (i == n2) && ( signalLength_ % 2 == 0 ) )
        (*(this->dftOutData_))[2*i+1] = 0.0;
      else      
        (*(this->dftOutData_))[2*i+1] = (*outResultTmp_)[signalLength_-i];
    }
  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_FFTW_Interface<std::vector<double> >::calculateIFT()
  {
    // If the plan needs to be constructed do that first, then execute the plan

    if (firstInverseFFT_)
    {
      // We need to create a temp vector for converting the storage format 
      inDataTmp_ = Teuchos::rcp( new std::vector<double>((this->iftInData_)->size(),0.0) );

      inversePlan_ = fftw_plan_r2r_1d(signalLength_, &(*inDataTmp_)[0], &(*this->iftOutData_)[0], 
                                     FFTW_HC2R, FFTW_ESTIMATE );
      firstInverseFFT_ = false;
    }

    // Copy the data into the appropriate place in iftInData_. 
    int n2 = (int)(signalLength_/2);
    (*inDataTmp_)[0] = (*(this->iftInData_))[0];
    for(int i=1; i<=n2; ++i)
    {
      (*inDataTmp_)[i] = (*(this->iftInData_))[2*i];

      if ( !( (i == n2) && (signalLength_ % 2 == 0 ) ) )
        (*inDataTmp_)[signalLength_-i] = (*(this->iftInData_))[2*i+1];
    }

    // Execute the IFT.
    fftw_execute(inversePlan_);

    // Scale the output by "n"
    for (unsigned int i=0; i<this->iftOutData_->size(); ++i)
    {
      (*(this->iftOutData_))[i] /= signalLength_;
    }
  }

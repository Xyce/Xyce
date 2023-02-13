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
// Purpose        : This file contains specializations for the Intel FFT interface
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

#include <N_UTL_IntelFFT_Interface.hpp>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Intel FFT Interface specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_IntelFFT_Interface<std::vector<double> >::calculateDFT()
  {
    // Although we used a const on input to show that we aren't changing dftInData_
    // we need to cast that away as the FFT library takes non-const pointers.
    std::vector<double>::const_iterator inDataItr = (this->dftInData_)->begin();
    double * inDataPtr = const_cast< double * >( &(*inDataItr) );
    std::vector<double>::iterator outResultItr = (this->dftOutData_)->begin();
    double * outResultPtr = &(*outResultItr);

    long status = DftiComputeForward( fftDescriptor, inDataPtr, outResultPtr);
    checkAndTrapErrors( status );
  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_IntelFFT_Interface<std::vector<double> >::calculateIFT()
  {
    // Although we used a const on input to show that we aren't changing iftInData_
    // we need to cast that away as the FFT library takes non-const pointers.
    std::vector<double>::const_iterator inDataItr = (this->iftInData_)->begin();
    double * inDataPtr = const_cast< double * >( &(*inDataItr) );
    std::vector<double>::iterator outResultItr = (this->iftOutData_)->begin();
    double * outResultPtr = &(*outResultItr);

    long status = DftiComputeBackward( fftDescriptor, inDataPtr, outResultPtr);
    checkAndTrapErrors( status );
  }

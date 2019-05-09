//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : This file contains specializations for the APFT 
//                  for various vector types.
//
// Special Notes  : 
//
// Creator        : Ting Mei
//
// Creation Date  : 4/9/14
//
//-------------------------------------------------------------------------
// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <N_LAS_BlockVector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <N_UTL_APFT.h>

// ----------   Other Includes   ----------

#include <iostream>
#include <vector>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
//#include <Teuchos_SerialDenseHelpers.hpp>
//#include <Teuchos_SerialDenseSolver.hpp>
// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// APFT specialization for std::vector
//-----------------------------------------------------------------------------

  template<>
  void N_UTL_APFT<std::vector<double> >::calculateDFT()
  {
    int signalLength_ = dftInData_->size();

    Teuchos::SerialDenseVector<int,double> ftInDataVector ( Teuchos::View, &(*dftInData_)[0], signalLength_);
    Teuchos::SerialDenseVector<int,double> ftOutDataVector( Teuchos::View, &(*dftOutData_)[1], signalLength_);

    ftOutDataVector.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, dftMatrix_, ftInDataVector, 0.0 );

     (*dftOutData_)[0] = (*dftOutData_)[1];
     (*dftOutData_)[1] = 0.0;

     ftOutDataVector.scale(1.0/2.0);

     for ( int i=1; i< (signalLength_ + 1)/2; i++)
     {
       (*dftOutData_)[2*i+1] = -(*dftOutData_)[2*i+1]; 
     }
//     std::cout << "The output of the calculateFFT() is:" << std::endl;
//     ftOutDataVector.print(std::cout);
//     std::cout << "The dftOutData_ is:"  << std::endl;
//     for (int i=0; i<= signalLength_; i++)
//     {
//       std::cout << (*dftOutData_)[i] << std::endl;
//     }

  }

  // Calculate IFT with the vectors that have been registered.
  template<>
  void N_UTL_APFT<std::vector<double> >::calculateIFT()
  {

    int inSignalLength_ = iftInData_->size() - 1;
    Teuchos::SerialDenseVector<int,double> iftInDataVector ( Teuchos::View, &(*iftInData_)[1],  inSignalLength_);
    int outSignalLength_ = iftOutData_->size();
    Teuchos::SerialDenseVector<int,double> iftOutDataVector( Teuchos::View, &(*iftOutData_)[0],  outSignalLength_);

    iftInDataVector.scale(2.0);

    for ( int i=1; i< ( inSignalLength_ + 1)/2; i++)
    {
      (*iftInData_)[2*i+1] = - (*iftInData_)[2*i+1]; 
    }

    iftInDataVector[0] = (*iftInData_)[0];

    iftOutDataVector.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, idftMatrix_, iftInDataVector, 0.0 );

 /*    std::cout << "The output of the calculateIFT() is:" << std::endl;
    iftOutDataVector.print(std::cout);
    std::cout << "The iftOutData_ is:"  << std::endl;

    for (int i=0; i<= signalLength_; i++)
    {
      std::cout << (*iftOutData_)[i] << std::endl;
    } */

  }

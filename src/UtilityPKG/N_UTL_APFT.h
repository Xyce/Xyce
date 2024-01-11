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
// Purpose        : This class performs DFT and IFT for APFT algorithm
//
// Special Notes  : 

//
// Creator        : Ting Mei
//
// Creation Date  : 4/09/14
//
//-------------------------------------------------------------------------
#ifndef N_UTL_APFT_H
#define N_UTL_APFT_H


// ---------- Standard Includes ----------

#include <N_UTL_FFTInterfaceDecl.hpp>

// ----------   Other Includes   ----------

#include <Teuchos_RCP.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>

// ---------- Structure definitions ----------

//-----------------------------------------------------------------------------
// Class         : N_UTL_APFT
// Purpose       : This class performs DFT and IFT 
// Special Notes :
// Creator       : Ting Mei
// Creation Date : 4/15/14
// Last Modified : 4/15/14
//-----------------------------------------------------------------------------
template<typename VectorType>
class N_UTL_APFT: public N_UTL_DFTInterfaceDecl<VectorType>
{
  public:
    // Basic constructor without passing in vector structures, which need to be registered later or
    // passed into the calculate[FFT/IFT] methods.

    N_UTL_APFT(const Teuchos::SerialDenseMatrix<int,double>& idftMatrix, const Teuchos::SerialDenseMatrix<int,double>& dftMatrix )
    {

      idftMatrix_ = idftMatrix;
      dftMatrix_ = dftMatrix;
     
//      std::cout << "checking the IDFT matrix in FT interface"  << std::endl;
//      idftMatrix_.print(std::cout);

//      std::cout << "checking the dftmatrix in FT interface"  << std::endl;
//      dftMatrix_.print(std::cout); 

    }

    // Basic destructor 
    virtual ~N_UTL_APFT() {}

    // Calculate FFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
    void calculateDFT();
    // Calculate IFT with the vectors that have been registered.
    // NOTE:  This method must be specialized for each type of vector used by this class,
    //        or the lack of method definition will result in a build failure.
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

    int getScalar()
    {return 1;}

    void setSignalLength(int length){ }
    
  private:

  // Fourier matrices
  Teuchos::SerialDenseMatrix<int,double> idftMatrix_, dftMatrix_;

};

#endif

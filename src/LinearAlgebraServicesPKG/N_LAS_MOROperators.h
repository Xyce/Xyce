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

//-----------------------------------------------------------------------------
//
// Purpose        : Specification file for operators that are necessary for 
//                  performing model-order reduction.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL
//
// Creation Date  : 06/04/12
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_MOROperators_h
#define Xyce_N_LAS_MOROperators_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_Solver.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Matrix.h>
#include <N_PDS_ParMap.h>

// ----------  Other Includes   ----------

// Include header for Epetra compressed-row storage matrix and linear problem
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h" 

// Include selected communicator class and map required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : MORGenOp
// Purpose       : Implementation of the Epetra_Operator interface to define
//                 inv(A)*B, where inv(A) is computed by Amesos.
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/04/12
//-----------------------------------------------------------------------------
class MORGenOp : public virtual Epetra_Operator
{
public:
  // Basic constructor
  MORGenOp( const Teuchos::RCP<Linear::Solver>& solver,
            const Teuchos::RCP<Linear::Matrix>& B,
            bool reuseFactor = true, bool useTranspose = false );
  // Destructor
  ~MORGenOp() {};

  // Methods for supporting Epetra_Operator interface
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;
  const char* Label() const { return "Amesos direct solver for applying A^{-1}B"; }
  bool UseTranspose() const { return useTranspose_; }
  int SetUseTranspose( bool useTranspose ) { useTranspose_ = useTranspose; return 0; }
  const Epetra_Comm& Comm() const { return B_->Comm(); };
  const Epetra_Map& OperatorDomainMap() const { return B_->OperatorDomainMap(); }
  const Epetra_Map& OperatorRangeMap() const { return B_->OperatorRangeMap(); }

  // Epetra_Operator interface methods that are not supported.
  // Note:  ApplyInverse not defined because M not guaranteed to have an inverse.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const { return -1; }
  bool HasNormInf() const { return false; };
  double NormInf() const { return -1.0; };
   
private:  
  // Default constructor
  MORGenOp () {};
  
  // Copy constructor 
  MORGenOp ( const MORGenOp& genOp ) {};
  
  // Epetra_LinearProblem contained in the Amesos_BaseSolver
  bool reuseFactors_, useTranspose_;
  Teuchos::RCP<Linear::Solver> solver_;
  Teuchos::RCP<Epetra_Operator> B_;
  Teuchos::RCP<Epetra_LinearProblem> problem_;
  
};

// Generate basis vectors for K(inv(G + s0*C)*C, R), where Op = inv(G + s0*C)*C
Teuchos::RCP<const Linear::MultiVector> createKrylovBasis( const Teuchos::RCP<Linear::MORGenOp>& Op,
                                                           const Teuchos::RCP<Linear::MultiVector>& R,
                                                           int numBlocks, int blockSize );

// Get a multivector that views the first k columns of the input multivector V.
Linear::MultiVector* cloneView( Linear::MultiVector* V, int k );

const Linear::MultiVector* cloneView( const Linear::MultiVector* V, int k );

Linear::MultiVector* transferSDMtoMV( Parallel::ParMap& sdm_map, Teuchos::SerialDenseMatrix<int, double>& sdm );

// Perform V'*W
void blockDotProduct( const Linear::MultiVector& V, const Linear::MultiVector& W, Teuchos::SerialDenseMatrix<int,double>& result );

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_MOROperators_h

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
// Purpose       : Implementation of the Epetra_Operator interface to define
//                 inv(A)*B, where inv(A) is computed by Amesos.
// Special Notes :
//
// Creator       : Heidi Thornquist, SNL
//
// Creation Date : 06/04/12
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <algorithm>

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_MOROperators.h>
#include <N_LAS_Problem.h>

#include <N_LAS_EpetraProblem.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_EpetraMultiVector.h>
#include <N_PDS_EpetraParMap.h>

// ----------  Other Includes   ----------

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresIter.hpp>
#include <BelosDGKSOrthoManager.hpp>
#include <BelosStatusTestMaxIters.hpp>
#include <BelosOutputManager.hpp>
#include <BelosEpetraAdapter.hpp>

#include <Epetra_MultiVector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : MORGenOp::MORGenOp
// Purpose       : Constructor 
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 06/04/12
//-----------------------------------------------------------------------------
MORGenOp::MORGenOp( const Teuchos::RCP<Linear::Solver>& solver,
                    const Teuchos::RCP<Linear::Matrix>& B,
                    bool reuseFactors, bool useTranspose )
  : reuseFactors_(reuseFactors),
    useTranspose_(useTranspose),
    solver_(solver)
{
  Linear::Matrix& lasMatrix = const_cast<Linear::Matrix&>( *B );
  B_ = Teuchos::rcp( &((dynamic_cast<EpetraMatrix&>(lasMatrix)).epetraObj()), false );

  Linear::Problem& lasProblem = const_cast<Linear::Problem&>( solver->getProblem() );
  problem_ = Teuchos::rcp( &((dynamic_cast<EpetraProblem&>(lasProblem)).epetraObj()), false );

  if ( B_->UseTranspose() )
    B_->SetUseTranspose(!useTranspose);
  else
    B_->SetUseTranspose(useTranspose);
}

//-----------------------------------------------------------------------------
// Function      : MORGenOp::Apply()
// Purpose       : Applies the operator inv(A)*B*X = Y
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL 
// Creation Date : 06/04/12
//-----------------------------------------------------------------------------
int MORGenOp::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{ 
  if (!useTranspose_) {

    // Storage for B*X 
    Epetra_MultiVector BX(X.Map(),X.NumVectors());

    // Apply B*X
    B_->Apply(X, BX);
    Y.PutScalar(0.0);

    // Set the LHS and RHS
    problem_->SetRHS(&BX);
    problem_->SetLHS(&Y);

    // Solve the linear system A*Y = BX
    solver_->solve( reuseFactors_ );
  }
  else {
    // Storage for A^{-T}*X
    Epetra_MultiVector ATX(X.Map(),X.NumVectors());
    Epetra_MultiVector tmpX = const_cast<Epetra_MultiVector&>(X);

    // Set the LHS and RHS
    problem_->SetRHS(&tmpX);
    problem_->SetLHS(&ATX);

    // Solve the linear system A^T*Y = X 
    solver_->solveTranspose( reuseFactors_ );

    // Apply B*ATX
    B_->Apply(ATX, Y);
  }
  
  return 0;
}

// Generate basis vectors for K(inv(G + s0*C)*C, R), where Op = inv(G + s0*C)*C
Teuchos::RCP<const Linear::MultiVector> createKrylovBasis( const Teuchos::RCP<Linear::MORGenOp>& Op,
                                                           const Teuchos::RCP<Linear::MultiVector>& R,
                                                           int numBlocks, int blockSize )
{
  // ---------------------------------------------------------------------
  // Now use Belos to compute the basis vectors for K_k(inv(G + s0*C)*C, R)
  // ---------------------------------------------------------------------

  // Helpful typedefs for the templates
  typedef double                            ST;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>      MVT;

  // Output manager.
  Belos::OutputManager<ST> printer;

  // Status test.
  Belos::StatusTestMaxIters<ST, MV, OP> maxIterTest( numBlocks );

  // Orthogonalization manager.
  Belos::DGKSOrthoManager<ST, MV, OP> orthoMgr;

  // Linear Problem.

  // Reuse R.  We need the basis vectors for K(inv(G + s0*C)*C, R)
  Teuchos::RCP<Linear::EpetraMultiVector> EMV = Teuchos::rcp_dynamic_cast<Linear::EpetraMultiVector>( R );
  Linear::EpetraMultiVector temp( *(R->pmap()), blockSize );
  Belos::LinearProblem<ST, MV, OP > problem( Op,
                                             Teuchos::rcp( &temp.epetraObj(), false ),
                                             Teuchos::rcp( &EMV->epetraObj(), false ));
  problem.setProblem();

  // Create parameter list.
  Teuchos::ParameterList params;
  params.set("Num Blocks", numBlocks);
  params.set("Block Size", blockSize);

  // Create Krylov subspace iteration from Belos
  Belos::BlockGmresIter<ST, MV, OP> krylovIter( Teuchos::rcp( &problem, false ),
                                                Teuchos::rcp( &printer, false ),
                                                Teuchos::rcp( &maxIterTest, false ),
                                                Teuchos::rcp( &orthoMgr, false ),
                                                params );

  // Get a matrix to hold the orthonormalization coefficients.
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ST> > z
    = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ST>( blockSize, blockSize) );

  // Orthonormalize the Krylov kernel vectors V
  Teuchos::RCP<MV> V = Teuchos::rcp( new Epetra_MultiVector( EMV->epetraObj() ) );
  orthoMgr.normalize( *V, z );

  // Set the new state and initialize the solver to use R as the kernel vectors.
  Belos::GmresIterationState<ST,MV> initState;
  initState.V = V;
  initState.z = z;
  initState.curDim = 0;
  krylovIter.initializeGmres(initState);

  // Have the solver iterate until the basis size is computed
  try {
    krylovIter.iterate();
  }
  catch (const Belos::GmresIterationOrthoFailure &e) {
    // This might happen if the basis size is the same as the operator's dimension.
  }
  catch (const std::exception &e) {
    return Teuchos::null;
  }

  // Get the basis vectors back from the iteration object
  Belos::GmresIterationState<ST,MV> newState = krylovIter.getState();

  // Return a copy of the Krylov vectors
  Epetra_MultiVector* copyV = new Epetra_MultiVector( *(newState.V) );
  Teuchos::RCP<const Linear::MultiVector> retV = 
    Teuchos::rcp( new Linear::EpetraMultiVector( copyV, true ) );

  return retV;
}

Linear::MultiVector* cloneView( Linear::MultiVector* V, int k )
{
  // Cast to an Epetra-aware object
  Linear::EpetraMultiVector* EMV = dynamic_cast<Linear::EpetraMultiVector *>( V );

  std::vector<int> indices(k);
  for (int i=0; i<k; ++i) { indices[i] = i; }

  Linear::EpetraMultiVector* newV =  new Linear::EpetraMultiVector( new Epetra_MultiVector( View, EMV->epetraObj(), &indices[0], k), true );

  return newV;
}

const Linear::MultiVector* cloneView( const Linear::MultiVector* V, int k )
{
  // Cast to an Epetra-aware object
  const Linear::EpetraMultiVector* EMV = dynamic_cast<const Linear::EpetraMultiVector *>( V );

  std::vector<int> indices(k);
  for (int i=0; i<k; ++i) { indices[i] = i; }

  const Linear::EpetraMultiVector* newV =  new Linear::EpetraMultiVector( new Epetra_MultiVector( View, EMV->epetraObj(), &indices[0], k), true );

  return newV;
}

Linear::MultiVector* transferSDMtoMV( Parallel::ParMap& sdm_map, Teuchos::SerialDenseMatrix<int, double>& sdm )
{
  // Create a multivector with the SerialDenseMatrix and ParMap
  Parallel::EpetraParMap& e_mapPtr = dynamic_cast<Parallel::EpetraParMap&>( sdm_map );

  Epetra_MultiVector* new_eMV = new Epetra_MultiVector( View, *(e_mapPtr.petraMap()), sdm.values(), sdm.stride(), sdm.numCols() );
  Linear::EpetraMultiVector* newMV = new Linear::EpetraMultiVector( new_eMV, true );

  return newMV;
}

// Perform V'*W
void blockDotProduct( const Linear::MultiVector& V, const Linear::MultiVector& W, 
                      Teuchos::SerialDenseMatrix<int,double>& result )
{
  const Linear::EpetraMultiVector& eV = dynamic_cast<const Linear::EpetraMultiVector&>(V);
  const Linear::EpetraMultiVector& eW = dynamic_cast<const Linear::EpetraMultiVector&>(W);

  // Helpful typedefs for the templates
  typedef double                            ST;
  typedef Epetra_MultiVector                MV;
  typedef Belos::MultiVecTraits<ST,MV>      MVT;

  MVT::MvTransMv( 1.0, eV.epetraObj(), eW.epetraObj(), result );
}

} // namespace Linear
} // namespace Xyce

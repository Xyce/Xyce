//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

//
// test the N_LAS_BlockSystemHelpers
//
#include <N_LAS_EpetraGraph.h>
#include <N_LAS_EpetraMultiVector.h>
#include <N_LAS_EpetraBlockMatrix.h>
#include <N_LAS_EpetraBlockVector.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_EpetraVector.h>

#include <N_PDS_EpetraParMap.h>
#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_EpetraMPIComm.h>
#include <mpi.h>
#else
#include <N_PDS_EpetraSerialComm.h>
#endif

#include <N_UTL_Math.h>
#include <N_ERH_Message.h>

#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MultiVectorIn.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>

#include <Amesos.h>

#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  int myPID = 0;
#ifdef Xyce_PARALLEL_MPI  
  // Initialize MPI (inside the MPIComm creation) 
  Xyce::Parallel::EpetraMPIComm comm( argc, argv );
  myPID = comm.procID();
#else
  Xyce::Parallel::EpetraSerialComm comm;
#endif

  // -----------------------------------------------------
  // Read in example HB base matrix.
  // -----------------------------------------------------

  // The memory for these objects will be deleted by the Matrix, Vector objects.
  Epetra_CrsMatrix *A = 0;
  Epetra_MultiVector *b = 0;

  // Get the filenames
  std::string matrixName = "hb_matrix1.mm";
  std::string rhsName = "hb_rhs1.mm";

  // Read the first matrix in.
  int finfo = EpetraExt::MatrixMarketFileToCrsMatrix( matrixName.c_str(), *(comm.petraComm()), A );
  if (finfo!=0 && myPID==0)
    std::cout << "First matrix file could not be read in!!!, info = "<< finfo << std::endl;
  //A->SetTracebackMode(2); // Shutdown Epetra Warning tracebacks

  Epetra_Map const& petraBaseMap = A->RowMap();
  Epetra_BlockMap const& bMap = dynamic_cast<Epetra_BlockMap const&>(petraBaseMap);

  finfo = EpetraExt::MatrixMarketFileToMultiVector( rhsName.c_str(), bMap, b );
  if (finfo!=0 && myPID==0)
    std::cout << "First rhs file could not be read in!!!, info = "<< finfo << std::endl;

  // -----------------------------------------------------
  // First test non-overlapping (ghosted) maps.
  // -----------------------------------------------------

  int numBlocks = 2;
  int numMyElements = petraBaseMap.NumMyElements();
  Xyce::Parallel::EpetraParMap baseMap( const_cast<Epetra_Map*>(&petraBaseMap), comm );

  // Create a block map for a map not including ground nodes.
  Teuchos::RCP<Xyce::Parallel::ParMap> blockMap = Xyce::Linear::createBlockParMap( numBlocks, baseMap ); 

  // Create a block vector
  Xyce::Linear::EpetraBlockVector blockVector( numBlocks, blockMap, Teuchos::rcpFromRef(baseMap) );

  // Create a block of vectors using the baseMap
  std::vector<Teuchos::RCP<Xyce::Linear::Vector> > blockVectors( numBlocks );
  for (int i=0; i<numBlocks; ++i) {
    blockVectors[i] = Teuchos::rcp( new Xyce::Linear::EpetraVector( baseMap ) );
    blockVectors[i]->putScalar( i+1.0 );  // Initialize each vector
  }

  // Copy the block of vectors into the BlockVector object.
  Xyce::Linear::copyToBlockVector( blockVectors, blockVector );

  blockVector.putScalar( numBlocks+1.0 );

  // Copy the BlockVector back into the block of vectors.
  Xyce::Linear::copyFromBlockVector( blockVector, blockVectors );
 
  // -----------------------------------------------------
  // Now test block graphs.
  // -----------------------------------------------------

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  int offset = Xyce::Linear::generateOffset( baseMap );

  Xyce::Linear::EpetraGraph baseGraph( Teuchos::rcp( const_cast<Epetra_CrsGraph*>(&(A->Graph())), false ) );

  Teuchos::RCP<Xyce::Linear::Graph> blockGraph = 
    Xyce::Linear::createBlockGraph( offset, blockPattern, *blockMap, baseGraph );
  
  // -----------------------------------------------------
  // Now test block matrices.
  // -----------------------------------------------------

  Teuchos::RCP<Xyce::Linear::BlockMatrix> blockMatrix = 
    Teuchos::rcp ( new Xyce::Linear::EpetraBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), &baseGraph ) );
  Teuchos::RCP<Xyce::Linear::EpetraBlockMatrix> eBlockMatrix =
    Teuchos::rcp_dynamic_cast<Xyce::Linear::EpetraBlockMatrix>( blockMatrix );

  // Insert the linear system on the diagonals and zeros on the off diagonal blocks
  Xyce::Linear::EpetraMatrix origMatrix( A );

  // First diagonal block
  blockMatrix->put( 0.0 ); // Zero out whole matrix
  blockMatrix->block( 0, 0 ).add( origMatrix );
  // Second diagonal block
  blockMatrix->block( 1, 1 ).add( origMatrix );

  // Create block vector RHS.
  Xyce::Linear::EpetraVector origRHS( dynamic_cast<Epetra_Vector *>(b), false );
  Xyce::Linear::EpetraVector origSoln( baseMap );

  // Call update directly on the returned block
  // Don't save to temporary Vector object, because it would only be a copy.
  blockVector.putScalar( 0.0 );
  blockVector.block( 0 ).update( 1.0, origRHS );
  blockVector.block( 1 ).update(-1.0, origRHS );;

  // -----------------------------------------------------
  // Solve the block system using Amesos.
  // -----------------------------------------------------
 
  Amesos amesosFactory;

  Xyce::Linear::EpetraBlockVector blockSoln( numBlocks, blockMap, Teuchos::rcpFromRef(baseMap) );
  blockSoln.putScalar( 0.0 );
  origSoln.putScalar( 0.0 );

  // Create the original problem and then the block problem
  Epetra_LinearProblem problem( &origMatrix.epetraObj(), &origSoln.epetraObj(), &origRHS.epetraObj() );
  Epetra_LinearProblem blockProblem( &eBlockMatrix->epetraObj(), &blockSoln.epetraObj(), &blockVector.epetraObj() );

  // Solve the original problem

  Teuchos::RCP<Amesos_BaseSolver> solver = Teuchos::rcp( amesosFactory.Create( "Klu", problem ) );

  int linearStatus = solver->SymbolicFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;

  linearStatus = solver->NumericFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
  
  linearStatus = solver->Solve();
  if (linearStatus != 0) Xyce::dout() << "Amesos solve exited with error: " << linearStatus;

  origSoln.print(Xyce::dout());

  // Solve the block problem

  Teuchos::RCP<Amesos_BaseSolver> blockSolver = Teuchos::rcp( amesosFactory.Create( "Klu", blockProblem ) );

  linearStatus = blockSolver->SymbolicFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;

  linearStatus = blockSolver->NumericFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
  
  linearStatus = blockSolver->Solve();
  if (linearStatus != 0) Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
 
  blockSoln.print(Xyce::dout());
  
  return 0;

}

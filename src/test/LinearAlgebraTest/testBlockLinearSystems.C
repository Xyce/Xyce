//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
#include <N_LAS_Graph.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_Vector.h>

#include <N_PDS_ParMap.h>
#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_MPIComm.h>
#include <mpi.h>
#else
#include <N_PDS_SerialComm.h>
#endif

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
#include <cmath>

int main(int argc, char* argv[])
{
  int myPID = 0;
#ifdef Xyce_PARALLEL_MPI  
  // Initialize MPI (inside the MPIComm creation) 
  N_PDS_MPIComm comm( argc, argv );
  myPID = comm.procID();
#else
  N_PDS_SerialComm comm;
#endif

  // -----------------------------------------------------
  // Read in example HB base matrix.
  // -----------------------------------------------------

  // The memory for these objects will be deleted by the N_LAS_Matrix, N_LAS_Vector objects.
  Epetra_CrsMatrix *A = 0;
  Epetra_MultiVector *b = 0;

  // Get the filenames
  std::string matrixName = "hb_matrix1.mm";
  std::string rhsName = "hb_rhs1.mm";

  // Read the first matrix in.
  int finfo = EpetraExt::MatrixMarketFileToCrsMatrix( matrixName.c_str(), *(comm.petraComm()), A );
  if (finfo!=0 && myPID==0)
    cout << "First matrix file could not be read in!!!, info = "<< finfo << endl;
  //A->SetTracebackMode(2); // Shutdown Epetra Warning tracebacks

  Epetra_Map const& petraBaseMap = A->RowMap();
  Epetra_BlockMap const& bMap = dynamic_cast<Epetra_BlockMap const&>(petraBaseMap);

  finfo = EpetraExt::MatrixMarketFileToMultiVector( rhsName.c_str(), bMap, b );
  if (finfo!=0 && myPID==0)
    cout << "First rhs file could not be read in!!!, info = "<< finfo << endl;

  // -----------------------------------------------------
  // First test non-overlapping (ghosted) maps.
  // -----------------------------------------------------

  int numBlocks = 2;
  int numMyElements = petraBaseMap.NumMyElements();
  N_PDS_ParMap baseMap( const_cast<Epetra_Map*>(&petraBaseMap), &comm );

  // Create a block map for a map not including ground nodes.
  Teuchos::RCP<N_PDS_ParMap> blockMap = createBlockParMap( numBlocks, baseMap ); 

  // Create a block vector
  N_LAS_BlockVector blockVector( numBlocks, *(blockMap->petraMap()), const_cast<Epetra_Map&>(petraBaseMap) );

  // Create a block of vectors using the baseMap
  std::vector<Teuchos::RCP<N_LAS_Vector> > blockVectors( numBlocks );
  for (int i=0; i<numBlocks; ++i) {
    blockVectors[i] = Teuchos::rcp( new N_LAS_Vector( baseMap ) );
    blockVectors[i]->putScalar( i+1.0 );  // Initialize each vector
  }

  // Copy the block of vectors into the N_LAS_BlockVector object.
  copyToBlockVector( blockVectors, blockVector );

  blockVector.putScalar( numBlocks+1.0 );

  // Copy the N_LAS_BlockVector back into the block of vectors.
  copyFromBlockVector( blockVector, blockVectors );
 
  // -----------------------------------------------------
  // First test overlapping (ghosted) maps.
  // -----------------------------------------------------

  // Determine overlap map using column map
  Epetra_Map const& petraColMap = A->ColMap();
  int numMyCols = petraColMap.NumMyElements();

  std::vector<int> baseGIDs( numMyCols+1 );
  petraColMap.MyGlobalElements( &baseGIDs[0] );
  baseGIDs[numMyCols] = -1;

  Epetra_Map oPetraBaseMap(
        -1, // Global elements
        numMyCols+1, // Local elements
        &baseGIDs[0],
        -1, // 0 or 1
        *(comm.petraComm()) // communicator
        );
  N_PDS_ParMap oBaseMap( &oPetraBaseMap, &comm );
  baseMap.print(Xyce::dout();
  oBaseMap.print(Xyce::dout());

  // Create block maps for a map including ground nodes.
  std::vector<Teuchos::RCP<N_PDS_ParMap> > blockMaps = 
    createBlockParMaps( numBlocks, baseMap, oBaseMap );
  blockMaps[0]->print(Xyce::dout());
  blockMaps[1]->print(Xyce::dout());

  Xyce::dout() << "CREATING NEW BLOCK MAPS!!!" << std::endl;
  std::vector<Teuchos::RCP<N_PDS_ParMap> > blockMaps2 = 
    createBlockParMaps2( numBlocks, baseMap, oBaseMap );
  blockMaps2[0]->print(Xyce::dout());
  blockMaps2[1]->print(Xyce::dout());

  // Create a block of vectors using the baseMap and oBaseMap
  std::vector<Teuchos::RCP<N_LAS_Vector> > blockVectors2( numBlocks );
  for (int i=0; i<numBlocks; ++i) {
    blockVectors2[i] = Teuchos::rcp( new N_LAS_Vector( baseMap, oBaseMap ) );
    blockVectors2[i]->putScalar( i+1.0 );  // Initialize each vector
  }

  // Create a block vector with ground node
  N_LAS_BlockVector blockVector2 = N_LAS_BlockVector( numBlocks,
                                                      blockMaps[0], 
                                                      blockMaps[1], 
                                                      Teuchos::rcp( &baseMap, false ),
                                                      Teuchos::rcp( &oBaseMap, false )
                                                    );

  // Copy the block of vectors into the N_LAS_BlockVector object.
  copyToBlockVector( blockVectors2, blockVector2 );

  blockVector2.putScalar( numBlocks+1.0 );

  // Copy the N_LAS_BlockVector back into the block of vectors.
  copyFromBlockVector( blockVector2, blockVectors2 );

  // -----------------------------------------------------
  // Now test block graphs.
  // -----------------------------------------------------

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  int MaxGID = baseMap.maxGlobalEntity();
  int offset=1;
  while ( offset <= MaxGID ) offset *= 10;

  Linear::Graph baseGraph( Teuchos::rcp( &(A->Graph()), false ) );
  Teuchos::RCP<Linear::Graph> blockGraph = createBlockGraph( offset, blockPattern, *blockMaps[0], baseGraph );
  
  // -----------------------------------------------------
  // Now test block matrices.
  // -----------------------------------------------------

  Teuchos::RCP<N_LAS_BlockMatrix> blockMatrix = Teuchos::rcp ( new N_LAS_BlockMatrix( numBlocks, blockPattern, *blockGraph, baseGraph );

  // Insert the linear system on the diagonals and zeros on the off diagonal blocks
  N_LAS_Matrix origMatrix( A );

  // First diagonal block
  blockMatrix->put( 0.0 ); // Zero out whole matrix
  blockMatrix->block( 0, 0 ).add( origMatrix );
  // Second diagonal block
  blockMatrix->block( 1, 1 ).add( origMatrix );

  // Create block vector RHS.
  N_LAS_Vector origRHS( dynamic_cast<Epetra_Vector *>(b) );
  N_LAS_Vector origSoln( origRHS );

  // Call addVec directly on the returned block
  // Don't save to temporary N_LAS_Vector object, because it would only be a copy.
  blockVector.putScalar( 0.0 );
  blockVector.block( 0 ).addVec( 1.0, origRHS );
  blockVector.block( 1 ).addVec(-1.0, origRHS );;

  // -----------------------------------------------------
  // Solve the block system using Amesos.
  // -----------------------------------------------------
 
  Amesos amesosFactory;

  N_LAS_BlockVector blockSoln( numBlocks, *(blockMap->petraMap()), const_cast<Epetra_Map&>(petraBaseMap) );
  blockSoln.putScalar( 0.0 );
  origSoln.putScalar( 0.0 );

  // Create the original problem and then the block problem
  Epetra_LinearProblem problem( &origMatrix.epetraObj(), &origSoln.epetraObj(), &origRHS.epetraObj() );
  Epetra_LinearProblem blockProblem( &blockMatrix->epetraObj(), &blockSoln.epetraObj(), &blockVector.epetraObj() );

  // Solve the original problem

  Teuchos::RCP<Amesos_BaseSolver> solver = Teuchos::rcp( amesosFactory.Create( "Klu", problem ) );

  int linearStatus = solver->SymbolicFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;

  linearStatus = solver->NumericFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
  
  linearStatus = solver->Solve();
  if (linearStatus != 0) Xyce::dout() << "Amesos solve exited with error: " << linearStatus;

  origSoln.printPetraObject();

  // Solve the block problem

  Teuchos::RCP<Amesos_BaseSolver> blockSolver = Teuchos::rcp( amesosFactory.Create( "Klu", blockProblem ) );

  linearStatus = blockSolver->SymbolicFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;

  linearStatus = blockSolver->NumericFactorization();
  if (linearStatus != 0) Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
  
  linearStatus = blockSolver->Solve();
  if (linearStatus != 0) Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
 
  blockSoln.printPetraObject();
  
  return 0;

}

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
// Purpose        : Implementation file for the Iterative linear solver
//                  interface.
//
// Special Notes  :
//
// Creator        : Heidi K. Thornquist, SNL, Electrical & Microsystem Modeling
//
// Creation Date  : 11/11/08
//
//
//
//
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

#include <Xyce_config.h>

#include <sstream>

// ----------   Xyce Includes   ----------


#include <N_LAS_HBBlockJacobiPrecond.h>
#include <N_LAS_HBBlockJacobiEpetraOperator.h>
#include <N_LAS_HBBuilder.h>
#include <N_LOA_HBLoader.h>
#include <N_MPDE_State.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Graph.h>
#include <N_LAS_FilteredMatrix.h>

#include <N_LAS_Problem.h>
#include <N_LAS_Builder.h>
#include <N_LAS_System.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LOA_Loader.h>

#include <N_UTL_Timer.h>
#include <N_UTL_Math.h>
#include <N_UTL_FeatureTest.h>

#include <N_ERH_ErrorMgr.h>

#include <N_PDS_EpetraParMap.h>
#include <N_PDS_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Comm.h>
#include <Amesos.h>

#ifdef Xyce_PARALLEL_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#endif

using Teuchos::RCP;
using Teuchos::rcp;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::HBBlockJacobiPrecond
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
HBBlockJacobiPrecond::HBBlockJacobiPrecond(
  Linear::Builder &builder)
  : Preconditioner(),
    isCorrected_(false),
    hbOsc_(false),
    builder_(builder)
{
  setDefaultOptions();
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::setDefaultOptions
// Purpose       : resets options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiPrecond::setDefaultOptions()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::setOptions
// Purpose       : sets options and params from modelblock
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiPrecond::setOptions( const Util::OptionBlock & OB )
{
  // Set the parameters from the list
  Util::ParamList::const_iterator it_tpL = OB.begin();
  Util::ParamList::const_iterator end_tpL = OB.end();
  for (; it_tpL != end_tpL; ++it_tpL)
    {
      this->setParam( *it_tpL );
    }

  // store for restart of solver_
  options_ = Teuchos::rcp( &OB, false );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::setParam
// Purpose       : sets options
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiPrecond::setParam( const Util::Param & param )
{
  std::string uTag = param.uTag();

  // Check if the one-step correction should be used.
  if (uTag=="BLOCK_JACOBI_CORRECTED")
  {
    isCorrected_ = static_cast<bool> (param.getImmutableValue<int>());
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::initGraph
// Purpose       : Initialize the graph using information from the System
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 12/11/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiPrecond::initGraph( const Teuchos::RCP<Problem> & problem )
{
  // Generate the graph of each real equivalent form and then generate
  // empty linear systems for each frequency.
  RCP<Matrix> appdQdx = rcp( builder_.createMatrix() );
  RCP<Matrix> appdFdx = rcp( builder_.createMatrix() );

  // Relate the conductance and inductance matrices:
  // G(t) = df/dx(x(t))
  // C(t) = dq/dx(x(t))
  // Compute:  (omega*i*j*C_bar + G_bar)^{-1} for i=0,1,...,M,-M,...,-1
  //           using real equivalent form K_1 = [G_bar -i*omega*C_bar; i*omega*C_bar G_bar]

  // Generate the new real equivalent map and graph
  RCP<Parallel::ParMap> origMap = builder_.getSolutionMap(); 
  const Xyce::Parallel::Communicator& comm = origMap->pdsComm();
  int origLocalRows = origMap->numLocalEntities();
  int origGlobalRows = origMap->numGlobalEntities();

  // Count up the row indices and number of nonzero entries for the 2x2 block matrix.
  int maxRefNNZs_ = 0;
  int refRows = 2*origLocalRows;
  std::vector<int> rowIdxs( refRows ), refNNZs(refRows);
  RCP<Parallel::EpetraParMap> e_origMap = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>(origMap);
  int * origIdxs = e_origMap->petraMap()->MyGlobalElements();
  for (int i=0; i<origLocalRows; ++i)
  {
    rowIdxs[i] = origIdxs[i];
    rowIdxs[origLocalRows+i] = origIdxs[i] + origGlobalRows;
    refNNZs[i] = appdQdx->getLocalRowLength(i) + appdFdx->getLocalRowLength(i);
    refNNZs[origLocalRows+i] = refNNZs[i];
    if (refNNZs[i] > maxRefNNZs_) maxRefNNZs_ = refNNZs[i];
  }

  // Create the real equivalent map
  epetraMap_ = rcp(new Epetra_Map( -1, refRows, &rowIdxs[0], 0, e_origMap->petraMap()->Comm() ) );

  // Epetra graph is being created with static constructor
  epetraGraph_ = rcp(new Epetra_CrsGraph( Copy, *epetraMap_, &refNNZs[0], true ));

  // Communicate the maximum number of nonzeros for all processors.
  int tmpMaxNNZs = maxRefNNZs_;
  comm.sumAll( &tmpMaxNNZs, &maxRefNNZs_, 1 );

  // Put together the indices for each row and insert them into the graph.
  int tmpNNZs=0, tmpNNZs2=0;
  std::vector<double> tmpCoeffs(maxRefNNZs_);
  std::vector<int> refIdxs(maxRefNNZs_), refIdxs2(maxRefNNZs_);

   for ( int i=0; i<origLocalRows; ++i ) {

    // Get the indices for the first block of the matrix (G_bar)
    appdFdx->getRowCopy( rowIdxs[i], maxRefNNZs_, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

    // Get the indices for the third block of the matrix (C_bar)
    appdQdx->getRowCopy( rowIdxs[i], maxRefNNZs_, tmpNNZs2, &tmpCoeffs[0], &refIdxs2[0] );

    // Insert the indices for the third block into refIdxs, as they are the indices of the second block
    for (int j=0; j<tmpNNZs2; ++j) {
      refIdxs[tmpNNZs+j] = refIdxs2[j]+origGlobalRows;
    }
    epetraGraph_->InsertGlobalIndices( rowIdxs[i], refNNZs[i], &refIdxs[0] );

    // Insert the indices for the first block into refIdxs2, as they are the indices of the fourth block
    for (int j=0; j<tmpNNZs; ++j) {
      refIdxs2[tmpNNZs2+j] = refIdxs[j]+origGlobalRows;
    }
    epetraGraph_->InsertGlobalIndices( rowIdxs[origLocalRows+i], refNNZs[origLocalRows+i], &refIdxs2[0] );
  }
  
  epetraGraph_->FillComplete();

  int numProcs = comm.numProc();
  int myPID = comm.procID();

  // Get the Fourier series information and generate the Epetra_LinearSystems.
  RCP<BlockVector> bXt = hbBuilderPtr_->createTimeDomainBlockVector();
  N_ = bXt->blockCount();
  M_ = (int)((N_-1)/2);

  if (numProcs > 1)
  {
    // Create the serial maps for all processors.
    // NOTE:  Yes, there must be as many maps/importers/graphs as processors.
    serialEpetraMap_.resize( numProcs );
    serialImporter_.resize( numProcs );
    serialGraph_.resize( numProcs );

    // Determine how many linear systems this processor will manage.
    int nOverP = (M_+1)/numProcs;
    if (nOverP > 0)
    {
      beginN_ = myPID*nOverP;
      endN_ = (myPID+1)*nOverP;
      if ( (M_+1)%numProcs && myPID==(numProcs-1))
        endN_ += (M_+1)%numProcs;
    }
    else
    {
      if (myPID < M_+1)
      {
        beginN_ = myPID;
        endN_ = myPID+1;
      }
      else
        beginN_ = endN_ = -1;
    }

    for (int proc = 0; proc < numProcs; ++proc )
    {
      serialEpetraMap_[proc] = Teuchos::rcp( new Epetra_Map( Epetra_Util::Create_Root_Map( *epetraMap_, proc ) ) );
      serialImporter_[proc] = Teuchos::rcp( new Epetra_Import( *serialEpetraMap_[proc], *epetraMap_ ) ); 
      serialGraph_[proc] = Teuchos::rcp( new Epetra_CrsGraph( Copy, *serialEpetraMap_[proc], 0 ) );
      serialGraph_[proc]->Import( *epetraGraph_, *serialImporter_[proc], Insert );
      serialGraph_[proc]->FillComplete(); 
      serialGraph_[proc]->OptimizeStorage(); 
    }

    // Get the GID list from the serial map, so the split map has the same order of GIDs in parallel.
    int numAll = serialEpetraMap_[myPID]->NumGlobalElements();
    gidList_.resize( numAll );
    for (int i=0; i<numAll; i++)
    {
      gidList_[i] = serialEpetraMap_[myPID]->GID(i);
    }

    // Now split the communicator so that each processor solves their own linear systems.
#ifdef Xyce_PARALLEL_MPI
    MPI_Comm splitComm;
    Epetra_MpiComm mpiComm( dynamic_cast<Epetra_MpiComm &>( const_cast<Epetra_Comm &>( epetraMap_->Comm() ) ) );
    MPI_Comm_split( mpiComm.GetMpiComm(), myPID, myPID, &splitComm); 

    int n=0, maxLen=serialGraph_[myPID]->MaxNumIndices();
    std::vector<int> lIdxs( maxLen ), gNnzs( numAll );
    std::vector< std::vector<int> > gIdxs( numAll );
    Epetra_MpiComm singleComm( splitComm );
    singleMap_ = Teuchos::rcp( new Epetra_Map( numAll, 0, singleComm ) );
    for (int i=0; i<numAll; i++)
    {
      serialGraph_[myPID]->ExtractMyRowCopy( i, maxLen, n, &lIdxs[0] );
      for (int j=0; j<n; j++)
      {
        gIdxs[i].push_back( gidList_[ lIdxs[j] ] );
      } 
      gNnzs[i] = n;
    }

    // Generate Epetra_CrsGraph using a static profile given by gNnzs.
    singleGraph_ = Teuchos::rcp( new Epetra_CrsGraph( Copy, *singleMap_, &gNnzs[0], true ) );

    // Insert values into graph that were collected in gIdxs
    for (int i=0; i<numAll; i++)
    {
      singleGraph_->InsertGlobalIndices( gidList_[i], gNnzs[i], &(gIdxs[i])[0] );
    }

    // Now the fill is complete.  No need to optimize storage, a static constructor was used.
    singleGraph_->FillComplete();
#endif
  } 
  else
  {
    beginN_ = 0;
    endN_ = M_+1;
  }
  
  Amesos amesosFactory;
  Teuchos::ParameterList params;
  // We are always moving the system to one processor, so the solver should not check inputs to reduce overhead.
  params.set( "TrustMe", true );
  
  // Generate the vectors and matrices for the N_ linear problems to be solved on the diagonal.
  // epetraMatrix_               : object  having the same distribution as the original system.
  // serialMatrix_/serialVector_ : objects having all rows on a single processor; target for exporting dist objects.
  // singleMatrix_/singleVector_ : objects having all rows on a single processor with a split communicator.
  if (numProcs > 1)
  {
    epetraMatrix_.resize(1);
    epetraMatrix_[0] = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *epetraGraph_ ) );
    epetraMatrix_[0]->FillComplete();

    serialVector_ = Teuchos::rcp( new Epetra_MultiVector( *serialEpetraMap_[myPID], 1 ) );

    // Update the importer for performing the matrix importing.
    for (int proc = 0; proc < numProcs; ++proc )
    {
      serialMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *serialGraph_[proc] ) );
      serialImporter_[proc] = Teuchos::rcp( new Epetra_Import( serialMatrix_->RowMap(), epetraMatrix_[0]->RowMap() ) ); 
    }
 
    singleMatrix_.resize(endN_-beginN_);
    epetraProblem_.resize(endN_-beginN_);
    amesosPtr_.resize(endN_-beginN_);

#ifdef Xyce_PARALLEL_MPI 
    singleRHS_ = Teuchos::rcp( new Epetra_MultiVector( *singleMap_, 1 ) );
    singleSoln_ = Teuchos::rcp( new Epetra_MultiVector( *singleMap_, 1 ) );

    for (int i=0; i<endN_-beginN_; ++i) {
      singleMatrix_[i] = rcp( new Epetra_CrsMatrix( Copy, *singleGraph_ ) );
      singleMatrix_[i]->FillComplete();
      epetraProblem_[i] = rcp( new Epetra_LinearProblem( &*singleMatrix_[i], &*singleRHS_, &*singleSoln_ ) );
      amesosPtr_[i] = rcp( amesosFactory.Create( "Klu", *epetraProblem_[i] ) );
      amesosPtr_[i]->SetParameters( params );
      amesosPtr_[i]->SymbolicFactorization();
    }
#endif
  }
  else 
  {
    epetraRHS_ = Teuchos::rcp( new Epetra_MultiVector( *epetraMap_, 1 ) );
    epetraSoln_ = Teuchos::rcp( new Epetra_MultiVector( *epetraMap_, 1 ) );
    epetraMatrix_.resize(M_+1);
    epetraProblem_.resize(M_+1);
    amesosPtr_.resize(M_+1);

    for (int i=0; i<M_+1; ++i) {
      epetraMatrix_[i] = rcp( new Epetra_CrsMatrix( Copy, *epetraGraph_ ) );
      epetraMatrix_[i]->FillComplete();
      epetraProblem_[i] = rcp( new Epetra_LinearProblem( &*epetraMatrix_[i], &*epetraRHS_, &*epetraSoln_ ) );
      amesosPtr_[i] = rcp( amesosFactory.Create( "Klu", *epetraProblem_[i] ) );
      amesosPtr_[i]->SetParameters( params );
      amesosPtr_[i]->SymbolicFactorization();
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::initValues
// Purpose       : Initialize the values of the Epetra_LinearSystem's
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 12/11/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiPrecond::initValues( const Teuchos::RCP<Problem> & problem )
{
  int size_ = freqs_.size();
  int posFreq = (size_-1)/2;

  RCP<Parallel::ParMap> origMap = builder_.getSolutionMap();
  const Xyce::Parallel::Communicator& comm = origMap->pdsComm();

  if (VERBOSE_LINEAR)
  {
    Xyce::dout() << "HBBlockJacobiPrecond::initValues: " << std::endl;
  }

  // Get matrix to contain the sum of all dQdx and dFdx
  // NOTE:  In parallel, this will contain only the nonlinear portion of dQdx and dFdx
  Teuchos::RCP<Matrix> appdQdxSum = rcp( builder_.createMatrix() );
  Teuchos::RCP<Matrix> appdFdxSum = rcp( builder_.createMatrix() );

  // Initialize memory.
  appdQdxSum->put(0.0);
  appdFdxSum->put(0.0);

  int numProcs = comm.numProc();
  int myPID = comm.procID();

  // Get the separated stored Jacobian matrices from the HB loader.
  Teuchos::RCP<FilteredMatrix>& linAppdQdx = hbLoaderPtr_->getStoreLindQdx();
  Teuchos::RCP<FilteredMatrix>& linAppdFdx = hbLoaderPtr_->getStoreLindFdx();

  std::vector<Teuchos::RCP<FilteredMatrix> >& vecNLAppdQdx = hbLoaderPtr_->getStoreNLdQdx();
  std::vector<Teuchos::RCP<FilteredMatrix> >& vecNLAppdFdx = hbLoaderPtr_->getStoreNLdFdx();

  // Get the frequency-domain Jacobian matrix entries from the HB loader.
  const std::vector< std::vector< Util::FreqMatEntry > >& linFreqdFdx =
    hbLoaderPtr_->getFreqDFDXMatrix();

  const std::vector<int>& linFreqNZRows = hbLoaderPtr_->getFreqNZLocalRows();
  std::map<int,int> linFreqNZRowMap = hbLoaderPtr_->getFreqNZLocalRowsMap();

  // Sum up dQdx and dFdx
  double periodLength = times_[N_] - times_[0];
  for( int i = 0; i < N_; ++i )
  {
    // Compute scale factor.
    double timeStep = times_[i+1] - times_[i];
    double scaleFactor = (timeStep/periodLength);      

    // Add into nonlinear dQdxSum and dFdxSum with scaling
    vecNLAppdQdx[i]->addToMatrix( *appdQdxSum, scaleFactor );
    vecNLAppdFdx[i]->addToMatrix( *appdFdxSum, scaleFactor );
  }

  // Compute Cdiff and Gdiff of nonlinear matrices, if necessary
  if (isCorrected_)
  {
    // Resize the storage for diffC and diffG.
    diffCMatrix_.resize(N_);
    diffGMatrix_.resize(N_);

    Teuchos::RCP<Matrix> tmpCMatrix = rcp( builder_.createMatrix() );
    Teuchos::RCP<Matrix> tmpGMatrix = rcp( builder_.createMatrix() );

    for( int i = 0; i < N_; ++i )
    {
      // Compute C_avg - C_i and G_avg - G_i
      // NOTE:  This is the negative of what we need, but this will be fixed in the operator.
      tmpCMatrix->put( 0.0 );
      tmpCMatrix->add( *appdQdxSum );
      vecNLAppdQdx[i]->addToMatrix( *tmpCMatrix, -1.0 );
      diffCMatrix_[i] = Teuchos::rcp( new Xyce::Linear::FilteredMatrix( tmpCMatrix.get(), &*origMap, false ) );

      tmpGMatrix->put( 0.0 );
      tmpGMatrix->add( *appdFdxSum );
      vecNLAppdFdx[i]->addToMatrix( *tmpGMatrix, -1.0 );
      diffGMatrix_[i] = Teuchos::rcp( new Xyce::Linear::FilteredMatrix( tmpGMatrix.get(), &*origMap, false ) );
    } 
  }

  // Add in linear device contributions, without scaling.
  linAppdQdx->addToMatrix( *appdQdxSum );

  // Add in linear device contributions, without scaling.
  linAppdFdx->addToMatrix( *appdFdxSum );

/*
  std::cout << "appdQdxSum = " << std::endl;
  appdQdxSum->print( std::cout );
  std::cout << "appdFdxSum = " << std::endl;
  appdFdxSum->print( std::cout );
*/
  int tmpNNZs=0;
  double cplxCoeff = 0.0;
  std::vector<double> tmpCoeffs(maxRefNNZs_), tmpCoeffs2(2*maxRefNNZs_);
  std::vector<int> refIdxs(maxRefNNZs_), refIdxs2(2*maxRefNNZs_);

  // Get the values for each row of appdFdxSum/appdQdxSum and insert them into the matrix.

  // Since we are reusing this object, reinitialize it. 
  if (numProcs > 1)
      epetraMatrix_[0]->PutScalar(0.0);

  int origLocalRows = origMap->numLocalEntities();
  int origGlobalRows = origMap->numGlobalEntities();

  // Only solving for the positive frequencies and copying over conjugate values.
  for ( int nB=0; nB<=M_; ++nB ) 
  {
    // Compute the coefficient on the complex matrix [0,1,...,M,-M,...,-1]
    cplxCoeff = 2.0 * M_PI * freqs_[posFreq+nB];

    // Insert the entries for all four blocks of the real-equivalent form
    for ( int i=0; i<origLocalRows; ++i ) 
    {
      // Get the global ID for this row.
      int gid = epetraMap_->GID( i );

      // Load [G_bar -i*omega*C_bar; i*omega*C_bar G_bar]

      // Get the indices for the first block of the matrix (G_bar)
      appdFdxSum->getRowCopy( gid, maxRefNNZs_, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      
      if (numProcs > 1)
        epetraMatrix_[0]->ReplaceGlobalValues( gid, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      else
        epetraMatrix_[nB]->ReplaceGlobalValues( gid, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

      // Modify the indices for the fourth block
      for ( int j=0; j<tmpNNZs; ++j ) 
      {
        refIdxs[j] += origGlobalRows;
      }
      if (numProcs > 1)
        epetraMatrix_[0]->ReplaceGlobalValues( gid+origGlobalRows, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      else
        epetraMatrix_[nB]->ReplaceGlobalValues( gid+origGlobalRows, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

      // Add in complex part of REF matrix if non-zero.
      if (cplxCoeff != 0.0) 
      {
        // Get the indices for the first block of the matrix (G_bar)
        appdQdxSum->getRowCopy( gid, maxRefNNZs_, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

        // Insert the second block
        for (int ii=0; ii<tmpNNZs; ++ii) { tmpCoeffs[ii] *= cplxCoeff; }
        if (numProcs > 1)
          epetraMatrix_[0]->ReplaceGlobalValues( gid+origGlobalRows, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
        else
          epetraMatrix_[nB]->ReplaceGlobalValues( gid+origGlobalRows, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );

        // Insert the third block
        for (int ii=0; ii<tmpNNZs; ++ii) {
          tmpCoeffs[ii] *= -1.0;
          refIdxs[ii] += origGlobalRows;
        }
        if (numProcs > 1)
          epetraMatrix_[0]->ReplaceGlobalValues( gid, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
        else 
          epetraMatrix_[nB]->ReplaceGlobalValues( gid, tmpNNZs, &tmpCoeffs[0], &refIdxs[0] );
      }
    }

    // Insert frequency domain contributions for this frequency.
    // First collect individual entries and convert local ids to global ids.
    int numRows = linFreqNZRows.size();
    std::vector< std::vector< double > > realEntries(numRows), imagEntries(numRows);
    std::vector< std::vector< int > > indices(numRows);
    std::vector< Util::FreqMatEntry >::const_iterator it = linFreqdFdx[nB].begin(); 
    for ( ; it != linFreqdFdx[nB].end(); it++ ) 
    {
       // Get the global ID for this row.
       int idx = linFreqNZRowMap[it->row_lid];
       int gid = appdFdxSum->getGraph()->localToGlobalColIndex( it->col_lid );
       if (gid != -1)
       {
         realEntries[idx].push_back( (*it).val.real() );
         imagEntries[idx].push_back( (*it).val.imag() );
         indices[idx].push_back( gid );
       }
    }

    // Second, sum into global values.
    for ( int i=0; i<numRows; i++ )
    {
      int gid = appdFdxSum->getGraph()->localToGlobalRowIndex( linFreqNZRows[i] );
      
      // Insert first quadrant, G
      if (numProcs > 1)
        epetraMatrix_[0]->SumIntoGlobalValues( gid, realEntries[i].size(), &realEntries[i][0], &indices[i][0] );
      else
        epetraMatrix_[nB]->SumIntoGlobalValues( gid, realEntries[i].size(), &realEntries[i][0], &indices[i][0] );

      // Insert second quadrant, positive C.
      if (cplxCoeff != 0.0)
      {
        if (numProcs > 1)
          epetraMatrix_[0]->SumIntoGlobalValues( gid+origGlobalRows, imagEntries[i].size(), &imagEntries[i][0], &indices[i][0] );
        else
          epetraMatrix_[nB]->SumIntoGlobalValues( gid+origGlobalRows, imagEntries[i].size(), &imagEntries[i][0], &indices[i][0] );
      }

      // Shift column indices for third and fourth quadrant loads and conjugate imagEntries[i].
      for ( unsigned j=0; j<indices[i].size(); ++j )
      {
        indices[i][j] += origGlobalRows;
        imagEntries[i][j] *= -1.0;
      }

      // Insert third quadrant, conjugate C.
      if (cplxCoeff != 0.0)
      {
        if (numProcs > 1)
          epetraMatrix_[0]->SumIntoGlobalValues( gid, imagEntries[i].size(), &imagEntries[i][0], &indices[i][0] );
        else
          epetraMatrix_[nB]->SumIntoGlobalValues( gid, imagEntries[i].size(), &imagEntries[i][0], &indices[i][0] );
      }

      // Insert fourth quadrant, G
      if (numProcs > 1)
        epetraMatrix_[0]->SumIntoGlobalValues( gid+origGlobalRows, realEntries[i].size(), &realEntries[i][0], &indices[i][0] );
      else
        epetraMatrix_[nB]->SumIntoGlobalValues( gid+origGlobalRows, realEntries[i].size(), &realEntries[i][0], &indices[i][0] );
    }

    // Wait for the other processors.
    epetraMatrix_[0]->Comm().Barrier();

    // In parallel, move the current matrix information to a local processor 
    if (numProcs > 1)
    {
      // Sum in any off processor contributions from the frequency-domain loading.
      epetraMatrix_[0]->FillComplete();

      // Migrate the matrix to one processor.
      int nOverP = (M_+1) / numProcs;
      int proc = nB;
      if (nOverP > 0)
        proc = nB / nOverP;
      if (proc >= numProcs)
        proc = numProcs-1;

      serialMatrix_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *serialGraph_[proc] ) );
      serialMatrix_->Import( *epetraMatrix_[0], *serialImporter_[proc], Insert );
      serialMatrix_->FillComplete();

      if ( nB >= beginN_ && nB < endN_ )
      {
        int idx = nB - beginN_;
        int myRows = serialEpetraMap_[myPID]->NumGlobalElements();
        for (int i=0; i<myRows; i++)
        {
          int gid_i=gidList_[i];
          serialMatrix_->ExtractGlobalRowCopy( gid_i, 2*maxRefNNZs_, tmpNNZs, &tmpCoeffs2[0], &refIdxs2[0] );
          singleMatrix_[idx]->ReplaceGlobalValues( gid_i, tmpNNZs, &tmpCoeffs2[0], &refIdxs2[0] );
        }
      } 
      
      // Wait for the other processors.
      epetraMatrix_[0]->Comm().Barrier();
    }

  } // for ( int nB=0; nB<=M_; ++nB ) 

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::compute
// Purpose       : Compute a preconditioner M such that M ~= A^{-1}.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiPrecond::compute()
{
  bool precStatus = true;
  
  if (VERBOSE_LINEAR)
  {
    Xyce::dout() << "HBBlockJacobiPrecond::compute: " << std::endl;
  }

  // Compute numeric factorization for each block.
  for ( int i=0; i<endN_-beginN_; ++i ) {
    amesosPtr_[i]->NumericFactorization();
  }

  if ( Teuchos::is_null( epetraPrec_ ) )
    epetraPrec_ = blockJacobiOperator( epetraProblem_, amesosPtr_, 
                                       diffCMatrix_, diffGMatrix_, 
                                       hbLoaderPtr_, hbBuilderPtr_, 
                                       freqs_, std::pair<int,int>(beginN_, endN_), 
                                       hbOsc_ );

  if ( Teuchos::is_null( epetraPrec_ ) )
    return false;

  return precStatus;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiPrecond::apply
// Purpose       : Calls the actual preconditioner to apply y = M*x.
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
int HBBlockJacobiPrecond::apply( MultiVector & x, MultiVector & y )
{
  int precStatus = 0;

  EpetraVectorAccess* e_x = dynamic_cast<EpetraVectorAccess *>( &x );
  EpetraVectorAccess* e_y = dynamic_cast<EpetraVectorAccess *>( &y );

  if (VERBOSE_LINEAR)
  {
    Xyce::dout() << "HBBlockJacobiPrecond::apply: " << std::endl;
  }

  // If there is no preconditioner to apply return a nonzero code
  if( Teuchos::is_null(epetraPrec_) )
    precStatus = -1;
  else
    precStatus = epetraPrec_->Apply( e_x->epetraObj(), e_y->epetraObj() );

  return precStatus;
}

} // namespace Linear
} // namespace Xyce

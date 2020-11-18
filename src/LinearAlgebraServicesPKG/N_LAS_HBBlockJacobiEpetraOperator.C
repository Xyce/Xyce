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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Creator        : Heidi Thornquist, 1437
//
// Creation Date  : 09/04/08
//
//
//
//
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_HBBlockJacobiEpetraOperator.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_EpetraParMap.h>
#include <N_PDS_Comm.h>
#include <N_UTL_Math.h>
#include <N_PDS_EpetraParMap.h>

// ----------   Trilinos Includes   ----------

#include <Epetra_LinearProblem.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Amesos_BaseSolver.h>

using Teuchos::RCP;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : blockJacobiOperator 
// Purpose       : non-member constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
RCP<HBBlockJacobiEpetraOperator> blockJacobiOperator(
    std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
    const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
    const std::vector<Teuchos::RCP<FilteredMatrix> >& diffCMatrix,
    const std::vector<Teuchos::RCP<FilteredMatrix> >& diffGMatrix,
    const Teuchos::RCP<Loader::HBLoader>& hbLoader,
    const Teuchos::RCP<HBBuilder>& hbBuilder,
    const std::vector<double>& freqs,
    const std::pair<int,int>& localRange,
    const bool hbOsc
    )
{
  RCP<HBBlockJacobiEpetraOperator> epetraOperator =
    rcp(new HBBlockJacobiEpetraOperator);
  epetraOperator->initialize(epetraProblems,
      amesosSolvers,
      diffCMatrix,
      diffGMatrix,
      hbLoader,
      hbBuilder,
      freqs,
      localRange,
      hbOsc
      );
  return epetraOperator;
}


//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::HBBlockJacobiEpetraOperator
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
HBBlockJacobiEpetraOperator::HBBlockJacobiEpetraOperator()
{
  isInitialized_ = false;
  isCorrected_ = false;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::HBBlockJacobiEpetraOperator
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
HBBlockJacobiEpetraOperator::~HBBlockJacobiEpetraOperator()
{
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::initialize
// Purpose       : Initialization
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
void HBBlockJacobiEpetraOperator::initialize(
      std::vector<Teuchos::RCP<Epetra_LinearProblem> >& epetraProblems,
      const std::vector<Teuchos::RCP<Amesos_BaseSolver> >& amesosSolvers,
      const std::vector<Teuchos::RCP<FilteredMatrix> >& diffCMatrix,
      const std::vector<Teuchos::RCP<FilteredMatrix> >& diffGMatrix,
      const Teuchos::RCP<Loader::HBLoader>& hbLoader,
      const Teuchos::RCP<HBBuilder>& hbBuilder,
      const std::vector<double>& freqs,
      const std::pair<int,int>& localRange,
      const bool hbOsc
    )
{
  epetraProblems_ = epetraProblems;
  amesosSolvers_ = amesosSolvers;
  diffCMatrix_ = diffCMatrix;
  diffGMatrix_ = diffGMatrix;
  hbLoader_ = hbLoader;
  hbBuilder_ = hbBuilder;
  freqs_ = freqs;
  myN_ = localRange;
  hbOsc_ = hbOsc;
  int globalUnk = hbBuilder_->getSolutionMap()->numGlobalEntities();
  int localUnk = epetraProblems[0]->GetMatrix()->NumGlobalRows();
  N_ = globalUnk / localUnk;
  numAugRows_ = (hbBuilder_->getAugmentedLIDs()).size();
#ifdef Xyce_PARALLEL_MPI
  int numProcs = hbBuilder_->getPDSComm()->numProc();
  if (numProcs > 1)
  {
    serialEpetraMap_.resize(numProcs);
    serialImporter_.resize(numProcs);
    Teuchos::RCP<N_PDS_EpetraParMap> e_map = Teuchos::rcp_dynamic_cast<N_PDS_EpetraParMap>(hbBuilder_->getSolutionMap()); 
    for (int proc = 0; proc < numProcs; ++proc )
    {
      serialEpetraMap_[proc] = Teuchos::rcp( new Epetra_Map( Epetra_Util::Create_Root_Map( *(e_map->petraMap()), proc ) ) );
    }
    
    // Get a sum of all the augmented rows.
    int tmpSize = numAugRows_;
    hbBuilder_->getPDSComm()->sumAll( &tmpSize, &numAugRows_, 1 );
  }
#endif
  M_ = (freqs.size()-1)/2;

  isInitialized_ = true;
  isCorrected_ = (diffCMatrix.size() > 0) ? true : false;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::SetUseTranspose
// Purpose       : Define if transpose Apply and ApplyInverse is to be used.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::SetUseTranspose(bool UseTranspose)
{
  // This is not supported for the HB load layers.
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::ApplyInverse
// Purpose       : Apply matrix free preconditioner with Epetra_MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::ApplyInverse(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  int status = 0;

  // Convert these to MultiVectors and call the other Apply
  // Cast away the const until the Apply which will enforce it.
  // This is necessary because there is no const view function in MultiVector
  std::vector<double> norm(X.NumVectors(), 0.0);
  X.NormInf( &norm[0] );
  double max = *std::max_element(norm.begin(), norm.end());
  if (max > 0.0)
  {
    Epetra_MultiVector* Xptr = const_cast<Epetra_MultiVector*>(&X);
    MultiVector las_X(Xptr, false);  
    MultiVector las_Y(&Y, false);

    status = ApplyInverse(las_X,las_Y);
  }
  else
  {
    Y.PutScalar( 0.0 );
  }

  return(status);
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::ApplyInverse
// Purpose       : Apply matrix free preconditioner with MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::ApplyInverse(
  const MultiVector& X,
  MultiVector& Y
  ) const
{
  if (!isInitialized_)
  {
    std::string msg = "HBBlockJacobiEpetraOperator::ApplyInverse:  I'm not initialized!";
    Report::DevelFatal0() << msg;
  }

  int size = (hbBuilder_->getAugmentedLIDs()).size();

  std::vector< std::vector <double> > Ytmp (X.numVectors(), std::vector < double>(size) );

  if (hbOsc_ && size)
  {
    const std::vector<int>& augLIDs = hbBuilder_->getAugmentedLIDs();
    for (int i=0; i<X.numVectors() ; ++i) 
    {

      int j = 0;
      for (std::vector<int>::const_iterator it = augLIDs.begin(); it != augLIDs.end(); it++)
      {
        Ytmp[i][j] = X[i][*it];
        j++;
      } 
    }
  }

  // Apply original block Jacobi preconditioner
  // Y = J_BD^{-1}*X 
  ApplyBlockJacobi( X, Y );

  // Apply correction to the block Jacobi preconditioner, if necessary.
  if (isCorrected_ && (diffCMatrix_.size()>0))
  {
     // Create temporary vectors for application of block Jacobi and correction.
     Teuchos::RCP<BlockVector> tmpX = hbBuilder_->createExpandedRealFormTransposeBlockVector();
     Teuchos::RCP<BlockVector> tmpY = hbBuilder_->createExpandedRealFormTransposeBlockVector();

     // tmpX = (omega*D*C_diff*D^{-1} + D*G_diff*D^{-1})*tmpY
     ApplyCorrection(Y, *tmpX);

     // X = -J_BD^{-1}*tmpX
     ApplyBlockJacobi(*tmpX, *tmpY);
     Y.update( -1.0, *tmpY, 1.0 );
  }

  if (hbOsc_ && size)
  {
    const std::vector<int>& augLIDs = hbBuilder_->getAugmentedLIDs();

    for (int i=0; i<X.numVectors() ; ++i) 
    {
      
      int j = 0;
      for (std::vector<int>::const_iterator it = augLIDs.begin(); it != augLIDs.end(); it++)
      {
        Y[i][*it] = Ytmp[i][j];
        j++; 
      } 
    }
  }

  // Copy augmented rows from X to Y, identity preconditioner.
/*  if (hbOsc_ && (hbBuilder_->getAugmentedLIDs()).size())
  {
    const std::vector<int>& augLIDs = hbBuilder_->getAugmentedLIDs();
    for (int i=0; i<X.numVectors() ; ++i) 
    {
      for (std::vector<int>::const_iterator it = augLIDs.begin(); it != augLIDs.end(); it++)
      {
        Y[i][*it] = X[i][*it];
      } 
    }
  } */

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::ApplyBlockJacobi
// Purpose       : Apply matrix free preconditioner with MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::ApplyBlockJacobi(
  const MultiVector& X,
  MultiVector& Y
  ) const
{
  // Determine number of time-domain variables.
  int n = (X.globalLength()-numAugRows_) / (2*N_);
  int numProcs = hbBuilder_->getPDSComm()->numProc();
  int myPID = hbBuilder_->getPDSComm()->procID();

  Epetra_MultiVector *nB_RHS=0;
  double *nB_Soln=0;

  // If this is being run on multiple processors, load into the RHS vectors owned using a local vector.
  if (numProcs > 1)
  {
    // Create the importer if we don't have one
    if ( Teuchos::is_null( serialImporter_[myPID] ) )
    {
      serialX_ = Teuchos::rcp( new Epetra_MultiVector( *serialEpetraMap_[myPID], X.numVectors() ) );
      serialY_ = Teuchos::rcp( new Epetra_MultiVector( *serialEpetraMap_[myPID], Y.numVectors() ) );
      serialImporter_[myPID] = Teuchos::rcp( new Epetra_Import( *(serialEpetraMap_[myPID]), X.epetraObj().Map() ) );
    }

    // Copy all the RHS vectors to X.
    serialX_->Import( X.epetraObj(), *serialImporter_[myPID], Insert );
  }

  Teuchos::RCP<const Vector> x;
  Teuchos::RCP<Vector> y;
 
  int size = freqs_.size();

  for (int nB=myN_.first; nB<myN_.second; ++nB) 
  {
    nB_RHS = epetraProblems_[nB-myN_.first]->GetRHS();
    nB_Soln = epetraProblems_[nB-myN_.first]->GetLHS()->Values();

    for (int i=0 ; i<X.numVectors() ; ++i) 
    {
      if (numProcs > 1)
      {
        x = Teuchos::rcp( new Vector((*serialX_)(i), false) );
        y = Teuchos::rcp( new Vector((*serialY_)(i), false) );
      }
      else
      {
        x = Teuchos::rcp( X.getVectorViewAssembled(i), true );
        y = Teuchos::rcp( Y.getNonConstVectorViewAssembled(i), true );
      }

      for (int j=0; j<n; ++j) 
      {
        nB_RHS->ReplaceMyValue(j, 0, (*x)[j*(2*N_)+2*nB]);        // real
        nB_RHS->ReplaceMyValue(n+j, 0, (*x)[j*(2*N_)+2*nB+1]);    // imaginary
      }

      amesosSolvers_[nB-myN_.first]->Solve();
   
      int nB2 = size-nB;
      // Copy the solutions back into y. 
      for (int j=0; j<n; ++j) 
      {
        // Copy positive frequency solution
        (*y)[j*(2*N_)+2*nB]   = nB_Soln[j];
        (*y)[j*(2*N_)+2*nB+1] = nB_Soln[n+j];
       
        if (nB)
        {
          // Copy negative frequency solution for any nonzero frequency (conjugate)
          (*y)[j*(2*N_)+2*nB2]   = nB_Soln[j];
          (*y)[j*(2*N_)+2*nB2+1] = -nB_Soln[n+j];
        }
      }
    }
  }

  // Return all the solutions back to Y, if running in parallel.
  if (numProcs > 1)
  {
    // Make sure the resulting vector is zeros, since this is the target.
    Y.putScalar( 0.0 );

    Y.epetraObj().Export( *serialY_, *serialImporter_[myPID], Add );

    // Wait for everyone else before moving on.
    Y.epetraObj().Comm().Barrier();
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::ApplyCorrection
// Purpose       : Apply matrix free preconditioner with MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::ApplyCorrection(
  const MultiVector& X,
  MultiVector& Y
  ) const
{
  // Create frequency and time domain work vectors.
  Teuchos::RCP<BlockVector> bXtPtr = hbBuilder_->createTimeDomainBlockVector();
  Teuchos::RCP<BlockVector> bCXtPtr = hbBuilder_->createTimeDomainBlockVector();
  Teuchos::RCP<BlockVector> bGXtPtr = hbBuilder_->createTimeDomainBlockVector();
  Teuchos::RCP<BlockVector> btmpYf = hbBuilder_->createExpandedRealFormTransposeBlockVector();
  Teuchos::RCP<BlockVector> bYf = hbBuilder_->createExpandedRealFormTransposeBlockVector();
  Teuchos::RCP<BlockVector> bYVec = hbBuilder_->createExpandedRealFormTransposeBlockVector();

  for (int col=0; col<X.numVectors() ; ++col)
  {
    // Loop over all the blocks and apply C_diff and G_diff to Xt.
    int blockCount = bXtPtr->blockCount();

    // Apply one column at a time to the multivector.
    Teuchos::RCP<const Vector> X_col = Teuchos::rcp( X.getVectorViewAssembled( col ) );
    BlockVector bXf( &*X_col, 2*N_ );

    // Permute the input vector from the frequency to time domain, since this
    // is a time domain preconditioner.
    hbLoader_->permutedIFT(bXf, &*bXtPtr);

    for( int i = 0; i < blockCount; ++i )
    {
      diffCMatrix_[ i ]->matvec(bXtPtr->block(i), bCXtPtr->block(i));
      diffGMatrix_[ i ]->matvec(bXtPtr->block(i), bGXtPtr->block(i));
    }

    // Scale the resulting vector from matvec since diffC/diffG are C_avg - Ci
    // and G_avg - Gi, the negative of what is needed.
    // NOTE:  This is due to the filtered matrices in HB loader.
    bCXtPtr->scale( -1.0 );
    bGXtPtr->scale( -1.0 );

    // Permute back from time to frequency domain.
    hbLoader_->permutedFFT(*bCXtPtr, btmpYf.get());
    hbLoader_->permutedFFT(*bGXtPtr, bYf.get());

    // Now scale btmpYf.
    int size = freqs_.size();
    int posFreq = (size-1)/2;
    double omega = 2.0 * M_PI * freqs_[posFreq];
 
    // Put block count in the frequency domain.
    blockCount = bYf->blockCount();
    int blockSize = bXf.block(0).globalLength();

    for( int i = 0; i < blockCount; ++i )
    {
      // QVec needs to be created here since only one processor owns each block
      // and we do not know which one it is in parallel.
      Vector& YVec = bYVec->block(i);
      Vector& freqVec = btmpYf->block(i);

      // Only one processor owns each block of the frequency-domain vector
      if (freqVec.localLength() > 0)
      {
        omega = 2.0 * M_PI * freqs_[posFreq];
    
        YVec[0] = -freqVec[1]*omega;
        YVec[1] = freqVec[0]*omega;
    
        for (int j=1; j < (blockSize/2+1)/2; ++j)
        {
          omega = 2.0 * M_PI * freqs_[posFreq+j];
          YVec[2*j] = -freqVec[2*j+1]*omega;
          YVec[2*(blockSize/2-j)] = -freqVec[2*j+1]*omega;
    
          YVec[2*j+1] = freqVec[2*j]*omega;
          YVec[2*(blockSize/2-j)+1] = -freqVec[2*j]*omega;
        }
      }
   
      // Correction := (omega*D^{-1}*C_diff*D + D^{-1}*G_diff*D) 
      bYf->block(i).update(1.0, YVec , 1.0);
    }

    // Assign the correction back to col of Y. 
    Teuchos::RCP<Vector> Y_col = Teuchos::rcp( Y.getNonConstVectorViewAssembled( col ) );
    Y_col->update( 1.0, *bYf, 0.0 );
  }

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::Apply
// Purpose       : Apply inverse of matrix free preconditioner with Epetra_MultiVectors
// Special Notes : Not supported!
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::Apply(
  const Epetra_MultiVector& X,
  Epetra_MultiVector& Y
  ) const
{
  std::string msg = "HBBlockJacobiEpetraOperator::Apply is not supported!";
  Report::DevelFatal0() << msg;
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::Apply
// Purpose       : Apply inverse of matrix free preconditioner with MultiVectors
// Special Notes : Not supported!
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
int HBBlockJacobiEpetraOperator::Apply(
  const MultiVector& X,
  MultiVector& Y
  ) const
{
  std::string msg = "HBBlockJacobiEpetraOperator::Apply is not supported!";
  Report::DevelFatal0() << msg;
  return -1;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::NormInf
// Purpose       : Norm Inf of matrix
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
double HBBlockJacobiEpetraOperator::NormInf() const
{
  std::string msg = "HBBlockJacobiEpetraOperator::NormInf is not supported!";
  Report::DevelFatal0() << msg;
  return -1.0;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::Label
// Purpose       : Label for operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const char * HBBlockJacobiEpetraOperator::Label() const
{
  return "Matrix Free Harmonic Balance Block Jacobi Preconditioner";
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::UseTranspose
// Purpose       : Query for useTranspose setting
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiEpetraOperator::UseTranspose() const
{
  // Use Transpose is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::HasNormInf
// Purpose       : Query for normInf support
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
bool HBBlockJacobiEpetraOperator::HasNormInf() const
{
  // Norm Inf is not supported, so always return false.
  return false;
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::Comm
// Purpose       : Return Epetra_Comm object
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Comm & HBBlockJacobiEpetraOperator::Comm() const
{
  if (!isInitialized_)
  {
    std::string msg = "HBBlockJacobiEpetraOperator::Comm:  I'm not initialized!";
    Report::DevelFatal0() << msg;
  }
  return(Teuchos::rcp_dynamic_cast<N_PDS_EpetraParMap>(hbBuilder_->getSolutionMap())->petraMap()->Comm());
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::OperatorDomainMap
// Purpose       : Return Epetra_Map corresponding to domain of operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Map & HBBlockJacobiEpetraOperator::OperatorDomainMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "HBBlockJacobiEpetraOperator::OperatorDomainMap:  I'm not initialized!";
    Report::DevelFatal0() << msg;
  }
  return(*Teuchos::rcp_dynamic_cast<N_PDS_EpetraParMap>(hbBuilder_->getSolutionMap())->petraMap());
}

//-----------------------------------------------------------------------------
// Function      : HBBlockJacobiEpetraOperator::OperatorRangeMap
// Purpose       : Return Epetra_Map corresponding to range of operator
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, 1437
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
const Epetra_Map & HBBlockJacobiEpetraOperator::OperatorRangeMap() const
{
  if (!isInitialized_)
  {
    std::string msg = "HBBlockJacobiEpetraOperator::OperatorRangeMap:  I'm not initialized!";
    Report::DevelFatal0() << msg;
  }
  return(*Teuchos::rcp_dynamic_cast<N_PDS_EpetraParMap>(hbBuilder_->getSolutionMap())->petraMap());
}

} // namespace Linear
} // namespace Xyce

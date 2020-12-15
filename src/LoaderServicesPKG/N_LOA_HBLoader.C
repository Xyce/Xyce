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
// Purpose       :
// Special Notes :
// Creator       :
// Creation Date :
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <numeric>

// ----------   Xyce Includes   ----------

#include <N_LOA_HBLoader.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Builder.h>
#include <N_LAS_HBBuilder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_MultiVector.h>

#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_UTL_Functors.h>

#include <N_PDS_ParMap.h>
#include <N_DEV_DeviceMgr.h>

#include <N_PDS_Comm.h>  

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace Xyce {
namespace Loader {

HBLoader::HBLoader(
  Device::DeviceMgr &                 device_manager,
  Linear::Builder &                   builder,
  const int refID,
  const bool hbOsc)
  : CktLoader(device_manager, builder),
    periodicTimesOffset_(0),
    period_(0),
    refID_(refID),
    hbOsc_(hbOsc),
    matrixFreeFlag_(false),
    deviceManager_(device_manager),
    freqLoadAnalysisDone_(false),
    numGlobalFreqRows_(0),
    totalNZOffProcRows_(0),
    totalOffProcBVecLIDs_(0),
    builder_(builder)
{
  // Now initialize all the time domain working vectors.
  appVecPtr_ = rcp(builder_.createVector());
  appNextStaVecPtr_ = rcp(builder_.createStateVector());
  appCurrStaVecPtr_ = rcp(builder_.createStateVector());
  appLastStaVecPtr_ = rcp(builder_.createStateVector());

  appdQdxPtr_ = rcp(builder_.createMatrix());
  appdFdxPtr_ = rcp(builder_.createMatrix());

  appNextStoVecPtr_ = rcp(builder_.createStoreVector());
  appCurrStoVecPtr_ = rcp(builder_.createStoreVector());
  
  appNextLeadFVecPtr_ = rcp(builder.createLeadCurrentVector());
  appLeadQVecPtr_     = rcp(builder.createLeadCurrentVector());
  appNextJunctionVVecPtr_ = rcp(builder.createLeadCurrentVector());
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::registerHBBuilder
// Purpose       : Registration method for the HB builder
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, Sandia Labs
// Creation Date : 11/08/13
//-----------------------------------------------------------------------------
void HBLoader::registerHBBuilder( Teuchos::RCP<Linear::HBBuilder> hbBuilderPtr )
{
  hbBuilderPtr_ = hbBuilderPtr;

  if (hbOsc_)
    bQPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

  // Now initialize all the frequency domain working vectors.
  bXtPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();
  bVtPtr_ = hbBuilderPtr_->createTimeDomainBlockVector();

  // Vectors related to lead currents
//  bStoreVecFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeStoreBlockVector();

  bLeadCurrentVecFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeLeadCurrentBlockVector();
  bLeadCurrentQVecFreqPtr_ = hbBuilderPtr_->createExpandedRealFormTransposeLeadCurrentBlockVector();
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::setFreq
// Purpose       : Assign times for fast time scale
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void HBLoader::setHBFreqs( const std::vector<double> & freqs)
{
   freqs_ = freqs;
}



//-----------------------------------------------------------------------------
// Function      : HBLoader::setFastTimes
// Purpose       : Assign times for fast time scale
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void HBLoader::setFastTimes( const std::vector<double> & times )
{
  times_ = times;
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool HBLoader::loadDAEMatrices( Linear::Vector * X,
                                Linear::Vector * S,
                                Linear::Vector * dSdt,
                                Linear::Vector * Store,
                                Linear::Matrix * dQdx,
                                Linear::Matrix * dFdx,
                                int loadType)
{
  if (matrixFreeFlag_)
  {
    if (DEBUG_HB)
    {
      Xyce::dout() << std::endl
                   << Xyce::section_divider << std::endl
                   << "  HBLoader::loadDAEMatrices:  matrixFree case" << std::endl;
    }

    dQdx->put(0.0);
    dFdx->put(0.0);
  }
  else
  {
    Report::DevelFatal0().in("HBLoader::loadDAEMatrices") << "This function actually was called in a non matrix free case";
  }

  return(true);
  
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::applyDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool HBLoader::applyDAEMatrices( Linear::Vector * Xf,
                                 Linear::Vector * S,
                                 Linear::Vector * dSdt,
                                 Linear::Vector * Store,
                                 const Linear::Vector & Vf,
                                 Linear::Vector * dQdxV,
                                 Linear::Vector * dFdxV )
{
  if (!matrixFreeFlag_)
  {
    Report::DevelFatal0().in("HBLoader::applyDAEMatrices") << "This function should only be called in the matrix free case";
  }
  if (DEBUG_HB) 
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  HBLoader::applyDAEMatrices" << std::endl;
  }

  const std::vector<int> & augmentedLIDs = hbBuilderPtr_->getAugmentedLIDs();
  double fScalar = 1.0;

  double vfScalar = 0.0;

  Linear::BlockVector & bXf = *dynamic_cast<Linear::BlockVector*>(Xf);

  // We have to do something special with Vf because AztecOO (or Belos)
  // probably used the LinearProblem's maps to create the input
  // vector here.  In this case, Vf is just an Linear::Vector and not a
  // Linear::BlockVector.
  const Linear::BlockVector bVf(&Vf, bXf.blockSize());

  if (hbOsc_)
  {
    double tmpfScalar = 0.0;   

    double tmpVf = 0.0; 

    if (augmentedLIDs.size() )
    {
      tmpfScalar = (*Xf)[(augmentedLIDs)[0 ]];
      tmpVf = (Vf)[(augmentedLIDs)[0 ]];
    }

    Xf->pmap()->pdsComm().sumAll( &tmpfScalar, &fScalar, 1 );
    Xf->pmap()->pdsComm().sumAll( &tmpVf, &vfScalar, 1 );

  }

  Linear::BlockVector * bdQdxV = dynamic_cast<Linear::BlockVector*>(dQdxV);
  Linear::BlockVector * bdFdxV = dynamic_cast<Linear::BlockVector*>(dFdxV);

  std::vector<double> norm(1, 0.0);
  Vf.infNorm( &norm[0] );
  if (norm[0] > 0.0)
  {
    // Calculate contributions of the portion of the Jacobian relating to the nonlinear devices.
    permutedIFT(bVf, &*bVtPtr_);

    // Initialize resulting vector, since operations on separated Jacobian may not touch every entry.
    bdFdxV->putScalar( 0.0 );
    bdQdxV->putScalar( 0.0 );

    Teuchos::RCP<Linear::BlockVector>  bdQdxVt = hbBuilderPtr_->createTimeDomainBlockVector();
    Teuchos::RCP<Linear::BlockVector>  bdFdxVt = hbBuilderPtr_->createTimeDomainBlockVector();

    int BlockCount = bVtPtr_->blockCount();
    std::vector< std::vector< Util::FreqVecEntry > > freqFVector( BlockCount );

    for( int i = 0; i < BlockCount; ++i )
    {
      if (DEBUG_HB)
      {
        Xyce::dout() << "Processing diagonal matrix block " << i << " of " << BlockCount-1 << std::endl;
      }

      // Get the already stored nonlinear time-domain Jacobian matrices
      vecNLAppdQdxPtr_[ i ]->matvec( bVtPtr_->block(i), bdQdxVt->block(i) );
      vecNLAppdFdxPtr_[ i ]->matvec( bVtPtr_->block(i), bdFdxVt->block(i) );

      if (DEBUG_HB)
      {
        Xyce::dout() << "bVtPtr block i = " << i << " : " << std::endl;
        bVtPtr_->block(i).print(dout());

        Xyce::dout() << "bdQdxVt block i = " << i << " : " << std::endl;
        bdQdxVt->block(i).print(dout());

        Xyce::dout() << "bdFdxVt block i = " << i << " : " << std::endl;
        bdFdxVt->block(i).print(dout());
      }
    }

    permutedFFT(*bdQdxVt, bdQdxV, &nonlinQNZRows_);
    permutedFFT(*bdFdxVt, bdFdxV, &nonlinFNZRows_);

    int blockCount = bXf.blockCount();
    int blockSize = bXf.block(0).globalLength();

    int size_ = freqs_.size();
    int posFreq = (size_-1)/2;
    double omega = 2.0 * M_PI * freqs_[posFreq] * fScalar;

    // Add in contributions from part of Jacobian related to the linear devices.
    // Act on bVf or Vf, see if we can put it into a MultiVector and 
    // apply linAppdQdxPtr_ and linAppdFdxPtr_ to the transpose.

    Teuchos::RCP<Linear::BlockVector> permlindQdxV = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();
    Teuchos::RCP<Linear::BlockVector> permlindFdxV = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

    applyLinearMatrices (Vf, *permlindQdxV, *permlindFdxV);

    if(!linAppdFdxPtr_->isEmpty() || numGlobalFreqRows_)
    {
      bdFdxV->update (1.0, *permlindFdxV, 1.0);
    } 
    if(!linAppdQdxPtr_->isEmpty()) 
    {
      bdQdxV->update (1.0, *permlindQdxV, 1.0);
    }

    Teuchos::RCP<Linear::BlockVector> bQVec = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();
    double omega_1 = 2.0 * M_PI * freqs_[posFreq];

    Teuchos::RCP<Linear::BlockVector> bcolVec;

    if (hbOsc_) 
      bcolVec = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

    for( int i = 0; i < blockCount; ++i )
    {
    // QVec needs to be created here since only one processor owns each block
    // and we do not know which one it is in parallel.
      Linear::Vector& QVec = bQVec->block(i);
      Linear::Vector& freqVec = bdQdxV->block(i);

    // Only one processor owns each block of the frequency-domain vector
      if (freqVec.localLength() > 0)
      {
        omega = 2.0 * M_PI * freqs_[posFreq] * fScalar;

        QVec[0] = -freqVec[1]*omega;
        QVec[1] = freqVec[0]*omega;
        

        for (int j=1; j < (blockSize/2+1)/2; ++j)
        {
          omega = 2.0 * M_PI * freqs_[posFreq+j] * fScalar;
          QVec[2*j] = -freqVec[2*j+1]*omega;
          QVec[2*(blockSize/2-j)] = -freqVec[2*j+1]*omega;

          QVec[2*j+1] = freqVec[2*j]*omega;
          QVec[2*(blockSize/2-j)+1] = -freqVec[2*j]*omega;
        }

        bdFdxV->block(i).update(1.0, QVec , 1.0);

      }

      if ( hbOsc_ && freqVec.localLength() > 0 )
      {

        Linear::Vector& fcolVec = bcolVec->block(i);
        Linear::Vector& fqVec = bQPtr_->block(i);

        omega_1 = 2.0 * M_PI * freqs_[posFreq];
        fcolVec[0] = -fqVec[1]*omega_1;
        fcolVec[1]  = fqVec[0]*omega_1;

        for (int j=1; j < (blockSize/2+1)/2; ++j)
        {
          omega_1 = 2.0 * M_PI * freqs_[posFreq+j];

          fcolVec[2*j] = -fqVec[2*j+1]*omega_1;
          fcolVec[2*(blockSize/2-j)] = -fqVec[2*j+1]*omega_1;

          fcolVec[2*j+1] = fqVec[2*j]*omega_1;
          fcolVec[2*(blockSize/2-j)+1] = -fqVec[2*j]*omega_1;
        }

        bdFdxV->block(i).update(vfScalar, fcolVec, 1.0);
      }

    }

  }
  else
  {
    // X is zero, so Y is zero.
    bdQdxV->putScalar( 0.0 );
    bdFdxV->putScalar( 0.0 );
  }

  if (hbOsc_)
  {

//    if (augmentedLIDs.size() )
//      (*bdFdxV)[(augmentedLIDs)[0 ]] = Vf[(augmentedLIDs)[0 ]];
    double refValue = 0.0;
    double tmpValue= 0.0;
    Linear::Vector & freqVec = bVf.block(refID_);

    if (freqVec.localLength() > 0)
    {
       tmpValue = freqVec[3];
    }

    bXf.pmap()->pdsComm().sumAll( &tmpValue, &refValue, 1 );

    if (augmentedLIDs.size() )
      (*bdFdxV)[(augmentedLIDs)[0 ]] = refValue;
  }

  if (DEBUG_HB)
  {
    Xyce::dout() << "HB bVf:" << std::endl;
    bVf.print(std::cout);
    Xyce::dout() << "HB bdQdxV:" << std::endl;
    bdQdxV->print(std::cout);
    Xyce::dout() << "HB bdFdxV:" << std::endl;
    bdFdxV->print(std::cout);

    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : HBLoader::updateState
// Purpose       :
// Special Notes : ERK.  This function needs to be a no-op.  The reason
//                 is that the state information needs to be the same
//                 at the time of updateState, loadDAEVectors and
//                 loadDAEMatrices.  Thus, they have to all happen inside
//                 of the same "fast time" loop.  So, this functionality
//                 has been moved up into loadDAEVectors.
//
//                 Note: for similar reasons, loadDAEMatrices is called from
//                 within that function as well.
//
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool HBLoader::updateState
 (Linear::Vector * nextSolVectorPtr,
  Linear::Vector * currSolVectorPtr,
  Linear::Vector * lastSolVectorPtr,
  Linear::Vector * nextStaVectorPtr,
  Linear::Vector * currStaVectorPtr,
  Linear::Vector * lastStaVectorPtr,
  Linear::Vector * nextStoVectorPtr,
  Linear::Vector * currStoVectorPtr,
  int loadType
  )
{
  bool bsuccess = true;

  // For HB case, this needs to be a no-op.

  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : HBLoader::applyLinearMatrices
// Purpose       : apply matrices related to linear nodes
// Special Notes : This is done in frequency domain
// Scope         : public
//
// Creator       : Ting Mei and Heidi Thornquist
// Creation Date : 04/10/2015
//-----------------------------------------------------------------------------
bool HBLoader::applyLinearMatrices( const Linear::Vector & Vf,
                                    Linear::BlockVector & permlindQdxV,
                                    Linear::BlockVector & permlindFdxV )
{
  int numharms = bVtPtr_->blockCount();
  const Linear::BlockVector bVf(&Vf, 2*numharms);
  int first = bVf.startBlock();

  Teuchos::RCP<const Linear::Vector> Vf_overlap;
  Linear::Vector freqDFDXtVf( *(hbBuilderPtr_->getSolutionMap()) );
  if ( overlapMap_ != Teuchos::null )
  {
    Teuchos::RCP<Linear::Vector> Vf_overlap_tmp = Teuchos::rcp( Xyce::Linear::createVector( *(hbBuilderPtr_->getSolutionMap()), *overlapMap_ ) ); 
    *Vf_overlap_tmp = Vf;
    Vf_overlap_tmp->importOverlap();
    Vf_overlap = Vf_overlap_tmp;
  }
  else
  {
    Vf_overlap = Teuchos::rcp( &Vf, false );
  }

  std::vector< Util::FreqVecEntry > freqFVector( freqNZLocalRows_.size() );
  if ( freqNZLocalRows_.size() )
  {
    for( int i = 0; i < (numharms + 1)/2; ++i )
    {
      // Perform matvec in the frequency domain.
      for (unsigned j = 0; j < freqNZLocalRows_.size(); j++)
      {
        freqFVector[j].val = std::complex<double>(0.0,0.0);
        freqFVector[j].lid = freqNZLocalRows_[j];
      }  
      for (unsigned j = 0; j < freqDFDXMatrix_[i].size(); j++)
      {
        int row_lid = freqNZLocalRowsMap_[freqDFDXMatrix_[i][j].row_lid];
        int col_lid = freqDFDXMatrix_[i][j].col_lid;
        freqFVector[row_lid].val += freqDFDXMatrix_[i][j].val * 
                                    std::complex<double>( (*Vf_overlap)[col_lid*2*numharms + 2*i],
                                                          (*Vf_overlap)[col_lid*2*numharms + 2*i+1] );
      }
      for (unsigned j = 0; j < freqFVector.size(); j++)
      {
        freqDFDXtVf[freqFVector[j].lid*2*numharms + 2*i] = freqFVector[j].val.real();
        freqDFDXtVf[freqFVector[j].lid*2*numharms + 2*i+1] = freqFVector[j].val.imag();
        if (i > 0)
        {
          freqDFDXtVf[freqFVector[j].lid*2*numharms + 2*(numharms-i)] = freqFVector[j].val.real();
          freqDFDXtVf[freqFVector[j].lid*2*numharms + 2*(numharms-i)+1] = -freqFVector[j].val.imag();
        } 
      }
    }
  }
  
  // Add the frequency domain contribution into permlindFdxV.
  // NOTE: update will add the assembled vector from freqDFDXtVf to 
  //       the underlying multivector of permlindFdxV, which all blocks
  //       view.
  permlindFdxV.update( 1.0, freqDFDXtVf, 0.0 );

  if ( !linAppdQdxPtr_->isEmpty() || !linAppdFdxPtr_->isEmpty())
  {
    Linear::MultiVector lindQdxV( *(bVtPtr_->blockPmap()), bVf.blockSize()/2 );
    Linear::MultiVector lindFdxV( *(bVtPtr_->blockPmap()), bVf.blockSize()/2 );

    int numMyRows = (bVtPtr_->blockPmap())->numLocalEntities();
    Linear::MultiVector permVf( *(bVtPtr_->blockPmap()), bVf.blockSize()/2 );

    // Now copy over data for the blocks that this processor owns.
    // NOTE:  Copy over all rows, some devices like VCVS may refer to columns
    //        which are not in the list of nonzero rows for the linear matrix storage.
    for (int row = 0; row < numMyRows; row++)
    {
      Linear::Vector & currBlock = bVf.block(first + row);

      // Insert zero-th value of the Fourier expansion, only need real value.
      permVf[0][row] = currBlock[0];
      for (int j=1; j<currBlock.localLength()/2; j++) 
      {
        permVf[j][row] = currBlock[j+1];
      } 
    }

    linAppdQdxPtr_->matvec(permVf, lindQdxV);
    linAppdFdxPtr_->matvec(permVf, lindFdxV);

    // Now copy over data for the blocks that this processor owns.
    for (std::vector<int>::iterator it = linNZRows_.begin(); it != linNZRows_.end(); it++)
    {
      int row = *it;

      if (!linAppdQdxPtr_->isEmpty())
      {
        Linear::Vector & currBlock = permlindQdxV.block(first + row);

        currBlock[0] = lindQdxV[0][row];
        currBlock[1] = 0.0;

        for (int j=1; j< (numharms + 1)/2; j++)
        {
          currBlock[2*j] = lindQdxV[2*j-1][row];
          currBlock[2*(numharms-j)] = lindQdxV[2*j-1][row];
     
          currBlock[2*j+1] = lindQdxV[2*j][row];
          currBlock[2*( numharms - j) + 1] = -lindQdxV[2*j][row];
        }
      }

      if (!linAppdFdxPtr_->isEmpty())
      {
        Linear::Vector & currBlockF = permlindFdxV.block(first + row);
        currBlockF[0] += lindFdxV[0][row];
        currBlockF[1] += 0.0;

        for (int j=1; j< (numharms + 1)/2; j++)
        {
          currBlockF[2*j] += lindFdxV[2*j-1][row];
          currBlockF[2*(numharms-j)] += lindFdxV[2*j-1][row];

          currBlockF[2*j+1] += lindFdxV[2*j][row];
          currBlockF[2*(numharms-j)+1] -= lindFdxV[2*j][row];
        }
      }
    }
  }

  // Put a barrier here for parallel.
  //(bVf.pmap()->pdsComm()).barrier();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : HBLoader::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool HBLoader::loadDAEVectors( Linear::Vector * Xf,
                               Linear::Vector * currX,
                               Linear::Vector * lastX,
                               Linear::Vector * S,
                               Linear::Vector * currS,
                               Linear::Vector * lastS,
                               Linear::Vector * dSdt,
                               Linear::Vector * Store,
                               Linear::Vector * currStore,
                               Linear::Vector * nextLeadFVectorPtr,
                               Linear::Vector * nextLeadQVectorPtr,
                               Linear::Vector * nextJunctionVVectorPtr,
                               Linear::Vector * Q,
                               Linear::Vector * F,
                               Linear::Vector * B,
                               Linear::Vector * dFdxdVp,
                               Linear::Vector * dQdxdVp,
                               int loadType )
{
  if (DEBUG_HB)
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  HBLoader::loadDAEVectors" << std::endl;
  }

  //Zero out vectors
  appVecPtr_->putScalar(0.0);
  appNextStaVecPtr_->putScalar(0.0);
  appCurrStaVecPtr_->putScalar(0.0);
  appLastStaVecPtr_->putScalar(0.0);
  Linear::Vector * appdSdt = builder_.createStateVector();

  appNextStoVecPtr_->putScalar(0.0);
  appCurrStoVecPtr_->putScalar(0.0);

  Linear::Vector * appQ = builder_.createVector();
  Linear::Vector * appF = builder_.createVector();
  Linear::Vector * appB = builder_.createVector();

  Linear::Vector * appdFdxdVp = builder_.createVector();
  Linear::Vector * appdQdxdVp = builder_.createVector();

  bXtPtr_->putScalar(0.0);

  Linear::Vector Xf_overlap( *(hbBuilderPtr_->getSolutionMap()), 
                             *(hbBuilderPtr_->getSolutionOverlapMap()) );
  Xf_overlap = *Xf;
  Xf_overlap.importOverlap();

  Linear::BlockVector & bXf = *dynamic_cast<Linear::BlockVector*>(Xf);
  permutedIFT(bXf, &*bXtPtr_);

  const std::vector<int> & augmentedLIDs = hbBuilderPtr_->getAugmentedLIDs();
  double fScalar = 1.0;

  if (hbOsc_)
  {
    double tmpfScalar = 0.0;   

    if (augmentedLIDs.size() )
      tmpfScalar = (*Xf)[(augmentedLIDs)[0 ]];

    Xf->pmap()->pdsComm().sumAll( &tmpfScalar, &fScalar, 1 );
  }


  // 12/8/06 tscoffe:   Note:  "b" at beginning of variable name means Linear::BlockVector
  Linear::BlockVector & bX = *bXtPtr_;
  Linear::BlockVector & bS = *dynamic_cast<Linear::BlockVector*>(S);
  Linear::BlockVector & bcurrS = *dynamic_cast<Linear::BlockVector*>(currS);
  Linear::BlockVector & blastS = *dynamic_cast<Linear::BlockVector*>(lastS);
  Linear::BlockVector & bdSdt = *dynamic_cast<Linear::BlockVector*>(dSdt);
  Linear::BlockVector & bStore = *dynamic_cast<Linear::BlockVector*>(Store);
  Linear::BlockVector & bcurrStore = *dynamic_cast<Linear::BlockVector*>(currStore);
  
  Linear::BlockVector & bNextLeadF = *dynamic_cast<Linear::BlockVector*>(nextLeadFVectorPtr);
  Linear::BlockVector & bLeadQ = *dynamic_cast<Linear::BlockVector*>(nextLeadQVectorPtr);
  Linear::BlockVector & bNextJunctionV = *dynamic_cast<Linear::BlockVector*>(nextJunctionVVectorPtr);
  
  Linear::BlockVector * bQ = dynamic_cast<Linear::BlockVector*>(Q);
  Linear::BlockVector * bF = dynamic_cast<Linear::BlockVector*>(F);
  Linear::BlockVector * bB = dynamic_cast<Linear::BlockVector*>(B);

  Linear::BlockVector * bdFdxdVp = dynamic_cast<Linear::BlockVector*>(dFdxdVp);
  Linear::BlockVector * bdQdxdVp = dynamic_cast<Linear::BlockVector*>(dQdxdVp);

  Teuchos::RCP<Linear::BlockVector> bQt = hbBuilderPtr_->createTimeDomainBlockVector();
  Teuchos::RCP<Linear::BlockVector> bFt = hbBuilderPtr_->createTimeDomainBlockVector();
  Teuchos::RCP<Linear::BlockVector> bBt = hbBuilderPtr_->createTimeDomainBlockVector();

  Teuchos::RCP<Linear::BlockVector> bdQdxdVpt = hbBuilderPtr_->createTimeDomainBlockVector();
  Teuchos::RCP<Linear::BlockVector> bdFdxdVpt = hbBuilderPtr_->createTimeDomainBlockVector();

  int BlockCount = bX.blockCount();

  if ((int)vecNLAppdQdxPtr_.size() != BlockCount)
  {
    vecNLAppdQdxPtr_.resize( BlockCount );
    vecNLAppdFdxPtr_.resize( BlockCount );
  }

  unsigned int nlQrows = nonlinQNZRows_.size();
  unsigned int nlFrows = nonlinFNZRows_.size();

  // Now perform implicit application of frequency domain Jacobian. 
  for( int i = 0; i < BlockCount; ++i )
  {
    if (DEBUG_HB)
    {
      Xyce::dout() << "Processing vectors for block " << i << " of " << BlockCount-1 << std::endl;
    }

    // Set Time for fast time scale somewhere
    deviceManager_.setFastTime( times_[i]/ fScalar);

    if (DEBUG_HB)
    {
      Xyce::dout() << "Calling updateSources on the appLoader" << std::endl;
    }

    //Update the sources
    appLoaderPtr_->updateSources();  // this is here to handle "fast" sources.

    *appVecPtr_ = bX.block(i);

    *appNextStaVecPtr_ = bS.block(i);
    *appCurrStaVecPtr_ = bcurrS.block(i);
    *appLastStaVecPtr_ = blastS.block(i);
    *appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);
    *appCurrStoVecPtr_ = bcurrStore.block(i);
    *appNextLeadFVecPtr_  = bNextLeadF.block(i);
    *appLeadQVecPtr_      =  bLeadQ.block(i);
    *appNextJunctionVVecPtr_  =  bNextJunctionV.block(i);
    
    if (DEBUG_HB)
    {
      Xyce::dout() << "Updating State for block " << i << " of " << BlockCount-1 << std::endl;
    }

    // Note: This updateState call is here (instead of in the
    // HBLoader::updateState function) because it has to be called
    // for the same fast time point.
    appLoaderPtr_->updateState
      ( &*appVecPtr_,
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_ , &*appLastStaVecPtr_,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_,
        Xyce::Device::NONLINEAR_FREQ );

    bS.block(i) = *appNextStaVecPtr_;
    bcurrS.block(i) = *appCurrStaVecPtr_;
    blastS.block(i) = *appLastStaVecPtr_;
    bStore.block(i) = *appNextStoVecPtr_;
    bcurrStore.block(i) = *appCurrStoVecPtr_;

    if (DEBUG_HB)
    {
      Xyce::dout() << "Calling loadDAEVectors on the appLoader" << std::endl;
    }

    // This has to be done because the app loader does NOT zero these vectors out.
    appQ->putScalar(0.0);
    appF->putScalar(0.0);
    appB->putScalar(0.0);
    appdFdxdVp->putScalar(0.0);
    appdQdxdVp->putScalar(0.0);

    appLoaderPtr_->loadDAEVectors
      ( &*appVecPtr_,
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appVecPtr_,  // note, this is a placeholder! ERK
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_, &*appLastStaVecPtr_, &*appdSdt,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_, 
        &*appNextLeadFVecPtr_, &*appLeadQVecPtr_, 
        &*appNextJunctionVVecPtr_, 
        appQ, appF, appB,
        appdFdxdVp, appdQdxdVp,
        Xyce::Device::NONLINEAR_FREQ );

    // Load the sources into appB.   
    appLoaderPtr_->loadBVectorsforSources();

    appB->fillComplete();

    bQt->block(i) = *appQ;
    bFt->block(i) = *appF;
    bBt->block(i) = *appB;

    bdQdxdVpt->block(i) = *appdQdxdVp;
    bdFdxdVpt->block(i) = *appdFdxdVp;


    bNextLeadF.block(i) = *appNextLeadFVecPtr_;  // lead currents loaded into lead current vector.
    bLeadQ.block(i) = *appLeadQVecPtr_;

    // Store the time domain Jacobian for future use.
    appdQdxPtr_->put(0.0);
    appdFdxPtr_->put(0.0);

    // Load dQdx and dFdx into the storage location for this time point, for linear devices.
    if (i == 0 && Teuchos::is_null(linAppdQdxPtr_) && Teuchos::is_null(linAppdFdxPtr_))
    {
      appLoaderPtr_->loadDAEMatrices( &*appVecPtr_, &*appNextStaVecPtr_, &*appdSdt, &*appNextStoVecPtr_, &*appdQdxPtr_, &*appdFdxPtr_, Xyce::Device::LINEAR_FREQ);

      // Copy over matrix values into linear storage
      linAppdQdxPtr_ = Teuchos::rcp( new Xyce::Linear::FilteredMatrix( &*appdQdxPtr_, appVecPtr_->pmap(), false ) );
      linAppdFdxPtr_ = Teuchos::rcp( new Xyce::Linear::FilteredMatrix( &*appdFdxPtr_, appVecPtr_->pmap(), false ) );

      // Create a vector containing the unique NZ rows in dQdx and dFdx.
      linNZRows_.clear();
      linNZRows_.insert( linNZRows_.end(), (linAppdQdxPtr_->getNZRows()).begin(), (linAppdQdxPtr_->getNZRows()).end() );
      linNZRows_.insert( linNZRows_.end(), (linAppdFdxPtr_->getNZRows()).begin(), (linAppdFdxPtr_->getNZRows()).end() );
      std::sort( linNZRows_.begin(), linNZRows_.end() );
      linNZRows_.erase( std::unique( linNZRows_.begin(), linNZRows_.end() ), linNZRows_.end() );

      appdQdxPtr_->put(0.0);
      appdFdxPtr_->put(0.0);
    }

    // Load dQdx and dFdx into the storage location for this time point, for nonlinear devices.
    appLoaderPtr_->loadDAEMatrices( &*appVecPtr_, &*appNextStaVecPtr_, &*appdSdt, &*appNextStoVecPtr_, &*appdQdxPtr_, &*appdFdxPtr_, Xyce::Device::NONLINEAR_FREQ);

    // Copy over matrix values into storage for nonlinear.
    // Try to reuse filtered matrix objects whenever possible, for efficiency.
    if (Teuchos::is_null(vecNLAppdQdxPtr_[i]))
    {
      vecNLAppdQdxPtr_[i] = Teuchos::rcp( new Xyce::Linear::FilteredMatrix( &*appdQdxPtr_, appVecPtr_->pmap(), false ) );

      // Add new rows to vector containing the unique NZ rows in dQdx.
      nonlinQNZRows_.insert( nonlinQNZRows_.end(), (vecNLAppdQdxPtr_[i]->getNZRows()).begin(), (vecNLAppdQdxPtr_[i]->getNZRows()).end() );
    }
    else
    {
      // Attempt to use previous NZ entries to filter new appdQdxPtr_, recreate filtered matrix if necessary
      bool reset = vecNLAppdQdxPtr_[i]->filterMatrix( &*appdQdxPtr_, appVecPtr_->pmap(), false );
      if (reset)
        nonlinQNZRows_.insert( nonlinQNZRows_.end(), (vecNLAppdQdxPtr_[i]->getNZRows()).begin(), (vecNLAppdQdxPtr_[i]->getNZRows()).end() );
    }

    if (Teuchos::is_null(vecNLAppdFdxPtr_[i]))
    {
      vecNLAppdFdxPtr_[i] = Teuchos::rcp( new Xyce::Linear::FilteredMatrix( &*appdFdxPtr_, appVecPtr_->pmap(), false ) );

      // Add new rows to vector containing the unique NZ rows in dFdx.
      nonlinFNZRows_.insert( nonlinFNZRows_.end(), (vecNLAppdFdxPtr_[i]->getNZRows()).begin(), (vecNLAppdFdxPtr_[i]->getNZRows()).end() );
    }
    else
    {
      // Attempt to use previous NZ entries to filter new appdFdxPtr_, recreate filtered matrix if necessary
      bool reset = vecNLAppdFdxPtr_[i]->filterMatrix( &*appdFdxPtr_, appVecPtr_->pmap(), false ); 
      if (reset)
        nonlinFNZRows_.insert( nonlinFNZRows_.end(), (vecNLAppdFdxPtr_[i]->getNZRows()).begin(), (vecNLAppdFdxPtr_[i]->getNZRows()).end() );
    }
  }

  // Generate unique list of NZ rows for dQdx and dFdx. 
  if (nlQrows < nonlinQNZRows_.size())
  {
    std::sort( nonlinQNZRows_.begin(), nonlinQNZRows_.end() );
    nonlinQNZRows_.erase( std::unique( nonlinQNZRows_.begin(), nonlinQNZRows_.end() ), nonlinQNZRows_.end() );
  }
  if (nlFrows < nonlinFNZRows_.size())
  {
    std::sort( nonlinFNZRows_.begin(), nonlinFNZRows_.end() );
    nonlinFNZRows_.erase( std::unique( nonlinFNZRows_.begin(), nonlinFNZRows_.end() ), nonlinFNZRows_.end() );
  }

  // Load frequency-domain devices.
  // NOTE:  Use overlapped vector size because some devices may reference ground or ghost nodes.
  int localOverlapGndN = appVecPtr_->omap()->numLocalEntities();
  int indexBase = appVecPtr_->omap()->indexBase();
  int M = (BlockCount - 1) / 2;

  // First convert real-equivalent form of the solution vector Xf to a complex-valued vector for loading.
  freqBVector_.resize(M+1);
  freqDFDXMatrix_.resize(M+1);

  if ( !freqLoadAnalysisDone_ )
  {
    for (int i = 0; i <= M; i++)
    {
      double freq = freqs_[M + i] * fScalar;

      // Copy values from Xf to Xf_complex vector for solver.
      std::vector< std::complex<double> > Xf_complex( localOverlapGndN );
      for (int nB = 0; nB < localOverlapGndN+indexBase; nB++)
      {
        std::complex<double> val( Xf_overlap[nB*2*BlockCount + 2*i], Xf_overlap[nB*2*BlockCount + 2*i+1] );
        Xf_complex[nB] = val;
      }
      if (indexBase == -1)
      Xf_complex[localOverlapGndN-1] = std::complex<double>(0.0, 0.0);

      // Now load freqBVector and freqDFDXMatrix_ for this frequency.
      std::vector< Util::FreqVecEntry > tmpFreqFVector;

      // Call updateFDIntermediate vars so the device can compute any
      // quantities it needs prior to loading vectors
      Teuchos::rcp_dynamic_cast<CktLoader>(appLoaderPtr_)->updateFDIntermediateVars( freq, &Xf_complex[0]);

      // Call application loader for frequency-domain contributions
      Teuchos::rcp_dynamic_cast<CktLoader>(appLoaderPtr_)->loadFreqDAEVectors( freq, &Xf_complex[0], 
                                                           tmpFreqFVector, freqBVector_[i] );

      // Perform some analysis of nonzero entries for the first harmonic, reuse for other harmonics.
      // This method will determine local and nonlocal nonzero rows, as well as communication PIDs.
      // NOTE:  This method will remove the row entries associated with the ground node. 
      if ( i == 0 )
      {
        compNZRowsAndCommPIDs( tmpFreqFVector, freqBVector_[i] );
      }

      // Load entries of DFDX matrix in frequency domain.
      Teuchos::rcp_dynamic_cast<CktLoader>(appLoaderPtr_)->loadFreqDAEMatrices( freq, &Xf_complex[0], freqDFDXMatrix_[i] );

      // Consolidate DFDX matrix entries not owned by this processor for this harmonic.
      std::vector< Util::FreqMatEntry > offProcFreqDFDXMatrix, newFreqDFDXMatrix;
      consolidateMatrixEntries( offProcLocalRows_, freqDFDXMatrix_[i], offProcFreqDFDXMatrix );

      // Send and/or receive off processor DFDX matrix entries for this harmonic.
      sendReceiveMatrixEntries( offProcFreqDFDXMatrix, newFreqDFDXMatrix );

      // Consolidate DFDX matrix entries owned by this processor for this harmonic.
      // NOTE:  If any entries were received from off processor in newFreqDFDXMatrix, 
      //        this will append the local matrix entries.
      // NOTE 2:  This method will remove the column entries associated with the ground node. 
      consolidateMatrixEntries( freqNZLocalRows_, freqDFDXMatrix_[i], newFreqDFDXMatrix, false );
      freqDFDXMatrix_[i] = newFreqDFDXMatrix;
    }

    // Create the permuted B vector to add to the transformed time-domain B vector.
    createPermFreqBVector( freqBVector_, permFreqBVector_ );
  }

  permutedFFT2(*bBt, bB);

  permutedFFT2(*bQt, bQ); 

  permutedFFT2(*bFt, bF);

  permutedFFT2(*bdQdxdVpt, bdQdxdVp);

  permutedFFT2(*bdFdxdVpt, bdFdxdVp);

  // Clean up temporary vectors before moving on.
  delete appdSdt;
  delete appQ; 
  delete appF;
  delete appB;
  delete appdFdxdVp;
  delete appdQdxdVp; 

  bLeadCurrentVecFreqPtr_->putScalar(0.0);
  bLeadCurrentQVecFreqPtr_->putScalar(0.0);
  permutedFFT2(bNextLeadF, &*bLeadCurrentVecFreqPtr_);
  permutedFFT2(bLeadQ, &*bLeadCurrentQVecFreqPtr_);

  // Add in linear contributions into vectors.
  Teuchos::RCP<Linear::BlockVector> permlindQdxX = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();
  Teuchos::RCP<Linear::BlockVector> permlindFdxX = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();

  applyLinearMatrices (*Xf, *permlindQdxX, *permlindFdxX);

  if(!linAppdFdxPtr_->isEmpty() || numGlobalFreqRows_)
  {
    bF->update (1.0, *permlindFdxX, 1.0);
  }
 
  if(!linAppdQdxPtr_->isEmpty())
  {
    bQ->update (1.0, *permlindQdxX, 1.0);
  }

  if(!Teuchos::is_null( permFreqBVector_ ))
  {
    bB->update (1.0, *permFreqBVector_, 1.0);
  }

/*
    Xyce::dout() << "HB X Vector" << std::endl;
    bX.print(std::cout);
    Xyce::dout() << "HB Q Vector" << std::endl;
    bQ->print(std::cout);
    Xyce::dout() << "HB F Vector" << std::endl;
    bF->print(std::cout);
    Xyce::dout() << "HB B Vector" << std::endl;
    bB->print(std::cout);
*/

  int blockCount = bXf.blockCount();
  int blockSize = bXf.block(0).globalLength();
  
  int posFreq = (freqs_.size()-1)/2;
  double omega = 2.0 * M_PI * freqs_[posFreq] * fScalar;

  if (hbOsc_) 
    bQPtr_->putScalar( 0.0);

  for( int i = 0; i < blockCount; ++i )
  {
    // Create work vectors from the current frequency block vector
    // NOTE:  This needs to be done for each block to make sure that the
    //        map is the same as the bF block.
    Linear::Vector * QVecPtr = (bQ->block(i)).clone();
    Linear::Vector& QVec = *QVecPtr;
    Linear::Vector& freqVec = bQ->block(i);

    Linear::Vector * dQdxdVpVecPtr = (bdQdxdVp->block(i)).clone();
    Linear::Vector& dQdxdVpVec = *dQdxdVpVecPtr;
    Linear::Vector& freqVec1 = bdQdxdVp->block(i);

    if (hbOsc_)
      bQPtr_->block(i) = bQ->block(i);

    // Only one processor owns each block of the frequency-domain vector
    if (freqVec.localLength() > 0)
    {
      omega = 2.0 * M_PI * freqs_[posFreq] * fScalar;

      QVec[0] = -freqVec[1]*omega;
      QVec[1] = freqVec[0]*omega;

      dQdxdVpVec[0] = -freqVec1[1]*omega;
      dQdxdVpVec[1] = freqVec1[0]*omega;

      for (int j=1; j < (blockSize/2+1)/2; ++j)
      {
        omega = 2.0 * M_PI * freqs_[posFreq + j] * fScalar;

        QVec[2*j] = -freqVec[2*j+1]*omega;
        QVec[2*(blockSize/2-j)] = -freqVec[2*j+1]*omega;

        QVec[2*j+1] = freqVec[2*j]*omega;
        QVec[2*(blockSize/2-j)+1] = -freqVec[2*j]*omega;

        dQdxdVpVec[2*j] = -freqVec1[2*j+1]*omega;
        dQdxdVpVec[2*(blockSize/2-j)] = -freqVec1[2*j+1]*omega;

        dQdxdVpVec[2*j+1] = freqVec1[2*j]*omega;
        dQdxdVpVec[2*(blockSize/2-j)+1] = -freqVec1[2*j]*omega;
      }
    }

    bF->block(i).update(1.0, QVec , 1.0);
    bdFdxdVp->block(i).update(1.0, dQdxdVpVec, 1.0);

    delete QVecPtr;
    delete dQdxdVpVecPtr;
  }

//  blockCount = bStoreVecFreqPtr_->blockCount();
//  blockSize =  bStoreVecFreqPtr_->blockSize();

  // 
  // take care of the lead current vector 
  blockCount = bLeadCurrentVecFreqPtr_->blockCount();
  blockSize =  bLeadCurrentVecFreqPtr_->blockSize();

  for( int i = 0; i < blockCount; ++i )
  {
    // Create work vectors from the current frequency block vector
    // NOTE:  This needs to be done for each block to make sure that the
    //        map is the same as the bF block.
    Linear::Vector * leadCurrdQdtVecPtr = (bLeadCurrentQVecFreqPtr_->block(i)).clone();
    Linear::Vector& leadCurrdQdtVec = *leadCurrdQdtVecPtr;
    Linear::Vector& leadCurrQVec = bLeadCurrentQVecFreqPtr_->block(i);

    // Only one processor owns each block of the frequency-domain vector
    if (leadCurrQVec.localLength() > 0)
    {

      omega = 2.0 * M_PI * freqs_[posFreq] * fScalar;
      leadCurrdQdtVec[0] = leadCurrQVec[1]*omega;
      leadCurrdQdtVec[1] = leadCurrQVec[0]*omega;

      for (int j=1; j < (blockSize/2+1)/2; ++j)
      {
        omega = 2.0 * M_PI * freqs_[posFreq+ j] * fScalar;

        leadCurrdQdtVec[2*j] = -leadCurrQVec[2*j+1]*omega;
        leadCurrdQdtVec[2*(blockSize/2-j)] = -leadCurrQVec[2*j+1]*omega;

        leadCurrdQdtVec[2*j+1] = leadCurrQVec[2*j]*omega;
        leadCurrdQdtVec[2*(blockSize/2-j)+1] = -leadCurrQVec[2*j]*omega;
      }
    }

    bLeadCurrentVecFreqPtr_->block(i).update(1.0, leadCurrdQdtVec, 1.0);

    delete leadCurrdQdtVecPtr;
  }

  if (DEBUG_HB)
  {
    Xyce::dout() << "HB X Vector" << std::endl;
    bX.print(std::cout);
    //  Xyce::dout() << "HB S Vector" << std::endl;
    //  bS.print(Xyce::dout());
    //  Xyce::dout() << "HB dSdt Vector" << std::endl;
    //  bdSdt.print(Xyce::dout());
    Xyce::dout() << "HB Store Vector" << std::endl;
    bStore.print(std::cout);
    Xyce::dout() << "HB Q Vector" << std::endl;
    bQ->print(std::cout);
    Xyce::dout() << "HB F Vector" << std::endl;
    bF->print(std::cout);
    Xyce::dout() << "HB bdFdxdVp Vector" << std::endl;
    bdFdxdVp->print(std::cout);
    Xyce::dout() << "HB bdQdxdVp Vector" << std::endl;
    bdQdxdVp->print(std::cout);
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  if (hbOsc_)
  {

//    if (augmentedLIDs.size() )
//      (*bF)[(augmentedLIDs)[0 ]] = fScalar - 1.0;
    double refValue = 0.0;
    double tmpValue= 0.0;
    Linear::Vector& freqVec = bXf.block(refID_);

    if (freqVec.localLength() > 0)
    {
       tmpValue = freqVec[3];
    }

    Xf->pmap()->pdsComm().sumAll( &tmpValue, &refValue, 1 );

    if (augmentedLIDs.size() )
      (*bF)[(augmentedLIDs)[0 ]] = refValue;
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : HBLoader::loadDeviceErrorWeightMask
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 11/1/2014
//-----------------------------------------------------------------------------
bool
HBLoader::loadDeviceErrorWeightMask(
  Linear::Vector *      deviceMask) const
{
#if 0
  Linear::BlockVector & bDevMask = *dynamic_cast<Linear::BlockVector*>(deviceMask);

#if 0
  Xyce::dout() << "HBLoader::loadDeviceErrorWeightMask.  Original (nonblock) deviceMask.size = "
    << appVecPtr_->globalLength() <<std::endl;
#endif

  appVecPtr_->putScalar(1.0);
  bool returnValue = deviceManager_.loadErrorWeightMask(&*appVecPtr_);

  int blockCount = bDevMask.blockCount();
  int blockSize =  bDevMask.blockSize();

#if 0
  Xyce::dout() << "bDevMask.blockCount = "<< blockCount <<std::endl;
  Xyce::dout() << "bDevMask.blockSize = "<< blockSize <<std::endl;
  appVecPtr_->print(Xyce::dout());
#endif

  //Teuchos::RCP<N_PDS_ParMap> baseMap = Teuchos::rcp_const_cast<N_PDS_ParMap>( hbBuilderPtr_->getBaseStoreMap() );

  for( int i = 0; i < blockCount; ++i )
  {
    // See if this variable is owned by the local processor.
    // If so, this processor owns the entire j-th block of the vector
    //int lid = baseMap->globalToLocalIndex( i );

    Linear::Vector& localVecRef =  bDevMask.block(i);

    for (int j=0;j<blockSize;++j)
    {
      localVecRef[j] = (*appVecPtr_)[i];
    }
  }

#if 0
  bDevMask.print(Xyce::dout());
#endif

  return returnValue;
#else
  return true;
#endif
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::permutedFFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/05/08
//---------------------------------------------------------------------------
void HBLoader::permutedFFT(const Linear::BlockVector & xt, Linear::BlockVector * xf, std::vector<int>* lids)
{
  // Call the function to compute the permuted FFT from the block system helper functions.
  computePermutedDFT( *dftInterface_, xt, xf, lids );
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::permutedFFT2
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/05/08
//---------------------------------------------------------------------------
void HBLoader::permutedFFT2(const Linear::BlockVector & xt, Linear::BlockVector * xf)
{
  // Call the function to compute the permuted FFT from the block system helper functions.
  computePermutedDFT2( *dftInterface_, xt, xf );
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::permutedIFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 09/05/08
//---------------------------------------------------------------------------
void HBLoader::permutedIFT(const Linear::BlockVector & xf, Linear::BlockVector * xt, int numTimePts_ )
{
  // Call the function to compute the permuted IFT from the block system helper functions.
  computePermutedIFT( *dftInterface_, xf, xt, numTimePts_ );  
} 

//-----------------------------------------------------------------------------
// Function      : HBLoader::getVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015
//---------------------------------------------------------------------------
bool HBLoader::getVoltageLimiterStatus()
{
  return appLoaderPtr_->getVoltageLimiterStatus();
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::setVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015 
//---------------------------------------------------------------------------
void HBLoader::setVoltageLimiterStatus(bool voltageLimterStatus)
{
  return appLoaderPtr_->setVoltageLimiterStatus(voltageLimterStatus);
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::consolidateMatrixEntries
// Purpose       :
// Special Notes : overlapIDs are necessary when processing IDs for off processor
//               : rows.  For local rows overlapIDs should be false.
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 8/21/2017
//---------------------------------------------------------------------------
void HBLoader::consolidateMatrixEntries( const std::vector<int>& nzRows,
                                         const std::vector< Util::FreqMatEntry >& matrixEntries,
                                         std::vector< Util::FreqMatEntry >& consolidatedEntries,
                                         bool overlapIDs )
{
  // Consolidate matrix entries for this harmonic.
  // Convert offset lids to column lids.
  int numEntries = 0;
  double * values = 0;
  int * indices = 0;
  N_PDS_ParMap * colmap = appdFdxPtr_->getColMap( *builder_.getPDSComm() );
  N_PDS_ParMap * ocolmap = appdFdxPtr_->getOverlapColMap( *builder_.getPDSComm() );

  std::vector<int>::const_iterator lid_it = nzRows.begin();

  for (int idx=0; lid_it != nzRows.end(); lid_it++, idx++)
  {
    // Get the current row id and extract row entries from the matrix to
    // convert the offset col_ids to column local ids.
    int curr_lid = *lid_it;
    appdFdxPtr_->extractLocalRowView( curr_lid, numEntries, values, indices);

    std::set<int> col_id_set;
    std::vector< Util::FreqMatEntry > tmpFreqDFDXMatrix;
    std::vector< Util::FreqMatEntry >::const_iterator it = matrixEntries.begin();
    for ( ; it != matrixEntries.end(); it++)
    {
      if ( (it->row_lid == curr_lid) && (it->col_lid != -1) )
      {
        int col_id = indices[ it->col_lid ];
        if (!overlapIDs)
        {
          int gid = ocolmap->localToGlobalIndex( col_id );
          col_id = colmap->globalToLocalIndex( gid );
        }
        std::pair<std::set<int>::iterator, bool> ret = col_id_set.insert( col_id );
        if (ret.second == false)
        {
          std::vector< Util::FreqMatEntry >::iterator it2 = tmpFreqDFDXMatrix.begin();
          std::vector< Util::FreqMatEntry >::iterator it2_end = tmpFreqDFDXMatrix.end();
          for ( ; it2 != it2_end; it2++ )
          {
            if ( it2->col_lid == col_id )
            {
              it2->val += it->val;
            }
          }
        }
        else
        {
          tmpFreqDFDXMatrix.push_back( *it );
          tmpFreqDFDXMatrix.back().col_lid = col_id;
        }
      }
    }
    consolidatedEntries.insert( consolidatedEntries.end(), tmpFreqDFDXMatrix.begin(),
                                tmpFreqDFDXMatrix.end() );
  }
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::sendReceiveMatrixEntries
// Purpose       :
// Special Notes : Input vector contains only entries that need to be
//                 communicated to other processors.
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 8/21/2017
//---------------------------------------------------------------------------
void HBLoader::sendReceiveMatrixEntries( const std::vector< Util::FreqMatEntry >& sendMatrixEntries,
                                         std::vector< Util::FreqMatEntry >& recvMatrixEntries )
{
  int numProcs = appVecPtr_->pmap()->pdsComm().numProc();
  int myProc = appVecPtr_->pmap()->pdsComm().procID();
  N_PDS_ParMap * colmap = appdFdxPtr_->getColMap( *builder_.getPDSComm() );
  N_PDS_ParMap * ocolmap = appdFdxPtr_->getOverlapColMap( *builder_.getPDSComm() );

  if ( totalNZOffProcRows_ )
  {
    for (int sendProc=0; sendProc<numProcs; sendProc++)
    {
      // This processor is packing up entries for the other processors.
      if ( sendProc == myProc )
      {
        int prevProc = -1;
        std::vector<int> row_col_gids;
        std::vector<double> vals;
        std::vector<int>::iterator send_it = offProcLocalRows_.begin();
        std::vector<int>::iterator send_pid_it = offProcLocalRowsRecvPIDs_.begin();
        for ( ; send_it != offProcLocalRows_.end(); send_it++, send_pid_it++ )
        {
          // Initialize the receiving processor.
          if (prevProc == -1)
            prevProc = *send_pid_it;

          // If the receiving processor has changed, send the current data.
          if (prevProc != *send_pid_it)
          {
            // Send the data that is in the row_col_gids and vals vector to prevProc.
            int len = row_col_gids.size();
            appVecPtr_->pmap()->pdsComm().send( &len, 1, prevProc );
            if ( len )
            {
              appVecPtr_->pmap()->pdsComm().send( &row_col_gids[0], len, prevProc );
              appVecPtr_->pmap()->pdsComm().send( &vals[0], len, prevProc );
            }
            row_col_gids.clear();
            vals.clear();
            prevProc = *send_pid_it;
          }
          int globalRow = appVecPtr_->omap()->localToGlobalIndex( *send_it );
          std::vector< Util::FreqMatEntry >::const_iterator mat_it = sendMatrixEntries.begin();
          for ( ; mat_it != sendMatrixEntries.end(); mat_it++ )
          {
            if ( mat_it->row_lid == *send_it )
            {
              // Don't send row entries that are associated with the ground node.
              if ( mat_it->col_lid != -1 )
              {
                // convert lids to gids and store in appropriate vectors.
                int globalCol = ocolmap->localToGlobalIndex( mat_it->col_lid );
                row_col_gids.push_back( globalRow );
                row_col_gids.push_back( globalCol );
                vals.push_back( mat_it->val.real() );
                vals.push_back( mat_it->val.imag() );
              }
            }
          }
        }
        // Final send of the data that is in the row_col_gids and vals vector to prevProc.
        int len = row_col_gids.size();
        if ( prevProc != -1 )
        {
          appVecPtr_->pmap()->pdsComm().send( &len, 1, prevProc );
          if ( len )
          {
            appVecPtr_->pmap()->pdsComm().send( &row_col_gids[0], len, prevProc );
            appVecPtr_->pmap()->pdsComm().send( &vals[0], len, prevProc );
          }
        }
      }
      else
      {
        // Wait for data if you are supposed to receive something.
        std::vector<int>::iterator recv_it = std::find( offProcNonlocalRowsSendPIDs_.begin(),
                                                        offProcNonlocalRowsSendPIDs_.end(), sendProc );
        if ( recv_it != offProcNonlocalRowsSendPIDs_.end() )
        {
          // Wait to receive data from sendProc.
          int len = 0;
          appVecPtr_->pmap()->pdsComm().recv( &len, 1, sendProc );
          if ( len )
          {  
            std::vector<int> row_col_gids( len );
            std::vector<double> vals( len );
            appVecPtr_->pmap()->pdsComm().recv( &row_col_gids[0], len, sendProc );
            appVecPtr_->pmap()->pdsComm().recv( &vals[0], len, sendProc );
            
            for (int i=0; i<len; i+=2)
            {
              // If the column is not the ground node, insert this entry.
              Util::FreqMatEntry newEntry;
              newEntry.row_lid = appVecPtr_->pmap()->globalToLocalIndex( row_col_gids[i] );
              newEntry.col_lid = colmap->globalToLocalIndex( row_col_gids[i+1] );
              newEntry.val = std::complex<double>(vals[i], vals[i+1]);
              recvMatrixEntries.push_back( newEntry );
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBLoader::compNZRowsAndCommPIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 8/21/2017
//---------------------------------------------------------------------------
void HBLoader::compNZRowsAndCommPIDs( const std::vector< Util::FreqVecEntry >& vectorEntries,
                                      const std::vector< Util::FreqVecEntry >& bVecEntries )
{
  int numProcs = appVecPtr_->pmap()->pdsComm().numProc();
  int myProc = appVecPtr_->pmap()->pdsComm().procID();

  int row_groundID = appVecPtr_->omap()->globalToLocalIndex( -1 );

  // Collect all of the off processor ids in the b vector
  std::vector< Util::FreqVecEntry >::const_iterator vE_it = bVecEntries.begin();
  for (; vE_it != bVecEntries.end(); vE_it++)
  {
    if ( (appVecPtr_->pmap()->localToGlobalIndex( vE_it->lid ) == -1) && (vE_it->lid != row_groundID) )
      offProcBVecLIDs_.push_back( vE_it->lid ); 
  }
  std::sort( offProcBVecLIDs_.begin(), offProcBVecLIDs_.end() );
  offProcBVecLIDs_.erase( std::unique( offProcBVecLIDs_.begin(), offProcBVecLIDs_.end() ), offProcBVecLIDs_.end() ); 

  std::vector<int> numOffProcBVecLIDs( numProcs ), sumOffProcBVecLIDs( numProcs );
  numOffProcBVecLIDs[ myProc ] = offProcBVecLIDs_.size();
  appVecPtr_->pmap()->pdsComm().sumAll( &numOffProcBVecLIDs[0], &sumOffProcBVecLIDs[0], numProcs );

  // Collect all the GIDs being communicated for the b vector load 
  totalOffProcBVecLIDs_ = std::accumulate( sumOffProcBVecLIDs.begin(), sumOffProcBVecLIDs.end(), 0 );
  std::vector<int> tmpOffProcGIDs( totalOffProcBVecLIDs_, 0 ), allOffProcGIDs( totalOffProcBVecLIDs_, 0 );
  for ( int p=0, ptr=0; p<numProcs; ++p )
  {
    if ( myProc == p )
    {
      for ( int idx=0; idx<sumOffProcBVecLIDs[p]; ++idx, ++ptr )
        tmpOffProcGIDs[ ptr ] = appVecPtr_->omap()->localToGlobalIndex( offProcBVecLIDs_[idx] );
    }
    else
      ptr += sumOffProcBVecLIDs[p];
   
    if ( p > 0 )
      sumOffProcBVecLIDs[p] += sumOffProcBVecLIDs[p-1];
  } 
  appVecPtr_->pmap()->pdsComm().sumAll( &tmpOffProcGIDs[0], &allOffProcGIDs[0], totalOffProcBVecLIDs_ );

  // Convert the GIDs being communicated to the PIDs that own them
  // Store the PIDs that b vector information needs to be received from.
  std::vector<int> tmpOffProcPIDs( totalOffProcBVecLIDs_, 0 ), allOffProcPIDs( totalOffProcBVecLIDs_, 0 );

  for ( int p=0; p<numProcs; ++p )
  {
    int idx = 0;
    if ( p > 0 )
      idx = sumOffProcBVecLIDs[ p-1 ];
    for ( ; idx < sumOffProcBVecLIDs[ p ]; ++idx )
    {
      int lid = appVecPtr_->pmap()->globalToLocalIndex( allOffProcGIDs[idx] );
      if ( lid > -1 )  
      {
        tmpOffProcPIDs[ idx ] = myProc;
        offProcBVecSendPIDs_.push_back( p );
        offProcBVecSendLIDs_.push_back( lid );
      }
    } 
  }
  appVecPtr_->pmap()->pdsComm().sumAll( &tmpOffProcPIDs[0], &allOffProcPIDs[0], totalOffProcBVecLIDs_ );

  // Resize the processors we are receiving the data from.
  offProcBVecPIDs_.resize( offProcBVecLIDs_.size() );

  // Store the PIDs that need to be sent to.
  int gidPtr = 0;
  if ( myProc > 0 )
    gidPtr = sumOffProcBVecLIDs[ myProc-1 ];

  for ( int idx=0; idx<offProcBVecLIDs_.size(); ++idx )
  { 
    offProcBVecPIDs_[ idx ] = allOffProcPIDs[ gidPtr + idx ];
  }

  // Now set up off-processor communication for the f vector and Jacobian. 
  vE_it = vectorEntries.begin();
  for (; vE_it != vectorEntries.end(); vE_it++)
  {
    if ( appVecPtr_->pmap()->localToGlobalIndex( vE_it->lid ) != -1 )
      freqNZLocalRows_.push_back( vE_it->lid );
    else if ( vE_it->lid != row_groundID )
      offProcLocalRows_.push_back( vE_it->lid );
  }
  std::sort( offProcLocalRows_.begin(), offProcLocalRows_.end() );
  offProcLocalRows_.erase( std::unique( offProcLocalRows_.begin(), offProcLocalRows_.end() ), offProcLocalRows_.end() ); 

  // Now communicate the number off processor rows that need to be sent after loadFreqDAEMatrices
  std::vector<int> numNZOffProcRows( numProcs ), sumNZOffProcRows( numProcs );
  numNZOffProcRows[ myProc ] = offProcLocalRows_.size(); 
  appVecPtr_->pmap()->pdsComm().sumAll( &numNZOffProcRows[0], &sumNZOffProcRows[0], numProcs );

  // Determine which global rows are being communicated and to which processor
  totalNZOffProcRows_ = std::accumulate( sumNZOffProcRows.begin(), sumNZOffProcRows.end(), 0 );

  if (totalNZOffProcRows_)
  {
    std::vector<int> tmpOffProcRows( totalNZOffProcRows_ ), allOffProcRows( totalNZOffProcRows_ );
    std::vector<int> tmpOffProcRowsSendPIDs( totalNZOffProcRows_ ), allOffProcRowsSendPIDs( totalNZOffProcRows_ );
    std::vector<int> tmpOffProcRowsRecvPIDs( totalNZOffProcRows_ ), allOffProcRowsRecvPIDs( totalNZOffProcRows_ );
    for (int idx=0, proc=0; proc<numProcs; proc++)
    {
      if ( proc == myProc )
      {
        for ( std::vector<int>::iterator it = offProcLocalRows_.begin(); it != offProcLocalRows_.end(); it++, idx++ )
        {
          tmpOffProcRows[idx] = appVecPtr_->omap()->localToGlobalIndex( *it );
          tmpOffProcRowsSendPIDs[idx] = proc;
        }
      }
      else
        idx += sumNZOffProcRows[proc];
    }
    appVecPtr_->pmap()->pdsComm().sumAll( &tmpOffProcRows[0], &allOffProcRows[0], totalNZOffProcRows_ );
    appVecPtr_->pmap()->pdsComm().sumAll( &tmpOffProcRowsSendPIDs[0], &allOffProcRowsSendPIDs[0], totalNZOffProcRows_ );

    std::vector<int>::iterator it=allOffProcRows.begin();
    for (int idx=0; it != allOffProcRows.end(); it++, idx++)
    {
      // Check if this processor owns this row, by looking at the non-overlapped map.
      // Claim the row, so other processors know to send it to this processor, save the local row.
      int lid = appVecPtr_->pmap()->globalToLocalIndex( *it );
      if ( lid != -1 )
      {
        tmpOffProcRowsRecvPIDs[ idx ] = myProc;
        freqNZLocalRows_.push_back( lid );
      }
    }
    appVecPtr_->pmap()->pdsComm().sumAll( &tmpOffProcRowsRecvPIDs[0], &allOffProcRowsRecvPIDs[0], totalNZOffProcRows_ );

    // Determine send and receive processors.
    for (int idx=0, proc=0; proc<numProcs; proc++)
    {
      if ( proc == myProc )
      {
        for ( std::vector<int>::iterator it = offProcLocalRows_.begin(); it != offProcLocalRows_.end(); it++, idx++ )
        {
          offProcLocalRowsRecvPIDs_.push_back( allOffProcRowsRecvPIDs[idx] );
        }
        SortContainer2( offProcLocalRowsRecvPIDs_, offProcLocalRows_ );
      }
      else
      {
        for ( int i=0; i<sumNZOffProcRows[proc]; i++, idx++ )
        {
          // Check if this processor owns this row, by looking at the non-overlapped map.
          int lid = appVecPtr_->pmap()->globalToLocalIndex( allOffProcRows[idx] );
          if ( lid != -1 )
          {
            offProcNonlocalRows_.push_back( lid );
            offProcNonlocalRowsSendPIDs_.push_back( allOffProcRowsSendPIDs[idx] );
          }
        }
      }
    }
  }

  // Order nonlocal rows and pids.
  SortContainer2( offProcNonlocalRowsSendPIDs_, offProcNonlocalRows_ );

  // Order and create maps for local nonzero rows
  std::sort( freqNZLocalRows_.begin(), freqNZLocalRows_.end() );
  freqNZLocalRows_.erase( std::unique( freqNZLocalRows_.begin(), freqNZLocalRows_.end() ), freqNZLocalRows_.end() ); 
  for (unsigned idx=0; idx < freqNZLocalRows_.size(); idx++)
  {
    freqNZLocalRowsMap_[freqNZLocalRows_[idx]] = idx;  /* Add lid to map for matvec later */
  }


  // Create an overlap map to be used for matrix-vector products with the frequency-domain matrix entries.
  // NOTE:  This map could be reduced to only those columns for which we have entries to reduce parallel communication.
  int numharms = bVtPtr_->blockCount();
  N_PDS_ParMap * colmap = appdFdxPtr_->getColMap( *builder_.getPDSComm() );

  if ( numProcs > 1 )
  {
    if (hbOsc_)
    {
      std::vector<int> augLIDs( 1 );
      overlapMap_ = Linear::createBlockFreqERFParMap( numharms, *appVecPtr_->pmap(), *colmap, 1, &augLIDs );
    }
    else
    {
      overlapMap_ = Linear::createBlockFreqERFParMap( numharms, *appVecPtr_->pmap(), *colmap );
    }
  }

  // Now collect a global number of nzRows
  int numNZLocalRows = freqNZLocalRows_.size();
  appVecPtr_->pmap()->pdsComm().sumAll( &numNZLocalRows, &numGlobalFreqRows_, 1 );

  // The frequency-domain loading analysis is done, don't perform again.
  freqLoadAnalysisDone_ = true;

}

void HBLoader::createPermFreqBVector( std::vector< std::vector< Util::FreqVecEntry > >& vectorEntries,
                                      Teuchos::RCP<Linear::BlockVector>& blockVector )
{
  int numharms = bVtPtr_->blockCount();
  Linear::Vector freqB( *(hbBuilderPtr_->getSolutionMap()) );

  int myProc = appVecPtr_->pmap()->pdsComm().procID();
  int numProcs = appVecPtr_->pmap()->pdsComm().numProc();

  // Filtering out local ground node entries.
  // In parallel, entries to nonlocal nodes must be migrated to the necessary processor.
  int local_gnd_id = appVecPtr_->omap()->globalToLocalIndex( -1 );

  if ( totalOffProcBVecLIDs_ )
  {
    // Send / receive vector entries
    for ( int p = 0; p < numProcs; ++p )
    {
      // Sending processor is p
      if ( p == myProc && offProcBVecPIDs_.size() )
      {
        for ( int recvp = 0; recvp < numProcs; ++recvp )
        {
          int totalData = 0;
          std::vector<int> gidVecs;
          std::vector<double> valVecs;

          std::vector<int>::iterator it = std::find( offProcBVecPIDs_.begin(), offProcBVecPIDs_.end(), recvp );
          while ( it != offProcBVecPIDs_.end() )
          {
            int idx = std::distance( offProcBVecPIDs_.begin(), it );
            for( int i = 0; i < (numharms + 1)/2; ++i )
            {
              for (unsigned j = 0; j < vectorEntries[i].size(); j++)
              {
                const Util::FreqVecEntry& currEntry = vectorEntries[i][j];

                // Check if this entry is on this processor.
                if ( currEntry.lid == offProcBVecLIDs_[idx] )
                {
                  totalData++;
                  gidVecs.push_back( i );
                  gidVecs.push_back( appVecPtr_->omap()->localToGlobalIndex( currEntry.lid ) );
                  valVecs.push_back( currEntry.val.real() );       
                  valVecs.push_back( currEntry.val.imag() );       
                }
              }
            }

            it++;
            it = std::find( it, offProcBVecPIDs_.end(), recvp );
          }

          // Send vector entries to processor recvp 
          if (totalData)
          {
            appVecPtr_->pmap()->pdsComm().send( &totalData, 1, recvp);
            appVecPtr_->pmap()->pdsComm().send( &gidVecs[0], 2*totalData, recvp);
            appVecPtr_->pmap()->pdsComm().send( &valVecs[0], 2*totalData, recvp);
          }
        }
      }
      // Otherwise wait to receive data from p
      else if ( std::find(offProcBVecSendPIDs_.begin(), offProcBVecSendPIDs_.end(), p ) != offProcBVecSendPIDs_.end() )
      {
        Util::FreqVecEntry tmpEntry;

        int totalData = 0;
        appVecPtr_->pmap()->pdsComm().recv( &totalData, 1, p);
        std::vector<int> gidVecs( 2*totalData );
        std::vector<double> valVecs( 2*totalData );
        appVecPtr_->pmap()->pdsComm().recv( &gidVecs[0], 2*totalData, p);
        appVecPtr_->pmap()->pdsComm().recv( &valVecs[0], 2*totalData, p);

        // Process the vector entries sent from p
        for (int i=0; i<totalData; ++i)
        {
          int lid = appVecPtr_->pmap()->globalToLocalIndex( gidVecs[2*i+1] );
          tmpEntry.lid = appVecPtr_->pmap()->globalToLocalIndex( gidVecs[2*i+1] );
          tmpEntry.val = std::complex<double>( valVecs[2*i], valVecs[2*i+1] );
          vectorEntries[ gidVecs[2*i] ].push_back( tmpEntry );        
        }
      }
    }
  }

  int addFreqB = 0;
  for( int i = 0; i < (numharms + 1)/2; ++i )
  {
    for (unsigned j = 0; j < vectorEntries[i].size(); j++)
    {
      const Util::FreqVecEntry& currEntry = vectorEntries[i][j];
      if (currEntry.lid != local_gnd_id)
      {
        // Check if this entry is on this processor.
        std::vector<int>::iterator it = std::find(offProcBVecLIDs_.begin(),offProcBVecLIDs_.end(),currEntry.lid);
        if ( it  == offProcBVecLIDs_.end() )
        {
          freqB[currEntry.lid*2*numharms + 2*i] += currEntry.val.real();
          freqB[currEntry.lid*2*numharms + 2*i+1] += currEntry.val.imag();
          if (i > 0)
          {
            freqB[currEntry.lid*2*numharms + 2*(numharms-i)] += currEntry.val.real();
            freqB[currEntry.lid*2*numharms + 2*(numharms-i)+1] += -currEntry.val.imag();
          }
          addFreqB = 1;
        }
      }
    }
  }

  int globalFreqB = 0;
  appVecPtr_->pmap()->pdsComm().maxAll( &addFreqB, &globalFreqB, 1 );

  if (globalFreqB)
  {
    blockVector = hbBuilderPtr_->createExpandedRealFormTransposeBlockVector();
    blockVector->update( 1.0, freqB, 0.0 );
  }
}

} // namespace Loader
} // namespace Xyce

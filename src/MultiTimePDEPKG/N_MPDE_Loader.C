//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_DeviceMgr.h>
#include <N_MPDE_Loader.h>
#include <N_MPDE_Builder.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_Manager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>

#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

using Xyce::DEBUG_MPDE;

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::~N_MPDE_Loader
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_MPDE_Loader::~N_MPDE_Loader()
{
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::setFastTimes
// Purpose       : Assign times for fast time scale
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
void N_MPDE_Loader::setFastTimes( const std::vector<double> & times )
{ 
  times_ = times;
  constructPeriodicTimes();
}


void N_MPDE_Loader::setPeriodFlags( const std::vector<bool> & periodicFlags )
{ 
  nonPeriodic_ = periodicFlags;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::constructPeriodicTimes
// Purpose       : Make a copy of times_ vector with padding
//                 at the beginning and end to make the calculation
//                 of derivatives easier around the periodicy condition
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
void N_MPDE_Loader::constructPeriodicTimes()
{
  // we will pad our array of times with 2 times the width, one
  // at the top and one at the bottom of the array
  periodicTimesOffset_ = fastTimeDiscretizer_.Width();
  int timesSize = times_.size();
  period_ = times_[timesSize - 1];
  periodicTimes_.resize( timesSize + 2*periodicTimesOffset_ );
  for( int i=0; i< periodicTimesOffset_; i++ )
  {
    periodicTimes_[i] = times_[i + timesSize - periodicTimesOffset_ - 1] - period_;
  }
  for( int i=periodicTimesOffset_; i< (timesSize + periodicTimesOffset_); i++ )
  {
    periodicTimes_[i] = times_[i - periodicTimesOffset_];
  }
  for( int i=(timesSize+periodicTimesOffset_); i< (timesSize + 2*periodicTimesOffset_); i++ )
  {
    periodicTimes_[i] = period_ - times_[i - (timesSize+periodicTimesOffset_) + 1];
  }
  
  if (DEBUG_MPDE)
  {
    Xyce::dout() << "periodicTimes_ array is" << std::endl;
    for( int i=0; i<(int)periodicTimes_.size(); i++)
    {
      Xyce::dout() << "  periodicTimes_[ " << i << " ] = " << periodicTimes_[i];
      if( (i >= periodicTimesOffset_) && (i < (timesSize + periodicTimesOffset_) ) )
      {
        Xyce::dout() << "\ttimes_[ " << (i-periodicTimesOffset_) << " ] = " << times_[i-periodicTimesOffset_];
      } 
      Xyce::dout() << std::endl;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::loadDAEMatrices( Xyce::Linear::Vector * X,
                                     Xyce::Linear::Vector * S,
                                     Xyce::Linear::Vector * dSdt,
                                     Xyce::Linear::Vector * Store,
                                     Xyce::Linear::Matrix * dQdx,
                                     Xyce::Linear::Matrix * dFdx,
                                     int loadType )
{
  if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_MPDE_Loader::loadDAEMatrices" << std::endl;
    
  }

  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  Xyce::Linear::Vector * appdSdt = appNextStaVecPtr_->cloneVector();

  Xyce::Linear::BlockMatrix & bdQdx = *dynamic_cast<Xyce::Linear::BlockMatrix*>(dQdx);
  Xyce::Linear::BlockMatrix & bdFdx = *dynamic_cast<Xyce::Linear::BlockMatrix*>(dFdx);
  Xyce::Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>(X);

#ifdef Xyce_FLEXIBLE_DAE_LOADS
  Xyce::Linear::BlockVector & bS = *dynamic_cast<Xyce::Linear::BlockVector*>(S);
  Xyce::Linear::BlockVector & bdSdt = *dynamic_cast<Xyce::Linear::BlockVector*>(dSdt);
  Xyce::Linear::BlockVector & bStore = *dynamic_cast<Xyce::Linear::BlockVector*>(Store);
#endif // Xyce_FLEXIBLE_DAE_LOADS
  
  int BlockCount = bX.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Processing diagonal matrix block " << i << " of " << BlockCount-1 << std::endl;
    }

#ifdef Xyce_FLEXIBLE_DAE_LOADS
    //Set Time for fast time scale somewhere
    state_.fastTime = times_[i];
    deviceManager_.setFastTime( times_[i] );

    //Update the sources
    loader_.updateSources();

    *appNextVecPtr_ = bX.block(i);
    *appNextStaVecPtr_ = bS.block(i);
    *appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);

    loader_.loadDAEMatrices( appNextVecPtr_, appNextStaVecPtr_, appdSdt, 
        appNextStoVecPtr_, appdQdxPtr_, appdFdxPtr_);

    bdQdx.block(i,i).add( *appdQdxPtr_ );
    bdFdx.block(i,i).add( *appdFdxPtr_ );
#else
    //For now, the matrices are loaded during the loadDAEVectors method
    //Just copied here
    bdQdx.block(i,i).add( bmdQdxPtr_->block(i,i) );
    bdFdx.block(i,i).add( bmdFdxPtr_->block(i,i) );

#endif // Xyce_FLEXIBLE_DAE_LOADS
  }
  if (warpMPDEPhasePtr_)
  {
    // Add \dot{phi} = omega dQdx contribution here
    int phiGID = warpMPDEPhasePtr_->getPhiGID();
    int phiLID = bX.pmap()->globalToLocalIndex(phiGID);
    if (phiLID >= 0)
    {
      std::vector<int> colIndices;
      std::vector<double> coeffs;
      colIndices.push_back(phiGID);
      coeffs.push_back(1.0);
      bdQdx.replaceAugmentedRow( phiGID, colIndices.size(), &coeffs[0], &colIndices[0] );
    
      // tscoffe/tmei 08/10/05:  Add omega dependency for warped MPDE to every dq/dt2 equation
      // tscoffe 12/12/06:  For WaMPDE, we're augmenting two columns, omega & phi, and
      // the indexing starts at zero, so the index of omega is 0.
      // tscoffe 01/11/07:  Add \dot{phi} = omega dFdx derivative term here:
      (*bOmegadQdt2Ptr_)[phiLID] = -1.0;
    }
    bdFdx.replaceAugmentedColumn(warpMPDEPhasePtr_->getOmegaGID(),*bOmegadQdt2Ptr_);
  }

  // Fast Time scale discretization terms:
  // These are d(dQ/dt1)/dx terms, but go into bdFdx.  
  // For this procedure, need to re-use the app matrix, appdQdx.
  Xyce::Linear::Matrix & dQdxdt = *appdQdxPtr_;

  int Start = fastTimeDiscretizer_.Start();
  int Width = fastTimeDiscretizer_.Width();
  
  const std::vector<double> & Coeffs = fastTimeDiscretizer_.Coeffs();
  
  for( int i = 0; i < BlockCount; ++i )
  {
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Processing off diagonal matrix blocks on row " << i << " of " << BlockCount-1 << std::endl;
    }

    int Loc;
    int indexT1 = i + Start + periodicTimesOffset_;
    int indexT2 = indexT1 + Width - 1;
    double invh2 = 1.0 / (periodicTimes_[indexT2] - periodicTimes_[indexT1]);
    
    for( int j = 0; j < Width; ++j )
    {
      Loc = i + (j+Start);
      
      if( Loc < 0 )
      {
        Loc += BlockCount;
      }
      else if( Loc > (BlockCount-1) )
      {
        Loc -= BlockCount;
      }
      
      dQdxdt.put(0.0);
      dQdxdt.add( bdQdx.block(Loc,Loc) );
      dQdxdt.scale( Coeffs[j]*invh2 );
      bdFdx.block(i,Loc).add( dQdxdt );
    }
  }

  // tscoffe/tmei 08/10/05:  Add omega equation for warped MPDE
  if ( warpMPDEPhasePtr_ ) 
  {
    std::vector<int> colIndices;
    std::vector<double> coeffs;
    int omegaGID = warpMPDEPhasePtr_->getOmegaGID();
    warpMPDEPhasePtr_->getPhaseConditionDerivative(&bX,times_,&colIndices,&coeffs);
    bdFdx.replaceAugmentedRow( omegaGID, colIndices.size(), &coeffs[0], &colIndices[0] );
  }

  // Now that the matrix loading is finished, call fillComplete().
  dQdx->fillComplete();
  dFdx->fillComplete();

  // For BlockMatrix objects, synchronize the global copy of the block matrix.
  bdQdx.assembleGlobalMatrix();
  bdFdx.assembleGlobalMatrix();
 
  if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PRINT_VECTORS))
  {
    Xyce::dout() << "MPDE bX:" << std::endl;
    bX.print(std::cout);
    Xyce::dout() << "MPDE bdQdx:" << std::endl;
    bdQdx.print(std::cout);
    Xyce::dout() << "MPDE bdFdx:" << std::endl;
    bdFdx.print(std::cout);
#ifdef Xyce_FLEXIBLE_DAE_LOADS
    Xyce::dout() << "MPDE bS:" << std::endl;
    bS.print(std::cout);
    Xyce::dout() << "MPDE dSdt:" << std::endl;
    bdSdt.print(std::cout);
    Xyce::dout() << "MPDE bStore:" << std::endl;
    bStore.print(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS
  
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  delete appdSdt;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::allDevicesConverged
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/08/2019
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::allDevicesConverged(Xyce::Parallel::Machine comm)
{
  return allDevicesAllTimePointsConverged_;
}

//-----------------------------------------------------------------------------
// Function      : N_LOA_MPDELoader::updateState
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
bool N_MPDE_Loader::updateState      
 (Xyce::Linear::Vector * nextSolVectorPtr, 
  Xyce::Linear::Vector * currSolVectorPtr,
  Xyce::Linear::Vector * lastSolVectorPtr,
  Xyce::Linear::Vector * nextStaVectorPtr,
  Xyce::Linear::Vector * currStaVectorPtr,
  Xyce::Linear::Vector * lastStaVectorPtr,
  Xyce::Linear::Vector * nextStoVectorPtr,
  Xyce::Linear::Vector * currStoVectorPtr,
  Xyce::Linear::Vector * lastStoVectorPtr,
  int loadType )
{
  bool bsuccess = true;

  // For MPDE case, this needs to be a no-op.

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Loader::loadDAEVectors( Xyce::Linear::Vector * X,
                                    Xyce::Linear::Vector * currX,
                                    Xyce::Linear::Vector * lastX,
                                    Xyce::Linear::Vector * S,
                                    Xyce::Linear::Vector * currS,
                                    Xyce::Linear::Vector * lastS,
                                    Xyce::Linear::Vector * dSdt,
                                    Xyce::Linear::Vector * Store,
                                    Xyce::Linear::Vector * currStore,
                                    Xyce::Linear::Vector * lastStore,
                                    Xyce::Linear::Vector * nextLeadFVectorPtr,
                                    Xyce::Linear::Vector * nextLeadQVectorPtr,
                                    Xyce::Linear::Vector * nextJunctionVVectorPtr,
                                    Xyce::Linear::Vector * Q,
                                    Xyce::Linear::Vector * F,
                                    Xyce::Linear::Vector * B,
                                    Xyce::Linear::Vector * dFdxdVp,
                                    Xyce::Linear::Vector * dQdxdVp,
                                    int loadType )
{
  if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
  {
    Xyce::dout() << std::endl
           << Xyce::section_divider << std::endl
           << "  N_MPDE_Loader::loadDAEVectors" << std::endl
           << "warpMPDE flag = " << (warpMPDEPhasePtr_ ? "true" : "false") << std::endl;
  }

  //Zero out vectors
  appNextVecPtr_->putScalar(0.0);
  appCurrVecPtr_->putScalar(0.0);
  appLastVecPtr_->putScalar(0.0);

  appNextStaVecPtr_->putScalar(0.0);
  appCurrStaVecPtr_->putScalar(0.0);
  appLastStaVecPtr_->putScalar(0.0);
  Xyce::Linear::Vector* appdSdt = appNextStaVecPtr_->cloneVector();
  appNextStoVecPtr_->putScalar(0.0);
  appCurrStoVecPtr_->putScalar(0.0);
  appLastStoVecPtr_->putScalar(0.0);

  Xyce::Linear::Vector * appQ = appNextVecPtr_->cloneVector();
  Xyce::Linear::Vector * appF = appNextVecPtr_->cloneVector();
  Xyce::Linear::Vector * appB = appNextVecPtr_->cloneVector();
  Xyce::Linear::Vector * appdFdxdVp = appNextVecPtr_->cloneVector();
  Xyce::Linear::Vector * appdQdxdVp = appNextVecPtr_->cloneVector();

  // This is a temporary load storage vector.
  Xyce::Linear::Vector * dQdt2 = appNextVecPtr_->cloneVector(); 

  // 12/8/06 tscoffe:   Note:  "b" at beginning of variable name means Xyce::Linear::BlockVector
  Xyce::Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>(X);
  Xyce::Linear::BlockVector & bcurrX = *dynamic_cast<Xyce::Linear::BlockVector*>(currX);
  Xyce::Linear::BlockVector & blastX = *dynamic_cast<Xyce::Linear::BlockVector*>(lastX);
  Xyce::Linear::BlockVector & bS = *dynamic_cast<Xyce::Linear::BlockVector*>(S);
  Xyce::Linear::BlockVector & bcurrS = *dynamic_cast<Xyce::Linear::BlockVector*>(currS);
  Xyce::Linear::BlockVector & blastS = *dynamic_cast<Xyce::Linear::BlockVector*>(lastS);
  Xyce::Linear::BlockVector & bdSdt = *dynamic_cast<Xyce::Linear::BlockVector*>(dSdt);
  Xyce::Linear::BlockVector & bStore = *dynamic_cast<Xyce::Linear::BlockVector*>(Store);
  Xyce::Linear::BlockVector & bcurrStore = *dynamic_cast<Xyce::Linear::BlockVector*>(currStore);
  Xyce::Linear::BlockVector & blastStore = *dynamic_cast<Xyce::Linear::BlockVector*>(lastStore);
 
  Xyce::Linear::BlockVector & bNextLeadF = *dynamic_cast<Xyce::Linear::BlockVector*>(nextLeadFVectorPtr);
  Xyce::Linear::BlockVector & bLeadQ = *dynamic_cast<Xyce::Linear::BlockVector*>(nextLeadQVectorPtr);
  Xyce::Linear::BlockVector & bNextJunctionV = *dynamic_cast<Xyce::Linear::BlockVector*>(nextJunctionVVectorPtr);
  
  Xyce::Linear::BlockVector & bQ = *dynamic_cast<Xyce::Linear::BlockVector*>(Q);
  Xyce::Linear::BlockVector & bF = *dynamic_cast<Xyce::Linear::BlockVector*>(F);
  Xyce::Linear::BlockVector & bB = *dynamic_cast<Xyce::Linear::BlockVector*>(B);

  Xyce::Linear::BlockVector & bdFdxdVp = *dynamic_cast<Xyce::Linear::BlockVector*>(dFdxdVp);
  Xyce::Linear::BlockVector & bdQdxdVp = *dynamic_cast<Xyce::Linear::BlockVector*>(dQdxdVp);

#ifndef Xyce_FLEXIBLE_DAE_LOADS
  bmdQdxPtr_->put(0.0);
  bmdFdxPtr_->put(0.0);
#endif
    
  allDevicesAllTimePointsConverged_ = true;
  int BlockCount = bQ.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Processing vectors for block " << i << " of " << BlockCount-1 << std::endl;
    }

    //Set Time for fast time scale somewhere
    state_.fastTime = times_[i];
    deviceManager_.setFastTime( times_[i] );
    
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Calling updateSources on the appLoader" << std::endl;
    }

    //Update the sources
    loader_.updateSources();  // this is here to handle "fast" sources.

    *appNextVecPtr_ = bX.block(i);
    *appCurrVecPtr_ = bcurrX.block(i);
    *appLastVecPtr_ = blastX.block(i);

    *appNextStaVecPtr_ = bS.block(i);
    *appCurrStaVecPtr_ = bcurrS.block(i);
    *appLastStaVecPtr_ = blastS.block(i);
    *appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);
    *appCurrStoVecPtr_ = bcurrStore.block(i);
    *appLastStoVecPtr_ = blastStore.block(i);
    
    *appNextLeadFVecPtr_  = bNextLeadF.block(i);
    *appLeadQVecPtr_      =  bLeadQ.block(i);
    *appNextJunctionVVecPtr_  =  bNextJunctionV.block(i);
    
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Updating State for block " << i << " of " << BlockCount-1 << std::endl;
    }

    // Note: This updateState call is here (instead of in the 
    // N_MPDE_Loader::updateState function) because it has to be called
    // for the same fast time point.
    loader_.updateState 
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_,
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_ , &*appLastStaVecPtr_ ,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_ , &*appLastStoVecPtr_ );

    bS.block(i) = *appNextStaVecPtr_;
    bcurrS.block(i) = *appCurrStaVecPtr_;
    blastS.block(i) = *appLastStaVecPtr_;
    bStore.block(i) = *appNextStoVecPtr_;
    bcurrStore.block(i) = *appCurrStoVecPtr_;
    blastStore.block(i) = *appLastStoVecPtr_;

    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Calling loadDAEVectors on the appLoader" << std::endl;
    }

    // This has to be done because the app loader does NOT zero these vectors out.
    appQ->putScalar(0.0);
    appF->putScalar(0.0);
    appB->putScalar(0.0);
    appdFdxdVp->putScalar(0.0);
    appdQdxdVp->putScalar(0.0);

    // RLS need to fix Store vector passage with lead current junction voltage equivalent
    loader_.loadDAEVectors
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_, 
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_, &*appLastStaVecPtr_, 
        &*appdSdt, &*appNextStoVecPtr_, &*appCurrStoVecPtr_, &*appLastStoVecPtr_, 
        &*appNextLeadFVecPtr_, &*appLeadQVecPtr_, 
        &*appNextJunctionVVecPtr_, 
        &*appQ, &*appF, &*appB,
        &*appdFdxdVp, &*appdQdxdVp );

    // get the device convergence status
    bool allDevsConv = loader_.allDevicesConverged(appQ->pmap()->pdsComm().comm());
    bool tmpVal = allDevicesAllTimePointsConverged_;
    allDevicesAllTimePointsConverged_ = tmpVal && allDevsConv;

    bQ.block(i) = *appQ;
    bF.block(i) = *appF;

    bB.block(i) = *appB;
    bdFdxdVp.block(i) = *appdFdxdVp;
    bdQdxdVp.block(i) = *appdQdxdVp;

#ifndef Xyce_FLEXIBLE_DAE_LOADS
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Processing matrices for block " << i << " of " << BlockCount-1 << std::endl;
    }

    // This has to be done because the app loader does NOT zero these out.
    appdQdxPtr_->put(0.0);
    appdFdxPtr_->put(0.0);

    loader_.loadDAEMatrices( &*appNextVecPtr_, &*appNextStaVecPtr_, &*appdSdt, &*appNextStoVecPtr_, 
                                    &*appdQdxPtr_, &*appdFdxPtr_);

    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Copying diagonal block into bmdQdx" << std::endl;
    }
    bmdQdxPtr_->block(i,i).add( *appdQdxPtr_ );
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Copying diagonal block into bmdFdx" << std::endl;
    }
    bmdFdxPtr_->block(i,i).add( *appdFdxPtr_ );
#endif
  }
  
  int phiGID = -1;
  int phiLID = -1;
  if (warpMPDEPhasePtr_)
  {
    // Add \dot{phi(t_1)} = omega(t_1) term to Q
    phiGID = warpMPDEPhasePtr_->getPhiGID();
    phiLID = bQ.pmap()->globalToLocalIndex(phiGID);
    if (phiLID >= 0)
    {
      double phiValue = bX[phiLID];

      if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
      {
        Xyce::dout() << "Inserting phi with value = " << phiValue << " into q" << std::endl;
      }

      bQ.setElementByGlobalIndex( phiGID, phiValue );
    }
  }

  //-------------------------------------------------------------
  // Now do the fast time scale time derivative -----------------
  //-------------------------------------------------------------

  int Start = fastTimeDiscretizer_.Start();
  int Width = fastTimeDiscretizer_.Width();
  const std::vector<double> & Coeffs = fastTimeDiscretizer_.Coeffs();
    
  double omega = 1.0;
  int omegaGID = -1;
  if (warpMPDEPhasePtr_)
  {
    omegaGID = warpMPDEPhasePtr_->getOmegaGID();

    // tscoffe/tmei 08/04/05:  zero out temporary vector for use in Jacobian
    bOmegadQdt2Ptr_->putScalar(0.0);
    // tscoffe/tmei 08/11/05:  Convert global indices to local indices:

    int omegaLID = bX.pmap()->globalToLocalIndex(omegaGID);

    double tmpOmega = 0.0; 
    if (omegaLID >= 0)
    {
      tmpOmega = bX[omegaLID];
    }
    bX.pmap()->pdsComm().sumAll( &tmpOmega, &omega, 1 );
  }

  for( int i = 0; i < BlockCount; ++i )
  {
    if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
    {
      Xyce::dout() << "Processing dQdt2 block " << i << " of " << BlockCount-1 << std::endl;
    }
    dQdt2->putScalar(0.0);
    
    int Loc;
    int indexT1 = i + Start + periodicTimesOffset_;
    int indexT2 = indexT1 + Width - 1;
    double invh2 = 1.0 / (periodicTimes_[indexT2] - periodicTimes_[indexT1]);
    
    for( int j = 0; j < Width; ++j )
    {
      Loc = i + (j+Start);
      
      if( Loc < 0 )
      {
        Loc += BlockCount;
      }
      else if( Loc > (BlockCount-1) )
      {
        Loc -= BlockCount;
      }
      dQdt2->update( Coeffs[j]*invh2, bQ.block(Loc), 1.0 );
    }
    
    if (warpMPDEPhasePtr_)
    {
      // 12/8/06 tscoffe:  save dQdt2 term for Jacobian evaluation
      bOmegadQdt2Ptr_->block(i) = *dQdt2; // This is not scaled by omega
    }
    // Update F with omega*dq/dt2 term
    bF.block(i).update( omega, *dQdt2, 1.0 );
  }

  //  tscoffe/tmei 08/04/05:  Add omega equation to end of f:
  if ( warpMPDEPhasePtr_ )
  {
    double phaseValue = warpMPDEPhasePtr_->getPhaseCondition(&bX,times_);
    if ( bF.pmap()->globalToLocalIndex(omegaGID) >= 0 )
    {
      if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
      {
        Xyce::dout() << "Inserting phase condition = " << phaseValue << " into f" << std::endl;
        Xyce::dout() << "Inserting omega = " << omega << " into f" << std::endl;
      }

      bF.setElementByGlobalIndex( omegaGID, phaseValue );
      // Add \dot{phi(t_1)} = omega(t_1) term to F
      bF.setElementByGlobalIndex( phiGID, -omega );
    }
  }

  if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PRINT_VECTORS))
  {
    Xyce::dout() << "MPDE X Vector" << std::endl;
    bX.print(std::cout);
    Xyce::dout() << "MPDE S Vector" << std::endl;
    bS.print(std::cout);
    Xyce::dout() << "MPDE dSdt Vector" << std::endl;
    bdSdt.print(std::cout);
    Xyce::dout() << "MPDE Store Vector" << std::endl;
    bStore.print(std::cout);
    Xyce::dout() << "MPDE Q Vector" << std::endl;
    bQ.print(std::cout);
    Xyce::dout() << "MPDE F Vector" << std::endl;
    bF.print(std::cout);

#ifndef Xyce_FLEXIBLE_DAE_LOADS
    bmdQdxPtr_->assembleGlobalMatrix();
    Xyce::dout() << "MPDE bmdQdx_" << std::endl;
    bmdQdxPtr_->print(std::cout);

    bmdFdxPtr_->assembleGlobalMatrix();
    Xyce::dout() << "MPDE bmdFdx_" << std::endl;
    bmdFdxPtr_->print(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS

    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  delete appQ;
  delete appF;
  delete appB;
  delete appdFdxdVp;
  delete appdQdxdVp;
  delete dQdt2;
  delete appdSdt;
 
  return true;
}


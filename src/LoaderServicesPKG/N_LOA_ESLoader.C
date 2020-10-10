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
// Creator       : Eric Keiter
// Creation Date :
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <numeric>

// ----------   Xyce Includes   ----------

#include <N_LOA_ESLoader.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Builder.h>
#include <N_LAS_ESBuilder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_MultiVector.h>

#include <N_ANP_UQSupport.h>

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

//-----------------------------------------------------------------------------
// Function      : ESLoader::ESLoader
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
ESLoader::ESLoader(
  Device::DeviceMgr &                 device_manager,
  Linear::Builder &                   builder,
  int numSamples,
  Analysis::SweepVector & samplingVector,
  const std::vector<double> & Y,
  bool useExprSamples
  )
  : CktLoader(device_manager, builder),
    deviceManager_(device_manager),
    builder_(builder),
    numSamples_(numSamples),
    samplingVector_(samplingVector),
    Y_(Y),
    allDevicesAllBlocksConverged_(true),
    useExpressionSamples_(useExprSamples)
{
  // Now initialize all the working vectors, size of the original system
  appNextVecPtr_ = rcp(builder_.createVector());
  appCurrVecPtr_ = rcp(builder_.createVector());
  appLastVecPtr_ = rcp(builder_.createVector());

  appQPtr_ = rcp(builder_.createVector());
  appFPtr_ = rcp(builder_.createVector());
  appBPtr_ = rcp(builder_.createVector());
  appdFdxdVpPtr_ = rcp(builder_.createVector());
  appdQdxdVpPtr_ = rcp(builder_.createVector());

  appNextStaVecPtr_ = rcp(builder_.createStateVector());
  appCurrStaVecPtr_ = rcp(builder_.createStateVector());
  appLastStaVecPtr_ = rcp(builder_.createStateVector());

  appdSdtPtr_ = rcp(builder_.createStateVector());

  appdQdxPtr_ = rcp(builder_.createMatrix());
  appdFdxPtr_ = rcp(builder_.createMatrix());

  appNextStoVecPtr_ = rcp(builder_.createStoreVector());
  appCurrStoVecPtr_ = rcp(builder_.createStoreVector());
  
  appNextLeadFVecPtr_ = rcp(builder.createLeadCurrentVector());
  appLeadQVecPtr_     = rcp(builder.createLeadCurrentVector());
  appNextJunctionVVecPtr_ = rcp(builder.createLeadCurrentVector());
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::registerESBuilder
// Purpose       : Registration method for the ES builder
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void ESLoader::registerESBuilder( Teuchos::RCP<Linear::ESBuilder> esBuilderPtr )
{
  esBuilderPtr_ = esBuilderPtr;
  bmdQdxPtr_ = esBuilderPtr_->createBlockMatrix();
  bmdFdxPtr_ = esBuilderPtr_->createBlockMatrix();
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::allDevicesConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/08/2019
//-----------------------------------------------------------------------------
bool ESLoader::allDevicesConverged(Xyce::Parallel::Machine comm)
{
  return allDevicesAllBlocksConverged_ ;
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool ESLoader::loadDAEMatrices( Linear::Vector * X,
                                Linear::Vector * S,
                                Linear::Vector * dSdt,
                                Linear::Vector * Store,
                                Linear::Matrix * dQdx,
                                Linear::Matrix * dFdx,
                                int loadType)
{

  if (DEBUG_ES)
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_LOA ESLoader::loadDAEMatrices" << std::endl;
    
  }

  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  //Xyce::Linear::Vector appdSdt( *appNextStaVecPtr_ );
  //Xyce::Linear::Vector & appdSdt = &*appdSdtPtr_;
  Xyce::Linear::Vector & appdSdt = *appdSdtPtr_;

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

    Xyce::Loader::Loader &loader_ = *(appLoaderPtr_);
    bool reset = false;
    if (useExpressionSamples_)
    {
      reset = Xyce::Analysis::UQ::updateExpressionSamplingTerms(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numSamples_, false);
    }
    else
    {
      reset = Xyce::Analysis::UQ::updateSamplingParams(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numSamples_, false);
    }

    if (DEBUG_ES)
    {
      Xyce::dout() << "Processing diagonal matrix block " << i << " of " << BlockCount-1 << std::endl;
    }

#ifdef Xyce_FLEXIBLE_DAE_LOADS
    // set params!
    //Set Time for fast time scale somewhere
    //state_.fastTime = times_[i];
    //deviceManager_.setFastTime( times_[i] );

    //Update the sources
    //loader_.updateSources();

    *appNextVecPtr_ = bX.block(i);
    *appNextStaVecPtr_ = bS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);

    //loader_.loadDAEMatrices( appNextVecPtr_, appNextStaVecPtr_, &appdSdt, 
    appLoaderPtr_->loadDAEMatrices( appNextVecPtr_, appNextStaVecPtr_, &appdSdt, 
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

  // Now that the matrix loading is finished, call fillComplete().
  dQdx->fillComplete();
  dFdx->fillComplete();

  // For BlockMatrix objects, synchronize the global copy of the block matrix.
  bdQdx.assembleGlobalMatrix();
  bdFdx.assembleGlobalMatrix();
 
  if (DEBUG_ES)
  {
    Xyce::dout() << "ES bX:" << std::endl;
    bX.printPetraObject(std::cout);
    Xyce::dout() << "ES bdQdx:" << std::endl;
    bdQdx.printPetraObject(std::cout);
    Xyce::dout() << "ES bdFdx:" << std::endl;
    bdFdx.printPetraObject(std::cout);
#ifdef Xyce_FLEXIBLE_DAE_LOADS
    Xyce::dout() << "ES bS:" << std::endl;
    bS.printPetraObject(std::cout);
    Xyce::dout() << "ES dSdt:" << std::endl;
    bdSdt.printPetraObject(std::cout);
    Xyce::dout() << "ES bStore:" << std::endl;
    bStore.printPetraObject(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS
  
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool ESLoader::loadDAEVectors( Linear::Vector * X,
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
  if (DEBUG_ES)
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  ESLoader::loadDAEVectors" << std::endl;
  }

  //Zero out vectors
  appNextVecPtr_->putScalar(0.0);
  appCurrVecPtr_->putScalar(0.0);
  appLastVecPtr_->putScalar(0.0);

  appNextStaVecPtr_->putScalar(0.0);
  appCurrStaVecPtr_->putScalar(0.0);
  appLastStaVecPtr_->putScalar(0.0);
  appNextStoVecPtr_->putScalar(0.0);
  appCurrStoVecPtr_->putScalar(0.0);

  appdSdtPtr_->putScalar(0.0);
  Xyce::Linear::Vector & appdSdt = *appdSdtPtr_;

  appQPtr_->putScalar(0.0);
  appFPtr_->putScalar(0.0);
  appBPtr_->putScalar(0.0);
  appdFdxdVpPtr_->putScalar(0.0);
  appdQdxdVpPtr_->putScalar(0.0);

  Xyce::Linear::Vector & appQ = *appQPtr_;
  Xyce::Linear::Vector & appF = *appFPtr_;
  Xyce::Linear::Vector & appB = *appBPtr_;
  Xyce::Linear::Vector & appdFdxdVp = *appdFdxdVpPtr_;
  Xyce::Linear::Vector & appdQdxdVp = *appdQdxdVpPtr_;

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
    
  allDevicesAllBlocksConverged_ = true;
  int BlockCount = bQ.blockCount();
  for( int i = 0; i < BlockCount; ++i )
  {
    Xyce::Loader::Loader &loader_ = *(appLoaderPtr_);
    bool reset = false;
    if (useExpressionSamples_)
    {
      reset = Xyce::Analysis::UQ::updateExpressionSamplingTerms(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numSamples_, false);
    }
    else
    {
      reset = Xyce::Analysis::UQ::updateSamplingParams(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numSamples_, false);
    }

    if (DEBUG_ES)
    {
      Xyce::dout() << "Processing vectors for block " << i << " of " << BlockCount-1 << std::endl;
    }
    
    if (DEBUG_ES)
    {
      Xyce::dout() << "Calling updateSources on the appLoader" << std::endl;
    }

    // note: should these be views rather than copies?
    *appNextVecPtr_ = bX.block(i);
    *appCurrVecPtr_ = bcurrX.block(i);
    *appLastVecPtr_ = blastX.block(i);

    *appNextStaVecPtr_ = bS.block(i);
    *appCurrStaVecPtr_ = bcurrS.block(i);
    *appLastStaVecPtr_ = blastS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bStore.block(i);
    *appCurrStoVecPtr_ = bcurrStore.block(i);
    
    *appNextLeadFVecPtr_  = bNextLeadF.block(i);
    *appLeadQVecPtr_      =  bLeadQ.block(i);
    *appNextJunctionVVecPtr_  =  bNextJunctionV.block(i);
    
    if (DEBUG_ES)
    {
      Xyce::dout() << "Updating State for block " << i << " of " << BlockCount-1 << std::endl;
    }

    // Note: This updateState call is here (instead of in the 
    // N_LOA_ESLoader::updateState function) because it has to be called
    // for the same fast time point.
    appLoaderPtr_->updateState 
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_,
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_ , &*appLastStaVecPtr_ ,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_ );

    bS.block(i) = *appNextStaVecPtr_;
    bcurrS.block(i) = *appCurrStaVecPtr_;
    blastS.block(i) = *appLastStaVecPtr_;
    bStore.block(i) = *appNextStoVecPtr_;
    bcurrStore.block(i) = *appCurrStoVecPtr_;

    if (DEBUG_ES)
    {
      Xyce::dout() << "Calling loadDAEVectors on the appLoader" << std::endl;
    }

    // This has to be done because the app loader does NOT zero these vectors out.
    appQ.putScalar(0.0);
    appF.putScalar(0.0);
    appB.putScalar(0.0);
    appdFdxdVp.putScalar(0.0);
    appdQdxdVp.putScalar(0.0);

    // RLS need to fix Store vector passage with lead current junction voltage equivalent
    appLoaderPtr_->loadDAEVectors
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_, 
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_, &*appLastStaVecPtr_, 
        &appdSdt, &*appNextStoVecPtr_, &*appCurrStoVecPtr_, 
        &*appNextLeadFVecPtr_, &*appLeadQVecPtr_, 
        &*appNextJunctionVVecPtr_, 
        &appQ, &appF, &appB,
        &appdFdxdVp, &appdQdxdVp );

    // get the device convergence status
    bool allDevsConv = appLoaderPtr_->allDevicesConverged(appQ.pmap()->pdsComm().comm());
    bool tmpVal = allDevicesAllBlocksConverged_;
    allDevicesAllBlocksConverged_ = tmpVal && allDevsConv;

    bQ.block(i) = appQ;
    bF.block(i) = appF;

    bB.block(i) = appB;
    bdFdxdVp.block(i) = appdFdxdVp;
    bdQdxdVp.block(i) = appdQdxdVp;

#ifndef Xyce_FLEXIBLE_DAE_LOADS
    if (DEBUG_ES)
    {
      Xyce::dout() << "Processing matrices for block " << i << " of " << BlockCount-1 << std::endl;
    }

    // This has to be done because the app loader does NOT zero these out.
    appdQdxPtr_->put(0.0);
    appdFdxPtr_->put(0.0);

    appLoaderPtr_->loadDAEMatrices( &*appNextVecPtr_, &*appNextStaVecPtr_, &appdSdt, &*appNextStoVecPtr_, 
                                    &*appdQdxPtr_, &*appdFdxPtr_);

    if (DEBUG_ES)
    {
      Xyce::dout() << "Copying diagonal block into bmdQdx" << std::endl;
    }

    bmdQdxPtr_->block(i,i).add( *appdQdxPtr_ );

    if (DEBUG_ES)
    {
      Xyce::dout() << "Copying diagonal block into bmdFdx" << std::endl;
    }
    bmdFdxPtr_->block(i,i).add( *appdFdxPtr_ );
#endif
  }
  
  // Now that the vector loading is finished, synchronize the global copy of the block vector
  bX.assembleGlobalVector();
  bS.assembleGlobalVector();
  bdSdt.assembleGlobalVector();
  bStore.assembleGlobalVector();
  bQ.assembleGlobalVector();
  bF.assembleGlobalVector();
  bdFdxdVp.assembleGlobalVector();
  bdQdxdVp.assembleGlobalVector();
  bcurrS.assembleGlobalVector();
  blastS.assembleGlobalVector();
  bcurrStore.assembleGlobalVector();
  bNextLeadF.assembleGlobalVector();
  bLeadQ.assembleGlobalVector();
  bNextJunctionV.assembleGlobalVector();

  if (DEBUG_ES)
  {
    Xyce::dout() << "ES X Vector" << std::endl;
    bX.printPetraObject(std::cout);
    Xyce::dout() << "ES S Vector" << std::endl;
    bS.printPetraObject(std::cout);
    Xyce::dout() << "ES dSdt Vector" << std::endl;
    bdSdt.printPetraObject(std::cout);
    Xyce::dout() << "ES Store Vector" << std::endl;
    bStore.printPetraObject(std::cout);
    Xyce::dout() << "ES Q Vector" << std::endl;
    bQ.printPetraObject(std::cout);
    Xyce::dout() << "ES F Vector" << std::endl;
    bF.printPetraObject(std::cout);

#ifndef Xyce_FLEXIBLE_DAE_LOADS
    bmdQdxPtr_->assembleGlobalMatrix();
    Xyce::dout() << "ES bmdQdx_" << std::endl;
    bmdQdxPtr_->printPetraObject(std::cout);

    bmdFdxPtr_->assembleGlobalMatrix();
    Xyce::dout() << "ES bmdFdx_" << std::endl;
    bmdFdxPtr_->printPetraObject(std::cout);
#endif // Xyce_FLEXIBLE_DAE_LOADS

    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::loadDeviceErrorWeightMask
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 11/1/2014
//-----------------------------------------------------------------------------
bool ESLoader::loadDeviceErrorWeightMask(Linear::Vector * deviceMask) const
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::getVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015
//---------------------------------------------------------------------------
bool ESLoader::getVoltageLimiterStatus()
{
  return appLoaderPtr_->getVoltageLimiterStatus();
}

//-----------------------------------------------------------------------------
// Function      : ESLoader::setVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015 
//---------------------------------------------------------------------------
void ESLoader::setVoltageLimiterStatus(bool voltageLimterStatus)
{
  return appLoaderPtr_->setVoltageLimiterStatus(voltageLimterStatus);
}

} // namespace Loader
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#ifdef Xyce_STOKHOS_ENABLE

// ---------- Standard Includes ----------
#include <iostream>
#include <numeric>

// ----------   Xyce Includes   ----------
#include <N_LOA_PCELoader.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Builder.h>
#include <N_LAS_PCEBuilder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_MultiVector.h>

#include <N_LAS_Solver.h>
#include <N_LAS_TranSolverFactory.h>
#include <N_LAS_Problem.h>

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
// Function      : PCELoader::PCELoader
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
PCELoader::PCELoader(
  Device::DeviceMgr &                 device_manager,
  Linear::Builder &                   builder,
  int numQuadPoints,
  int numBlockRows,
  Analysis::SweepVector & samplingVector,
  const std::vector<double> & Y,
  const Xyce::IO::CmdParse & cp,
  int voltLimAlg,
  bool useExprSamples
  )
  : CktLoader(device_manager, builder),
    deviceManager_(device_manager),
    builder_(builder),
    numQuadPoints_(numQuadPoints),
    numBlockRows_(numBlockRows),
    samplingVector_(samplingVector),
    Y_(Y),
    commandLine_(cp),
    voltLimAlgorithm_(voltLimAlg),
    allDevicesAllQuadPointsConverged_(true),
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

  app_dV_voltlim_Ptr_ = rcp(builder_.createVector());

  appNextStaVecPtr_ = rcp(builder_.createStateVector());
  appCurrStaVecPtr_ = rcp(builder_.createStateVector());
  appLastStaVecPtr_ = rcp(builder_.createStateVector());

  appdSdtPtr_ = rcp(builder_.createStateVector());

  appdQdxPtr_ = rcp(builder_.createMatrix());
  appdFdxPtr_ = rcp(builder_.createMatrix());

  appNextStoVecPtr_ = rcp(builder_.createStoreVector());
  appCurrStoVecPtr_ = rcp(builder_.createStoreVector());
  appLastStoVecPtr_ = rcp(builder_.createStoreVector());
  
  appNextLeadFVecPtr_ = rcp(builder.createLeadCurrentVector());
  appLeadQVecPtr_     = rcp(builder.createLeadCurrentVector());
  appNextJunctionVVecPtr_ = rcp(builder.createLeadCurrentVector());

  // Create Linear::Problem
  lasProblemPtr_ = rcp(Xyce::Linear::createProblem(appdFdxPtr_.get(), app_dV_voltlim_Ptr_.get(), appdFdxdVpPtr_.get()));
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::~PCELoader
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
PCELoader::~PCELoader()
{
  delete bmdQdx_ptr_;
  delete bmdFdx_ptr_;
  delete bmdQdx_quad_ptr_;
  delete bmdFdx_quad_ptr_;
  delete b_dV_voltlim_coef_Ptr_;
  delete b_dV_voltlim_quad_Ptr_;

  delete bQ_quad_ptr_;
  delete bF_quad_ptr_;
  delete bB_quad_ptr_;
  delete bdFdxdVp_quad_ptr_;
  delete bdQdxdVp_quad_ptr_;

  delete bXNext_quad_ptr_;
  delete bXCurr_quad_ptr_;
  delete bXLast_quad_ptr_;
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::registerPCEBuilder
// Purpose       : Registration method for the PCE builder
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
void PCELoader::registerPCEBuilder( Teuchos::RCP<Linear::PCEBuilder> pceBuilderPtr )
{
  pceBuilderPtr_ = pceBuilderPtr;

  bmdQdx_ptr_ = dynamic_cast<Linear::BlockMatrix *>(pceBuilderPtr_->createMatrix());
  bmdFdx_ptr_ = dynamic_cast<Linear::BlockMatrix *>(pceBuilderPtr_->createMatrix());

  bmdQdx_quad_ptr_ = pceBuilderPtr_->createQuadMatrix();
  bmdFdx_quad_ptr_ = pceBuilderPtr_->createQuadMatrix();

  bQ_quad_ptr_ = pceBuilderPtr_->createQuadVector();
  bF_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 
  bB_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 
  bdFdxdVp_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 
  bdQdxdVp_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 

  bXNext_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 
  bXCurr_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 
  bXLast_quad_ptr_ = pceBuilderPtr_->createQuadVector(); 

  b_dV_voltlim_quad_Ptr_ = pceBuilderPtr_->createQuadVector();
  b_dV_voltlim_coef_Ptr_ = dynamic_cast<Linear::BlockVector *>(pceBuilderPtr_->createVector());
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::registerSolverFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/04/2019
//-----------------------------------------------------------------------------
void PCELoader::registerSolverFactory (Xyce::Linear::SolverFactory *tmpLasSolverPtr) 
{ 
  lasSolverFactoryPtr_ = tmpLasSolverPtr; 
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::allDevicesConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/08/2019
//-----------------------------------------------------------------------------
bool PCELoader::allDevicesConverged(Xyce::Parallel::Machine comm)
{
  return allDevicesAllQuadPointsConverged_ ;
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
bool PCELoader::loadDAEMatrices( Linear::Vector * X,
                                 Linear::Vector * S,
                                 Linear::Vector * dSdt,
                                 Linear::Vector * Store,
                                 Linear::Matrix * dQdx,
                                 Linear::Matrix * dFdx,
                                 int loadType)
{

  if (DEBUG_PCE)
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  N_LOA PCELoader::loadDAEMatrices" << std::endl;
    
  }

  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  Xyce::Linear::BlockMatrix & bdQdx = *dynamic_cast<Xyce::Linear::BlockMatrix*>(dQdx);
  Xyce::Linear::BlockMatrix & bdFdx = *dynamic_cast<Xyce::Linear::BlockMatrix*>(dFdx);
  Xyce::Linear::BlockVector & bnextX = *dynamic_cast<Xyce::Linear::BlockVector*>(X);

  int basisSize = basis_->size();
  for( int i = 0; i < basisSize; ++i )
  {
    for( int j = 0; j < basisSize; ++j )
    {
      //The matrices are loaded during the loadDAEVectors method, and are copied here
      bdQdx.block(i,j).add( bmdQdx_ptr_->block(i,j) );
      bdFdx.block(i,j).add( bmdFdx_ptr_->block(i,j) );
    }
  }

  // Now that the matrix loading is finished, call fillComplete().
  dQdx->fillComplete();
  dFdx->fillComplete();

  // For BlockMatrix objects, synchronize the global copy of the block matrix.
  bdQdx.assembleGlobalMatrix();
  bdFdx.assembleGlobalMatrix();
 
  if (DEBUG_PCE)
  {
    Xyce::dout() << "PCE bnextX:" << std::endl;
    bnextX.print(std::cout);
    Xyce::dout() << "PCE bdQdx:" << std::endl;
    bdQdx.print(std::cout);
    Xyce::dout() << "PCE bdFdx:" << std::endl;
    bdFdx.print(std::cout);
  
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : returnDenseMatrixEntry
// Purpose       : This function takes a PCE expansion of a single Jacobian 
//                 entry and returns the fully Stochastic Jacobian entry that 
//                 results from performing a Kronecker product between that 
//                 expansion and the triple product tensor (Cjik).
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 8/29/2019
//-----------------------------------------------------------------------------
void returnDenseMatrixEntry(
    Sacado::PCE::OrthogPoly<double,Stokhos::StandardStorage<int,double> > const &inval,
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double> > & denseEntry,
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > const &Cijk,
    const Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & basis
    )
{
  Stokhos::OrthogPolyApprox<int, double> val= inval.getOrthogPolyApprox();
  typedef Stokhos::Sparse3Tensor<int, double> Cijk_type;
  int pb = val.size();
  const double* cv = val.coeff();

  const Teuchos::Array<double>& norms = basis->norm_squared(); 

  denseEntry->putScalar(0.0);
  typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
  typename Cijk_type::k_iterator k_end = Cijk->k_end();
  if (pb < Cijk->num_k())
    k_end = Cijk->find_k(pb);
  double cijk;
  int i,j,k;
  for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) 
  {
    k = index(k_it);
    for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); j_it != Cijk->j_end(k_it); ++j_it) 
    {
       j = index(j_it);
       for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it); i_it != Cijk->i_end(j_it); ++i_it) 
       {
         i = index(i_it);
         cijk = value(i_it);
         (*denseEntry)(i,j) += (cijk/norms[i])*cv[k];
       }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::allocateVoltageLimitingSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
void PCELoader::allocateVoltageLimitingSolver ()
{
  if ( Teuchos::is_null(lasSolverRCPtr_) )
  {
    if (lasSolverFactoryPtr_)
    {
      lasSolverRCPtr_ = Teuchos::rcp( lasSolverFactoryPtr_->create( saved_lsOB_, *lasProblemPtr_ , commandLine_) );
    }
    else
    {
      Linear::TranSolverFactory lasFactory;
      lasSolverRCPtr_ = Teuchos::rcp( lasFactory.create( saved_lsOB_, *lasProblemPtr_ , commandLine_) );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
bool PCELoader::loadDAEVectors( Linear::Vector * X,
                                Linear::Vector * currX,
                                Linear::Vector * lastX,
                                Linear::Vector * S,
                                Linear::Vector * currS,
                                Linear::Vector * lastS,
                                Linear::Vector * dSdt,
                                Linear::Vector * Store,
                                Linear::Vector * currStore,
                                Linear::Vector * lastStore,
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
  if (DEBUG_PCE)
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  PCELoader::loadDAEVectors" << std::endl;
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
  appLastStoVecPtr_->putScalar(0.0);

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

  // Note:  "b" at beginning of variable name means Xyce::Linear::BlockVector
  Xyce::Linear::BlockVector & bnextX = *dynamic_cast<Xyce::Linear::BlockVector*>(X);
  Xyce::Linear::BlockVector & bcurrX = *dynamic_cast<Xyce::Linear::BlockVector*>(currX);
  Xyce::Linear::BlockVector & blastX = *dynamic_cast<Xyce::Linear::BlockVector*>(lastX);
  Xyce::Linear::BlockVector & bnextS = *dynamic_cast<Xyce::Linear::BlockVector*>(S);
  Xyce::Linear::BlockVector & bcurrS = *dynamic_cast<Xyce::Linear::BlockVector*>(currS);
  Xyce::Linear::BlockVector & blastS = *dynamic_cast<Xyce::Linear::BlockVector*>(lastS);
  Xyce::Linear::BlockVector & bdSdt = *dynamic_cast<Xyce::Linear::BlockVector*>(dSdt);
  Xyce::Linear::BlockVector & bnextStore = *dynamic_cast<Xyce::Linear::BlockVector*>(Store);
  Xyce::Linear::BlockVector & bcurrStore = *dynamic_cast<Xyce::Linear::BlockVector*>(currStore);
  Xyce::Linear::BlockVector & blastStore = *dynamic_cast<Xyce::Linear::BlockVector*>(lastStore);
 
  Xyce::Linear::BlockVector & bNextLeadF = *dynamic_cast<Xyce::Linear::BlockVector*>(nextLeadFVectorPtr);
  Xyce::Linear::BlockVector & bLeadQ = *dynamic_cast<Xyce::Linear::BlockVector*>(nextLeadQVectorPtr);
  Xyce::Linear::BlockVector & bNextJunctionV = *dynamic_cast<Xyce::Linear::BlockVector*>(nextJunctionVVectorPtr);
  
  Xyce::Linear::BlockVector & bQ = *dynamic_cast<Xyce::Linear::BlockVector*>(Q);
  Xyce::Linear::BlockVector & bF = *dynamic_cast<Xyce::Linear::BlockVector*>(F);
  Xyce::Linear::BlockVector & bB = *dynamic_cast<Xyce::Linear::BlockVector*>(B);

  Xyce::Linear::BlockVector & bDV = *(b_dV_voltlim_coef_Ptr_);

  Xyce::Linear::BlockVector & bdFdxdVp = *dynamic_cast<Xyce::Linear::BlockVector*>(dFdxdVp);
  Xyce::Linear::BlockVector & bdQdxdVp = *dynamic_cast<Xyce::Linear::BlockVector*>(dQdxdVp);

  bmdQdx_ptr_->put(0.0);
  bmdFdx_ptr_->put(0.0);

  bmdQdx_quad_ptr_->put(0.0);
  bmdFdx_quad_ptr_->put(0.0);

  // ERK.  Convert the solution vector bnextX from coefficients to quadrature variable values.
  //
  // Use the UQ helper function for evaluating the PCE expansion.
  // Evaluate the PCE approximation at the sample points, then use this as input 
  // to the device models load function calls
  //
  // ERK.  Note:  to be quicker, I am converting all 3 solution vectors: next, curr and last.  
  // But of course, once a transient is up and running, the curr and last should have been 
  // converted already.  So I am doing some redundant work.
  {
  std::vector< Stokhos::OrthogPolyApprox<int,double> > pceVec(3);
  int solutionSize = bnextX.block(0).localLength();  // get local length.  SERIAL ONLY HERE!!!  sigh, fix later.

  for (int isol=0;isol<solutionSize;++isol)
  {
    pceVec[0].reset(basis_);
    pceVec[1].reset(basis_);
    pceVec[2].reset(basis_);
    int basisSize = basis_->size();
    for (int icoef=0;icoef<basisSize;++icoef)
    {
      *appNextVecPtr_ = bnextX.block(icoef);
      *appCurrVecPtr_ = bcurrX.block(icoef);
      *appLastVecPtr_ = blastX.block(icoef);

      pceVec[0][icoef] = (*appNextVecPtr_)[isol];
      pceVec[1][icoef] = (*appCurrVecPtr_)[isol];
      pceVec[2][icoef] = (*appLastVecPtr_)[isol];
    }

    std::vector < std::vector<double> > xSamples(pceVec.size(), std::vector<double>(numQuadPoints_,0.0) );
    Xyce::Analysis::UQ::evaluateApproximationPCE(samplingVector_, Y_, numQuadPoints_, pceVec, xSamples);
    for (int iquad=0;iquad<numQuadPoints_;++iquad)
    {
      (bXNext_quad_ptr_->block(iquad))[isol] = xSamples[0][iquad];
      (bXCurr_quad_ptr_->block(iquad))[isol] = xSamples[1][iquad];
      (bXLast_quad_ptr_->block(iquad))[isol] = xSamples[2][iquad];
    }
  }
  }

  // This loop is over the number of quadrature points
  bool applyLimit=true;// ERK note.  This is (for now) unconditionally true. 
  if (voltLimAlgorithm_==1)
  {
    b_dV_voltlim_quad_Ptr_->putScalar(0.0);
  }

#if 0
  {
  std::cout << "Block Count of bnextS = " << bnextS.blockCount() << std::endl;
  std::cout << "Block Size  of bnextS = " << bnextS.blockSize () << std::endl;
  std::cout << "numQuadPoints         = " << numQuadPoints_ <<std::endl;
  int solutionSize = bnextX.block(0).localLength();  // get local length.  SERIAL ONLY HERE!!!  sigh, fix later.
  std::cout << "solutionSize          = " << solutionSize <<std::endl;
  std::cout << "Basis size            = " << basis_->size() <<std::endl;
  }
#endif

  allDevicesAllQuadPointsConverged_ = true;

  for( int i = 0; i < numQuadPoints_; ++i )
  {
    Xyce::Loader::Loader &loader_ = *(appLoaderPtr_);
    bool reset = false;

    if (useExpressionSamples_)
    {
      reset = Xyce::Analysis::UQ::updateExpressionSamplingTerms2(loader_, i, samplingVector_, Y_, numQuadPoints_, false);
    }
    else
    {
      reset = Xyce::Analysis::UQ::updateSamplingParams(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numQuadPoints_, false);
    }

    if (DEBUG_PCE)
    {
      Xyce::dout() << "Processing vectors for block " << i << " of " << numQuadPoints_-1 << std::endl;
    }

    // pull the various vectors out of the block objects
    //
    // these 3 solution vectors were transformed from PCE coefs to sample (quad) values
    *appNextVecPtr_ = bXNext_quad_ptr_->block(i);
    *appCurrVecPtr_ = bXCurr_quad_ptr_->block(i);
    *appLastVecPtr_ = bXLast_quad_ptr_->block(i);

    // all the other vectors (state, store, lead current, etc) don't need to go thru a transformation between
    // the PCE representation and sample values.  The reason is that they aren't part of the PCE linear system
    // solve.  So they can always have the values produced by the device models, and not go back and forth.
    //
    // These vectors, below all need to have their number of blocks = # quad points
    *appNextStaVecPtr_ = bnextS.block(i);
    *appCurrStaVecPtr_ = bcurrS.block(i);
    *appLastStaVecPtr_ = blastS.block(i);
    appdSdt = bdSdt.block(i);
    *appNextStoVecPtr_ = bnextStore.block(i);
    *appCurrStoVecPtr_ = bcurrStore.block(i);
    *appLastStoVecPtr_ = blastStore.block(i);
    
    *appNextLeadFVecPtr_  = bNextLeadF.block(i);
    *appLeadQVecPtr_      =  bLeadQ.block(i);
    *appNextJunctionVVecPtr_  =  bNextJunctionV.block(i);
    
    if (DEBUG_PCE)
    {
      Xyce::dout() << "Updating State for block " << i << " of " << numQuadPoints_-1 << std::endl;
    }

    // Note: updateState call is here (instead of N_LOA_PCELoader::updateState function) 
    // because it has to be called for the same fast time point. Same for the matrix load, below
    appLoaderPtr_->updateState 
      ( &*appNextVecPtr_, &*appCurrVecPtr_, &*appLastVecPtr_,
        &*appNextStaVecPtr_, &*appCurrStaVecPtr_ , &*appLastStaVecPtr_ ,
        &*appNextStoVecPtr_, &*appCurrStoVecPtr_, &*appLastStoVecPtr_ );

    bnextS.block(i) = *appNextStaVecPtr_;
    bcurrS.block(i) = *appCurrStaVecPtr_;
    blastS.block(i) = *appLastStaVecPtr_;
    bnextStore.block(i) = *appNextStoVecPtr_;
    bcurrStore.block(i) = *appCurrStoVecPtr_;
    blastStore.block(i) = *appLastStoVecPtr_;

    if (DEBUG_PCE)
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
        &appdSdt, &*appNextStoVecPtr_, &*appCurrStoVecPtr_, &*appLastStoVecPtr_, 
        &*appNextLeadFVecPtr_, &*appLeadQVecPtr_, 
        &*appNextJunctionVVecPtr_, 
        &appQ, &appF, &appB,
        &appdFdxdVp, &appdQdxdVp );

    // get the device convergence status
    bool allDevsConv = appLoaderPtr_->allDevicesConverged(appQ.pmap()->pdsComm().comm());
    bool tmpVal = allDevicesAllQuadPointsConverged_;
    allDevicesAllQuadPointsConverged_ = tmpVal && allDevsConv;

    // all of these bQ, bF, bB, bdFdxdVp, bdQdxdVp need to be in temporary structures that have num blocks = # quad points
    bQ_quad_ptr_->block(i) = appQ;
    bF_quad_ptr_->block(i) = appF;
    bB_quad_ptr_->block(i) = appB;
    bdFdxdVp_quad_ptr_->block(i) = appdFdxdVp;
    bdQdxdVp_quad_ptr_->block(i) = appdQdxdVp;

    if (DEBUG_PCE)
    {
      Xyce::dout() << "Processing matrices for block " << i << " of " << numQuadPoints_-1 << std::endl;
    }

    // This has to be done because the app loader does NOT zero these out.
    appdQdxPtr_->put(0.0);
    appdFdxPtr_->put(0.0);

    appLoaderPtr_->loadDAEMatrices( &*appNextVecPtr_, &*appNextStaVecPtr_, &appdSdt, &*appNextStoVecPtr_, 
                                    &*appdQdxPtr_, &*appdFdxPtr_);

    bmdQdx_quad_ptr_->block(i,i).add( *appdQdxPtr_ );
    bmdFdx_quad_ptr_->block(i,i).add( *appdFdxPtr_ );

    // solve volt lim problem, if necessary.  
    // check max norm. if tiny don't bother
    double maxNormFlimiter=0.0; 
    appdFdxdVp.infNorm(&maxNormFlimiter);
    double maxNormQlimiter=0.0; 
    appdQdxdVp.infNorm(&maxNormQlimiter);

    if (appLoaderPtr_->getVoltageLimiterStatus() && !allDevsConv)
    {
      if (voltLimAlgorithm_==1)
      {
        // get a solution for dV for this quad point
        allocateVoltageLimitingSolver (); // only allocates if null
        bool reuseFactors_ = false;
        lasSolverRCPtr_->solve(reuseFactors_);
        b_dV_voltlim_quad_Ptr_->block(i) = *app_dV_voltlim_Ptr_;
      }
      else if (voltLimAlgorithm_==2)
      {
        // do nothing.  No linear solve necessary.  Doing projections of volt lim vectors instead
      }
      else
      { // do nothing
      }

      applyLimit=true; // this is unconditionally true for now.  I made a mistake with it
    }
#if 0
    else
    {
      // this is probably redundant
      b_dV_voltlim_quad_Ptr_->block(i).putScalar(0.0);
    }
#endif
  }
  
  // Now that the vector loading is finished, synchronize the global copy of the relevant block vectors
  // For the solution-sized vectors (x,f,q,b), use the PCE block vectors
  // For the state, store and lead current-sized vectors, use the quad versions
  
  // matrices
  bmdQdx_quad_ptr_->assembleGlobalMatrix();
  bmdFdx_quad_ptr_->assembleGlobalMatrix();
  bmdQdx_quad_ptr_->fillComplete();
  bmdFdx_quad_ptr_->fillComplete();

#if 0
  std::cout << "Printing bF_quad_ptr_ (quadrature points):" <<std::endl;
  bF_quad_ptr_->print(std::cout);
  double maxNormFquad=0.0; 
  bF_quad_ptr_->infNorm(&maxNormFquad);
  std::cout << "Max norm of bF_quad_ptr_ = " << maxNormFquad <<std::endl;
#endif

  {
  // obtain the PCE coefficients of both f and q, and volt lim if necessary
  // ERK.  Fix this
  // the std::vectors f and q need to contain all the quadrature points for a given variable
  
  int solutionSize = bnextX.block(0).localLength();  // get local length.  SERIAL ONLY HERE!!!  sigh, fix later.
  for (int isol=0;isol<solutionSize;++isol)
  {
    std::vector<double> f(numQuadPoints_,0.0);
    std::vector<double> q(numQuadPoints_,0.0);
    std::vector<double> b(numQuadPoints_,0.0);
    std::vector<double> dv, dFdxdVp, dQdxdVp;

    if (voltLimAlgorithm_==1)
    {
      dv.resize(numQuadPoints_,0.0);
    }
    else if (voltLimAlgorithm_==2)
    {
      dFdxdVp.resize(numQuadPoints_,0.0);
      dQdxdVp.resize(numQuadPoints_,0.0);
    }
    else
    {
    }

    for (int iquad=0;iquad<numQuadPoints_;++iquad)
    {
      f[iquad] = (bF_quad_ptr_->block(iquad))[isol];
      q[iquad] = (bQ_quad_ptr_->block(iquad))[isol];
      b[iquad] = (bB_quad_ptr_->block(iquad))[isol];
      if (applyLimit)
      {
        if (voltLimAlgorithm_==1)
        {
          dv[iquad] = (b_dV_voltlim_quad_Ptr_->block(iquad))[isol];
        }
        else if (voltLimAlgorithm_==2)
        {
          dFdxdVp[iquad] = (bdFdxdVp_quad_ptr_->block(iquad))[isol];
          dQdxdVp[iquad] = (bdQdxdVp_quad_ptr_->block(iquad))[isol];
        }
        else
        {
          // do nothing
        }
      }
    }

    pceF.init(0.0); pceF.reset(expnMethod_);
    pceQ.init(0.0); pceQ.reset(expnMethod_);
    pceB.init(0.0); pceB.reset(expnMethod_);
    Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, f, pceF);
    Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, q, pceQ);
    Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, b, pceB);
    if (applyLimit)
    {
      if (voltLimAlgorithm_==1)
      {
        pceDV.init(0.0); pceDV.reset(expnMethod_); 
        Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, dv, pceDV);
      }
      else if (voltLimAlgorithm_==2)
      {
        pce_dFdxdVp.init(0.0); pce_dFdxdVp.reset(expnMethod_); 
        Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, dFdxdVp, pce_dFdxdVp);

        pce_dQdxdVp.init(0.0); pce_dQdxdVp.reset(expnMethod_); 
        Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, dQdxdVp, pce_dQdxdVp);
      }
      else 
      { // do nothing
      }
    }

    int basisSize = basis_->size();
    for (int icoef=0;icoef<basisSize;++icoef)
    {
      (bF.block(icoef))[isol] = pceF.coeff(icoef);
      (bQ.block(icoef))[isol] = pceQ.coeff(icoef);
      (bB.block(icoef))[isol] = pceB.coeff(icoef);
      if (applyLimit)
      {
        if (voltLimAlgorithm_==1)
        {
          (bDV.block(icoef))[isol] = pceDV.coeff(icoef);
        }
        else if (voltLimAlgorithm_==2)
        {
          (bdFdxdVp.block(icoef))[isol] = pce_dFdxdVp.coeff(icoef);
          (bdQdxdVp.block(icoef))[isol] = pce_dQdxdVp.coeff(icoef); 
        }
        else 
        { // do nothing
        }
      }
    }
  }

#if 0
  std::cout << "Printing bF (PCE coefs):" <<std::endl;
  bF.print(std::cout);
  double maxNormF=0.0; 
  bF.infNorm(&maxNormF);
  std::cout << "Max norm of bF = " << maxNormF <<std::endl;
#endif
  }

  {
  // solve for the PCE coefficients of both dfdx and dqdx
  int solutionSize = bnextX.block(0).localLength();  // get local length.  SERIAL ONLY HERE!!!  sigh, fix later.

  Xyce::Linear::Matrix & subMatRef = bmdFdx_quad_ptr_->block(0,0); // use this to get the structure
  int numLocalRowsRef = subMatRef.getLocalNumRows(); // num ckt vars = n_

  int basisSize = basis_->size();
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double> > denseEntryF = 
    Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,double>( basisSize, basisSize));
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,double> > denseEntryQ = 
    Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,double>( basisSize, basisSize));

  for (int irow=0;irow<numLocalRowsRef;++irow)
  {
    std::vector<double> dfdx(numQuadPoints_,0.0);
    std::vector<double> dqdx(numQuadPoints_,0.0);

    int length; double * coeffs; int * colIndices;
    subMatRef.getLocalRowView(irow, length, coeffs, colIndices);
    for (int icol=0;icol<length;++icol)
    {
      for (int iquad=0;iquad<numQuadPoints_;++iquad)
      {
        Xyce::Linear::Matrix & subMatF = bmdFdx_quad_ptr_->block(iquad,iquad); 
        Xyce::Linear::Matrix & subMatQ = bmdQdx_quad_ptr_->block(iquad,iquad); 
        double fval = subMatF[irow][icol];
        double qval = subMatQ[irow][icol];
        dfdx[iquad] = fval;
        dqdx[iquad] = qval;
#if 0
        std::cout << "dfdx["<<iquad<<"] = " << dfdx[iquad] << std::endl;
#endif
      }

      pceF.init(0.0); pceF.reset(expnMethod_);
      pceQ.init(0.0); pceQ.reset(expnMethod_);
      Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, dfdx, pceF);
      Xyce::Analysis::UQ::solveProjectionPCE(basis_, quadMethod_, dqdx, pceQ);

      returnDenseMatrixEntry( pceF, denseEntryF, Cijk_,basis_);
      returnDenseMatrixEntry( pceQ, denseEntryQ, Cijk_,basis_);
#if 0
      {
      std::cout.setf(std::ios::scientific);
      std::cout.setf(std::ios::right);
      std::cout.width(16);
      std::cout << "--------------------------------------------------------------" <<std::endl;
      std::cout << "PCE expansion of dfdx for row="<<irow<<" col="<<icol<<std::endl;
      pceF.print(std::cout);
      std::cout << "dFdx: dense block for row="<<irow<<" col="<<icol<<std::endl;
      denseEntryF->print(std::cout);
      std::cout << "--------------------------------------------------------------" <<std::endl;

#if 0
      std::cout << "--------------------------------------------------------------" <<std::endl;
      std::cout << "PCE expansion of dqdx for row="<<irow<<" col="<<icol<<std::endl;
      pceQ.print(std::cout);
      std::cout << "dQdx: dense block for row="<<irow<<" col="<<icol<<std::endl;
      denseEntryQ->print(std::cout);
      std::cout << "--------------------------------------------------------------" <<std::endl;
#endif
      }
#endif

      for (int icoefRow=0;icoefRow<basisSize;++icoefRow)
      {
        for (int icoefCol=0;icoefCol<basisSize;++icoefCol)
        {
          Xyce::Linear::Matrix & subMatF = bmdFdx_ptr_->block(icoefRow,icoefCol);
          Xyce::Linear::Matrix & subMatQ = bmdQdx_ptr_->block(icoefRow,icoefCol);
          double fval = (*denseEntryF)(icoefRow,icoefCol);
          double qval = (*denseEntryQ)(icoefRow,icoefCol);
          subMatF[irow][icol] += fval;
          subMatQ[irow][icol] += qval;
        }
      }
    }
  }

  bmdFdx_ptr_->assembleGlobalMatrix();
  bmdQdx_ptr_->assembleGlobalMatrix();
  bmdFdx_ptr_->fillComplete();
  bmdQdx_ptr_->fillComplete();

#if 0
  std::cout << "--------------------------------------------------------------" <<std::endl;
  std::cout << "Full Jacobian" << std::endl;
  bmdFdx_ptr_->print(std::cout);
  std::cout << "--------------------------------------------------------------" <<std::endl;
#endif
  }

  if (applyLimit)
  {
    if (voltLimAlgorithm_==1)
    {
      // perform matvecs to convert bDV to dFdxdVp and dQdxdVp
      bool Transpose = false;
      dFdxdVp->putScalar(0.0);
      bmdFdx_ptr_->matvec( Transpose , bDV, *dFdxdVp );

      dQdxdVp->putScalar(0.0);
      bmdQdx_ptr_->matvec( Transpose , bDV, *dQdxdVp );
    }
  }

  if (DEBUG_PCE)
  {
    Xyce::dout() << "PCE X Vector" << std::endl;
    bnextX.print(std::cout);
    Xyce::dout() << "PCE S Vector" << std::endl;
    bnextS.print(std::cout);
    Xyce::dout() << "PCE dSdt Vector" << std::endl;
    bdSdt.print(std::cout);
    Xyce::dout() << "PCE Store Vector" << std::endl;
    bnextStore.print(std::cout);
    Xyce::dout() << "PCE Q Vector" << std::endl;
    bQ.print(std::cout);
    Xyce::dout() << "PCE F Vector" << std::endl;
    bF.print(std::cout);

    bmdQdx_ptr_->assembleGlobalMatrix();
    Xyce::dout() << "PCE bmdQdx_" << std::endl;
    bmdQdx_ptr_->print(std::cout);

    bmdFdx_ptr_->assembleGlobalMatrix();
    Xyce::dout() << "PCE bmdFdx_" << std::endl;
    bmdFdx_ptr_->print(std::cout);

    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::loadDeviceErrorWeightMask
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
bool PCELoader::loadDeviceErrorWeightMask(Linear::Vector * deviceMask) const
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::getVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//---------------------------------------------------------------------------
bool PCELoader::getVoltageLimiterStatus()
{
  return appLoaderPtr_->getVoltageLimiterStatus();
}

//-----------------------------------------------------------------------------
// Function      : PCELoader::setVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//---------------------------------------------------------------------------
void PCELoader::setVoltageLimiterStatus(bool voltageLimterStatus)
{
  return appLoaderPtr_->setVoltageLimiterStatus(voltageLimterStatus);
}

} // namespace Loader
} // namespace Xyce

#endif // stokhos

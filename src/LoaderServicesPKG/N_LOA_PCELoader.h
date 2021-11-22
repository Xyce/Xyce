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
//
// Purpose        : PCE Specific Loader
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 07/27/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_PCELoader_H
#define Xyce_LOA_PCELoader_H

#ifdef Xyce_STOKHOS_ENABLE
#include <vector>

#include <Teuchos_RCP.hpp>

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>

#include <N_LAS_TranSolverFactory.h>

#include <N_LOA_CktLoader.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_UTL_OptionBlock.h>

#ifdef Xyce_STOKHOS_ENABLE
#include <Stokhos_Sacado.hpp>
#include <Sacado_No_Kokkos.hpp>
#include <Stokhos_Sparse3TensorUtilities.hpp>
#endif

// ---------- Forward declarations --------

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Class         : PCELoader
// Purpose       : PCE specific CktLoader interface
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 07/27/2019
//-----------------------------------------------------------------------------
class PCELoader : public CktLoader
{
public:
  PCELoader(
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    int numQuadPoints, 
    int numBlockRows, 
    Analysis::SweepVector & samplingVector,
    const std::vector<double> & Y,
    const Xyce::IO::CmdParse & commandLine,
    int voltLimAlg,
    bool useExprSamples
    ); 

  ~PCELoader();

  // Get convergence info from devices
  bool allDevicesConverged(Xyce::Parallel::Machine comm);

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( Linear::Vector * X,
                        Linear::Vector * S,
                        Linear::Vector * dSdt,
                        Linear::Vector * Store,
                        Linear::Matrix * dQdx,
                        Linear::Matrix * dFdx,
                        int loadType = Xyce::Device::ALL );

  // Method is called to load the mask to be used in calculating error norms.
  bool loadDeviceErrorWeightMask(Linear::Vector * deviceMask) const;

  // Method which is called to load the new-DAE vectors
  bool loadDAEVectors( Linear::Vector * X,
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
                       int loadType = Xyce::Device::ALL );

  bool updateState(    Linear::Vector * nextSolVectorPtr,
                       Linear::Vector * currSolVectorPtr,
                       Linear::Vector * lastSolVectorPtr,
                       Linear::Vector * nextStaVectorPtr,
                       Linear::Vector * currStaVectorPtr,
                       Linear::Vector * lastStaVectorPtr,
                       Linear::Vector * nextStoVectorPtr,
                       Linear::Vector * currStoVectorPtr,
                       Linear::Vector * lastStoVectorPtr,
                       int loadType = Xyce::Device::ALL
                       )
  { return true; }

    // Virtual method which initializes the nonlinear problem.
  virtual bool initializeProblem( Linear::Vector * nextSolVectorPtr,
                          Linear::Vector * currSolVectorPtr,
                          Linear::Vector * lastSolVectorPtr,
                          Linear::Vector * nextStaVectorPtr,
                          Linear::Vector * currStaVectorPtr,
                          Linear::Vector * lastStaVectorPtr,
                          Linear::Vector * StateDerivVectorPtr,
                          Linear::Vector * nextStoVectorPtr,
                          Linear::Vector * currStoVectorPtr,
                          Linear::Vector * lastStoVectorPtr,
                          Linear::Vector * QVectorPtr,
                          Linear::Vector * FVectorPtr,
                          Linear::Vector * BVectorPtr,
                          Linear::Vector * dFdxdVpVectorPtr,
                          Linear::Vector * dQdxdVpVectorPtr) const
  {
    return false;
  }

  // Get the voltage limiter flag:
  bool getLimiterFlag () { return PCELoader::appLoaderPtr_->getLimiterFlag (); }

  // voltage limiter solver
  void allocateVoltageLimitingSolver ();

  // Registration method for the device packaage
  void registerAppLoader( Teuchos::RCP<Loader> appLoaderPtr )
  { appLoaderPtr_ = appLoaderPtr; }

  void registerPCEBuilder(Teuchos::RCP<Linear::PCEBuilder> pceBuilderPtr);

  void registerPCEbasis (Teuchos::RCP<const Stokhos::ProductBasis<int,double> > & tmpBasis) { basis_ = tmpBasis; }

  void registerPCEquadMethod ( Teuchos::RCP<const Stokhos::Quadrature<int,double> > & tmpQuadMethod) { quadMethod_ = tmpQuadMethod; }

  void registerPCEexpnMethod ( Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > & tmpExpnMethod) { expnMethod_ = tmpExpnMethod; }

  void registerPCEtripleProductTensor ( Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > & tmpCijk) { Cijk_ = tmpCijk; }

  void registerSolverFactory (Xyce::Linear::SolverFactory *tmpLasSolverPtr);

  //void registerLinearSystem (Linear::System * linear_system_ptr) { lasSysPtr_ = linear_system_ptr; }

  void setLinSolOptions(const Util::OptionBlock & OB) { saved_lsOB_ = OB; }

  virtual bool analyticSensitivitiesAvailable (std::string & name) { return false; }
  virtual void getAnalyticSensitivities(
      std::string & name, 
      std::vector<double> & dfdpVec, 
      std::vector<double> & dqdpVec,
      std::vector<double> & dbdpVec,
      std::vector<int> & FindicesVec,
      std::vector<int> & QindicesVec,
      std::vector<int> & BindicesVec) const
  {}

  virtual bool setParam (std::string & name, double val, bool overrideOriginal=false) { return false; }
  virtual double getParamAndReduce(Parallel::Machine comm, const std::string & name) const { return 0.0; }

  // voltage limiter toggle functions
  virtual bool getVoltageLimiterStatus();
  virtual void setVoltageLimiterStatus(bool voltageLimterStatus);

private:
 
  // Base Application loader
  Teuchos::RCP<Loader>          appLoaderPtr_;          ///< Actually a CktLoader
  Device::DeviceMgr &           deviceManager_;

  // Application Linear Objects
  Teuchos::RCP<Linear::Vector> appNextVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrVecPtr_;
  Teuchos::RCP<Linear::Vector> appLastVecPtr_;

  Teuchos::RCP<Linear::Vector> appNextStaVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrStaVecPtr_;
  Teuchos::RCP<Linear::Vector> appLastStaVecPtr_;

  Teuchos::RCP<Linear::Matrix> appdQdxPtr_;
  Teuchos::RCP<Linear::Matrix> appdFdxPtr_;

  // Time domain vectors for loading.  
  Teuchos::RCP<Linear::Vector> appNextStoVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrStoVecPtr_;
  Teuchos::RCP<Linear::Vector> appLastStoVecPtr_;
  
  Teuchos::RCP<Linear::Vector> appNextLeadFVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrLeadFVecPtr_;
  Teuchos::RCP<Linear::Vector> appLeadQVecPtr_;
  Teuchos::RCP<Linear::Vector> appNextJunctionVVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrJunctionVVecPtr_;

  Teuchos::RCP<Linear::Vector> appdSdtPtr_;

  Teuchos::RCP<Linear::Vector> appQPtr_;
  Teuchos::RCP<Linear::Vector> appFPtr_;
  Teuchos::RCP<Linear::Vector> appBPtr_;
  Teuchos::RCP<Linear::Vector> appdFdxdVpPtr_;
  Teuchos::RCP<Linear::Vector> appdQdxdVpPtr_;

  Teuchos::RCP<Linear::Vector> app_dV_voltlim_Ptr_;

  // PCE Builder 
  Teuchos::RCP<Linear::PCEBuilder> pceBuilderPtr_;

  // App Builder
  Linear::Builder &             builder_;

  // Tmp storage block matrices 
  Linear::BlockMatrix* bmdQdx_ptr_;
  Linear::BlockMatrix* bmdFdx_ptr_;
  Linear::BlockMatrix* bmdQdx_quad_ptr_;
  Linear::BlockMatrix* bmdFdx_quad_ptr_;

  // Tmp storage block vectors 
  Linear::BlockVector* bQ_quad_ptr_;
  Linear::BlockVector* bF_quad_ptr_;
  Linear::BlockVector* bB_quad_ptr_;
  Linear::BlockVector* bdFdxdVp_quad_ptr_;
  Linear::BlockVector* bdQdxdVp_quad_ptr_;

  Linear::BlockVector* bXNext_quad_ptr_;
  Linear::BlockVector* bXCurr_quad_ptr_;
  Linear::BlockVector* bXLast_quad_ptr_;

  Linear::BlockVector* b_dV_voltlim_quad_Ptr_;
  Linear::BlockVector* b_dV_voltlim_coef_Ptr_;

  // PCE stuff:
  int numQuadPoints_;
  int numBlockRows_; // this is the size of the PCE expansion, not number of quad points
  Analysis::SweepVector & samplingVector_;
  const std::vector<double> & Y_; 

#ifdef Xyce_STOKHOS_ENABLE
  // many of these objects are copied from the N_ANP_PCE.h header; 
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis_;

  // Quadrature method
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quadMethod_;

  // Triple product tensor
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk_;

  // Expansion method
  Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > expnMethod_;

  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pceF;
  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pceQ;
  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pceB;
  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pceDV;

  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_dFdxdVp;
  Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_dQdxdVp;

#endif

  // these are for solving voltage limiting.  
  // They are the size of the original non-block circuit linear problem
  Xyce::Linear::SolverFactory * lasSolverFactoryPtr_;
  Teuchos::RCP<Xyce::Linear::Problem> lasProblemPtr_;
  Xyce::Linear::Matrix* jacobianMatrixPtr_;
  Xyce::Linear::Vector* rhsVectorPtr_;
  Xyce::Linear::Vector* NewtonVectorPtr_;
  Teuchos::RCP<Xyce::Linear::Solver> lasSolverRCPtr_;
  Util::OptionBlock  saved_lsOB_;
  const Xyce::IO::CmdParse & commandLine_;
  int voltLimAlgorithm_;
  bool allDevicesAllQuadPointsConverged_;
  bool useExpressionSamples_;
};

} // namespace Loader
} // namespace Xyce

#endif // Stokhos

#endif // Xyce_LOA_PCELoader_H

//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : ES Specific Loader
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 06/01/2018
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_ESLoader_H
#define Xyce_LOA_ESLoader_H

#include <vector>

#include <Teuchos_RCP.hpp>

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>

#include <N_LOA_CktLoader.h>
#include <N_UTL_AssemblyTypes.h>

// ---------- Forward declarations --------

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Class         : ESLoader
// Purpose       : ES specific CktLoader interface
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 06/01/2018
//-----------------------------------------------------------------------------
class ESLoader : public CktLoader
{
public:
  ESLoader(
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    int numSamples, 
    Analysis::SweepVector & samplingVector,
    const std::vector<double> & Y); 

  ~ESLoader()
  {}

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
                       int loadType = Xyce::Device::ALL
                       ) {return true; }

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
                          Linear::Vector * QVectorPtr,
                          Linear::Vector * FVectorPtr,
                          Linear::Vector * BVectorPtr,
                          Linear::Vector * dFdxdVpVectorPtr,
                          Linear::Vector * dQdxdVpVectorPtr) const
  {
    return false;
  }

  // Get the voltage limiter flag:
  bool getLimiterFlag () { return ESLoader::appLoaderPtr_->getLimiterFlag (); }

  // Registration method for the device packaage
  void registerAppLoader( Teuchos::RCP<Loader> appLoaderPtr )
  { appLoaderPtr_ = appLoaderPtr; }

  void registerESBuilder(Teuchos::RCP<Linear::ESBuilder> esBuilderPtr);

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

  bool getNumSamples() { return numSamples_; }
  bool setNumSamples(int numS) { numSamples_ = numS; }

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

  Teuchos::RCP<Linear::FilteredMatrix> linAppdQdxPtr_;
  std::vector<Teuchos::RCP<Linear::FilteredMatrix> > vecNLAppdQdxPtr_;
  Teuchos::RCP<Linear::FilteredMatrix> linAppdFdxPtr_;
  std::vector<Teuchos::RCP<Linear::FilteredMatrix> > vecNLAppdFdxPtr_;

  std::vector<int> linNZRows_, nonlinQNZRows_, nonlinFNZRows_;

  // Time domain vectors for loading.  
  Teuchos::RCP<Linear::Vector> appNextStoVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrStoVecPtr_;
  
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

  // ES Builder:
  Teuchos::RCP<Linear::ESBuilder> esBuilderPtr_;

  // App Builder:
  Linear::Builder &             builder_;

  // stuff copied from MPDE loader:
  // Tmp storage block matrices 
  Teuchos::RCP<Xyce::Linear::BlockMatrix> bmdQdxPtr_;
  Teuchos::RCP<Xyce::Linear::BlockMatrix> bmdFdxPtr_;

// sampling stuff:
  int numSamples_;
  Analysis::SweepVector & samplingVector_;
  const std::vector<double> & Y_; 
};

} // namespace Loader
} // namespace Xyce

#endif // Xyce_LOA_ESLoader_H

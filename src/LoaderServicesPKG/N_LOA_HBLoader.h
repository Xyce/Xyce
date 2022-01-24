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

//-----------------------------------------------------------------------------
//
// Purpose        : HB Specific Loader
//
// Special Notes  :
//
// Creator        : Todd Coffey, Ting Mei
//
// Creation Date  : 07/28/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_HBLoader_H
#define Xyce_LOA_HBLoader_H

#include <vector>

#include <Teuchos_RCP.hpp>

#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_DFTInterfaceDecl.hpp>
#include <N_LOA_CktLoader.h>
#include <N_UTL_AssemblyTypes.h>

// ---------- Forward declarations --------

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Class         : HBLoader
// Purpose       : HB specific CktLoader interface
// Special Notes :
// Creator       : Todd Coffey, Ting Mei, Rich Schiek
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
class HBLoader : public CktLoader
{
public:
  HBLoader(
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,

    const int refID,
    const bool hbOsc);

  ~HBLoader()
  {}

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( Linear::Vector * X,
                        Linear::Vector * S,
                        Linear::Vector * dSdt,
                        Linear::Vector * Store,
                        Linear::Matrix * dQdx,
                        Linear::Matrix * dFdx,
                        int loadType = Xyce::Device::ALL );

  // Method for matrix-free 
  bool applyDAEMatrices( Linear::Vector * X,
                         Linear::Vector * S,
                         Linear::Vector * dSdt,
                         Linear::Vector * Store,
                         const Linear::Vector & V,
                         Linear::Vector * dQdxV,
                         Linear::Vector * dFdxV );

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
                       );

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
  bool getLimiterFlag () { return HBLoader::appLoaderPtr_->getLimiterFlag (); }

  // Get the stored time-domain Jacobians from the HB loader.
  Teuchos::RCP<Linear::FilteredMatrix>& getStoreLindQdx() { return linAppdQdxPtr_; }
  Teuchos::RCP<Linear::FilteredMatrix>& getStoreLindFdx() { return linAppdFdxPtr_; }
  std::vector<Teuchos::RCP<Linear::FilteredMatrix> >& getStoreNLdQdx() { return vecNLAppdQdxPtr_; }
  std::vector<Teuchos::RCP<Linear::FilteredMatrix> >& getStoreNLdFdx() { return vecNLAppdFdxPtr_; }

  // Get the stored frequency-domain Jacobian entries from the HB loader.
  int getNumFreqGlobalRows()
  {
    return numGlobalFreqRows_;
  }

  const std::vector< std::vector< Util::FreqMatEntry > >& getFreqDFDXMatrix()
  {
    return freqDFDXMatrix_;
  }

  const std::vector<int>& getFreqNZLocalRows() 
  {
    return freqNZLocalRows_;
  }

  const std::map<int,int>& getFreqNZLocalRowsMap()
  {
    return freqNZLocalRowsMap_;
  }

  bool applyLinearMatrices( const Linear::Vector & Vf,
                            Linear::BlockVector & permlindQdxV,
                            Linear::BlockVector & permlindFdxV );

  // New functions for LOA_HBLoader
  // Assumption:  xt is in the block format:
  // xt.block(i) = { x_n(t_i) }_{n=0..N} N = number of solution components, t_i = time point i
  // P takes xt's block format and converts it to: 
  // (P*xt).block(n) = { x_n(t_i) }_{i=0..T} T = number of time points in block vector.

   
  void setFastTimes( const std::vector<double> & times );

  void setHBFreqs( const std::vector<double> & freqs); 

//  Teuchos::RCP<Linear::BlockVector> & getStoreVecFreqPtr()  { return bStoreVecFreqPtr_;} 
  Teuchos::RCP<Linear::BlockVector> & getLeadCurrentVecFreqPtr()  { return bLeadCurrentVecFreqPtr_;} 

  // Is the problem being applied as a matrix free operator.
  void setMatrixFreeFlag(bool matrixFreeFlag)
  {
    matrixFreeFlag_ = matrixFreeFlag;
  }
 
  void setLoadTimeBFlag(bool loadTimeB )
  {
    loadTimeB_ = loadTimeB;
  }

  // xf = D*P*xt, xf has the same block format as (P*xt), computation limited to input IDs.
  void permutedFFT(const Linear::BlockVector & xt, Linear::BlockVector * xf, std::vector<int>* lids = 0 ); 
  // xf = D*P*xt, xf has the same block format as (P*xt), computation limited to input IDs.
  // NOTE: This method computes 1 norm of input waveform and then performs FFT if nonzero.
  void permutedFFT2(const Linear::BlockVector & xt, Linear::BlockVector * xf); 
  // xt = P^{-1}D^{-1}*xf
  void permutedIFT(const Linear::BlockVector & xf, Linear::BlockVector * xt, int numTimePts_= 0); 

  // Registration method for the device packaage
  void registerAppLoader( Teuchos::RCP<Loader> appLoaderPtr )
  { appLoaderPtr_ = appLoaderPtr; }

  void registerHBBuilder(Teuchos::RCP<Linear::HBBuilder> hbBuilderPtr);

  void registerDFTInterface( const Teuchos::RCP<N_UTL_DFTInterfaceDecl<std::vector<double> > >& dftInterface )
  { dftInterface_ = dftInterface; }

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

  // Helper functions for frequency-domain loading
  void compNZRowsAndCommPIDs( const std::vector< Util::FreqVecEntry >& vectorEntries,
                              const std::vector< Util::FreqVecEntry >& bVecEntries );

  void consolidateMatrixEntries( const std::vector<int>& nzRows,
                                 const std::vector< Util::FreqMatEntry >& matrixEntries,
                                 std::vector< Util::FreqMatEntry >& consolidatedEntries,
                                 bool overlapIDs = true );

  void sendReceiveMatrixEntries( const std::vector< Util::FreqMatEntry >& sendMatrixEntries,
                                 std::vector< Util::FreqMatEntry >& recvMatrixEntries );

  void createPermFreqBVector( std::vector< std::vector< Util::FreqVecEntry > >& vectorEntries,
                              Teuchos::RCP<Linear::BlockVector>& blockVector );

  //Fast Time Scale Points
  std::vector<double> times_;
  std::vector<double> freqs_ ;
  int periodicTimesOffset_;
  std::vector<double> periodicTimes_;
  double period_;
 
  const int refID_;
  const bool hbOsc_;

  bool loadTimeB_;

  // Matrix free flag, operator is being applied not loaded
  bool matrixFreeFlag_;

  // Base Application loader
  Teuchos::RCP<Loader>          appLoaderPtr_;          ///< Actually a CktLoader
  Device::DeviceMgr &           deviceManager_;

  // Application Linear Objects
  Teuchos::RCP<Linear::Vector> appVecPtr_;
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

  // Frequency domain loading objects.
  bool freqLoadAnalysisDone_;
  int numGlobalFreqRows_, totalNZOffProcRows_, totalOffProcBVecLIDs_; 
  std::vector<int> freqNZLocalRows_;
  std::map<int,int> freqNZLocalRowsMap_;
  std::vector< std::vector< Util::FreqVecEntry > > freqBVector_;
  std::vector< std::vector< Util::FreqMatEntry > > freqDFDXMatrix_;
  std::vector<int> offProcBVecLIDs_, offProcBVecSendLIDs_, offProcBVecSendPIDs_, offProcBVecPIDs_;
  std::vector<int> offProcLocalRows_, offProcLocalRowsRecvPIDs_;
  std::vector<int> offProcNonlocalRows_, offProcNonlocalRowsSendPIDs_;
  Teuchos::RCP<Parallel::ParMap> overlapMap_;
  Teuchos::RCP<Linear::BlockVector> permFreqBVector_;

  // Time domain vectors for loading.  
  Teuchos::RCP<Linear::Vector> appNextStoVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrStoVecPtr_;
  Teuchos::RCP<Linear::Vector> appLastStoVecPtr_;
  
  Teuchos::RCP<Linear::Vector> appNextLeadFVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrLeadFVecPtr_;
  Teuchos::RCP<Linear::Vector> appLeadQVecPtr_;
  Teuchos::RCP<Linear::Vector> appNextJunctionVVecPtr_;
  Teuchos::RCP<Linear::Vector> appCurrJunctionVVecPtr_;

  // HB Builder:  (needed to convert AztecOO created Linear::Vectors into Linear::BlockVectors
  Teuchos::RCP<Linear::BlockVector> bQPtr_;
  Teuchos::RCP<Linear::HBBuilder> hbBuilderPtr_;

  // App Builder:  (needed to load time domain vectors and matrices)
  Linear::Builder &             builder_;

  Teuchos::RCP<Linear::BlockVector> bXtPtr_;
  Teuchos::RCP<Linear::BlockVector> bVtPtr_;

//  Teuchos::RCP<Linear::BlockVector> bStoreVecFreqPtr_;
  Teuchos::RCP<Linear::BlockVector> bLeadCurrentVecFreqPtr_; 
  Teuchos::RCP<Linear::BlockVector> bLeadCurrentQVecFreqPtr_; 

  // DFT Interface
  Teuchos::RCP<N_UTL_DFTInterfaceDecl<std::vector<double> > > dftInterface_;
};

} // namespace Loader
} // namespace Xyce

#endif // Xyce_LOA_HBLoader_H

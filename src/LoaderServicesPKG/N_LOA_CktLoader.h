//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Filename       : $RCSfile: N_LOA_CktLoader.h,v $
//
// Purpose        : This file contains class definitions for the loader
//                  services package.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.131 $
//
// Revision Date  : $Date: 2016/03/04 00:34:49 $
//
// Current Owner  : $Author: hkthorn $
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_CktLoader_H
#define Xyce_LOA_CktLoader_H

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_LOA_Loader.h>

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Class         : CktLoader
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
class CktLoader : public Loader
{
public:
  CktLoader(
    Device::DeviceMgr & device_manager,
    Linear::Builder &   builder);

  virtual ~CktLoader();

  // Method which is called to load the new-DAE contributions to the Jacobian matrix.
  bool loadDAEMatrices(
    Linear::Vector *    tmpSolVectorPtr,
    Linear::Vector *    tmpStaVectorPtr,
    Linear::Vector *    tmpStaDerivVectorPtr,
    Linear::Vector *    tmpStoVectorPtr,
    Linear::Matrix *    tmpdQdxMatrixPtr,
    Linear::Matrix *    tmpdFdxMatrixPtr,
    int loadType = Xyce::Device::ALL);

  // Method which is called to load the new-DAE vectors, which contribute
  // to the residual (RHS) vector.
  bool loadDAEVectors(
    Linear::Vector * nextSolVectorPtr,
    Linear::Vector * currSolVectorPtr,
    Linear::Vector * lastSolVectorPtr,
    Linear::Vector * nextStaVectorPtr,
    Linear::Vector * currStaVectorPtr,
    Linear::Vector * lastStaVectorPtr,
    Linear::Vector * StaDerivVectorPtr,
    Linear::Vector * nextStoVectorPtr,
    Linear::Vector * currStoVectorPtr,
    Linear::Vector * lastStoVectorPtr,
    Linear::Vector * nextLeadFVectorPtr,
    Linear::Vector * nextLeadQVectorPtr,
    Linear::Vector * nextJunctionVVectorPtr,
    Linear::Vector * QVectorPtr,
    Linear::Vector * FVectorPtr,
    Linear::Vector * BVectorPtr,
    Linear::Vector * dFdxdVpVectorPtr,
    Linear::Vector * dQdxdVpVectorPtr,
    int loadType = Xyce::Device::ALL);

  /// @brief initializes devices' frequency domain intermediate variables prior to vector/matrix loads
  ///
  /// @param[in] frequency   Frequency of input source
  /// @param[in] freqSolVec  Solution vector at this frequency
  ///
  bool updateFDIntermediateVars(
    double frequency,
    std::complex<double>* freqSolVec);

  /// \brief Obtains DAE contributions to the Jacobian matrix for frequency-based analyses
  /// \param[in]  frequency         The frequency of the input source
  /// \param[in]  freqSolVec        The solution vector in the frequency domain
  /// \param[out] dFdxEntries       Vector of entries into the Jacobian
  bool loadFreqDAEMatrices(
    double frequency,
    std::complex<double>* freqSolVec,
    std::vector<Util::FreqMatEntry>& dFdxEntries);

  /// \brief Obtains contributions to the `F` and `B` vectors for frequency-based analyses
  /// \param[in]  frequency         The frequency of the input source
  /// \param[in]  freqSolVec        The solution vector in the frequency domain
  /// \param[out] FVecEntries        Vector of entries into the `F` vector
  /// \param[out] BVecEntries        Vector of entries into the `B` vector
  bool loadFreqDAEVectors(
    double frequency,
    std::complex<double>* freqSolVec,
    std::vector<Util::FreqVecEntry>& FVecEntries,
    std::vector<Util::FreqVecEntry>& BVecEntries);

  bool loadFreqBVectorsforSources( 
    double frequency,
//    std::complex<double>* freqSolVec,
    std::vector<Util::FreqVecEntry>& BVecEntries);

  // Method is called to load the mask to be used in calculating error norms.
  bool loadDeviceErrorWeightMask(Linear::Vector * deviceMask) const;

  // Initializes the nonlinear problem.
  bool initializeProblem(
    Linear::Vector * nextSolVectorPtr,
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
    Linear::Vector * dQdxdVpVectorPtr) const;

  bool updateState(
    Linear::Vector * nextSolVectorPtr,
    Linear::Vector * currSolVectorPtr,
    Linear::Vector * lastSolVectorPtr,
    Linear::Vector * nextStaVectorPtr,
    Linear::Vector * currStaVectorPtr,
    Linear::Vector * lastStaVectorPtr,
    Linear::Vector * nextStoVectorPtr,
    Linear::Vector * currStoVectorPtr,
    Linear::Vector * lastStoVectorPtr,
    int loadType = Xyce::Device::ALL);

  bool loadBVectorsforAC (Linear::Vector * bVecRealPtr,
                          Linear::Vector * bVecImagPtr);

  bool loadBVectorsforSources ();

  int getNumNoiseSources();

  int getNumNoiseDevices();

  void setupNoiseSources(std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec);

  void getNoiseSources(std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec);

  bool getBMatrixEntries(std::vector<int>& bMatEntriesVec, std::vector<int>& portVec, std::vector<double> * Z0sVec  = NULL);

  // Function for setting the initial guess.
  bool setInitialGuess(Linear::Vector * solVectorPtr);

  // Function for setting a single parameter value.
  bool setParam(std::string & name, double val, bool overrideOriginal=false);

  virtual bool setParamRandomExpressionTerms2(
      const std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
      bool overrideOriginal = false);

  // Function for getting a single parameter value.
  virtual double getParamAndReduce(Xyce::Parallel::Machine comm, const std::string & name) const;

  virtual void getRandomParams(std::vector<Xyce::Analysis::SweepParam> & SamplingParams, Parallel::Communicator & parallel_comm);

  // Method which is called to update the sources.
  bool updateSources();

  // Get the voltage limiter flag:
  bool getLimiterFlag ();

  // Gets the double DC Operating Point flag - used for PDE devices.
  bool isPDESystem() const;
  bool outputPlotFiles() const;
  bool finishOutput() const;

  // two-level newton functions:
  int  enablePDEContinuation();
  bool disablePDEContinuation ();

  void getNumInterfaceNodes (std::vector<int> & numINodes);
  bool loadCouplingRHS(int iSubProblem, int iCouple, Linear::Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iCouple, const Linear::Vector * dxdvPtr);
  // Gets the time integration required breakpoint times (in a vector).
  void resetBreakPoints();// needed for 2-level

  bool getBreakPoints(
      std::vector< Util::BreakPoint > & breakPointTimes,
      std::vector< Util::BreakPoint > & pauseBreakPointTimes
      ) const;

  // Accessor which returns the maximum time step size (in seconds).
  double getMaxTimeStepSize();

  // Get convergence info from devices
  virtual bool allDevicesConverged(Xyce::Parallel::Machine comm);

  // Get convergence info from inner-solves
  bool innerDevicesConverged(Xyce::Parallel::Machine comm);

  void stepSuccess(Xyce::Analysis::TwoLevelMode analysis);
  void stepFailure(Xyce::Analysis::TwoLevelMode analysis);

  void acceptStep();

  virtual bool getInitialQnorm (std::vector<TimeIntg::TwoLevelError> & tleVec );

  virtual bool getInnerLoopErrorSums (std::vector<TimeIntg::TwoLevelError> & tleVec) const;

  bool updateStateArrays ();
  bool startTimeStep(
    bool                          beginIntegrationFlag,
    double                        nextTimeStep,
    double                        nextTime,
    int                           currentOrder);
  void setExternalSolverState(bool external_initJctFlag);

    // Function for determining if an analytic sensitivity (df/dp or dqdp) is available.
  virtual bool analyticSensitivitiesAvailable(const std::string & name);
  virtual bool numericalSensitivitiesAvailable (const std::string & name);

  virtual void getAnalyticSensitivities(
    std::string & name, 
    std::vector<double> & dfdpVec, 
    std::vector<double> & dqdpVec,
    std::vector<double> & dbdpVec,
    std::vector<int> & FindicesVec,
    std::vector<int> & QindicesVec,
    std::vector<int> & BindicesVec) const;

  virtual void getNumericalSensitivities(
      std::string &             name, 
      std::vector<double> &     dfdpVec, 
      std::vector<double> &     dqdpVec,
      std::vector<double> &     dbdpVec,
      std::vector<int> &        FindicesVec,
      std::vector<int> &        QindicesVec,
      std::vector<int> &        BindicesVec) const;

  bool analyticBVecSensAvailable(const std::string & name);
  bool numericalBVecSensAvailable(const std::string & name);

  bool analyticMatrixSensitivitiesAvailable(const std::string & name);
  bool numericalMatrixSensitivitiesAvailable(const std::string & name);

  void getAnalyticBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const;

  void getNumericalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const;

  void getAnalyticMatrixSensitivities(
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const;

  void getNumericalMatrixSensitivities(
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const;

  // voltage limiter toggle functions
  virtual bool getVoltageLimiterStatus();
  virtual void setVoltageLimiterStatus(bool voltageLimterStatus);

#if 0
  virtual void updateDependentParams ();
#endif
  virtual void resetScaledParams();

public:
  Device::DeviceMgr &   deviceManager_;         ///< Device manager
  Linear::Builder &     builder_;               ///< Matrix and vector builder

  // Detect if this is the first step out of a DCOP solve.
  bool dcopState_;

  // Pointers to the linear portion of the Jacobian matrix, if loading is separated.
  Linear::Matrix *      lindQdxMatrixPtr_;
  Linear::Matrix *      lindFdxMatrixPtr_;
  Linear::FilteredMatrix * filtered_lindQdxMatrixPtr_;
  Linear::FilteredMatrix * filtered_lindFdxMatrixPtr_;
};

} // namespace Loader
} // namespace Xyce

#endif

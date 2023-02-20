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
//
// Purpose        : This file contains class definitions for the base
//                  loader class.
//
// Special Notes  : This is a potentially confusing set of classes.  As of
//                  this writing(2/11/2007), loaders are used in 2 different
//                  interface layers:
//
//                  1) between the nonlinear solver and time integrator:
//                      Nonlinear Equation loader.
//
//                  2) between the time integrator and the device package:
//                      Ckt loader, MPDE, and SawTooth loaders
//
//                  This gets confusing, because 2 different purposes are being
//                  served.
//
//                  The Nonlinear Equation loader exists to insulate the nonlinear
//                  solver from having to know what kind of problem(in a math
//                  sense) this is.  For DCOP, the rhs consists of F(x) and B(t),
//                  but for transient F(x,t) , dQdt(x,t) and B(t). The Nonlinear Equation
//                  loader hides this decision from the nonlinear solver.
//
//                  The other layer(Ckt, MPDE and sawtooth) does 2 things.
//                 (1) it insulates the solvers from the specific physics.
//                  ie, the solvers(theoretically) don't know this is a
//                  circuit code.(2) We can hide MPDE details behind this
//                  interface, so the rest of the code can mostly be unaware
//                  if we are running an MPDE simulation or not.
//
//                  This structure is not perfect, and because it serves several
//                  different purposes, confusing.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_LOA_Loader_H
#define Xyce_LOA_Loader_H

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TIA_fwd.h>
#include <N_UTL_fwd.h>

#include <astRandEnum.h>

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Class         : Loader
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
class Loader
{
public:
  Loader()
  {}

  virtual ~Loader()
  {}

  // Get convergence info from devices
  virtual bool allDevicesConverged(Xyce::Parallel::Machine comm) = 0;

  // return value indicates whether mask is nontrivial
  // Virtual nonlinear new-DAE Matrices  load method.
  virtual bool loadDAEMatrices(
    Linear::Vector * tmpSolVectorPtr,
    Linear::Vector * tmpStaVectorPtr,
    Linear::Vector * tmpStaDerivVectorPtr,
    Linear::Vector * tmpStoVectorPtr,
    Linear::Matrix * tmpdQdxMatrixPtr,
    Linear::Matrix * tmpdFdxMatrixPtr,
    int loadType = Xyce::Device::ALL)
  {
    return false;
  }

  // Virtual nonlinear new-DAE vectors load method.
  virtual bool loadDAEVectors(
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
    int loadType = Xyce::Device::ALL)
  {
    return false;
  }

  // base class method does nothing.
  virtual bool loadDeviceErrorWeightMask(Linear::Vector * deviceMask) const
  {
    return false;
  }

  // Virtual method which initializes the nonlinear problem.
  virtual bool initializeProblem(
    Linear::Vector *    nextSolVectorPtr,
    Linear::Vector *    currSolVectorPtr,
    Linear::Vector *    lastSolVectorPtr,
    Linear::Vector *    nextStaVectorPtr,
    Linear::Vector *    currStaVectorPtr,
    Linear::Vector *    lastStaVectorPtr,
    Linear::Vector *    StateDerivVectorPtr,
    Linear::Vector *    nextStoVectorPtr,
    Linear::Vector *    currStoVectorPtr,
    Linear::Vector *    lastStoVectorPtr,
    Linear::Vector *    QVectorPtr,
    Linear::Vector *    FVectorPtr,
    Linear::Vector *    BVectorPtr,
    Linear::Vector *    dFdxdVpVectorPtr,
    Linear::Vector *    dQdxdVpVectorPtr) const = 0;

  // Virtual nonlinear new-DAE Matrices  apply method.
  // tmpdQdxVecVector = tmpdQdxMatrix * tmpVecVector
  // tmpdFdxVecVector = tmpdFdxMatrix * tmpVecVector
  virtual bool applyDAEMatrices(
    Linear::Vector * tmpSolVectorPtr,
    Linear::Vector * tmpStaVectorPtr,
    Linear::Vector * tmpStaDerivVectorPtr,
    Linear::Vector * tmpStoVectorPtr,
    const Linear::Vector & tmpVecVectorPtr,
    Linear::Vector * tmpdQdxVecVectorPtr,
    Linear::Vector * tmpdFdxVecVectorPtr)
  {
    return false;
  }

  virtual bool updateState(
    Linear::Vector * nextSolVectorPtr,
    Linear::Vector * currSolVectorPtr,
    Linear::Vector * lastSolVectorPtr,
    Linear::Vector * nextStaVectorPtr,
    Linear::Vector * currStaVectorPtr,
    Linear::Vector * lastStaVectorPtr,
    Linear::Vector * nextStoVectorPtr,
    Linear::Vector * currStoVectorPtr,
    Linear::Vector * lastStoVectorPtr,
    int loadType = Xyce::Device::ALL)
  {
    return false;
  }

  virtual bool loadBVectorsforAC(
    Linear::Vector * bVecRealPtr,
    Linear::Vector * bVecImagPtr)
  {
    return false;
  }

  virtual bool loadBVectorsforSources ()
  {
    return false;
  }

  virtual int getNumNoiseSources()
  {
    return 0;
  }

  virtual int getNumNoiseDevices()
  {
    return 0;
  }

  virtual void setupNoiseSources(std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec) 
  {}

  virtual void getNoiseSources(std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec) 
  {}

  virtual bool getBMatrixEntries(std::vector<int>& bMatEntriesVec, std::vector<int>& portVec, std::vector<double> * Z0sVec  = NULL )
  {
    return false;
  }

  // Virtual function for setting the initial guess.
  virtual bool setInitialGuess(Linear::Vector * solVectorPtr)
  {
    return false;
  }

  // Virtual function for setting a single parameter value.
  virtual bool setParam(std::string & name, double val, bool overrideOriginal = false) = 0; 

  virtual bool setParamRandomExpressionTerms2(
      const std::vector<Xyce::Analysis::SweepParam> & SamplingParams,
      bool overrideOriginal = false) { return true; };

  // Virtual function for getting a single parameter value.
  virtual double getParamAndReduce(Xyce::Parallel::Machine comm, const std::string & name) const = 0; 

  virtual void getRandomParams(std::vector<Xyce::Analysis::SweepParam> & SamplingParams, Parallel::Communicator & parallel_comm) {};

  // Virtual method which is called to update the sources.
  virtual bool updateSources()
  {
    return false;
  }
  // This is somewhat circuit-specific, unfortunately.
  virtual bool getLimiterFlag()
  {
    return false;
  }

  // Virtual method which gets the double DC Operating Point flag - used for
  // PDE devices.
  virtual bool isPDESystem() const
  {
    return false;
  }

  virtual bool outputPlotFiles() const
  {
    return false;
  }
  
  virtual bool finishOutput() const
  {
    return false;
  }

  // Virtual method which gets the time integration required breakpoint times
  //(in a vector).
  virtual bool getBreakPoints(
      std::vector< Util::BreakPoint > & breakPointTimes,
      std::vector< Util::BreakPoint > & pauseBreakPointTimes
      ) const
  {
    return false;
  }

  // Virtual accessor which returns the maximum time step size(in seconds).
  virtual double getMaxTimeStepSize()
  {
    return 0.0;
  }

  virtual void stepSuccess(Analysis::TwoLevelMode analysis)
  {}

  virtual void stepFailure(Analysis::TwoLevelMode analysis)
  {}

  //TVR: Method to be called when time integrator accepts a step, before any
  // tinkering with times or vectors
  virtual void acceptStep()
  {}

  virtual bool getInitialQnorm(std::vector<TimeIntg::TwoLevelError> & tleVec )
  {
    return false;
  }

  virtual bool getInnerLoopErrorSums(std::vector<TimeIntg::TwoLevelError> & tleVec ) const
  {
    return false;
  }

  virtual bool startTimeStep(
    bool                          beginIntegrationFlag,
    double                        nextTimeStep,
    double                        nextTime,
    int                           currentOrder)
  {
    return true;
  }

  virtual void setExternalSolverState(bool external_initJctFlag)
  {}

  virtual bool analyticBVecSensAvailable(const std::string & name) {return true;}
  virtual bool numericalBVecSensAvailable(const std::string & name) {return true;}

  virtual bool analyticMatrixSensitivitiesAvailable(const std::string & name) {return true;}
  virtual bool numericalMatrixSensitivitiesAvailable(const std::string & name) {return true;}

  virtual void getAnalyticBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const {}

  virtual void getNumericalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const {}

  virtual void getAnalyticMatrixSensitivities(
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const {}

  virtual void getNumericalMatrixSensitivities(
      const std::string & name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const {}

  // voltage limiter toggle functions
  virtual bool getVoltageLimiterStatus() = 0;
  virtual void setVoltageLimiterStatus(bool voltageLimterStatus) = 0;

  virtual void updateDependentParams () { return; }
  virtual void resetScaledParams() { return; }
};

} // namespace Loader
} // namespace Xyce

#endif

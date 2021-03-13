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
// Purpose        : MPDE Sawtooth IC Specific Loader
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/23/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_SawtoothLoader_H
#define Xyce_MPDE_SawtoothLoader_H

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------

#include <N_ANP_fwd.h>
#include <N_LAS_fwd.h>

#include <N_LOA_Loader.h>

#include <N_MPDE_State.h>

//-----------------------------------------------------------------------------
// Class         : N_MPDE_SawtoothLoader
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/23/04
//-----------------------------------------------------------------------------
class N_MPDE_SawtoothLoader : public Xyce::Loader::Loader
{
public:

  // Default constructor
  N_MPDE_SawtoothLoader(
    N_MPDE_State &                      state,
    Xyce::Analysis::AnalysisManager &   analysis_manager,
    Xyce::Loader::Loader &              loader)
  : state_(state),
    loader_(loader),
    analysisManager_(analysis_manager),
    timeShift_(0.0)
  {}

  // Destructor
  ~N_MPDE_SawtoothLoader() {}

  bool allDevicesConverged(Xyce::Parallel::Machine comm) { return loader_.allDevicesConverged(comm); } // ERK CHECK THIS

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( Xyce::Linear::Vector * X,
                        Xyce::Linear::Vector * S,
                        Xyce::Linear::Vector * dSdt,
                        Xyce::Linear::Vector * Store,
                        Xyce::Linear::Matrix * dQdx,
                        Xyce::Linear::Matrix * dFdx);

  // Method which is called to load the new-DAE vectors
  bool loadDAEVectors( Xyce::Linear::Vector * X,
                       Xyce::Linear::Vector * currX,
                       Xyce::Linear::Vector * lastX,
                       Xyce::Linear::Vector * nextS,
                       Xyce::Linear::Vector * currS,
                       Xyce::Linear::Vector * lastS,
                       Xyce::Linear::Vector * dSdt,
                       Xyce::Linear::Vector * nextStore,
                       Xyce::Linear::Vector * currStore,
                       Xyce::Linear::Vector * nextLeadFVectorPtr,
                       Xyce::Linear::Vector * nextLeadQVectorPtr,
                       Xyce::Linear::Vector * nextJunctionVVectorPtr,
                       Xyce::Linear::Vector * Q,
                       Xyce::Linear::Vector * F,
                       Xyce::Linear::Vector * B,
                       Xyce::Linear::Vector * dFdxdVp,
                       Xyce::Linear::Vector * dQdxdVp );

  bool updateState      (Xyce::Linear::Vector * nextSolVectorPtr,
                         Xyce::Linear::Vector * currSolVectorPtr,
                         Xyce::Linear::Vector * lastSolVectorPtr,
                         Xyce::Linear::Vector * nextStaVectorPtr,
                         Xyce::Linear::Vector * currStaVectorPtr,
                         Xyce::Linear::Vector * lastStaVectorPtr,
                         Xyce::Linear::Vector * nextStoVectorPtr,
                         Xyce::Linear::Vector * currStoVectorPtr
                         );

  // Virtual method which initializes the nonlinear problem.
  virtual bool initializeProblem( Xyce::Linear::Vector * nextSolVectorPtr,
                          Xyce::Linear::Vector * currSolVectorPtr,
                          Xyce::Linear::Vector * lastSolVectorPtr,
                          Xyce::Linear::Vector * nextStaVectorPtr,
                          Xyce::Linear::Vector * currStaVectorPtr,
                          Xyce::Linear::Vector * lastStaVectorPtr,
                          Xyce::Linear::Vector * StateDerivVectorPtr,
                          Xyce::Linear::Vector * nextStoVectorPtr,
                          Xyce::Linear::Vector * currStoVectorPtr,
                          Xyce::Linear::Vector * QVectorPtr,
                          Xyce::Linear::Vector * FVectorPtr,
                          Xyce::Linear::Vector * BVectorPtr,
                          Xyce::Linear::Vector * dFdxdVpVectorPtr,
                          Xyce::Linear::Vector * dQdxdVpVectorPtr) const  { return false; }

  void getAnalyticSensitivities(
      std::string & name, 
      std::vector<double> & dfdpVec, 
      std::vector<double> & dqdpVec,
      std::vector<double> & dbdpVec,
      std::vector<int> & FindicesVec,
      std::vector<int> & QindicesVec,
      std::vector<int> & BindicesVec) const
  {}

  bool analyticSensitivitiesAvailable (std::string & name) { return false; }

  virtual bool setParam (std::string & name, double val, bool overrideOriginal=false) { return false; }
  virtual double getParamAndReduce(Xyce::Parallel::Machine comm, const std::string & name) const { return 0.0; }

  // Assign times for fast time scale
  void setTimeShift( double timeShift )
  { timeShift_ = timeShift; }

  // Registration method for the system loader 
  // void registerAppLoader( Xyce::Loader::Loader *appLoader )
  // { loader_ = appLoader; }

  // // Registration method for the TIA
  // void registerTIA( Xyce::Analysis::AnalysisManager *anaInt )
  // { anaInt_ = anaInt; }

  // voltage limiter toggle functions (do nothing for now)
  bool getVoltageLimiterStatus() {return true;}
  void setVoltageLimiterStatus(bool voltageLimterStatus) {}

private :
  N_MPDE_State &                        state_;                 ///< MPDE State
  Xyce::Loader::Loader &                loader_;                ///< Base Application loader
  Xyce::Analysis::AnalysisManager &     analysisManager_;       ///< Analysis Manager
  double                                timeShift_;             ///< Time Shift for Fast Src
};

#endif


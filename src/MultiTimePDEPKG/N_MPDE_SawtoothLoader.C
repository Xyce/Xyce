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
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/23/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_MPDE_SawtoothLoader.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_ANP_AnalysisManager.h>

//-----------------------------------------------------------------------------
// Function      : N_MPDE_SawtoothLoader::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/23/04
//-----------------------------------------------------------------------------
bool N_MPDE_SawtoothLoader::loadDAEMatrices( Xyce::Linear::Vector * X,
                                             Xyce::Linear::Vector * S,
                                             Xyce::Linear::Vector * dSdt,
                                             Xyce::Linear::Vector * Store,
                                             Xyce::Linear::Matrix * dQdx,
                                             Xyce::Linear::Matrix * dFdx)
{
  //Zero out matrices
  dQdx->put(0.0);
  dFdx->put(0.0);

  double fastTime = analysisManager_.getTime() + timeShift_;
  //Set Time for fast time scale somewhere
  state_.fastTime = fastTime;

  return loader_.loadDAEMatrices( X, S, dSdt, Store, dQdx, dFdx);
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_SawtoothLoader::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_SawtoothLoader::loadDAEVectors( Xyce::Linear::Vector * X,
                                            Xyce::Linear::Vector * currX,
                                            Xyce::Linear::Vector * lastX,
                                            Xyce::Linear::Vector * nextS,
                                            Xyce::Linear::Vector * currS,
                                            Xyce::Linear::Vector * lastS,
                                            Xyce::Linear::Vector * dSdt,
                                            Xyce::Linear::Vector * nextStore,
                                            Xyce::Linear::Vector * currStore,
                                            Xyce::Linear::Vector * lastStore,
                                            Xyce::Linear::Vector * nextLeadFVectorPtr,
                                            Xyce::Linear::Vector * nextLeadQVectorPtr,
                                            Xyce::Linear::Vector * nextJunctionVVectorPtr,
                                            Xyce::Linear::Vector * Q,
                                            Xyce::Linear::Vector * F,
                                            Xyce::Linear::Vector * B,
                                            Xyce::Linear::Vector * dFdxdVp,
                                            Xyce::Linear::Vector * dQdxdVp )
{
  double fastTime = analysisManager_.getTime() + timeShift_;

  //Set Time for fast time scale somewhere
  state_.fastTime = fastTime;

  loader_.updateSources();

  return loader_.loadDAEVectors
    ( X, currX, lastX, nextS, currS, lastS, dSdt, 
      nextStore, currStore, lastStore, 
      nextLeadFVectorPtr, nextLeadQVectorPtr, 
      nextJunctionVVectorPtr, 
      Q, F, B, 
      dFdxdVp, dFdxdVp );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_SawtoothLoader::updateState
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_MPDE_SawtoothLoader::updateState (Xyce::Linear::Vector * nextSolVectorPtr,
                                         Xyce::Linear::Vector * currSolVectorPtr,
                                         Xyce::Linear::Vector * lastSolVectorPtr,
                                         Xyce::Linear::Vector * nextStaVectorPtr,
                                         Xyce::Linear::Vector * currStaVectorPtr,
                                         Xyce::Linear::Vector * lastStaVectorPtr,
                                         Xyce::Linear::Vector * nextStoVectorPtr,
                                         Xyce::Linear::Vector * currStoVectorPtr,
                                         Xyce::Linear::Vector * lastStoVectorPtr
                                         )
{
  return loader_.updateState
    ( nextSolVectorPtr, currSolVectorPtr, lastSolVectorPtr,
      nextStaVectorPtr, currStaVectorPtr, lastStaVectorPtr,
      nextStoVectorPtr, currStoVectorPtr, lastStoVectorPtr );
}


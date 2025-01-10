//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Builder for MPDE specific linear objects
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/11/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_MPDE_BUILDER_H
#define  Xyce_MPDE_BUILDER_H

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_LAS_Builder.h>
#include <N_MPDE_WarpedPhaseCondition.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;

// ---------- Forward Declarations ----------

class N_MPDE_Discretization;
class N_MPDE_Manager;

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Builder
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_MPDE_Builder : public Xyce::Linear::Builder
{

 public:

  // Default Constructor
  N_MPDE_Builder( N_MPDE_Manager &mgr,
                  const int Size,
                  const N_MPDE_Discretization &Disc,
                  const bool warpMPDE)
  : Size_(Size),
    Disc_(Disc),
    mpdeMgr_(mgr),
    warpMPDE_(warpMPDE),
    offset_(0),
    stateOffset_(0),
    storeOffset_(0),
    omegaGID_(0),
    phiGID_(0),
    augProcID_(-1)
  {}

  // Destructor
  virtual ~N_MPDE_Builder() {}

  void setWarpedPhaseCondition(const N_MPDE_WarpedPhaseCondition *warpMPDEPhasePtr)
  {
    warpMPDEPhasePtr_ = warpMPDEPhasePtr;
  }

  // Vector and Matrix creators

  // Vector factory 
  Xyce::Linear::Vector * createVector() const;

  // State Vector factory 
  Xyce::Linear::Vector * createStateVector() const;

  // Store Vector factory 
  Xyce::Linear::Vector * createStoreVector() const;

  // Lead Current Vector factory 
  Xyce::Linear::Vector * createLeadCurrentVector() const;

  // Matrix factory
  Xyce::Linear::Matrix * createMatrix() const;

  bool generateMaps( const RCP<Xyce::Parallel::ParMap>& BaseMap );

  bool generateStateMaps( const RCP<Xyce::Parallel::ParMap>& BaseStateMap );
  bool generateStoreMaps( const RCP<Xyce::Parallel::ParMap>& BaseStoreMap );
  bool generateLeadCurrentMaps( const RCP<Xyce::Parallel::ParMap>& BaseLeadCurrentMap );

  bool generateGraphs( const Xyce::Linear::Graph& BaseGraph );

  Teuchos::RCP<const Xyce::Parallel::ParMap> getSolutionMap() const;
  Teuchos::RCP<Xyce::Parallel::ParMap> getSolutionMap();

  // Return GID offset for blocks to Manager for construction of Loader
  int getMPDEOffset()
  { return offset_; }

  // Return omega GID to Manager for construction of Loader
  int getMPDEomegaGID()
  { return omegaGID_; }

  // Return phi GID to Manager for construction of Loader
  int getMPDEphiGID()
  { return phiGID_; }

  // Return processor id for augmented systems (i.e. warped MPDE)
  int getMPDEaugProcID()
  { return augProcID_; }

private:

  const int Size_;
  const N_MPDE_Discretization &   Disc_;
  
  // mpde manager:
  N_MPDE_Manager &              mpdeMgr_;
  
  const N_MPDE_WarpedPhaseCondition * warpMPDEPhasePtr_;

  bool warpMPDE_;
  int offset_, stateOffset_, storeOffset_, leadCurrentOffset_;
  int omegaGID_;
  int phiGID_;
  int augProcID_;

  RCP<Xyce::Parallel::ParMap> BaseMap_, BaseStateMap_;
  RCP<Xyce::Parallel::ParMap> MPDEMap_, MPDEStateMap_;
  RCP<Xyce::Parallel::ParMap> BaseStoreMap_;
  RCP<Xyce::Parallel::ParMap> BaseLeadCurrentMap_;
  RCP<Xyce::Parallel::ParMap> MPDEStoreMap_;
  RCP<Xyce::Parallel::ParMap> MPDELeadCurrentMap_;

  RCP<Xyce::Linear::Graph> BaseGraph_;
  RCP<Xyce::Linear::Graph> MPDEGraph_;
};

#endif



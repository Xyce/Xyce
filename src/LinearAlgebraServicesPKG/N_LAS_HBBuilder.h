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
// Purpose        : Builder for HB specific linear objects
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_HBBUILDER_H
#define  Xyce_LAS_HBBUILDER_H

// ---------- Standard Includes ----------

#include <string>
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_LAS_Builder.h>

#include <Teuchos_RCP.hpp>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : HBBuilder
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
class HBBuilder : public Builder
{

 public:

  // Default Constructor
  
  HBBuilder( const int Size, const bool hbOsc);

  // Destructor
  ~HBBuilder() {}

  // Vector and Matrix creators

  // Vector factory 
  Vector * createVector() const;

  // State Vector factory 
  Vector * createStateVector() const;

  // Store Vector factory 
  Vector * createStoreVector() const;

  // Lead Current Vector factory 
  Vector * createLeadCurrentVector() const;

  // HB time-domain block vector creation:
  Teuchos::RCP<BlockVector> createTimeDomainBlockVector() const;
  Teuchos::RCP<BlockVector> createTimeDomainStateBlockVector() const;
  Teuchos::RCP<BlockVector> createTimeDomainStoreBlockVector() const;
  Teuchos::RCP<BlockVector> createTimeDomainLeadCurrentBlockVector() const;
  
  Teuchos::RCP<BlockVector> createExpandedRealFormBlockVector() const;
  Teuchos::RCP<BlockVector> createExpandedRealFormTransposeBlockVector() const;
  Teuchos::RCP<BlockVector> createExpandedRealFormTransposeStateBlockVector() const;
  Teuchos::RCP<BlockVector> createExpandedRealFormTransposeStoreBlockVector() const;
  Teuchos::RCP<BlockVector> createExpandedRealFormTransposeLeadCurrentBlockVector() const;

  // Matrix factory
  Matrix * createMatrix() const { return 0; }

  bool generateMaps( const Teuchos::RCP<Parallel::ParMap>& BaseMap, 
                     const Teuchos::RCP<Parallel::ParMap>& oBaseMap );

  bool generateStateMaps( const Teuchos::RCP<Parallel::ParMap>& BaseStateMap );
  bool generateStoreMaps( const Teuchos::RCP<Parallel::ParMap>& BaseStoreMap );
  bool generateLeadCurrentMaps( const Teuchos::RCP<Parallel::ParMap>& BaseLeadCurrentMap );

  bool generateGraphs( const Graph& baseFullGraph );

  // Return maps for Harmonic Balance linear system.
  Teuchos::RCP<const Parallel::ParMap> getSolutionMap() const;
  Teuchos::RCP<Parallel::ParMap> getSolutionMap();

  Teuchos::RCP<Parallel::ParMap> getSolutionOverlapMap() const
  { return HBExpandedRealFormBVOverlapMap_; }

  // Return the base map for each block in the expanded maps (a.k.a. time-domain maps)
  Teuchos::RCP<const Parallel::ParMap> getBaseSolutionMap() const
  { return BaseMap_; }

  Teuchos::RCP<const Parallel::ParMap> getBaseStateMap() const
  { return BaseStateMap_; }

  Teuchos::RCP<const Parallel::ParMap> getBaseStoreMap() const
  { return BaseStoreMap_; }
  
  Teuchos::RCP<const Parallel::ParMap> getBaseLeadCurrentMap() const
  { return BaseLeadCurrentMap_; }

  const std::vector<int>& getAugmentedLIDs() const
  { return augmentedLIDs_; }

  // Return GID offset for blocks for construction of Loader
  int getHBOffset()
  { return offset_; }

  int getHBStateOffset()
  { return stateOffset_; }

  int getHBStoreOffset()
  { return storeOffset_; }
  
  int getHBLeadCurrentOffset()
  { return leadCurrentOffset_; }

  int getNumHarmonics()
  { return numHarmonics_; }

private:

  const int numHarmonics_;
  int numSolVariables_, numStateVariables_;
  int numStoreVariables_;
  int numLeadCurrentVariables_;

  int offset_, stateOffset_;
  int storeOffset_;

  bool hbOsc_;
  std::vector<int> augmentedLIDs_;
  std::vector<int> augmentedOverlapLIDs_;
//  std::vector<int>* augmentedLIDs_;
  int leadCurrentOffset_;

  // HB maps for block vectors:
  // numBlocks = 2*(number of harmonics), numElem = number of solution variables
  Teuchos::RCP<Parallel::ParMap> HBExpandedRealFormBVMap_; 
  Teuchos::RCP<Parallel::ParMap> HBExpandedRealFormBVOverlapMap_; 
  Teuchos::RCP<Parallel::ParMap> HBExpandedRealFormStateBVMap_;
  Teuchos::RCP<Parallel::ParMap> HBExpandedRealFormStoreBVMap_;
  Teuchos::RCP<Parallel::ParMap> HBExpandedRealFormLeadCurrentBVMap_;
  
  // numBlocks = number of solution variables, numElem = 2*(number of harmonics) 
  // We don't need a special map here, its the same as the non-tranpose

  Teuchos::RCP<Parallel::ParMap> BaseMap_, oBaseMap_;

  Teuchos::RCP<Parallel::ParMap> BaseStateMap_;
  Teuchos::RCP<Parallel::ParMap> BaseStoreMap_;
  Teuchos::RCP<Parallel::ParMap> BaseLeadCurrentMap_;

  Teuchos::RCP<Graph> baseFullGraph_;

  Teuchos::RCP<Parallel::ParMap> HBMap_, oHBMap_;

  Teuchos::RCP<Parallel::ParMap> HBStateMap_;
  Teuchos::RCP<Parallel::ParMap> HBStoreMap_;
  Teuchos::RCP<Parallel::ParMap> HBLeadCurrentMap_;

};

} // namespace Linear
} // namespace Xyce

#endif

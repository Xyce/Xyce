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
// Purpose        : Builder for ES specific linear objects
//
// Special Notes  : This version produces a block structure in which the 
//                  outer loop is parameters, and the inner loop is the original 
//                  (non-block) circuit matrix structure.
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 05/31/2018
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_ESBUILDER_H
#define  Xyce_LAS_ESBUILDER_H

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
// Class         : ESBuilder
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
class ESBuilder : public Builder
{
 public:

  // Default Constructor
  
  ESBuilder( const int Size );

  // Destructor
  ~ESBuilder() {}

  // Vector and Matrix creators

  // Vector factory 
  Vector * createVector() const;

  // MultiVector factory 
  MultiVector * createMultiVector( int numVectors = 1 ) const;

  // State Vector factory 
  Vector * createStateVector() const;

  // Store Vector factory 
  Vector * createStoreVector() const;

  // Lead Current Vector factory 
  Vector * createLeadCurrentVector() const;

  // Matrix factory
  Matrix * createMatrix() const;

  //Coloring Assoc with Variable Types in Solution Vector
  const std::vector<int> & createSolnColoring() const;

  //Coloring needed for imposing .IC and .NODESET
  const std::vector<int> & createInitialConditionColoring() const;

  // Convert topology op data to analysis specific op data
  bool createInitialConditionOp( std::map<int,double> & op ) const;

  bool createInitialConditionOp( std::vector<int> & op ) const;

  bool generateMaps( const Teuchos::RCP<Parallel::ParMap>& BaseMap, 
                     const Teuchos::RCP<Parallel::ParMap>& oBaseMap );

  bool generateStateMaps( const Teuchos::RCP<Parallel::ParMap>& BaseStateMap );
  bool generateStoreMaps( const Teuchos::RCP<Parallel::ParMap>& BaseStoreMap );
  bool generateLeadCurrentMaps( const Teuchos::RCP<Parallel::ParMap>& BaseLeadCurrentMap );

  bool generateGraphs( const Graph& BaseFullGraph );

  // Return maps for sampling linear system.
  Teuchos::RCP<const Parallel::ParMap> getSolutionMap() const
  { return( ESMap_ ); }

  Teuchos::RCP<Parallel::ParMap> getSolutionMap()
  { return( ESMap_ ); }

  Teuchos::RCP<Parallel::ParMap> getSolutionOverlapMap() const
  { return oESMap_; }

  const std::vector<int> & vnodeGIDVec() const;

  // Return the base map for each block in the expanded maps (a.k.a. time-domain maps)
  Teuchos::RCP<const Parallel::ParMap> getBaseSolutionMap() const
  { return BaseMap_; }

  Teuchos::RCP<const Parallel::ParMap> getBaseStateMap() const
  { return BaseStateMap_; }

  Teuchos::RCP<const Parallel::ParMap> getBaseStoreMap() const
  { return BaseStoreMap_; }
  
  Teuchos::RCP<const Parallel::ParMap> getBaseLeadCurrentMap() const
  { return BaseLeadCurrentMap_; }

  // Return GID offset for blocks for construction of Loader
  int getESOffset()
  { return offset_; }

  int getESStateOffset()
  { return stateOffset_; }

  int getESStoreOffset()
  { return storeOffset_; }
  
  int getESLeadCurrentOffset()
  { return leadCurrentOffset_; }

  int getNumSamples()
  { return numSamples_; }

private:
  const int numSamples_;
  int numSolVariables_, numStateVariables_;
  int numStoreVariables_;
  int numLeadCurrentVariables_;

  int offset_, stateOffset_;
  int storeOffset_;
  int leadCurrentOffset_;

  mutable std::vector<int> vnodeVec_;

  // ES maps for block vectors (BV):
 // numBlocks = number of samples, numElem = number of solution variables
  
  // numBlocks = number of solution variables, numElem = number of samples
  Teuchos::RCP<Parallel::ParMap> BaseMap_, oBaseMap_;

  Teuchos::RCP<Parallel::ParMap> BaseStateMap_;
  Teuchos::RCP<Parallel::ParMap> BaseStoreMap_;
  Teuchos::RCP<Parallel::ParMap> BaseLeadCurrentMap_;

  Teuchos::RCP<Graph> baseFullGraph_;
  Teuchos::RCP<Graph> blockGraph_;

  std::vector<std::vector<int> > blockPattern_;

  Teuchos::RCP<Parallel::ParMap> ESMap_, oESMap_;
  Teuchos::RCP<Parallel::ParMap> ESStateMap_;
  Teuchos::RCP<Parallel::ParMap> ESStoreMap_;
  Teuchos::RCP<Parallel::ParMap> ESLeadCurrentMap_;
};

} // namespace Linear
} // namespace Xyce

#endif

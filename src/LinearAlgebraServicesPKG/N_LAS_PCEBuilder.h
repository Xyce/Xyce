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
// Purpose        : Builder for PCE specific linear objects
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 06/27/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_PCEBUILDER_H
#define  Xyce_LAS_PCEBUILDER_H

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
// Class         : PCEBuilder
// Purpose       :
//
// Special Notes : Unlike some other block-analysis builders, this builder 
//                 supports two different block structures.  
//
//                 (1) The PCE expansion.  This is the main one, which 
//                 creates the block system objects used by the nonlinear 
//                 solver.  This structure is used by all the default 
//                 functions such as "createMatrix", "createVector" etc.
//
//                 (2) The quadrature points, which are needed in order 
//                 to evaluate the integrals that are computed to set up 
//                 the matrix and vector entries into the PCE block system 
//                 described by (1).  This structure is used by functions 
//                 and data with the word "quad" in the name, such as 
//                 "createQuadVector".  The only part of the code that 
//                 will use these objects is the PCELoader.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 06/27/2019
//-----------------------------------------------------------------------------
class PCEBuilder : public Builder
{
 public:

  // Default Constructor
  PCEBuilder( const int Size, const int quadPointsSize);

  // Destructor
  ~PCEBuilder() {}

  // Vector and Matrix creators

  // Vector factory 
  Vector * createVector() const;
  BlockVector * createQuadVector() const;

  // State Vector factory 
  Vector * createStateVector() const;

  // Store Vector factory 
  Vector * createStoreVector() const;

  // Lead Current Vector factory 
  Vector * createLeadCurrentVector() const;

  // Matrix factory
  Matrix * createMatrix() const;
  BlockMatrix * createQuadMatrix() const;

#if 0
  // ERK. The following functions need to be implemented to support 
  // various solver methods with intrusive PCE.  Without them, things
  // like .IC and gmin stepping will not work.

  //Coloring Assoc with Variable Types in Solution Vector
  const std::vector<int> & createSolnColoring() const;

  //Coloring needed for imposing .IC and .NODESET
  const std::vector<int> & createInitialConditionColoring() const;

  // Convert topology op data to analysis specific op data
  bool createInitialConditionOp( std::map<int,double> & op ) const;

  bool createInitialConditionOp( std::vector<int> & op ) const;

  const std::vector<int> & vnodeGIDVec() const;
#else
  const std::vector<int> & vnodeGIDVec() const;
#endif

  bool generateMaps( const Teuchos::RCP<Parallel::ParMap>& BaseMap, 
                     const Teuchos::RCP<Parallel::ParMap>& oBaseMap );

  bool generateStateMaps( const Teuchos::RCP<Parallel::ParMap>& BaseStateMap );
  bool generateStoreMaps( const Teuchos::RCP<Parallel::ParMap>& BaseStoreMap );
  bool generateLeadCurrentMaps( const Teuchos::RCP<Parallel::ParMap>& BaseLeadCurrentMap );

  bool generateGraphs( 
    const Graph& pceGraph,
    const Graph& baseFullGraph );

  // Return maps for sampling linear system.
  Teuchos::RCP<const Parallel::ParMap> getSolutionMap() const
  { return( PCEMap_ ); }

  Teuchos::RCP<Parallel::ParMap> getSolutionMap()
  { return( PCEMap_ ); }

  Teuchos::RCP<Parallel::ParMap> getSolutionOverlapMap() const
  { return oPCEMap_; }

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
  int getPCEOffset()
  { return offset_; }

  int getPCEStateOffset()
  { return stateOffset_; }

  int getPCEStoreOffset()
  { return storeOffset_; }
  
  int getPCELeadCurrentOffset()
  { return leadCurrentOffset_; }

  int getNumBlockRows ()
  { return numBlockRows_; }

private:
  const int numBlockRows_;
  const int numQuadPoints_;
  int numSolVariables_, numStateVariables_;
  int numStoreVariables_;
  int numLeadCurrentVariables_;

  int offset_, stateOffset_;
  int storeOffset_;
  int leadCurrentOffset_;

  mutable std::vector<int> vnodeVec_;

  // PCE maps for block vectors (BV):
  Teuchos::RCP<Parallel::ParMap> BaseMap_, oBaseMap_;

  Teuchos::RCP<Parallel::ParMap> BaseStateMap_;
  Teuchos::RCP<Parallel::ParMap> BaseStoreMap_;
  Teuchos::RCP<Parallel::ParMap> BaseLeadCurrentMap_;

  Teuchos::RCP<Graph> pceGraph_;
  Teuchos::RCP<Graph> baseFullGraph_;
  Teuchos::RCP<Graph> blockGraph_;
  Teuchos::RCP<Graph> quadBlockGraph_;

  std::vector<std::vector<int> > blockPattern_;
  std::vector<std::vector<int> > quadBlockPattern_;

  Teuchos::RCP<Parallel::ParMap> PCEMap_, oPCEMap_;
  Teuchos::RCP<Parallel::ParMap> quadMap_, oquadMap_;

  Teuchos::RCP<Parallel::ParMap> PCEStateMap_;
  Teuchos::RCP<Parallel::ParMap> quadStateMap_;

  Teuchos::RCP<Parallel::ParMap> PCEStoreMap_;
  Teuchos::RCP<Parallel::ParMap> quadStoreMap_;

  Teuchos::RCP<Parallel::ParMap> PCELeadCurrentMap_;
  Teuchos::RCP<Parallel::ParMap> quadLeadCurrentMap_;
};

} // namespace Linear
} // namespace Xyce

#endif

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
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_PCEBuilder.h>

#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_QueryUtil.h>
#include <N_LAS_Graph.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_fwd.h>
#include <N_UTL_FeatureTest.h>

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Utils.hpp>

#include <N_ANP_UQSupport.h>

using Teuchos::rcp;
using Teuchos::RCP;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::PCEBuilder
// Purpose       : Default constructor 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
PCEBuilder::PCEBuilder( const int Size, const int quadPointsSize )
: numBlockRows_(Size),
  numQuadPoints_(quadPointsSize), // ERK. FIX THIS
  numSolVariables_(0),
  numStateVariables_(0),
  numStoreVariables_(0),
  offset_(0),
  stateOffset_(0),
  storeOffset_(0),
  leadCurrentOffset_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createVector() const
{
  return Xyce::Linear::createBlockVector( numBlockRows_, PCEMap_, BaseMap_ );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createQuadVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
BlockVector * PCEBuilder::createQuadVector() const
{
  return Xyce::Linear::createBlockVector( numQuadPoints_, quadMap_, BaseMap_ );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Matrix * PCEBuilder::createMatrix() const
{
  return Xyce::Linear::createBlockMatrix( numBlockRows_, offset_, blockPattern_, blockGraph_.get(), baseFullGraph_.get() );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createQuadMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
BlockMatrix * PCEBuilder::createQuadMatrix() const
{
  return Xyce::Linear::createBlockMatrix( numQuadPoints_, offset_, quadBlockPattern_, quadBlockGraph_.get(), baseFullGraph_.get() );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createStateVector() const
{
  return Xyce::Linear::createBlockVector( numQuadPoints_, quadStateMap_, BaseStateMap_ );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createStoreVector() const
{
  return Xyce::Linear::createBlockVector( numQuadPoints_, quadStoreMap_, BaseStoreMap_ );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createLeadCurrentVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createLeadCurrentVector() const
{
  return Xyce::Linear::createBlockVector( numQuadPoints_, quadLeadCurrentMap_, BaseLeadCurrentMap_ );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
bool PCEBuilder::generateMaps( const RCP<Parallel::ParMap>& BaseMap, 
                               const RCP<Parallel::ParMap>& oBaseMap )
{
  //Save copy of base map
  BaseMap_ = BaseMap;
  oBaseMap_ = oBaseMap;

  //Determine block offset
  offset_ = BaseMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  PCEMap_ = Linear::createBlockParMap(numBlockRows_, *BaseMap, 0, 0, offset_);

  // Use the block linear system helper to create the block parallel maps
  quadMap_ = Linear::createBlockParMap(numQuadPoints_, *BaseMap, 0, 0, offset_);

  // Helpful names for various sizes (subtract 1 for ground node):
  numSolVariables_ = oBaseMap_->numLocalEntities()-1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::generateStateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
bool PCEBuilder::generateStateMaps( const RCP<Parallel::ParMap>& BaseStateMap )
{
  //Save copy of base map
  BaseStateMap_ = BaseStateMap;

  //determine block offset
  stateOffset_ = BaseStateMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  PCEStateMap_ = createBlockParMap(numBlockRows_, *BaseStateMap);

  quadStateMap_ = createBlockParMap(numQuadPoints_, *BaseStateMap);

  // Helpful names for various sizes:
  numStateVariables_ = BaseStateMap_->numGlobalEntities();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::generateStoreMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
bool PCEBuilder::generateStoreMaps( const RCP<Parallel::ParMap>& BaseStoreMap )
{
  //Save copy of base map
  BaseStoreMap_ = BaseStoreMap;

  //determine block offset
  storeOffset_ = BaseStoreMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  PCEStoreMap_ = createBlockParMap(numBlockRows_, *BaseStoreMap);

  quadStoreMap_ = createBlockParMap(numQuadPoints_, *BaseStoreMap);

  // Helpful names for various sizes:
  numStoreVariables_ = BaseStoreMap_->numGlobalEntities();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::generateLeadCurrentMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
bool PCEBuilder::generateLeadCurrentMaps( const RCP<Parallel::ParMap>& BaseLeadCurrentMap )
{
  //Save copy of base map
  BaseLeadCurrentMap_ = BaseLeadCurrentMap;

  //determine block offset
  leadCurrentOffset_ = BaseLeadCurrentMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  PCELeadCurrentMap_ = createBlockParMap(numBlockRows_, *BaseLeadCurrentMap);

  quadLeadCurrentMap_ = createBlockParMap(numQuadPoints_, *BaseLeadCurrentMap);

  // Helpful names for various sizes:
  numLeadCurrentVariables_ = BaseLeadCurrentMap_->numGlobalEntities();

  return true;
}
//-----------------------------------------------------------------------------
// Function      : PCEBuilder::generateGraphs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
bool PCEBuilder::generateGraphs( 
    const Graph& pceGraph,
    const Graph& baseFullGraph 
    )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("PCEBuilder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of graphs
  pceGraph_ = rcp( pceGraph.cloneCopy() );
  baseFullGraph_ = rcp( baseFullGraph.cloneCopy() );

  int numBlockRows = numBlockRows_;
  blockPattern_.clear();
  blockPattern_.resize(numBlockRows);

  // the coefficient (PCE) pattern is dense.  At first I tried to use a 
  // sparsity problem but ran into problems.  So for now, it is dense.
  for (int i=0;i<numBlockRows;++i)
  {
    blockPattern_[i].clear();
    blockPattern_[i].resize(numBlockRows);

#if 0
    // sparse version
    int maxIndices = pceGraph.maxNumIndices();
    std::vector<int> indices(maxIndices);
    int numIndices=0;
    int pceRow = pceGraph.localToGlobalIndex(i);
    pceGraph.extractGlobalRowCopy( pceRow, maxIndices, numIndices, &indices[0] );
    blockPattern_[i].resize(numIndices,0);
    for (int j=0;j<numIndices;++j)
    {
      int col = indices[j];
      blockPattern_[i][col] = j; 
    }
#else
    // dense version
    for (int j=0;j<numBlockRows;++j)
    {
      blockPattern_[i][j] = j; 
    }
#endif
  }

  // the quad pattern is a block diagonal
  quadBlockPattern_.clear();
  quadBlockPattern_.resize(numQuadPoints_);
  for (int i=0;i<numQuadPoints_;++i)
  {
    quadBlockPattern_[i].resize(1);
    quadBlockPattern_[i][0] = i; 
  }

  blockGraph_ = Linear::createBlockGraph(offset_, blockPattern_, *PCEMap_, *baseFullGraph_ );

  quadBlockGraph_ = Linear::createBlockGraph(offset_, quadBlockPattern_, *quadMap_, *baseFullGraph_ );

  return true;
}

//-----------------------------------------------------------------------------
// ERK.  4/16/2021. This function is only here to serve as an error trap and
// prevent a seg fault.
//
// The intrusive PCE method, which this class supports, cannot support gmin 
// stepping without this function.  But, if this function isn't implemented
// present, then gmin stepping will cause a seg fault.  This is a bit of a
// problem because gmin stepping is attempted automatically in Xyce when 
// standard Newton fails the DCOP.
//
// There is an implemented vnodeGILDVec function for in the ESBuilder,
// but naively copying that function to this class won't work.  Doing it right
// will require a little thought.
//-----------------------------------------------------------------------------
const std::vector<int> & PCEBuilder::vnodeGIDVec() const
{
  Report::DevelFatal().in("PCEBuilder::vnodeGIDVec")
  << "This function is not implemented yet.";

  return vnodeVec_;
}

} // namespace Linear
} // namespace Xyce

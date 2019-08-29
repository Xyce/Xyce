//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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

#include <N_PDS_ParMap.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_fwd.h>
#include <N_UTL_FeatureTest.h>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
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
Vector * PCEBuilder::createVector( double initialValue ) const
{
  RCP<Vector> vector = createBlockVector(); 
  vector.release(); // Release ownership of the object.
  return(&*vector);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numBlockRows_, PCEMap_, BaseMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeBlockVector() const
{
  RCP<BlockVector> vec = rcp(
      new BlockVector( numBlockRows_, PCEMap_ )
      );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeStateBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeStateBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numBlockRows_, PCEStateMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeStoreBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeStoreBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numBlockRows_, PCEStoreMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeLeadCurrentBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeLeadCurrentBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numBlockRows_, PCELeadCurrentMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createQuadVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createQuadVector( double initialValue ) const
{
  RCP<Vector> vector = createQuadBlockVector(); 
  vector.release(); // Release ownership of the object.
  return(&*vector);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createQuadBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createQuadBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numQuadPoints_, quadMap_, BaseMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeQuadBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeQuadBlockVector() const
{
  RCP<BlockVector> vec = rcp(
      new BlockVector( numQuadPoints_, quadMap_ )
      );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeStateQuadBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeStateQuadBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numQuadPoints_, quadStateMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeStoreQuadBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeStoreQuadBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numQuadPoints_, quadStoreMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createTransposeLeadCurrentQuadBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
RCP<BlockVector> PCEBuilder::createTransposeLeadCurrentQuadBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numQuadPoints_, quadLeadCurrentMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Matrix * PCEBuilder::createMatrix( double initialValue ) const
{
  RCP<Matrix> matrix = createBlockMatrix( initialValue );
  matrix.release(); // Release ownership of the object.
  return(&*matrix);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createBlockMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Teuchos::RCP<BlockMatrix> PCEBuilder::createBlockMatrix( double initialValue ) const
{
  return rcp (new Linear::BlockMatrix( numBlockRows_, offset_, blockPattern_, *blockGraph_, *BaseFullGraph_) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createQuadMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
Matrix * PCEBuilder::createQuadMatrix( double initialValue ) const
{
  RCP<Matrix> matrix = createQuadBlockMatrix( initialValue );
  matrix.release(); // Release ownership of the object.
  return(&*matrix);
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createQuadBlockMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
Teuchos::RCP<BlockMatrix> PCEBuilder::createQuadBlockMatrix( double initialValue ) const
{
  return rcp (new Linear::BlockMatrix( numQuadPoints_, offset_, quadBlockPattern_, *quadBlockGraph_, *BaseFullGraph_) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createStateVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numBlockRows_, PCEStateMap_, BaseStateMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createStoreVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numBlockRows_, PCEStoreMap_, BaseStoreMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createLeadCurrentVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createLeadCurrentVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numBlockRows_, PCELeadCurrentMap_, BaseLeadCurrentMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createStateQuadVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createStateQuadVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numQuadPoints_, quadStateMap_, BaseStateMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createStoreQuadVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createStoreQuadVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numQuadPoints_, quadStoreMap_, BaseStoreMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::createLeadCurrentQuadVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 8/25/2019
//-----------------------------------------------------------------------------
Vector * PCEBuilder::createLeadCurrentQuadVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numQuadPoints_, quadLeadCurrentMap_, BaseLeadCurrentMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : PCEBuilder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/27/2019
//-----------------------------------------------------------------------------
bool PCEBuilder::generateMaps( const RCP<N_PDS_ParMap>& BaseMap, 
                               const RCP<N_PDS_ParMap>& oBaseMap )
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
bool PCEBuilder::generateStateMaps( const RCP<N_PDS_ParMap>& BaseStateMap )
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
bool PCEBuilder::generateStoreMaps( const RCP<N_PDS_ParMap>& BaseStoreMap )
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
bool PCEBuilder::generateLeadCurrentMaps( const RCP<N_PDS_ParMap>& BaseLeadCurrentMap )
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
    const Epetra_CrsGraph & pceGraph,
    const Epetra_CrsGraph & BaseFullGraph 
    )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("PCEBuilder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of graphs
  pceGraph_ = rcp(new Epetra_CrsGraph( pceGraph ));
  BaseFullGraph_ = rcp(new Epetra_CrsGraph( BaseFullGraph ));

  int numBlockRows = numBlockRows_;
  blockPattern_.clear();
  blockPattern_.resize(numBlockRows);

  for (int i=0;i<numBlockRows;++i)
  {
    blockPattern_[i].resize(numBlockRows);

#if 0
    int maxIndices = pceGraph.MaxNumIndices();
    std::vector<int> indices(maxIndices);
    int numIndices=0;
    int pceRow = pceGraph.GRID(i);
    pceGraph.ExtractGlobalRowCopy( pceRow, maxIndices, numIndices, &indices[0] );
    blockPattern_[i].resize(numIndices,0);
    for (int j=0;j<numIndices;++j)
    {
      int col = indices[j];
      blockPattern_[i][col] = j; 
    }
#else
    // making this dense for now
    for (int j=0;j<numBlockRows;++j)
    {
      blockPattern_[i][j] = j; 
    }
#endif
  }

  quadBlockPattern_.clear();
  quadBlockPattern_.resize(numQuadPoints_);
  for (int i=0;i<numQuadPoints_;++i)
  {
    quadBlockPattern_[i].resize(1);
    quadBlockPattern_[i][0] = i; 
  }

  blockGraph_ = Linear::createBlockGraph(offset_, blockPattern_, *PCEMap_, *BaseFullGraph_ );

  quadBlockGraph_ = Linear::createBlockGraph(offset_, quadBlockPattern_, *quadMap_, *BaseFullGraph_ );

  return true;
}

} // namespace Linear
} // namespace Xyce

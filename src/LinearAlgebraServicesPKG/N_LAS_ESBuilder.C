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
// Creation Date : 5/31/2018
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_ESBuilder.h>

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

using Teuchos::rcp;
using Teuchos::RCP;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : ESBuilder::ESBuilder
// Purpose       : Default constructor 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
ESBuilder::ESBuilder( const int Size )
: numSamples_(Size),
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
// Function      : ESBuilder::createVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
Vector * ESBuilder::createVector( double initialValue ) const
{
  RCP<Vector> vector = createBlockVector(); 
  vector.release(); // Release ownership of the object.
  return(&*vector);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
RCP<BlockVector> ESBuilder::createBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numSamples_, ESMap_, BaseMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createTransposeBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
RCP<BlockVector> ESBuilder::createTransposeBlockVector() const
{
  RCP<BlockVector> vec = rcp(
      new BlockVector( numSamples_, ESMap_ )
      );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createTransposeStateBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
RCP<BlockVector> ESBuilder::createTransposeStateBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numSamples_, ESStateMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createTransposeStoreBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
RCP<BlockVector> ESBuilder::createTransposeStoreBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numSamples_, ESStoreMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createTransposeLeadCurrentBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
RCP<BlockVector> ESBuilder::createTransposeLeadCurrentBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( numSamples_, ESLeadCurrentMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/14/2018
//-----------------------------------------------------------------------------
Matrix * ESBuilder::createMatrix( double initialValue ) const
{
  RCP<Matrix> matrix = createBlockMatrix( initialValue );
  matrix.release(); // Release ownership of the object.
  return(&*matrix);
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createBlockMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/14/2018
//-----------------------------------------------------------------------------
Teuchos::RCP<BlockMatrix> ESBuilder::createBlockMatrix( double initialValue ) const
{
  return rcp (new Linear::BlockMatrix( numSamples_, offset_, blockPattern_, *blockGraph_, *BaseFullGraph_) );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
Vector * ESBuilder::createStateVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numSamples_, ESStateMap_, BaseStateMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
Vector * ESBuilder::createStoreVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numSamples_, ESStoreMap_, BaseStoreMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createLeadCurrentVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
Vector * ESBuilder::createLeadCurrentVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numSamples_, ESLeadCurrentMap_, BaseLeadCurrentMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
bool ESBuilder::generateMaps( const RCP<N_PDS_ParMap>& BaseMap, 
                                    const RCP<N_PDS_ParMap>& oBaseMap )
{
  //Save copy of base map
  BaseMap_ = BaseMap;
  oBaseMap_ = oBaseMap;

  //Determine block offset
  offset_ = BaseMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  ESMap_ = Linear::createBlockParMap(numSamples_, *BaseMap, 0, 0, offset_);

  // Helpful names for various sizes (subtract 1 for ground node):
  numSolVariables_ = oBaseMap_->numLocalEntities()-1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::generateStateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
bool ESBuilder::generateStateMaps( const RCP<N_PDS_ParMap>& BaseStateMap )
{
  //Save copy of base map
  BaseStateMap_ = BaseStateMap;

  //determine block offset
  stateOffset_ = BaseStateMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  ESStateMap_ = createBlockParMap(numSamples_, *BaseStateMap);

  // Helpful names for various sizes:
  numStateVariables_ = BaseStateMap_->numGlobalEntities();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::generateStoreMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
bool ESBuilder::generateStoreMaps( const RCP<N_PDS_ParMap>& BaseStoreMap )
{
  //Save copy of base map
  BaseStoreMap_ = BaseStoreMap;

  //determine block offset
  storeOffset_ = BaseStoreMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  ESStoreMap_ = createBlockParMap(numSamples_, *BaseStoreMap);

  // Helpful names for various sizes:
  numStoreVariables_ = BaseStoreMap_->numGlobalEntities();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::generateLeadCurrentMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
bool ESBuilder::generateLeadCurrentMaps( const RCP<N_PDS_ParMap>& BaseLeadCurrentMap )
{
  //Save copy of base map
  BaseLeadCurrentMap_ = BaseLeadCurrentMap;

  //determine block offset
  leadCurrentOffset_ = BaseLeadCurrentMap_->maxGlobalEntity() + 1;  // Use this offset to create a contiguous gid map for direct solvers.

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  ESLeadCurrentMap_ = createBlockParMap(numSamples_, *BaseLeadCurrentMap);

  // Helpful names for various sizes:
  numLeadCurrentVariables_ = BaseLeadCurrentMap_->numGlobalEntities();

  return true;
}
//-----------------------------------------------------------------------------
// Function      : ESBuilder::generateGraphs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
bool ESBuilder::generateGraphs( const Epetra_CrsGraph & BaseFullGraph )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("ESBuilder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of base graphs
  BaseFullGraph_ = rcp(new Epetra_CrsGraph( BaseFullGraph ));

  int numBlocks = numSamples_;
  blockPattern_.clear();
  blockPattern_.resize(numBlocks);
  for (int i=0;i<numBlocks;++i)
  {
    blockPattern_[i].resize(1);
    blockPattern_[i][0] = i; 
  }

  blockGraph_ = Linear::createBlockGraph(offset_, blockPattern_, *ESMap_, *BaseFullGraph_ );

  return true;
}

} // namespace Linear
} // namespace Xyce

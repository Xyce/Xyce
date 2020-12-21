//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : Builder for ES specific linear objects
//
// Special Notes  : 
//                  This version produces a block structure in which the 
//                  outer loop is parameters, and the inner loop is the original 
//                  (non-block) circuit matrix structure.
//
// Creator       : Eric Keiter, SNL
// Creation Date : 5/31/2018
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_ESBuilder.h>

#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMultiVector.h>
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
Vector * ESBuilder::createVector() const
{
  return Xyce::Linear::createBlockVector( numSamples_, ESMap_, BaseMap_ );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createMultiVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
MultiVector * ESBuilder::createMultiVector( int numVectors ) const
{
  return Xyce::Linear::createBlockMultiVector( numSamples_, numVectors, ESMap_, BaseMap_ );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/14/2018
//-----------------------------------------------------------------------------
Matrix * ESBuilder::createMatrix() const
{
  return Xyce::Linear::createBlockMatrix( numSamples_, offset_, blockPattern_, blockGraph_.get(), baseFullGraph_.get() );
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/31/2018
//-----------------------------------------------------------------------------
Vector * ESBuilder::createStateVector() const
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
Vector * ESBuilder::createStoreVector() const
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
Vector * ESBuilder::createLeadCurrentVector() const
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
bool ESBuilder::generateMaps( const RCP<Parallel::ParMap>& BaseMap, 
                              const RCP<Parallel::ParMap>& oBaseMap )
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
bool ESBuilder::generateStateMaps( const RCP<Parallel::ParMap>& BaseStateMap )
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
bool ESBuilder::generateStoreMaps( const RCP<Parallel::ParMap>& BaseStoreMap )
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
bool ESBuilder::generateLeadCurrentMaps( const RCP<Parallel::ParMap>& BaseLeadCurrentMap )
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
bool ESBuilder::generateGraphs( const Graph& baseFullGraph )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("ESBuilder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of base graphs
  baseFullGraph_ = rcp( baseFullGraph.cloneCopy() );

  int numBlocks = numSamples_;
  blockPattern_.clear();
  blockPattern_.resize(numBlocks);
  for (int i=0;i<numBlocks;++i)
  {
    blockPattern_[i].resize(1);
    blockPattern_[i][0] = i; 
  }

  blockGraph_ = Linear::createBlockGraph(offset_, blockPattern_, *ESMap_, *baseFullGraph_ );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::createSolnColoring
// Purpose       : Color Map representing variable types in solution vector
// Special Notes : This is the block version
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 03/08/04
//-----------------------------------------------------------------------------
const std::vector<int> & ESBuilder::createSolnColoring() const
{
  if (solnColoring_.empty())
  {
    const std::vector<char> & charColors = lasQueryUtil_->rowList_VarType();

    int size = charColors.size();
    solnColoring_.resize( size * numSamples_ );
    for( int i = 0; i < size; ++i )
    {
      switch( charColors[i] )
      {
        case 'V': 
        {
          for (int j=0; j<numSamples_; ++j)
            solnColoring_[j*size + i] = 0;
          break;
        }
        case 'I': 
        {
          for (int j=0; j<numSamples_; ++j)
            solnColoring_[j*size + i] = 1;
          break;
        }
        default : 
        {
          for (int j=0; j<numSamples_; ++j)
            solnColoring_[j*size + i] = 2;
          break;
        }
      }
    }
  }

  return solnColoring_;
}


//-----------------------------------------------------------------------------
// Function      : ESBuilder::createInitialConditionColoring
// Purpose       :
// Special Notes : The .IC and .NODESET capabilities will use the variables
//                 which are colored 0.  This will be all voltage nodes not
//                 connected to independent sources.
// Scope         : Public
// Creator       : Eric R. Keiter,  SNL
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
const std::vector<int> & ESBuilder::createInitialConditionColoring() const
{
  if (icColoring_.empty())
  {
    const std::vector<char> & charColors = lasQueryUtil_->rowList_VarType();
    const std::vector<int> & vsrcGIDColors = lasQueryUtil_->vsrcGIDVec();

    int size = charColors.size();
    icColoring_.resize( size * numSamples_ );
    for( int i = 0; i < size; ++i )
    {
      switch( charColors[i] )
      {
        case 'V': 
        {
          for (int j=0; j<numSamples_; ++j)
            icColoring_[j*size + i] = 0;
          break;
        }
        case 'I': 
        {
          for (int j=0; j<numSamples_; ++j)
            icColoring_[j*size + i] = 1;
          break;
        }
        default : 
        {
          for (int j=0; j<numSamples_; ++j)
            icColoring_[j*size + i] = 2;
          break;
        }
      }
    }  

    int vsrcSize = vsrcGIDColors.size();
    for( int i=0; i < vsrcSize; ++i )
    {  
      int vsrcID = vsrcGIDColors[i];
      // Convert the ID from local to global if it is valid and the build is parallel.
      if (vsrcID >= 0)
      {
        if (!pdsMgr_->getPDSComm()->isSerial())
          vsrcID = BaseMap_->globalToLocalIndex( vsrcGIDColors[i] );

        if (vsrcID < size && vsrcID >= 0)
        {
          for (int j=0; j<numSamples_; ++j)
            icColoring_[j*size + vsrcID] = 1;
        }
      }
    }
  }

  return icColoring_;
}

//-----------------------------------------------------------------------------
// Function      : ESBuilder::vnodeGIDVec()
// Purpose       :
// Special Notes : This is overridden for blockAnalysis types (like MPDE & HB)
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/11/2019
//-----------------------------------------------------------------------------
const std::vector<int> & ESBuilder::vnodeGIDVec() const
{
  if (vnodeVec_.empty())
  {
    const std::vector<int>& vnodeTmp = lasQueryUtil_->vnodeGIDVec();
    int size = vnodeTmp.size();
    vnodeVec_.resize( size * numSamples_ );
   
    for( int i = 0; i < size; ++i )
    {
      for (int j = 0; j < numSamples_; ++j)
      {
        // The GID should be a multiple of the offset plus the current GID. 
        vnodeVec_[j*size + i] = j*offset_ + vnodeTmp[i];  
      }
    } 
  }
  return vnodeVec_;

}

} // namespace Linear
} // namespace Xyce

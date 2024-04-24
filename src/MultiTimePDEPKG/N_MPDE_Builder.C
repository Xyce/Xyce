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
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_MPDE_Builder.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_Manager.h>

#include <N_LAS_Graph.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_SystemHelpers.h>
#include <N_PDS_Comm.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>

using Xyce::DEBUG_MPDE;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
Xyce::Linear::Vector * N_MPDE_Builder::createVector() const
{
  if (warpMPDE_)
  {
    // tscoffe/tmei 08/11/05:  Appending an extra row for omega and phi
    return Xyce::Linear::createBlockVector( Size_, MPDEMap_, BaseMap_, 2 );
  }
  else
  {
    return Xyce::Linear::createBlockVector( Size_, MPDEMap_, BaseMap_ );
  }
}



//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
Xyce::Linear::Vector * N_MPDE_Builder::createStateVector() const
{
  return Xyce::Linear::createBlockVector( Size_, MPDEStateMap_, BaseStateMap_ );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
Xyce::Linear::Vector * N_MPDE_Builder::createStoreVector() const
{
  return Xyce::Linear::createBlockVector( Size_, MPDEStoreMap_, BaseStoreMap_ );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createLeadCurrentVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
Xyce::Linear::Vector * N_MPDE_Builder::createLeadCurrentVector() const
{
  return dynamic_cast<Xyce::Linear::Vector*>(
        Xyce::Linear::createBlockVector( Size_, MPDELeadCurrentMap_, BaseLeadCurrentMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::createMatrix
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
Xyce::Linear::Matrix * N_MPDE_Builder::createMatrix() const
{
  std::vector< std::vector<int> > Cols(Size_);
  int Start = Disc_.Start();
  int Width = Disc_.Width();

  for( int i = 0; i < Size_; ++i )
  {
    Cols[i].resize(Width);
    for( int j = 0; j < Width; ++j )
    {
      int Loc = i+(j+Start);
      if( Loc < 0 )           Loc += Size_;
      else if( Loc >= Size_ ) Loc -= Size_;
      Cols[i][j] = Loc;
    }
    sort(Cols[i].begin(),Cols[i].end());
  }

  if (warpMPDE_)
  {
    // tscoffe/tmei 08/11/05:  Appending an extra row & column for omega and phi
    return Xyce::Linear::createBlockMatrix( Size_,
                                            offset_, 
                                            Cols,
                                            MPDEGraph_.get(),
                                            BaseGraph_.get(),
                                            2 );
  }
  else
  {
    return Xyce::Linear::createBlockMatrix( Size_,
                                            offset_,
                                            Cols,
                                            MPDEGraph_.get(),
                                            BaseGraph_.get() );
  }
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateMaps( const RCP<Xyce::Parallel::ParMap>& BaseMap )
{
  //Save copy of base map
  BaseMap_ = BaseMap;
  
  //determine block offset
  offset_ = Xyce::Linear::generateOffset( *BaseMap );

  std::vector<int> augGIDs(2);
  if (warpMPDE_)
  {
    // tscoffe/tmei 08/02/05:  Added two to size of map for omega and phi
    MPDEMap_ = Xyce::Linear::createBlockParMap( Size_, *BaseMap, 2, &augGIDs );
    omegaGID_ = augGIDs[0];
    phiGID_ = augGIDs[1];
    
    // Figure out which processor owns the augmented rows.
    int omegaLID = MPDEMap_->globalToLocalIndex( omegaGID_ );
    int augProc = -1;
    if (omegaLID >= 0)
      augProc = BaseMap->pdsComm().procID();
    BaseMap->pdsComm().maxAll( &augProc, &augProcID_, 1 );
  }
  else
  {
    MPDEMap_ = Xyce::Linear::createBlockParMap( Size_, *BaseMap );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateStateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateStateMaps( const RCP<Xyce::Parallel::ParMap>& BaseStateMap )
{
  //Save copy of base map
  BaseStateMap_ = BaseStateMap;

  //determine block offset
  stateOffset_= Xyce::Linear::generateOffset( *BaseStateMap );

  //Setup Block Maps for state
  MPDEStateMap_ = Xyce::Linear::createBlockParMap( Size_, *BaseStateMap );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateStoreMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateStoreMaps( const RCP<Xyce::Parallel::ParMap>& BaseStoreMap )
{
  //Save copy of base map
  BaseStoreMap_ = BaseStoreMap;

  //determine block offset
  storeOffset_= Xyce::Linear::generateOffset( *BaseStoreMap );

  //Setup Block Maps for store
  MPDEStoreMap_ = Xyce::Linear::createBlockParMap( Size_, *BaseStoreMap );

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateLeadCurrentMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateLeadCurrentMaps( const RCP<Xyce::Parallel::ParMap>& BaseLeadCurrentMap )
{
  //Save copy of base map
  BaseLeadCurrentMap_ = BaseLeadCurrentMap;

  //determine block offset
  leadCurrentOffset_ = Xyce::Linear::generateOffset( *BaseLeadCurrentMap );

  //Setup Block Maps for store
  MPDELeadCurrentMap_ = Xyce::Linear::createBlockParMap( Size_, *BaseLeadCurrentMap );

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::generateGraphs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Builder::generateGraphs( const Xyce::Linear::Graph & BaseGraph )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("N_MPDE_Builder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of base graphs
  BaseGraph_ = rcp( BaseGraph.cloneCopy() );

  int BlockSize = BaseMap_->numLocalEntities();
  int MaxIndices = BaseGraph_->maxNumIndices();
  int numLocalRows = MPDEMap_->numLocalEntities();
  std::vector<int> Indices(MaxIndices);
  std::vector<int> nnzs(numLocalRows, 0);
  std::vector<std::vector<int> > allIndices(numLocalRows);
  int DiscStart = Disc_.Start();
  int DiscWidth = Disc_.Width();
  std::vector<int> Cols(DiscWidth);
  int NumIndices;
  int BaseRow;
  int MPDERow, localMPDERow;
  for( int i = 0; i < Size_; ++i )
  {
    for( int j = 0; j < DiscWidth; ++j )
    {
      Cols[j] = i + (j+DiscStart);
      if( Cols[j] < 0 ) Cols[j] += Size_;
      else if( Cols[j] > (Size_-1) ) Cols[j] -= Size_;
    }

    for( int j = 0; j < BlockSize; ++j )
    {
      BaseRow = BaseMap_->localToGlobalIndex(j);
      BaseGraph_->extractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      MPDERow = BaseRow + offset_*i;
      localMPDERow = MPDEMap_->globalToLocalIndex( MPDERow );

      for( int k = 0; k < DiscWidth; ++k )
      {
        int Shift = Cols[k]*offset_;
        for( int kk = 0; kk < NumIndices; ++kk )
          allIndices[localMPDERow].push_back( Indices[kk] + Shift );
      }

      //Diagonal Block
      for( int k = 0; k < NumIndices; ++k )
        allIndices[localMPDERow].push_back( Indices[k] + offset_*i );
    }
  }

  if (warpMPDE_)
  {
    // tscoffe 01/15/07 This block adds dependence on omega in the dFdx matrix everywhere that q is nonzero.
    // dqdt1 + omega dqdt2 + f + b
    // q = dqdt1   f = omega dqdt2 + f   b = b
    for( int i = 0; i < Size_; ++i )
    {
      for( int j = 0; j < BlockSize; ++j )
      {
        BaseRow = BaseMap_->localToGlobalIndex(j);
        MPDERow = BaseRow + offset_*i;
        localMPDERow = MPDEMap_->globalToLocalIndex( MPDERow );
        allIndices[localMPDERow].push_back( omegaGID_ );
      }
    }

    if ( BaseMap_->pdsComm().procID() == augProcID_ )
    {
      // Enter graph details for phi:  \dot{phi(t_1)} = omega(t_1)
      localMPDERow = MPDEMap_->globalToLocalIndex( phiGID_ );
      allIndices[localMPDERow].push_back( phiGID_ );
      allIndices[localMPDERow].push_back( omegaGID_ );

      Teuchos::RCP<std::vector<int> > phaseGraph = warpMPDEPhasePtr_->getPhaseGraph();
      localMPDERow = MPDEMap_->globalToLocalIndex( omegaGID_ );
      allIndices[localMPDERow].insert( allIndices[localMPDERow].end(), phaseGraph->begin(), phaseGraph->end() );

      // An (omegaGID, omegaGID) entry must be inserted if not done by the phase graph.
      // This is because when the augmented column is loaded, an entry is expected for each
      // row of column omegaGID.
      bool isOmegaCol = false;
      for (int i=0; i<phaseGraph->size(); ++i)
      {
        if ((*phaseGraph)[i] == omegaGID_)
          isOmegaCol = true;
      }

      // Add (omegaGID, omegaGID) entry if one doesn't already exist.
      if (!isOmegaCol) 
      {
        localMPDERow = MPDEMap_->globalToLocalIndex( omegaGID_ );
        allIndices[localMPDERow].push_back( omegaGID_ );
      }
    }
  }
  for (int i=0; i<numLocalRows; ++i)
  {
    std::sort( allIndices[i].begin(), allIndices[i].end() );
    allIndices[i].erase( std::unique( allIndices[i].begin(), allIndices[i].end() ), allIndices[i].end() );
    nnzs[i] = allIndices[i].size();
  }

  MPDEGraph_ = rcp( Xyce::Linear::createGraph( *MPDEMap_, *MPDEMap_, nnzs, allIndices ) );
  
  if (DEBUG_MPDE && Xyce::isActive(Xyce::Diag::MPDE_PARAMETERS))
  {  
    Xyce::dout() << "Final MPDEGraph = " << std::endl;
    MPDEGraph_->print(Xyce::dout());
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::getSolutionMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
Teuchos::RCP<const Xyce::Parallel::ParMap> N_MPDE_Builder::getSolutionMap() const
{
  return( MPDEMap_ );
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Builder::getSolutionMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
Teuchos::RCP<Xyce::Parallel::ParMap> N_MPDE_Builder::getSolutionMap() 
{
  return( MPDEMap_ );
}


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
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_HBBuilder.h>

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
// Function      : HBBuilder::HBBuilder
// Purpose       : Default constructor 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
HBBuilder::HBBuilder( const int Size, const bool hbOsc)
: numHarmonics_(Size),
  numSolVariables_(0),
  numStateVariables_(0),
  numStoreVariables_(0),
  offset_(0),
  stateOffset_(0),
  storeOffset_(0),
  hbOsc_(hbOsc)
{
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
Vector * HBBuilder::createVector( double initialValue ) const
{
  RCP<Vector> vector =
    createExpandedRealFormTransposeBlockVector(); 
  vector.release(); // Release ownership of the object.
  return(&*vector);
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createTimeDomainBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createTimeDomainBlockVector( ) const
{
  RCP<BlockVector> vector = rcp(
        new BlockVector( numHarmonics_, HBMap_, BaseMap_ ) 
        );
  return(vector);
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createTimeDomainStateBlockVector( ) 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Ting Mei 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createTimeDomainStateBlockVector( ) const
{
  RCP<BlockVector> vector = rcp(
        new BlockVector( numHarmonics_, HBStateMap_, BaseStateMap_ ) 
        );
  return(vector);
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createTimeDomainStoreBlockVector( ) 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createTimeDomainStoreBlockVector( ) const
{
  RCP<BlockVector> vector = rcp(
        new BlockVector( numHarmonics_, HBStoreMap_, BaseStoreMap_ ) 
        );
  return(vector);
}


//-----------------------------------------------------------------------------
// Function      : HBBuilder::createTimeDomainLeadCurrentBlockVector( ) 
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createTimeDomainLeadCurrentBlockVector() const
{
  RCP<BlockVector> vector = rcp(
        new BlockVector( numHarmonics_, HBLeadCurrentMap_, BaseLeadCurrentMap_) 
        );
  return(vector);
}


//-----------------------------------------------------------------------------
// Function      : HBBuilder::createExpandedRealFormBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Rich Schiek 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createExpandedRealFormBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( 2*numHarmonics_, HBExpandedRealFormBVMap_, BaseMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createExpandedRealFormTransposeBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Rich Schiek 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createExpandedRealFormTransposeBlockVector() const
{

  if(hbOsc_)
  {
    RCP<BlockVector> vec = rcp(
        new BlockVector( 2*numHarmonics_, HBExpandedRealFormBVMap_, 1)
        );
    return(vec);
  }
  else
  {
    RCP<BlockVector> vec = rcp(
        new BlockVector( 2*numHarmonics_, HBExpandedRealFormBVMap_ )
        );
    return(vec);
  }

}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createExpandedRealFormTransposeStateBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey 1414, Rich Schiek 1437
// Creation Date : 9/9/08
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createExpandedRealFormTransposeStateBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( 2*numHarmonics_, HBExpandedRealFormStateBVMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createExpandedRealFormTransposeStoreBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createExpandedRealFormTransposeStoreBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( 2*numHarmonics_, HBExpandedRealFormStoreBVMap_ )
        );
  return(vec);
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createExpandedRealFormTransposeLeadCurrentBlockVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
RCP<BlockVector> HBBuilder::createExpandedRealFormTransposeLeadCurrentBlockVector() const
{
  RCP<BlockVector> vec = rcp(
        new BlockVector( 2*numHarmonics_, HBExpandedRealFormLeadCurrentBVMap_ )
        );
  return(vec);
}


//-----------------------------------------------------------------------------
// Function      : HBBuilder::createStateVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
Vector * HBBuilder::createStateVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numHarmonics_, HBStateMap_, BaseStateMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createStoreVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
Vector * HBBuilder::createStoreVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numHarmonics_, HBStoreMap_, BaseStoreMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::createLeadCurrentVector
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
Vector * HBBuilder::createLeadCurrentVector( double initialValue ) const
{
  return dynamic_cast<Vector*>(
        new BlockVector( numHarmonics_, HBLeadCurrentMap_, BaseLeadCurrentMap_ ) );
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::generateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool HBBuilder::generateMaps( const RCP<N_PDS_ParMap>& BaseMap, 
                                    const RCP<N_PDS_ParMap>& oBaseMap )
{
  //Save copy of base map
  BaseMap_ = BaseMap;
  oBaseMap_ = oBaseMap;
  //Xyce::dout() << "BaseMap" << std::endl;
  //BaseMap->petraMap()->Print(std::cout);
  //Xyce::dout() << "oBaseMap" << std::endl;
  //oBaseMap->petraMap()->Print(std::cout);

  //Determine block offset
  offset_ = generateOffset( *BaseMap );

  //Check to see if a graph can be generated using 32-bit integers for the 
  //specified number of harmonics.
  int MaxGID = BaseMap->maxGlobalEntity();
  long int hbMaxGID = MaxGID + offset_*(numHarmonics_-1);
  if ( hbMaxGID > Teuchos::OrdinalTraits<int>::max() )
  {
      int allowableHarmonics = (Teuchos::OrdinalTraits<int>::max() - MaxGID + offset_)/offset_ - 1;
      std::string msg = "HBBuilder::generateMaps():  Harmonic Balance map cannot be constructed, reduce number of harmonics to " 
                      + Teuchos::Utils::toString( allowableHarmonics ) + ".";
      Report::UserError0() <<  msg;
  }                          

  // Use the block linear system helper to create the block parallel maps
  std::vector<RCP<N_PDS_ParMap> > blockMaps = createBlockParMaps(numHarmonics_, *BaseMap, *oBaseMap);

  HBMap_ = blockMaps[0];
  oHBMap_ = blockMaps[1];
  
  //Xyce::dout() << "HBMap_" << std::endl;
  //HBMap_->petraMap()->Print(std::cout);
  //Xyce::dout() << "oHBMap_" << std::endl;
  //oHBMap_->petraMap()->Print(std::cout); 

  int indexBase = oBaseMap_->indexBase();

  // Helpful names for various sizes (subtract indexBase to remove ground node, if there is one):
  numSolVariables_ = oBaseMap_->numLocalEntities() + indexBase;

  if (hbOsc_)
  {
    HBExpandedRealFormBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseMap, 1, &augmentedLIDs_);
    HBExpandedRealFormBVOverlapMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseMap, *oBaseMap, 1, &augmentedOverlapLIDs_);
  }
  else
  {
    HBExpandedRealFormBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseMap );
    HBExpandedRealFormBVOverlapMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseMap, *oBaseMap );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::generateStateMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
bool HBBuilder::generateStateMaps( const RCP<N_PDS_ParMap>& BaseStateMap )
{
  //Save copy of base map
  BaseStateMap_ = BaseStateMap;

  //determine block offset
  stateOffset_ = generateOffset( *BaseStateMap );

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  HBStateMap_ = createBlockParMap(numHarmonics_, *BaseStateMap);

  // Helpful names for various sizes:
  numStateVariables_ = BaseStateMap_->numGlobalEntities();

  // Expanded Real Form for Block Vector map:
  HBExpandedRealFormStateBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseStateMap );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::generateStoreMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool HBBuilder::generateStoreMaps( const RCP<N_PDS_ParMap>& BaseStoreMap )
{
  //Save copy of base map
  BaseStoreMap_ = BaseStoreMap;

  //determine block offset
  storeOffset_ = generateOffset( *BaseStoreMap );

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  HBStoreMap_ = createBlockParMap(numHarmonics_, *BaseStoreMap);

  // Helpful names for various sizes:
  numStoreVariables_ = BaseStoreMap_->numGlobalEntities();

  // Expanded Real Form for Block Vector map:
  HBExpandedRealFormStoreBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseStoreMap );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::generateLeadCurrentMaps
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
bool HBBuilder::generateLeadCurrentMaps( const RCP<N_PDS_ParMap>& BaseLeadCurrentMap )
{
  //Save copy of base map
  BaseLeadCurrentMap_ = BaseLeadCurrentMap;

  //determine block offset
  leadCurrentOffset_ = generateOffset( *BaseLeadCurrentMap );

  // Use the block linear system helper to create the block parallel maps
  // NOTE:  At this time augmented parallel maps are not supported.
  HBLeadCurrentMap_ = createBlockParMap(numHarmonics_, *BaseLeadCurrentMap);

  // Helpful names for various sizes:
  numLeadCurrentVariables_ = BaseLeadCurrentMap_->numGlobalEntities();

  // Expanded Real Form for Block Vector map:
  HBExpandedRealFormLeadCurrentBVMap_ = createBlockFreqERFParMap( numHarmonics_, *BaseLeadCurrentMap );

  return true;
}
//-----------------------------------------------------------------------------
// Function      : HBBuilder::generateGraphs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool HBBuilder::generateGraphs( const Epetra_CrsGraph & BaseFullGraph )
{
  if( Teuchos::is_null(BaseMap_) )
    Xyce::Report::DevelFatal0().in("HBBuilder::generateGraphs")
      << "Need to setup Maps first";

  //Copies of base graphs
  BaseFullGraph_ = rcp(new Epetra_CrsGraph( BaseFullGraph ));

  return true;
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::getSolutionMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<const N_PDS_ParMap> HBBuilder::getSolutionMap() const
{
    return( HBExpandedRealFormBVMap_ );
}

//-----------------------------------------------------------------------------
// Function      : HBBuilder::getSolutionMap
// Purpose       : 
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<N_PDS_ParMap> HBBuilder::getSolutionMap() 
{
    return( HBExpandedRealFormBVMap_ );
}

} // namespace Linear
} // namespace Xyce

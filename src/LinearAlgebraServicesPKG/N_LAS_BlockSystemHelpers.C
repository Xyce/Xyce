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
//
// Purpose        : This is collection of non-member functions that help
//                  in the construction of block linear systems, like those
//                  found in AC analysis.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 06/22/11
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_BlockSystemHelpers.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParHelpers.h>
#include <N_PDS_EpetraParMap.h>

#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>

#include <N_LAS_Graph.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_BlockMatrix.h>

#include <N_LAS_System.h>
#include <N_LAS_QueryUtil.h>

#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : generateOffset
// Purpose       : A helper function that standardizes how offsets are computed 
// Special Notes : Block maps, graphs, vectors, and matrices require a global
//               : numbering scheme.  Xyce uses an offset index to space the 
//               : global ids apart so that they are unique.  The computation
//               : of this offset is performed by this function using the
//               : N_PDS_ParMap from the base block object.
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 10/8/13
//-----------------------------------------------------------------------------
int generateOffset( const N_PDS_ParMap& baseMap )
{
   // Compute the offset needed for global indexing.
   int offset = baseMap.maxGlobalEntity();
   if (baseMap.indexBase() == 0)
     offset++;

   if (offset <= 0)
     offset = 1;

   return offset;
}

//-----------------------------------------------------------------------------
// Function      : copyToBlockVector
// Purpose       : A helper function that copies the array of Vectors
//               : into an BlockVector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyToBlockVector( std::vector<Teuchos::RCP<Vector> >& inputVectors, BlockVector& blockVector )
{
  int inputVecLength = inputVectors.size();
  int blockVecLength = blockVector.blockCount();

  // If the number of blocks is not the same, throw an error.
  if (inputVecLength != blockVecLength) {}
   
  // Loop over the vectors and copy each one.  
  for (int i=0 ; i<blockVecLength ; ++i)
  {
    blockVector.block(i) = *(inputVectors[i]);
  }
}

//-----------------------------------------------------------------------------
// Function      : copyFromBlockVector
// Purpose       : A helper function that copies a BlockVector to an 
//               : array of Vectors.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyFromBlockVector( BlockVector& blockVector, std::vector<Teuchos::RCP<Vector> >& outputVectors )
{
  int outputVecLength = outputVectors.size();
  int blockVecLength = blockVector.blockCount();

  // If the number of blocks is not the same, throw an error.
  if (outputVecLength != blockVecLength) {}
 
  // Loop over the vectors and copy each one. 
  for (int i=0 ; i<blockVecLength ; ++i)
  {
    // Copy over the vector and then import any off-processor contributions.
    *(outputVectors[i]) = blockVector.block(i);
    outputVectors[i]->importOverlap();
  }
}

//-----------------------------------------------------------------------------
// Function      : copyFromBlockVector
// Purpose       : A helper function that copies a set of BlockVectors to a
//               : MultiVector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyFromBlockVectors( std::vector<Teuchos::RCP<BlockVector> >& blockVectors,
                           MultiVector& outputMultiVector )
{
  int blockVecNumVecs = blockVectors.size();
  int blockVecLength = blockVectors[0]->blockCount();
  int outputNumVecs = outputMultiVector.numVectors();

  // If the number of blocks is not the same, throw an error.
  if (outputNumVecs != blockVecLength*blockVecNumVecs) {}
 
  // Loop over the vectors and copy each one. 
  for (int j=0 ; j<blockVecNumVecs; ++j)
  {
    for (int i=0 ; i<blockVecLength ; ++i)
    {
      Teuchos::RCP<Vector> mvCol = 
        Teuchos::rcp( outputMultiVector.getNonConstVectorView(j*blockVecLength+i) );

      // Copy over the vector and then import any off-processor contributions.
      *(mvCol) = blockVectors[j]->block(i);
    }
  }

  outputMultiVector.importOverlap();
}

//-----------------------------------------------------------------------------
// Function      : createBlockParMaps
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns both the map and overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//----------------------------------------------------------------------------- 
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap )
{
   // The information about the current maps, given by pmap and omap (overlap map)
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int olocalBlockSize = omap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();
   int oglobalBlockSize = omap.numGlobalEntities();

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();
   int oBaseIndex = omap.indexBase();

   // Compute the offset needed for global indexing.
   int offset = generateOffset( pmap );

   // Determine size of block maps
   int numGlobalElements = numBlocks*globalBlockSize;
   int onumGlobalElements = numBlocks*oglobalBlockSize;
   int numLocalElements = numBlocks*localBlockSize;
   int onumLocalElements = numBlocks*olocalBlockSize;

   // Initialize vectors to hold the GIDs for the original and block map.
   std::vector<int> BaseGIDs(localBlockSize), oBaseGIDs(olocalBlockSize);
   std::vector<int> GIDs(numLocalElements), oGIDs(onumLocalElements);

   // Extract the global indices.
   N_PDS_EpetraParMap& e_pmap = dynamic_cast<N_PDS_EpetraParMap&>(pmap);
   e_pmap.petraMap()->MyGlobalElements( &BaseGIDs[0] );
   N_PDS_EpetraParMap& e_omap = dynamic_cast<N_PDS_EpetraParMap&>(omap);
   e_omap.petraMap()->MyGlobalElements( &oBaseGIDs[0] );
   
   int gnd_node = 0;  // Will be decremented before first use.

   for( int i = 0; i < numBlocks; ++i )
   {
     // Setting up GIDs for the map without overlap and ground nodes 
     for( int j = 0; j < localBlockSize; ++j )
     {
       GIDs[i*localBlockSize+j] = BaseGIDs[j] + offset*i;
     }

     // Setting up GIDs for the map with overlap and ground nodes
     for( int j = 0; j < (olocalBlockSize+oBaseIndex); ++j )
     {
       oGIDs[i*olocalBlockSize+j] = oBaseGIDs[j] + offset*i;
     }

     if (oBaseIndex == -1)
     {
       // The last GID is the ground node (-1)
       gnd_node--;
       oGIDs[(i+1)*olocalBlockSize-1] = gnd_node;
     }
   }

   // Adapt base index for unique ground node numbering.
   oBaseIndex = std::min( oBaseIndex, gnd_node );

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap, oBlockMap;
   blockMap = Teuchos::rcp(Parallel::createPDSParMap(numGlobalElements, numLocalElements, GIDs, BaseIndex, pmap.pdsComm()));
   oBlockMap = Teuchos::rcp(Parallel::createPDSParMap(onumGlobalElements, onumLocalElements, oGIDs, oBaseIndex, pmap.pdsComm()));

   std::vector<Teuchos::RCP<N_PDS_ParMap> > allMaps;
   allMaps.push_back(blockMap);
   allMaps.push_back(oBlockMap);

   return allMaps;
}

//-----------------------------------------------------------------------------
// Function      : createBlockParMaps
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns both the map and overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//----------------------------------------------------------------------------- 
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps2( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap )
{
   // The information about the current maps, given by pmap and omap (overlap map)
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int olocalBlockSize = omap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();
   int oglobalBlockSize = omap.numGlobalEntities();
   int overlapSize = (olocalBlockSize - 1) - localBlockSize;

   // If overlapSize == -1, then the overlap map is the same size as the parallel map.

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();
   int oBaseIndex = omap.indexBase();

   // Compute the offset needed for global indexing.
   int offset = generateOffset( pmap );

   // Determine size of block maps
   int numProcs = pmap.pdsComm().numProc();
   int numGlobalElements = numBlocks*globalBlockSize;
   int onumGlobalElements = numBlocks*(oglobalBlockSize-numProcs) + numProcs;
   int numLocalElements = numBlocks*localBlockSize;
   int onumLocalElements = numBlocks*(olocalBlockSize-1) + 1;

   // Initialize vectors to hold the GIDs for the original and block map.
   std::vector<int> BaseGIDs(localBlockSize), oBaseGIDs(olocalBlockSize);
   std::vector<int> GIDs(numLocalElements), oGIDs(onumLocalElements);

   // Extract the global indices.
   N_PDS_EpetraParMap& e_pmap = dynamic_cast<N_PDS_EpetraParMap&>(pmap);
   e_pmap.petraMap()->MyGlobalElements( &BaseGIDs[0] );
   N_PDS_EpetraParMap& e_omap = dynamic_cast<N_PDS_EpetraParMap&>(omap);
   e_omap.petraMap()->MyGlobalElements( &oBaseGIDs[0] );
   
   for( int i = 0; i < numBlocks; ++i )
   {
     // Setting up GIDs for the map without overlap and ground nodes 
     for( int j = 0; j < localBlockSize; ++j )
     {
       GIDs[i*localBlockSize+j] = BaseGIDs[j] + offset*i;
       
       // Load all the local elements first in the block map.
       // All the external nodes should be inserted at the end of the GID list.
       oGIDs[i*localBlockSize+j] = oBaseGIDs[j] + offset*i;
     }

     // Now insert the external (overlap) nodes at the end of the GID list.
     for ( int j = localBlockSize, jj=0; j < olocalBlockSize-1; ++j, ++jj )
     {
       oGIDs[numLocalElements+ i*overlapSize + jj] = oBaseGIDs[j] + offset*i;
     }
   }
   // Insert ground node.
   oGIDs[onumLocalElements-1] = -1;

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap, oBlockMap;
   blockMap = Teuchos::rcp(Parallel::createPDSParMap(numGlobalElements, numLocalElements, GIDs, BaseIndex, pmap.pdsComm()));
   oBlockMap = Teuchos::rcp(Parallel::createPDSParMap(onumGlobalElements, onumLocalElements, oGIDs, oBaseIndex, pmap.pdsComm()));

   std::vector<Teuchos::RCP<N_PDS_ParMap> > allMaps;
   allMaps.push_back(blockMap);
   allMaps.push_back(oBlockMap);

   return allMaps;
}


//-----------------------------------------------------------------------------
// Function      : createBlockParMap
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns only the map and not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//----------------------------------------------------------------------------- 
Teuchos::RCP<N_PDS_ParMap> createBlockParMap( int numBlocks, N_PDS_ParMap& pmap, 
                                              int augmentRows, std::vector<int>* augmentedGIDs,
                                              int offset )
{
   // The information about the current maps, given by pmap 
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();

   // Compute the offset needed for global indexing.
   if (offset < 0)
     offset = generateOffset( pmap );

   // Determine size of block maps
   int numGlobalElements = numBlocks*globalBlockSize + augmentRows;
   int numLocalElements = numBlocks*localBlockSize;

   // Find which processor has the maximum global ID.
   // NOTE:  In some cases the "last processor" does not have the largest, or any, IDs assigned to it.
   int maxGID = pmap.maxGlobalEntity();
   int maxProc = -1;
   if ( pmap.globalToLocalIndex( maxGID ) >= 0 ) 
     maxProc = pmap.pdsComm().procID();
 
   // Add the augmented rows to the final processor (assume there aren't too many of these rows)
   if (augmentRows && maxProc >= 0)
   {
     numLocalElements += augmentRows;
   }

   // Initialize vectors to hold the GIDs for the original and block map.
   std::vector<int> BaseGIDs(localBlockSize);
   std::vector<int> GIDs(numLocalElements);

   // Extract the global indices.
   N_PDS_EpetraParMap& e_pmap = dynamic_cast<N_PDS_EpetraParMap&>(pmap);
   e_pmap.petraMap()->MyGlobalElements( &BaseGIDs[0] );

   for( int i = 0; i < numBlocks; ++i )
   {
     // Setting up GIDs for the map without overlap and ground nodes 
     for( int j = 0; j < localBlockSize; ++j )
     {
       GIDs[i*localBlockSize+j] = BaseGIDs[j] + offset*i;
     }
   }

   // Add the augmented rows to the final processor (assume there aren't too many of these rows)
   // All processors will be returned the augmented GIDs.  They will have to determine via the map or by
   // checking if they are the last processor in the communicator if they own the GID.
   if (augmentRows && augmentedGIDs)
   {
     std::vector<int> tmpAugGIDs( augmentRows, -1 );
     augmentedGIDs->resize( augmentRows );

     // Add the augmented GIDs to the processor that has the maximum GID.
     if ( maxProc >= 0 )
     {
       for ( int i=0; i<augmentRows; i++ )
       {
         GIDs[numLocalElements-augmentRows+i] = GIDs[numLocalElements-augmentRows-1+i] + 1;
         tmpAugGIDs[i] = GIDs[numLocalElements-augmentRows+i];
       }
     }

     // Now communicate the GIDs to all the processors.  Users of this map will have to check if they own the
     // GID or just check if they are the last processor to whom the GID was assigned.  
     pmap.pdsComm().maxAll( &tmpAugGIDs[0], &(*augmentedGIDs)[0], augmentRows );
   }

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap = 
     Teuchos::rcp(Parallel::createPDSParMap(numGlobalElements, numLocalElements, GIDs, BaseIndex, pmap.pdsComm()));

   return blockMap;
}

//-----------------------------------------------------------------------------
// Function      : createBlockGraph
// Purpose       : A helper function for creating block parallel graphs.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<Graph> createBlockGraph( int offset, std::vector<std::vector<int> >& blockPattern, 
                                      N_PDS_ParMap& blockMap, const Graph& baseGraph )
{
  Teuchos::RCP<const Epetra_CrsGraph> epetraGraph = baseGraph.epetraObj();

  int numBlockRows = blockPattern.size();
  int numMyBaseRows = epetraGraph->NumMyRows();
  int maxIndices = epetraGraph->MaxNumIndices();
 
  int maxBlockCols = blockPattern[0].size();
  for (int i=1; i<numBlockRows; ++i) { 
    int cols=blockPattern[i].size();
    if (cols > maxBlockCols)
      maxBlockCols = cols;
  }
 
  //Construct block graph based on  [All graphs are the same, so only one needs to be made]
  N_PDS_EpetraParMap& e_blockMap = dynamic_cast<N_PDS_EpetraParMap&>(blockMap);
 
  Teuchos::RCP<Epetra_CrsGraph> newEpetraGraph = rcp(new Epetra_CrsGraph( Copy, *dynamic_cast<Epetra_BlockMap*>(e_blockMap.petraMap()), 0 ));
  
  std::vector<int> indices(maxIndices);
  int shift=0, index=0, baseRow=0, blockRow=0, numIndices=0;
  int maxNNZs = maxIndices*maxBlockCols;
  std::vector<int> newIndices(maxNNZs);  // Make as large as the combined maximum of indices and column blocks

  for( int j = 0; j < numMyBaseRows; ++j )
  {
    // Extract the base entries from the base row.
    baseRow = epetraGraph->GRID(j);
    epetraGraph->ExtractGlobalRowCopy( baseRow, maxIndices, numIndices, &indices[0] );

    for( int i = 0; i < numBlockRows; ++i )
    {
      // For this harmonic, which row will be inserted.
      blockRow = baseRow + offset*i;

      int numBlockCols = blockPattern[i].size();

      // Find all entries from a row before inserting it.
      for( int k = 0; k < numBlockCols; ++k )
      {
        // Find which block column to start at.
        shift = blockPattern[i][k]*offset;  // Actual column index.
        index = k*numIndices;  // Pointer to next block of column indices.
        for( int kk = 0; kk < numIndices; ++kk ) newIndices[index+kk] = indices[kk] + shift;
      }

      // Insert entire row for all blocks.
      newEpetraGraph->InsertGlobalIndices( blockRow, numBlockCols*numIndices, &newIndices[0] );
    }
  }
  newEpetraGraph->FillComplete();
  newEpetraGraph->OptimizeStorage();

  Teuchos::RCP<Graph> newGraph = Teuchos::rcp( new Graph( newEpetraGraph ) );
 
  return newGraph;
}
//-----------------------------------------------------------------------------
// Function      : createBlockFreqERFParMap
// Purpose       : A helper function for creating block parallel maps for
//               : the frequency domain.  The map generated here has all the
//               : harmonics for one time point grouped together in expanded
//               : real form.
//               : This function returns only the map, not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 12/3/13
//-----------------------------------------------------------------------------
Teuchos::RCP<N_PDS_ParMap> createBlockFreqERFParMap( int numHarmonics, N_PDS_ParMap& pmap, 
                                                     int augmentRows, std::vector<int>* augmentedLIDs )
{
   // The information about the current maps, given by pmap
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   int localBlockSize = pmap.numLocalEntities();
   int globalBlockSize = pmap.numGlobalEntities();

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();

   // Determine size of block maps
   int erfNumHarms = 2*numHarmonics;
   int numGlobalElements = erfNumHarms*globalBlockSize;
   int numLocalElements = erfNumHarms*localBlockSize;
   std::vector<int> newGIDs( numLocalElements );

   // Determine GIDs for block maps
   for (int i=0; i<localBlockSize; i++)
   {
     int gid = pmap.localToGlobalIndex( i );

     for (int numharms=0; numharms<erfNumHarms; numharms++)
     {
       newGIDs[i*erfNumHarms+numharms] = gid*erfNumHarms + numharms;
     }
   } 

   // Find which processor has the maximum global ID.
   // NOTE:  In some cases the "last processor" does not have the largest, or any, IDs assigned to it.
   int maxProc = -1;
   std::vector<int> augmentedGIDs;

   // Add the augmented rows to the final processor (assume there aren't too many of these rows)
   if (augmentRows)
   {
     int maxGID = pmap.maxGlobalEntity();
     if ( pmap.globalToLocalIndex( maxGID ) >= 0 ) 
       maxProc = pmap.pdsComm().procID();
    
     // Increment the local count for the processor with the maximum global ID
     if (maxProc >= 0)
     {
       augmentedGIDs.resize( augmentRows );
       numLocalElements += augmentRows;
       
       for (int i=0; i<augmentRows; i++)
         augmentedGIDs[i] = numGlobalElements+i;

       newGIDs.insert( newGIDs.end(), augmentedGIDs.begin(), augmentedGIDs.end() );
     }
   }

   // Increment the global count.
   numGlobalElements += augmentRows;

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap = 
     Teuchos::rcp(Parallel::createPDSParMap(numGlobalElements, numLocalElements, newGIDs, BaseIndex, pmap.pdsComm()));

   // Now capture the LIDs of the augmented rows.
   if ( maxProc >= 0 )
   {
     augmentedLIDs->resize( augmentRows );

     for ( int i=0; i<augmentRows; i++ )
     {
       int localID = blockMap->globalToLocalIndex( augmentedGIDs[i] );
        
       // Add the augmented LIDs to the processor that has the maximum GID.
       (*augmentedLIDs)[i] = localID;
     }
   }

   return blockMap;
}

//-----------------------------------------------------------------------------
// Function      : createBlockFreqERFParMap
// Purpose       : A helper function for creating block parallel maps for
//               : the frequency domain.  The map generated here has all the
//               : harmonics for one time point grouped together in expanded
//               : real form.
//               : This function returns only the map, not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 12/3/13
//-----------------------------------------------------------------------------
Teuchos::RCP<N_PDS_ParMap> createBlockFreqERFParMap( int numHarmonics, N_PDS_ParMap& pmap, 
                                                     N_PDS_ParMap& omap,
                                                     int augmentRows, std::vector<int>* augmentedLIDs )
{
   // The information about the current maps, given by pmap
   // will be used to generate new parallel maps for the block system.

   // Get the current number of entries owned by this processor and globally.
   // Subtract 1 from the overlap maps for the ground node.
   int oLocalBlockSize = omap.numLocalEntities();
   if ( omap.indexBase() == -1 )
     oLocalBlockSize--;

   // Get the index base from the original maps
   int BaseIndex = pmap.indexBase();

   // Determine size of block maps
   int erfNumHarms = 2*numHarmonics;
   int numGlobalElements = -1;
   int numLocalElements = erfNumHarms*oLocalBlockSize;
   std::vector<int> newGIDs( numLocalElements );

   // Determine GIDs for block maps
   for (int i=0; i<oLocalBlockSize; i++)
   {
     int gid = omap.localToGlobalIndex( i );

     for (int numharms=0; numharms<erfNumHarms; numharms++)
     {
       newGIDs[i*erfNumHarms+numharms] = gid*erfNumHarms + numharms;
     }
   } 

   // Find which processor has the maximum global ID.
   // NOTE:  In some cases the "last processor" does not have the largest, or any, IDs assigned to it.
   int maxProc = -1;
   std::vector<int> augmentedGIDs;
   int maxGlobalElement = erfNumHarms*pmap.numGlobalEntities();

   // Add the augmented rows to the final processor (assume there aren't too many of these rows)
   if (augmentRows)
   {
     int maxGID = pmap.maxGlobalEntity();
     if ( pmap.globalToLocalIndex( maxGID ) >= 0 ) 
       maxProc = pmap.pdsComm().procID();
 
     // Increment the local count for the processor with the maximum global ID
     if (maxProc >= 0)
     {
       augmentedGIDs.resize( augmentRows );
       numLocalElements += augmentRows;
       
       for (int i=0; i<augmentRows; i++)
         augmentedGIDs[i] = maxGlobalElement+i;

       newGIDs.insert( newGIDs.end(), augmentedGIDs.begin(), augmentedGIDs.end() );
     }
   }

   // Create new maps for the block system 
   Teuchos::RCP<N_PDS_ParMap> blockMap = 
     Teuchos::rcp(Parallel::createPDSParMap(numGlobalElements, numLocalElements, newGIDs, BaseIndex, pmap.pdsComm()));

   // Now capture the GIDs of the augmented rows.
   if (augmentRows)
   {
     if ( maxProc >= 0 )
     {
       augmentedLIDs->resize( augmentRows );

       for ( int i=0; i<augmentRows; i++ )
       {
         int localID = blockMap->globalToLocalIndex( augmentedGIDs[i] );
        
         // Add the augmented LIDs to the processor that has the maximum GID.
         (*augmentedLIDs)[i] = localID;
       }
     }
   }

   return blockMap;
}


//-----------------------------------------------------------------------------
// Function      : computePermutedDFT
// Purpose       : A helper function for applying the DFT to a block vector.
// Special Notes : xf = D*P*xt, xf has the same block format as (P*xt).
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 3/28/14
//-----------------------------------------------------------------------------
void computePermutedDFT(N_UTL_DFTInterfaceDecl<std::vector<double> > & dft, 
                        const BlockVector & xt, BlockVector * xf,
                        const std::vector<int>* lids)
{
  int blockCount = xt.blockCount();
  int localN = xt.block(0).localLength();

  // It's necessary to get the blockmap from Epetra because the N_PDS_ParMap is not always guaranteed to be valid.
  Epetra_BlockMap blockMap = xt.block(0).epetraObj().Map();

  // Obtain registered vectors with the DFT interface.
  Teuchos::RCP< std::vector<double> > inputSignal, outputSignal;
  dft.getDFTVectors( inputSignal, outputSignal );

  int dftScalar = dft.getScalar(); 

  // Loop through all the solution variables and compute the DFT of each variable.
  if ( lids )
    localN = lids->size();

  for (int j=0; j<localN; j++)
  {
    int lid = j;
    if ( lids )
      lid = (*lids)[j];

    // Get the global id for this variable.
    int gid = blockMap.GID(lid);

    Vector& freqVecRef = xf->block(gid);

    for (int i=0; i<blockCount; ++i)
    {
      Vector& timeVecRef = xt.block(i);
      (*inputSignal)[i] = timeVecRef[lid];
    }

    // Calculate the DFT for the inputSignal.
    dft.calculateDFT();

    freqVecRef[0] =  (*outputSignal)[0]/dftScalar;
    freqVecRef[1] =  (*outputSignal)[1]/dftScalar;

    for (int i=1; i<(blockCount+1)/2; ++i)
    {
      freqVecRef[2*i] =  (*outputSignal)[2*i]/dftScalar;
      freqVecRef[2*(blockCount-i)] =  (*outputSignal)[2*i]/dftScalar;

      freqVecRef[2*i+1] =  (*outputSignal)[2*i+1]/dftScalar;
      freqVecRef[2*(blockCount-i)+1] = -(*outputSignal)[2*i+1]/dftScalar;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : computePermutedDFT2
// Purpose       : A helper function for applying the DFT to a block vector.
// Special Notes : xf = D*P*xt, xf has the same block format as (P*xt).
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 3/28/14
//-----------------------------------------------------------------------------
void computePermutedDFT2(N_UTL_DFTInterfaceDecl<std::vector<double> > & dft, 
                        const BlockVector & xt, BlockVector * xf)
{
  int blockCount = xt.blockCount();
  int localN = xt.block(0).localLength();

  // It's necessary to get the blockmap from Epetra because the N_PDS_ParMap is not always guaranteed to be valid.
  Epetra_BlockMap blockMap = xt.block(0).epetraObj().Map();

  // Obtain registered vectors with the DFT interface.
  Teuchos::RCP< std::vector<double> > inputSignal, outputSignal;
  dft.getDFTVectors( inputSignal, outputSignal );

  int dftScalar = dft.getScalar(); 

  for (int j=0; j<localN; j++)
  {
    // Get the global id for this variable.
    int gid = blockMap.GID(j);

    Vector& freqVecRef = xf->block(gid);

    double norm1 = 0.0;
    for (int i=0; i<blockCount; ++i)
    {
      Vector& timeVecRef = xt.block(i);
      double val = timeVecRef[j];
      (*inputSignal)[i] = val;
      if (std::abs<double>(val) > norm1)
        norm1 = std::abs<double>(val);
    }

    // Calculate the DFT for the inputSignal.
    if (norm1 > 0.0)
    {
      dft.calculateDFT();

      freqVecRef[0] =  (*outputSignal)[0]/dftScalar;
      freqVecRef[1] =  (*outputSignal)[1]/dftScalar;

      for (int i=1; i<(blockCount+1)/2; ++i)
      {
        freqVecRef[2*i] =  (*outputSignal)[2*i]/dftScalar;
        freqVecRef[2*(blockCount-i)] =  (*outputSignal)[2*i]/dftScalar;

        freqVecRef[2*i+1] =  (*outputSignal)[2*i+1]/dftScalar;
        freqVecRef[2*(blockCount-i)+1] = -(*outputSignal)[2*i+1]/dftScalar;
      }
    }
  }
}
//----------------------------------------------------------------------------- 
// Function      : computePermutedIFT
// Purpose       : A helper function for applying the IFT to a block vector.
// Special Notes : xt = P^{-1}D^{-1}*xf, xf has the same block format as (P*xt).
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 3/28/14
//-----------------------------------------------------------------------------
void computePermutedIFT(N_UTL_DFTInterfaceDecl<std::vector<double> > & dft, 
                        const BlockVector & xf, BlockVector * xt,
                        int numTimePts_, const std::vector<int>* lids)
{
  int N = xf.block(0).globalLength();
  int localBC = xt->block(0).localLength();

  // It's necessary to get the blockmap from Epetra because the N_PDS_ParMap is not always guaranteed to be valid.
  Epetra_BlockMap blockMap = (xt->block(0)).epetraObj().Map();

  // Register input and output signals for the IFT.
  Teuchos::RCP< std::vector<double> > inputSignal, outputSignal;
  dft.getIFTVectors( inputSignal, outputSignal );

  if ( numTimePts_ == 0 )
    numTimePts_ = N/2;

  if (numTimePts_ != N/2 && (inputSignal->size() == (numTimePts_ +1)) ) 
    inputSignal->assign( numTimePts_ +1, 0.0); 

  int dftScalar = dft.getScalar();

  if ( lids )
    localBC = lids->size();
   
  for (int j=0; j<localBC; j++)
  {
    int lid = j;
    if ( lids )
      lid = (*lids)[j];

    // Get the global id for this variable.
    int gid = blockMap.GID(lid);

    Vector& freqVecRef = xf.block(gid);

    for (int i=0; i<(N/2+1); ++i)
    {
      (*inputSignal)[i] = freqVecRef[i];
    }

    // Calculate the inverse FFT.
    dft.calculateIFT();

    for (int i=0; i<numTimePts_;  ++i)
    {
      Vector& timeVecRef = xt->block(i);
      timeVecRef[j] =  (*outputSignal)[i]*dftScalar;
    }
  }
}

} // namespace Linear
} // namespace Xyce

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
//                  found in AC or HB analysis.
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

#ifndef  Xyce_LAS_BLOCKSYSTEMHELPERS_H
#define  Xyce_LAS_BLOCKSYSTEMHELPERS_H

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>
#include <N_UTL_DFTInterfaceDecl.hpp>

#include <Teuchos_RCP.hpp>

// ---------- Forward Declarations ----------

class N_PDS_ParMap;

namespace Xyce {
namespace Linear {

BlockVector* createBlockVector( int numBlocks,
                                const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                                const Teuchos::RCP<N_PDS_ParMap> & subBlockMap,
                                int augmentRows = 0 );

BlockMultiVector* createBlockMultiVector( int numBlocks, int numVectors,
                                          const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                                          const Teuchos::RCP<N_PDS_ParMap> & subBlockMap );

BlockMatrix* createBlockMatrix( int size,
                                int offset,
                                const std::vector< std::vector<int> > & blockColumns,
                                const Graph* globalGraph,
                                const Graph* subBlockGraph,
                                int augmentCount = 0 );

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
int generateOffset( const N_PDS_ParMap& baseMap );

//-----------------------------------------------------------------------------
// Function      : createBlockParMaps
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns both the map and overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap );
std::vector<Teuchos::RCP<N_PDS_ParMap> > createBlockParMaps2( int numBlocks, N_PDS_ParMap& pmap, N_PDS_ParMap& omap );

//-----------------------------------------------------------------------------
// Function      : createBlockParMap
// Purpose       : A helper function for creating block parallel maps.
//               : This function returns only the map, not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<N_PDS_ParMap> createBlockParMap( int numBlocks, N_PDS_ParMap& pmap,
                                              int augmentRows=0, std::vector<int>* augmentedGIDs = 0,
                                              int offset = -1 );

//-----------------------------------------------------------------------------
// Function      : createBlockGraph
// Purpose       : A helper function for creating block parallel graphs.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<Graph> createBlockGraph( int offset, std::vector<std::vector<int> >& blockPattern,
                                      N_PDS_ParMap& blockMap, const Graph& baseGraph );

//-----------------------------------------------------------------------------
// Function      : createBlockFreqERFParMap
// Purpose       : A helper function for creating block parallel maps for 
//               : the frequency domain.  The map generated here has all the
//               : harmonics for one time point grouped together in expanded
//               : real form.
//               : This function returns only the map, not the overlap map.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
Teuchos::RCP<N_PDS_ParMap> createBlockFreqERFParMap( int numHarmonics, N_PDS_ParMap& pmap,
                                                     int augmentRows = 0, std::vector<int>* augmentedLIDs = 0 );

Teuchos::RCP<N_PDS_ParMap> createBlockFreqERFParMap( int numHarmonics, N_PDS_ParMap& pmap, N_PDS_ParMap& omap,
                                                     int augmentRows = 0, std::vector<int>* augmentedLIDs = 0 );
//-----------------------------------------------------------------------------
// Function      : copyToBlockVector
// Purpose       : A helper function that copies the array of Vectors
//               : into a BlockVector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyToBlockVector( std::vector<Teuchos::RCP<Vector> >& inputVectors, BlockVector& blockVector ); 

//-----------------------------------------------------------------------------
// Function      : copyFromBlockVector
// Purpose       : A helper function that copies a BlockVector to an 
//               : array of Vectors.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyFromBlockVector( BlockVector& blockVector, std::vector<Teuchos::RCP<Vector> >& outputVectors ); 

//-----------------------------------------------------------------------------
// Function      : copyFromBlockVector
// Purpose       : A helper function that copies a BlockVector to a 
//               : MultiVector.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 6/22/11
//-----------------------------------------------------------------------------
void copyFromBlockVectors( std::vector<Teuchos::RCP<BlockVector> >& blockVectors, 
                           MultiVector& outputMultiVector );

//-----------------------------------------------------------------------------
// Function      : computePermutedDFT
// Purpose       : A helper function for applying the DFT to a block vector.
// Special Notes : xf = D*P*xt, xf has the same block format as (P*xt).
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 3/28/14
//-----------------------------------------------------------------------------
void computePermutedDFT(N_UTL_DFTInterfaceDecl<std::vector<double> > & dft, 
                        const BlockVector & xt, BlockVector * xf,
                        const std::vector<int>* lids = 0);

//-----------------------------------------------------------------------------
// Function      : computePermutedDFT2
// Purpose       : A helper function for applying the DFT to a block vector.
// Special Notes : xf = D*P*xt, xf has the same block format as (P*xt).
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 3/28/14
//-----------------------------------------------------------------------------
void computePermutedDFT2(N_UTL_DFTInterfaceDecl<std::vector<double> > & dft, 
                        const BlockVector & xt, BlockVector * xf);

//----------------------------------------------------------------------------- 
// Function      : computePermutedIFT
// Purpose       : A helper function for applying the IFT to a block vector.
// Special Notes : xt = P^{-1}D^{-1}*xf, xf has the same block format as (P*xt).
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 3/28/14
//-----------------------------------------------------------------------------
void computePermutedIFT(N_UTL_DFTInterfaceDecl<std::vector<double> > & dft, 
                        const BlockVector & xf, BlockVector * xt,
                        int numTimePts_ =0, const std::vector<int>* lids = 0);
} // namespace Linear
} // namespace Xyce

#endif

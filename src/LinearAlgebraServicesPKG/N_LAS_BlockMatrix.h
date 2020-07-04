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
// Purpose        : Specification file for the Abstract interface to sparse
//                  block matrix type.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/12/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_BlockMatrix_h
#define Xyce_N_LAS_BlockMatrix_h

// ---------- Standard Includes ----------
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Graph.h>

// ----------  Other Includes   ----------
#include <Teuchos_RCP.hpp>
using Teuchos::RCP;

class Epetra_CrsGraph;
class Epetra_Map;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : BlockMatrix
// Purpose       : Abstract interface to sparse block matrix type.
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
class BlockMatrix : public Matrix
{
 public:

  BlockMatrix( int size,
                     int offset,
                     const std::vector< std::vector<int> > & blockColumns,
                     const Graph* globalGraph,
                     const Graph* subBlockGraph,
                     int augmentCount = 0 );

  //Destructor
  ~BlockMatrix() {}

  //Block Access
  Matrix & block( int row, int col );

  int blockSize()
  { return blockSize_; }
  
  int numBlockRows()
  { return numBlockRows_; }

  // Put function for the block sparse-matrix.
  void put(double s);

  // Replace the entries of an augmented row using the row GID.
  void replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices); 

  void replaceAugmentedColumn(int augmentedColumn, const BlockVector & vec);

  // Assemble global matrix with blocks
  // NOTE:  The global matrix is not always a view of the local matrix, so this function ensures
  // that the values are sync'ed up.  Call this before using the global matrix for computations.
  void assembleGlobalMatrix();
 
  void fillComplete();
 
  void printPetraObject(std::ostream &os) const;

 private:

  bool blocksViewGlobalMat_;
  const int blockSize_;
  const int offset_;
  const int numBlockRows_;
  const int augmentCount_;

  std::vector<int> augmentGIDs_, baseNumCols_, baseIndices_;
  const std::vector< std::vector<int> > cols_;
  std::vector< std::vector<Teuchos::RCP<Matrix> > > blocks_;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_BlockMatrix_h

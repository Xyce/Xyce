//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_LAS_Matrix.h>

// ----------  Other Includes   ----------

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

  BlockMatrix() {}

  //Destructor
  virtual ~BlockMatrix() {}

  //Block Access
  virtual Matrix & block( int row, int col ) = 0;
  virtual const Matrix & block( int row, int col ) const = 0;

  virtual int blockSize() const = 0;
  
  virtual int numBlockRows() const = 0;

  // Replace the entries of an augmented row using the row GID.
  virtual void replaceAugmentedRow(int rowGID, int length, double * coeffs, int * colIndices) = 0; 

  virtual void replaceAugmentedColumn(int augmentedColumn, const BlockVector & vec) = 0;

  // Assemble global matrix with blocks
  // NOTE:  The global matrix is not always a view of the local matrix, so this function ensures
  // that the values are sync'ed up.  Call this before using the global matrix for computations.
  virtual void assembleGlobalMatrix() = 0;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_BlockMatrix_h

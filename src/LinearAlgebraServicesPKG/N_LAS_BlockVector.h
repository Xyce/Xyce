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
//
// Purpose        : Block Vector access
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Computational Sciences
//
// Creation Date  : 3/12/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_BlockVector_h
#define Xyce_N_LAS_BlockVector_h

#include <N_LAS_Vector.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : BlockVector
// Purpose       : Provides an abstract interface for block vectors
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 3/12/04
//-----------------------------------------------------------------------------
class BlockVector : public Vector
{
 public:

  // Default constructor
  BlockVector() {}

  // Destructor
  virtual ~BlockVector() {}

  // Assignment operator
  virtual BlockVector & operator=(const BlockVector & right) = 0;
  virtual BlockVector & operator=(const Vector & right) = 0;

  // Block accessors
  virtual Vector & block( int Loc ) const = 0;

  virtual int blockSize() const = 0;

  virtual int blockCount() const = 0;

  virtual int startBlock() const = 0;

  virtual int endBlock() const = 0;

  // Get the ParMap objects for each BLOCK in this block vector.
  virtual const Parallel::ParMap * blockPmap() const = 0;

};

} // namespace Linear
} // namespace Xyce

#endif


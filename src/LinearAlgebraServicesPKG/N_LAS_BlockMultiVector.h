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
// Purpose        : Block MultiVector access
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

#ifndef Xyce_N_LAS_BlockMultiVector_h
#define Xyce_N_LAS_BlockMultiVector_h

#include <vector>

#include <N_LAS_MultiVector.h>

#include <Teuchos_RCP.hpp>

class N_PDS_ParMap;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : BlockMultiVector
// Purpose       : Provides an abstract interface for block vectors
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 3/12/04
//-----------------------------------------------------------------------------
class BlockMultiVector : public MultiVector
{
 public:
  BlockMultiVector( int numBlocks, int numVectors,
                    const Teuchos::RCP<N_PDS_ParMap> & globalMap,
                    const Teuchos::RCP<N_PDS_ParMap> & subBlockMap
                  );

  //Copy constructor
  BlockMultiVector( const BlockMultiVector & right );

  // Destructor
  virtual ~BlockMultiVector() {};

  // Block accessors
  MultiVector & block( int Loc ) const
  { return *blocks_[Loc]; }

  int blockSize() const
  { return globalBlockSize_; }

  int blockCount() const
  { return numBlocks_; }

  int startBlock() const
  { return startBlock_; }

  int endBlock() const
  { return endBlock_; }

  // Return whether the local vector is a view of the global vector.
  bool isBlockView()
  { return blocksViewGlobalVec_; }

  // Assemble global vector with blocks
  // NOTE:  The global vector is not always a view of the local vector, so this function ensures
  // that the values are sync'ed up.  Call this before using the global vector for computations.
  void assembleGlobalVector();

  // Get the ParMap objects for each BLOCK in this block vector.
  N_PDS_ParMap * blockPmap() { return newBlockMap_.get(); }
  const N_PDS_ParMap * blockPmap() const { return newBlockMap_.get(); }

  // Print out the underlying data in this object.
  void print(std::ostream &os) const;

 private:

  bool blocksViewGlobalVec_;
  const int globalBlockSize_;
  const int localBlockSize_;
  const int numBlocks_;

  // In frequency domain, whole blocks may be owned by one processor.
  // NOTE:  This assumes they are contiguous.  By default these routines
  //        will return 0 and numBlocks_ (which is sane for the time domain specs).
  int startBlock_, endBlock_;

  Teuchos::RCP<N_PDS_ParMap> newBlockMap_; 

  std::vector<Teuchos::RCP<MultiVector> > blocks_;

};

} // namespace Linear
} // namespace Xyce

#endif


//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
#include <N_PDS_fwd.h>

#include <Teuchos_RCP.hpp>

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
                    const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                    const Teuchos::RCP<const Parallel::ParMap> & subBlockMap
                  );

  // Destructor
  virtual ~BlockMultiVector() {};

  // Block accessors
  virtual MultiVector & block( int Loc ) const
  { return *blocks_[Loc]; }

  virtual int blockSize() const
  { return globalBlockSize_; }

  virtual int blockCount() const
  { return numBlocks_; }

  virtual int startBlock() const
  { return startBlock_; }

  virtual int endBlock() const
  { return endBlock_; }

  // Return whether the local vector is a view of the global vector.
  virtual bool isBlockView()
  { return blocksViewGlobalVec_; }

  // Get the ParMap objects for each BLOCK in this block vector.
  virtual const Parallel::ParMap * blockPmap() const { return newBlockMap_.get(); }

  // Print out the underlying data in this object.
  virtual void print(std::ostream &os) const;

 private:

  bool blocksViewGlobalVec_;
  const int globalBlockSize_;
  const int localBlockSize_;
  const int numBlocks_;

  // In frequency domain, whole blocks may be owned by one processor.
  // NOTE:  This assumes they are contiguous.  By default these routines
  //        will return 0 and numBlocks_ (which is sane for the time domain specs).
  int startBlock_, endBlock_;

  Teuchos::RCP<const Parallel::ParMap> newBlockMap_; 

  std::vector<Teuchos::RCP<MultiVector> > blocks_;

};

} // namespace Linear
} // namespace Xyce

#endif


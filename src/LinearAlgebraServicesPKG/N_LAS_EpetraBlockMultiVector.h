//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Block MultiVector implementation for Epetra
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Computational Sciences
//
// Creation Date  : 3/12/04
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_EpetraBlockMultiVector_h
#define Xyce_N_LAS_EpetraBlockMultiVector_h

#include <vector>

#include <N_LAS_BlockMultiVector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_EpetraHelpers.h>

#include <N_PDS_fwd.h>

#include <Teuchos_RCP.hpp>
#include <Epetra_MultiVector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : BlockMultiVector
// Purpose       : Provides an abstract interface for block vectors
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 3/12/04
//-----------------------------------------------------------------------------
class EpetraBlockMultiVector : public BlockMultiVector, public EpetraVectorAccess 
{
 public:

  EpetraBlockMultiVector( int numBlocks, int numVectors,
                          const Teuchos::RCP<const Parallel::ParMap> & globalMap,
                          const Teuchos::RCP<const Parallel::ParMap> & subBlockMap
                        );

  // Destructor
  virtual ~EpetraBlockMultiVector();

  // Assignment operator
  BlockMultiVector& operator=(const BlockMultiVector& right);
  BlockMultiVector& operator=(const MultiVector& right);

  // Index operator
  double * operator() (int row_lid, int col_lid)
  {
    if (row_lid >= 0 && col_lid >= 0)
      return ((*aMultiVector_)[col_lid]+row_lid);
    else
      return &groundNode_;
  }

  // Index operator
  const double * operator() (int row_lid, int col_lid) const
  {
    if (row_lid >= 0 && col_lid >= 0)
      return ((*aMultiVector_)[col_lid]+row_lid);
    else
      return &groundNode_;
  }

  // Clone operation:
  MultiVector* clone() const;

  // Clone operation:
  MultiVector* cloneCopy() const;

  // Returns the dot product of "this" vector and another.
  // NOTE:  If *this or y is a single vector, it will return the dot product of each column.
  //        Otherwise, *this and y should have the same number of columns and d[i] = (*this)[i]'*y[i].
  void dotProduct(const MultiVector & y, std::vector<double>& d) const;

  // Scale every entry in the multi-vector by "a"
  void scale(const double a);

  // Matrix-Matrix multiplication.  this[i] = this[i]*x[i] for each vector
  void multiply(const MultiVector & x);

  // Standard blas AXPY operations
  void update(double a, const MultiVector & A, double s = 1.0);

  void update(double a, const MultiVector & A, double b,
              const MultiVector & B, double s = 1.0);

  // Compute the l_p norm (e.g., 2-norm is l_2)
  int lpNorm(const int p, double * result) const;

  // Infinity norm
  int infNorm(double * result, int * index = 0) const;

  // Weighted root-mean-square norm
  int wRMSNorm(const MultiVector & weights, double * result) const;

  // Weighted max-norm (with index if allocated by the caller)
  int wMaxNorm(const MultiVector & weights, double * result, int * index = 0) const;

  // Generate random number
  void random()
  { aMultiVector_->Random(); }

  // Fill vector with constant value.
  void putScalar(const double scalar);

  // Add to vector with constant value.
  void addScalar(const double scalar);

  // Absolute value element-wise for vector
  void absValue(const MultiVector & A);

  // Reciprocal of elements in vector
  void reciprocal(const MultiVector & A);

  // Get the global (across all processors) length of the multi-vector
  int globalLength() const
  { return aMultiVector_->GlobalLength(); }
  // Get the local (on processor component) length of the multi-vector
  int localLength() const
  { return aMultiVector_->MyLength(); }
  // Get the number of individual vectors in the multi-vector
  int numVectors() const
  { return aMultiVector_->NumVectors(); }

  // Get the external vector size
  int externVectorSize() const { return 0; }

  // Import/Export capability
  bool vectorImport(const MultiVector * vec, Importer * importer);
  bool importOverlap() { return true; }

  // Vector access function
  const Vector* getVectorView(int index) const;
  Vector* getNonConstVectorView(int index);
  const Vector* getVectorViewAssembled(int index) const;
  Vector* getNonConstVectorViewAssembled(int index);

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

  // Clear the external vector map
  void clearExternVectorMap() {}

  // Add an element to the external vector map
  void addElementToExternVectorMap(const int & global_index,
                                   const double & value) {}

  // Get the parallel map associated with this multi-vector
  const Parallel::ParMap * pmap() const { return parallelMap_; }
  const Parallel::ParMap * omap() const { return parallelMap_; }

  // Get the parallel communicator associated with this multi-vector
  const Parallel::Communicator* pdsComm() const { return pdsComm_.get(); }

  Epetra_MultiVector & epetraObj() { return *aMultiVector_; }
  const Epetra_MultiVector & epetraObj() const { return *aMultiVector_; }
  Epetra_MultiVector & epetraOverlapObj() { return *aMultiVector_; }
  const Epetra_MultiVector & epetraOverlapObj() const { return *aMultiVector_; }

  // Get for vector elements by their global index 
  const double & getElementByGlobalIndex(const int & global_index, const int & vec_index = 0) const;

  // Set for vector elements by their global index
  bool setElementByGlobalIndex(const int & global_index, const double & val,
                                       const int & vec_index = 0);

  // Sum vector elements by their global index
  bool sumElementByGlobalIndex(const int & global_index, const double & val,
                                       const int & vec_index = 0);

  // Accumulate off processor fill contributions if necessary
  void fillComplete() {}

  const Parallel::ParMap * blockPmap() const { return subBlockMap_.get(); }

  // Print out the underlying data in this object.
  void print(std::ostream &os) const;

  // Dump vector entries to file.
  void writeToFile( const char * filename, bool useLIDs=false, bool mmFormat=false ) const 
  { Xyce::Linear::writeToFile( *aMultiVector_, filename, useLIDs, mmFormat ); }

 private:

  // Copy constructor
  EpetraBlockMultiVector(const EpetraBlockMultiVector & right);

  // Pointer to the multi-vector's parallel map object
  const Parallel::ParMap* parallelMap_;

  // Pointer the Petra multi-vector object.
  Epetra_MultiVector * aMultiVector_;

  // Communicator object, if one is needed.
  Teuchos::RCP<const Parallel::Communicator> pdsComm_;

  // isOwned flags
  bool vecOwned_, mapOwned_;

  // Dummy variable for loading ground node contributions.
  double groundNode_;

  bool blocksViewGlobalVec_;
  const int globalBlockSize_;
  const int localBlockSize_;
  const int numBlocks_;

  // In frequency domain, whole blocks may be owned by one processor.
  // NOTE:  This assumes they are contiguous.  By default these routines
  //        will return 0 and numBlocks_ (which is sane for the time domain specs).
  int startBlock_, endBlock_;

  Teuchos::RCP<const Parallel::ParMap> subBlockMap_; 

  std::vector<Teuchos::RCP<MultiVector> > blocks_;

  // Process library error codes.
  void processError(const char *methodMsg, int error) const;

};

} // namespace Linear
} // namespace Xyce

#endif


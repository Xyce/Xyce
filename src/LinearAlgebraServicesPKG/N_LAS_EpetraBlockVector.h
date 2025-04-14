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

#ifndef Xyce_N_LAS_EpetraBlockVector_h
#define Xyce_N_LAS_EpetraBlockVector_h

#include <vector>

#include <N_LAS_BlockVector.h>
#include <N_LAS_Vector.h>
#include <N_PDS_fwd.h>

#include <Teuchos_RCP.hpp>

#include <N_LAS_EpetraHelpers.h>
#include <Epetra_Vector.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : EpetraBlockVector
// Purpose       : Provides an Epetra implementation for block vectors
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 3/12/04
//-----------------------------------------------------------------------------
class EpetraBlockVector : public BlockVector, public EpetraVectorAccess
{
 public:
  EpetraBlockVector( int numBlocks,
               const Teuchos::RCP<const Parallel::ParMap> & globalMap,
               const Teuchos::RCP<const Parallel::ParMap> & subBlockMap,
               int augmentRows = 0 );

  // Constructor that uses the block size to divide up the number of elements on
  // each processor into vectors whose values are all "owned" by one processor.
  // NOTE:  This constructor is handy for frequency-domain representations of time-domain vectors.
  EpetraBlockVector( int blockSize,
               const Teuchos::RCP<const Parallel::ParMap> & globalMap,
               int augmentRows = 0 );

  // View constructor
  //NOTE:  This constructor assumes that the Vector is divided up into blockSize subvectors,
  //       whose values are solely owned by one of the processors.
  EpetraBlockVector( const Vector * right, int blockSize );

  // Destructor
  virtual ~EpetraBlockVector();

  // Assignment operator
  BlockVector & operator=(const BlockVector & right);
  BlockVector & operator=(const Vector & right);

  // Clone operations:
  Vector* cloneVector() const;
  Vector* cloneCopyVector() const;

  MultiVector* clone() const { return cloneVector(); }
  MultiVector* cloneCopy() const { return cloneCopyVector(); }

  // Vector view methods
  const Vector* getVectorView(int index) const;
  Vector* getNonConstVectorView(int index);
  const Vector* getVectorViewAssembled(int index) const { return getVectorView( index ); }
  Vector* getNonConstVectorViewAssembled(int index) { return getNonConstVectorView( index ); }

  // Block accessors
  Vector & block( int Loc ) const
  { return *blocks_[Loc]; }

  int blockSize() const
  { return globalBlockSize_; }

  int blockCount() const
  { return numBlocks_; }

  int startBlock() const
  { return startBlock_; }

  int endBlock() const
  { return endBlock_; }

  // Index operator
  double * operator() (int row_lid, int col_lid);
  const double * operator() (int row_lid, int col_lid) const;

  // Operation: operator []
  double & operator[] (int index);
  const double & operator[] (int index) const;

  // Dot product methods
  double dotProduct( const Vector & y ) const;
  void dotProduct(const MultiVector & y, std::vector<double>& d) const;

  // Scale every entry in the multi-vector by "a"
  void scale(const double a);

  // Matrix-Matrix multiplication.  this[i] = this[i]*x[i] for each vector
  void multiply(const MultiVector & x); 

  // Linear combination with one and two constants and vectors
  void update(double a, const MultiVector & A, double s = 1.0);

  void update(double a, const MultiVector & A, double b,
                      const MultiVector & B, double s = 1.0);

  // Compute the l_p norm (e.g., 2-norm is l_2)
  int lpNorm(const int p, double * result) const;

  // Weighted root-mean-square norm
  int wRMSNorm(const MultiVector & weights, double * result) const;

  // Infinity norm (with index if allocated by the caller)
  int infNorm(double * result, int * index = 0) const;

  // Weighted max-norm (with index if allocated by the caller)
  int wMaxNorm(const MultiVector & weights, double * result, int * index = 0) const;

  // Generate random number
  void random() { aMultiVector_->Random(); }

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
  int numVectors() const { return 1; }
  // Get the external vector size
  int externVectorSize() const { return 0; }

  // Get the ParMap objects for each BLOCK in this block vector.
  const Parallel::ParMap * blockPmap() const 
  { return newBlockMap_.get(); }

  // Dump vector entries to file.
  void writeToFile( const char * filename, bool useLIDs=false, bool mmFormat=false ) const
  { Xyce::Linear::writeToFile( *aMultiVector_, filename, useLIDs, mmFormat ); }

  // Get for vector elements by their global index (const version)
  const double & getElementByGlobalIndex(const int & global_index, const int & vec_index = 0) const;

  // Set for vector elements by their global index
  bool setElementByGlobalIndex(const int & global_index, const double & val,
                               const int & vec_index = 0);

  // Sum vector elements by their global index
  bool sumElementByGlobalIndex(const int & global_index, const double & val,
                               const int & vec_index = 0);

  // Import/Export capability
  bool vectorImport(const MultiVector * vec, Importer * importer)
  { return true; }
  bool importOverlap() 
  { return true; }

  // Clear the external vector map
  void clearExternVectorMap() {}
  // Add an element to the external vector map
  void addElementToExternVectorMap(const int & global_index, const double & value) {}

  // Get the parallel map associated with this multi-vector
  const Parallel::ParMap * pmap() const { return parallelMap_; } 
  const Parallel::ParMap * omap() const { return parallelMap_; }

  // Get the parallel communicator associated with this multi-vector
  const Parallel::Communicator* pdsComm() const { return pdsComm_.get(); }

  Epetra_MultiVector & epetraObj() { return *aMultiVector_; }
  const Epetra_MultiVector & epetraObj() const { return *aMultiVector_; }

  Epetra_MultiVector & epetraOverlapObj() { return *aMultiVector_; }
  const Epetra_MultiVector & epetraOverlapObj() const { return *aMultiVector_; }

  //Accumulate off processor fill contributions if necessary
  void fillComplete() {}

  // Print out the underlying data in this object.
  void print(std::ostream &os) const;

 private:

  // Pointer to the multi-vector's parallel map object
  const Parallel::ParMap* parallelMap_;

  // Pointer the Petra multi-vector object.
  Epetra_MultiVector * aMultiVector_;

  // isOwned flags
  bool vecOwned_, mapOwned_;

  // Communicator object, if one is needed.
  Teuchos::RCP<const Parallel::Communicator> pdsComm_;

  // Dummy variable for loading ground node contributions.
  double groundNode_;

  bool blocksViewGlobalVec_;
  int globalBlockSize_;
  int localBlockSize_;
  int overlapBlockSize_;
  int numBlocks_;
  int augmentCount_;

  // In frequency domain, whole blocks may be owned by one processor.
  // NOTE:  This assumes they are contiguous.  By default these routines
  //        will return 0 and numBlocks_ (which is sane for the time domain specs).
  int startBlock_, endBlock_;

  Teuchos::RCP<const Parallel::ParMap> newBlockMap_;

  std::vector<Teuchos::RCP<Vector> > blocks_;

};

// Index operator
inline double * EpetraBlockVector::operator() (int row_lid, int col_lid)
{
  if (row_lid >= 0 && col_lid >= 0)
    return (*aMultiVector_)[col_lid]+row_lid;
  else
    return &groundNode_;
}

// Index operator
inline const double * EpetraBlockVector::operator() (int row_lid, int col_lid) const
{
  if (row_lid >= 0 && col_lid >= 0)
    return (*aMultiVector_)[col_lid]+row_lid;
  else
    return &groundNode_;
}

// Operation: operator []
inline double & EpetraBlockVector::operator[] (int index)
{
  if (index >= 0)
    return (*aMultiVector_)[0][index];
  else
    return groundNode_;
}

// Operation: operator []
inline const double & EpetraBlockVector::operator[] (int index) const
{
  if (index >= 0)
    return (*aMultiVector_)[0][index];
  else
    return groundNode_;
}

inline void EpetraBlockVector::addScalar(const double scalar)
{
  int length  = localLength();

  for (int j = 0; j < length; ++j)
    (*this)[j] += scalar;
}

} // namespace Linear
} // namespace Xyce

#endif


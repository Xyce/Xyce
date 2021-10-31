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
// Purpose        : Specification file for the Abstract interface to the
//                  multi-vector types (RDP, RSP, CDP or CSP).
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_EpetraVector_h
#define Xyce_N_LAS_EpetraVector_h

// ---------- Standard Includes ----------
#include <string>
#include <map>
#include <vector>

#include <N_PDS_fwd.h>
#include <N_NLS_fwd.h>
#include <N_LAS_fwd.h>

#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_Vector.h>

// --------  Forward Declarations --------

class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Export;
class Epetra_Import;
class Epetra_Map;

namespace EpetraExt {

 class MultiVector_View;

}

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : EpetraVector
// Purpose       : Provides an Epetra interface for the vector type
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class EpetraVector : public Vector, public EpetraVectorAccess
{
public:

  // Constructors to map to Petra constructors.
  EpetraVector( const Parallel::ParMap & map );

  EpetraVector( const Parallel::ParMap & map, const Parallel::ParMap & ol_map );

  // Constructor that wraps an Epetra vector inside a Linear::Vector.
  // This is used in the nonlinear solver and linear solver interface.
  EpetraVector( Epetra_Vector * origV, bool isOwned);

  // Constructor takes the overlap Epetra vector and generates the assembled vector.
  EpetraVector( Epetra_Vector * overlapMV, const Epetra_BlockMap& parMap, bool isOwned = true );

  // Assignment operator
  Vector & operator=(const Vector & right);
  
  //Destructor
  virtual ~EpetraVector();

  // Clone operation:
  Vector* cloneVector() const;

  // Clone operation:
  Vector* cloneCopyVector() const;

  // Clone operation:
  virtual MultiVector* clone() const
  { return cloneVector(); }

  // Clone operation:
  virtual MultiVector* cloneCopy() const
  { return cloneCopyVector(); }

  // Dot product with another vector.
  double dotProduct(const Vector & y) const;

  // Returns the dot product of "this" vector and another.
  // NOTE:  If *this or y is a single vector, it will return the dot product of each column.
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

  // Infinity norm (with index if allocated by the caller)
  int infNorm(double * result, int * index = 0) const;

  // Weighted root-mean-square norm
  int wRMSNorm(const MultiVector & weights, double * result) const;

  // Weighted max-norm (with index if allocated by the caller)
  int wMaxNorm(const MultiVector & weights, double * result, int * index = 0) const;

  // Generate random number
  void random();

  // Fill vector with constant value.
  void putScalar(const double scalar);

  // Add to vector with constant value.
  void addScalar(const double scalar);

  // Absolute value element-wise for vector
  void absValue(const MultiVector & A);

  // Reciprocal of elements in vector
  void reciprocal(const MultiVector & A);

  // Index operator
  double * operator() (int row_lid, int col_lid)
  {
    if (row_lid >= 0 && col_lid >= 0)
      return (epetraOverlapObj()[col_lid]+row_lid);
    else
      return &groundNode_;
  }

  // Index operator
  const double * operator() (int row_lid, int col_lid) const
  {
    if (row_lid >= 0 && col_lid >= 0)
      return ((*oMultiVector_)[col_lid]+row_lid);
    else
      return &groundNode_;
  }

  // Operation: operator []
  virtual double & operator[] (int index)
  {
    if (index >= 0)
      return (*oMultiVector_)[0][index];
    else
      return groundNode_;
  }

  // Operation: operator []
  virtual const double & operator[] (int index) const
  {
    if (index >= 0)
      return (*oMultiVector_)[0][index];
    else
      return groundNode_;
  }

  // Vector access function
  const Vector* getVectorView(int index) const;
  Vector* getNonConstVectorView(int index);
  const Vector* getVectorViewAssembled(int index) const;
  Vector* getNonConstVectorViewAssembled(int index);

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
  int externVectorSize() const 
  { return externVectorMap_.size(); }

  // Import/Export capability
  bool vectorImport(const MultiVector * vec, Importer * importer);
  bool importOverlap();

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

  // Clear the external vector map
  void clearExternVectorMap() { externVectorMap_.clear(); }

  // Add an element to the external vector map
  void addElementToExternVectorMap(const int & global_index,
  	                           const double & value);

  // Print the underlying object.
  void print(std::ostream &os) const;

  // Get the parallel map associated with this multi-vector
  const Parallel::ParMap * pmap() const { return parallelMap_; }
  const Parallel::ParMap * omap() const { return overlapMap_; }

  // Get the parallel communicator associated with this multi-vector
  const Parallel::Communicator* pdsComm() const { return pdsComm_.get(); }

  Epetra_MultiVector & epetraObj() { return *aMultiVector_; }
  const Epetra_MultiVector & epetraObj() const { return *aMultiVector_; }
  Epetra_MultiVector & epetraOverlapObj() { return *oMultiVector_; }
  const Epetra_MultiVector & epetraOverlapObj() const { return *oMultiVector_; }

  //Accumulate off processor fill contributions if necessary
  void fillComplete();

protected:

  // Copy constructor
  EpetraVector(const EpetraVector & right);

  // Pointer to the multi-vector's parallel map object
  const Parallel::ParMap* parallelMap_;

  // Parallel Map for overlapped data
  const Parallel::ParMap* overlapMap_;

  // Pointer the Petra multi-vector object.
  Epetra_MultiVector * aMultiVector_;

  // Overlapped view of multi-vector
  Epetra_MultiVector * oMultiVector_;

  //Used to distribute dependent data to processors needing it for fills
  Epetra_Import * importer_;

  //Used to accumulate off processor partial contributions to fills
  Epetra_Export * exporter_;

  // Subset View Transform
  EpetraExt::MultiVector_View * viewTransform_;

  // Communicator object, if one is needed.
  Teuchos::RCP<const Parallel::Communicator> pdsComm_;

  // isOwned flags
  bool vecOwned_, mapOwned_;

  // Map containing extern elements from migration
  std::map<int,double> externVectorMap_;

  // Dummy variable for loading ground node contributions.
  double groundNode_;

  // Process library error codes.
  void processError(const char *methodMsg, int error) const;

};

//-----------------------------------------------------------------------------
inline void EpetraVector::addElementToExternVectorMap(const int &
	global_index, const double & value)
{
  if( !externVectorMap_.count(global_index) )
   externVectorMap_[global_index] = value;
}

} // namespace Linear
} // namespace Xyce

#endif


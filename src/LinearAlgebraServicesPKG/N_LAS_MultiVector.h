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

#ifndef Xyce_N_LAS_MultiVector_h
#define Xyce_N_LAS_MultiVector_h

// ---------- Standard Includes ----------
#include <string>
#include <map>
#include <vector>

#include <N_PDS_fwd.h>
#include <N_NLS_fwd.h>
#include <N_LAS_fwd.h>

#include <Epetra_MultiVector.h>

// --------  Forward Declarations --------

class N_PDS_ParMap;

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
// Class         : MultiVector
// Purpose       : Provides an interface for the multi-vector type
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class MultiVector
{
public:

  // Constructors to map to Petra constructors.
  MultiVector( N_PDS_ParMap & map,
               int numVectors = 1 );

  MultiVector( N_PDS_ParMap & map,
               N_PDS_ParMap & ol_map,
               int numVectors = 1 );

  // Constructor that wraps an Epetra multivector inside a Linear::MultiVector.
  // This is used in the nonlinear solver and linear solver interface.
  MultiVector( Epetra_MultiVector * origMV, bool isOwned);

  // Constructor takes the oMultiVector and generates the aMultiVector
  MultiVector( Epetra_MultiVector * overlapMV, const Epetra_BlockMap& parMap, bool isOwned = true );
  
  // Assignment operator
  MultiVector & operator=(const MultiVector & right);

  // Copy constructor
  MultiVector(const MultiVector & right);

  //Destructor
  virtual ~MultiVector();

  // Returns the dot product of "this" vector and another.
  // NOTE:  If *this or y is a single vector, it will return the dot product of each column.
  //        Otherwise, *this and y should have the same number of columns and d[i] = (*this)[i]'*y[i].
  void dotProduct(const MultiVector & y, std::vector<double>& d) const;

  // Scale every entry in the multi-vector by "a"
  void scale(const double a);

  // Scale every entry ([i]) in the multi-vector by "a*x[i]"
  void scale(const double a, const MultiVector & x);

  // Matrix-Matrix multiplication.  this[i] = this[i]*x[i] for each vector
  void multiply(const MultiVector & x);

  // Standard blas DAXPY operation
  void daxpy(const MultiVector & y, const double a,
  	const MultiVector & x);

  // Linear combination with two constants and vectors
  void linearCombo(const double a, const MultiVector & x,
  	const double b, const MultiVector & y);

  void update(double a, const MultiVector & A, double s);

  void update(double a, const MultiVector & A, double b,
  	const MultiVector & B, double s);

  // Compute the l_p norm (e.g., 2-norm is l_2)
  int lpNorm(const int p, double * result) const;

  // Infinity norm
  int infNorm(double * result) const;

  // Infinity norm index (the index pointer should be allocated by the caller)
  int infNormIndex(int * index) const;

  // Weighted root-mean-square norm
  int wRMSNorm(const MultiVector & weights, double * result) const;

  // Weighted max-norm
  int wMaxNorm(const MultiVector & weights, double * result) const;

  // maximum value
  int maxValue(double * result) const;

  // minimum value
  int minValue(double * result) const;

  // Add a multiple of a multi-vector (this = a*y).
  void addVec(const double a, const MultiVector & y);

  // Generate random number
  void random();

  // Fill vector with constant value.
  virtual void putScalar(const double scalar);

  // Add to vector with constant value.
  virtual void addScalar(const double scalar);

  // Absolute value element-wise for vector
  void absValue(const MultiVector & A);

  // Reciprocal of elements in vector
  void reciprocal(const MultiVector & A);

  // Index operator
  double *& operator[] (int index)
  {
    return (*oMultiVector_)[index];
  }

  // Index operator
  double * const & operator[] (int index) const
  {
    return (*oMultiVector_)[index];
  }

  // Index operator
  double * operator() (int row_lid, int col_lid)
  {
    if (row_lid >= 0 && col_lid >= 0)
      return ((*oMultiVector_)[col_lid]+row_lid);
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

  // Vector access function
  const Vector* getVectorView(int index) const;
  Vector* getNonConstVectorView(int index);
  const Vector* getVectorViewAssembled(int index) const;
  Vector* getNonConstVectorViewAssembled(int index);

  // Get the global (across all processors) length of the multi-vector
  int globalLength() const;
  // Get the local (on processor component) length of the multi-vector
  int localLength() const;

  // Get the number of individual vectors in the multi-vector
  int numVectors() const;
  // Get the external vector size
  int externVectorSize() const { return externVectorMap_.size(); }

  // Import/Export capability
  bool vectorImport(const MultiVector * vec, Epetra_Import * importer);
  bool importOverlap();

  // Dump vector entries to file.
  virtual void writeToFile( const char * filename, bool useLIDs=false, bool mmFormat=false ) const;

  // Friend of the Matrix and IterativeSolver classes so their
  // member functions can access our private members.
  friend class Matrix;
  friend class FilteredMatrix;
  friend class Nonlinear::DampedNewton;

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
  virtual void print(std::ostream &os) const;

  // Get the parallel map associated with this multi-vector
  N_PDS_ParMap * pmap() { return parallelMap_; }
  N_PDS_ParMap * omap() { return overlapMap_; }
  const N_PDS_ParMap * pmap() const { return parallelMap_; }
  const N_PDS_ParMap * omap() const { return overlapMap_; }

  // Get the parallel communicator associated with this multi-vector
  N_PDS_Comm * pdsComm() { return pdsComm_.get(); }
  const N_PDS_Comm * pdsComm() const { return pdsComm_.get(); }

  Epetra_MultiVector & epetraObj() { return *aMultiVector_; }
  
  const Epetra_MultiVector & epetraObj() const { return *aMultiVector_; }

  //Accumulate off processor fill contributions if necessary
  void fillComplete();

protected:

  // Pointer to the multi-vector's parallel map object
  N_PDS_ParMap* parallelMap_;

  // Parallel Map for overlapped data
  N_PDS_ParMap* overlapMap_;

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
  Teuchos::RCP<N_PDS_Comm> pdsComm_;

  // isOwned flag
  bool isOwned_;

  // Map containing extern elements from migration
  std::map<int,double> externVectorMap_;

  // Process library error codes.
  void processError(const char *methodMsg, int error) const;

  // Dummy variable for loading ground node contributions.
  double groundNode_;

};

//-----------------------------------------------------------------------------
inline void MultiVector::addElementToExternVectorMap(const int &
	global_index, const double & value)
{
  if( !externVectorMap_.count(global_index) )
   externVectorMap_[global_index] = value;
}

} // namespace Linear
} // namespace Xyce

#endif

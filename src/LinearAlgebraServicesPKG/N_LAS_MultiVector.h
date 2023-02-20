//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : MultiVector
// Purpose       : Provides an abstract interface for the multi-vector type
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class MultiVector 
{
public:

  // Default constructor
  MultiVector() {}

  //Destructor
  virtual ~MultiVector() {}

  // Assignment operator
  virtual MultiVector & operator=(const MultiVector & right) = 0;

  // Clone operation:
  virtual MultiVector* clone() const = 0;

  // Clone operation:
  virtual MultiVector* cloneCopy() const = 0;

  // Returns the dot product of "this" vector and another.
  // NOTE:  If *this or y is a single vector, it will return the dot product of each column.
  //        Otherwise, *this and y should have the same number of columns and d[i] = (*this)[i]'*y[i].
  virtual void dotProduct(const MultiVector & y, std::vector<double>& d) const = 0;

  // Scale every entry in the multi-vector by "a"
  virtual void scale(const double a) = 0;

  // Matrix-Matrix multiplication.  this[i] = this[i]*x[i] for each vector
  virtual void multiply(const MultiVector & x) = 0;

  // Linear combination with one and two constants and vectors
  virtual void update(double a, const MultiVector & A, double s = 1.0) = 0;

  virtual void update(double a, const MultiVector & A, double b,
  	              const MultiVector & B, double s = 1.0) = 0;

  // Compute the l_p norm (e.g., 2-norm is l_2)
  virtual int lpNorm(const int p, double * result) const = 0;

  // Infinity norm (with index if allocated by the caller)
  virtual int infNorm(double * result, int * index = 0) const = 0;

  // Weighted root-mean-square norm
  virtual int wRMSNorm(const MultiVector & weights, double * result) const = 0;

  // Weighted max-norm (with index if allocated by the caller)
  virtual int wMaxNorm(const MultiVector & weights, double * result, int * index = 0) const = 0;

  // Generate random number
  virtual void random() = 0;

  // Fill vector with constant value.
  virtual void putScalar(const double scalar) = 0;

  // Add to vector with constant value.
  virtual void addScalar(const double scalar) = 0;

  // Absolute value element-wise for vector
  virtual void absValue(const MultiVector & A) = 0;

  // Reciprocal of elements in vector
  virtual void reciprocal(const MultiVector & A) = 0;

  // Index operator
  virtual double * operator() (int row_lid, int col_lid) = 0;

  // Index operator
  virtual const double * operator() (int row_lid, int col_lid) const = 0;

  // Vector access function
  virtual const Vector* getVectorView(int index) const = 0;
  virtual Vector* getNonConstVectorView(int index) = 0;
  virtual const Vector* getVectorViewAssembled(int index) const = 0;
  virtual Vector* getNonConstVectorViewAssembled(int index) = 0;

  // Get the global (across all processors) length of the multi-vector
  virtual int globalLength() const = 0;
  // Get the local (on processor component) length of the multi-vector
  virtual int localLength() const = 0;

  // Get the number of individual vectors in the multi-vector
  virtual int numVectors() const = 0;
  // Get the external vector size
  virtual int externVectorSize() const = 0;

  // Import/Export capability
  virtual bool vectorImport(const MultiVector * vec, Importer * importer) = 0;
  virtual bool importOverlap() = 0;

  // Dump vector entries to file.
  virtual void writeToFile( const char * filename, bool useLIDs=false, bool mmFormat=false ) const = 0;

  // Get for vector elements by their global index (const version)
  virtual const double & getElementByGlobalIndex(const int & global_index, const int & vec_index = 0) const = 0;

  // Set for vector elements by their global index
  virtual bool setElementByGlobalIndex(const int & global_index, const double & val,
  	                               const int & vec_index = 0) = 0;

  // Sum vector elements by their global index
  virtual bool sumElementByGlobalIndex(const int & global_index, const double & val,
  	                               const int & vec_index = 0) = 0;

  // Clear the external vector map
  virtual void clearExternVectorMap() = 0;
  // Add an element to the external vector map
  virtual void addElementToExternVectorMap(const int & global_index,
  	                                   const double & value) = 0;

  // Print the underlying object.
  virtual void print(std::ostream &os) const = 0;

  // Get the parallel map associated with this multi-vector
  virtual const Parallel::ParMap * pmap() const = 0;
  virtual const Parallel::ParMap * omap() const = 0;

  // Get the parallel communicator associated with this multi-vector
  virtual const Parallel::Communicator* pdsComm() const = 0;

  //Accumulate off processor fill contributions if necessary
  virtual void fillComplete() = 0;
};

} // namespace Linear
} // namespace Xyce

#endif

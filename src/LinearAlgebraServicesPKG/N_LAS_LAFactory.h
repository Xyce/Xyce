//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Specification file for the interface for creating linear
//                  algebra objects using the GoF Abstract Factory design
//                  pattern.  It must be used with a compatible linear algebra
//                  package (e.g., Trilinos/Petra to provide the concrete
//                  implementations.
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

#ifndef Xyce_N_LAS_LAFactory_h
#define Xyce_N_LAS_LAFactory_h

// ---------- Standard Includes ----------
#include <vector>


#include <N_LAS_fwd.h>

class N_PDS_ParMap;

namespace Xyce {
namespace Linear {

typedef double DataType;

//-----------------------------------------------------------------------------
// Class         : LAFactory
// Purpose       : Provides an interface for creating linear algebra objects
//                 using the GoF Abstract Factory design pattern.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/20/00
//-----------------------------------------------------------------------------
class LAFactory
{

 public:

  // Creates a new LAS vector
  static Vector * newVector( DataType type, N_PDS_ParMap & map );
  static Vector * newVector( DataType type,
                                   N_PDS_ParMap & map,
                                   N_PDS_ParMap & ol_map );

  // Creates a new LAS multi-vector
  static MultiVector * newMultiVector( DataType type,
                                             N_PDS_ParMap & map,
                                             int numVectors );
  static MultiVector * newMultiVector( DataType type,
                                             N_PDS_ParMap & map,
                                             int numVectors,
                                             N_PDS_ParMap & ol_map );

  // Creates a new LAS matrix
  static Matrix * newMatrix( DataType type,
                             N_PDS_ParMap & map,
                             std::vector<int> & diagArray);

 private:

  // The Linear Algebra Services default factory
  LAFactory();

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_LAFactory_h

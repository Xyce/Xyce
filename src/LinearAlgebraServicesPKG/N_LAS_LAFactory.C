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

//-------------------------------------------------------------------------
//
// Purpose        : Implementation file for the interface for creating linear
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_LAS_LAFactory.h>

#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : LAFactory::newVector
// Purpose       : Concrete implementation of the "newVector" function for
//                 creating a single vector.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
Vector * LAFactory::newVector( DataType type,
                                           N_PDS_ParMap & map )
{
  Vector * tmp = new Vector( map );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : LAFactory::newVector
// Purpose       : Concrete implementation of the "newVector" function for
//                 creating a single vector.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 10/16/00
//-----------------------------------------------------------------------------
Vector * LAFactory::newVector( DataType type,
                                           N_PDS_ParMap & map,
                                           N_PDS_ParMap & ol_map )
{
  Vector * tmp = new Vector( map, ol_map );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : LAFactory::newMultiVector
// Purpose       : Concrete implementation of the "newMultiVector" function.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
MultiVector * LAFactory::newMultiVector( DataType type,
                                                     N_PDS_ParMap & map,
                                                     int numVectors)
{
  MultiVector * tmp = new MultiVector( map, numVectors );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : LAFactory::newMultiVector
// Purpose       : Concrete implementation of the "newMultiVector" function.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
MultiVector * LAFactory::newMultiVector( DataType type,
                                                     N_PDS_ParMap & map,
                                                     int numVectors,
                                                     N_PDS_ParMap & ol_map )
{
  MultiVector * tmp = new MultiVector( map, ol_map, numVectors );
//  tmp->putScalar( type );
  return tmp;
}

//-----------------------------------------------------------------------------
// Function      : LAFactory::newMatrix
// Purpose       : Concrete implementation of the "newMatrix" function.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/22/00
//-----------------------------------------------------------------------------
Matrix * LAFactory::newMatrix( DataType type,
                                           N_PDS_ParMap & map,
					   std::vector<int> & diagArray)
{
  Matrix * tmp = new Matrix( map, diagArray );
//  tmp->put( type );
  return tmp;
}

} // namespace Linear
} // namespace Xyce

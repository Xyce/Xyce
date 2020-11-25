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
// Purpose        : This is collection of non-member functions that help
//                  in the construction of linear systems.  These methods
//                  were developed for the seperation of systems into variables
//                  related to linear and nonlinear devices.
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 04/06/15
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_SYSTEMHELPERS_H
#define  Xyce_LAS_SYSTEMHELPERS_H

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>
#include <Teuchos_RCP.hpp>

// ---------- Forward Declarations ----------

class N_PDS_ParMap;

namespace Xyce {
namespace Linear {

  // Non-member creation methods.
MultiVector* createMultiVector( N_PDS_ParMap & map,
                                int numVectors = 1 );

MultiVector* createMultiVector( N_PDS_ParMap & map,
                                      N_PDS_ParMap & ol_map,
                                int numVectors = 1 );

Vector* createVector( N_PDS_ParMap & map );

Vector* createVector( N_PDS_ParMap & map, N_PDS_ParMap & ol_map );

Matrix* createMatrix( const Graph* overlapGraph,
                      const Graph* baseGraph );

Graph* createGraph( N_PDS_ParMap & map,
                    const std::vector<int>& numIndicesPerRow );

Graph* createGraph( N_PDS_ParMap & map,
                    int maxNumIndicesPerRow );

//-----------------------------------------------------------------------------
// Function      : extractValues
// Purpose       : A helper function that extracts the values into another matrix
// Special Notes : When a matrix is being separated into various components,
//               : we need to have a way to extract out the subset of the entries
//               : into different matrices.  It is assumed that the matrix has
//               : the same parallel distribution as the original, just a subset
//               : of the entries.
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 04/06/15
//-----------------------------------------------------------------------------
void extractValues( const Matrix& inputMatrix, 
                    std::vector< Teuchos::RCP<Matrix> >& outputMatrices );

} // namespace Linear
} // namespace Xyce

#endif

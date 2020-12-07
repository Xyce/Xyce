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
//                  in the construction of block linear systems, like those
//                  found in AC analysis.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 06/22/11
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_LAS_SystemHelpers.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Operator.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_EpetraProblem.h>

#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Linear {

  // Non-member creation methods.
MultiVector* createMultiVector( N_PDS_ParMap & map,
                                int numVectors )
{
  return new MultiVector( map, numVectors );
}

MultiVector* createMultiVector( N_PDS_ParMap & map,
                                N_PDS_ParMap & ol_map,
                                int numVectors )
{
  return new MultiVector( map, ol_map, numVectors );
}

Vector* createVector( N_PDS_ParMap & map ) 
{
  return new Vector( map );
}

Vector* createVector( N_PDS_ParMap & map, N_PDS_ParMap & ol_map )
{
  return new Vector( map, ol_map );
}

Matrix* createMatrix( const Graph* overlapGraph,
                      const Graph* baseGraph )
{
  return new Matrix( overlapGraph, baseGraph );
}

Graph* createGraph( N_PDS_ParMap & map, 
                    const std::vector<int>& numIndicesPerRow )
{
  return new Graph( map, numIndicesPerRow );
}

Graph* createGraph( N_PDS_ParMap & map,
                    int maxNumIndicesPerRow )
{
  return new Graph( map, maxNumIndicesPerRow );
}

Problem* createProblem( Matrix* A, MultiVector* x, MultiVector* b )
{
  return new EpetraProblem( A, x, b );
}
               
Problem* createProblem( Operator* Op, MultiVector* x, MultiVector* b )
{
  return new EpetraProblem( Op, x, b );
}
               
// ///////////////////////////////////////////////////////////////////
//
// Implementation of the Xyce::Linear::EpetraTransOp class.
//
// ///////////////////////////////////////////////////////////////////

EpetraTransOp::EpetraTransOp (const Teuchos::RCP<Epetra_Operator> &Op)
  : Epetra_Op(Op)
{
}

// The version of Apply() that takes two arguments and returns int
// implements the Epetra_Operator interface.
int
EpetraTransOp::Apply (const Epetra_MultiVector &X,
                     Epetra_MultiVector &Y) const
{
  // This operation computes Y = A^{-1}*X.
  const int info = Epetra_Op->Apply( X, Y );

  return info;
}

// This implements Epetra_Operator::ApplyInverse().
int
EpetraTransOp::ApplyInverse (const Epetra_MultiVector &X,
                            Epetra_MultiVector &Y) const
{
  // This operation computes Y = A*X.
  const int info = Epetra_Op->ApplyInverse( X, Y );

  return info;
}


} // namespace Linear
} // namespace Xyce

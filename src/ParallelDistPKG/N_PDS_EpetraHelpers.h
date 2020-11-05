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
//                  in the construction of parallel distribution objects
//                  in Epetra
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 09/29/20
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_PDS_EPETRAHELPERS_H
#define  Xyce_PDS_EPETRAHELPERS_H

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------

#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>

// ---------- Forward Declarations ----------

class Epetra_Comm;

namespace Xyce {
namespace Parallel {

// Return a new ParMap
N_PDS_ParMap * createPDSParMap( int & numGlobalEntities,
                                int numLocalEntities,
                                const std::vector<int> & lbMap,
                                const int index_base,
                                N_PDS_Comm & aComm );

// Return a new ParMap
//  -> let the underlying linear algebra determine the IDs.
N_PDS_ParMap * createPDSParMap( int & numGlobalEntities,
                                int numLocalEntities,
                                const int index_base,
                                N_PDS_Comm & aComm );

// Return a new Comm
N_PDS_Comm * createPDSComm(int iargs = 0, char * cargs[] = 0, Xyce::Parallel::Machine comm = MPI_COMM_NULL );
N_PDS_Comm * createPDSComm( Epetra_Comm* comm );

const Epetra_Comm* getEpetraComm( const N_PDS_Comm* comm );
Epetra_Comm* getEpetraComm( N_PDS_Comm* comm );

} // namespace Parallel 
} // namespace Xyce

#endif

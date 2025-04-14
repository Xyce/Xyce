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

#include <Xyce_config.h>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#endif

#include <N_PDS_ParHelpers.h>
#include <N_PDS_EpetraHelpers.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_EpetraParMap.h>

#include <N_PDS_Comm.h>
#include <N_PDS_EpetraMPIComm.h>
#include <N_PDS_EpetraSerialComm.h>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Function      : createPDSParMap
// Purpose       : Creates an ParMap object based on linear algebra.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
ParMap * createPDSParMap( int & numGlobalEntities,
                          int numLocalEntities,
                          const std::vector<int> & lbMap,
                          const int index_base,
                          const Communicator & aComm )
{
  const int * mArray = lbMap.empty() ? 0 : static_cast<const int *>(&lbMap[0]);

  // fix for empty maps
  int nGE = std::max( -1, numGlobalEntities );
  int nLE = std::max( 0, numLocalEntities );
  // Call the Petra constructor for the true Petra map.
  const Epetra_Comm* petraComm = Xyce::Parallel::getEpetraComm( &aComm );
  Epetra_Map*  petraMap = new Epetra_Map( nGE,
                              nLE,
                              mArray,
                              index_base,
                              *petraComm );
  Communicator& nonconst_comm = const_cast<Communicator&>(aComm);

  return ( new EpetraParMap( petraMap, nonconst_comm, true ) );
}

//-----------------------------------------------------------------------------
// Function      : createPDSParMap
// Purpose       : Creates an ParMap object based on linear algebra.
// Special Notes : let the underlying linear algebra determine the IDs
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
ParMap * createPDSParMap( int & numGlobalEntities,
                          int numLocalEntities,
                          const int index_base,
                          const Communicator & aComm )
{
  const Epetra_Comm* petraComm = Xyce::Parallel::getEpetraComm( &aComm );
  Epetra_Map*  petraMap = new Epetra_Map( numGlobalEntities,
                                          numLocalEntities,
                                          index_base,
                                          *petraComm );
  Communicator& nonconst_comm = const_cast<Communicator&>(aComm);

  return ( new EpetraParMap( petraMap, nonconst_comm, true ) );
}

//-----------------------------------------------------------------------------
// Function      : createPDSComm
// Purpose       : Creates an Epetra serial or parallel comm object based on
//                 system.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
Communicator * createPDSComm( int iargs, char * cargs[], Machine comm)
{
  Communicator * theComm = NULL;

#ifdef Xyce_PARALLEL_MPI
  if (comm != MPI_COMM_NULL)
    theComm = new EpetraMPIComm( comm );
  else
    theComm = new EpetraMPIComm( iargs, cargs );
#else
  theComm = new EpetraSerialComm();
#endif
  return theComm;
}

//-----------------------------------------------------------------------------
// Function      : createPDSComm
// Purpose       : Creates an Communicator based on either the Epetra serial or 
//               : parallel comm object.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
Communicator * createPDSComm( Epetra_Comm* comm )
{
  Communicator * pdsComm = NULL;

#ifdef Xyce_PARALLEL_MPI
  Epetra_MpiComm * mpicomm = dynamic_cast<Epetra_MpiComm *>( comm );

  if (mpicomm)
    pdsComm = new EpetraMPIComm( mpicomm->Comm() );
  else
    pdsComm = new EpetraMPIComm( MPI_COMM_WORLD );
#else
  pdsComm = new EpetraSerialComm();
#endif

  return pdsComm;
}

//-----------------------------------------------------------------------------
// Function      : createPDSComm
// Purpose       : Creates an Communicator based on either the Epetra serial or 
//               : parallel comm object.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
const Communicator * createPDSComm( const Epetra_Comm* comm )
{
  Communicator * pdsComm = NULL;

#ifdef Xyce_PARALLEL_MPI
  const Epetra_MpiComm * mpicomm = dynamic_cast<const Epetra_MpiComm *>( comm );

  if (mpicomm)
    pdsComm = new EpetraMPIComm( mpicomm->Comm() );
  else
    pdsComm = new EpetraMPIComm( MPI_COMM_WORLD );
#else
  pdsComm = new EpetraSerialComm();
#endif

  return pdsComm;
}

//-----------------------------------------------------------------------------
// Function      : getEpetraComm
// Purpose       : Gets an Epetra serial or parallel comm object based on
//                 system.
// Special Notes : Non-const communicator
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/29/20
//-----------------------------------------------------------------------------
Epetra_Comm* getEpetraComm( Communicator* comm )
{ 
  Epetra_Comm* petraComm = 0;

#ifdef Xyce_PARALLEL_MPI
  EpetraMPIComm* mpiComm = dynamic_cast<EpetraMPIComm*>( comm );
  if ( mpiComm )
    petraComm = mpiComm->petraComm();
#endif
  
  if (!petraComm)
  { 
    EpetraSerialComm* serialComm = dynamic_cast<EpetraSerialComm*>( comm );
    if ( serialComm )
      petraComm = serialComm->petraComm();
  }
  
  return petraComm;
}

//-----------------------------------------------------------------------------
// Function      : getEpetraComm
// Purpose       : Gets an Epetra serial or parallel comm object based on
//                 system.
// Special Notes : Const communicator
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/29/20
//-----------------------------------------------------------------------------
const Epetra_Comm* getEpetraComm( const Communicator* comm )
{ 
  const Epetra_Comm* petraComm = 0;

#ifdef Xyce_PARALLEL_MPI
  const EpetraMPIComm* mpiComm = dynamic_cast<const EpetraMPIComm*>( comm );
  if ( mpiComm )
    petraComm = mpiComm->petraComm();
#endif
  
  if (!petraComm)
  { 
    const EpetraSerialComm* serialComm = dynamic_cast<const EpetraSerialComm*>( comm );
    if ( serialComm )
      petraComm = serialComm->petraComm();
  }
  
  return petraComm;
}

} // namespace Parallel
} // namespace Xyce



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

#include <Xyce_config.h>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <N_PDS_EpetraHelpers.h>

#include <N_PDS_Comm.h>
#include <N_PDS_EpetraMPIComm.h>
#include <N_PDS_EpetraSerialComm.h>

#include <Epetra_Comm.h>

//-----------------------------------------------------------------------------
// Function      : createPDSComm
// Purpose       : Creates an Epetra serial or parallel comm object based on
//                 system.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_Comm *Xyce::Parallel::createPDSComm( int iargs, char * cargs[], Xyce::Parallel::Machine comm)
{
  N_PDS_Comm * theComm = NULL;

#ifdef Xyce_PARALLEL_MPI
  if (comm != MPI_COMM_NULL)
    theComm = new N_PDS_EpetraMPIComm( comm );
  else
    theComm = new N_PDS_EpetraMPIComm( iargs, cargs );
#else
  theComm = new N_PDS_EpetraSerialComm();
#endif
  return theComm;
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
Epetra_Comm* Xyce::Parallel::getEpetraComm( N_PDS_Comm* comm )
{ 
  Epetra_Comm* petraComm = 0;

#ifdef Xyce_PARALLEL_MPI
  N_PDS_EpetraMPIComm* mpiComm = dynamic_cast<N_PDS_EpetraMPIComm*>( comm );
  if ( mpiComm )
    petraComm = mpiComm->petraComm();
#endif
  
  if (!petraComm)
  { 
    N_PDS_EpetraSerialComm* serialComm = dynamic_cast<N_PDS_EpetraSerialComm*>( comm );
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
const Epetra_Comm* Xyce::Parallel::getEpetraComm( const N_PDS_Comm* comm )
{ 
  const Epetra_Comm* petraComm = 0;

#ifdef Xyce_PARALLEL_MPI
  const N_PDS_EpetraMPIComm* mpiComm = dynamic_cast<const N_PDS_EpetraMPIComm*>( comm );
  if ( mpiComm )
    petraComm = mpiComm->petraComm();
#endif
  
  if (!comm)
  { 
    const N_PDS_EpetraSerialComm* serialComm = dynamic_cast<const N_PDS_EpetraSerialComm*>( comm );
    if ( serialComm )
      petraComm = serialComm->petraComm();
  }
  
  return petraComm;
}


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


//-------------------------------------------------------------------------
//
// Purpose        : Specification file for testing the ParallelDist package.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/22/00
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Test_h
#define Xyce_N_PDS_Test_h

// ---------- Standard Includes ----------

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_PDS_Manager.h>

//-----------------------------------------------------------------------------
// Class         : N_PDS_ParTest
// Purpose       : This is the top level class for the parallel testing
//                 program.  The member function, RunTests, is the "main"
//                 function, essentially.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

class N_PDS_ParTest
{

public:

  // Default constructor & destructor.
  N_PDS_ParTest();
  ~N_PDS_ParTest();

  int RunTests(int iargs, char *cargs[], N_PDS_ParTest *ParTest);

protected:

private:
  N_PDS_Manager* parMgrPtr;

};

//-----------------------------------------------------------------------------

#endif

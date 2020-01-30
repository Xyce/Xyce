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
// Purpose        : This function is the test driver for the ParallelDist
//                  package.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/22/00
//
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

#include <iostream>
#include <vector>
#include <list>
#include <string>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_PDS_Test.h>
#include <N_PDS_Manager.h>
#include <N_ERH_ErrorMgr.h>

// Some of these are temporary includes - for testing purposes only!

// Class N_PDS_ParTest

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParTest::N_PDS_ParTest
// Purpose       : Default constructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

N_PDS_ParTest::N_PDS_ParTest()
{

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParTest::~N_PDS_ParTest
// Purpose       : Default destructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

N_PDS_ParTest::~N_PDS_ParTest()
{

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParTest::RunTests
// Purpose       : Performs tests on ParallelDist package.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

int N_PDS_ParTest::RunTests(int iargs, char *cargs[], N_PDS_ParTest *ParTest)

{
  static int iSuccess = false;

#ifdef Xyce_PARALLEL_MPI
  static string              lbMethod = ("PARMETIS");
  static list<string_params> lbParams;
#endif

  cout  << endl <<
    "Welcome to the Xyce(TM) ParallelDist testing program." << endl <<
    endl;

#ifdef Xyce_PARALLEL_MPI

  if (iargs > 3)
  {
    if (iargs%2 != 0)
      Xyce::Report::DevelFatal0().in("N_PDS_ParTest::RunTests") << "wrong number of arguments.";

    // Setup the arguments for the manager constructor.
    for (int i = 3; i < iargs; i += 2)
    {
      string_params sp(cargs[i], cargs[i+1]);
      lbParams.push_back(sp);
    }
  }

#endif

  N_PDS_Manager *parMgr = new N_PDS_Manager(atoi(cargs[1])
#ifdef Xyce_PARALLEL_MPI
                                            , lbMethod, lbParams
#endif
                                            );

  iSuccess = (parMgr != NULL);

  parMgr->reportLoadBalance();

  if (iSuccess == false)
    cout << "Test of ParallelDist NOT completed successfully." <<
      endl;
  else
    cout << "Test of ParallelDist completed successfully." << endl;

  delete parMgr;

  return iSuccess;

}

//-----------------------------------------------------------------------------
// Function      : main
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/12/00
//-----------------------------------------------------------------------------

int main(int iargs, char *cargs[])

{
  N_PDS_ParTest *ParTest = new N_PDS_ParTest();

  if (iargs < 2)
  {
    Xyce::Report::DevelFatal0().in("main") << "Wrong number of arguments.";
  }

  ParTest->RunTests(iargs, cargs, ParTest);

  delete ParTest;
  exit(0);
}

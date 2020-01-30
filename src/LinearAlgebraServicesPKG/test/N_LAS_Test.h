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
// Purpose        : Specification file for testing the LinearAlgebraServices
//                  package.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/23/00
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_LAS_Test_h
#define Xyce_N_LAS_Test_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_Manager.h>
#include <N_LAS_LAFactory.h>

//-----------------------------------------------------------------------------
// Class         : N_LAS_LATest
// Purpose       : This is the top level class for the linear-algebra testing
//                 program.  The member function, RunTests, is the "main"
//                 function, essentially.
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Compuational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------

class N_LAS_LATest
{

public:

  // Default constructor & destructor.
  N_LAS_LATest();
  ~N_LAS_LATest();

  int RunTests(int iargs, char *cargs[], N_LAS_LATest *LATest);

  int vectorTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                  int numVectors);
  int matrixVectorTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                        int numVectors);

  int solverTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                  int numVectors);

protected:

private:

};

//-----------------------------------------------------------------------------

#endif

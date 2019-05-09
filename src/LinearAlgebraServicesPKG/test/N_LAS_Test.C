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
// Purpose        : Implementation file for testing the LinearAlgebraServices
//                  package.
//
// Special Notes  :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/23/00
//
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------

#include <iostream>
#include <string>
#include <limits>

// ----------   Xyce Includes   ----------

#include <N_LAS_Test.h>
#include <N_PDS_Manager.h>
#include <N_ERH_ErrorMgr.h>

// Some of these are temporary includes - for testing purposes only!

// Class N_LAS_LATest

//-----------------------------------------------------------------------------
// Function      : N_LAS_LATest::N_LAS_LATest
// Purpose       : Default constructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 05/23/00
//-----------------------------------------------------------------------------

N_LAS_LATest::N_LAS_LATest()
{

}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LATest::~N_LAS_LATest
// Purpose       : Default destructor.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

N_LAS_LATest::~N_LAS_LATest()

{
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LATest::vectorTests
// Purpose       : Performs tests on N_LAS_MultiVector objects.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/01/00
//-----------------------------------------------------------------------------

int N_LAS_LATest::vectorTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                              int numVectors)

{

  static const string msg("N_LAS_LATest::vectorTests - ");

  bool                iSuccess = true;
  double              result = 0.0, alpha = 1.0, beta = 2.0;
  N_LAS_MultiVector  *a = NULL, *b = NULL;

  // Initially, only tests for single vectors.
  if (numVectors > 1)
  {
    static const string errorMsg("numVectors > 1");
    Xyce::Report::DevelFatal0() << msg + errorMsg;
  }

  a = factory->newVector(result, *parMap, 1);
  b = factory->newVector(result, *parMap, 1);

  a->putScalar(alpha);
  b->putScalar(beta);

  result = a->dotProduct(*b);

  if (result != (beta * parMap->getNumGlobalEntities()))
  {
    static const string errorMsg("One or more of the following "
                                      "vector function tests failed:"
                                      "\n\t* dotProduct\n\t* putScalar");
    cout << "Error - result: " << result << endl;
    Xyce::Report::DevelFatal0() <<  msg + errorMsg;
  }
  else
  {
    Xyce::lout() << Xyce::section_divider << std::endl
           << msg + "All of the following vector function tests passed:\n\t* dotProduct\n\t* putScalar" << std::endl
           << Xyce::section_divider << std::endl;
  }

  // Randomize
  a->random();

  // MultiVector copy "constructor"
  N_LAS_MultiVector  *alphaA = factory->newVector(*a);

  // Testing scale, addVec and lpNorm
  alphaA->scale(beta);
  alphaA->addVec(-beta, *a);
  alphaA->lpNorm(2, &result);

//  if (result > numeric_limits<double>::epsilon())  // I don't know why this
  // doesn't work - it should be part of STL.
  if (result > DBL_EPSILON)
  {
    static const string
      errorMsg("One or more of the following vector function tests failed:"
               "\n\t* random\n\t* newVector copy constructor\n\t* scale"
               "\n\t* addVec\n\t* lpNorm");
    cout << "Error - result: " << result << endl;
    Xyce::Report::DevelFatal0() << msg + errorMsg;
  }
  else
  {
    Xyce::lout() << Xyce::section_divider << std::endl
                 << msg << "All of the following vector function tests passed:" << std::endl
                 << "\t* newVector copy constructor\n\t* scale\n\t* addVec\n\t* lpNorm" << std::endl
                 << Xyce::section_divider << std::endl;
  }

  // Testing daxpy.
  *b = *a;
  b->scale(-beta);
  alphaA->daxpy(*b, beta, *a);
  alphaA->lpNorm(2, &result);

//  if (result > numeric_limits<double>::epsilon())  // I don't know why this
  // doesn't work - it should be part of STL.
  if (result > DBL_EPSILON)
  {
    static const string
      errorMsg("One or more of the following vector function tests failed:"
               "\n\t* = operator\n\t* scale\n\t* daxpy\n\t* lpNorm");

    cout << "Error - result: " << result << endl;
    Xyce::Report::DevelFatal0() <<  msg + errorMsg;
  }
  else
  {
    Xyce::lout() << Xyce::section_divider << std::endl
                 << msg
                 << "All of the following vector function tests passed:\n" << std::endl
                 << "\t* = operator\n\t* scale\n\t* daxpy\n\t* lpNorm" << std::endl
                 << Xyce::section_divider << std::endl;
  }

  // Testing indexing...
//   for (int i = 0; i < parMap->getNumGlobalEntities(); i++)
//   {
//     result = a[i];
//     cout << "a[" << i << "] = " << result << endl;
//   }

  // Clean up
  delete a;
  delete b;
  delete alphaA;

  return iSuccess;

}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LATest::matrixVectorTests
// Purpose       : Performs tests on N_LAS_Matrix and N_LAS_MultiVector
//                 objects using a tridiagonal matrix and the Power Method.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/01/00
//-----------------------------------------------------------------------------

int N_LAS_LATest::matrixVectorTests(N_PDS_ParMap *parMap,
                                    N_LAS_LAFactory *factory, int numVectors)

{

  static const string msg("N_LAS_LATest::matrixVectorTests - ");
  bool                iSuccess = true;
  double              result   = 0.0;
  N_LAS_Matrix       *A        = NULL;
  N_LAS_MultiVector  *q = NULL, *z = NULL, *resid = NULL;

  // Initially, only tests for single vectors.
  if (numVectors > 1)
  {
    static const string errorMsg("numVectors > 1");
    Xyce::Report::DevelFatal0() <<  msg + errorMsg;
  }

  A = factory->newMatrix(result, *parMap);

  // Create an integer vector numNz that is used to build the matrix.  numNz[i]
  // is the number of OFF-DIAGONAL term for the ith global equation on this
  // processor.  This code is largely borrowed from Mike Heroux and Petra.
  int  numLocalEquations  = parMap->getNumLocalEntities();
  int  numGlobalEquations = parMap->getNumGlobalEntities();
  int *numNz              = new int[numLocalEquations];
  int *updateList         = parMap->getParMap();

  for (int i = 0; i < numLocalEquations ; i++)
    if (updateList[i]==0 || updateList[i] == numGlobalEquations-1)
      numNz[i] = 1;
    else
      numNz[i] = 2;

  // Allocate space using numNz
  A->allocate(numNz);

  // Add  rows one-at-a-time.  Need some vectors to help off diagonal values
  // will always be -1.

  double *values  = new double[2];
  int    *indices = new int[2];
  double  two     = 2.0;
  int     numEntries;

  values[0] = -1.0; values[1] = -1.0;

  for (int i = 0; i < numLocalEquations; i++)
  {
    if (updateList[i] == 0)
    {
      indices[0] = 1;
      numEntries = 1;
    }
    else if (updateList[i] == numGlobalEquations - 1)
    {
      indices[0] = numGlobalEquations - 2;
      numEntries = 1;
    }
    else
    {
      indices[0] = updateList[i] - 1;
      indices[1] = updateList[i] + 1;
      numEntries = 2;
    }

    A->putRow(updateList[i], numEntries, values, indices);

    // Put in the diagonal entry.
    A->putRow(updateList[i], 1, &two, updateList+i);
  }

  // Finish up.
  A->fillComplete();

  // Create vectors and variables for Power method.
  q     = factory->newVector(result, *parMap, 1);
  z     = factory->newVector(result, *parMap, 1);
  resid = factory->newVector(result, *parMap, 1);

  double normz, lambda, residual;

  // Fill z with random numbers.
  z->random();

  // Iterate
  int    nIters    = 100 * numGlobalEquations;
  double tolerance = pow(10, -10.0 + 2.0*log10((double) numGlobalEquations));

  for (int iter = 0; iter < nIters; iter++)
  {
    z->lpNorm(2, &normz);       // Compute 2-norm of z
    *q    = *z;
    q->scale(1.0/normz);
    A->matvec(*q, *z);          // Compute z = A*q
    lambda = q->dotProduct(*z); // Approximate maximum eigenvalue

    if (iter % 500 == 0 || iter + 1 == nIters)
    {
      resid->daxpy(*z, -lambda, *q); // Compute A*q - lambda*q
      resid->lpNorm(2, &residual));

      if (parMap->getPDSComm()->getProcID() == 0)
        cout << "Iter = " << iter << "  normz = " << normz << "  Lambda = " <<
          lambda << "  Residual of A*q - lambda*q = " << residual << endl;
    }

    if (residual < tolerance) break;
  }

  if (residual > tolerance)
    iSuccess = false;
  else
  {
    Xyce::lout() << Xyce::section_divider << std::endl
                 << msg  <<"Power method converged\n\n" << std::endl
                 << Xyce::section_divider << std::endl;
  }

  // Clean up
  delete [] numNz;
  delete [] values;
  delete [] indices;

  delete A;
  delete q;
  delete z;
  delete resid;

  return iSuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LATest::solverTests
// Purpose       : Performs tests on N_LAS solvers.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 06/12/00
//-----------------------------------------------------------------------------

int N_LAS_LATest::solverTests(N_PDS_ParMap *parMap, N_LAS_LAFactory *factory,
                              int numVectors)

{

  static const string msg("N_LAS_LATest::solverTests - ");
  bool                iSuccess = true;
  double              result   = 0.0;
  N_LAS_Matrix       *A        = NULL;
  N_LAS_MultiVector  *x = NULL, *b = NULL;
  N_LAS_IterativeSolver *solver = NULL;

  // Initially, only tests for single vectors.
  if (numVectors > 1)
  {
    static const string errorMsg("numVectors > 1");
    Xyce::Report::DevelFatal0() << msg + errorMsg;
  }

  A = factory->newMatrix(result, *parMap);
  x = factory->newVector(result, *parMap, 1);
  b = factory->newVector(result, *parMap, 1);

  // Setup a test problem in A and b.

  // Create an integer vector numNz that is used to build the matrix.  numNz[i]
  // is the number of OFF-DIAGONAL term for the ith global equation on this
  // processor.  This code is largely borrowed from Mike Heroux and Petra.
  int  numLocalEquations  = parMap->getNumLocalEntities();
  int  numGlobalEquations = parMap->getNumGlobalEntities();
  int *numNz              = new int[numLocalEquations];
  int *updateList         = parMap->getParMap();

  for (int i = 0; i < numLocalEquations ; i++)
    if (updateList[i]==0 || updateList[i] == numGlobalEquations-1)
      numNz[i] = 1;
    else
      numNz[i] = 2;

  // Allocate space using numNz
  A->allocate(numNz);

  // Add  rows one-at-a-time.  Need some vectors to help off diagonal values
  // will always be -1.

  double *values  = new double[2];
  double *zeroes  = new double[2];
  int    *indices = new int[2];
  double  two     = 2.0, bc = 10.0, one = 1.0;
  int     numEntries;

  values[0] = -1.0; values[1] = -1.0;
  zeroes[0] = 0.0;  zeroes[1] = 0.0;

  for (int i = 0; i < numLocalEquations; i++)
  {
    if (updateList[i] == 0)
    {
      indices[0] = 1;
      numEntries = 1;
    }
    else if (updateList[i] == numGlobalEquations - 1)
    {
      indices[0] = numGlobalEquations - 2;
      numEntries = 1;
    }
    else
    {
      indices[0] = updateList[i] - 1;
      indices[1] = updateList[i] + 1;
      numEntries = 2;
    }

    if (updateList[i] == 0)
      A->putRow(updateList[i], numEntries, zeroes, indices);
    else
      A->putRow(updateList[i], numEntries, values, indices);

    // Put in the diagonal entry.
    if (updateList[i] == 0)
      A->putRow(updateList[i], 1, &one, updateList+i);
    else
      A->putRow(updateList[i], 1, &two, updateList+i);
  }

  // Finish up.
  A->fillComplete();

  // Fill a RHS vector b and an initial guess for x.
  b->putScalar(bc);
  x->random();

  // Create solver.
  solver = factory->newIterativeSolver(*A, *x, *b);

  // Setup the parameters
  solver->setParams(200, 1.0e-08, 3);

  // Solve the problem.
  iSuccess = solver->solve();

  if (iSuccess)
  {
    Xyce::lout() << Xyce::section_divider << std::endl
           << msg + "Solver converged\n\n" << std::endl
           << Xyce::section_divider << std::endl;
  }

  // Clean up
  delete [] numNz;
  delete [] values;
  delete [] zeroes;
  delete [] indices;

  delete A;
  delete x;
  delete b;

  return iSuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_LAS_LATest::RunTests
// Purpose       : Performs tests on LinearAlgebraServices package.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/22/00
//-----------------------------------------------------------------------------

int N_LAS_LATest::RunTests(int iargs, char *cargs[], N_LAS_LATest *LATest)

{

  static int iSuccess = 1, vectorSuccess = 1, matrixVectorSuccess = 1,
    solverSuccess = 1;
  static int vectorLength = atoi(cargs[1]);

  // Error string
  static const string msg("N_LAS_LATest::RunTests - ");

#ifdef Xyce_PARALLEL_MPI
  static string             lbMethod = ("PARMETIS");
  static list<string_param> lbParams;
#endif

  N_PDS_Manager *parMgr = new N_PDS_Manager(vectorLength
#ifdef Xyce_PARALLEL_MPI
                                            , lbMethod, lbParams
#endif
                                            );

  iSuccess = (parMgr != NULL);

  Xyce::lout() << "\n\n\tWelcome to the "
                         "Xyce(TM << std::endl LinearAlgebraServices testing "
                         "program.\n\n" << std::endl;

#ifdef Xyce_PARALLEL_MPI

  if (iargs > 3)
  {
    if (iargs%2 != 0)
      Xyce::Report::DevelFatal0() <<  msg + "wrong number of arguments.";

    // Setup the arguments for the manager constructor.
    for (int i = 3; i < iargs; i += 2)
    {
      string_param sp(cargs[i], cargs[i+1]);
      lbParams.push_back(sp);
    }
  }

#endif

  parMgr->reportLoadBalance();

  int                 numVectors = 1;

  N_PDS_LoadBalance  *LAS_LB  = parMgr->getLAS_LB();
  N_PDS_LoadBalance  *DEV_LB  = parMgr->getDEV_LB();
  N_PDS_ParMap       *parMap  = LAS_LB->getParallelMap();
  N_LAS_LAFactory    *factory = new N_LAS_LAFactory();

  vectorSuccess = vectorTests(parMap, factory, numVectors);

  matrixVectorSuccess = matrixVectorTests(parMap, factory, numVectors);

  solverSuccess = solverTests(parMap, factory, numVectors);

  iSuccess = (iSuccess && vectorSuccess && matrixVectorSuccess &&
              solverSuccess);

  if (!iSuccess)
  {
    static string
      errorMsg("Test of LinearAlgebraServices NOT completed successfully.");

    if (!vectorSuccess)
      errorMsg += "\nFailed in vectorTests.";
    else if (!matrixVectorSuccess)
      errorMsg += "\nFailed in matrixVectorTests.";
    else if (!solverSuccess)
      errorMsg += "\nFailed in solverTests.";

    Xyce::Report::DevelFatal0() << msg + errorMsg;
  }

  Xyce::lout() << "Test of "
                         "LinearAlgebraServices completed successfully." << std::endl;

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

  N_LAS_LATest *LATest = new N_LAS_LATest();

  if (iargs < 2)
  {
    // Error string
    static const string error_msg("main - Wrong number of arguments.");
    cout << error_msg << endl;
    exit(-1);
  }

  LATest->RunTests(iargs, cargs, LATest);

  delete LATest;
  exit(0);
}

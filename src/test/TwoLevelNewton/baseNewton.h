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

//-----------------------------------------------------------------------------
// Filename       : baseNewton.h
//
// Purpose        : base class
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 12/5/2018
//
//-----------------------------------------------------------------------------

#ifndef Xyce_baseNewton_h
#define Xyce_baseNewton_h

#include <Xyce_config.h>

#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>
#include <N_TIA_TwoLevelError.h>
#include <N_ANP_fwd.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include <vector>

// Function to be called if memory runs out:
void _new_handler (void)
{
  Xyce::Report::UserFatal0() << "OUT OF MEMORY (error in 'new')";
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class baseNewton
{
  public:
    // variables: 
     std::string netlistFilename_;

  public:

     // functions:
     // constructor
  private:
     baseNewton();

  public:
     baseNewton(std::string & netlist): netlistFilename_ (netlist) { };

     // destructor
    virtual ~baseNewton() {};

    void allocateLinearSystem(int size);
    bool solveLinearSystem();

    virtual void updateLinearSystem() = 0;
    virtual bool checkConvergence() = 0;

    bool newtonSolve();

    // top level linear system:
    Teuchos::SerialDenseMatrix< int, double > JacobianMatrix;
    Teuchos::SerialDenseMatrix< int, double > RHS;
    Teuchos::SerialDenseMatrix< int, double > dX;
    Teuchos::SerialDenseMatrix< int, double > X;
};

//-----------------------------------------------------------------------------
inline void baseNewton::allocateLinearSystem(int size)
{
  // "shape" is equivalent to resize for these objects
  JacobianMatrix.shape(size,size);

  RHS.shape(size,1);
  dX.shape(size,1);
  X.shape(size,1);

  // initialize everything to zero
  JacobianMatrix.putScalar(0.0);
  RHS.putScalar(0.0);
  dX.putScalar(0.0);
  X.putScalar(0.0);
}

//-----------------------------------------------------------------------------
inline bool baseNewton::solveLinearSystem()
{
  Teuchos::LAPACK< int, double > lapack;
  int size = JacobianMatrix.numRows();
  std::vector<int> ipiv(size,0);
  int info1, info2;

  //! Computes an LU factorization of a general \c m by \c n matrix \c A using partial pivoting with row interchanges.
  lapack.GETRF( JacobianMatrix.numRows(), JacobianMatrix.numCols(), JacobianMatrix.values(), JacobianMatrix.numRows(), &ipiv[0], &info1) ;

  //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general \c n by \c n matrix \c A using the LU factorization computed by GETRF.
  lapack.GETRS( 'N', JacobianMatrix.numRows(), 1, JacobianMatrix.values(), JacobianMatrix.numRows(), &ipiv[0], RHS.values(), JacobianMatrix.numRows(), &info2) ;

#if 0
  Xyce::lout() << "LAPACK linear solve result" << std::endl;
  Xyce::lout() << "info1 = " << info1 <<std::endl;
  Xyce::lout() << "info2 = " << info1 <<std::endl;
  for (int i=0;i<size;++i)
  {
    Xyce::lout() << "RHS("<<i<<") = " << RHS(i,0) <<std::endl;
  }
#endif

  // check this later
  return true;
}

//-----------------------------------------------------------------------------
//
// This is the function executing the top level Newton solve.  
//
//-----------------------------------------------------------------------------
inline bool baseNewton::newtonSolve()
{
  int iNewt=0;
  bool converged = false;
  bool bsuccess  = false;

  int maxSteps = 20;

  Xyce::lout() << "--------------\nX vector:";
  X.print(Xyce::lout());

  while (!converged)
  {
    Xyce::lout() << "Outer loop Newton iteration : " << iNewt <<std::endl;

    // sum the contributions into the linear system, which has the form of J * dx = -f
    // solve linear system for dX and update X. 
    updateLinearSystem();
    solveLinearSystem();

    dX = RHS;
    X += dX;

#if 0
    Xyce::lout() << "--------------\ndX vector:";
    dX.print(Xyce::lout());
    Xyce::lout() << "--------------\nX vector:";
    X.print(Xyce::lout());
#endif
  
    converged = checkConvergence();

    if (iNewt >= maxSteps) 
    {
      Xyce::lout() << "maximum number of Outer loop Newton steps reached.  iNewt = " << iNewt <<std::endl;
      break;
    }

    ++iNewt;
  }

  return converged;
}

#endif



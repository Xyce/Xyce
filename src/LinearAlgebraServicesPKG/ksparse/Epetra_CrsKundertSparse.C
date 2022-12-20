//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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


#include "Epetra_CrsKundertSparse.h"

extern "C" {
#include "spmatrix.h"
  int spFactorAndSolve(char *eMatrix, double *RHS); // Sparse has no prototype for this function

}

#include "Epetra_LinearProblem.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"


Epetra_CrsKundertSparse::Epetra_CrsKundertSparse( Epetra_LinearProblem * Problem,
			      const double RelThreshold,
			      const double AbsThreshold,
			      const int DiagPivoting)

  : RelThreshold_(RelThreshold),
    AbsThreshold_(AbsThreshold),
    DiagPivoting_(DiagPivoting),
    Problem_(Problem),
    FirstSolve_(1)
{

  Epetra_CrsMatrix * A = dynamic_cast<Epetra_CrsMatrix *> (Problem->GetOperator());
  MyPID_ = A->Comm().MyPID();
  NumMyRows_ = A->NumMyRows();
  NumMyCols_ = A->NumMyCols();
  NumGlobalRows_ = A->NumGlobalRows();
  NumGlobalCols_ = A->NumGlobalCols();

  if (A->Comm().NumProc()>1 && ((MyPID_==0)&&(NumMyRows_!=NumGlobalRows_))) 
    throw A->ReportError("Can only use Spice Sparse in serial on processor 0", -1);

  if (NumGlobalRows_ != NumGlobalCols_) 
    throw A->ReportError("Matrix must be square", -2);

  // Create a Sparse matrix
  int err = 0;
  if (MyPID_ == 0) {
    Matrix_ = (char *) spCreate (NumGlobalRows_, 0, &err);
    if (err!=0) 
      throw A->ReportError("Error occurred in Spice Sparse spCreate", err);
  } 
  else 
    Matrix_ = 0;

  int NumEntries;
  int * Indices;
  double * Values;
    
  int NumMyNonzeros = A->NumMyNonzeros();
  
  int curValue = 0;
  addr_list_ = new double *[NumMyNonzeros];
  for (int i=0; i<NumMyRows_; i++) {
    // int curGRID = A->RowMap().GID(i); // Needed for parallel (later)
    // View of current row
    int ierr = A->ExtractMyRowView(i, NumEntries, Values, Indices); 
    if (ierr!=0) 
      throw A->ReportError("Error occurred in ExtractMyRowView", ierr);
    for (int j=0; j<NumEntries; j++) {
      // int columnIndex = A->ImportMap().GID(Indices[j]);// parallel (later)
      int columnIndex = Indices[j];
      // Register this entry into Sparse matrix
      double * p = (double *) spGetElement(Matrix_, i+1, columnIndex+1);
      *p = Values[j];
      addr_list_[curValue++] = p;
    }
  }

}

Epetra_CrsKundertSparse::~Epetra_CrsKundertSparse() {
  deleteArrays();
}

void Epetra_CrsKundertSparse::deleteArrays() {

  if (Matrix_!=0) {
    spDestroy(Matrix_);
    Matrix_ = 0;
  }
  if (addr_list_!=0) {
    delete [] addr_list_;
    addr_list_ = 0;
  }

}

int Epetra_CrsKundertSparse::Solve(const bool ComputeFactor, const bool Transpose) {

  // Do some sanity checks and make some local pointers
  int orderStatus=0, factorStatus=0, solveStatus=0;
  EPETRA_CHK_ERR(Problem_->CheckInput());  // Check to make sure all problem parameters are well-defined.

  Epetra_CrsMatrix * A = dynamic_cast<Epetra_CrsMatrix *> (Problem_->GetOperator());
  if (A==0) EPETRA_CHK_ERR(-6); // Couldn't cast Operator to a CrsMatrix
  Epetra_MultiVector * X = Problem_->GetLHS();
  Epetra_MultiVector * B = Problem_->GetRHS();
  // If not first call to solver, we need to copy values to solver matrix.
  // NOTE: We are proceeding through the matrix in the same order as it was
  //       constructed.  As a result, we do not need to access index information.

  // We are only solving the linear system on processor 0.
  if (MyPID_ == 0) { 
    if ((FirstSolve_==0) && ComputeFactor) {
      spClear (Matrix_); // Clear previous factorization and matrix values
      int curValue = 0;
      int NumEntries;
      double * Values;
      for (int i=0; i<NumMyRows_; i++) {
      // int curGRID = A->RowMap().GID(i); // Needed for parallel (later)
      // View of current row
        EPETRA_CHK_ERR(A->ExtractMyRowView(i, NumEntries, Values)); 
        for (int j=0; j<NumEntries; j++)
	  *(addr_list_[curValue++]) = Values[j];
      }
    }

    /* Create right-hand side matrix B. */
    double ** rhsptrs;
    double ** solutionptrs;
    double * rhs;
    double * solution;
    B->ExtractView(&rhsptrs);
    X->ExtractView(&solutionptrs);
    rhs = rhsptrs[0];
    solution = solutionptrs[0];
    rhs--; solution--; // adjust for 1-based indexing

    if (ComputeFactor)
      orderStatus = spOrderAndFactor (Matrix_, NULL, RelThreshold_, AbsThreshold_, DiagPivoting_, FirstSolve_);

    if (!orderStatus)
    {
      if (Transpose)
        solveStatus = spSolveTransposed (Matrix_, rhs, solution, NULL, NULL);
      else
        solveStatus = spSolve (Matrix_, rhs, solution, NULL, NULL);
    }
    else
    {
      solveStatus = orderStatus;
    }

    if (FirstSolve_)
      FirstSolve_ = 0;

    // Check if there are more RHS to solve
    if (B->NumVectors()>1) {
      for (int i=1; i<B->NumVectors();i++) {
        rhs = rhsptrs[i];
        solution = solutionptrs[i];
        rhs--; solution--; // adjust for 1-based indexing
        if (Transpose)
          solveStatus = spSolveTransposed (Matrix_, rhs, solution, NULL, NULL);
        else
          solveStatus = spSolve (Matrix_, rhs, solution, NULL, NULL);
      }
    }
  }

  // Communicate failures to all processors.
  int tmpOrderStatus = orderStatus;
  int tmpFactorStatus = factorStatus;
  int tmpSolveStatus = solveStatus;
  A->Comm().MaxAll(&tmpOrderStatus, &orderStatus, 1);
  A->Comm().MaxAll(&tmpFactorStatus, &factorStatus, 1);
  A->Comm().MaxAll(&tmpSolveStatus, &solveStatus, 1);
  solveStatus += (orderStatus + factorStatus);

  return solveStatus;
}


/*
The original code in this file is Copyright (c) 1985-1991 The Regents of the
University of California and is under the Spice 3f5 BSD Copyright.

All additions and changes are under the following:
//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
*/

void timer_clear(), timer_start(), timer_stop();
double timer_get();

/* #define FINGERPRINT */
/* #define CHECK_ANS */
/*
#define WRITE_MAT
*/

#ifdef FINGERPRINT
void FingerPrint (double *, int *, int *, int);
#endif

/*
#define DEBUG_OVER
*/

#include <sys/types.h>

/*
 *  MATRIX FACTORIZATION MODULE
 *
 *  Author:                     Advising Professor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains the routines to factor the matrix into LU form.
 *
 *  >>> User accessible functions contained in this file:
 *  spOrderAndFactor
 *  spFactor
 *  spPartition
 *
 *  >>> Other functions contained in this file:
 *  FactorComplexMatrix         spcCreateInternalVectors
 *  CountMarkowitz              MarkowitzProducts
 *  SearchForPivot              SearchForSingleton
 *  QuicklySearchDiagonal       SearchDiagonal
 *  SearchEntireMatrix          FindLargestInCol
 *  FindBiggestInColExclude     ExchangeRowsAndCols
 *  spcRowExchange              spcColExchange
 *  ExchangeColElements         ExchangeRowElements
 *  RealRowColElimination       ComplexRowColElimination
 *  UpdateMarkowitzNumbers      CreateFillin
 *  MatrixIsSingular            ZeroPivot
 *  WriteStatus
 */


/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88,89,90
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any purpose and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the authors and the University of
 *  California are properly credited.  The authors and the University of
 *  California make no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without express
 *  or implied warranty.
 */

#ifdef notdef
static char copyright[] =
    "Sparse1.3: Copyright (c) 1985,86,87,88,89,90 by Kenneth S. Kundert";
static char RCSid[] =
    "@@(#)$Header$";
#endif


/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spConfig.h
 *    Macros that customize the sparse matrix routines.
 *  spMatrix.h
 *    Macros and declarations to be imported by the user.
 *  spDefs.h
 *    Matrix type and macro definitions for the sparse matrix routines.
 */

#define spINSIDE_SPARSE
#include "spconfig.h"
#include "spmatrix.h"
#include "spdefs.h"

/*
 * Function declarations
 */

static int  FactorComplexMatrix( MatrixPtr );
static void CountMarkowitz( MatrixPtr, RealVector, int );
static void MarkowitzProducts( MatrixPtr, int );
static ElementPtr SearchForPivot( MatrixPtr, int, int );
static ElementPtr SearchForSingleton( MatrixPtr, int );
static ElementPtr QuicklySearchDiagonal( MatrixPtr, int );
ElementPtr SearchDiagonal( MatrixPtr, int );
ElementPtr SearchEntireMatrix( MatrixPtr, int );
static RealNumber FindLargestInCol( ElementPtr );
static RealNumber FindBiggestInColExclude( MatrixPtr, ElementPtr, int );
static void ExchangeRowsAndCols( MatrixPtr, ElementPtr, int );
static void ExchangeColElements( MatrixPtr, int, ElementPtr, int,
                                 ElementPtr, int );
static void ExchangeRowElements( MatrixPtr, int, ElementPtr, int,
                                 ElementPtr, int );
void RealRowColElimination( MatrixPtr, ElementPtr, int );
static void ComplexRowColElimination( MatrixPtr, ElementPtr, int );
static void UpdateMarkowitzNumbers( MatrixPtr, ElementPtr );
static ElementPtr CreateFillin( MatrixPtr, int, int, int, ElementPtr *);
static int  MatrixIsSingular( MatrixPtr, int );
static int  ZeroPivot( MatrixPtr, int );
static void WriteStatus( MatrixPtr, int );
void add_fast_col_index (MatrixPtr, int, int, ElementPtr);
void add_fast_row_index (MatrixPtr, int, int, ElementPtr);
static void remove_fast_col_index (MatrixPtr, int, int, ElementPtr, ElementPtr);
static void remove_fast_row_index (MatrixPtr, int, int, ElementPtr, ElementPtr);
extern int f_ind(MatrixPtr, int, int);
void spCheckInd (MatrixPtr, char *);
void spSetIndex (MatrixPtr);
void spColInd (MatrixPtr, int);
void spRowInd (MatrixPtr, int);
void spExpandFormat (MatrixPtr Matrix);

#ifndef DEBUG
#define DEBUG
#endif
#undef DEBUG

/* #define CHECK_IND */
#define WIDTH 100
#define DEL_RED 1

/* routine to break optimization */
void dummy_call_x(int *i, int *j)
{
    if (*i+10 == -*j)
      printf ("*i+10=-*j\n");
    return;
}

void check_col(MatrixPtr Matrix, int Column, char* errmsg)
{
    ElementPtr p;
    int Row;

    p = Matrix->FirstInCol[Column];
    Row = 0;
    while (p != NULL) {
      if (p->Col != Column || p->Row <= Row) {
        printf ("Error found in column %d links: %s\n",Column,errmsg);
        break;
      }
      Row = p->Row;
      p = p->NextInCol;
    }
    return;
}

void print_col(MatrixPtr Matrix, int Column)
{
    ElementPtr p;
    int Row;

    p = Matrix->FirstInCol[Column];
    Row = 0;
    while (p != NULL) {
      printf ("Column entry: %d, Row = %d\n",p->Col,p->Row);
      if (p->Col != Column || p->Row <= Row) {
        printf ("Error found in column %d\n",Column);
        break;
      }
      Row = p->Row;
      p = p->NextInCol;
    }
    return;
}

void print_row(MatrixPtr Matrix, int Row)
{
    ElementPtr p;
    int Col;

    p = Matrix->FirstInRow[Row];
    Col = 0;
    while (p != NULL) {
      printf ("Row entry: %d, Col = %d\n",p->Row,p->Col);
      if (p->Row != Row || p->Col <= Col) {
        printf ("Error found in row %d\n",Row);
        break;
      }
      Col = p->Col;
      p = p->NextInRow;
    }
    return;
}

#ifdef CHECK_ANS
void make_copy (MatrixPtr Matrix)
{
    ElementPtr  pElement;
    int i, Size;

    Size = Matrix->Size;
    for (i=1 ; i<=Size ; i++) {
      pElement = Matrix->FirstInRow[i];
      while(pElement) {
        pElement->RealCopy = pElement->Real;
        pElement = pElement->NextInRow;
      }
    }

    return;
}

void swap_copy (MatrixPtr Matrix)
{
    ElementPtr  pElement;
    int i, Size;
    double temp;

    Size = Matrix->Size;
    for (i=1 ; i<=Size ; i++) {
      pElement = Matrix->FirstInRow[i];
      while(pElement) {
        temp = pElement->RealCopy;
        pElement->RealCopy = pElement->Real;
        pElement->Real = temp;
        pElement = pElement->NextInRow;
      }
    }

    return;
}
#endif


static int curr_width, *w_goal;


/*
 *  ORDER AND FACTOR MATRIX
 *
 *  This routine chooses a pivot order for the matrix and factors it
 *  into LU form.  It handles both the initial factorization and subsequent
 *  factorizations when a reordering is desired.  This is handled in a manner
 *  that is transparent to the user.  The routine uses a variation of
 *  Gauss's method where the pivots are associated with L and the
 *  diagonal terms of U are one.
 *
 *  >>> Returned:
 *  The error code is returned.  Possible errors are listed below.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  RHS  <input>  (RealVector)
 *      Representative right-hand side vector that is used to determine
 *      pivoting order when the right hand side vector is sparse.  If
 *      RHS is a NULL pointer then the RHS vector is assumed to
 *      be full and it is not used when determining the pivoting
 *      order.
 *  RelThreshold  <input>  (RealNumber)
 *      This number determines what the pivot relative threshold will
 *      be.  It should be between zero and one.  If it is one then the
 *      pivoting method becomes complete pivoting, which is very slow
 *      and tends to fill up the matrix.  If it is set close to zero
 *      the pivoting method becomes strict Markowitz with no
 *      threshold.  The pivot threshold is used to eliminate pivot
 *      candidates that would cause excessive element growth if they
 *      were used.  Element growth is the cause of roundoff error.
 *      Element growth occurs even in well-conditioned matrices.
 *      Setting the RelThreshold large will reduce element growth and
 *      roundoff error, but setting it too large will cause execution
 *      time to be excessive and will result in a large number of
 *      fill-ins.  If this occurs, accuracy can actually be degraded
 *      because of the large number of operations required on the
 *      matrix due to the large number of fill-ins.  A good value seems
 *      to be 0.001.  The default is chosen by giving a value larger
 *      than one or less than or equal to zero.  This value should be
 *      increased and the matrix resolved if growth is found to be
 *      excessive.  Changing the pivot threshold does not improve
 *      performance on matrices where growth is low, as is often the
 *      case with ill-conditioned matrices.  Once a valid threshold is
 *      given, it becomes the new default.  The default value of
 *      RelThreshold was choosen for use with nearly diagonally
 *      dominant matrices such as node- and modified-node admittance
 *      matrices.  For these matrices it is usually best to use
 *      diagonal pivoting.  For matrices without a strong diagonal, it
 *      is usually best to use a larger threshold, such as 0.01 or
 *      0.1.
 *  AbsThreshold  <input>  (RealNumber)
 *      The absolute magnitude an element must have to be considered
 *      as a pivot candidate, except as a last resort.  This number
 *      should be set significantly smaller than the smallest diagonal
 *      element that is is expected to be placed in the matrix.  If
 *      there is no reasonable prediction for the lower bound on these
 *      elements, then AbsThreshold should be set to zero.
 *      AbsThreshold is used to reduce the possibility of choosing as a
 *      pivot an element that has suffered heavy cancellation and as a
 *      result mainly consists of roundoff error.  Once a valid
 *      threshold is given, it becomes the new default.
 *  DiagPivoting  <input>  (BOOLEAN)
 *      A flag indicating that pivot selection should be confined to the
 *      diagonal if possible.  If DiagPivoting is nonzero and if
 *      DIAGONAL_PIVOTING is enabled pivots will be chosen only from
 *      the diagonal unless there are no diagonal elements that satisfy
 *      the threshold criteria.  Otherwise, the entire reduced
 *      submatrix is searched when looking for a pivot.  The diagonal
 *      pivoting in Sparse is efficient and well refined, while the
 *      off-diagonal pivoting is not.  For symmetric and near symmetric
 *      matrices, it is best to use diagonal pivoting because it
 *      results in the best performance when reordering the matrix and
 *      when factoring the matrix without ordering.  If there is a
 *      considerable amount of nonsymmetry in the matrix, then
 *      off-diagonal pivoting may result in a better equation ordering
 *      simply because there are more pivot candidates to choose from.
 *      A better ordering results in faster subsequent factorizations.
 *      However, the initial pivot selection process takes considerably
 *      longer for off-diagonal pivoting.
 *
 *  >>> Local variables:
 *  pPivot  (ElementPtr)
 *      Pointer to the element being used as a pivot.
 *  ReorderingRequired  (BOOLEAN)
 *      Flag that indicates whether reordering is required.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 *  spSINGULAR
 *  spSMALL_PIVOT
 *  Error is cleared in this function.
 */

int
spOrderAndFactor( char* eMatrix, RealNumber RHS[], RealNumber RelThreshold, RealNumber AbsThreshold, BOOLEAN DiagPivoting, int FirstTime )
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
ElementPtr  pPivot;
int  Step, N_Col, Size, ReorderingRequired;
RealNumber LargestInCol;
int num_elem, num_fill;
ElementPtr  pElement;
    int i, *dist, num_el, del, p_acct, p_accum, row1;
    double f_goal, emax, emaxinv;
#ifdef WRITE_MAT
int N_Row;
FILE *fmat;
#endif

/* Begin `spOrderAndFactor'. */
    ASSERT( IS_VALID(Matrix) AND NOT Matrix->Factored);

/*
    printf ("In spOrderAndFactor\n");
    if (FirstTime == YES)
      printf ("Matrix size = %d\n",Matrix->Size);
*/

    if (Matrix->Format != FORMAT_SPARSE)
      spExpandFormat (Matrix);
    Matrix->Error = spOKAY;
    Matrix->Format = FORMAT_SPARSE;
    Matrix->Updated = 0;
    Matrix->NeedsScale = YES; 
    Size = Matrix->Size;

#ifdef WRITE_MAT
    fmat = fopen ("a.mat","w");
    fprintf (fmat,"%d\n",Size);
    for (N_Row=1 ; N_Row <= Size; N_Row++) {
       pElement = Matrix->FirstInRow[N_Row];
       while (pElement != NULL) {
         if (pElement->Fillin == 0) {
           fprintf (fmat,"%d %d %.12g\n",Matrix->IntToExtRowMap[pElement->Row],
              Matrix->IntToExtColMap[pElement->Col],pElement->Real);
         }
         pElement = pElement->NextInRow;
       }
     }
     fclose(fmat);
#endif

    if (Matrix->Pivots_d < Size) {
      if (Matrix->Pivots_d == 0) {
        Matrix->Pivots = (double *) tmalloc((Size+1)*sizeof(double));
      }
      else {
        Matrix->Pivots = (double *) trealloc(Matrix->Pivots, (Size+1)*sizeof(double));
      }
      Matrix->Pivots_d = Size;
    }
    if (RelThreshold <= 0.0) RelThreshold = Matrix->RelThreshold;
    if (RelThreshold > 1.0) RelThreshold = Matrix->RelThreshold;
    if (FirstTime == YES)
      Matrix->RelThreshold = RelThreshold;
    if (AbsThreshold < 0.0) AbsThreshold = Matrix->AbsThreshold;
    if (FirstTime == YES)
      Matrix->AbsThreshold = AbsThreshold;
    ReorderingRequired = NO;
    Matrix->has_scale_factors = 0;

#if REORDER_SCALING
    if (Matrix->NeedsScale) {
      if (NOT Matrix->RowsLinked)
          spcLinkRows( Matrix );
      Matrix->NeedsScale = NO;
      if (Matrix->scale_factors_d < Size) {
        Matrix->scale_factors_d = Size;
        Matrix->row_scale_factors = (double *) trealloc(Matrix->row_scale_factors,
                                (Matrix->scale_factors_d+1)*sizeof(double));
        Matrix->col_scale_factors = (double *) trealloc(Matrix->col_scale_factors,
                                (Matrix->scale_factors_d+1)*sizeof(double));
      }
      for (i= 1 ;i<= Size ; i++) {
        pElement = Matrix->FirstInCol[i];
        if (!pElement) {
#ifdef CHILE
          fprintf (stderr,"Fatal error in factorization, unpopulated column number %d (%s) in matrix\n",
                   i,CKTnodName(Matrix->Ckt,Matrix->IntToExtColMap[i]));
#endif
          return MatrixIsSingular( Matrix, i );
        }
        emax = fabs(pElement->Real);
        while ((pElement = pElement->NextInCol)) {
          if (fabs(pElement->Real) > emax) {
            emax = fabs(pElement->Real);
          }
        }
        if (emax == 0) {
/* not needed, but shows how to print node corresponding to column:
          printf ("zero Column found in scaling\n");
          printf ("Node is: %s\n",CKTnodName(Matrix->Ckt,Matrix->IntToExtColMap[i]));
*/
          return MatrixIsSingular( Matrix, i );
        }
        emaxinv = 1./emax;
        pElement = Matrix->FirstInCol[i];
        while (pElement) {
          pElement->Real *= emaxinv;
          pElement = pElement->NextInCol;
        }
        Matrix->col_scale_factors[Matrix->IntToExtColMap[i]] = emaxinv;
      }
      for (i= 1 ;i<= Size ; i++) {
        pElement = Matrix->FirstInRow[i];
        emax = fabs(pElement->Real);
        while ((pElement = pElement->NextInRow)) {
          if (fabs(pElement->Real) > emax) {
            emax = fabs(pElement->Real);
          }
        }
        if (emax == 0) {
          printf ("zero Row found in scaling\n");
          return MatrixIsSingular( Matrix, i );
        }
        emaxinv = 1./emax;
        pElement = Matrix->FirstInRow[i];
        while (pElement) {
          pElement->Real *= emaxinv;
          pElement = pElement->NextInRow;
        }
        Matrix->row_scale_factors[Matrix->IntToExtRowMap[i]] = emaxinv;
      }
      Matrix->has_scale_factors = 1;
    }
#endif

    if (NOT Matrix->NeedsOrdering) {
/* Matrix has been factored before and reordering is not required. */
        for (Step = 1; Step <= Size; Step++)
        {   pPivot = Matrix->Diag[Step];
            LargestInCol = FindLargestInCol(pPivot->NextInCol);
            if ((LargestInCol * RelThreshold < ELEMENT_MAG(pPivot))) {
                Matrix->Pivots[Step] = pPivot->Real;
                if (Matrix->Complex)
                    ComplexRowColElimination( Matrix, pPivot, Step );
                else
                    RealRowColElimination( Matrix, pPivot, Step );
            }
            else {
                ReorderingRequired = YES;
                break; /* for loop */
            }
        }
        if (NOT ReorderingRequired)
            goto Done;
        else
        {
/*
 * A pivot was not large enough to maintain accuracy,
 * so a partial reordering is required.
 */

#if (ANNOTATE >= ON_STRANGE_BEHAVIOR)
            printf("Reordering,  Step = %1d\n", Step);
#endif
        }
    } /* End of if(NOT Matrix->NeedsOrdering) */
    else
    {
/*
 * This is the first time the matrix has been factored.  These few statements
 * indicate to the rest of the code that a full reodering is required rather
 * than a partial reordering, which occurs during a failure of a fast
 * factorization.
 */
        Step = 1;
        if (NOT Matrix->RowsLinked)
            spcLinkRows( Matrix );
        if (NOT Matrix->InternalVectorsAllocated)
            spcCreateInternalVectors( Matrix );
        if (Matrix->Error >= spFATAL)
            return Matrix->Error;
    }

/* First step is to eliminate Fillin elements as these tend to
   accumulate when multiple reorderings are performed */

#ifdef DEBUG
    num_elem = 0;
    num_fill = 0;
    for (N_Col=1 ; N_Col <= Size; N_Col++) {
       pElement = Matrix->FirstInCol[N_Col];
       while (pElement != NULL)
       {
           if (pElement->Fillin > 0) {
             num_fill++;
           }
           else {
             num_elem++;
           }
           pElement = pElement->NextInCol;
       }
    }
    printf ("Before StripFills, elems = %d, fills = %d, Step = %d\n",
            num_elem,num_fill,Step);
#endif

    spStripFills((char *)Matrix, Step);
    memset (Matrix->Intermediate4, 0, (Size+1)*sizeof(double));
#ifdef DEBUG
    num_elem = 0;
    num_fill = 0;
    for (N_Col=1 ; N_Col <= Size; N_Col++) {
       pElement = Matrix->FirstInCol[N_Col];
       while (pElement != NULL)
       {
           if (pElement->Fillin > 0)
             num_fill++;
           else
             num_elem++;
           pElement = pElement->NextInCol;
       }
    }
    printf ("After StripFills, elems = %d, fills %d\n",num_elem,num_fill);
#endif

/* Form initial Markowitz products. */
    CountMarkowitz( Matrix, RHS, Step );
    MarkowitzProducts( Matrix, Step );
    Matrix->MaxRowCountInLowerTri = -1;

/* Perform reordering and factorization. */
    for (; Step <= Size; Step++) {
        pPivot = SearchForPivot( Matrix, Step, DiagPivoting );
/*
        printf ("At step = %d, pivot method = %c, pivot = %g\n",Step,Matrix->PivotSelectionMethod, pPivot->Real);
*/
        if (pPivot == NULL) {
          printf ("Null pivot returned from Search\n");
          return MatrixIsSingular( Matrix, Step );
        }
        Matrix->Pivots[Step] = pPivot->Real;
        Matrix->Intermediate4[pPivot->Row] = 0;
        ExchangeRowsAndCols( Matrix, pPivot, Step );
        pElement = pPivot;
        while(pElement) {
          Matrix->Intermediate4[pElement->Col] = 0;
          pElement = pElement->NextInRow;
        }
        if (Step > 1) {
          pElement = Matrix->Diag[Step-1];
          while(pElement) {
            if (pElement->Real == Matrix->Intermediate4[pElement->Col])
              Matrix->Intermediate4[pElement->Col] = 0;
            pElement = pElement->NextInRow;
          }
        }

        if (Matrix->Complex)
            ComplexRowColElimination( Matrix, pPivot, Step );
        else
            RealRowColElimination( Matrix, pPivot, Step );

        if (Matrix->Error >= spFATAL) {
          printf ("Error returned from elimination\n");
          return Matrix->Error;
        }
        UpdateMarkowitzNumbers( Matrix, pPivot );

#ifdef DEBUG
        if (Step%100 == 0) printf ("At: %d/%d  %d/%d\n",Step,Size,Matrix->Elements, Matrix->Fillins);
#endif
#if (ANNOTATE == FULL)
        WriteStatus( Matrix, Step );
#endif
    }
    num_elem = 0;
    num_fill = 0;
    for (N_Col=1 ; N_Col <= Size; N_Col++) {
       pElement = Matrix->FirstInCol[N_Col];
       while (pElement != NULL)
       {
           if (pElement->Fillin <= 0)
             num_fill++;
           else
             num_elem++;
           pElement = pElement->NextInCol;
       }
    }
#ifdef DEBUG
    printf ("At end of spOrderAndFactor, elems = %d, fills %d\n",num_elem,num_fill);
#endif

Done:
    Matrix->NeedsOrdering = NO;
    Matrix->Reordered = YES;
    Matrix->Factored = YES;
    Matrix->Format = FORMAT_SPARSE;
    Matrix->DensePointers = 0;
    Matrix->Partitioned = NO;
    Matrix->OverflowDanger = 0;
/*
    spPrint(eMatrix,1,1,1);
*/
    return Matrix->Error;
}

void spExpand (MatrixPtr Matrix)
{
    return;
}

void spExpandFormat (MatrixPtr Matrix)
{
    ElementPtr pElement;
    int i, my_p;
    struct context_m *my_context;

    if (Matrix->Format == FORMAT_SPARSE) {
      return;
    }
    if (Matrix->Format == FORMAT_DENSE) {
      my_context = &Matrix->MyStuff[0];
      for (i=1 ; i<=Matrix->Size ; i++) {
        pElement = Matrix->FirstInCol[i];
        my_p = my_context->ColStart_s[i];
        while (pElement) {
          pElement->Real = my_context->MyD[my_p++];
          pElement = pElement->NextInCol;
        }
      }
      Matrix->Format = FORMAT_SPARSE;
      return;
    }
    fprintf (stderr, "Internal error: Unknown factored format\n");
    return;
}

/* static void *BAD_ADD= (void *) 0x405d2acc80; */
void check_mat(MatrixPtr Matrix, int *num, int tag)
{
    ElementPtr pElement;
    int i, Size;
    return;
}


/*
 *  FACTOR MATRIX
 *
 *  This routine is the companion routine to spOrderAndFactor().
 *  Unlike spOrderAndFactor(), spFactor() cannot change the ordering.
 *  It is also faster than spOrderAndFactor().  The standard way of
 *  using these two routines is to first use spOrderAndFactor() for the
 *  initial factorization.  For subsequent factorizations, spFactor()
 *  is used if there is some assurance that little growth will occur
 *  (say for example, that the matrix is diagonally dominant).  If
 *  spFactor() is called for the initial factorization of the matrix,
 *  then spOrderAndFactor() is automatically called with the default
 *  threshold.  This routine uses "row at a time" LU factorization.
 *  Pivots are associated with the lower triangular matrix and the
 *  diagonals of the upper triangular matrix are ones.
 *
 *  >>> Returned:
 *  The error code is returned.  Possible errors are listed below.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 *  spSINGULAR
 *  spZERO_DIAG
 *  spSMALL_PIVOT
 *  Error is cleared in this function.
 */

int
spFactor( char *eMatrix )
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;

#if spCOMPLEX
    if (Matrix->Complex) {
      return(FactorComplexMatrix( Matrix ));
    }
#endif
    return (spFactorAndSolve(eMatrix, (double *) NULL));
}
   


int spFactorAndSolve(char *eMatrix, double *RHS)

{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
ElementPtr  pElement, pColumn;
int  Step, Size, i, ii, k_start, k_end, row_lo, j, k, n, m;
int lo_row, lo_col, hi_col, hi_lim, c_c, hi_col_lim;
double time_interval;
long tstart,tend;
int *work, n_lo, n_hi, num_bins, coverage, my_width, my_ops;

int num_ent;
int row_step, ind, load_bal, n_max, n_min;
int load_bal_done;
double x_max, x_avg, x_min;
int low_ind, wait, shift;
long t1, t2, time_start, time_end;
int first_pass;
int num_chunks, stat, mid_pe, d_avg;
long c_sum, target, local_target, total_clocks, nm, nm_goal;
double avgpiv_ratio, minpiv, minpiv_ratio, piv, del_d;
double t_chunk, t_goal, t_accum, t_total, t_rat, t_avg, t_max, t_min_p;
int r_chunk, row_curr, n_min_p, ni, num_ind, n_d, me;
int c_up, s_num, repartitioned, num_strips, update, my_p, need_solve, num_report;
long goal, tot, tot1, tot2, num_s, curr;
#ifdef CHECK_ANS
static double *Dest;
static int Dest_d;
#endif

ElementPtr *f_col;
static double *my_rhs;
int use_scr = 0;

static int *valid_f, *valid_b;
static int num_call, my_rhs_d;
  int num_els, num_fills,jj;

double *scr, *value, diag_value;
int *index, *cs, *cd, my_i, n_ind, n_ind_i, last_my_i;
int *ind_list_i, numi;
unsigned char *ind_list;
unsigned char n_char, maxchar=0xff;

long min_rd,max_rd;

struct context_m *my_context;
struct strip_out *strip, *old_strip, *my_strip, *prev_strip;

/* spPrint(eMatrix,1,1,1); */
/*   printf ("RHS: %g %g\n",RHS[1], RHS[2]); */
#ifdef DEBUG
    printf ("In spFactorAndSolve\n");
#endif

#if spCOMPLEX
      if (Matrix->Complex) {
        printf ("Encountered factor/solve call for Complex\n");
      }
#endif

    Matrix->has_scale_factors = 0;
    Size = Matrix->Size;
#ifdef CHECK_ANS
    if (Dest_d < Size) {
      Dest_d = Size;
      Dest = trealloc(Dest, (Dest_d+1)*sizeof(double));
    }
#endif
    if (my_rhs_d < Size) {
      my_rhs_d = Size;
      my_rhs = (double *) tmalloc((my_rhs_d+1)*sizeof(double));
    }
    spSetReal( (char *)Matrix );
    
    if (Matrix->NewFlags || !Matrix->Updated) {
      update = 1;
    }
    else {
      update = 0;
    }
/*
    printf ("PE: %d, update = %d :: %d %d\n",pe_in_smp,update,Matrix->NewFlags,Matrix->Updated);
*/
    my_context = &Matrix->MyStuff[0];

    minpiv = Matrix->AbsThreshold + 1;
    minpiv_ratio = 1.;
    avgpiv_ratio = 0.;
/* Begin `spFactor'. */
    ASSERT( IS_VALID(Matrix) AND NOT Matrix->Factored);

    repartitioned = 0;
    if (my_context->Dsize < Size) {
      if (my_context->Dsize == 0) {
        memset(my_context,0,sizeof(struct context_m));
      }
      my_context->Dsize = Size;
      my_context->ColStart_s = (int *) trealloc(my_context->ColStart_s,
                               (my_context->Dsize+2)*sizeof(int));
      my_context->ColStart[0] = my_context->ColStart_s;
      my_context->ColDiag = (int *) trealloc(my_context->ColDiag,
                               (my_context->Dsize+2)*sizeof(int));
      my_context->Dest = (double *) trealloc(my_context->Dest,
                               (my_context->Dsize+1)*sizeof(double));
    }
      Matrix->Error = spOKAY;
      if (Matrix->NeedsOrdering || Matrix->OverflowDanger>=OF_THRESHOLD || Matrix->SmallTimeStep) {

#ifdef DEBUG_OVER
        printf ("In spFactorAndSolve, finding new pivots (with format = %d) due to ",
                Matrix->Format);
        if (Matrix->NeedsOrdering) {
          printf ("new ordering required"); 
        }
        else if (Matrix->OverflowDanger) {
          printf ("small pivot"); 
        }
        else if (Matrix->SmallTimeStep) {
          printf ("small time step");
        }
        printf ("\n");
#endif

/*    check_mat(Matrix, &num_call, 4); */
        stat = spOrderAndFactor( eMatrix, (RealVector)NULL,
             Matrix->RelThreshold, Matrix->AbsThreshold,
             DIAG_PIVOTING_AS_DEFAULT, 0 );
        Matrix->Max_TS = 0;
        Matrix->OverflowDanger = 0;
        Matrix->SmallTimeStep = 0;
        Matrix->NeedsOrdering = 0;
        if (stat == spOKAY) {
          if (RHS) {
            stat = spSolve (eMatrix, RHS, RHS, NULL, NULL );
          }
        }
        return stat;
      }

/*
        num_els = 0;
        num_fills = 0;
        for (i=1 ; i<=Size ; i++) {
          pElement = Matrix->FirstInCol[i];
          while (pElement) {
            if (pElement->Fillin)
              num_fills++;
            else
              num_els++;
            pElement = pElement->NextInCol;
          }
        }
        if (num_els+num_fills != Matrix->Elements + Matrix->Fillins)
          printf ("Element count: %d/%d  Fill count: %d/%d\n",num_els,
              Matrix->Elements, num_fills, Matrix->Fillins);
        my_context->BufUsed = num_els+num_fills;
*/
        my_context->BufUsed = Matrix->Elements + Matrix->Fillins;
        if (my_context->BufUsed > my_context->BufDim) {
          my_context->BufDim = my_context->BufUsed;
          my_context->MyD = (double *) trealloc(my_context->MyD, my_context->BufDim*sizeof(int)+
                                    my_context->BufDim*sizeof(double));
          my_context->MyI = (int *) (my_context->MyD + my_context->BufDim);
        }

      if (Matrix->Format != FORMAT_SPARSE)
        Matrix->Diag[1]->Real = my_context->MyD[0];
      if (minpiv > fabs(Matrix->Diag[1]->Real)) {
        minpiv = fabs(Matrix->Diag[1]->Real);
      }
      if (Matrix->Diag[1]->Real == 0.0) {
        printf ("Zero pivot in spFactor with format = %d\n",Matrix->Format);
        Matrix->OverflowDanger = 2*OF_THRESHOLD;
        return ZeroPivot( Matrix, 1 );
      }
      if (minpiv_ratio*fabs(Matrix->Pivots[1]) > fabs(Matrix->Diag[1]->Real)) {
        minpiv_ratio = fabs(Matrix->Diag[1]->Real/Matrix->Pivots[1]);
      }
      if (fabs(Matrix->Pivots[1]) < fabs(Matrix->Diag[1]->Real)) {
        avgpiv_ratio += fabs(Matrix->Diag[1]->Real/Matrix->Pivots[1]);
      }
      else {
        avgpiv_ratio += 1;
      }
      Matrix->Diag[1]->Real = 1.0 / Matrix->Diag[1]->Real;
      if (Matrix->Format != FORMAT_SPARSE)
        Matrix->MyStuff->MyD[0] = Matrix->Diag[1]->Real;

#ifdef CHECK_ANS
      make_copy(Matrix);
      Matrix->Factored = YES;
      Matrix->Format = FORMAT_SPARSE;

      for (Step = 2; Step <= Size; Step++) {
/* Scatter. */
        pElement = Matrix->FirstInCol[Step];
        while (pElement != NULL)
        {   Dest[pElement->Row] = pElement->Real;
            pElement = pElement->NextInCol;
        }

/* Update column. */
        pColumn = Matrix->FirstInCol[Step];
        while (pColumn->Row < Step)
        {   pElement = Matrix->Diag[pColumn->Row];
            pColumn->Real = Dest[pColumn->Row] * pElement->Real;
            while ((pElement = pElement->NextInCol) != NULL) {
                Dest[pElement->Row] -= pColumn->Real * pElement->Real;
            }
            pColumn = pColumn->NextInCol;
        }

/* Gather. */
        pElement = Matrix->Diag[Step]->NextInCol;
        while (pElement != NULL)
        {   pElement->Real = Dest[pElement->Row];
            pElement = pElement->NextInCol;
        }

/* Check for singular matrix. */
        if (Dest[Step] == 0.0) {
#ifdef CHILE
          printf ("Zero pivot at Step = %d,(%s)\n", Step, CKTnodName(Matrix->Ckt,Matrix->IntToExtColMap[Step]));
#endif
          return ZeroPivot( Matrix, Step );
        }
        if (minpiv > fabs(Dest[Step]))
            minpiv = fabs(Dest[Step]);
        if (minpiv_ratio*fabs(Matrix->Pivots[Step]) > fabs(Dest[Step]))
            minpiv_ratio = fabs(Dest[Step]/Matrix->Pivots[Step]);
        if (fabs(Matrix->Pivots[Step]) < fabs(Dest[Step])) {
          avgpiv_ratio += fabs(Dest[Step]/Matrix->Pivots[Step]);
        }
        else {
          avgpiv_ratio += 1;
        }
        Matrix->Diag[Step]->Real = 1.0 / Dest[Step];
      }

      if (RHS) {
        memcpy (Matrix->Intermediate4, RHS, (Size+1)*sizeof(double));
        stat = spSolve (eMatrix, Matrix->Intermediate4, Matrix->Intermediate4, NULL, NULL );
      }
      swap_copy(Matrix);
#endif /* CHECK_ANS */

/* Start serial factorization. */
        scr = my_context->Dest;
        value =  my_context->MyD;
        index = my_context->MyI;
        cd = my_context->ColDiag;
        cs = my_context->ColStart_s;
        if (update)
          Matrix->Updated = 1;
        if (RHS)
          need_solve = 1;
        pElement = Matrix->FirstInCol[1];
        my_p = 0;
        if (update) {
          cs[1] = my_p;
          while (pElement != NULL) {
            if (pElement->Row == 1)
              cd[1] = my_p;
#ifdef CHILE
            pElement->RealDense = &value[my_p];
#endif
            value[my_p] = pElement->Real;
            index[my_p++] = pElement->Row;
            pElement = pElement->NextInCol;
          }
          cs[2] = my_p;
        }
        else {
          if (Matrix->Format == FORMAT_SPARSE) {
            while (pElement != NULL) {
              value[my_p++] = pElement->Real;
              pElement = pElement->NextInCol;
            }
          }
        }
        n_ind = 0;
        n_ind_i = 0;
        ind_list = my_context->ind_list;
        ind_list_i = my_context->ind_list_i;
        for (Step = 2; Step <= Size; Step++) {
/* Scatter. */
          if (Matrix->Format == FORMAT_SPARSE) {
            pElement = Matrix->FirstInCol[Step];
            my_p = cs[Step];
            while (pElement != NULL)
            {   if (update) {
                  scr[pElement->Row] = pElement->Real;
#ifdef CHILE
                  pElement->RealDense = &value[my_p];
#endif
                  index[my_p] = pElement->Row;
                  if (Step == pElement->Row)
                    cd[Step] = my_p;
                }
                else {
                  if (use_scr)
                    scr[pElement->Row] = pElement->Real;
                  else
                    value[my_p] = pElement->Real;
                }
                my_p++;
                pElement = pElement->NextInCol;
            }
          }
          else {
            if (use_scr) {
              for (i=cs[Step] ; i<cs[Step+1] ; i++) {
                scr[index[i]] = value[i];
              }
            }
          }

/* Update column. */
          if (update) {
            cs[Step+1] = my_p;
            for (i=cs[Step] ; i<cd[Step] ; i++) {
              piv = scr[index[i]] * value[cd[index[i]]];
              value[i] = piv;
              my_i = i;
              last_my_i = my_i;
              if (my_context->ind_list_d < n_ind+Size) {
                my_context->ind_list_d = 2*n_ind+Size;
                my_context->ind_list = (unsigned char *) trealloc(my_context->ind_list,
                    my_context->ind_list_d);
                ind_list = my_context->ind_list;
              }
              if (my_context->ind_list_i_d < n_ind_i+Size) {
                my_context->ind_list_i_d = 2*n_ind_i+Size;
                my_context->ind_list_i = (int *) trealloc(my_context->ind_list_i,
                    my_context->ind_list_i_d*sizeof(int));
                ind_list_i = my_context->ind_list_i;
              }
              for (j=cd[index[i]]+1 ; j<cs[index[i]+1] ; j++) {
                scr[index[j]] -= piv * value[j];
                while (index[++my_i] < index[j]) { }
                n = my_i-last_my_i;
                if (n <= maxchar) {
                  ind_list[n_ind++] = (unsigned char) my_i-last_my_i;
                }
                else {
                  ind_list[n_ind++] = 0;
                  ind_list_i[n_ind_i++] = n;
                }
                last_my_i = my_i;
              }
            }
            diag_value = scr[Step];
          }
          else {
            if (use_scr) {
              for (i=cs[Step] ; i<cd[Step] ; i++) {
                piv = scr[index[i]] * value[cd[index[i]]];
                value[i] = piv;
                for (j=cd[index[i]]+1 ; j<cs[index[i]+1] ; j++) {
                  scr[index[j]] -= piv * value[j];
                }
              }
              diag_value = scr[Step];
            }
            else {
              for (i=cs[Step] ; i<cd[Step] ; i++) {
                j = cd[index[i]];
                piv = value[i] * value[j++];
                if (fabs(piv) > HUGE) {
                  printf ("excessively big pivot at Step = %d\n", Step);
                  Matrix->OverflowDanger = 2*OF_THRESHOLD;
                  return ZeroPivot( Matrix, Step );
                }
                value[i] = piv;
                last_my_i = i;
                for ( ; j<cs[index[i]+1] ; j++) {
                  n_char = ind_list[n_ind++];
                  last_my_i += (int) n_char;
                  if (n_char == 0)
                    last_my_i += ind_list_i[n_ind_i++];
                  value[last_my_i] -= piv * value[j];
                }
              }
              diag_value = value[cd[Step]];
            }
          }

/* Check for singular matrix. */
          if (fabs(diag_value) < 0.001*Matrix->AbsThreshold) {
#ifdef CHILE
            printf ("Zero pivot at Step = %d,(%s)\n", Step, CKTnodName(Matrix->Ckt,Matrix->IntToExtColMap[Step]));
#endif
            Matrix->OverflowDanger = 2*OF_THRESHOLD;
            return ZeroPivot( Matrix, Step );
          }
          if (minpiv > fabs(diag_value))
              minpiv = fabs(diag_value);
          if (minpiv_ratio*fabs(Matrix->Pivots[Step]) > fabs(diag_value))
              minpiv_ratio = fabs(diag_value/Matrix->Pivots[Step]);
          if (fabs(Matrix->Pivots[Step]) < fabs(diag_value)) {
            avgpiv_ratio += fabs(diag_value/Matrix->Pivots[Step]);
          }
          else {
            avgpiv_ratio += 1;
          }
          value[cd[Step]] = 1.0 / diag_value;

/* Gather. */
          if (update | use_scr) {
            for (i=cd[Step]+1 ; i<cs[Step+1] ; i++)
              value[i] = scr[index[i]];
          }

          if (!RHS) {
            pElement = Matrix->FirstInCol[Step];
            my_p = cs[Step];
            while (pElement) {
              pElement->Real = value[my_p++];
              pElement = pElement->NextInCol;
            }
          }
        }
/*
void *buf;
int ndimi, ndimd, ndimc;
        if (update) {
          ndimi = (my_context->Bufdim+4-1)/4;
          ndimd = my_context->Bufdim;
          ndimc = (n_ind+8-1)/8;
          buf = tmalloc(ndimi*sizeof(int)+ndimd*sizeof(double)+ndimc+n_ind_i*sizeof(int));
          memcpy (buf, my_context->MyI, my_context->Bufdim*sizeof(int));
          tfree(my_context->MyI);
          my_context->MyI = (int *) buf;
        }
*/
        if (RHS) {
          need_solve = 0;
          for (i=1 ; i<=Size ; i++ )
            my_rhs[i] = RHS[Matrix->IntToExtRowMap[i]];
          for (Step=1 ; Step<=Size ; Step++ ) {
            piv = value[cd[Step]];
            if (piv != 0) {
              my_rhs[Step] = (piv *= my_rhs[Step]);;
              for (j=cd[Step]+1 ; j<cs[Step+1] ; j++) {
                my_rhs[index[j]] -= piv*value[j];
              }
            }
          }
          for (Step=Size ; Step>0 ; Step-- ) {
            piv = my_rhs[Step];
            for (j=cs[Step] ; j<cd[Step] ; j++) {
              my_rhs[index[j]] -= piv*value[j];
            }
          }
          for (i=1 ; i<=Size ; i++ )
            RHS[Matrix->IntToExtColMap[i]] = my_rhs[i];
          Matrix->Format = FORMAT_DENSE;
        }
        else {
          spExpandFormat(Matrix);
        }
        Matrix->Factored = YES;
/*
        min_rd = (long) Matrix->FirstInCol[1]->RealDense;
        max_rd = (long) Matrix->FirstInCol[1]->RealDense;
        for (Step = 1; Step <= Size; Step++) {
          pElement = Matrix->FirstInCol[Step];
          while(pElement) {
            if (min_rd > (long) pElement->RealDense)
              min_rd = (long) pElement->RealDense;
            if (max_rd < (long) pElement->RealDense)
              max_rd = (long) pElement->RealDense;
            pElement = pElement->NextInCol;
          }
        }
        printf ("Min,Max RealDense: %lx, %lx\n",min_rd,max_rd);
*/
          
      avgpiv_ratio /= Size;
      Matrix->OverflowDanger = Matrix->OverflowDanger*19/20;
#ifdef DEBUG_OVER
      if (Matrix->OverflowDanger > 0)
        printf ("revised Matrix->OverflowDanger = %d\n",Matrix->OverflowDanger);
#endif
      if (minpiv < Matrix->AbsThreshold) {
        del_d = OF_THRESHOLD;
      }
      else {
        del_d = OF_THRESHOLD*Matrix->AbsThreshold/minpiv;
      }
      del_d /= DEL_RED;
      Matrix->OverflowDanger += del_d;
#ifdef DEBUG_OVER
      if (Matrix->OverflowDanger > 0)
        printf ("After Matrix->AbsThreshold test, revised Matrix->OverflowDanger = %d\n",
               Matrix->OverflowDanger);
#endif
/* These are only used when into transient analysis, state is set in CKTload */
      if (Matrix->PivotRefining) {
/*
        if (avgpiv_ratio > 10000) {
          del_d = avgpiv_ratio/1000;
          if (del_d > OF_THRESHOLD) del_d = OF_THRESHOLD*DEL_RED;
*/
        if (avgpiv_ratio > 100) {
          del_d = avgpiv_ratio;
          if (del_d > OF_THRESHOLD) del_d = OF_THRESHOLD;
          del_d /= DEL_RED;
          Matrix->OverflowDanger += del_d;
#ifdef DEBUG_OVER
          if (Matrix->OverflowDanger > 0)
            printf ("After avgpiv_ratio test, revised Matrix->OverflowDanger = %d\n",
                   Matrix->OverflowDanger);
#endif
        }
        if (minpiv_ratio < 1.e-8 && minpiv_ratio > 0) {
          del_d = 1.e-8/minpiv_ratio;
          if (del_d > OF_THRESHOLD) del_d = OF_THRESHOLD*DEL_RED;
          del_d /= DEL_RED;
          Matrix->OverflowDanger += del_d;
#ifdef DEBUG_OVER
          if (Matrix->OverflowDanger > 0)
            printf ("After minpiv_ratio test, revised Matrix->OverflowDanger = %d\n",
                   Matrix->OverflowDanger);
#endif
        }
      }

#ifdef CHECK_IND
    spCheckInd(Matrix, "spFactor bottom");
#endif
#ifdef TIMERS
    timer_stop(10);
    timer_wait1 = wait1*(*timer_calibration);
    timer_wait2 = wait2*(*timer_calibration);
    timer_wait3 = wait3*(*timer_calibration);
    timer_wait4 = wait4*(*timer_calibration);
    total_time = timer_get(10);
    timer_wait1 *= 100./total_time;
    timer_wait2 *= 100./total_time;
    timer_wait3 *= 100./total_time;
    timer_wait4 *= 100./total_time;
    printf ("PE: %d spFactorAndSolve time: %.4f, Wait times: factor: %5.2f%% forward: %5.2f%% back: %5.2f%%\n",
            pe_in_smp, total_time, timer_wait2, timer_wait3, timer_wait4);
#endif

    { if (need_solve)
        stat = spSolve (eMatrix, RHS, RHS, NULL, NULL );
    }

#ifdef CHECK_ANS
      if (RHS) {
        printf ("Check of RHS answers:");
        num_report = 0;
        for (i=Size ; i>=1 ; i--) {
          if (Matrix->Intermediate4[i] != RHS[i]) {
            if (num_report++ < 10)
              printf ("\nDifferent RHS at i = %d,  %g :: %g\n",i,Matrix->Intermediate4[i],RHS[i]);
            RHS[i] = Matrix->Intermediate4[i];
          }
        }
        if (num_report == 0)
          printf (" Identical\n");
        else
          printf ("\nTotal of %d differences found\n",num_report);
      }
      num_report = 0;
      printf ("Check of factorized matrix:");

          me = 0;
          num_report = 0;
          my_context->ColStart[me] = my_context->ColStart_s;
          f_col = Matrix->FirstInCol;

          for (j=1 ; j<=Size ; j++) {
            pElement = f_col[j];
            if (pElement) {
              for (k=my_context->ColStart[me][j]; k<my_context->ColStart[me][j+1] ; k++) {
                if (my_context->MyI[k] != pElement->Row) {
                  if (num_report++ < 10)
                    printf ("\nError in row in spFactor copy: %d %d :: %d %d",
                      my_context->MyI[k], pElement->Row, j, pElement->Col);
                }
                if (pElement->Row == pElement->Col) {
                  if (k != my_context->ColDiag[j]) {
                    if (num_report++ < 10)
                      printf ("\nError in diag in spFactor copy: %d %d :: %d %d",
                        my_context->MyI[k], pElement->Row, j, pElement->Col);
                  }
                }
                if (pElement->RealCopy != my_context->MyD[k]) {
                  if (num_report++ < 10)
                    printf ("\nMismatch in factor: (%d,%d) %g != %g",pElement->Row,pElement->Col,
                      pElement->RealCopy, my_context->MyD[k]);
                }
                pElement = pElement->NextInCol;
              }
            }
          }
      if (num_report == 0)
        printf (" Identical\n");
      else
        printf ("\nTotal of %d differences found\n",num_report);

#endif /* CHECK_ANS */

#ifdef FINGERPRINT
      FingerPrint (RHS, Matrix->IntToExtRowMap, Matrix->IntToExtColMap, Size);
#endif
/*
    printf ("Exiting spFactorAndSolve, with Matrix->Factored = %d\n",Matrix->Factored);
*/
    return (Matrix->Error);
}






#if spCOMPLEX
/*
 *  FACTOR COMPLEX MATRIX
 *
 *  This routine is the companion routine to spFactor(), it
 *  handles complex matrices.  It is otherwise identical.
 *
 *  >>> Returned:
 *  The error code is returned.  Possible errors are listed below.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *
 *  >>> Possible errors:
 *  spSINGULAR
 *  Error is cleared in this function.
 */

static int
FactorComplexMatrix( MatrixPtr Matrix )
{
ElementPtr  pElement;
ElementPtr  pColumn;
int  Step, Size;
ComplexNumber Mult, Pivot;

/* Begin `FactorComplexMatrix'. */
    ASSERT(Matrix->Complex);

    Size = Matrix->Size;
    pElement = Matrix->Diag[1];
    if (ELEMENT_MAG(pElement) == 0.0) return ZeroPivot( Matrix, 1 );
/* Cmplx expr: *pPivot = 1.0 / *pPivot. */
    CMPLX_RECIPROCAL( *pElement, *pElement );

/* Start factorization. */
    for (Step = 2; Step <= Size; Step++)
    {   if (Matrix->DoCmplxDirect[Step])
        {   /* Update column using direct addressing scatter-gather. */
            ComplexNumber  *Dest;
            Dest = (ComplexNumber *)Matrix->Intermediate;

/* Scatter. */
            pElement = Matrix->FirstInCol[Step];
            while (pElement != NULL)
            {   Dest[pElement->Row] = *(ComplexNumber *)pElement;
                pElement = pElement->NextInCol;
            }

/* Update column. */
            pColumn = Matrix->FirstInCol[Step];
            while (pColumn->Row < Step)
            {   pElement = Matrix->Diag[pColumn->Row];
                /* Cmplx expr: Mult = Dest[pColumn->Row] * (1.0 / *pPivot). */
                CMPLX_MULT(Mult, Dest[pColumn->Row], *pElement);
                CMPLX_ASSIGN(*pColumn, Mult);
                while ((pElement = pElement->NextInCol) != NULL)
                {   /* Cmplx expr: Dest[pElement->Row] -= Mult * pElement */
                    CMPLX_MULT_SUBT_ASSIGN(Dest[pElement->Row],Mult,*pElement);
                }
                pColumn = pColumn->NextInCol;
            }

/* Gather. */
            pElement = Matrix->Diag[Step]->NextInCol;
            while (pElement != NULL)
            {   *(ComplexNumber *)pElement = Dest[pElement->Row];
                pElement = pElement->NextInCol;
            }

/* Check for singular matrix. */
            Pivot = Dest[Step];
            if (CMPLX_1_NORM(Pivot) == 0.0) return ZeroPivot( Matrix, Step );
            CMPLX_RECIPROCAL( *Matrix->Diag[Step], Pivot );  
        }
        else
        {   /* Update column using direct addressing scatter-gather. */
            ComplexNumber  **pDest;
            pDest = (ComplexNumber **)Matrix->Intermediate;

/* Scatter. */
            pElement = Matrix->FirstInCol[Step];
            while (pElement != NULL)
            {   pDest[pElement->Row] = (ComplexNumber *)pElement;
                pElement = pElement->NextInCol;
            }

/* Update column. */
            pColumn = Matrix->FirstInCol[Step];
            while (pColumn->Row < Step)
            {   pElement = Matrix->Diag[pColumn->Row];
                /* Cmplx expr: Mult = *pDest[pColumn->Row] * (1.0 / *pPivot). */
                CMPLX_MULT(Mult, *pDest[pColumn->Row], *pElement);
                CMPLX_ASSIGN(*pDest[pColumn->Row], Mult);
                while ((pElement = pElement->NextInCol) != NULL)
                {  /* Cmplx expr: *pDest[pElement->Row] -= Mult * pElement */
                   CMPLX_MULT_SUBT_ASSIGN(*pDest[pElement->Row],Mult,*pElement);
                }
                pColumn = pColumn->NextInCol;
            }

/* Check for singular matrix. */
            pElement = Matrix->Diag[Step];
            if (ELEMENT_MAG(pElement) == 0.0) return ZeroPivot( Matrix, Step );
            CMPLX_RECIPROCAL( *pElement, *pElement );  
        }
    }

    Matrix->Factored = YES;
    Matrix->Format = FORMAT_SPARSE;
    return (Matrix->Error = spOKAY);
}
#endif /* spCOMPLEX */






/*
 *  PARTITION MATRIX
 *
 *  This routine determines the cost to factor each row using both
 *  direct and indirect addressing and decides, on a row-by-row basis,
 *  which addressing mode is fastest.  This information is used in
 *  spFactor() to speed the factorization.
 *
 *  When factoring a previously ordered matrix using spFactor(), Sparse
 *  operates on a row-at-a-time basis.  For speed, on each step, the
 *  row being updated is copied into a full vector and the operations
 *  are performed on that vector.  This can be done one of two ways,
 *  either using direct addressing or indirect addressing.  Direct
 *  addressing is fastest when the matrix is relatively dense and
 *  indirect addressing is best when the matrix is quite sparse.  The
 *  user selects the type of partition used with Mode.  If Mode is set
 *  to spDIRECT_PARTITION, then the all rows are placed in the direct
 *  addressing partition.  Similarly, if Mode is set to
 *  spINDIRECT_PARTITION, then the all rows are placed in the indirect
 *  addressing partition.  By setting Mode to spAUTO_PARTITION, the
 *  user allows Sparse to select the partition for each row
 *  individually.  spFactor() generally runs faster if Sparse is
 *  allowed to choose its own partitioning, however choosing a
 *  partition is expensive.  The time required to choose a partition is
 *  of the same order of the cost to factor the matrix.  If you plan to
 *  factor a large number of matrices with the same structure, it is
 *  best to let Sparse choose the partition.  Otherwise, you should
 *  choose the partition based on the predicted density of the matrix.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  Mode  <input>  (int)
 *      Mode must be one of three special codes: spDIRECT_PARTITION,
 *      spINDIRECT_PARTITION, or spAUTO_PARTITION.
 */

#define C_SIZE 20

void
spPartition( char* eMatrix, int Mode )
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
ElementPtr  pElement, pColumn;
int  Step, Size;
int  *Nc, *No, *Nm, *Nb;
BOOLEAN *DoRealDirect, *DoCmplxDirect;
int start, stride, i, j, n_diag;
static int *scr, scr_d;
long flops;

/*
   Gather time estimate and carryover columns to solve for each row.  These are
   used in spFactorAndSolve for allocating arrays, and estimating load balance.
   Nc is time estimate
   No is number of carryover columns from above row
*/
/* Begin `spPartition'. */
      ASSERT( IS_SPARSE( Matrix ) );
      if (Matrix->Partitioned) return;
      Size = Matrix->Size;
      DoRealDirect = Matrix->DoRealDirect;
      DoCmplxDirect = Matrix->DoCmplxDirect;
      Matrix->Partitioned = YES;
    return;
}

/*
 *  CREATE INTERNAL VECTORS
 *
 *  Creates the Markowitz and Intermediate vectors.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *
 *  >>> Local variables:
 *  SizePlusOne  (unsigned)
 *      Size of the arrays to be allocated.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

void
spcCreateInternalVectors( MatrixPtr Matrix )
{
int  Size;

/* Begin `spcCreateInternalVectors'. */
/* Create Markowitz arrays. */
    Size= Matrix->Size;

    if (Matrix->MarkowitzRow == NULL)
    {   if (( Matrix->MarkowitzRow = ALLOC(int, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
    if (Matrix->MarkowitzCol == NULL)
    {   if (( Matrix->MarkowitzCol = ALLOC(int, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
    if (Matrix->MarkowitzProd == NULL)
    {
        if (( Matrix->MarkowitzProd = ALLOC(long, Size+2)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
    if (Matrix->Nc == NULL)
    {   if (( Matrix->Nc = ALLOC(int, Size+2)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
    if (Matrix->Nm == NULL)
    {   if (( Matrix->Nm = ALLOC(int, Size+2)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
    if (Matrix->No == NULL)
    {   if (( Matrix->No = ALLOC(int, Size+2)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }

/* Create DoDirect vectors for use in spFactor(). */
#if REAL
    if (Matrix->DoRealDirect == NULL)
    {   if (( Matrix->DoRealDirect = ALLOC(BOOLEAN, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
#endif
#if spCOMPLEX
    if (Matrix->DoCmplxDirect == NULL)
    {   if (( Matrix->DoCmplxDirect = ALLOC(BOOLEAN, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
#endif

/* Create Intermediate vectors for use in MatrixSolve. */
#if spCOMPLEX
    if (Matrix->Intermediate == NULL) {
      if ((Matrix->Intermediate = ALLOC(RealNumber,2*(Size+1))) == NULL)
            Matrix->Error = spNO_MEMORY;
      if ((Matrix->Intermediate2 = ALLOC(RealNumber,2*(Size+1))) == NULL)
            Matrix->Error = spNO_MEMORY;
      if ((Matrix->Intermediate3 = ALLOC(RealNumber,2*(Size+1))) == NULL)
            Matrix->Error = spNO_MEMORY;
      if ((Matrix->Intermediate4 = ALLOC(RealNumber,2*(Size+1))) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
#else
    if (Matrix->Intermediate == NULL) {
      if ((Matrix->Intermediate = ALLOC(RealNumber, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
      if ((Matrix->Intermediate2 = ALLOC(RealNumber, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
      if ((Matrix->Intermediate3 = ALLOC(RealNumber, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
      if ((Matrix->Intermediate4 = ALLOC(RealNumber, Size+1)) == NULL)
            Matrix->Error = spNO_MEMORY;
    }
#endif

    if (Matrix->Error != spNO_MEMORY)
        Matrix->InternalVectorsAllocated = YES;
    return;
}







/*
 *  COUNT MARKOWITZ
 *
 *  Scans Matrix to determine the Markowitz counts for each row and column.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  RHS  <input>  (RealVector)
 *      Representative right-hand side vector that is used to determine
 *      pivoting order when the right hand side vector is sparse.  If
 *      RHS is a NULL pointer then the RHS vector is assumed to be full
 *      and it is not used when determining the pivoting order.
 *  Step  <input>  (int)
 *     Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  Count  (int)
 *     Temporary counting variable.
 *  ExtRow  (int)
 *     The external row number that corresponds to I.
 *  pElement  (ElementPtr)
 *     Pointer to matrix elements.
 *  Size  (int)
 *     The size of the matrix.
 */

static void
CountMarkowitz( MatrixPtr Matrix, RealVector RHS, int Step )
{
int  Count, I, Size = Matrix->Size;
ElementPtr  pElement;
int  ExtRow;

/* Begin `CountMarkowitz'. */

/* Correct array pointer for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS OR NOT spCOMPLEX
        if (RHS != NULL) --RHS;
#else
        if (RHS != NULL)
        {   if (Matrix->Complex) RHS -= 2;
            else --RHS;
        }
#endif
#endif

/* Generate MarkowitzRow Count for each row. */
    for (I = Step; I <= Size; I++)
    {
/* Set Count to -1 initially to remove count due to pivot element. */
        Count = -1;
        pElement = Matrix->FirstInRow[I];
        while (pElement != NULL AND pElement->Col < Step)
            pElement = pElement->NextInRow;
        while (pElement != NULL)
        {   Count++;
            pElement = pElement->NextInRow;
        }

/* Include nonzero elements in the RHS vector. */
        ExtRow = Matrix->IntToExtRowMap[I];

#if spSEPARATED_COMPLEX_VECTORS OR NOT spCOMPLEX
        if (RHS != NULL)
            if (RHS[ExtRow] != 0.0)  Count++;
#else
        if (RHS != NULL)
        {   if (Matrix->Complex)
            {   if ((RHS[2*ExtRow] != 0.0) OR (RHS[2*ExtRow+1] != 0.0))
                    Count++;
            }
            else if (RHS[I] != 0.0) Count++;
        }
#endif
        Matrix->MarkowitzRow[I] = Count;
    }

/* Generate the MarkowitzCol count for each column. */
    for (I = Step; I <= Size; I++)
    {
/* Set Count to -1 initially to remove count due to pivot element. */
        Count = -1;
        pElement = Matrix->FirstInCol[I];
        while (pElement != NULL AND pElement->Row < Step)
            pElement = pElement->NextInCol;
        while (pElement != NULL)
        {   Count++;
            pElement = pElement->NextInCol;
        }
        Matrix->MarkowitzCol[I] = Count;
    }
    return;
}










/*
 *  MARKOWITZ PRODUCTS
 *
 *  Calculates MarkowitzProduct for each diagonal element from the Markowitz
 *  counts.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local Variables:
 *  pMarkowitzProduct  (long *)
 *      Pointer that points into MarkowitzProduct array. Is used to
 *      sequentially access entries quickly.
 *  pMarkowitzRow  (int *)
 *      Pointer that points into MarkowitzRow array. Is used to sequentially
 *      access entries quickly.
 *  pMarkowitzCol  (int *)
 *      Pointer that points into MarkowitzCol array. Is used to sequentially
 *      access entries quickly.
 *  Product  (long)
 *      Temporary storage for Markowitz product./
 *  Size  (int)
 *      The size of the matrix.
 */

static void
MarkowitzProducts( MatrixPtr Matrix, int Step )
{
int  I, *pMarkowitzRow, *pMarkowitzCol;
long  Product, *pMarkowitzProduct;
int  Size = Matrix->Size;
double fProduct;

/* Begin `MarkowitzProducts'. */
    Matrix->Singletons = 0;

    pMarkowitzProduct = &(Matrix->MarkowitzProd[Step]);
    pMarkowitzRow = &(Matrix->MarkowitzRow[Step]);
    pMarkowitzCol = &(Matrix->MarkowitzCol[Step]);

    for (I = Step; I <= Size; I++)
    {
/* If chance of overflow, use real numbers. */
        if ((*pMarkowitzRow > LARGEST_SHORT_INTEGER AND *pMarkowitzCol != 0) OR
            (*pMarkowitzCol > LARGEST_SHORT_INTEGER AND *pMarkowitzRow != 0))
        {   fProduct = (double)(*pMarkowitzRow++) * (double)(*pMarkowitzCol++);
            if ((long)fProduct >= LARGEST_LONG_INTEGER)
                *pMarkowitzProduct++ = LARGEST_LONG_INTEGER;
            else
                *pMarkowitzProduct++ = fProduct;
        }
        else
        {   Product = *pMarkowitzRow++ * *pMarkowitzCol++;
            if ((*pMarkowitzProduct++ = Product) == 0)
                Matrix->Singletons++;
        }
    }
    return;
}











/*
 *  SEARCH FOR BEST PIVOT
 *
 *  Performs a search to determine the element with the lowest Markowitz
 *  Product that is also acceptable.  An acceptable element is one that is
 *  larger than the AbsThreshold and at least as large as RelThreshold times
 *  the largest element in the same column.  The first step is to look for
 *  singletons if any exist.  If none are found, then all the diagonals are
 *  searched. The diagonal is searched once quickly using the assumption that
 *  elements on the diagonal are large compared to other elements in their
 *  column, and so the pivot can be chosen only on the basis of the Markowitz
 *  criterion.  After a element has been chosen to be pivot on the basis of
 *  its Markowitz product, it is checked to see if it is large enough.
 *  Waiting to the end of the Markowitz search to check the size of a pivot
 *  candidate saves considerable time, but is not guaranteed to find an
 *  acceptable pivot.  Thus if unsuccessful a second pass of the diagonal is
 *  made.  This second pass checks to see if an element is large enough during
 *  the search, not after it.  If still no acceptable pivot candidate has
 *  been found, the search expands to cover the entire matrix.
 *
 *  >>> Returned:
 *  A pointer to the element chosen to be pivot.  If every element in the
 *  matrix is zero, then NULL is returned.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      The row and column number of the beginning of the reduced submatrix.
 *
 *  >>> Local variables:
 *  ChosenPivot  (ElementPtr)
 *      Pointer to element that has been chosen to be the pivot.
 *
 *  >>> Possible errors:
 *  spSINGULAR
 *  spSMALL_PIVOT
 */

static ElementPtr
SearchForPivot( MatrixPtr Matrix, int Step, int DiagPivoting )
{
ElementPtr  ChosenPivot;

/* Begin `SearchForPivot'. */

/* If singletons exist, look for an acceptable one to use as pivot. */
    if (Matrix->Singletons)
    {   ChosenPivot = SearchForSingleton( Matrix, Step );
        if (ChosenPivot != NULL)
        {   Matrix->PivotSelectionMethod = 's';
            return ChosenPivot;
        }
    }

#if DIAGONAL_PIVOTING
    if (DiagPivoting)
    {
/*
 * Either no singletons exist or they weren't acceptable.  Take quick first
 * pass at searching diagonal.  First search for element on diagonal of 
 * remaining submatrix with smallest Markowitz product, then check to see
 * if it okay numerically.  If not, QuicklySearchDiagonal fails.
 */
/* This can lead to stability problems, so it is not used:
        ChosenPivot = QuicklySearchDiagonal( Matrix, Step );
        if (ChosenPivot != NULL)
        {   Matrix->PivotSelectionMethod = 'q';
            return ChosenPivot;
        }
*/

/*
 * Quick search of diagonal failed, carefully search diagonal and check each
 * pivot candidate numerically before even tentatively accepting it.
 */
        ChosenPivot = SearchDiagonal( Matrix, Step );
        if (ChosenPivot != NULL)
        {   Matrix->PivotSelectionMethod = 'd';
            return ChosenPivot;
        }
    }
#endif /* DIAGONAL_PIVOTING */

/* No acceptable pivot found yet, search entire matrix. */
    ChosenPivot = SearchEntireMatrix( Matrix, Step );
    Matrix->PivotSelectionMethod = 'e';
    if (ChosenPivot == NULL) {
      printf ("NULL returned from SearchEntireMatrix on step = %d\n",Step);
    }

    return ChosenPivot;
}









/*
 *  SEARCH FOR SINGLETON TO USE AS PIVOT
 *
 *  Performs a search to find a singleton to use as the pivot.  The
 *  first acceptable singleton is used.  A singleton is acceptable if
 *  it is larger in magnitude than the AbsThreshold and larger
 *  than RelThreshold times the largest of any other elements in the same
 *  column.  It may seem that a singleton need not satisfy the
 *  relative threshold criterion, however it is necessary to prevent
 *  excessive growth in the RHS from resulting in overflow during the
 *  forward and backward substitution.  A singleton does not need to
 *  be on the diagonal to be selected.
 *
 *  >>> Returned:
 *  A pointer to the singleton chosen to be pivot.  In no singleton is
 *  acceptable, return NULL.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  ChosenPivot  (ElementPtr)
 *      Pointer to element that has been chosen to be the pivot.
 *  PivotMag  (RealNumber)
 *      Magnitude of ChosenPivot.
 *  Singletons  (int)
 *      The count of the number of singletons that can be used as pivots.
 *      A local version of Matrix->Singletons.
 *  pMarkowitzProduct  (long *)
 *      Pointer that points into MarkowitzProduct array. It is used to quickly
 *      access successive Markowitz products.
 */

static ElementPtr
SearchForSingleton( MatrixPtr Matrix, int Step )
{
ElementPtr  ChosenPivot;
int  I;
long  *pMarkowitzProduct;
int  Singletons;
RealNumber  PivotMag;

/* Begin `SearchForSingleton'. */
/* Initialize pointer that is to scan through MarkowitzProduct vector. */
    pMarkowitzProduct = &(Matrix->MarkowitzProd[Matrix->Size+1]);
    Matrix->MarkowitzProd[Matrix->Size+1] = Matrix->MarkowitzProd[Step];

/* Decrement the count of available singletons, on the assumption that an
 * acceptable one will be found. */
    Singletons = Matrix->Singletons--;

/*
 * Assure that following while loop will always terminate, this is just
 * preventive medicine, if things are working right this should never
 * be needed.
 */
    Matrix->MarkowitzProd[Step-1] = 0;

    while (Singletons-- > 0)
    {
/* Singletons exist, find them. */

/*
 * This is tricky.  Am using a pointer to sequentially step through the
 * MarkowitzProduct array.  Search terminates when singleton (Product = 0)
 * is found.  Note that the conditional in the while statement
 * ( *pMarkowitzProduct ) is true as long as the MarkowitzProduct is not
 * equal to zero.  The row (and column) index on the diagonal is then
 * calculated by subtracting the pointer to the Markowitz product of
 * the first diagonal from the pointer to the Markowitz product of the
 * desired element, the singleton.
 *
 * Search proceeds from the end (high row and column numbers) to the
 * beginning (low row and column numbers) so that rows and columns with
 * large Markowitz products will tend to be move to the bottom of the
 * matrix.  However, choosing Diag[Step] is desirable because it would
 * require no row and column interchanges, so inspect it first by
 * putting its Markowitz product at the end of the MarkowitzProd
 * vector.
 */

        while ( *pMarkowitzProduct-- )
        {   /*
             * N bottles of beer on the wall;
             * N bottles of beer.
             * you take one down and pass it around;
             * N-1 bottles of beer on the wall.
             */
        }
        I = pMarkowitzProduct - Matrix->MarkowitzProd + 1;

/* Assure that I is valid. */
        if (I < Step) break;  /* while (Singletons-- > 0) */
        if (I > Matrix->Size) I = Step;

/* Singleton has been found in either/both row or/and column I. */
        if ((ChosenPivot = Matrix->Diag[I]) != NULL)
        {
/* Singleton lies on the diagonal. */
            PivotMag = ELEMENT_MAG(ChosenPivot);
            if
            (    PivotMag > Matrix->AbsThreshold AND
                 PivotMag > Matrix->RelThreshold *
                            FindBiggestInColExclude( Matrix, ChosenPivot, Step )
            ) return ChosenPivot;
        }
        else
        {
/* Singleton does not lie on diagonal, find it. */
            if (Matrix->MarkowitzCol[I] == 0)
            {   ChosenPivot = Matrix->FirstInCol[I];
                while ((ChosenPivot != NULL) AND (ChosenPivot->Row < Step))
                    ChosenPivot = ChosenPivot->NextInCol;
		if (ChosenPivot != NULL)
		{   /* Reduced column has no elements, matrix is singular. */
		    break;
		}
                PivotMag = ELEMENT_MAG(ChosenPivot);
                if
                (    PivotMag > Matrix->AbsThreshold AND
                     PivotMag > Matrix->RelThreshold *
                                FindBiggestInColExclude( Matrix, ChosenPivot,
                                                         Step )
                ) return ChosenPivot;
                else
                {   if (Matrix->MarkowitzRow[I] == 0)
                    {   ChosenPivot = Matrix->FirstInRow[I];
                        while((ChosenPivot != NULL) AND (ChosenPivot->Col<Step))
                            ChosenPivot = ChosenPivot->NextInRow;
			if (ChosenPivot != NULL)
			{ /* Reduced row has no elements, matrix is singular. */
			    break;
			}
                        PivotMag = ELEMENT_MAG(ChosenPivot);
                        if
                        (    PivotMag > Matrix->AbsThreshold AND
                             PivotMag > Matrix->RelThreshold *
                                        FindBiggestInColExclude( Matrix,
                                                                 ChosenPivot,
                                                                 Step )
                        ) return ChosenPivot;
                    }
                }
            }
            else
            {   ChosenPivot = Matrix->FirstInRow[I];
                while ((ChosenPivot != NULL) AND (ChosenPivot->Col < Step))
                    ChosenPivot = ChosenPivot->NextInRow;
		if (ChosenPivot != NULL)
		{   /* Reduced row has no elements, matrix is singular. */
		    break;
		}
                PivotMag = ELEMENT_MAG(ChosenPivot);
                if
                (    PivotMag > Matrix->AbsThreshold AND
                     PivotMag > Matrix->RelThreshold *
                                FindBiggestInColExclude( Matrix, ChosenPivot,
                                                         Step )
                ) return ChosenPivot;
            }
        }
/* Singleton not acceptable (too small), try another. */
    } /* end of while(lSingletons>0) */

/*
 * All singletons were unacceptable.  Restore Matrix->Singletons count.
 * Initial assumption that an acceptable singleton would be found was wrong.
 */
    Matrix->Singletons++;
    return NULL;
}












#if DIAGONAL_PIVOTING
#if MODIFIED_MARKOWITZ
/*
 *  QUICK SEARCH OF DIAGONAL FOR PIVOT WITH MODIFIED MARKOWITZ CRITERION
 *
 *  Searches the diagonal looking for the best pivot.  For a pivot to be
 *  acceptable it must be larger than the pivot RelThreshold times the largest
 *  element in its reduced column.  Among the acceptable diagonals, the
 *  one with the smallest MarkowitzProduct is sought.  Search terminates
 *  early if a diagonal is found with a MarkowitzProduct of one and its
 *  magnitude is larger than the other elements in its row and column.
 *  Since its MarkowitzProduct is one, there is only one other element in
 *  both its row and column, and, as a condition for early termination,
 *  these elements must be located symmetricly in the matrix.  If a tie
 *  occurs between elements of equal MarkowitzProduct, then the element with
 *  the largest ratio between its magnitude and the largest element in its
 *  column is used.  The search will be terminated after a given number of
 *  ties have occurred and the best (largest ratio) of the tied element will
 *  be used as the pivot.  The number of ties that will trigger an early
 *  termination is MinMarkowitzProduct * TIES_MULTIPLIER.
 *
 *  >>> Returned:
 *  A pointer to the diagonal element chosen to be pivot.  If no diagonal is
 *  acceptable, a NULL is returned.
 *
 *  >>> Arguments:
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  ChosenPivot  (ElementPtr)
 *      Pointer to the element that has been chosen to be the pivot.
 *  LargestOffDiagonal  (RealNumber)
 *      Magnitude of the largest of the off-diagonal terms associated with
 *      a diagonal with MarkowitzProduct equal to one.
 *  Magnitude  (RealNumber)
 *      Absolute value of diagonal element.
 *  MaxRatio  (RealNumber)
 *      Among the elements tied with the smallest Markowitz product, MaxRatio
 *      is the best (smallest) ratio of LargestInCol to the diagonal Magnitude
 *      found so far.  The smaller the ratio, the better numerically the
 *      element will be as pivot.
 *  MinMarkowitzProduct  (long)
 *      Smallest Markowitz product found of pivot candidates that lie along
 *      diagonal.
 *  NumberOfTies  (int)
 *      A count of the number of Markowitz ties that have occurred at current
 *      MarkowitzProduct.
 *  pDiag  (ElementPtr)
 *      Pointer to current diagonal element.
 *  pMarkowitzProduct  (long *)
 *      Pointer that points into MarkowitzProduct array. It is used to quickly
 *      access successive Markowitz products.
 *  Ratio  (RealNumber)
 *      For the current pivot candidate, Ratio is the ratio of the largest
 *      element in its column (excluding itself) to its magnitude.
 *  TiedElements  (ElementPtr[])
 *      Array of pointers to the elements with the minimum Markowitz
 *      product.
 *  pOtherInCol  (ElementPtr)
 *      When there is only one other element in a column other than the
 *      diagonal, pOtherInCol is used to point to it.  Used when Markowitz
 *      product is to determine if off diagonals are placed symmetricly.
 *  pOtherInRow  (ElementPtr)
 *      When there is only one other element in a row other than the diagonal,
 *      pOtherInRow is used to point to it.  Used when Markowitz product is
 *      to determine if off diagonals are placed symmetricly.
 */

static ElementPtr
QuicklySearchDiagonal( Matrix, Step )

MatrixPtr Matrix;
int Step;
{
long  MinMarkowitzProduct, *pMarkowitzProduct;
ElementPtr  pDiag, pOtherInRow, pOtherInCol;
int  I, NumberOfTies;
ElementPtr  ChosenPivot, TiedElements[MAX_MARKOWITZ_TIES + 1];
RealNumber  Magnitude, LargestInCol, Ratio, MaxRatio;
RealNumber  LargestOffDiagonal;

/* Begin `QuicklySearchDiagonal'. */
    NumberOfTies = -1;
    MinMarkowitzProduct = LARGEST_LONG_INTEGER;
    pMarkowitzProduct = &(Matrix->MarkowitzProd[Matrix->Size+2]);
    Matrix->MarkowitzProd[Matrix->Size+1] = Matrix->MarkowitzProd[Step];

/* Assure that following while loop will always terminate. */
    Matrix->MarkowitzProd[Step-1] = -1;

/*
 * This is tricky.  Am using a pointer in the inner while loop to
 * sequentially step through the MarkowitzProduct array.  Search
 * terminates when the Markowitz product of zero placed at location
 * Step-1 is found.  The row (and column) index on the diagonal is then
 * calculated by subtracting the pointer to the Markowitz product of
 * the first diagonal from the pointer to the Markowitz product of the
 * desired element. The outer for loop is infinite, broken by using
 * break.
 *
 * Search proceeds from the end (high row and column numbers) to the
 * beginning (low row and column numbers) so that rows and columns with
 * large Markowitz products will tend to be move to the bottom of the
 * matrix.  However, choosing Diag[Step] is desirable because it would
 * require no row and column interchanges, so inspect it first by
 * putting its Markowitz product at the end of the MarkowitzProd
 * vector.
 */

    for(;;)  /* Endless for loop. */
    {   while (MinMarkowitzProduct < *(--pMarkowitzProduct))
        {   /*
             * N bottles of beer on the wall;
             * N bottles of beer.
             * You take one down and pass it around;
             * N-1 bottles of beer on the wall.
             */
        }

        I = pMarkowitzProduct - Matrix->MarkowitzProd;

/* Assure that I is valid; if I < Step, terminate search. */
        if (I < Step) break; /* Endless for loop */
        if (I > Matrix->Size) I = Step;

        if ((pDiag = Matrix->Diag[I]) == NULL)
            continue; /* Endless for loop */
        if ((Magnitude = ELEMENT_MAG(pDiag)) <= Matrix->AbsThreshold)
            continue; /* Endless for loop */

        if (*pMarkowitzProduct == 1)
        {
/* Case where only one element exists in row and column other than diagonal. */

/* Find off diagonal elements. */
            pOtherInRow = pDiag->NextInRow;
            pOtherInCol = pDiag->NextInCol;
            if (pOtherInRow == NULL AND pOtherInCol == NULL)
            {    pOtherInRow = Matrix->FirstInRow[I];
                 while(pOtherInRow != NULL)
                 {   if (pOtherInRow->Col >= Step AND pOtherInRow->Col != I)
                         break;
                     pOtherInRow = pOtherInRow->NextInRow;
                 }
                 pOtherInCol = Matrix->FirstInCol[I];
                 while(pOtherInCol != NULL)
                 {   if (pOtherInCol->Row >= Step AND pOtherInCol->Row != I)
                         break;
                     pOtherInCol = pOtherInCol->NextInCol;
                 }
            }

/* Accept diagonal as pivot if diagonal is larger than off diagonals and the
 * off diagonals are placed symmetricly. */
            if (pOtherInRow != NULL  AND  pOtherInCol != NULL)
            {   if (pOtherInRow->Col == pOtherInCol->Row)
                {   LargestOffDiagonal = MAX(ELEMENT_MAG(pOtherInRow),
                                                      ELEMENT_MAG(pOtherInCol));
                    if (Magnitude >= LargestOffDiagonal)
                    {
/* Accept pivot, it is unlikely to contribute excess error. */
                        return pDiag;
                    }
                }
            }
        }

        if (*pMarkowitzProduct < MinMarkowitzProduct)
        {
/* Notice strict inequality in test. This is a new smallest MarkowitzProduct. */
            TiedElements[0] = pDiag;
            MinMarkowitzProduct = *pMarkowitzProduct;
            NumberOfTies = 0;
        }
        else
        {
/* This case handles Markowitz ties. */
            if (NumberOfTies < MAX_MARKOWITZ_TIES)
            {   TiedElements[++NumberOfTies] = pDiag;
                if (NumberOfTies >= MinMarkowitzProduct * TIES_MULTIPLIER)
                    break; /* Endless for loop */
            }
        }
    } /* End of endless for loop. */

/* Test to see if any element was chosen as a pivot candidate. */
    if (NumberOfTies < 0)
        return NULL;

/* Determine which of tied elements is best numerically. */
    ChosenPivot = NULL;
    MaxRatio = 1.0 / Matrix->RelThreshold;

    for (I = 0; I <= NumberOfTies; I++)
    {   pDiag = TiedElements[I];
        Magnitude = ELEMENT_MAG(pDiag);
        LargestInCol = FindBiggestInColExclude( Matrix, pDiag, Step );
        Ratio = LargestInCol / Magnitude;
        if (Ratio < MaxRatio)
        {   ChosenPivot = pDiag;
            MaxRatio = Ratio;
        }
    }
    return ChosenPivot;
}










#else /* Not MODIFIED_MARKOWITZ */
/*
 *  QUICK SEARCH OF DIAGONAL FOR PIVOT WITH CONVENTIONAL MARKOWITZ
 *  CRITERION
 *
 *  Searches the diagonal looking for the best pivot.  For a pivot to be
 *  acceptable it must be larger than the pivot RelThreshold times the largest
 *  element in its reduced column.  Among the acceptable diagonals, the
 *  one with the smallest MarkowitzProduct is sought.  Search terminates
 *  early if a diagonal is found with a MarkowitzProduct of one and its
 *  magnitude is larger than the other elements in its row and column.
 *  Since its MarkowitzProduct is one, there is only one other element in
 *  both its row and column, and, as a condition for early termination,
 *  these elements must be located symmetricly in the matrix.
 *
 *  >>> Returned:
 *  A pointer to the diagonal element chosen to be pivot.  If no diagonal is
 *  acceptable, a NULL is returned.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  ChosenPivot  (ElementPtr)
 *      Pointer to the element that has been chosen to be the pivot.
 *  LargestOffDiagonal  (RealNumber)
 *      Magnitude of the largest of the off-diagonal terms associated with
 *      a diagonal with MarkowitzProduct equal to one.
 *  Magnitude  (RealNumber)
 *      Absolute value of diagonal element.
 *  MinMarkowitzProduct  (long)
 *      Smallest Markowitz product found of pivot candidates which are
 *      acceptable.
 *  pDiag  (ElementPtr)
 *      Pointer to current diagonal element.
 *  pMarkowitzProduct  (long *)
 *      Pointer that points into MarkowitzProduct array. It is used to quickly
 *      access successive Markowitz products.
 *  pOtherInCol  (ElementPtr)
 *      When there is only one other element in a column other than the
 *      diagonal, pOtherInCol is used to point to it.  Used when Markowitz
 *      product is to determine if off diagonals are placed symmetricly.
 *  pOtherInRow  (ElementPtr)
 *      When there is only one other element in a row other than the diagonal,
 *      pOtherInRow is used to point to it.  Used when Markowitz product is
 *      to determine if off diagonals are placed symmetricly.
 */

static ElementPtr
QuicklySearchDiagonal( MatrixPtr Matrix, int Step )
{
long  MinMarkowitzProduct, *pMarkowitzProduct;
ElementPtr  pDiag;
int  I;
ElementPtr  ChosenPivot, pOtherInRow, pOtherInCol;
RealNumber  Magnitude, LargestInCol, LargestOffDiagonal;

/* Begin `QuicklySearchDiagonal'. */
    ChosenPivot = NULL;
    MinMarkowitzProduct = LARGEST_LONG_INTEGER;
    pMarkowitzProduct = &(Matrix->MarkowitzProd[Matrix->Size+2]);
    Matrix->MarkowitzProd[Matrix->Size+1] = Matrix->MarkowitzProd[Step];

/* Assure that following while loop will always terminate. */
    Matrix->MarkowitzProd[Step-1] = -1;

/*
 * This is tricky.  Am using a pointer in the inner while loop to
 * sequentially step through the MarkowitzProduct array.  Search
 * terminates when the Markowitz product of zero placed at location
 * Step-1 is found.  The row (and column) index on the diagonal is then
 * calculated by subtracting the pointer to the Markowitz product of
 * the first diagonal from the pointer to the Markowitz product of the
 * desired element. The outer for loop is infinite, broken by using
 * break.
 *
 * Search proceeds from the end (high row and column numbers) to the
 * beginning (low row and column numbers) so that rows and columns with
 * large Markowitz products will tend to be move to the bottom of the
 * matrix.  However, choosing Diag[Step] is desirable because it would
 * require no row and column interchanges, so inspect it first by
 * putting its Markowitz product at the end of the MarkowitzProd
 * vector.
 */

    for (;;)  /* Endless for loop. */
    {   while (*(--pMarkowitzProduct) >= MinMarkowitzProduct)
        {   /* Just passing through. */
        }

        I = pMarkowitzProduct - Matrix->MarkowitzProd;

/* Assure that I is valid; if I < Step, terminate search. */
        if (I < Step) break; /* Endless for loop */
        if (I > Matrix->Size) I = Step;

        if ((pDiag = Matrix->Diag[I]) == NULL)
            continue; /* Endless for loop */
        if ((Magnitude = ELEMENT_MAG(pDiag)) <= Matrix->AbsThreshold)
            continue; /* Endless for loop */

        if (*pMarkowitzProduct == 1)
        {
/* Case where only one element exists in row and column other than diagonal. */

/* Find off-diagonal elements. */
            pOtherInRow = pDiag->NextInRow;
            pOtherInCol = pDiag->NextInCol;
            if (pOtherInRow == NULL AND pOtherInCol == NULL)
            {    pOtherInRow = Matrix->FirstInRow[I];
                 while(pOtherInRow != NULL)
                 {   if (pOtherInRow->Col >= Step AND pOtherInRow->Col != I)
                         break;
                     pOtherInRow = pOtherInRow->NextInRow;
                 }
                 pOtherInCol = Matrix->FirstInCol[I];
                 while(pOtherInCol != NULL)
                 {   if (pOtherInCol->Row >= Step AND pOtherInCol->Row != I)
                         break;
                     pOtherInCol = pOtherInCol->NextInCol;
                 }
            }

/* Accept diagonal as pivot if diagonal is larger than off-diagonals and the
 * off-diagonals are placed symmetricly. */
            if (pOtherInRow != NULL  AND  pOtherInCol != NULL)
            {   if (pOtherInRow->Col == pOtherInCol->Row)
                {   LargestOffDiagonal = MAX(ELEMENT_MAG(pOtherInRow),
                                                      ELEMENT_MAG(pOtherInCol));
                    if (Magnitude >= LargestOffDiagonal)
                    {
/* Accept pivot, it is unlikely to contribute excess error. */
                        return pDiag;
                    }
                }
            }
        }

        MinMarkowitzProduct = *pMarkowitzProduct;
        ChosenPivot = pDiag;
    }  /* End of endless for loop. */

    if (ChosenPivot != NULL)
    {   LargestInCol = FindBiggestInColExclude( Matrix, ChosenPivot, Step );
        if( ELEMENT_MAG(ChosenPivot) <= Matrix->RelThreshold * LargestInCol )
            ChosenPivot = NULL;
    }
    return ChosenPivot;
}
#endif /* Not MODIFIED_MARKOWITZ */









/*
 *  SEARCH DIAGONAL FOR PIVOT WITH MODIFIED MARKOWITZ CRITERION
 *
 *  Searches the diagonal looking for the best pivot.  For a pivot to be
 *  acceptable it must be larger than the pivot RelThreshold times the largest
 *  element in its reduced column.  Among the acceptable diagonals, the
 *  one with the smallest MarkowitzProduct is sought.  If a tie occurs
 *  between elements of equal MarkowitzProduct, then the element with
 *  the largest ratio between its magnitude and the largest element in its
 *  column is used.  The search will be terminated after a given number of
 *  ties have occurred and the best (smallest ratio) of the tied element will
 *  be used as the pivot.  The number of ties that will trigger an early
 *  termination is MinMarkowitzProduct * TIES_MULTIPLIER.
 *
 *  >>> Returned:
 *  A pointer to the diagonal element chosen to be pivot.  If no diagonal is
 *  acceptable, a NULL is returned.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  ChosenPivot  (ElementPtr)
 *      Pointer to the element that has been chosen to be the pivot.
 *  Size  (int)
 *      Local version of size which is placed in a register to increase speed.
 *  Magnitude  (RealNumber)
 *      Absolute value of diagonal element.
 *  MinMarkowitzProduct  (long)
 *      Smallest Markowitz product found of those pivot candidates which are
 *      acceptable.
 *  NumberOfTies  (int)
 *      A count of the number of Markowitz ties that have occurred at current
 *      MarkowitzProduct.
 *  pDiag  (ElementPtr)
 *      Pointer to current diagonal element.
 *  pMarkowitzProduct  (long*)
 *      Pointer that points into MarkowitzProduct array. It is used to quickly
 *      access successive Markowitz products.
 *  Ratio  (RealNumber)
 *      For the current pivot candidate, Ratio is the
 *      Ratio of the largest element in its column to its magnitude.
 *  RatioOfAccepted  (RealNumber)
 *      For the best pivot candidate found so far, RatioOfAccepted is the
 *      Ratio of the largest element in its column to its magnitude.
 */
ElementPtr
SearchDiagonal( MatrixPtr Matrix, int Step )
{
int  J;
long  *pMarkowitzProduct, MarkOfAccepted;
int  I;
ElementPtr  pDiag;
int  NumberOfTies, Size = Matrix->Size;
ElementPtr  ChosenPivot;
RealNumber  Magnitude, Ratio, RatioOfAccepted, LargestInCol;

int start,step,iteration;
static long MinMarkowitzProduct, SecondMinMarkowitzProduct;

/* Begin `SearchDiagonal'. */
    start = Size+1;
    step = 1;

    ChosenPivot = NULL;
    iteration = 0;
    while (iteration++ < 2 && ChosenPivot == NULL) {
      RatioOfAccepted = 2/Matrix->RelThreshold;
      MarkOfAccepted = MinMarkowitzProduct;
      Matrix->MarkowitzProd[start] = Matrix->MarkowitzProd[Step];
      pMarkowitzProduct = &(Matrix->MarkowitzProd[start+step]);
/* Start search of diagonal. */
      for (J = start; J > Step; J-=step)
      {
          pMarkowitzProduct -= step;
/* qqq this is other part of selection criterion */
          if (*pMarkowitzProduct > MinMarkowitzProduct)
            continue;
          if (J > Matrix->Size) {
            I = Step;
          }
          else {
            I = J;
          }
          if ((pDiag = Matrix->Diag[I]) == NULL)
              continue; /* for loop */
          if ((Magnitude = ELEMENT_MAG(pDiag)) <= Matrix->AbsThreshold)
              continue; /* for loop */

/* Test to see if diagonal's magnitude is acceptable. */
/*
          LargestInCol = FindBiggestInColExclude( Matrix, pDiag, Step );
*/
          LargestInCol = Matrix->Intermediate4[I];
          if (LargestInCol == 0) {
            LargestInCol = FindBiggestInColExclude( Matrix, pDiag, Step );
            Matrix->Intermediate4[I] = LargestInCol;
          }
          if (Magnitude <= Matrix->RelThreshold * LargestInCol)
              continue; /* for loop */

          if (*pMarkowitzProduct < MinMarkowitzProduct)
          {
/* Notice strict inequality in test. This is a new smallest MarkowitzProduct. */
              if (MinMarkowitzProduct == LARGEST_LONG_INTEGER) {
                ChosenPivot = pDiag;
                RatioOfAccepted = LargestInCol / Magnitude;
                MarkOfAccepted = *pMarkowitzProduct;
              }
              NumberOfTies = 0;
              SecondMinMarkowitzProduct = MinMarkowitzProduct;
              MinMarkowitzProduct = *pMarkowitzProduct;
          }
/* This case handles Markowitz ties. */
          NumberOfTies++;
          Ratio = LargestInCol / Magnitude;
/* qqq This condition is critical to stability!!
   best for ken's problem:
          if (Ratio < RatioOfAccepted || *pMarkowitzProduct*2 <= MarkOfAccepted) { }
   best for test suite (sandler4):
          if (Ratio < RatioOfAccepted || *pMarkowitzProduct < MarkOfAccepted) { }
   gives consistent parallel and serial runs
          if (*pMarkowitzProduct < MarkOfAccepted ||
              *pMarkowitzProduct == MarkOfAccepted && Ratio < RatioOfAccepted) { }
*/
          if (*pMarkowitzProduct < MarkOfAccepted ||
              *pMarkowitzProduct == MarkOfAccepted && Ratio < RatioOfAccepted) {
              ChosenPivot = pDiag;
              RatioOfAccepted = Ratio;
              MarkOfAccepted = *pMarkowitzProduct;
          }
      } /* End of for(Step) */

/*
        if (ChosenPivot)
          printf ("Final pivot at Step: %d = %d, mark = %ld, ratio = %g\n",Step, ChosenPivot->Col,
                   MinMarkowitzProduct, RatioOfAccepted);
*/

      if (ChosenPivot == NULL) {
        MinMarkowitzProduct = LARGEST_LONG_INTEGER;
      }
    }
    return ChosenPivot;
}
#endif /* DIAGONAL_PIVOTING */










/*
 *  SEARCH ENTIRE MATRIX FOR BEST PIVOT
 *
 *  Performs a search over the entire matrix looking for the acceptable
 *  element with the lowest MarkowitzProduct.  If there are several that
 *  are tied for the smallest MarkowitzProduct, the tie is broken by using
 *  the ratio of the magnitude of the element being considered to the largest
 *  element in the same column.  If no element is acceptable then the largest
 *  element in the reduced submatrix is used as the pivot and the
 *  matrix is declared to be spSMALL_PIVOT.  If the largest element is
 *  zero, the matrix is declared to be spSINGULAR.
 *
 *  >>> Returned:
 *  A pointer to the diagonal element chosen to be pivot.  If no element is
 *  found, then NULL is returned and the matrix is spSINGULAR.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  ChosenPivot  (ElementPtr)
 *      Pointer to the element that has been chosen to be the pivot.
 *  LargestElementMag  (RealNumber)
 *      Magnitude of the largest element yet found in the reduced submatrix.
 *  Size  (int)
 *      Local version of Size; placed in a register for speed.
 *  Magnitude  (RealNumber)
 *      Absolute value of diagonal element.
 *  MinMarkowitzProduct  (long)
 *      Smallest Markowitz product found of pivot candidates which are
 *      acceptable.
 *  NumberOfTies  (int)
 *      A count of the number of Markowitz ties that have occurred at current
 *      MarkowitzProduct.
 *  pElement  (ElementPtr)
 *      Pointer to current element.
 *  pLargestElement  (ElementPtr)
 *      Pointer to the largest element yet found in the reduced submatrix.
 *  Product  (long)
 *      Markowitz product for the current row and column.
 *  Ratio  (RealNumber)
 *      For the current pivot candidate, Ratio is the
 *      Ratio of the largest element in its column to its magnitude.
 *  RatioOfAccepted  (RealNumber)
 *      For the best pivot candidate found so far, RatioOfAccepted is the
 *      Ratio of the largest element in its column to its magnitude.
 *
 *  >>> Possible errors:
 *  spSINGULAR
 *  spSMALL_PIVOT
 */

ElementPtr
SearchEntireMatrix( MatrixPtr Matrix, int Step )
{
int  I, Size = Matrix->Size;
ElementPtr  pElement;
int  NumberOfTies;
long  Product, MinMarkowitzProduct, MarkOfAccepted;
ElementPtr  ChosenPivot, pLargestElement;
RealNumber  Magnitude, LargestElementMag, Ratio, RatioOfAccepted, LargestInCol;

/* Begin `SearchEntireMatrix'. */
    ChosenPivot = NULL;
    LargestElementMag = 0.0;
    MinMarkowitzProduct = LARGEST_LONG_INTEGER;

/* Start search of matrix on column by column basis. */
    for (I = Step; I <= Size; I++)
    {   pElement = Matrix->FirstInCol[I];

        while (pElement != NULL AND pElement->Row < Step)
            pElement = pElement->NextInCol;

        if((LargestInCol = FindLargestInCol(pElement)) == 0.0)
            continue; /* for loop */

        while (pElement != NULL)
        {
/* Check to see if element is the largest encountered so far.  If so, record
   its magnitude and address. */
            if ((Magnitude = ELEMENT_MAG(pElement)) > LargestElementMag)
            {   LargestElementMag = Magnitude;
                pLargestElement = pElement;
            }
/* Calculate element's MarkowitzProduct. */
            Product = Matrix->MarkowitzRow[pElement->Row] *
                      Matrix->MarkowitzCol[pElement->Col];

/* Test to see if element is acceptable as a pivot candidate. */
            if ((Product <= MinMarkowitzProduct) AND
                (Magnitude > Matrix->RelThreshold * LargestInCol) AND
                (Magnitude > Matrix->AbsThreshold))
            {
/* Test to see if element has lowest MarkowitzProduct yet found, or whether it
   is tied with an element found earlier. */
                if (Product < MinMarkowitzProduct)
                {
/* Notice strict inequality in test. This is a new smallest MarkowitzProduct. */
                    if (MinMarkowitzProduct == LARGEST_LONG_INTEGER) {
                      ChosenPivot = pElement;
                      RatioOfAccepted = LargestInCol / Magnitude;
                      MarkOfAccepted = Product;
                    }
                    NumberOfTies = 0;
                    MinMarkowitzProduct = Product;
                }
/* This case handles Markowitz ties. */
                NumberOfTies++;
                Ratio = LargestInCol / Magnitude;
                if (Ratio < RatioOfAccepted || Product < MarkOfAccepted)
                {   ChosenPivot = pElement;
                    RatioOfAccepted = Ratio;
                    MarkOfAccepted = Product;
                }
                if (NumberOfTies >= MinMarkowitzProduct * TIES_MULTIPLIER)
                    return ChosenPivot;
            }
            pElement = pElement->NextInCol;
        }  /* End of while(pElement != NULL) */
    } /* End of for(Step) */

    if (ChosenPivot != NULL) return ChosenPivot;

    if (LargestElementMag == 0.0)
    {   Matrix->Error = spSINGULAR;
        return NULL;
    }

    Matrix->Error = spSMALL_PIVOT;
    return pLargestElement;
}











/*
 *  DETERMINE THE MAGNITUDE OF THE LARGEST ELEMENT IN A COLUMN
 *
 *  This routine searches a column and returns the magnitude of the largest
 *  element.  This routine begins the search at the element pointed to by
 *  pElement, the parameter.
 *
 *  The search is conducted by starting at the element specified by a pointer,
 *  which should be one below the diagonal, and moving down the column.  On
 *  the way down the column, the magnitudes of the elements are tested to see
 *  if they are the largest yet found.
 *
 *  >>> Returned:
 *  The magnitude of the largest element in the column below and including
 *  the one pointed to by the input parameter.
 *
 *  >>> Arguments:
 *  pElement  <input>  (ElementPtr)
 *      The pointer to the first element to be tested.  Also, used by the
 *      routine to access all lower elements in the column.
 *
 *  >>> Local variables:
 *  Largest  (RealNumber)
 *      The magnitude of the largest element.
 *  Magnitude  (RealNumber)
 *      The magnitude of the currently active element.
 */

static RealNumber
FindLargestInCol( ElementPtr pElement )
{
RealNumber  Magnitude, Largest = 0.0;

/* Begin `FindLargestInCol'. */
/* Search column for largest element beginning at Element. */
    while (pElement != NULL)
    {   if ((Magnitude = ELEMENT_MAG(pElement)) > Largest)
            Largest = Magnitude;
        pElement = pElement->NextInCol;
    }

    return Largest;
}










/*
 *  DETERMINE THE MAGNITUDE OF THE LARGEST ELEMENT IN A COLUMN
 *  EXCLUDING AN ELEMENT
 *
 *  This routine searches a column and returns the magnitude of the largest
 *  element.  One given element is specifically excluded from the search.
 *
 *  The search is conducted by starting at the first element in the column
 *  and moving down the column until the active part of the matrix is entered,
 *  i.e. the reduced submatrix.  The rest of the column is then traversed
 *  looking for the largest element.
 *
 *  >>> Returned:
 *  The magnitude of the largest element in the active portion of the column,
 *  excluding the specified element, is returned.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  pElement  <input>  (ElementPtr)
 *      The pointer to the element that is to be excluded from search. Column
 *      to be searched is one that contains this element.  Also used to
 *      access the elements in the column.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.  Indicates where
 *      the active part of the matrix begins.
 *
 *  >>> Local variables:
 *  Col  (int)
 *      The number of the column to be searched.  Also the column number of
 *      the element to be avoided in the search.
 *  Largest  (RealNumber)
 *      The magnitude of the largest element.
 *  Magnitude  (RealNumber)
 *      The magnitude of the currently active element.
 *  Row  (int)
 *      The row number of element to be excluded from the search.
 */

static RealNumber
FindBiggestInColExclude( MatrixPtr Matrix, ElementPtr pElement, int Step )
{
int  Row;
int  Col;
RealNumber  Largest, Magnitude;

/* Begin `FindBiggestInColExclude'. */
    Row = pElement->Row;
    Col = pElement->Col;

    pElement = Matrix->Col_fast[Col][f_ind(Matrix, Col,Step)];
    if (pElement == NULL || pElement->Col != Col || pElement->Row > Step) {
        pElement = Matrix->FirstInCol[Col];
    }
/* Travel down column until reduced submatrix is entered. */
    while ((pElement != NULL) AND (pElement->Row < Step))
        pElement = pElement->NextInCol;

/* Initialize the variable Largest. */
    if (pElement->Row != Row)
        Largest = ELEMENT_MAG(pElement);
    else
        Largest = 0.0;

/* Search rest of column for largest element, avoiding excluded element. */
    while ((pElement = pElement->NextInCol) != NULL)
    {   if ((Magnitude = ELEMENT_MAG(pElement)) > Largest)
        {   if (pElement->Row != Row)
                Largest = Magnitude;
        }
    }

    return Largest;
}










/*
 *  EXCHANGE ROWS AND COLUMNS
 *
 *  Exchanges two rows and two columns so that the selected pivot is moved to
 *  the upper left corner of the remaining submatrix.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  pPivot  <input>  (ElementPtr)
 *      Pointer to the current pivot.
 *  Step  <input>  (int)
 *      Index of the diagonal currently being eliminated.
 *
 *  >>> Local variables:
 *  Col  (int)
 *      Column where the pivot was found.
 *  Row  (int)
 *      Row where the pivot was found.
 *  OldMarkowitzProd_Col  (long)
 *      Markowitz product associated with the diagonal element in the row
 *      the pivot was found in.
 *  OldMarkowitzProd_Row  (long)
 *      Markowitz product associated with the diagonal element in the column
 *      the pivot was found in.
 *  OldMarkowitzProd_Step  (long)
 *      Markowitz product associated with the diagonal element that is being
 *      moved so that the pivot can be placed in the upper left-hand corner
 *      of the reduced submatrix.
 */

static void
ExchangeRowsAndCols( MatrixPtr Matrix, ElementPtr pPivot, int Step )
{
int   Row, Col;
long  OldMarkowitzProd_Step, OldMarkowitzProd_Row, OldMarkowitzProd_Col;

/* Begin `ExchangeRowsAndCols'. */
#ifdef CHECK_IND
    spCheckInd(Matrix,"Top ExchangeRowsAndCols");
#endif
    Row = pPivot->Row;
    Col = pPivot->Col;
    Matrix->PivotsOriginalRow = Row;
    Matrix->PivotsOriginalCol = Col;

    if ((Row == Step) AND (Col == Step)) return;

/* Exchange rows and columns. */
    if (Row == Col) {
        spcRowExchange( Matrix, Step, Row );
        spcColExchange( Matrix, Step, Col );
        SWAP( long, Matrix->MarkowitzProd[Step], Matrix->MarkowitzProd[Row] );
        SWAP( ElementPtr, Matrix->Diag[Row], Matrix->Diag[Step] );
    }
    else
    {

/* Initialize variables that hold old Markowitz products. */
        OldMarkowitzProd_Step = Matrix->MarkowitzProd[Step];
        OldMarkowitzProd_Row = Matrix->MarkowitzProd[Row];
        OldMarkowitzProd_Col = Matrix->MarkowitzProd[Col];

/* Exchange rows. */
        if (Row != Step) {
            spcRowExchange( Matrix, Step, Row );
            Matrix->NumberOfInterchangesIsOdd =
                                       NOT Matrix->NumberOfInterchangesIsOdd;
            Matrix->MarkowitzProd[Row] = Matrix->MarkowitzRow[Row] *
                                                   Matrix->MarkowitzCol[Row];

/* Update singleton count. */
            if ((Matrix->MarkowitzProd[Row]==0) != (OldMarkowitzProd_Row==0))
            {   if (OldMarkowitzProd_Row == 0)
                    Matrix->Singletons--;
                else
                    Matrix->Singletons++;
            }
        }

/* Exchange columns. */
        if (Col != Step) {
            spcColExchange( Matrix, Step, Col );
            Matrix->NumberOfInterchangesIsOdd =
                                       NOT Matrix->NumberOfInterchangesIsOdd;
            Matrix->MarkowitzProd[Col] = Matrix->MarkowitzCol[Col] *
                                                   Matrix->MarkowitzRow[Col];

/* Update singleton count. */
            if ((Matrix->MarkowitzProd[Col]==0) != (OldMarkowitzProd_Col==0))
            {   if (OldMarkowitzProd_Col == 0)
                    Matrix->Singletons--;
                else
                    Matrix->Singletons++;
            }

            Matrix->Diag[Col] = spcFindElementInCol( Matrix,
                                                     Matrix->FirstInCol+Col,
                                                     Col, Col, NO );
        }
        if (Row != Step)
        {   Matrix->Diag[Row] = spcFindElementInCol( Matrix,
                                                     Matrix->FirstInCol+Row,
                                                     Row, Row, NO );
        }
        Matrix->Diag[Step] = spcFindElementInCol( Matrix,
                                                  Matrix->FirstInCol+Step,
                                                  Step, Step, NO );

/* Update singleton count. */
        Matrix->MarkowitzProd[Step] = Matrix->MarkowitzCol[Step] *
                                                    Matrix->MarkowitzRow[Step];
        if ((Matrix->MarkowitzProd[Step]==0) != (OldMarkowitzProd_Step==0))
        {   if (OldMarkowitzProd_Step == 0)
                Matrix->Singletons--;
            else
                Matrix->Singletons++;
        }
    }
#ifdef CHECK_IND
    spCheckInd(Matrix,"Bottom ExchangeRowsAndCols");
#endif
    return;
}









/*
 *  EXCHANGE ROWS
 *
 *  Performs all required operations to exchange two rows. Those operations
 *  include: swap FirstInRow pointers, fixing up the NextInCol pointers,
 *  swapping row indexes in MatrixElements, and swapping Markowitz row
 *  counts.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  Row1  <input>  (int)
 *      Row index of one of the rows, becomes the smallest index.
 *  Row2  <input>  (int)
 *      Row index of the other row, becomes the largest index.
 *
 *  Local variables:
 *  Column  (int)
 *      Column in which row elements are currently being exchanged.
 *  Row1Ptr  (ElementPtr)
 *      Pointer to an element in Row1.
 *  Row2Ptr  (ElementPtr)
 *      Pointer to an element in Row2.
 *  Element1  (ElementPtr)
 *      Pointer to the element in Row1 to be exchanged.
 *  Element2  (ElementPtr)
 *      Pointer to the element in Row2 to be exchanged.
 */

/* The following would disable parallelization of row/col swapping

*/

void
spcRowExchange( MatrixPtr Matrix, int Row1, int Row2 )
{
ElementPtr  Row1Ptr, Row2Ptr;
int  Column;
ElementPtr  Element1, Element2;


/* Begin `spcRowExchange'. */
    if (Row1 > Row2)  SWAP(int, Row1, Row2);

    Row1Ptr = Matrix->FirstInRow[Row1];
    Row2Ptr = Matrix->FirstInRow[Row2];
    while (Row1Ptr != NULL OR Row2Ptr != NULL)
    {
/* Exchange elements in rows while traveling from left to right. */
        if (Row1Ptr == NULL)
        {   Column = Row2Ptr->Col;
            Element1 = NULL;
            Element2 = Row2Ptr;
            Row2Ptr = Row2Ptr->NextInRow;
        }
        else if (Row2Ptr == NULL)
        {   Column = Row1Ptr->Col;
            Element1 = Row1Ptr;
            Element2 = NULL;
            Row1Ptr = Row1Ptr->NextInRow;
        }
        else if (Row1Ptr->Col < Row2Ptr->Col)
        {   Column = Row1Ptr->Col;
            Element1 = Row1Ptr;
            Element2 = NULL;
            Row1Ptr = Row1Ptr->NextInRow;
        }
        else if (Row1Ptr->Col > Row2Ptr->Col)
        {   Column = Row2Ptr->Col;
            Element1 = NULL;
            Element2 = Row2Ptr;
            Row2Ptr = Row2Ptr->NextInRow;
        }
        else   /* Row1Ptr->Col == Row2Ptr->Col */
        {   Column = Row1Ptr->Col;
            Element1 = Row1Ptr;
            Element2 = Row2Ptr;
            Row1Ptr = Row1Ptr->NextInRow;
            Row2Ptr = Row2Ptr->NextInRow;
        }

          ExchangeColElements( Matrix, Row1, Element1, Row2, Element2, Column);

    }  /* end of while(Row1Ptr != NULL OR Row2Ptr != NULL) */

    if (Matrix->InternalVectorsAllocated)
        SWAP( int, Matrix->MarkowitzRow[Row1], Matrix->MarkowitzRow[Row2]);
    SWAP( ElementPtr, Matrix->FirstInRow[Row1], Matrix->FirstInRow[Row2]);
    SWAP( int, Matrix->IntToExtRowMap[Row1], Matrix->IntToExtRowMap[Row2]);
    SWAP( ElementPtr *, Matrix->Row_fast[Row1], Matrix->Row_fast[Row2]);
#if TRANSLATE
    Matrix->ExtToIntRowMap[ Matrix->IntToExtRowMap[Row1] ] = Row1;
    Matrix->ExtToIntRowMap[ Matrix->IntToExtRowMap[Row2] ] = Row2;
#endif

    spRowInd(Matrix, Row1);
    spRowInd(Matrix, Row2);
    return;
}





/*
 *  EXCHANGE COLUMNS
 *
 *  Performs all required operations to exchange two columns. Those operations
 *  include: swap FirstInCol pointers, fixing up the NextInRow pointers,
 *  swapping column indexes in MatrixElements, and swapping Markowitz 
 *  column counts.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  Col1  <input>  (int)
 *      Column index of one of the columns, becomes the smallest index.
 *  Col2  <input>  (int)
 *      Column index of the other column, becomes the largest index
 *
 *  Local variables:
 *  Row  (int)
 *      Row in which column elements are currently being exchanged.
 *  Col1Ptr  (ElementPtr)
 *      Pointer to an element in Col1.
 *  Col2Ptr  (ElementPtr)
 *      Pointer to an element in Col2.
 *  Element1  (ElementPtr)
 *      Pointer to the element in Col1 to be exchanged.
 *  Element2  (ElementPtr)
 *      Pointer to the element in Col2 to be exchanged.
 */


void
spcColExchange( MatrixPtr Matrix, int Col1, int Col2 )
{
ElementPtr  Col1Ptr, Col2Ptr;
int  Row;
ElementPtr  Element1, Element2;

/* Begin `spcColExchange'. */
    if (Col1 > Col2)  SWAP(int, Col1, Col2);

    Col1Ptr = Matrix->FirstInCol[Col1];
    Col2Ptr = Matrix->FirstInCol[Col2];
    while (Col1Ptr != NULL OR Col2Ptr != NULL)
    {
/* Exchange elements in rows while traveling from top to bottom. */
        if (Col1Ptr == NULL)
        {   Row = Col2Ptr->Row;
            Element1 = NULL;
            Element2 = Col2Ptr;
            Col2Ptr = Col2Ptr->NextInCol;
        }
        else if (Col2Ptr == NULL)
        {   Row = Col1Ptr->Row;
            Element1 = Col1Ptr;
            Element2 = NULL;
            Col1Ptr = Col1Ptr->NextInCol;
        }
        else if (Col1Ptr->Row < Col2Ptr->Row)
        {   Row = Col1Ptr->Row;
            Element1 = Col1Ptr;
            Element2 = NULL;
            Col1Ptr = Col1Ptr->NextInCol;
        }
        else if (Col1Ptr->Row > Col2Ptr->Row)
        {   Row = Col2Ptr->Row;
            Element1 = NULL;
            Element2 = Col2Ptr;
            Col2Ptr = Col2Ptr->NextInCol;
        }
        else   /* Col1Ptr->Row == Col2Ptr->Row */
        {   Row = Col1Ptr->Row;
            Element1 = Col1Ptr;
            Element2 = Col2Ptr;
            Col1Ptr = Col1Ptr->NextInCol;
            Col2Ptr = Col2Ptr->NextInCol;
        }

          ExchangeRowElements( Matrix, Col1, Element1, Col2, Element2, Row);
    }  /* end of while(Col1Ptr != NULL OR Col2Ptr != NULL) */

    if (Matrix->InternalVectorsAllocated)
        SWAP( int, Matrix->MarkowitzCol[Col1], Matrix->MarkowitzCol[Col2]);
    SWAP( ElementPtr, Matrix->FirstInCol[Col1], Matrix->FirstInCol[Col2]);
    SWAP( int, Matrix->IntToExtColMap[Col1], Matrix->IntToExtColMap[Col2]);
    SWAP( ElementPtr *, Matrix->Col_fast[Col1], Matrix->Col_fast[Col2]);
#if TRANSLATE
    Matrix->ExtToIntColMap[ Matrix->IntToExtColMap[Col1] ] = Col1;
    Matrix->ExtToIntColMap[ Matrix->IntToExtColMap[Col2] ] = Col2;
#endif

    spColInd(Matrix, Col1);
    spColInd(Matrix, Col2);
    return;
}


/* This flag turns on the new search method, but it does not work.  Can't
   figure out why */
/* #define NEW_M */



/*
 *  EXCHANGE TWO ELEMENTS IN A COLUMN
 *
 *  Performs all required operations to exchange two elements in a column.
 *  Those operations are: restring NextInCol pointers and swapping row indexes
 *  in the MatrixElements.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  Row1  <input>  (int)
 *      Row of top element to be exchanged.
 *  Element1  <input>  (ElementPtr)
 *      Pointer to top element to be exchanged.
 *  Row2  <input>  (int)
 *      Row of bottom element to be exchanged.
 *  Element2  <input>  (ElementPtr)
 *      Pointer to bottom element to be exchanged.
 *  Column <input>  (int)
 *      Column that exchange is to take place in.
 *
 *  >>> Local variables:
 *  ElementAboveRow1  (ElementPtr *)
 *      Location of pointer which points to the element above Element1. This
 *      pointer is modified so that it points to correct element on exit.
 *  ElementAboveRow2  (ElementPtr *)
 *      Location of pointer which points to the element above Element2. This
 *      pointer is modified so that it points to correct element on exit.
 *  ElementBelowRow1  (ElementPtr)
 *      Pointer to element below Element1.
 *  ElementBelowRow2  (ElementPtr)
 *      Pointer to element below Element2.
 *  pElement  (ElementPtr)
 *      Pointer used to traverse the column.
 */

static void
ExchangeColElements( MatrixPtr Matrix, int Row1, ElementPtr Element1, int Row2, ElementPtr Element2, int Column )
{
ElementPtr  *ElementAboveRow1, *ElementAboveRow2;
ElementPtr  ElementAbove1, ElementAbove2;
ElementPtr  ElementBelowRow1, ElementBelowRow2, ElementRight, ElementAbove;
ElementPtr  pElement;

/* Begin `ExchangeColElements'. */
/* Search to find the ElementAboveRow1. */

    if (Row1 > Row2) SWAP (int, Row1, Row2);
    ElementAbove = Matrix->Col_fast[Column][f_ind(Matrix, Column, Row1)];
    if (ElementAbove == NULL || ElementAbove->Col != Column ||
        ElementAbove->Row >= Row1) {
          ElementAboveRow1 = &(Matrix->FirstInCol[Column]);
          ElementAbove1 = NULL;
    }
    else {
      ElementAboveRow1 = &(ElementAbove->NextInCol);
      ElementAbove1 = ElementAbove;
    }
    pElement = *ElementAboveRow1;
#ifdef NEW_M
    if (pElement != NULL) {
      while (pElement->NextInCol != NULL && pElement->NextInCol->Row < Row1) 
        pElement = pElement->NextInCol;
      ElementAbove1 = pElement;
      ElementAboveRow1 = &(pElement->NextInCol);
    }
#else
    while (pElement != NULL && pElement->Row < Row1)
    {   ElementAboveRow1 = &(pElement->NextInCol);
        ElementAbove1 = pElement;
        pElement = *ElementAboveRow1;
    }
#endif


    ElementAbove = Matrix->Col_fast[Column][f_ind(Matrix, Column, Row2)];
    if (ElementAbove == NULL || ElementAbove->Col != Column ||
      ElementAbove->Row >= Row2) {
        ElementAboveRow2 = ElementAboveRow1;
        ElementAbove2 = ElementAbove1;
    }
    else {
      ElementAboveRow2 = &(ElementAbove->NextInCol);
      ElementAbove2 = ElementAbove;
    }
    pElement = *ElementAboveRow2;
#ifdef NEW_M
    if (pElement != NULL) {
      while (pElement->NextInCol != NULL && pElement->NextInCol->Row < Row2) 
        pElement = pElement->NextInCol;
      ElementAbove2 = pElement;
      ElementAboveRow2 = &(pElement->NextInCol);
    }
#else
    while (pElement != NULL && pElement->Row < Row2)
    {   ElementAboveRow2 = &(pElement->NextInCol);
        ElementAbove2 = pElement;
        pElement = *ElementAboveRow2;
    }
#endif

    if (Element1 != NULL)
    {
        remove_fast_col_index (Matrix, Row1, Column, ElementAbove1, *ElementAboveRow1);
        ElementBelowRow1 = Element1->NextInCol;
        if (Element2 == NULL)
        {
/* Element2 does not exist, move Element1 down to Row2. */
            if ( ElementBelowRow1 != NULL AND ElementBelowRow1->Row < Row2 )
            {
/* Element1 must be removed from linked list and moved. */
                *ElementAboveRow1 = ElementBelowRow1;
                Element1->NextInCol = *ElementAboveRow2;
                *ElementAboveRow2 = Element1;
            }
        }
        else
        {
            remove_fast_col_index (Matrix, Row2, Column, ElementAbove2, *ElementAboveRow2);
/* Element2 does exist, and the two elements must be exchanged. */
            if ( ElementBelowRow1->Row == Row2)
            {
/* Element2 is just below Element1, exchange them. */
                Element1->NextInCol = Element2->NextInCol;
                Element2->NextInCol = Element1;
                *ElementAboveRow1 = Element2;
            }
            else
            {
                ElementBelowRow2 = Element2->NextInCol;
/* Switch Element1 and Element2. */
                *ElementAboveRow1 = Element2;
                Element2->NextInCol = ElementBelowRow1;
                *ElementAboveRow2 = Element1;
                Element1->NextInCol = ElementBelowRow2;
            }
            Element2->Row = Row1;
            add_fast_col_index (Matrix, Row1, Column, Element2);
        }
        Element1->Row = Row2;
        add_fast_col_index (Matrix, Row2, Column, Element1);
    }
    else
    {
/* Element1 does not exist. */
        remove_fast_col_index (Matrix, Row2, Column, ElementAbove2, *ElementAboveRow2);

        if (*ElementAboveRow1 != NULL && (*ElementAboveRow1)->Row < Row1)
          ElementBelowRow1 = (*ElementAboveRow1)->NextInCol;
        else
          ElementBelowRow1 = *ElementAboveRow1;

        if (ElementBelowRow1 != NULL && ElementBelowRow1->Row != Row2)
        {
          *ElementAboveRow2 = Element2->NextInCol;
          *ElementAboveRow1 = Element2;
          Element2->NextInCol = ElementBelowRow1;
        }
        Element2->Row = Row1;
        add_fast_col_index (Matrix, Row1, Column, Element2);
    }
    return;
}







/*
 *  EXCHANGE TWO ELEMENTS IN A ROW
 *
 *  Performs all required operations to exchange two elements in a row.
 *  Those operations are: restring NextInRow pointers and swapping column
 *  indexes in the MatrixElements.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  Col1  <input>  (int)
 *      Col of left-most element to be exchanged.
 *  Element1  <input>  (ElementPtr)
 *      Pointer to left-most element to be exchanged.
 *  Col2  <input>  (int)
 *      Col of right-most element to be exchanged.
 *  Element2  <input>  (ElementPtr)
 *      Pointer to right-most element to be exchanged.
 *  Row <input>  (int)
 *      Row that exchange is to take place in.
 *
 *  >>> Local variables:
 *  ElementLeftOfCol1  (ElementPtr *)
 *      Location of pointer which points to the element to the left of
 *      Element1. This pointer is modified so that it points to correct
 *      element on exit.
 *  ElementLeftOfCol2  (ElementPtr *)
 *      Location of pointer which points to the element to the left of
 *      Element2. This pointer is modified so that it points to correct
 *      element on exit.
 *  ElementRightOfCol1  (ElementPtr)
 *      Pointer to element right of Element1.
 *  ElementRightOfCol2  (ElementPtr)
 *      Pointer to element right of Element2.
 *  pElement  (ElementPtr)
 *      Pointer used to traverse the row.
 */

static void
ExchangeRowElements( MatrixPtr Matrix, int Col1, ElementPtr Element1, int Col2, ElementPtr Element2, int Row )
{
    ElementPtr  *ElementLeftOfCol1, *ElementLeftOfCol2;
    ElementPtr  ElementLeftOf1, ElementLeftOf2;
    ElementPtr  ElementRightOfCol1, ElementRightOfCol2, ElementLeft;
    ElementPtr  pElement;

/* Begin `ExchangeRowElements'. */
/* Search to find the ElementLeftOfCol1. */

    if (Col1 > Col2) SWAP (int, Col1, Col2);
    ElementLeft = Matrix->Row_fast[Row][f_ind(Matrix, Row, Col1)];
    if (ElementLeft == NULL || ElementLeft->Row != Row ||
        ElementLeft->Col >= Col1) {
          ElementLeftOfCol1 = &(Matrix->FirstInRow[Row]);
          ElementLeftOf1 = NULL;
    }
    else {
      ElementLeftOfCol1 = &(ElementLeft->NextInRow);
      ElementLeftOf1 = ElementLeft;
    }
    pElement = *ElementLeftOfCol1;
#ifdef NEW_M
    if (pElement != NULL) {
      while (pElement->NextInRow != NULL && pElement->NextInRow->Col < Col1) 
        pElement = pElement->NextInRow;
      ElementLeftOf1 = pElement;
      ElementLeftOfCol1 = &(pElement->NextInRow);
    }
#else
    while (pElement != NULL && pElement->Col < Col1)
    {   ElementLeftOfCol1 = &(pElement->NextInRow);
        ElementLeftOf1 = pElement;
        pElement = *ElementLeftOfCol1;
    }
#endif

    ElementLeft = Matrix->Row_fast[Row][f_ind(Matrix, Row, Col2)];
    if (ElementLeft == NULL || ElementLeft->Row != Row ||
      ElementLeft->Col >= Col2) {
        ElementLeftOfCol2 = ElementLeftOfCol1;
          ElementLeftOf2 =  ElementLeftOf1;
    }
    else {
      ElementLeftOfCol2 = &(ElementLeft->NextInRow);
      ElementLeftOf2 = ElementLeft;
    }
    pElement = *ElementLeftOfCol2;
#ifdef NEW_M
    if (pElement != NULL) {
      while (pElement->NextInRow != NULL && pElement->NextInRow->Col < Col2)
        pElement = pElement->NextInRow;
      ElementLeftOf2 = pElement;
      ElementLeftOfCol2 = &(pElement->NextInRow);
    }
#else
    while (pElement != NULL && pElement->Col < Col2)
    {   ElementLeftOfCol2 = &(pElement->NextInRow);
        ElementLeftOf2 = pElement;
        pElement = *ElementLeftOfCol2;
    }
#endif

    if (Element1 != NULL)
    {
        remove_fast_row_index (Matrix, Row, Col1, ElementLeftOf1, *ElementLeftOfCol1);
        ElementRightOfCol1 = Element1->NextInRow;
        if (Element2 == NULL)
        {
/* Element2 does not exist, move Element1 to right to Col2. */
            if ( ElementRightOfCol1 != NULL AND ElementRightOfCol1->Col < Col2 )
            {
/* Element1 must be removed from linked list and moved. */
                *ElementLeftOfCol1 = ElementRightOfCol1;
                Element1->NextInRow = *ElementLeftOfCol2;
                *ElementLeftOfCol2 = Element1;
            }
        }
        else
        {
/* Element2 does exist, and the two elements must be exchanged. */
            remove_fast_row_index (Matrix, Row, Col2, ElementLeftOf2, *ElementLeftOfCol2);
            if ( ElementRightOfCol1->Col == Col2)
            {
/* Element2 is just right of Element1, exchange them. */
                Element1->NextInRow = Element2->NextInRow;
                Element2->NextInRow = Element1;
                *ElementLeftOfCol1 = Element2;
            }
            else
            {
                ElementRightOfCol2 = Element2->NextInRow;
/* Switch Element1 and Element2. */
                *ElementLeftOfCol1 = Element2;
                Element2->NextInRow = ElementRightOfCol1;
                *ElementLeftOfCol2 = Element1;
                Element1->NextInRow = ElementRightOfCol2;
            }
            Element2->Col = Col1;
            add_fast_row_index (Matrix, Row, Col1, Element2);
        }
        Element1->Col = Col2;
        add_fast_row_index (Matrix, Row, Col2, Element1);
    }
    else
    {
/* Element1 does not exist. */
        remove_fast_row_index (Matrix, Row, Col2, ElementLeftOf2, *ElementLeftOfCol2);
        if (*ElementLeftOfCol1 != NULL && (*ElementLeftOfCol1)->Col < Col1)
          ElementRightOfCol1 = (*ElementLeftOfCol1)->NextInRow;
        else
          ElementRightOfCol1 = *ElementLeftOfCol1;

        if (ElementRightOfCol1 != NULL && ElementRightOfCol1->Col != Col2) {
          *ElementLeftOfCol2 = Element2->NextInRow;
          *ElementLeftOfCol1 = Element2;
          Element2->NextInRow = ElementRightOfCol1;
        }
        Element2->Col = Col1;
        add_fast_row_index (Matrix, Row, Col1, Element2);
    }
    return;
}




/*
 *  PERFORM ROW AND COLUMN ELIMINATION ON REAL MATRIX
 *
 *  Eliminates a single row and column of the matrix and leaves single row of
 *  the upper triangular matrix and a single column of the lower triangular
 *  matrix in its wake.  Uses Gauss's method.
 *
 *  >>> Argument:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  pPivot  <input>  (ElementPtr)
 *      Pointer to the current pivot.
 *
 *  >>> Local variables:
 *  pLower  (ElementPtr)
 *      Points to matrix element in lower triangular column.
 *  pSub (ElementPtr)
 *      Points to elements in the reduced submatrix.
 *  Row  (int)
 *      Row index.
 *  pUpper  (ElementPtr)
 *      Points to matrix element in upper triangular row.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

void
RealRowColElimination( MatrixPtr Matrix, ElementPtr pPivot , int Step)
{
#if REAL
ElementPtr  pSub, pFast;
int  Row;
ElementPtr  pLower, pUpper;

/* Begin `RealRowColElimination'. */

/* Test for zero pivot. */
      if (ABS(pPivot->Real) == 0.0)
      {   (void)MatrixIsSingular( Matrix, pPivot->Row );
          return;
      }
      pPivot->Real = 1.0 / pPivot->Real;

    pUpper = pPivot->NextInRow;

    while (pUpper != NULL)
    {
/* Calculate upper triangular element. */
        pUpper->Real *= pPivot->Real;

        pSub = pUpper;
        pLower = pPivot->NextInCol;
        while (pLower != NULL)
        {
            Row = pLower->Row;

/* Find element in row that lines up with current lower triangular element. */
            if (pSub != NULL) {
              if (pSub->Row > Row) {
                pFast = Matrix->Col_fast[pSub->Col][f_ind(Matrix, pSub->Col,Row)];
                if (pFast != NULL && pFast->Col == pSub->Col && pFast->Row <=Row)
                  pSub = pFast;
                else
                  pSub = pUpper;
              }
              while (pSub->NextInCol != NULL && pSub->NextInCol->Row <= Row)
                pSub = pSub->NextInCol;
            }

/* Test to see if desired element was not found, if not, create fill-in. */
            if (pSub == NULL OR pSub->Row != Row)
            {   
                if (pSub == NULL)
                  pSub = CreateFillin( Matrix, Row, pUpper->Col, Step, NULL);
                else
                  pSub = CreateFillin( Matrix, Row, pUpper->Col, Step, &pSub->NextInCol);
                pSub->Real -= pUpper->Real * pLower->Real;
            }
            else {
              pSub->Real -= pUpper->Real * pLower->Real;
            }
            pLower = pLower->NextInCol;
        }
        pUpper = pUpper->NextInRow;
    }
    return;
#endif /* REAL */
}









/*
 *  PERFORM ROW AND COLUMN ELIMINATION ON COMPLEX MATRIX
 *
 *  Eliminates a single row and column of the matrix and leaves single row of
 *  the upper triangular matrix and a single column of the lower triangular
 *  matrix in its wake.  Uses Gauss's method.
 *
 *  >>> Argument:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  pPivot  <input>  (ElementPtr)
 *      Pointer to the current pivot.
 *
 *  >>> Local variables:
 *  pLower  (ElementPtr)
 *      Points to matrix element in lower triangular column.
 *  pSub (ElementPtr)
 *      Points to elements in the reduced submatrix.
 *  Row  (int)
 *      Row index.
 *  pUpper  (ElementPtr)
 *      Points to matrix element in upper triangular row.
 *
 *  Possible errors:
 *  spNO_MEMORY
 */

static void
ComplexRowColElimination( MatrixPtr Matrix, ElementPtr pPivot, int Step )
{
#if spCOMPLEX
ElementPtr  pSub;
int  Row;
ElementPtr  pLower, pUpper;

/* Begin `ComplexRowColElimination'. */

/* Test for zero pivot. */
    if (ELEMENT_MAG(pPivot) == 0.0)
    {   (void)MatrixIsSingular( Matrix, pPivot->Row );
        return;
    }
    CMPLX_RECIPROCAL(*pPivot, *pPivot);

    pUpper = pPivot->NextInRow;
    while (pUpper != NULL)
    {
/* Calculate upper triangular element. */
/* Cmplx expr: *pUpper = *pUpper * (1.0 / *pPivot). */
        CMPLX_MULT_ASSIGN(*pUpper, *pPivot);

        pSub = pUpper->NextInCol;
        pLower = pPivot->NextInCol;
        while (pLower != NULL)
        {   Row = pLower->Row;

/* Find element in row that lines up with current lower triangular element. */
            while (pSub != NULL AND pSub->Row < Row)
                pSub = pSub->NextInCol;

/* Test to see if desired element was not found, if not, create fill-in. */
            if (pSub == NULL OR pSub->Row > Row)
            {   pSub = CreateFillin( Matrix, Row, pUpper->Col, Step, NULL);
                if (pSub == NULL)
                {   Matrix->Error = spNO_MEMORY;
                    return;
                }
            }

/* Cmplx expr: pElement -= *pUpper * pLower. */
            CMPLX_MULT_SUBT_ASSIGN(*pSub, *pUpper, *pLower);
            pSub = pSub->NextInCol;
            pLower = pLower->NextInCol;
        }
        pUpper = pUpper->NextInRow;
    }
    return;
#endif /* spCOMPLEX */
}





/*
 *  UPDATE MARKOWITZ NUMBERS
 *
 *  Updates the Markowitz numbers after a row and column have been eliminated.
 *  Also updates singleton count.
 *
 *  >>> Argument:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  pPivot  <input>  (ElementPtr)
 *      Pointer to the current pivot.
 *
 *  >>> Local variables:
 *  Row  (int)
 *      Row index.
 *  Col  (int)
 *      Column index.
 *  ColPtr  (ElementPtr)
 *      Points to matrix element in upper triangular column.
 *  RowPtr  (ElementPtr)
 *      Points to matrix element in lower triangular row.
 */

static void
UpdateMarkowitzNumbers( MatrixPtr Matrix, ElementPtr pPivot )
{
int  Row, Col;
ElementPtr  ColPtr, RowPtr;
int *MarkoRow = Matrix->MarkowitzRow, *MarkoCol = Matrix->MarkowitzCol;
double Product;

/* Begin `UpdateMarkowitzNumbers'. */

/* Update Markowitz numbers. */
    for (ColPtr = pPivot->NextInCol; ColPtr != NULL; ColPtr = ColPtr->NextInCol)
    {   Row = ColPtr->Row;
        --MarkoRow[Row];

/* Form Markowitz product while being cautious of overflows. */
        if ((MarkoRow[Row] > LARGEST_SHORT_INTEGER AND MarkoCol[Row] != 0) OR
            (MarkoCol[Row] > LARGEST_SHORT_INTEGER AND MarkoRow[Row] != 0))
        {   Product = MarkoCol[Row] * MarkoRow[Row];
            if (Product >= (double)LARGEST_LONG_INTEGER)
                Matrix->MarkowitzProd[Row] = LARGEST_LONG_INTEGER;
            else
                Matrix->MarkowitzProd[Row] = Product;
        }
        else Matrix->MarkowitzProd[Row] = MarkoRow[Row] * MarkoCol[Row];
        if (MarkoRow[Row] == 0)
            Matrix->Singletons++;
    }

    for (RowPtr = pPivot->NextInRow; RowPtr != NULL; RowPtr = RowPtr->NextInRow)
    {   Col = RowPtr->Col;
        --MarkoCol[Col];

/* Form Markowitz product while being cautious of overflows. */
        if ((MarkoRow[Col] > LARGEST_SHORT_INTEGER AND MarkoCol[Col] != 0) OR
            (MarkoCol[Col] > LARGEST_SHORT_INTEGER AND MarkoRow[Col] != 0))
        {   Product = MarkoCol[Col] * MarkoRow[Col];
            if (Product >= (double)LARGEST_LONG_INTEGER)
                Matrix->MarkowitzProd[Col] = LARGEST_LONG_INTEGER;
            else
                Matrix->MarkowitzProd[Col] = Product;
        }
        else Matrix->MarkowitzProd[Col] = MarkoRow[Col] * MarkoCol[Col];
        if ((MarkoCol[Col] == 0) AND (MarkoRow[Col] != 0))
            Matrix->Singletons++;
    }
    return;
}








/*
 *  CREATE FILL-IN
 *
 *  This routine is used to create fill-ins and splice them into the
 *  matrix.
 *
 *  >>> Returns:
 *  Pointer to fill-in.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  Col  <input>  (int)
 *      Column index for element.
 *  Row  <input>  (int)
 *      Row index for element.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.
 *  ppElementAbove  (ElementPtr *)
 *      This contains the address of the pointer to the element just above the
 *      one being created. It is used to speed the search and it is updated
 *      with address of the created element.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

static ElementPtr
CreateFillin( MatrixPtr Matrix, int Row, int Col, int Step, ElementPtr* ab)
{
ElementPtr  pElement, *ppElementAbove;

/* Begin `CreateFillin'. */

/* Find Element above fill-in. */
#ifdef CHECK_IND
    spCheckInd(Matrix, "Top CreateFillin");
#endif
    if (ab)
      ppElementAbove = ab;
    else
      ppElementAbove = &Matrix->FirstInCol[Col];
    pElement = *ppElementAbove;
    while (pElement != NULL)
    {   if (pElement->Row < Row)
        {   ppElementAbove = &pElement->NextInCol;
            pElement = *ppElementAbove;
        }
        else break; /* while loop */
    }
    if (pElement != NULL && pElement->Row == Row) return pElement;

/* End of search, create the element. */
    pElement = spcCreateElement( Matrix, Row, Col, ppElementAbove, YES );
    pElement->Fillin = Step;
    if (pElement->Col < Step) {
      printf ("Assumption about Col >= Step false\n");
      exit(-1);
    }
    add_fast_row_index(Matrix, Row, Col, pElement);
    add_fast_col_index(Matrix, Row, Col, pElement);

/* Update Markowitz counts and products. */
    Matrix->MarkowitzProd[Row] = ++Matrix->MarkowitzRow[Row] *
                                   Matrix->MarkowitzCol[Row];
    if ((Matrix->MarkowitzRow[Row] == 1) AND (Matrix->MarkowitzCol[Row] != 0))
        Matrix->Singletons--;
    Matrix->MarkowitzProd[Col] = ++Matrix->MarkowitzCol[Col] *
                                   Matrix->MarkowitzRow[Col];
    if ((Matrix->MarkowitzRow[Col] != 0) AND (Matrix->MarkowitzCol[Col] == 1))
        Matrix->Singletons--;

    return pElement;
}

/* Fast matrix index update routines: */

void add_fast_col_index (MatrixPtr Matrix, int Row, int Col,
                                ElementPtr pElement)
{
    int ind;

    ind = f_ind(Matrix, Col, Row)+1;
    while (ind < Matrix->Indsize &&
          (Matrix->Col_fast[Col][ind] == NULL ||
           Matrix->Col_fast[Col][ind]->Row < Row)) {
      Matrix->Col_fast[Col][ind++] = pElement;
    }
    return;
}

void add_fast_row_index (MatrixPtr Matrix, int Row, int Col,
                                ElementPtr pElement)
{
    int ind;

    ind = f_ind(Matrix, Row, Col)+1;
    while (ind < Matrix->Indsize &&
          (Matrix->Row_fast[Row][ind] == NULL ||
           Matrix->Row_fast[Row][ind]->Col < Col)) {
      Matrix->Row_fast[Row][ind++] = pElement;
    }
    return;
}

static void remove_fast_col_index (MatrixPtr Matrix, int Row, int Col,
                                   ElementPtr pAbove, ElementPtr pElement)
{
    int ind;

    ind = f_ind(Matrix, Col, Row);
    while (ind < Matrix->Indsize && (!Matrix->Col_fast[Col][ind] ||
           Matrix->Col_fast[Col][ind]->Row <= Row)) {
      if (Matrix->Col_fast[Col][ind] == pElement)
        Matrix->Col_fast[Col][ind] = pAbove;
      ind++;
    }

    return;
}

static void remove_fast_row_index (MatrixPtr Matrix, int Row, int Col,
                                   ElementPtr pLeft, ElementPtr pElement)
{
    int ind;

    ind = f_ind(Matrix, Row, Col);
    while (ind < Matrix->Indsize && (!Matrix->Row_fast[Row][ind] ||
           Matrix->Row_fast[Row][ind]->Col <= Col)) {
      if (Matrix->Row_fast[Row][ind] == pElement)
        Matrix->Row_fast[Row][ind] = pLeft;
      ind++;
    }

    return;
}



/*
 *  ZERO PIVOT ENCOUNTERED
 *
 *  This routine is called when a singular matrix is found.  It then
 *  records the current row and column and exits.
 *
 *  >>> Returned:
 *  The error code spSINGULAR or spZERO_DIAG is returned.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Step  <input>  (int)
 *      Index of diagonal that is zero.
 */

static int
MatrixIsSingular( MatrixPtr Matrix, int Step )
{
/* Begin `MatrixIsSingular'. */

    Matrix->SingularRow = Matrix->IntToExtRowMap[ Step ];
    Matrix->SingularCol = Matrix->IntToExtColMap[ Step ];
    return (Matrix->Error = spSINGULAR);
}


static int
ZeroPivot( MatrixPtr Matrix, int Step )
{
/* Begin `ZeroPivot'. */

    Matrix->SingularRow = Matrix->IntToExtRowMap[ Step ];
    Matrix->SingularCol = Matrix->IntToExtColMap[ Step ];
    return (Matrix->Error = spZERO_DIAG);
}


int
SMPmergeConstraint (MatrixPtr Matrix, double *rhs, int Row, double Value)
{
    ElementPtr cElem, rElem, baseElem, curr_E;
    double pivot, scale;
    int currRow, r_val, R;

    r_val = 0;
    return r_val;
    if (Matrix->RowsLinked == NO)
      spcLinkRows(Matrix);

    pivot = Matrix->Diag[Row]->Real;
/*    printf ("Merging constraint for row: %d\n",Row); */
    if (pivot != 0) {
      for (cElem = Matrix->FirstInCol[Row] ; cElem ; cElem = cElem->NextInCol) {
/*              printf ("Doing constraint for row: %d\n",cElem->Row); */
        currRow = cElem->Row;
        if (currRow != Row) {
          rhs[currRow] -= cElem->Real*Value;
          cElem->Real = 0;
        }
      }
    }
    return r_val;
}



#if (ANNOTATE == FULL)

/*
 *
 *  WRITE STATUS
 *
 *  Write a summary of important variables to standard output.
 */

static void
WriteStatus( Matrix, Step )

MatrixPtr Matrix;
int Step;
{
int  I;

/* Begin `WriteStatus'. */

    printf("Step = %1d   ", Step);
    printf("Pivot found at %1d,%1d using ", Matrix->PivotsOriginalRow,
                                            Matrix->PivotsOriginalCol);
    switch(Matrix->PivotSelectionMethod)
    {   case 's': printf("SearchForSingleton\n");  break;
        case 'q': printf("QuicklySearchDiagonal\n");  break;
        case 'd': printf("SearchDiagonal\n");  break;
        case 'e': printf("SearchEntireMatrix\n");  break;
    }

    printf("MarkowitzRow     = ");
    for (I = 1; I <= Matrix->Size; I++)
        printf("%2d  ", Matrix->MarkowitzRow[I]);
    printf("\n");

    printf("MarkowitzCol     = ");
    for (I = 1; I <= Matrix->Size; I++)
        printf("%2d  ", Matrix->MarkowitzCol[I]);
    printf("\n");

    printf("MarkowitzProduct = ");
    for (I = 1; I <= Matrix->Size; I++)
        printf("%2d  ", Matrix->MarkowitzProd[I]);
    printf("\n");

    printf("Singletons = %2d\n", Matrix->Singletons);

    printf("IntToExtRowMap     = ");
    for (I = 1; I <= Matrix->Size; I++)
        printf("%2d  ", Matrix->IntToExtRowMap[I]);
    printf("\n");

    printf("IntToExtColMap     = ");
    for (I = 1; I <= Matrix->Size; I++)
        printf("%2d  ", Matrix->IntToExtColMap[I]);
    printf("\n");

    printf("ExtToIntRowMap     = ");
    for (I = 1; I <= Matrix->ExtSize; I++)
        printf("%2d  ", Matrix->ExtToIntRowMap[I]);
    printf("\n");

    printf("ExtToIntColMap     = ");
    for (I = 1; I <= Matrix->ExtSize; I++)
        printf("%2d  ", Matrix->ExtToIntColMap[I]);
    printf("\n\n");

/*  spPrint((char *)Matrix, NO, YES);  */

    return;

}
#endif /* ANNOTATE == FULL */

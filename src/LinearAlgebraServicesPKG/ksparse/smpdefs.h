/*
The original code in this file is Copyright (c) 1985-1991 The Regents of the
University of California and is under the Spice 3f5 BSD Copyright.

All additions and changes are under the following:
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
*/

#ifndef KSPARSE_SMP
#define KSPARSE_SMP

typedef  char SMPmatrix;
typedef  struct MatrixElement  *SMPelement;

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "complex.h"
#include <stdio.h>

#ifdef __STDC__
int SMPaddElt( SMPmatrix *, int , int , double );
void SMPcClear( SMPmatrix *);
int SMPcLUfac( SMPmatrix *, double );
int SMPcProdDiag( SMPmatrix *, SPcomplex *, int *);
int SMPcReorder( SMPmatrix * , double , double , int *);
int SMPcSolve( SMPmatrix *, double [], double [], double [], double []);
void SMPclear( SMPmatrix *);
void SMPcolSwap( SMPmatrix * , int , int );
void SMPdestroy( SMPmatrix *);
int SMPfillup( SMPmatrix * );
SMPelement * SMPfindElt( SMPmatrix *, int , int , int );
void SMPgetError( SMPmatrix *, int *, int *);
int SMPluFac( SMPmatrix *, double , double );
double * SMPmakeElt( SMPmatrix * , int , int );
int SMPmatSize( SMPmatrix *);
int SMPnewMatrix( SMPmatrix ** );
int SMPnewNode( int , SMPmatrix *);
int SMPpreOrder( SMPmatrix *);
void SMPprint( SMPmatrix * , char *);
int SMPreorder( SMPmatrix * , double , double , double );
void SMProwSwap( SMPmatrix * , int , int );
int SMPsolve( SMPmatrix *, double [], double []);
#else /* stdc */
int SMPaddElt();
void SMPcClear();
int SMPcLUfac();
int SMPcProdDiag();
int SMPcReorder();
int SMPcSolve();
void SMPclear();
void SMPcolSwap();
void SMPdestroy();
int SMPfillup();
SMPelement * SMPfindElt();
void SMPgetError();
int SMPluFac();
double * SMPmakeElt();
int SMPmatSize();
int SMPnewMatrix();
int SMPnewNode();
int SMPpreOrder();
void SMPprint();
int SMPreorder();
void SMProwSwap();
int SMPsolve();
#endif /* stdc */

#endif /*KSPARSE_SMP*/

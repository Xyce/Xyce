/*
The original code in this file is Copyright (c) 1985-1991 The Regents of the
University of California and is under the Spice 3f5 BSD Copyright.

All additions and changes are under the following:
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
*/

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#ifndef KSPARSE_UTIL
#define KSPARSE_UTIL

#ifdef CHILE
/* #define FREE(x) {if (x) {free((char *)(x));(x) = 0;}} */
#define FREE(x) {if (x) {free((char *)(x));}}
#endif

#ifdef CHILE
#define MALLOC(x) tmalloc((unsigned)(x))
#define REALLOC(x,y) trealloc((char *)(x),(unsigned)(y))
#else
char *MALLOC(unsigned);
char *REALLOC (char *, unsigned);
void FREE(char *);
#endif

#define ZERO(PTR,TYPE)	(bzero((PTR),sizeof(TYPE)))

#include <stdlib.h>

extern char *trealloc();
extern char *tmalloc();

#define TRUE 1
#define FALSE 0

#ifdef DEBUG
#define DEBUGMSG(textargs) printf(textargs)
#else
#define DEBUGMSG(testargs)
#endif

#ifdef HAS_NOINLINE
#define FABS(a) fabs(a)
double fabs();
#else
#define FABS(a) ( ((a)<0) ? -(a) : (a) )
#endif

/* XXX Move these into the above ifdef someday */
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ( b >= 0 ? (a >= 0 ? a : - a) : (a >= 0 ? - a : a))

#define ABORT() fflush(stderr);fflush(stdout);abort();

#define ERROR(CODE,MESSAGE)	{					      \
	errMsg = MALLOC(strlen(MESSAGE) + 1);				      \
	strcpy(errMsg, (MESSAGE));					      \
	return (CODE);							      \
	}

#define	NEW(TYPE)	((TYPE *) MALLOC(sizeof(TYPE)))
#define	NEWN(TYPE,COUNT) ((TYPE *) MALLOC(sizeof(TYPE) * (COUNT)))
#define	NEW_SM(TYPE)	((TYPE *) MALLOC_SM(sizeof(TYPE)))
#define	NEWN_SM(TYPE,COUNT) ((TYPE *) MALLOC_SM(sizeof(TYPE) * (COUNT)))

#endif /*KSPARSE_UTIL*/

#define	R_NORM(A,B) {							      \
	if ((A) == 0.0) {						      \
	    (B) = 0;							      \
	} else {							      \
	    while (FABS(A) > 1.0) {					      \
		(B) += 1;						      \
		(A) /= 2.0;						      \
	    }								      \
	    while (FABS(A) < 0.5) {					      \
		(B) -= 1;						      \
		(A) *= 2.0;						      \
	    }								      \
	}								      \
    }

#define MOD_01 here->method.curr ^= (long) 00000000000000000000001;
#define MOD_02 here->method.curr ^= (long) 00000000000000000000002;
#define MOD_03 here->method.curr ^= (long) 00000000000000000000004;
#define MOD_04 here->method.curr ^= (long) 00000000000000000000010;
#define MOD_05 here->method.curr ^= (long) 00000000000000000000020;
#define MOD_06 here->method.curr ^= (long) 00000000000000000000040;
#define MOD_07 here->method.curr ^= (long) 00000000000000000000100;
#define MOD_08 here->method.curr ^= (long) 00000000000000000000200;
#define MOD_09 here->method.curr ^= (long) 00000000000000000000400;
#define MOD_10 here->method.curr ^= (long) 00000000000000000001000;
#define MOD_11 here->method.curr ^= (long) 00000000000000000002000;
#define MOD_12 here->method.curr ^= (long) 00000000000000000004000;
#define MOD_13 here->method.curr ^= (long) 00000000000000000010000;
#define MOD_14 here->method.curr ^= (long) 00000000000000000020000;
#define MOD_15 here->method.curr ^= (long) 00000000000000000040000;
#define MOD_16 here->method.curr ^= (long) 00000000000000000100000;
#define MOD_17 here->method.curr ^= (long) 00000000000000000200000;
#define MOD_18 here->method.curr ^= (long) 00000000000000000400000;
#define MOD_19 here->method.curr ^= (long) 00000000000000001000000;
#define MOD_20 here->method.curr ^= (long) 00000000000000002000000;
#define MOD_21 here->method.curr ^= (long) 00000000000000004000000;
#define MOD_22 here->method.curr ^= (long) 00000000000000010000000;
#define MOD_23 here->method.curr ^= (long) 00000000000000020000000;
#define MOD_24 here->method.curr ^= (long) 00000000000000040000000;
#define MOD_25 here->method.curr ^= (long) 00000000000000100000000;
#define MOD_26 here->method.curr ^= (long) 00000000000000200000000;
#define MOD_27 here->method.curr ^= (long) 00000000000000400000000;
#define MOD_28 here->method.curr ^= (long) 00000000000001000000000;
#define MOD_29 here->method.curr ^= (long) 00000000000002000000000;
#define MOD_30 here->method.curr ^= (long) 00000000000004000000000;
#define MOD_31 here->method.curr ^= (long) 00000000000010000000000;
#define MOD_32 here->method.curr ^= (long) 00000000000020000000000;
#define MOD_33 here->method.curr ^= (long) 00000000000040000000000;
#define MOD_34 here->method.curr ^= (long) 00000000000100000000000;
#define MOD_35 here->method.curr ^= (long) 00000000000200000000000;
#define MOD_36 here->method.curr ^= (long) 00000000000400000000000;
#define MOD_37 here->method.curr ^= (long) 00000000001000000000000;
#define MOD_38 here->method.curr ^= (long) 00000000002000000000000;
#define MOD_39 here->method.curr ^= (long) 00000000004000000000000;
#define MOD_40 here->method.curr ^= (long) 00000000010000000000000;
#define MOD_41 here->method.curr ^= (long) 00000000020000000000000;
#define MOD_42 here->method.curr ^= (long) 00000000040000000000000;
#define MOD_43 here->method.curr ^= (long) 00000000100000000000000;
#define MOD_44 here->method.curr ^= (long) 00000000200000000000000;
#define MOD_45 here->method.curr ^= (long) 00000000400000000000000;
#define MOD_46 here->method.curr ^= (long) 00000001000000000000000;
#define MOD_47 here->method.curr ^= (long) 00000002000000000000000;
#define MOD_48 here->method.curr ^= (long) 00000004000000000000000;
#define MOD_49 here->method.curr ^= (long) 00000010000000000000000;
#define MOD_50 here->method.curr ^= (long) 00000020000000000000000;
#define MOD_51 here->method.curr ^= (long) 00000040000000000000000;
#define MOD_52 here->method.curr ^= (long) 00000100000000000000000;
#define MOD_53 here->method.curr ^= (long) 00000200000000000000000;
#define MOD_54 here->method.curr ^= (long) 00000400000000000000000;
#define MOD_55 here->method.curr ^= (long) 00001000000000000000000;
#define MOD_56 here->method.curr ^= (long) 00002000000000000000000;
#define MOD_57 here->method.curr ^= (long) 00004000000000000000000;
#define MOD_58 here->method.curr ^= (long) 00010000000000000000000;
#define MOD_59 here->method.curr ^= (long) 00020000000000000000000;
#define MOD_60 here->method.curr ^= (long) 00040000000000000000000;
#define MOD_61 here->method.curr ^= (long) 00100000000000000000000;
#define MOD_62 here->method.curr ^= (long) 00200000000000000000000;
#define MOD_63 here->method.curr ^= (long) 00400000000000000000000;
#define MOD_64 here->method.curr ^= (long) 01000000000000000000000;

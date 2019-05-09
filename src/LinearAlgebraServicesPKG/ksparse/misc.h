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

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

#ifndef KSPARSE_MISC_H
#define KSPARSE_MISC_H

#define BSIZE_SP      512

#ifdef HAS_EXIT1
#  define EXIT_NORMAL 1
#  define EXIT_BAD    0
#else
#  define EXIT_NORMAL 0
#  define EXIT_BAD    1
#endif

#ifdef HAS_CTYPE
#  ifndef isalpha
#    include <ctype.h>
#  endif
#endif

#define eq(a,b)  (!strcmp((a), (b)))
#define eqc(a,b)  (cieq((a), (b)))
#define isalphanum(c)   (isalpha(c) || isdigit(c))
#define hexnum(c) ((((c) >= '0') && ((c) <= '9')) ? ((c) - '0') : ((((c) >= \
        'a') && ((c) <= 'f')) ? ((c) - 'a' + 10) : ((((c) >= 'A') && \
        ((c) <= 'F')) ? ((c) - 'A' + 10) : 0)))

#include "strext.h"

extern char *tmalloc();
extern char *trealloc();
extern void txfree();

#define tfree(x)	(txfree(x), x = 0)

#define	alloc(TYPE)	((TYPE *) tmalloc(sizeof(TYPE)))

extern char *copy();
extern char *gettok();
extern void appendc();
extern int scannum();
extern int prefix();
extern int ciprefix();
extern int cieq();
extern void strtolower();
extern int substring();
extern char *tilde_expand( );
extern void cp_printword();

extern char *datestring();
extern char *date_only();
extern char *time_only();
extern double seconds();

extern char *smktemp();

/* Externs from libc */

#ifdef HAS_STDLIB

#  ifndef _STDLIB_INCLUDED
#    define _STDLIB_INCLUDED
#    include <stdlib.h>
#  endif
#  ifndef HAS_BSDRAND
#    define random	rand
#    define srandom	srand
#  endif
#  ifdef HAS_DOSDIRS
#include <unistd.h>
#  endif

#else

#  ifdef HAS_BSDRAND
extern long random();
extern void srandom();
#  else
#    define random	rand
#    define srandom	srand
#  endif

extern void *calloc();
extern void *malloc();
extern void *realloc();
extern char *getenv();
extern int errno;
extern char *getenv();
extern char *getwd();
extern int rand();
extern void srand();
extern int atoi();
extern int kill();
extern int getpid();
extern void qsort();
#  ifdef notdef
extern void exit();
#  endif

#  ifdef HAS_GETCWD
char *getcwd(char *, size_t);
#  endif

#  ifdef HAS_CLEARERR
#    ifndef clearerr
extern void clearerr();
#    endif /* clearerr */
#  endif /* HAS_CLEARERR */

#  ifndef index
#    ifdef HAS_INDEX
extern char *rindex();
extern char *index();
#    else
#      ifdef HAS_STRCHR
/* For some strange reason these lines screw up the compile on linux:
extern char *strchr();
extern char *strrchr();
*/
#      else
#      endif
#    endif
#  endif

#endif	/* else STDLIB */

#ifndef HAS_INDEX
#  ifndef index
#    ifdef HAS_STRCHR
#      define	index	strchr
#      define	rindex	strrchr
#    endif
#  endif
#endif

#ifdef HAS_VPERROR
extern void perror();
#endif

#ifdef HAS_TIME_
#  ifdef HAS_BSDTIME
extern char *timezone();
#  endif
extern char *asctime();
extern struct tm *localtime();
#endif

#ifndef HAS_MEMAVL
#  ifdef HAS_RLIMIT_
extern char *sbrk();
#  endif
#endif

#define false 0
#define true 1

#ifdef HAS_DOSDIRS
typedef	int	*DIR;
struct direct {
	int	d_reclen;
	short	d_ino;
	short	d_namelen;
	char	d_name[20];
	};

#  ifdef __STDC__
extern DIR *opendir(char *);
extern struct direct *readdir(DIR *);
#  else
extern DIR *opendir( );
extern struct direct *readdir( );
#  endif

#endif

#endif /* KSPARSE_MISC_H */

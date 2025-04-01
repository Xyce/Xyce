/*
The original code in this file is Copyright (c) 1985-1991 The Regents of the
University of California and is under the Spice 3f5 BSD Copyright.

All additions and changes are under the following:
//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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

#ifndef KSPARSE_CAPABIL_H
#define KSPARSE_CAPABIL_H

/*
 *  #define-s that are always on
 */

#define CAPZEROBYPASS
#define NEWCONV
/* #define CAPBYPASS	Internal use only */

/*
 *  #define-s to identify common capabilities
 */

#ifdef WANT_X11
#  define HAS_X11
#else
#  ifdef WANT_X10
#    define HAS_X10
#    ifdef WANT_XT
#      define HAS_XT
#    endif
#  endif	/* X10 */
#endif	/* !X11 */

#ifdef WANT_MFB
#  define HAS_MFB
#endif

#ifdef HAS_BSDDIRS
#  define HAS_DIRS_
#endif
#ifdef HAS_SYSVDIRS
#  define HAS_DIRS_
#endif
#ifdef HAS_DOSDIRS
#  define HAS_DIRS_
#endif

#ifdef HAS_BSDTTY
#  define HAS_TTY_
#endif
#ifdef HAS_SYSVTTY
#  define HAS_TTY_
#endif

#ifdef HAS_BSDTIME
#  define HAS_TIME_
#  define HAS_LOCALTIME
#endif
#ifdef HAS_SYSVTIME
#  define HAS_TIME_
#  define HAS_LOCALTIME
#endif

#ifdef HAS_BSDRLIMIT
#  define HAS_RLIMIT_
#endif
#ifdef HAS_SYSVRLIMIT
#  define HAS_RLIMIT_
#endif
#ifdef HAS_MEMAVL
#  define HAS_RLIMIT_
#endif

#ifdef HAS_BSDRUSAGE
#  define HAS_RUSAGE_
#endif
#ifdef HAS_SYSVRUSAGE
#  define HAS_RUSAGE_
#endif

#ifdef HAS_X10
#  define HAS_X_
#endif
#ifdef HAS_X11
#  define HAS_X_
#endif

#ifndef DIR_TERM
#  define DIR_TERM	0
#endif

#ifndef DIR_PATHSEP
#  define DIR_PATHSEP	0
#endif

#ifndef DIR_CWD
#  define DIR_CWD		0
#endif

#ifdef HAS_NO_ERFC
#  ifndef HAS_NO_ERFC_DECL
#    define HAS_NO_ERFC_DECL
#  endif
#endif

#ifdef HAS_NO_IEEE_LOGB
#  ifndef HAS_NO_IEEE_LOGB_DECL
#    define HAS_NO_IEEE_LOGB_DECL
#  endif
#endif

#ifdef HAS_IEEE_SCALBN
#  define scalb scalbn
#endif

#ifdef HAS_NO_IEEE_LOGB_DECL
#  ifdef __STDC__
extern double logb(double), scalb(double, int);
#  else
extern double logb( ), scalb( );
#  endif
#endif

#ifdef HAS_NO_ERFC_DECL
#  ifdef __STDC__
extern double erfc(double);
#  else
extern double erfc( );
#  endif
#endif

#ifndef SIGNAL_TYPE
#  define SIGNAL_TYPE void
#endif

#ifndef SIGNAL_FUNCTION
#  define SIGNAL_FUNCTION SIGNAL_TYPE (*)( )
#endif

#endif /* KSPARSE_CAPABIL_H */

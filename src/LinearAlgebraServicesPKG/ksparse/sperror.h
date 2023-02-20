/*
The original code in this file is Copyright (c) 1985-1991 The Regents of the
University of California and is under the Spice 3f5 BSD Copyright.

All additions and changes are under the following:
//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

#ifndef ERRORS
#define ERRORS

#include "iferrmsg.h"
#include "strext.h"

        /*
         * definitions for error codes returned by SPICE3 routines.
         */

#define E_INTERN E_PANIC
#define E_BADMATRIX (E_PRIVATE+1)/* ill-formed matrix can't be decomposed */
#define E_SINGULAR (E_PRIVATE+2) /* matrix is singular */
#define E_ITERLIM (E_PRIVATE+3)  /* iteration limit reached,operation aborted */
#define E_ORDER (E_PRIVATE+4)    /* integration order not supported */
#define E_METHOD (E_PRIVATE+5)   /* integration method not supported */
#define E_TIMESTEP (E_PRIVATE+6) /* timestep too small */
#define E_XMISSIONLINE (E_PRIVATE+7)    /* transmission line in pz analysis */
#define E_MAGEXCEEDED (E_PRIVATE+8) /* pole-zero magnitude too large */
#define E_SHORT (E_PRIVATE+9)   /* pole-zero input or output shorted */
#define E_INISOUT (E_PRIVATE+10)    /* pole-zero input is output */
#define E_ASKCURRENT (E_PRIVATE+11) /* ac currents cannot be ASKed */
#define E_ASKPOWER (E_PRIVATE+12)   /* ac powers cannot be ASKed */
#define E_NODUNDEF (E_PRIVATE+13) /* node not defined in noise anal */
#define E_NOACINPUT (E_PRIVATE+14) /* no ac input src specified for noise */
#define E_NOF2SRC (E_PRIVATE+15) /* no source at F2 for IM disto analysis */
#define E_NODISTO (E_PRIVATE+16) /* no distortion analysis - NODISTO defined */
#define E_NONOISE (E_PRIVATE+17) /* no noise analysis - NONOISE defined */
#define E_NOSOL   (E_PRIVATE+18) /* no solution of linear system */
#define E_OVERFLOW (E_PRIVATE+19) /* overflow in RHS of matrix -- danger of exception */

char *SPerror();
#endif

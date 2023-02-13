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

/*
 *	Portability, global externs, kitchen sink (yuk).
 *	In future releases this will include things from misc.h and util.h,
 *	   which duplicate each other in places
 */
#ifndef KSPARSE_SPICE_H
#define KSPARSE_SPICE_H

#include <Xyce_config.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifndef	M_PI
#  define M_PI		3.14159265358979323846
#endif
#ifndef	M_E
#  define M_E  	   2.7182818284590452354
#endif
#ifndef	M_LOG2E
#  define M_LOG2E		1.4426950408889634074
#endif
#ifndef	M_LOG10E
#  define M_LOG10E        0.43429448190325182765
#endif
#define TS_LIMIT 5

#include "hw.h"
#include "config.h"
#include "capabil.h"

#define	NUMELEMS(ARRAY)	(sizeof(ARRAY)/sizeof(*ARRAY))

extern char *Spice_Exec_Dir;
extern char *Spice_Lib_Dir;
extern char *Spice_Help_Dir;
extern char *Spice_Model_Dir;
extern char Spice_OptChar;
extern char *Def_Editor;
extern char *Bug_Addr;
extern int AsciiRawFile;
extern char *Spice_Host;
extern char *Spiced_Log;

extern char Spice_Notice[ ];
extern char Spice_Version[ ];
extern char Spice_Build_Date[ ];

extern char *News_File;
extern char *Default_MFB_Cap;
extern char *Spice_Path;
extern char *Help_Path;
extern char *Lib_Path;
extern int  Patch_Level;

#ifdef MAIN_PROGRAM
    int report_interval, new_raw_head, load_mode, mat_dense;
    int device_error, model_error, *timer_flag, dev_math_error;
    double *simulation_time, *logic_break, *timer_calibration;
    double *accepted_simulation_time;
    int exp_num, exp_d;
    void **exp_list, **exp_values;

    char *xfile;
#else
    extern int report_interval, new_raw_head, load_mode, mat_dense;
    extern int device_error, model_error, *timer_flag, dev_math_error;
    extern double *simulation_time, *logic_break, *timer_calibration;
    extern double *accepted_simulation_time;
    extern char *xfile;
    extern int exp_num, exp_d;
    extern void **exp_list, **exp_values;
#endif

#define LOAD_NORMAL 1
#define LOAD_ENERGY 2

#endif

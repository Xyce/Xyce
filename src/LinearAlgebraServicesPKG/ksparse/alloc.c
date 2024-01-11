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

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

/*
 * Memory alloction functions
 */

#include "spice.h"
#include "stdio.h"
#include "misc.h"
#include "util.h"
#ifdef DEBUG_MALLOC
#include "sm.h"
#endif

#include <stdio.h>

#ifndef HAS_BCOPY
#define bzero(s,n) memset(s,0,n)
#endif

/* Malloc num bytes and initialize to zero. Fatal error if the space can't
 * be malloc'd.   Return NULL for a request for 0 bytes.
 */

void bye_bye(i)
    int i;
{
    printf ("inv = %d\n",1/i);
}
/*
*/
char *
tmalloc(num)
    int num;
{
    char *s;

    if (!num)
	return NULL;

#ifdef DEBUG_MALLOC
    s = sm_malloc((unsigned) num);
#else
    s = malloc((unsigned) num);
#endif
    if (!s) {
        fprintf(stderr, 
		"malloc: Internal Error: can't allocate %d bytes.\n", num);
        exit(EXIT_BAD);
    }

    bzero(s, num);

    return(s);
}

char *
trealloc(str, num)
    char *str;
    int num;
{
    char *s;

    if (!num) {
	if (str)

#ifdef DEBUG_MALLOC
		sm_free(str);
#else
		free(str);
#endif

	return NULL;
    }
    if (!str)
	s = tmalloc(num);
    else
#ifdef DEBUG_MALLOC
        s = sm_realloc(str, (unsigned) num);
#else
        s = realloc(str, (unsigned) num);
#endif

    if (!s) {
        fprintf(stderr, 
		"realloc: Internal Error: can't allocate %d bytes.\n", num);
        perror ("realloc");
        s = malloc((unsigned) num);
        bye_bye(0);
        fprintf (stderr, "From malloc of %d bytes: %lx\n",num,s);
        perror ("malloc");
        exit(EXIT_BAD);
    }
    return(s);
}

void
txfree(ptr)
	char	*ptr;
{
	if (ptr)
#ifdef DEBUG_MALLOC
		sm_free(ptr);
#else
		free(ptr);
#endif
}

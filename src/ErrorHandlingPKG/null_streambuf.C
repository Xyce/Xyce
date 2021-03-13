//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/02/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <null_streambuf.h>

/*--------------------------------------------------------------------*/

null_streambuf::null_streambuf() : streambuf()
{
  setp( buf , buf + sizeof(buf) );
}

null_streambuf::~null_streambuf() {}

/*--------------------------------------------------------------------*/
/* Overflow */

int null_streambuf::overflow( int c )
{
  setp( buf , buf + sizeof(buf) );

  return c ;
}

/*--------------------------------------------------------------------*/

int null_streambuf::sync()
{
  return 0 ;
}

streambuf * null_streambuf::setbuf( char * s , streamsize n )
{
  return this ;
}

/*--------------------------------------------------------------------*/

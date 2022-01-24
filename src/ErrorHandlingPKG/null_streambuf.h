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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/05/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef null_streambuf_h
#define null_streambuf_h

#include <iostream>

// Specialize the ANSI Standard C++ streambuf class that throws away everything
// given to it without generating an error.

class null_streambuf : public streambuf {
public:

  // Constructor
  null_streambuf();

  // Destructor
  virtual ~null_streambuf();

protected:

  // Called when output buffer is filled
  virtual int overflow(int c = EOF);

  // Sync is a no-op
  virtual int sync();

  // Setbuf is a no-op
  virtual streambuf * setbuf(char * s , streamsize n);

private:

  null_streambuf(const null_streambuf & ); // Not allowed
  null_streambuf & operator = (const null_streambuf & ); // Not allowed

  char buf[64]; // Throw away buffer
};

/*--------------------------------------------------------------------*/

#endif

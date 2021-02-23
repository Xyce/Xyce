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
// Creation Date  : 05/05/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef mpi_filebuf_h
#define mpi_filebuf_h

#ifdef Xyce_PARALLEL_MPI

#include <iostream>
#include <cstdio>
#include <mpi.h>

//: Specialize the ANSI Standard C++ streambuf class
//: for a parallel file buffer.  The actual file is
//: only touched by the root processor.
//
//  READ MODE: The file is read on the root processor and
//  broadcast one buffer at a time to the remaining processors.
//
//  WRITE MODE: Each processor has a buffer that is locally
//  filled.  When the buffer is full on the root processor the
//  buffer is written to the output file.  When the buffer is
//  full on any other processor the size of the buffer is doubled.
//  The 'mpi_filebuf::flush' method gathers all buffers on the
//  root processor and writes the buffers to the output file.
//
//  GLOBAL: Calls to the 'open', 'flush', 'close', destructor,
//  and 'underflow' methods are global; these calls must be
//  made on all processors.  The 'underflow' method is called
//  by the 'istream' that uses the 'mpi_filebuf' object when
//  ever the input buffer is empty.  Thus reading from an 'istream'
//  that uses an 'mpi_filebuf' must be globally consistent.


class mpi_filebuf : public streambuf {
public:

  //: Construct an MPI-parallel input/output file buffer
  mpi_filebuf();

  //: GLOBAL: Open a file.
  // The file name is only significant on the root processsor.
  // May only be opened as ios::in, ios::out, ios:app.
  mpi_filebuf * open(
    MPI_Comm       communicator ,       /* All processors */
    const int            root_processor ,     /* All processors */
    const std::ios::open_mode file_mode ,          /* All processors */
    const char * const   file_name = NULL );  /* Root processor */

  //: GLOBAL: Close the file.
  // If output mode then flush the output.
  mpi_filebuf * close();

  //: GLOBAL: Flush the buffered output to the file.
  // Sends all buffers to the root processor,
  // write to the file, and flushes the file.
  mpi_filebuf * flush();

  //: GLOBAL: Destructor
  //  Close and then reclaim memory.
  virtual ~mpi_filebuf();

  //: Query if open, a local operations
  int is_open() const ;

  //: When the file buffer is in the 'closed' state set the buffer length,
  //  The input argument must be consistent on all processors; however,
  //  this condition is not checked until the next 'open' operation.
  mpi_filebuf * set_buffer_length( const size_t buffer_length );

  //: Query the current buffer
  void get_buffer( const char * & , size_t & ) const ;

protected:

  //: Called to refill the input buffer
  virtual int underflow();

  //: Called when output buffer is filled
  virtual int overflow( int c = EOF );

  //: Sync is a no-op
  virtual int sync();

  //: Setbuf is a no-op
  virtual streambuf * setbuf( char * s , streamsize n );

private:

  mpi_filebuf( const mpi_filebuf & ); // Not allowed
  mpi_filebuf & operator = ( const mpi_filebuf & ); // Not allowed

  MPI_Comm comm ;            // Communicator
  int      comm_root ;       // Rank of root processor
  FILE *   comm_root_fp ;    // Root processor's file
  int      comm_output ;     // Output file
  char *   comm_buffer ;     // local buffer
  size_t   comm_buffer_len ; // length of buffer
};

/*--------------------------------------------------------------------*/

inline int  mpi_filebuf::is_open() const { return NULL != comm_buffer ; }

/* The SUN has the 'streambuf::pptr()' as a non-const method,
    which violates the ISO/ANSI standard specification.
    Therefore, must cast away the const. */

inline void mpi_filebuf::get_buffer( const char * & b , size_t & n ) const
{ b = comm_buffer ; n = ((mpi_filebuf*)this)->pptr() - comm_buffer ; }

#endif

#endif

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
// Purpose        :This file was stolen from SIERRA to help support the error
//                 reporting in Xyce which was modeled, in part, on that in
//                 SIERRA.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/05/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_UTL_MPIStreambuf_H
#define Xyce_UTL_MPIStreambuf_H

#ifdef Xyce_PARALLEL_MPI

#include <iostream>
#include <stdio.h>
#include <mpi.h>

namespace Xyce {
namespace MPI {

//: Specialize the ANSI Standard C++ std::streambuf class
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
//  The 'streambuf::flush' method gathers all buffers on the
//  root processor and writes the buffers to the output file.
//
//  GLOBAL: Calls to the 'open', 'flush', 'close', destructor,
//  and 'underflow' methods are global; these calls must be
//  made on all processors.  The 'underflow' method is called
//  by the 'istream' that uses the 'streambuf' object when
//  ever the input buffer is empty.  Thus reading from an 'istream'
//  that uses an 'streambuf' must be globally consistent.

class streambuf : public std::streambuf {
public:

  //: Construct an MPI-parallel input/output file buffer
  streambuf();

  //: GLOBAL: Open a file.
  // The file name is only significant on the root processsor.
  // May only be opened as ios::in, ios::out, ios:app.
  streambuf * open(
    MPI_Comm                    communicator ,       /* All processors */
    const int                   root_processor ,     /* All processors */
    const std::ios::open_mode   file_mode ,          /* All processors */
    const char * const          file_name = NULL );  /* Root processor */

  //: GLOBAL: Close the file.
  // If output mode then flush the output.
  streambuf * close();

  //: GLOBAL: Flush the buffered output to the file.
  // Sends all buffers to the root processor,
  // write to the file, and flushes the file.
  streambuf * flush();

  //: GLOBAL: Destructor
  //  Close and then reclaim memory.
  virtual ~streambuf();

  //: Query if open, a local operations
  int is_open() const ;

  //: When the file buffer is in the 'closed' state set the buffer length,
  //  The input argument must be consistent on all processors; however,
  //  this condition is not checked until the next 'open' operation.
  streambuf * set_buffer_length( const size_t buffer_length );

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
  virtual std::streambuf *setbuf( char * s , streamsize n );

private:
  streambuf( const streambuf & ); // Not allowed
  streambuf & operator = ( const streambuf & ); // Not allowed

  MPI_Comm comm ;            // Communicator
  int      comm_root ;       // Rank of root processor
  FILE *   comm_root_fp ;    // Root processor's file
  int      comm_output ;     // Output file
  char *   comm_buffer ;     // local buffer
  size_t   comm_buffer_len ; // length of buffer
};

/*--------------------------------------------------------------------*/

inline int  streambuf::is_open() const { return NULL != comm_buffer ; }

/* The SUN has the 'streambuf::pptr()' as a non-const method,
    which violates the ISO/ANSI standard specification.
    Therefore, must cast away the const. */

inline void streambuf::get_buffer( const char * & b , size_t & n ) const
{ b = comm_buffer ; n = ((streambuf*)this)->pptr() - comm_buffer ; }

#endif

} // namespace MPI
} // namespace Xyce

#endif // Xyce_UTL_MPIStreambuf_H

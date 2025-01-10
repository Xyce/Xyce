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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Creator        : Eric Keiter
//
// Creation Date  : 06/02/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#ifdef Xyce_PARALLEL_MPI
#include <stdio.h>
#include <stdlib.h>
#include <mpi_filebuf.h>

enum { buffer_default_length = 4096 };
enum { buffer_putback_length =   16 };

/*--------------------------------------------------------------------*/

mpi_filebuf::mpi_filebuf()
  : streambuf(),
    comm( MPI_COMM_NULL ),
    comm_root( -1 ),
    comm_root_fp( NULL ),
    comm_output( 0 ),
    comm_buffer( NULL ),
    comm_buffer_len( buffer_default_length )
{}

mpi_filebuf::~mpi_filebuf()
{
  close();
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::set_buffer_length( const size_t len )
{
  // If already open then abort
  if ( NULL != comm_buffer ) return (mpi_filebuf *) NULL ;

   // Wait and verify upon the attempt to open
  comm_buffer_len = buffer_putback_length < len ? len : buffer_putback_length ;

  return this ;
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::open(
  MPI_Comm       communicator ,
  const int            root_processor ,
  const std::ios::open_mode file_mode ,
  const char * const   file_name )
{
  // If already open then abort
  if ( NULL != comm_buffer ) return (mpi_filebuf *) NULL ;

  const int mode =
    ( std::ios::in  == file_mode ) ? 'r' : (
      ( std::ios::out == file_mode ) ? 'w' : (
        ( std::ios::app == file_mode ) ? 'a' : -1 ) );

  int err ;
  int rank ;
  int local, global ;
  int data[3] ;

  // Broadcast the selected root processor and 'C' file mode

  data[0] = root_processor ;
  data[1] = mode ;
  data[2] = comm_buffer_len ;

  if ( MPI_SUCCESS != ( err = MPI_Bcast(data,3,MPI_INT,0,communicator) ) )
    MPI_Abort( communicator , err );

   // Verify that all processors have the same root, mode, and buffer length:

  local = data[0] != root_processor ||
    data[1] != mode ||
    data[2] != comm_buffer_len ;

  if ( MPI_SUCCESS != ( err =
                        MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,communicator) ) )
    MPI_Abort( communicator , err );

  if ( global ) return (mpi_filebuf *) NULL ;

  //--------------------------------------------------------------------
  // Root processor and mode are consistent.
  // All processors try to allocate buffers and the
  // root processor tries to open the file.

  if ( MPI_SUCCESS != ( err =  MPI_Comm_rank( communicator , &rank ) ) )
    MPI_Abort( communicator , err );

  char * const tmp_buf = (char *) malloc( comm_buffer_len );
  FILE *       tmp_fp  = NULL ;

  local = tmp_buf == NULL ; // Failed allocation ?

  if ( root_processor == rank && ! local )
  {
    tmp_fp = fopen( file_name , ( ( ( mode == 'r' ) ? "r" :
                                    ( mode == 'w' ) ? "w" : "a" ) ) );
    local = NULL == tmp_fp ;
  }

  if ( MPI_SUCCESS != ( err =
                        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_BOR,
                                      communicator) ) )
    MPI_Abort( communicator , err );

  if ( global )
  {
    if ( NULL != tmp_buf ) free(   tmp_buf ); // Deallocate
    if ( NULL != tmp_fp  ) fclose( tmp_fp );  // Close the file
    return (mpi_filebuf *) NULL ;
  }

  //--------------------------------------------------------------------
  // All memory allocated and root processor openned the file
  // Update the internal members accordingly.

  comm         = communicator ;
  comm_root    = root_processor ;
  comm_root_fp = tmp_fp ;
  comm_buffer  = tmp_buf ;
  comm_output  = mode != 'r' ;

  // If output then set up put-buffer

  if ( comm_output ) setp( comm_buffer, comm_buffer + comm_buffer_len );

  return this ;
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::close()
{
  mpi_filebuf * tmp = NULL ;

  if ( NULL != comm_buffer )
  {

    flush(); // Flush the buffers

    if ( NULL != comm_root_fp ) fclose( comm_root_fp ); // Close the file

    free( comm_buffer ); // Free the buffer

    if ( comm_output ) setp(NULL,NULL);
    else               setg(NULL,NULL,NULL);

    // Reset the members:

    comm         = MPI_COMM_NULL ;
    comm_root    = -1 ;
    comm_root_fp = NULL ;
    comm_output  = 0 ;
    comm_buffer  = NULL ;

    tmp = this ;
  }

  return tmp ;
}

/*--------------------------------------------------------------------*/
/* Underflow, a global call.
   Read more data from the root processor's file and
   broadcast it to all processors.
 */

int mpi_filebuf::underflow()
{
  if ( NULL != comm_buffer && ! comm_output &&     // Open for read
       ( gptr() == NULL || gptr() >= egptr() ) )  // valid get buffer
  {
    // Length of the buffer, consistent on all processors
    // Entire buffer is offset to accomodate putbacks

    const size_t size = comm_buffer_len - buffer_putback_length ;
    char * const buf  = comm_buffer     + buffer_putback_length ;

    int nread ;
    int err ;

    // Root processor reads from the file and broadcasts the result

    if ( NULL != comm_root_fp ) nread = fread(buf,1,size,comm_root_fp);

    if ( MPI_SUCCESS != ( err =
                          MPI_Bcast( &nread, 1, MPI_INT, comm_root, comm ) ) )
      MPI_Abort(comm,err);

    // If the read is successfull then update the get buffer pointers:

    if ( 0 < nread )
    {

      // Broadcast the read buffer to all processors:

      if ( MPI_SUCCESS != ( err =
                            MPI_Bcast( buf, nread, MPI_BYTE, comm_root, comm )
                            ) )
        MPI_Abort(comm,err);

      // Set the get buffer:

      setg( comm_buffer, buf, buf + nread );

      // Return the next character from the file:
      return *buf ;
    }
  }

  // Failed: set the get buffer to NULL and return EOF
  setg(NULL, NULL, NULL);

  return EOF;
}

/*--------------------------------------------------------------------*/
/* Overflow, a local call.
   Output complete lines of data on the root processor.
   Increase the buffer size on all other processors.
*/

int mpi_filebuf::overflow( int c )
{
  if ( NULL != comm_buffer && comm_output )  // open for write
  {
    // Determine current offset and length:
    char * cur_buffer = comm_buffer ;
    size_t cur_offset = pptr()  - cur_buffer ;
    size_t cur_length = epptr() - cur_buffer ;

    if ( NULL != comm_root_fp )
    {
      if ( fwrite(cur_buffer,1,cur_offset,comm_root_fp) != cur_offset )
      {
        return EOF ; // Write failed
      }
      cur_offset = 0 ;
    }
    else
    {
      // Not root processor, cannot write so increase the buffer size:
      cur_buffer = (char *) realloc( cur_buffer , cur_length *= 2 );
    }

    // If buffer is still good then reset the put-buffer

    if ( NULL != cur_buffer )
    {

      comm_buffer = cur_buffer ;

      setp( cur_buffer + cur_offset, cur_buffer + cur_length );

      if ( c != EOF ) sputc(c);

      return c ;
    }
  }
  return EOF ;
}

/*--------------------------------------------------------------------*/
/* Send output buffers to root processor and
   write them to the output file.
 */

mpi_filebuf * mpi_filebuf::flush()
{
  int result = -1 ; // Failure return value

  if ( NULL != comm_buffer && comm_output )  // Open for write
  {
    int err ;

    result = 0 ;

    // Determine the local length:

    char * cur_buf = comm_buffer ;
    int    cur_len = pptr() - cur_buf ;

    // Determine the global lengths

    char * recv_buf  = NULL ;
    int  * recv_len  = NULL ;
    int  * recv_disp = NULL ;

    int nproc = 0 ;

    if ( NULL != comm_root_fp )
    {

      if ( MPI_SUCCESS != ( err = MPI_Comm_size(comm,&nproc) ) )
        MPI_Abort( comm , err );

      recv_len = (int*) malloc( sizeof(int) * nproc );

      if ( NULL == recv_len ) MPI_Abort( comm , MPI_ERR_UNKNOWN );
    }

    // Gather buffer lengths on the root processor

    if ( MPI_SUCCESS != ( err =
                          MPI_Gather(&cur_len, 1, MPI_INT, recv_len, 1,
                          MPI_INT, comm_root, comm)))
      MPI_Abort( comm , err );

    // Root processor must allocate enough buffer space:

    if ( NULL != comm_root_fp )
    {

      recv_len[ comm_root ] = 0 ; // Don't send to self

      int i ;

      if ( NULL == ( recv_disp = (int*) malloc( sizeof(int) * (nproc + 1) ) ) )
        result = -1 ;

      if ( 0 == result )  // Allocation succeeded
      {
        recv_disp[0] = 0 ;

        for ( i = 0 ; i < nproc ; ++i )
          recv_disp[i+1] = recv_disp[i] + recv_len[i] ;

        if ( 0 < recv_disp[nproc] )
        {
          if ( NULL == ( recv_buf = (char*) malloc( recv_disp[nproc] ) ) )
            result = -1 ;
        }
        else
        {
          result = 1 ; // No need to gather!
        }

        if ( -1 != result )
        {

          // Write the root processor's buffer

          if ( 0 < cur_len )
          {
            if ( fwrite(cur_buf,1,cur_len,comm_root_fp) != cur_len )
              result = -1 ; // Write failed

            cur_len = 0 ; // Wrote this buffer
          }
        }
      }
    }

    // Root process broadcasts that all is well with the allocation

    if ( MPI_SUCCESS != ( err = MPI_Bcast(&result,1,MPI_INT,comm_root,comm)))
      MPI_Abort( comm , err );

    if ( 0 == result ) // All-is-well, need to gather and write
    {
      // Gather the buffers to the root processor

      if ( MPI_SUCCESS != ( err =
                            MPI_Gatherv(cur_buf, cur_len, MPI_BYTE, recv_buf,
                            recv_len, recv_disp, MPI_BYTE, comm_root, comm ) )
           )
        MPI_Abort( comm , err );

      // Output the buffers, beginning with 'comm_root'

      if ( NULL != comm_root_fp )
      {

        int i ;

        for ( i = 1 ; i < nproc && 0 == result ; ++i )
        {
          const int j   = ( i + comm_root ) % nproc ;
          const int len = recv_len[j] ;

          if ( 0 < len )
            if ( fwrite(recv_buf+recv_disp[j],1,len,comm_root_fp) != len )
              result = -1 ; // Write failed
        }
      }

      // Broadcast that the write succeeded

      if ( MPI_SUCCESS != ( err = MPI_Bcast(&result,1,MPI_INT,comm_root,comm)))
        MPI_Abort( comm , err );
    }
    else if ( 1 == result )
    {
      // Did not need to gather

      result = 0 ;
    }

    // Reset the output buffer

    setp( comm_buffer , epptr() );

    // Clean up allocated memory

    if ( NULL != recv_buf  ) free( recv_buf );
    if ( NULL != recv_len  ) free( recv_len );
    if ( NULL != recv_disp ) free( recv_disp );

    // Flush the output file

    if ( NULL != comm_root_fp ) fflush( comm_root_fp );
  }

  return -1 == result ? (mpi_filebuf *) NULL : this ;
                                                       }

/*--------------------------------------------------------------------*/

int mpi_filebuf::sync()
{
  // The root processor will push to file, all others ignore

  if ( NULL != comm_root_fp )
  {

    // Determine the local length:

    char * cur_buf = comm_buffer ;
    int    cur_len = pptr() - cur_buf ;

    if ( 0 < cur_len ) fwrite(cur_buf,1,cur_len,comm_root_fp);

    fflush( comm_root_fp );

    setp( comm_buffer , epptr() );
  }

  return 0 ;
}

streambuf * mpi_filebuf::setbuf( char * s , streamsize n )
{
  return this ;
}

#endif


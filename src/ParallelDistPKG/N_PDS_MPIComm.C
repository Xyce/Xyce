//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose        : Implementation file for the parallel communication
//                  class for Xyce.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/26/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef Xyce_PARALLEL_MPI
 #include <mpi.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
// ----------   Xyce Includes   ----------

#include <N_PDS_MPIComm.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>

using Xyce::DEBUG_PARALLEL;

// ----------  Other Includes   ----------

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::N_PDS_MPIComm( int iargs,
                              char * cargs[] )
  :
  petraCommOwned_(true),
  petraComm_(0)
{
#ifdef Xyce_PARALLEL_MPI
  createMPIComm( iargs, cargs );
  petraComm_ = new Epetra_MpiComm( mpiComm_ );
#else
  petraComm_ = new Epetra_SerialComm();
#endif

  isSerial_ = ( numProc() == 1 );
}

#ifdef Xyce_PARALLEL_MPI
//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/16/06
//-----------------------------------------------------------------------------
N_PDS_MPIComm::N_PDS_MPIComm( MPI_Comm comm )
  :
  petraCommOwned_(true),
  mpiCommOwned_(false)
{
  mpiComm_ = comm;
  petraComm_ = new Epetra_MpiComm( mpiComm_ );

  isSerial_ = ( numProc() == 1 );
}
#endif

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::N_PDS_MPIComm(const N_PDS_MPIComm &right)
  :
  isSerial_(right.isSerial_) ,
  petraCommOwned_(false),
  petraComm_(right.petraComm_)
#ifdef Xyce_PARALLEL_MPI
  ,
  mpiCommOwned_(false),
  mpiComm_(right.mpiComm_)
#endif
{
}

Xyce::Parallel::Machine N_PDS_MPIComm::comm() const 
{
#ifdef Xyce_PARALLEL_MPI
      return mpiComm_;
#else
      return 0;
#endif
}
    
//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::operator=
// Purpose       : Assignment operator.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm & N_PDS_MPIComm::operator=(const N_PDS_MPIComm &right)
{
  if (this != &right)
  {
    isSerial_ = right.isSerial_;
    petraComm_  = right.petraComm_;
    petraCommOwned_ = false;
#ifdef Xyce_PARALLEL_MPI
    mpiComm_ = right.mpiComm_;
    mpiCommOwned_ = false;
#endif
  }

  return *this;
}
//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::~N_PDS_MPIComm
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoeksta, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::~N_PDS_MPIComm()
{
  if( petraCommOwned_ )
    if( petraComm_ ) delete petraComm_;

#ifdef Xyce_PARALLEL_MPI
  if( mpiCommOwned_ )
    MPI_Finalize();
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::numProc
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_MPIComm::numProc() const
{
  int size = 1;
#ifdef Xyce_PARALLEL_MPI
  // Get the machine size.
  if (MPI_SUCCESS != MPI_Comm_size(mpiComm_, &size) )
    Xyce::Report::DevelFatal0().in("N_PDS_MPIComm::numProc")
      << "MPI_Comm_size failed.";
#endif
  return size;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::procID
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_MPIComm::procID() const
{
  int id = 0;
#ifdef Xyce_PARALLEL_MPI
  // Get the machine size.
  if (MPI_SUCCESS != MPI_Comm_rank(mpiComm_, &id) )
    Xyce::Report::DevelFatal0().in("N_PDS_MPIComm::procID")
      << "MPI_Comm_rank failed.";
#endif
  return id;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::scanSum( const double * vals, double * sums,
		const int & count ) const
{
  return ( petraComm_->ScanSum( const_cast<double *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::sumAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::sumAll( const double * vals, double * sums,
		const int & count ) const
{
  return ( petraComm_->SumAll( const_cast<double *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::maxAll( const double * vals, double * maxs,
		const int & count ) const
{
  return ( petraComm_->MaxAll( const_cast<double *> (vals), maxs,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::minAll( const double * vals, double * mins,
		const int & count ) const
{
  return ( petraComm_->MinAll( const_cast<double *> (vals), mins,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::scanSum( const int * vals, int * sums,
		const int & count ) const
{
  return ( petraComm_->ScanSum( const_cast<int *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::sumAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::sumAll( const int * vals, int * sums,
		const int & count ) const
{
  return ( petraComm_->SumAll( const_cast<int *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::maxAll( const int * vals, int * maxs,
		const int & count ) const
{
  return ( petraComm_->MaxAll( const_cast<int *> (vals), maxs,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::minAll( const int * vals, int * mins,
		const int & count ) const
{
  return ( petraComm_->MinAll( const_cast<int *> (vals), mins,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( int * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Bcast( val, const_cast<int &> (count), MPI_INT,
             const_cast<int &> (root), mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( char * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Bcast( val, const_cast<int &> (count), MPI_CHAR,
             const_cast<int &> (root), mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( double * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Bcast( val, const_cast<int &> (count), MPI_DOUBLE,
             const_cast<int &> (root), mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const int * val, const int & count, const int & dest)
									const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Send( const_cast<int *> (val), const_cast<int &> (count), MPI_INT,
            const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : LNGs 
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const long * val, const int & count, const int & dest)
                                                                        const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Send( const_cast<long *> (val), const_cast<int &> (count), MPI_LONG,
            const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const char * val, const int & count, const int & dest)
									const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Send( const_cast<char *> (val), const_cast<int &> (count), MPI_CHAR,
            const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const double * val, const int & count, const int & dest)									const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Send( const_cast<double *> (val), const_cast<int &> (count), MPI_DOUBLE,
            const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( int * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Recv( val, const_cast<int &> (count), MPI_INT, const_cast<int &> (src),
            0, mpiComm_, &status_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( long * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Recv( val, const_cast<int &> (count), MPI_LONG, const_cast<int &> (src),
            0, mpiComm_, &status_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( char * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Recv( val, const_cast<int &> (count), MPI_CHAR, const_cast<int &> (src),
            0, mpiComm_, &status_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( double * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Recv( val, const_cast<int &> (count), MPI_DOUBLE, const_cast<int &> (src),
            0, mpiComm_, &status_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const int * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Rsend( const_cast<int *> (val), const_cast<int &> (count), MPI_INT,
             const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const char * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Rsend( const_cast<char *> (val), const_cast<int &> (count), MPI_CHAR,
             const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const double * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Rsend( const_cast<double *> (val), const_cast<int &> (count), MPI_DOUBLE,
             const_cast<int &> (dest), 0, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::iRecv
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( int * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_INT, const_cast<int &> (src),
             0, mpiComm_, &request_.front() );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::iRecv
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( char * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_CHAR, const_cast<int &> (src),
             0, mpiComm_, &request_.front() );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::iRecv
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( double * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_DOUBLE, const_cast<int &> (src),
             0, mpiComm_, &request_.front() );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::waitAll
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::waitAll()
{
#ifdef Xyce_PARALLEL_MPI
  std::list<MPI_Request>::iterator r=request_.begin();
  std::list<MPI_Request>::iterator r_end = request_.end();

  for ( ; r!=r_end ; ++r)
  {
    MPI_Wait( &(*r), &status_ );
  }
  request_.clear();

  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const int * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  if (DEBUG_PARALLEL)
    MPI_Comm_set_errhandler(mpiComm_, MPI_ERRORS_RETURN);

  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<int *> (val), count, MPI_INT,
                  buf, size, &pos, mpiComm_ );
  if (DEBUG_PARALLEL && err != MPI_SUCCESS)
  {
    std::cerr << "Processor " << procID() << " encountered an error ("<< err << ") calling MPI_Pack on an MPI_INT" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++)
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }

  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const char * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  if (DEBUG_PARALLEL)
    MPI_Comm_set_errhandler(mpiComm_, MPI_ERRORS_RETURN);

  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<char *> (val), count, MPI_CHAR,
                buf, size, &pos, mpiComm_ );
  if (DEBUG_PARALLEL && err != MPI_SUCCESS)
  {
    std::cerr << "Processor " << procID() << " encountered an error (" << err << ") calling MPI_Pack on an MPI_CHAR" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++)
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    for (int i=0; i<count; i++)
      std::cerr << val[i];
    std::cerr << std::endl;
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }

  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const double * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  if (DEBUG_PARALLEL)
    MPI_Comm_set_errhandler(mpiComm_, MPI_ERRORS_RETURN);

  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<double *> (val), count, MPI_DOUBLE,
                buf, size, &pos, mpiComm_ );
  if (DEBUG_PARALLEL && err != MPI_SUCCESS)
  {
    std::cerr << "Processor " << procID() << " encountered an error (" << err << ") calling MPI_Pack on an MPI_DOUBLE" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++)
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }

  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const long * val, const int count, char * buf,
                const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  if (DEBUG_PARALLEL)
    MPI_Comm_set_errhandler(mpiComm_, MPI_ERRORS_RETURN);
  
  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<long *> (val), count, MPI_LONG,
                buf, size, &pos, mpiComm_ );
  if (DEBUG_PARALLEL && err != MPI_SUCCESS)
  { 
    std::cerr << "Processor " << procID() << " encountered an error ("<< err << ") calling MPI_Pack on an MPI_LONG_INT" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;
    
    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++)
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }
  
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos,
		int * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
                count, MPI_INT, mpiComm_ );
  return true;
#else
  return false;
#endif
}


//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos,
		char * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
                count, MPI_CHAR, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos,
		double * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
                count, MPI_DOUBLE, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos,
                long * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
                count, MPI_LONG, mpiComm_ );
  return true;
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::barrier
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
void N_PDS_MPIComm::barrier() const
{
#ifdef Xyce_PARALLEL_MPI
  MPI_Barrier(mpiComm_);
#endif
}

#ifdef Xyce_PARALLEL_MPI
//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::createMPIComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 05/01/2018
//-----------------------------------------------------------------------------
void N_PDS_MPIComm::createMPIComm( int iargs, char * cargs[] )
{
  // Set the MPI Communicator
  mpiCommOwned_ = true;  
  mpiComm_ = MPI_COMM_WORLD;

#ifdef HAVE_UNISTD_H
  //
  // In some implementations of MPICH using the ch_p4 communicator,
  // calling mpirun -np 1 Xyce , or just running Xyce (which implies -np 1)
  // causes MPI_Init() to reset the current working directory to the location
  // of the Xyce binary.  Here we record the current working directory
  // and reset it if it has changed by calling MPI_Init().  This requires
  // unistd.h so I'm enclosing this test in ifdefs so it doesn't break
  // other builds unnecessarily
  //
  //                
  const int maxPathLen = 4096;
  char originalCWD[maxPathLen];
  char * getcwdRV1 = getcwd( originalCWD, maxPathLen );

  // Check if MPI has already been initialized.  If it has, we don't own the MPI_Comm.
  int initialized = false;
  MPI_Initialized(&initialized);
 
  if ( !initialized )
  {
    if ( MPI_SUCCESS != MPI_Init( &iargs, &cargs ) )
      Xyce::Report::DevelFatal().in("N_PDS_MPIComm::N_PDS_MPIComm")
        << "N_PDS_MPIComm::initMPI - MPI_Init failed.";
  }
  else
  {
    mpiCommOwned_ = false;
  }

  char currentCWD[maxPathLen];
  char * getcwdRV2 = getcwd( currentCWD, maxPathLen );

  if( (getcwdRV1 != NULL) && (getcwdRV2 != NULL) )
  {
    if( strcmp( originalCWD, currentCWD ) != 0 )
    {
      if( chdir( originalCWD ) != 0 )
      {
        // changing the directory back may have failed.  Issue a warning
        // and try and continue.
        Xyce::Report::UserWarning()
          << "Resetting working directory failed.  Trying to continue.";
      }
    }
  }
  else
  {
    Xyce::Report::UserWarning()
      << "Could not get working directory.  Trying to continue.";
  }
#endif
 
  if (Xyce::DAKOTA)
  {
    // if -ppp <number> was specified, then split the communicator
    int procPerProblem = 0;
    for( int i=0; i<iargs; i++ )
    {
      std::string anArg( cargs[i] );
      if( anArg.compare( "-ppp" ) == 0 )
      {
        // next arg should be number of processors per problem
        i++;
        int nextI = i;
        if( nextI < iargs )
        {
          std::string argValue( cargs[ nextI ] );
          std::stringstream iost;
          iost << argValue;
          iost >> procPerProblem;
        }
        else
        {
          // issue warning that "-ppp " was not followed by a number
          Xyce::Report::UserWarning()
            << "Could not split communicator, -ppp option must be followed by a number.";

        }
        break;
      }
    }

    if( procPerProblem > 0 )
    {
      // split mpi communicator
      int processID = procID();
      MPI_Comm myNewMpiComm;
      Xyce::dout() << "N_PDS_MPIComm on procID, " << processID << ", procPerProblem = " << procPerProblem << std::endl; 
        const std::string errorMsgForSplit( "N_PDS_MPIComm::initMPI - MPI_Comm_split failed.");
        if( MPI_SUCCESS != MPI_Comm_split( MPI_COMM_WORLD, processID, procPerProblem, &myNewMpiComm ) )
          Xyce::Report::DevelFatal() << errorMsgForSplit;
        mpiComm_ = myNewMpiComm;
    }
  }
}
#endif

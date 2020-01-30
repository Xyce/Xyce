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

//-----------------------------------------------------------------------------
//
// Purpose        : Specification file for the abstract parallel communication
//                  class for Xyce.  This class will contain parallel data and
//                  functions.
//
// Special Notes  : It assumes that all parallel communication is being handled
//                  through MPI.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Comm_h
#define Xyce_N_PDS_Comm_h

#include <list>

#include <N_PDS_ParallelMachine.h>


class Epetra_Comm;

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : Communicator
// Purpose       : Abstract parallel communication class for Xyce.  This class
//                 will contain parallel data and functions.
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
class Communicator
{

public:
  Communicator()
  {}

  virtual ~Communicator()
  {}

private:
  Communicator(const Communicator & right);
  Communicator & operator=(const Communicator & right);

  bool operator==(const Communicator & right) const;
  bool operator!=(const Communicator & right) const;

public:
  virtual Communicator *clone() const = 0;

  // Get my processor ID.
  virtual int procID() const = 0;

  // Get the total number of processors in this communicator.
  virtual int numProc() const = 0;

  // Get the serial flag (true => serial, not parallel operation).
  virtual bool isSerial() const = 0;

  // Get the last processor flag (by default true in serial)
  // NOTE:  This is used for generating augmented linear systems.
  virtual bool isLastProc() const = 0;

  // Wrappers for Petra_Comm functionality
  virtual bool scanSum(const double * vals, double * sums, const int & count) const = 0;
  virtual bool sumAll(const double * vals, double * sums, const int & count) const = 0;
  virtual bool maxAll(const double * vals, double * maxs, const int & count) const = 0;
  virtual bool minAll(const double * vals, double * mins, const int & count) const = 0;

  virtual bool scanSum(const int * vals, int * sums, const int & count) const = 0;
  virtual bool sumAll(const int * vals, int * sums, const int & count) const = 0;
  virtual bool maxAll(const int * vals, int * maxs, const int & count) const = 0;
  virtual bool minAll(const int * vals, int * mins, const int & count) const = 0;

  // Wrapper for MPI broadcasts
  virtual bool bcast(int * val, const int & count, const int & root) const = 0;
  virtual bool bcast(char * val, const int & count, const int & root) const = 0;
  virtual bool bcast(double * val, const int & count, const int & root) const = 0;

  // Wrappers for MPI Sends
  virtual bool send(const int * val, const int & count, const int & dest) const = 0;
  virtual bool send(const char * val, const int & count, const int & dest) const = 0;
  virtual bool send(const double * val, const int & count, const int & dest) const = 0;
  virtual bool send(const long * val, const int & count, const int & dest) const = 0;

  // Wrappers for MPI Recvs
  virtual bool recv(int * val, const int & count, const int & src) const = 0;
  virtual bool recv(char * val, const int & count, const int & src) const = 0;
  virtual bool recv(double * val, const int & count, const int & src) const = 0;
  virtual bool recv(long * val, const int & count, const int & src) const = 0;

  // Wrappers for MPI RSends
  virtual bool rSend(const int * val, const int & count, const int & dest) const = 0;
  virtual bool rSend(const char * val, const int & count, const int & dest) const = 0;
  virtual bool rSend(const double * val, const int & count, const int & dest) const = 0;

  // Wrappers for MPI IRecvs
  virtual bool iRecv(int * val, const int & count, const int & src) const = 0;
  virtual bool iRecv(char * val, const int & count, const int & src) const = 0;
  virtual bool iRecv(double * val, const int & count, const int & src) const = 0;

  // Wrapper for MPI Waitall
  virtual bool waitAll() = 0;

  // Wrappers for MPI pack and unpack functionality
  virtual bool pack(const int * val, const int count, char * buf,
                    const int size, int & pos) const = 0;
  virtual bool pack(const char * val, const int count, char * buf,
                    const int size, int & pos) const = 0;
  virtual bool pack(const double * val, const int count, char * buf,
                    const int size, int & pos) const = 0;
  virtual bool pack(const long * val, const int count, char * buf,
                    const int size, int & pos) const = 0;

  virtual bool unpack(const char * buf, const int size, int & pos, int * val,
                      const int count) const = 0;
  virtual bool unpack(const char * buf, const int size, int & pos, char * val,
                      const int count) const = 0;
  virtual bool unpack(const char * buf, const int size, int & pos,
                      double * val, const int count) const = 0;
  virtual bool unpack(const char * buf, const int size, int & pos, long * val,
                      const int count) const = 0; 
 
  virtual Epetra_Comm * petraComm() = 0;

  // Communicator Barrier function.
  // A no-op for a serial communicator.  For MPI, it causes each processor in
  // the communicator to wait until all processors have arrived.
  virtual void barrier() const = 0;

  virtual Parallel::Machine comm() const = 0;
};

// Return a new Comm
N_PDS_Comm * createPDSComm(int iargs = 0, char * cargs[] = 0, Xyce::Parallel::Machine comm = MPI_COMM_NULL );

} // namespace Parallel
} // namespace Xyce

typedef Xyce::Parallel::Communicator N_PDS_Comm;

#endif // Xyce_N_PDS_Comm_h

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
//                  
//
// Special Notes  : 
//                  
//
// Creator        : David Baur
//
// Creation Date  : 
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_PDS_Serial_h
#define Xyce_PDS_Serial_h

#include <stdexcept>
#include <string>
#include <vector>

#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_PDS_ParallelMachine.h>

#ifndef Xyce_PARALLEL_MPI

namespace {

static const int MPI_MAX = 0;
static const int MPI_MIN = 0;
static const int MPI_SUM = 0;
static const int MPI_PROD = 0;
static const int MPI_LAND = 0;
static const int MPI_BAND = 0;
static const int MPI_LOR = 0;
static const int MPI_BOR = 0;
static const int MPI_LXOR = 0;
static const int MPI_BXOR = 0;

}

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Function      : op_identifier_compare_op
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:01 2015
//-----------------------------------------------------------------------------
inline int op_identifier_compare_op() {
  return 0;
}

struct ReduceInterface
{};

struct Sum
{
  template <typename T>
  inline Sum(T * dest, const T *source) {
    *dest += *source;
  }
};

template<class Op, class LocalIt, class GlobalIt = LocalIt>
struct Reduce : public ReduceInterface
{
  typedef typename std::iterator_traits<LocalIt>::value_type value_type;
  typedef typename std::iterator_traits<LocalIt>::difference_type difference_type;

  Reduce(LocalIt local_begin, LocalIt local_end, GlobalIt global_begin, GlobalIt global_end)
  {}
  
};

//-----------------------------------------------------------------------------
// Function      : AllReduce
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:06 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
AllReduce(Machine comm, int op, T *src_dest, size_t size)
{}

//-----------------------------------------------------------------------------
// Function      : AllReduce
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:10 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
AllReduce(Machine comm, int op, const T *source, T *dest, size_t size)
{
  std::copy(source, source + size, dest);
}

//-----------------------------------------------------------------------------
// Function      : AllReduce
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:16 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
AllReduce(Machine comm, int op, std::vector<T> &src_dest)
{}

//-----------------------------------------------------------------------------
// Function      : AllGather
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:20 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
AllGather(Machine mpi_comm, const std::vector<T> &source, std::vector<T> &dest)
{
  if (source.size() != dest.size())
    throw std::runtime_error("Xyce::Serial::AllGather(MPI_Comm mpi_comm, std::vector<T> &source, std::vector<T> &dest) vector lengths not equal");

  dest = source;
}

//-----------------------------------------------------------------------------
// Function      : AllGather
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:24 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
AllGather(Machine mpi_comm, const T &source, std::vector<T> &dest)
{
  if (dest.size() != 1)
    throw std::runtime_error("Xyce::Serial::AllGather(MPI_Comm mpi_comm, const T &source, std::vector<T> &dest) vector lengths not equal");

  dest[0] = source;
}

//-----------------------------------------------------------------------------
// Function      : Broadcast
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:28 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
Broadcast(Machine comm, T *src_dest, size_t len, int root)
{}

//-----------------------------------------------------------------------------
// Function      : Broadcast
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:31 2015
//-----------------------------------------------------------------------------
inline void Broadcast(Machine mpi_comm, std::string &s, int root)
{}

//-----------------------------------------------------------------------------
// Function      : Broadcast
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:34 2015
//-----------------------------------------------------------------------------
inline void Broadcast(Machine mpi_comm, Util::Marshal &m, int root)
{}

//-----------------------------------------------------------------------------
// Function      : AllWriteString
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:40 2015
//-----------------------------------------------------------------------------
inline void
AllWriteString(
  Machine               comm,
  std::ostream &        os,
  const std::string &   message)
{
  os << message;
}

//-----------------------------------------------------------------------------
// Function      : Barrier
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:45 2015
//-----------------------------------------------------------------------------
inline void Barrier(Machine comm)
{}

//-----------------------------------------------------------------------------
// Function      : GatherV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:48 2015
//-----------------------------------------------------------------------------
inline void
GatherV(
  Machine                       comm,
  unsigned                      root,
  const std::string &           src,
  std::vector<std::string> &    dest)
{
  dest.resize(1);
  dest[0] = src;
}

//-----------------------------------------------------------------------------
// Function      : GatherV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:53 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
GatherV(
  Machine                       mpi_comm,
  unsigned                      root,
  const std::vector<T> &        src,
  std::vector<T> &              dest)
{
  dest = src;
}

//-----------------------------------------------------------------------------
// Function      : AllGatherV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:42:56 2015
//-----------------------------------------------------------------------------
inline void
AllGatherV(
  Machine                       comm,
  const std::string &           src,
  std::vector<std::string> &    dest)
{
  dest.resize(1);
  dest[0] = src;
}

//-----------------------------------------------------------------------------
// Function      : AllGatherV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Apr  6 10:43:00 2015
//-----------------------------------------------------------------------------
template<class T>
inline void
AllGatherV(
  Machine                       mpi_comm,
  const std::vector<T> &        src,
  std::vector<T> &              dest)
{
  dest = src;
}

} // namespace Parallel
} // namespace Xyce

#endif // Xyce_PARALLEL_MPI

#endif // Xyce_PDS_Serial_h

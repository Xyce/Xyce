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
// Purpose        :
//
// Special Notes  :
//
// Creator        : David G. Baur
//
// Creation Date  : 05/12/15
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_UTL_Pack_h
#define Xyce_N_UTL_Pack_h

#include <N_PDS_fwd.h>


namespace Xyce {

template <class T>
class Pack
{
public:
  // Counts bytes needed to pack block.
  static int packedByteCount(const T &t);

  // Packs OptionBlock into char buffer using MPI_PACK.
  static void pack(const T &t, char * buf, int bsize, int & pos, Parallel::Communicator * comm);

  // Unpacks OptionBlock from char buffer using MPI_UNPACK.
  static void unpack(T &t, char * pB, int bsize, int & pos, Parallel::Communicator * comm);
};

// Counts bytes needed to pack block.
template <class T>
inline int packedByteCount(const T &t) {
  return Pack<T>::packedByteCount(t);
}


// Packs OptionBlock into char buffer using MPI_PACK.
template <class T>
void pack(const T &t, char * buf, int bsize, int & pos, Parallel::Communicator * comm) {
  Pack<T>::pack(t, buf, bsize, pos, comm);
}

// Unpacks OptionBlock from char buffer using MPI_UNPACK.
template <class T>
void unpack(T &t, char * pB, int bsize, int & pos, Parallel::Communicator * comm) {
  Pack<T>::unpack(t, pB, bsize, pos, comm);
}

} // namespace Xyce

#endif // Xyce_N_UTL_Pack_h

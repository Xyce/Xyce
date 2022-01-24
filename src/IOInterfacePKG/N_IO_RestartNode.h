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

//-----------------------------------------------------------------------------
//
// Purpose        : Node storing restart info associated with an ID
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 8/22/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_RestartNode_h
#define Xyce_N_IO_RestartNode_h

#include <vector>
#include <iosfwd>

#include <N_DEV_fwd.h>
#include <N_TOP_fwd.h>

#include <N_UTL_Pack.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : RestartNode
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/22/01
//-----------------------------------------------------------------------------
class RestartNode
{
public:

  // Constructor
  RestartNode(const std::string & id = "", const int inType = _VNODE ) : ID(id), type(inType), devState(0) { }

  // Destructor
  ~RestartNode();

  // Copy constructor
  RestartNode(const RestartNode & right);

  // Assignment operator
  RestartNode & operator = (const RestartNode & right);

  void dump( std::ostream & os ) const;
  void restore( std::istream & is );

  std::string ID;
  int type;

  std::vector< std::vector< double > > solnVarData;
  std::vector< std::vector< double > > stateVarData;
  std::vector< std::vector< double > > storeVarData;

  Device::DeviceState * devState;

  friend std::ostream & operator << (std::ostream & os, const RestartNode & rn);
};

} // namespace IO
} // namespace Xyce

#endif

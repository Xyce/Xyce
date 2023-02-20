//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_TOP_CktNode_V.h>
#include <N_TOP_Topology.h>
#include <N_TOP_Indexor.h>
#include <N_PDS_Manager.h>

namespace Xyce {
namespace Topo {

void
CktNode_V::loadNodeSymbols(
  Topology &            topology) const
{
  Util::SymbolTable &node_symbols = topology.getNodeSymbols();

  Indexor indexor(topology.getPDSManager());

  int global_id = get_SolnVarGIDList().front();
  if (global_id >= 0)
  {
    std::vector<int> translate_vector(1, global_id);
    indexor.globalToLocal(Parallel::SOLUTION_OVERLAP_GND, translate_vector);

    node_symbols[Util::SOLUTION_SYMBOL][get_id()] = translate_vector[0];
    node_symbols[Util::EXTERN_SYMBOL][get_id()] = translate_vector[0];
  }
}


void
CktNode_V::varTypeList(
  std::vector<char> &   varTypeVec ) const
{
  varTypeVec.push_back('V');
}

//-----------------------------------------------------------------------------
// Function      : CktNode_V::put
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
std::ostream& CktNode_V::put(std::ostream& os) const
{
  os << "CN_V: " << get_id() << std::endl;
  os << "   GID: " << gID_ << "  Proc: " << procNum_ << std::endl;
  os << "   Owned: " << isOwned_ << std::endl;
  os << "   Soln Var GID List: ";
  for (std::vector<int>::const_iterator it_iL = solnVarGIDList_.begin(), it_iL_end = solnVarGIDList_.end(); it_iL != it_iL_end; ++it_iL )
    os << *it_iL << "  ";

  return os << std::endl;

}

} // namespace Topo
} // namespace Xyce

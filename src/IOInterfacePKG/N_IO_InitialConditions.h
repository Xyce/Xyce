//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_InitialConditions_h
#define Xyce_N_IO_InitialConditions_h

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace IO {

struct InitialConditionsData
{
  typedef std::map<std::string, std::pair<int, double>, LessNoCase > NodeNamePairMap;

  enum {IC_TYPE_UNDEFINED, IC_TYPE_IC, IC_TYPE_NODESET};

  InitialConditionsData()
    : icType_(IC_TYPE_UNDEFINED),
      ICflag_(false),
      nodesetflag_(false),
      saveFlag_(false),
      saveFileLevel_("ALL"),
      saveFileType_(".NODESET"),
      saveOutputFile_(""),
      saveTime_(0.0),
      total_soln_(0),
      op_found_(0)
  {}

  int icType_;
  bool ICflag_;
  bool nodesetflag_;
  bool saveFlag_;
  std::string saveFileLevel_;
  std::string saveFileType_;
  std::string saveOutputFile_;
  double saveTime_;
  int total_soln_;
  int op_found_;

  NodeNamePairMap opData_;
};

class InitialConditionsManager 
{
public:
  InitialConditionsManager(
    const std::string &         netlist_filename);

  InitialConditionsData::NodeNamePairMap &getICData( int &found, int &ic_type);

  double getSaveTime() const
  {
    return icData_.saveTime_;
  }

  bool getSaveFlag() const
  {
    return icData_.saveFlag_;
  }
  
  bool registerIC(const Util::OptionBlock & option_block);
  bool registerNodeSet(const Util::OptionBlock & option_block);
  bool registerSave(const Util::OptionBlock & option_block);

  bool setupInitialConditions(
    Parallel::Machine           comm,
    const NodeNameMap &         allNodes_,
    const AliasNodeMap &        alias_nodes,
    Linear::Vector &            solnVec,
    Linear::Vector &            flagVec);

  void outputDCOP(
    Parallel::Machine           comm,
    const NodeNameMap &         all_nodes, 
    const Linear::Vector &      solution);

private:
  const std::string             netlistFilename_;
  bool                          outputOnceAlreadyFlag_;

  InitialConditionsData         icData_;

  std::vector<Util::OptionBlock> ICblockVec_;
  std::vector<Util::OptionBlock> nodesetblockVec_;
};


bool registerPkgOptionsMgr(InitialConditionsManager & output_manager, PkgOptionsMgr &options_manager);

bool setupIC_or_NODESET(Parallel::Machine comm, const NodeNameMap &all_nodes, const AliasNodeMap &alias_nodes, Linear::Vector & solnVec, Linear::Vector & flagVec, int icType, std::vector<Util::OptionBlock> & initBlockVec, InitialConditionsData::NodeNamePairMap &opData_, int &op_found_, int &total_soln_);

void outputIC_or_NODESET(Parallel::Machine comm, std::ofstream &os, const std::string &saveFileType_, InitialConditionsData::NodeNamePairMap all_nodes, const Linear::Vector & solution);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_InitialConditions_h

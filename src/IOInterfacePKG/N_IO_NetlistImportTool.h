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
// Filename      :  NetlistImportTool.h
//
// Purpose       : Declare the interface to read and parse a netlist for
//                 an electrical circuit.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/00
//
//-------------------------------------------------------------------------

#ifndef NetlistImportTool_H
#define NetlistImportTool_H

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>

#include <N_IO_OutputTypes.h>

#include <N_IO_CircuitBlock.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CircuitMetadata.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_ParsingMgr.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Stats.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace IO {

class UndefinedName
{
public:
  UndefinedName()
    : name_(),
      netlistLocation_()
  {}

  UndefinedName(const std::string &name, const NetlistLocation &netlist_location)
    : name_(name),
      netlistLocation_(netlist_location)
  {}

  UndefinedName(const UndefinedName &s0)
    : name_(s0.name_),
      netlistLocation_(s0.netlistLocation_)
  {}

  UndefinedName &operator=(const UndefinedName &s0)
    {
      name_ = s0.name_;
      netlistLocation_ = s0.netlistLocation_;
      return *this;
    }

public:
  const std::string &getName() const
  {
    return name_;
  }

  void setName(const std::string &name)
  {
    name_ = name;
  }

  const NetlistLocation &getNetlistLocation() const
  {
    return netlistLocation_;
  }

  void setNetlistLocation(const NetlistLocation &netlist_location)
  {
    netlistLocation_ = netlist_location;
  }
  
private:
  std::string         name_;
  NetlistLocation     netlistLocation_;
};

inline bool operator<(const UndefinedName &s0, const UndefinedName &s1)
{
  return s0.getName() < s1.getName();
}

typedef std::set<UndefinedName> UndefinedNameSet;

//-----------------------------------------------------------------------------
// Class         : NetlistImportTool
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 07/28/00 ?
//-----------------------------------------------------------------------------
class NetlistImportTool
{
public:
  
  NetlistImportTool(Util::Op::BuilderManager &op_builder_manager,
                    const ParsingMgr &parsing_manager);

  ~NetlistImportTool();

  // R Result
  // R- The result of constructing a circuit from a netlist.
  // I netlistFile
  // I- The file containing the netlist describing the circuit to be
  // I- constructed.
  // This function performs three basic steps. First, it reads and parses a
  // netlist from a specified file and stores the parsed fields in a structure
  // as strings. Second, it interprets the parsed data and stores the data in
  // an appropriate structure. Finally, it builds the Xyce Circuit.
  int constructCircuitFromNetlist(
    const CmdParse &                                            command_line,
    HangingResistor &                                           hanging_resistor,
    const std::string &                                         netlistFilename,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
    Topo::Topology &                                            topology,
    Parallel::Communicator &                                    pds_comm,
    PkgOptionsMgr &                                             options_manager,
    OutputMgr &                                                 output_manager,
    Device::DeviceMgr &                                         device_manager,
    Measure::Manager &                                          measure_manager,
    FourierMgr &                                                fourier_manager);

  static void populateMetadata(IO::PkgOptionsMgr &   options_manager);

  bool registerDCOptions(Util::ParamList::const_iterator it, Util::ParamList::const_iterator end);
  bool registerSTEPOptions(Util::ParamList::const_iterator it, Util::ParamList::const_iterator end);
  bool setDISTOptions( const Util::OptionBlock& distOptions );
  bool registerSAMPLINGOptions(Util::ParamList::const_iterator it, Util::ParamList::const_iterator end);
  bool registerEMBEDDEDSAMPLINGOptions(Util::ParamList::const_iterator it, Util::ParamList::const_iterator end);
  bool registerPCEOptions(Util::ParamList::const_iterator it, Util::ParamList::const_iterator end);

  bool setParserOptions (const Util::OptionBlock & OB);

  const IO::AliasNodeMap &getAliasNodeMap() const
  {
    return aliasNodeMap_;
  }

  const UndefinedNameSet &getDeferredParameterCheck() const
  {
    return deferredUndefinedParameters_;
  }
  
  const Util::ParamMap &getMainContextFunctions() const
  {
    return circuitContext_.getFunctions();
  }

  const Util::UParamList &getMainContextParams() const
  {
    return circuitContext_.getParams();
  }

  const Util::ParamList &getMainContextGlobalParams() const
  {
    return circuitContext_.getGlobals();
  }

  bool getUseMOR() const
  {
    return useMOR_;
  }

  std::string getAnalysisName() const
  {
    return mainCircuitBlock_->getAnalysisName();
  }

public:

private:
  const ParsingMgr &                    parsing_manager;
  unordered_set<std::string>            nodeNames_;
  std::vector<std::string>              stepParams_;
  std::vector<std::string>              samplingParams_;
  std::vector<std::string>              embeddedSamplingParams_;
  std::vector<std::string>              pceParams_;
  std::vector<std::string>              dcParams_;
  Util::OptionBlock                     distOptions_;
  Util::OptionBlock                     parserOptions_;
  UndefinedNameSet                      deferredUndefinedParameters_;

  // Maps for parsing information and streams.
  std::map<std::string, FileSSFPair>    ssfMap_;
  std::map<std::string, IncludeFileInfo> iflMap_;


  CircuitBlock *                        mainCircuitBlock_;
  DistributionTool *                    distributionTool_;
  AliasNodeMap                          aliasNodeMap_;

  std::list<CircuitContext *>           contextList_;
  CircuitContext *                      currentContextPtr_;

  CircuitMetadata                       metadata_;
  CircuitContext                        circuitContext_;
  unordered_set<std::string>            modelNames_;
  bool                                  useMOR_;
};

bool registerPkgOptionsMgr(NetlistImportTool &netlist_import_tool, PkgOptionsMgr &options_manager);

// check for name collisions between nodes and devices
void checkNodeDevConflicts(const unordered_set<std::string> &device_names, N_PDS_Comm &pdsComm);

void getLeadCurrentDevices(const Util::ParamList &variable_list, std::set<std::string> &devicesNeedingLeadCurrents);

void getWildCardLeadCurrentDevices(
  const Util::ParamList &variable_list,
  const unordered_set<std::string> & device_names,
  std::set<std::string> &devicesNeedingLeadCurrents);

void printLineDiagnostics(
  Parallel::Machine                             comm,
  const OutputParameterMap &                    output_parameter_map,
  const Device::DeviceMgr &                     device_manager,
  const Measure::Manager &                      measure_manager,
  const unordered_set<std::string> &            node_names,
  const std::vector<std::string> &              step_params,
  const std::vector<std::string> &              sampling_params,
  const std::vector<std::string> &              embeddedSampling_params,
  const std::vector<std::string> &              pce_params,
  const std::vector<std::string> &              dc_params,
  const unordered_set<std::string> &            device_names,
  const IO::AliasNodeMap &                      alias_node_map,
  UndefinedNameSet &                            deferred_parameter_check,
  bool &                                        iStarRequested);

void deferredParameterDiagnostics(
  Parallel::Machine                             comm,
  const UndefinedNameSet &                      deferred_parameter_check,
  const Device::DeviceMgr &                     device_manager);

bool deferErrorCheckUntilOpCreation(std::string devStr);

template <class It>
void registerGlobalParams(Device::DeviceMgr &device_manager, It begin, It end)
{
  for (; begin != end; ++begin)
  {
    device_manager.addGlobalPar(*begin);
  }
}

bool registerCircuitOptions(PkgOptionsMgr &options_manager, std::list<Util::OptionBlock> &option_block_list);

bool registerDistOptions(PkgOptionsMgr &options_manager, std::list<Util::OptionBlock> &option_block_list);

} // namespace IO
} // namespace Xyce

#endif // NetlistImportTool_H

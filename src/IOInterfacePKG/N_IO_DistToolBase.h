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

//-----------------------------------------------------------------------------
//
// Purpose        : Declares the DistToolBase class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
//
// Special Notes  :
//
// Creator        : Eric Rankin, SNL
//
// Creation Date  : 03/12/2003
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_IO_DistToolBase_h
#define Xyce_N_IO_DistToolBase_h

#include <vector>
#include <string>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_ParsingMgr.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class          : DistToolBase
// Purpose        : Provide basic common functionality for all distribution tools.
//-----------------------------------------------------------------------------
class DistToolBase : public DistributionTool
{
public:
  DistToolBase(
    Parallel::Communicator *                 pdsCommPtr,
    CircuitBlock &                           circuit_block,
    std::map<std::string,FileSSFPair>      & ssfMap,
    const ParsingMgr &                       parsing_manager
    );

  virtual ~DistToolBase(); 

  // send options, metatdata, and context to all procs
  virtual bool broadcastGlobalData();

  std::list<Util::OptionBlock>& getAdditionalOptions() { return addOptions_; }

  Util::UParamList & getAdditionalGlobalParams() { return addResolvedGlobalParams_; }

protected:

  // Parse the given include file for 2nd pass
  virtual bool parseIncludeFile(std::string const& includeFiles,
                                const std::string &libSelect) = 0;

  // helper function for parseIncludeFile()
  virtual void restorePrevssfInfo(
    SpiceSeparatedFieldTool* oldssfPtr,
    const std::string& old_netlistFilename,
    int oldFilePos,
    int oldLineNumber) = 0;

  // Expand a subcircuit instance by adding the devices and
  // device models that compose the subcircuit to the main
  // (top level) circuit. Prepend device names and nodes with
  // subcircuitPrefix.
  virtual bool expandSubcircuitInstance(DeviceBlock & subcircuitInstance,
                                        const std::string &libSelect, 
                                        std::vector<std::string> &libInside) = 0;

  // get the next line from the SpiceSeparatedField pointer
  bool getLine(TokenVector& line,
               const std::string &libSelect, 
               std::vector<std::string> &libInside);

  // Process a device line on processor zero, or serial.
  bool handleDeviceLine(TokenVector const& deviceLine,
                        const std::string &libSelect, 
                        std::vector<std::string> &libInside);

  // Post process a mutual inductor in the current circuit.
  bool handleMutualInductance( DeviceBlock & device );

  // Post process an IBIS Buffer device in the current circuit.
  bool handleIBISdevice( DeviceBlock & device );

  // Fully parse and instantiate a single device.
  bool instantiateDevice( DeviceBlock & device,
      const std::string & prefix, const unordered_map<std::string,std::string>& nodeMap,
      const std::string &libSelect, std::vector<std::string> &libInside);

  // Find any .IC or .NODESET statements associated with the current 
  // subcircuit.  This will be called in expandSubcircuitInstance.
  void find_IC_NODESET_OptionBlock(const std::string& modelName,
                                   const std::string& subcircuitPrefix,
                                   const std::vector<std::string>& subcircuitNodes,
                                   const std::vector<std::string>& subcircuitInstanceNodes);

  // Return true if the option block needs to be checked during each subcircuit expansion.
  bool check_IC_NODESET_OptionBlock();

  // send circuit context to all procs
  void setCircuitContext();

  // stage options for tx
  void setCircuitOptions();

  // Parallel distribution data
  Parallel::Communicator *      pdsCommPtr_;            ///< parallel services
  int                           numProcs_;              ///< total number of available procs for distribution

  bool                          circuitContextReady_;   ///< flags for managing global data
  bool                          circuitOptionsReady_;   ///< flags for managing global data

  bool                          checkSubcktICNODESET_;
  bool                          resolveSubcktICNODESET_; 

  int                           charBufferSize_;        ///< length of buffer
  int                           charBufferPos_;         ///< length of data packed in buffer
  char *                        charBuffer_;

  // This list will be used to collect .IC and .NODESET statements from subcircuits
  // during device distribution.
  std::list<Util::OptionBlock>  addOptions_;

  // this container will contain global parameters from subcircuits.  
  // These aren't known until subcircuits are resolved in pass 2
  Util::UParamList addResolvedGlobalParams_; 

  // global data
  CircuitBlock &                circuitBlock_;
  CircuitContext *              circuitContext_;
  std::map<std::string,FileSSFPair>  & ssfMap_;
  std::list<Util::OptionBlock> &options_;
  std::string                   netlistFilename_;

  DeviceBlock                   device_;
  CircuitBlock*                 mainCircuitPtr_;
  CircuitBlock*                 parentCircuitPtr_;
  CircuitBlock*                 currentCircuitPtr_;
  std::vector<bool>             preprocessFilter_;
  SpiceSeparatedFieldTool*      ssfPtr_;
  const ParsingMgr &            parsingMgr_;

  bool                          remove_any_redundant_;  ///< check device nodes for removal
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_DistToolBase_h

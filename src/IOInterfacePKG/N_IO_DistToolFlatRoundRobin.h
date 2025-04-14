//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Declares the DistToolFlatRoundRobin class.  Distribution tool
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


#ifndef Xyce_N_IO_DistToolFlatRoundRobin_h
#define Xyce_N_IO_DistToolFlatRoundRobin_h

#include <vector>
#include <string>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_DistToolBase.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class          : DistToolFlatRoundRobin
// Purpose        : Buffers and distributes circuit blocks, device lines, and related data
//                  for/during parsing.  This class has the strategy of first-come-first-served
//                  distribution of devices over processors.
//-----------------------------------------------------------------------------
class DistToolFlatRoundRobin: public DistToolBase
{
public:
  DistToolFlatRoundRobin(
    Parallel::Communicator *                 pdsCommPtr,
    CircuitBlock &                           circuit_block,
    std::map<std::string,FileSSFPair>      & ssfMap, 
    std::map<std::string, IncludeFileInfo> & iflMap,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
    const ParsingMgr                       & parsing_manager
    );

  virtual ~DistToolFlatRoundRobin() {}

  // send options, metatdata, and context to all procs
  bool broadcastGlobalData();

  // Distribute devices using hierarchical context object.
  void distributeDevices();

  // change current subcircuit context after a new subcircuit is started
  bool circuitStart( std::string const & subcircuitName,
                     std::vector<std::string> const & nodes,
                     std::string const & prefix,
                     std::vector<Device::Param> const & params )
  { // Throw error here!
    return false;
  }

  // change current subcircuit context to previous context
  bool circuitEnd()
  { // Throw error here!
    return false;
  }

protected:

  // Parse the given include file for 2nd pass
  bool parseIncludeFile(std::string const& includeFiles,
                        const std::string &libSelect);

  // helper function for parseIncludeFile()
  void restorePrevssfInfo(
    SpiceSeparatedFieldTool* oldssfPtr,
    const std::string& old_netlistFilename,
    int oldFilePos,
    int oldLineNumber);

  bool expandSubcircuitInstance(DeviceBlock & subcircuitInstance,
                                const std::string &libSelect,
                                std::vector<std::string> &libInside)
  { // Throw error here! 
    return false;
  }

private:
  bool bufferDeviceLines(std::string fileName, 
                         std::string libSelect,
                         std::vector<std::string>& libInside,
                         long filePos, long lineNum,
                         int totalDevCount, int startIdx=0);

  // Tell neighbor processor that parsing has stopped.
  void endDeviceLines(int stopProc);

  void setFileName ( std::string const & fileNameIn );

  // wait in queue for instructions to retrieve device lines from the netlist.
  bool waitYourTurn();

  // process all the devices held in the buffer
  bool processDeviceBuffer();

private:

  int                           myProc_, fromProc_, toProc_; ///< round-robin procs

  std::map<std::string, IncludeFileInfo> & iflMap_;

  int                           procDeviceCount_;       ///< max number of devices to tx to each proc
  int                           blockSize_;             ///< size of blocks used in round robin
  int                           numBlocks_;             ///< number of blocks used in round robin
  int                           devices_;               ///< total number of devices;
  int                           myDeviceLines_;         ///< keep track of local number of devices processed
  std::vector<TokenVector>      bufs_;

  const std::vector< std::pair< std::string, std::string> >&  externalNetlistParams_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_DistToolFlatRoundRobin_h

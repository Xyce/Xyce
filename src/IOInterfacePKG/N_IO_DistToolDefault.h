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
// Purpose        : Declares the DistToolDefault class.  Distribution tool
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


#ifndef Xyce_N_IO_DistToolDefault_h
#define Xyce_N_IO_DistToolDefault_h

#include <vector>
#include <string>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_DistToolBase.h>
#include <N_IO_ParsingMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class          : DistToolDefault
// Purpose        : Buffers and distributes circuit blocks, device lines, and related data
//                  for/during parsing.  This class has the strategy of first-come-first-served
//                  distribution of devices over processors.
//-----------------------------------------------------------------------------
class DistToolDefault: public DistToolBase
{
public:
  DistToolDefault(
    Parallel::Communicator *                 pdsCommPtr,
    CircuitBlock &                           circuit_block,
    std::map<std::string,FileSSFPair>      & ssfMap, 
    std::map<std::string, IncludeFileInfo> & iflMap,
    const ParsingMgr                       & parsing_manager
    );

  virtual ~DistToolDefault() {}

  // send options, metatdata, and context to all procs
  bool broadcastGlobalData();

  // Distribute devices using hierarchical context object.
  void distributeDevices();

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

  // Expand a subcircuit instance by adding the devices and
  // device models that compose the subcircuit to the main
  // (top level) circuit. Prepend device names and nodes with
  // subcircuitPrefix.
  bool expandSubcircuitInstance(DeviceBlock & subcircuitInstance,
                                const std::string &libSelect, 
                                std::vector<std::string> &libInside);

private:

  // Send a circuit device line to current proc
  bool circuitDeviceLine(TokenVector & deviceLine );
  void endDeviceLines();

  // change current subcircuit context after a new subcircuit is started
  bool circuitStart( std::string const & subcircuitName,
                     std::vector<std::string> const & nodes,
                     std::string const & prefix,
                     std::vector<Device::Param> const & params );

  // change current subcircuit context to previous context
  bool circuitEnd();

  void setFileName ( std::string const & fileNameIn );

  // pack subcircuit data into buffer; helper for circuitStart()
  bool packSubcircuitData( std::string const & subcircuitName,
                           std::vector<std::string> const & nodes,
                           std::string const & prefix,
                           std::vector<Device::Param> const & params );

  // receive all data from proc 0
  bool receiveCircuitData();

  // send buffer from proc 0
  void send(int size = -1);

  // process all the devices held in the buffer
  bool processDeviceBuffer();

private:
  int                           currProc_;              ///< proc that will receive the next tx

  std::map<std::string, IncludeFileInfo> & iflMap_;

  int                           procDeviceCount_;       ///< max number of devices to tx to each proc
  int                           devices_;               ///< total number of devices;
  int                           deviceLinesSent_;       ///< keep track to txmitted devices for the current proc

  std::vector<char *>           bufs_;                  ///< tx/rx buffer
  std::vector<int>              bufSize_;               ///< buffer size

  // subcircuit data
  std::vector<std::string>                      subcircuitNames_;
  std::vector<std::vector<std::string> >        subcircuitNodes_;
  std::vector<std::string>                      subcircuitPrefixes_;
  std::vector<std::vector<Device::Param> >      subcircuitParams_;

};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_DistToolDefault_h

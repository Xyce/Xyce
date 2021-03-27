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

//-----------------------------------------------------------------------------
//
// Purpose        : Non-member functions that help the parser.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/06/2001
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_ParsingHelpers_h
#define Xyce_N_IO_ParsingHelpers_h

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>

#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_DeviceBlock.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace IO {

// translate a relative path for an include file into the correct relative path
void handleIncludeFilePath(
   const std::string& netlistFileName,
   std::string& includeFileName);

// helper functions for handleIncludeFilePath()
bool isAbsolutePath(const std::string& includeFile);
bool hasWinDriveLetter(const std::string& includeFile);

// Handle a netlist .include or .lib line, return the include file, and lib strings.
void handleIncludeLine(
    const std::string& netlistFileName,
    TokenVector const& parsedLine,
    const ExtendedString &, std::string& includefile, std::string &libSelect, std::string &libInside);

// Handle a netlist .endl line.
void handleEndlLine(
    const std::string& netlistFileName,
    TokenVector const& parsedLine,
    const std::string &libInside);

// Read a line of input from an istream
void readLine( std::istream & in, std::string& line );

// Split the input charLine into fields that are stored in splitLine.
// NOTE:  This method does NOT skip comments and blank lines.
void splitLine( const std::string& charLine, TokenVector & line );

// Split the input charLine into fields that are stored in splitLine.
// NOTE:  This method does NOT skip comments and blank lines, but no
//        whitespace characters are included.
void splitLineNoWS( const std::string& charLine, TokenVector & line );

// Remove two terminal devices if both nodes on the device are the same
bool removeThreeTerminalDevice(const std::vector<bool>& pFilter,
                               const char& linetype,
                               const ExtendedString & node1,
                               const ExtendedString & node2,
                               const ExtendedString & node3);

// Remove two terminal devices if both nodes on the device are the same
bool removeTwoTerminalDevice(const std::vector<bool>& pFilter,
                             const char& linetype,
                             const ExtendedString & node1,
                             const ExtendedString & node2);

bool removeDevice(const std::vector<bool>& pFilter,
                  const DeviceBlock& device);

void readExternalParamsFromFile( Parallel::Communicator& comm,
                                 std::string filename, 
                                 std::vector< std::pair< std::string, std::string > > & paramList );

bool extractParamData( CircuitBlock &            circuit_block,
                       const std::string &       netlist_filename,
                       const TokenVector &       parsed_line);

bool extractGlobalParamData( CircuitBlock &            circuit_block,
                             const std::string &       netlist_filename,
                             const TokenVector &       parsed_line);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_ParsingHelpers_h

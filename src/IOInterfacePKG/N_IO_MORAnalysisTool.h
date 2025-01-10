//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Declares the MORAnalysisTool class.  Finds linear reducible
//                  subcircuits and then performs model order reduction (MOR)
//                  on them.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 09/10/2014
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef N_IO_MORANALYSISTOOL_H
#define N_IO_MORANALYSISTOOL_H

#include <vector>
#include <string>

#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_IO_CircuitContext.h>
#include <N_UTL_NetlistLocation.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class          : MORAnalysisTool
// Purpose        : Analyzes CircuitContext object to find reducible linear subckts
//                  and then reduces them.
//-----------------------------------------------------------------------------
class MORAnalysisTool
{

public:
  // ctor
  MORAnalysisTool( CircuitContext & circuitContext,
                   std::list<Util::OptionBlock> & optionsTable,
                   const CmdParse & commandLine,
                   std::map<std::string, FileSSFPair> & ssfMap,
                   std::map<std::string, IncludeFileInfo> & iflMap )
  : ssfMap_(ssfMap),
    iflMap_(iflMap),
    circuitContext_(circuitContext),
    optionsTable_(optionsTable),
    commandLine_(commandLine),
    analysisOnly_(false)
  {
  }
    
  // dtor
  ~MORAnalysisTool() {}

  // Find linear subcircuits within the CircuitContext object.
  bool reduceLinearSubcircuits();

  // Remove the .MOR and .OPTIONS MOR_OPTS lines from the OptionsBlock object
  void removeMOROptions();

private:

  struct sort_by_line
  {
    bool operator()(const std::pair<int,std::string> &left, const std::pair<int,std::string> &right)
    {
        return left.first < right.first;
    }
  };

  void collectReducibleSubcircuits( const std::string& name );
  void generateMORFiles();
  void runFiles( const std::vector< std::string >& files );
  void collectReducedSubcircuits( const std::vector< std::string >& reducedNames );
  void collectIncludeFileLocations( std::vector< NetlistLocation >& includeFiles );
  void generateReintegratedFiles( std::string & newNetlistName );

  // Container for viable subckts and their children, providing the location of the subckt.
  std::map< std::string, std::map< std::string, NetlistLocation > > redSubcktNames_;

  // Container for node lists for viable subckts.
  std::map< std::string, std::vector< std::string > > redSubcktNodes_;

  // Container for models for viable subckts.
  std::map< std::string, std::list<ModelMap> > redSubcktModels_;

  // Container for user-specified subcircuit names.
  std::vector< std::string > desiredSubckts_;

  // Netlist file istream pairs.
  std::map< std::string, FileSSFPair > ssfMap_;

  // Include file location map.
  std::map< std::string, IncludeFileInfo > iflMap_;

  // Map of reduced subcircuit to locations in original circuit.
  std::map< std::string, std::vector< NetlistLocation > > usedSubcktLocations_;

  // circuit hierarchy
  CircuitContext & circuitContext_;

  // options table
  std::list<Util::OptionBlock> & optionsTable_;
 
  // command line
  const CmdParse & commandLine_;
 
  // MOR files
  std::vector< std::string > morFilenames_;

  // MOR options
  std::string morOptionsLine_;

  // Whether to perform reduction or only analysis (useful for debugging)
  bool analysisOnly_;

};

} // namespace IO
} // namespace Xyce

#endif //N_IO_MORANALYSISTOOL_H

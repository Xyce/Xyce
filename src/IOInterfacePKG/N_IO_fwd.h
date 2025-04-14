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

//-------------------------------------------------------------------------
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_IO_fwd_h
#define Xyce_N_IO_fwd_h

#include <unordered_map>
using std::unordered_map;

#include <map>
#include <string>
#include <vector>
#include <utility>

#include <N_UTL_NoCase.h>
#include <N_UTL_NetlistLocation.h>

namespace Xyce {
namespace IO {

namespace PreprocessType {
enum PreprocessType {REDUNDANT_C, REDUNDANT_D, REDUNDANT_I, REDUNDANT_L, REDUNDANT_M, REDUNDANT_Q, REDUNDANT_R, REDUNDANT_V, REPLACE_GROUND, NUM_PREPROCESS};
}

namespace DistStrategy {
enum DistStrategy {DEFAULT, FLAT_ROUND_ROBIN, DEVICE_BALANCED, NUM_STRATEGIES};
}

class ActiveOutput;
class CircuitBlock;
class CircuitContext;
class CircuitMetadata;
class CmdParse;
class DeviceBlock;
class DistributionTool;
class FourierMgr;
class FFTMgr;
class FFTAnalysis;
class FunctionBlock;
class InitialConditionsManager;
class HangingResistor;
class NetlistImportTool;
struct LoadManager;
class OutputMOR;
class OutputMOR;
class OutputMgr;
class OutputResponse;
class OutputResults;
class ParameterBlock;
class ParsingMgr;
class PkgOptionsMgr;
class RestartMgr;
class RestartNode;
class SpiceSeparatedFieldTool;
struct StringToken;
struct DeviceMetadata;
struct PkgOptionsReg;

class OutputFileBase;

struct PrintParameters;
struct Table;

typedef StringToken Token;
typedef std::vector<Token> TokenVector;

typedef bool (*ParseFunction)(PkgOptionsMgr &options_manager, CircuitBlock &circuit_block, const std::string &netlist_filename, const TokenVector &parsed_line);

typedef std::pair<double, double> Interval;
typedef std::vector<Interval> IntervalVector;
typedef unordered_map<std::string, std::string, HashNoCase, EqualNoCase> AliasNodeMap;
typedef std::pair<std::ifstream *, SpiceSeparatedFieldTool *> FileSSFPair;

typedef std::map<std::string, ParameterBlock *, LessNoCase> ModelMap;

namespace Measure {
class Manager;
}

class FourierMgr;

// Information stored about the include file during the first pass of the parser.
struct IncludeFileInfo
{
  IncludeFileInfo()
  :
  numDevices(0),
  numSubckts(0),
  numModels(0),
  numSUBCKTdefs(0),
  inSUBCKT(false),
  parentSUBCKT(""),
  location(NetlistLocation())
  {}

  virtual ~IncludeFileInfo()
  {}

  int numDevices;
  int numSubckts;
  int numModels;
  int numSUBCKTdefs;
  bool inSUBCKT;
  std::string parentSUBCKT;
  NetlistLocation location;
};

// Class used to wrap objects from user codes and the thing it wraps
class ExternalOutputWrapper;
class ExternalOutputInterface;
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_fwd_h


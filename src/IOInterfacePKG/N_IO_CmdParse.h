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
// Filename      : N_IO_CmdParse.h
//
// Purpose       : This file defines the command line parser class.
//
// Special Notes :
//
//
// Creator       : Eric Keiter
//
// Creation Date : 06/17/99
//
//
//-----------------------------------------------------------------------------

#ifndef _N_IO_CMDPARSE_
#define _N_IO_CMDPARSE_

#include <iosfwd>
#include <string>
#include <map>

#include <N_IO_fwd.h>
#include <N_DEV_fwd.h>
#include <N_PDS_fwd.h>

#include <N_IO_HangingResistor.h>

namespace Xyce {
namespace IO {

void usage(std::ostream &os);
void setTimeIntegratorDebugLevel(const CmdParse &command_line, int level);
void setDeviceDebugLevel(const CmdParse &command_line, int level);
void setSensitivityDebugLevel(const CmdParse &command_line, int level);

std::string makeOutputFileName(
  const CmdParse &      command_line,
  const std::string&    suffix);

std::string makeOutputFileNameWithStepNum(
  const CmdParse &      command_line,
  const std::string&    suffix,
  int                   stepNumber);

std::string removeExtensionsFromDashoFilename(std::string dashoFilename);

//-----------------------------------------------------------------------------
// Class         : CmdParse
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
class CmdParse
{
public:
  enum ParseState {ERROR, SUCCESS, DONE};

  CmdParse();

  virtual ~CmdParse();

  CmdParse(const CmdParse &);

private:
  CmdParse &operator=(const CmdParse &);

public:
  int argc() const
  {
    return argc_;
  }

  char **argv() const
  {
    return argv_;
  }

  ParseState parseCommandLine(Parallel::Machine comm, int argc, char **argv);

  bool argExists(const std::string & arg_tmp) const;

  std::string getArgumentValue(const std::string & argumentName) const;
  int getArgumentIntValue(const std::string & argument_name, int default_value) const;

  // Returns the number of cmd line args.
  int numArgs();

  void printArgMap();

  const char * getArg(int i);

  void setNetlist(const std::string &newNetlist);

private:
  void setCommandArgs();

  bool isSwitchArg(const std::string &arg);
  bool isStringValuedArg(const std::string &arg);

private :
  int           argc_;
  char **       argv_;
  char **       myArgv_;

  std::map<std::string, int>            swArgs;
  std::map<std::string, std::string>    stArgs;

  std::map<std::string, int>            argIndex;

  HangingResistor                       hangingResistor_;
};

} // namespace IO
} // namespace Xyce

#endif

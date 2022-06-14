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
// Filename      : N_IO_CmdParse.C
//
// Purpose       : This file contains the command line parser class.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 06/17/99
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <sstream>
#include <cctype>
#include <cstring>
#include <regex>

#include <N_IO_CmdParse.h>
#include <N_IO_PrintTypes.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_Print.h>
#include <N_DEV_RegisterDevices.h>
#include <N_ERH_Message.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Progress.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {

void
setTimeIntegratorDebugLevel(
  const CmdParse &      command_line,
  int                   level)
{
  Xyce::setTimeIntegratorDebugLevel(command_line.getArgumentIntValue("-tdl", DEBUG_TIME ? level : 0));
}

void
setDeviceDebugLevel(
  const CmdParse &      command_line,
  int                   level)
{
  Xyce::setDeviceDebugLevel(command_line.getArgumentIntValue("-ddl", DEBUG_TIME ? level : 0));
}

void
setSensitivityDebugLevel(
  const CmdParse &      command_line,
  int                   level)
{
  Xyce::setDeviceSensitivityDebugLevel(command_line.getArgumentIntValue("-sdl", DEBUG_TIME ? level : 0));
}

//-----------------------------------------------------------------------------
// Function      : makeOutputFileName
// Purpose       : make the output file name
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 05/28/2021
//-----------------------------------------------------------------------------
std::string makeOutputFileName(
  const CmdParse &      command_line,
  const std::string&    suffix)
{
  std::string netlist_filename = command_line.getArgumentValue("netlist");
  std::string dasho_filename = command_line.getArgumentValue("-o");
  std::string base_filename = (!dasho_filename.empty()) ? dasho_filename : netlist_filename;

  std::string file_name = base_filename + suffix;
  if (file_name == netlist_filename)
    file_name = file_name + suffix;

  return file_name;
}

//-----------------------------------------------------------------------------
// Function      : makeOutputFileNameWithStepNum
// Purpose       : make the output file mame, when that name includes the step
//                 number as the trailing part of the file suffix
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 05/28/2021
//-----------------------------------------------------------------------------
std::string makeOutputFileNameWithStepNum(
  const CmdParse &      command_line,
  const std::string&    suffix,
  int                   stepNumber)
{
  std::ostringstream converterBuff;
  converterBuff << stepNumber;

  return makeOutputFileName(command_line, suffix + converterBuff.str());
}

//-----------------------------------------------------------------------------
// Function      : removeExtensionsFromDashoFilename
// Purpose       : removes various "well-known .PRINT line extensions" from
//                the -o filename requested on the command line.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 06/28/2021
//-----------------------------------------------------------------------------
std::string removeExtensionsFromDashoFilename(std::string dashoFilename)
{
  std::regex e;

  for (int i =0; i<sizeof(dasho_print_extensions_regex)/sizeof(dasho_print_extensions_regex[0]); i++)
  {
    std::string regexString(dasho_print_extensions_regex[i]);

    // assume that at least one character must come before the "well-known extension".
    try
    {
      e.assign("(.+)" + regexString);
    }
    catch (std::regex_error& regexErr)
    {
      Report::DevelFatal().in("removeExtensionsFromDashoFilename") <<
      "Error converting " << "(.+)" + regexString << " into std::regex object";
    }

    if (std::regex_match(dashoFilename, e))
    {
      dashoFilename = std::regex_replace(dashoFilename, e, "$1");
      break;
    }
  }

  return dashoFilename;
}
//-----------------------------------------------------------------------------
// Function      : usage
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : David Baur
// Creation Date : 3/7/2014
//-----------------------------------------------------------------------------
///
/// Print Xyce usage
///
/// @param os   Stream to write usage information to
///
void usage(std::ostream &os)
{
  os << "Usage: Xyce [arguments] netlist\n\n"
     << "Arguments:\n"
     << "  -b                          batch mode flag for spice compatibility (ignored)\n"
     << "  -h                          print usage and exit\n"
     << "  -v                          print version info and exit\n"
     << "  -capabilities               print compiled-in options and exit\n"
     << "  -license                    print license and exit\n"
     << "  -param [device [level [<inst|mod>]]] print a terse summary of model and/or device parameters\n"
     << "  -doc [device [level [<inst|mod>]]] output latex tables of model and device parameters to files\n"
     << "  -doc_cat [device [level [<inst|mod>]]] output latex tables of model and device parameters to files\n"
     << "  -count                      device count without netlist syntax or topology check\n"
     << "  -syntax                     check netlist syntax and exit\n"
     << "  -norun                      netlist syntax and topology and exit\n"
     << "  -namesfile <path>           output internal names file to <path> and exit\n"
     << "  -noise_names_file <path>    output noise source names file to <path> and exit\n"
     << "  -quiet                      suppress some of the simulation-progress messages sent to stdout\n"
     << "  -jacobian_test              jacobian matrix diagnostic\n"
     << "  -hspice-ext  <option>       apply hspice compatibility features during parsing.  option=all applies them all\n"
     << "  -redefined_params <option>  set option for redefined .params as ignore (use last), usefirst, warn or error\n"
     << "  -subckt_multiplier <option> set option to true(default) or false to apply implicit subcircuit multipliers\n"
     << "  -delim <TAB|COMMA|string>   set the output file field delimiter\n"
     << "  -o <basename>               <basename> for the output file(s)\n"
     << "  -l <path>                   place the log output into <path>, \"cout\" to log to stdout\n"
     << "  -per-processor              create log file for each procesor, add .<n>.<r> to log path\n"
     << "  -remeasure [existing Xyce output file] recompute .measure() results with existing data\n"
     << "  -nox <on|off>               NOX nonlinear solver usage\n"
     << "  -linsolv <solver>           force usage of specific linear solver\n"
     << "  -maxord <1..5>              maximum time integration order\n"
     << "  -prf <param file name>      specify a file with simulation parameters\n"
     << "  -rsf <response file name>   specify a file to save simulation responses functions.\n"
     << "  -r <file>                   generate a rawfile named <file> in binary format\n"
     << "  -a                          use with -r <file> to output in ascii format\n"
     << "  -randseed <number>          seed random number generator used by expressions and sampling methods\n"

#ifdef HAVE_DLFCN_H
     << "  -plugin <plugin list>       load device plugin libraries (comma-separated list)\n"
#endif

#ifdef Xyce_Dakota
     << "\nDakota arguments:\n"
     << "  -dakota <dakota input file> dakota input file for this simulation\n"
#endif

#if defined Xyce_DEBUG_NONLINEAR || defined Xyce_DEBUG_TIME || defined Xyce_DEBUG_DEVICE || defined Xyce_TOTALVIEW_BOGON
     << "\nDebug arguments:\n"
#endif

#ifdef Xyce_DEBUG_NONLINEAR
     << "  -ndl <0..N>                 set the nonlinear solver debug level (overrides netlist)\n"
#endif

#ifdef Xyce_DEBUG_TIME
     << "  -tdl <0..N>                 set the time integration debug level\n"
#endif
     << (DEBUG_DEVICE ? "  -ddl <0..N>                 set the device debug level (overrides netlist)\n" : "")
     << (DEBUG_DEVICE ? "  -sdl <0..N>                 set the device sensitivity debug level\n" : "")

#ifdef Xyce_TOTALVIEW_BOGON
     << "  -mpichtv                    required for totalview \n"
#endif
     << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::CmdParse
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
CmdParse::CmdParse()
: argc_(0),
  argv_(0),
  myArgv_(0)
{
  setCommandArgs();
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::CmdParse
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 01/07/07
//-----------------------------------------------------------------------------
CmdParse::CmdParse(const CmdParse & right) :
  argc_(right.argc_),
  swArgs(right.swArgs),
  stArgs(right.stArgs),
  argIndex(right.argIndex)
{
  myArgv_ = argv_ = new char*[argc_];
  for (int i = 0; i < argc_; ++i)
  {
    argv_[i] = new char[std::strlen(right.argv_[i]) + 1];
    std::strcpy(argv_[i], right.argv_[i]);
  }
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::CmdParse
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
CmdParse::~CmdParse()
{
  if (myArgv_)
  {
    for (int i = 0; i < argc_; i++)
      delete[] argv_[i];

    delete[] myArgv_;
  }
}


//-----------------------------------------------------------------------------
// Function      : CmdParse::setCommandArgs
// Purpose       : Initialize command line arguments
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 01/24/11
//-----------------------------------------------------------------------------
void
CmdParse::setCommandArgs()
{
  stArgs.clear();
  swArgs.clear();
  argIndex.clear();

  stArgs[ "netlist" ] = "";     // The input netlist will be a special case.

  // Set the maps containing the expected arguments. There are two types
  // of arguments switch arguments and string valued arguments.
  swArgs[ "-b" ] = 0;
  swArgs[ "-h" ] = 0;           // Help option, print usage and exit.
  swArgs[ "-test" ] = 0;
  swArgs[ "-v" ] = 0;
  swArgs[ "-capabilities" ] = 0;
  swArgs[ "-license" ] = 0;
  swArgs[ "-count" ] = 0;
  swArgs[ "-syntax" ] = 0;
  swArgs[ "-norun" ] = 0;
  stArgs[ "-namesfile" ] = "";
  stArgs[ "-noise_names_file" ] = "";
  swArgs[ "-quiet" ] = 0;
  swArgs[ "-jacobian_test" ] = 0;
  swArgs[ "-error-test" ] = 0;
  swArgs[ "-hspice-ext"] = 0;
  swArgs[ "-redefined_params"] = 0;
  swArgs[ "-subckt_multiplier"] = 0;

  stArgs[ "-delim" ] = "";
  stArgs[ "-o" ] = "";
  stArgs[ "-l" ] = "";                  // Output log information to a file.
  stArgs[ "-verbose" ] = "";            // Output verbose information to seperate file.
  swArgs[ "-per-processor" ] = 0;
  swArgs[ "-messy-cout" ] = 0;
  stArgs[ "-remeasure" ] = "";          // specify a existing data file on which Xyce will recompute .measure functions.
  stArgs[ "-nox" ] = "";
  stArgs[ "-linsolv" ] = "";
  stArgs[ "-maxord" ] = "";
  stArgs[ "-prf" ] = "";                // specify a parameter input file to set runtime params from a file
  stArgs[ "-rsf" ] = "";                // specify a response output file to save results to a file
  stArgs[ "-r" ] = "";                  // Output binary rawfile.
  swArgs[ "-a" ] = 0;                   // Use ascii instead of binary in rawfile output
  stArgs[ "-randseed" ] = "";           // random number seed
  
#ifdef HAVE_DLFCN_H
  stArgs[ "-plugin" ] = "";
#endif

#ifdef Xyce_Dakota
  stArgs[ "-dakota" ] = "";     // specify a dakota input file to couple with this simulation
#endif

#ifdef Xyce_DEBUG_NONLINEAR
  stArgs[ "-ndl" ] = "";
#endif
#ifdef Xyce_DEBUG_TIME
  stArgs[ "-tdl" ] = "";
#endif
  if (DEBUG_DEVICE)
  {
    stArgs[ "-ddl" ] = "";
    stArgs[ "-sdl" ] = "";
  }
#ifdef Xyce_TOTALVIEW_BOGON
  swArgs[ "-mpichtv" ] = 0;
#endif
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::setNetlist
// Purpose       : This function is used for 2-level newton solves.
//
// Special Notes : For 2-level newton, the inner solve is handled by a new
//                 allocation of Circuit::Simulator, which needs to be passed in
//                 a set of command line args.  Those command line args
//                 need to be identical to the ones passed in by the user,
//                 except that the netlist needs to be different, as the
//                 inner problem is described in a different file.
//
//                 This function is called after a copy of the original
//                 CmdParse class has been made.
//
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 01/07/07
//-----------------------------------------------------------------------------
void
CmdParse::setNetlist(
  const std::string &   newNetlist)
{
  int netIndex = 0;

  if ( argIndex.find("netlist") != argIndex.end())
  {
    netIndex = argIndex["netlist"];
  }
  else
  {
    Report::DevelFatal0().in("CmdParse::setNetlist") << "Unable to find netlist argument.";
  }


#if 0
  //++netIndex;

  if ( netIndex >= argc_ )
  {
    // Unexectedly ran out of arguments on the command line.
    Report::DevelFatal0().in("CmdParse::setNetlist") << "Did not find previous netlist setting.";
  }
  else if (argv_[netIndex][0] == '-')
  {
    // Error if we ran into another option here.
    Report::DevelFatal0().in("CmdParse::setNetlist") << "Expected option value, but found option " << argv_[netIndex];
  }
  else // found it!
#endif

  {
    delete [] argv_[netIndex];

    int newSize = newNetlist.size()+2;
    argv_[netIndex] = new char[newSize];
    for (int i=0;i<newSize;++i)  argv_[netIndex][i] = 0;

    sprintf(argv_[netIndex], "%s", newNetlist.c_str());

    stArgs["netlist"] = newNetlist;
  }

}

//-----------------------------------------------------------------------------
// Function      : CmdParse::parseCommandLine
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters
// Creation Date : 03/13/02
//-----------------------------------------------------------------------------
CmdParse::ParseState
CmdParse::parseCommandLine(
  Parallel::Machine     comm,
  int                   argc,
  char **               argv)
{
  argc_ = argc;
  argv_ = argv;

  int parse_state = SUCCESS;

  setCommandArgs();

  if (Parallel::rank(comm) == 0)
  {
    for (int i = 1; i < argc; ++i)
    {
      // MPICH sometimes creates arguments that get replaced by NULLs by
      // MPI_INIT :-P.  They are always after the arguments that we care about,
      // so stop command line parsing at that point.
      if (argv[i] == NULL) {
        argc_ = i;
        break;
      }

      const std::string arg = argv[i];

      argIndex[arg] = i;

      if (arg[0] == '-')
      {
        if (arg == "-h" )
        {
          usage(Xyce::lout());
          parse_state = DONE;
          break;
        }
        else if (arg == "-v" )
        {
          Xyce::lout() << Util::Version::getFullVersionString() << std::endl;
          parse_state = DONE;
          break;
        }
        else if (arg == "-capabilities" )
        {
          Xyce::lout() << Util::Version::getCapabilities() << std::endl;
          parse_state = DONE;
          break;
        }
        else if (arg == "-license" )
        {
          Xyce::lout() << Util::Version::getLicense() << std::endl;
          parse_state = DONE;
          break;
        }
        else if (arg == "-plugin")
        {
          ++i;

          if ( i >= argc )
          {
            // Unexectedly ran out of arguments on the command line.
            Xyce::lout() << "Did not find required value for option " << arg << std::endl;
            parse_state = ERROR;
            continue;
          }
          else if (argv[i][0] == '-')
          {
            // Error if we ran into another option here.
            Xyce::lout() << "Expected option value, but found option " << argv[i] << std::endl;
            parse_state = ERROR;
            continue;
          }
          else
          {
            stArgs[arg] = argv[i];
          }
        }
        else if (arg == "-hspice-ext")
        {
          ++i;

          if ( i >= argc )
          {
            // Unexectedly ran out of arguments on the command line.
            Xyce::lout() << "Did not find required value for option " << arg << std::endl;
            parse_state = ERROR;
            continue;
          }
          else if (argv[i][0] == '-')
          {
            // Error if we ran into another option here.
            Xyce::lout() << "Expected option value, but found option " << argv[i] << std::endl;
            parse_state = ERROR;
            continue;
          }
          else
          {
            stArgs[arg] = argv[i];
          }
        }
        else if (arg == "-redefined_params")
        {
          ++i;

          if ( i >= argc )
          {
            // Unexectedly ran out of arguments on the command line.
            Xyce::lout() << "Did not find required value for option " << arg << std::endl;
            parse_state = ERROR;
            continue;
          }
          else if (argv[i][0] == '-')
          {
            // Error if we ran into another option here.
            Xyce::lout() << "Expected option value, but found option " << argv[i] << std::endl;
            parse_state = ERROR;
            continue;
          }
          else
          {
            stArgs[arg] = argv[i];
          }
        }
        else if (arg == "-subckt_multiplier")
        {
          ++i;

          if ( i >= argc )
          {
            // Unexectedly ran out of arguments on the command line.
            Xyce::lout() << "Did not find required value for option " << arg << std::endl;
            parse_state = ERROR;
            continue;
          }
          else if (argv[i][0] == '-')
          {
            // Error if we ran into another option here.
            Xyce::lout() << "Expected option value, but found option " << argv[i] << std::endl;
            parse_state = ERROR;
            continue;
          }
          else
          {
            stArgs[arg] = argv[i];
          }
        }
        else if (arg == "-o")
        {
          ++i;

          if ( i >= argc )
          {
            // Unexectedly ran out of arguments on the command line.
            Xyce::lout() << "Did not find required value for option " << arg << std::endl;
            parse_state = ERROR;
            continue;
          }
          else if (argv[i][0] == '-')
          {
            // Error if we ran into another option here.
            Xyce::lout() << "Expected option value, but found option " << argv[i] << std::endl;
            parse_state = ERROR;
            continue;
          }
          else
          {
            stArgs[arg] = removeExtensionsFromDashoFilename(argv[i]);
          }
        }
        else if (arg == "-param" || 
                 arg == "-doc"   || arg == "-doc_cat" )
        {
          std::string option_device_name="";
          int option_device_level = -1;
          std::string option_device_level_string = "-1";
          bool print_model = true;
          bool print_instance = true;
          std::string option_flags;
          if (i < argc - 1 && argv[i + 1][0] != '-')
          {
            option_device_name = argv[i + 1];
            ++i;
            if (i < argc - 1 && argv[i + 1][0] != '-')
            {
              option_device_level = atoi(argv[i + 1]);
              option_device_level_string = argv[i+1];
              if (option_device_level != 0)
              {
                ++i;
                if (i < argc - 1 && argv[i + 1][0] != '-')
                {
                  if (tolower(argv[i + 1][0]) == 'm')
                    print_instance = false;
                  else if (tolower(argv[i + 1][0]) == 'i')
                    print_model = false;
                  else
                    --i;
                  ++i;
                }
              }
            }
          }
          // Now, simply place what we found on the command line into the
          // args map.  Unfortunately, there is no such thing as an integer
          // args array, only a string and a flag array.  So we are doing
          // some goofy stuff here in order to allow the higher level code
          // to process -param and friends, which is necessary for them
          // to work with shared library plugins in parallel.
          
          // kludge: argExists will return false if we put "" in there.
          // this means we have to fix it later.  Ugh.
          if (option_device_name == "") option_device_name=" ";
          stArgs[arg] = option_device_name;

          stArgs[arg+"_level"] = option_device_level_string;

          // Ugh... trying to avoid wasteful use of stringstream here
          stArgs[arg+"_flags"] = (print_model&&print_instance)?std::string("3"):(print_model)?std::string("2"):(print_instance)?std::string("1"):std::string("0");
        }
        else if (argExists(arg))
        {
          Report::UserWarning0() << "More than one \"" << arg << "\" command line argument found, using last one";
        }
        else if (isSwitchArg(arg))
        {
          swArgs[arg] = i;
        }
        else if (isStringValuedArg(arg))
        {
          ++i;

          if ( i >= argc )
          {
            // Unexectedly ran out of arguments on the command line.
            Xyce::lout() << "Did not find required value for option " << arg << std::endl;
            parse_state = ERROR;
          }
          else if (argv[i][0] == '-')
          {
            // Error if we ran into another option here.
            Xyce::lout() << "Expected option value, but found option " << argv[i] << std::endl;
            parse_state = ERROR;
          }
          else
          {
            stArgs[arg] = argv[i];
          }
        }
        else
        {
          // Invalid option, stop here.
          Xyce::lout() << "Invalid option given: " << arg << std::endl;
          usage(Xyce::lout());
          parse_state = ERROR;
        }
      }
      else
      {
        if (stArgs["netlist"] == "")
        {
          // Assume this field is the input netlist.
          stArgs["netlist"] = arg;
          argIndex["netlist"] = i;
        }
        else
        {
          // Already found netlist, report error and terminate.
          // or this could be the case of a Dakota call to Xyce of
          // the format Xyce mycircuit.cir -prf params.in results.out
          // The more correct way would be to use the "-rsf <filename>"
          // to specify the results file.
          // This is a bit of a hack, but if the extra netlist starts with
          // res* (as in results* or response*) then assume it is the response
          // file and save it with the tag "-rsf" (provided that tag has not
          // been used).
          if ((stArgs["-rsf"]=="") &&
              ((arg.find_first_of("res") != std::string::npos)
               || (arg.find_first_of("RES") != std::string::npos)) )
          {
            // this is likely a response file name.  Store it as such
            stArgs["-rsf"]=arg;
          }
          else
          {
            Xyce::lout() << "Found second netlist on command line: " << arg << std::endl;
            parse_state = ERROR;
          }
        }
      }
    }
  }

  N_ERH_ErrorMgr::safeBarrier(comm);

  if (Parallel::is_parallel_run(comm)) {
    if (Parallel::rank(comm) == 0) {
      Util::Marshal mout; // (Util::Marshal::TYPE_CHECK_ALL);
      mout << parse_state << swArgs << stArgs;
      Parallel::Broadcast(comm, mout, 0);
    }
    else {
      Util::Marshal min("");
      Parallel::Broadcast(comm, min, 0);
      swArgs.clear();
      stArgs.clear();
      min >> parse_state >> swArgs >> stArgs;
    }
  }

  return (ParseState) parse_state;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::numArgs
// Purpose       : Returns the number of cmd line args.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
int CmdParse::numArgs()
{
  return argc_;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::printArgMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
void CmdParse::printArgMap()
{
  std::map<std::string,std::string>::iterator iter;
  std::map<std::string,std::string>::iterator begin = stArgs.begin();
  std::map<std::string,std::string>::iterator end   = stArgs.end();

  Xyce::dout() << std::endl << "Command Line Argument Map:" << std::endl;
  Xyce::dout() << std::endl;

  for (iter=begin;iter!=end;++iter)
  {
    Xyce::dout() << "   map[ ";
    Xyce::dout() << (iter->first);
    Xyce::dout() << " ] = ";
    Xyce::dout() << (iter->second) << std::endl;
  }
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::argExists
// Purpose       : This function returns true if the specified
//                 argument exists on the command line. It returns
//                 false if either the specified argument does not exist
//                 on the command line or there is no such option.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
bool CmdParse::argExists(const std::string & arg_tmp) const
{
  std::map<std::string,int>::const_iterator it = swArgs.find(arg_tmp);
  if (it != swArgs.end() && (*it).second != 0)
    return true;
  else {
    std::map<std::string,std::string>::const_iterator it = stArgs.find(arg_tmp);
    if (it == stArgs.end())
      return false;
    else
      return (*it).second != "";
  }
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::getArgumentValue
// Purpose       : This function returns the value of the specified
//                 string argument
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters
// Creation Date : 03/14/2002
//-----------------------------------------------------------------------------
std::string CmdParse::getArgumentValue(const std::string & argumentName) const
{
  std::map<std::string,std::string>::const_iterator it = stArgs.find(argumentName);

  return it == stArgs.end() ? "" : (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::getArgumentIntValue
// Purpose       : This function returns the value of the specified
//                 string argument
// Special Notes :
// Scope         : Public
// Creator       : David G. Baur
// Creation Date : 11/12/2014
//-----------------------------------------------------------------------------
int CmdParse::getArgumentIntValue(const std::string & argument_name, int default_value) const
{
  std::map<std::string,std::string>::const_iterator it = stArgs.find(argument_name);

  return it == stArgs.end() || (*it).second.empty() ? default_value : atoi((*it).second.c_str());
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::isSwitchArg
// Purpose       : This function returns true if the specified argument
//                 is a switch argument
// Special Notes :
// Scope         : Private
// Creator       : Lon Waters
// Creation Date : 03/13/2002
//-----------------------------------------------------------------------------
bool CmdParse::isSwitchArg(const std::string &arg)
{
  return swArgs.count(arg) >0;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::isStringValuedArg
// Purpose       : This function returns true if the specified argument
//                 is a string valued argument
// Special Notes :
// Scope         : Private
// Creator       : Lon Waters
// Creation Date : 03/13/2002
//-----------------------------------------------------------------------------
bool CmdParse::isStringValuedArg(const std::string &arg)
{
  return stArgs.count(arg) > 0;
}


//-----------------------------------------------------------------------------
// Function      : CmdParse::getArg
// Purpose       :
// Special Notes :
// Scope         : Private
// Creator       : Eric Keiter
// Creation Date : 10/24/2006
//-----------------------------------------------------------------------------
const char *CmdParse::getArg(int i)
{
  if (i < argc_)
    return argv_[i];
  else
    return 0;
}

} // namespace IO
} // namespace Xyce

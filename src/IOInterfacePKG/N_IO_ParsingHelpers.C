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

#include <Xyce_config.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include <ctime>
#include <cstring>

#include <N_DEV_Configuration.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_ParsingHelpers.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Expression.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

namespace Xyce {
namespace IO {

//--------------------------------------------------------------------------
// Function      : handleIncludeFilePath
// Purpose       : Translate a relative path for an include file into the
//                 correct relative path.  See SON Bug 1325 for more details.
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 11/04/2020
//--------------------------------------------------------------------------
void handleIncludeFilePath(
   const std::string& topLevelPath,
   const std::string& netlistFileName,
   std::string& includeFile)
{
  if (includeFile.empty())
    return;

  if (DEBUG_IO)
  {
    Xyce::dout() << "In IO::handlIncludeFilePath:" << std::endl;
    Xyce::dout() << "  netlistFileName = " << netlistFileName << std::endl;
    Xyce::dout() << "  includeFile name = " << includeFile << std::endl;
  }

  std::string netlistFilePath = getPathFromFileName(netlistFileName);

  if (isAbsolutePath(includeFile))
  {
    // Include file path is absolute.  So, use it as is.
  }
  else if (hasWinDriveLetter(includeFile))
  {
    // This catches Windows paths of the form C:dirname, which are relative
    // to the current execution directory.
  }
  else if (!netlistFilePath.empty())
  {
    // netlistFilename_ has a relative path.  So, make two file names.  The first priority
    // is the path that is relative to the subdirectory of netlistFilename_.  The second
    // priority is the path that is relative to the subdirectory of the top-level netlist.
    std::string includeFileWithRP = netlistFilePath + includeFile;
    std::string includeFileWithTLpath = topLevelPath + includeFile;

    // Modify the include file name if either path exists.  Otherwise, Xyce will
    // look in its execution subdirectory, by default.
    if ( Util::checkIfValidFile(includeFileWithRP) )
      includeFile = includeFileWithRP;
    else if ( !topLevelPath.empty() && Util::checkIfValidFile(includeFileWithTLpath) )
      includeFile = includeFileWithTLpath;
  }
  else
  {
    // no op, since netlistFilename_ is in the execution subdirectory, and will
    // be found by default.
  }

  if (DEBUG_IO)
    Xyce::dout() << "  Modified includeFile name = " << includeFile << std::endl;

  return;
}

//--------------------------------------------------------------------------
// Function      : getPathFromFileName
// Purpose       : Extracts the path portion from a file name. This path
//                 string may be absolute or relative. It may also be an
//                 empty string.  If non-empty, then it will be terminated
//                 with the file separator used in fileName.
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 04/05/2021
//--------------------------------------------------------------------------
std::string getPathFromFileName(const std::string& fileName)
{
  std::string path("");
  size_t posLast = fileName.find_last_of('/');

  // Also handle the Windows canonical file separator if running on Windows
  #ifdef HAVE_WINDOWS_H
  if (posLast == std::string::npos)
    posLast = fileName.find_last_of('\\');
  #endif

  if (posLast != std::string::npos)
    path = fileName.substr(0,posLast+1);

  return path;
}

//--------------------------------------------------------------------------
// Function      : isAbsolutePath
// Purpose       : Does the includeFile use an absolute path? On Windows,
//                 this also includes paths that start with C:\ and UNC paths
//                 that start with the characters `\\`
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 03/23/2021
//--------------------------------------------------------------------------
bool isAbsolutePath(const std::string& includeFile)
{
  bool retVal=false;

  size_t posFirst = includeFile.find_first_of('/');

  #ifdef HAVE_WINDOWS_H
  if (posFirst == std::string::npos)
    posFirst = includeFile.find_first_of('\\');
  #endif

  if (posFirst != std::string::npos)
  {
    if (posFirst == 0)
    {
      // include file path is absolute
      retVal=true;
    }
    #ifdef HAVE_WINDOWS_H
    else if ((includeFile.size() >=3) && (posFirst == 3) && hasWinDriveLetter(includeFile))
      retVal=true;
    #endif
  }

  return retVal;
}

//--------------------------------------------------------------------------
// Function      : hasWinDriveLetter
// Purpose       : Does the includeFile start with a Windows drive letter
//                 such as C:
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 03/23/2021
//--------------------------------------------------------------------------
bool hasWinDriveLetter(const std::string& includeFile)
{
  bool retVal=false;

  #ifdef HAVE_WINDOWS_H
  std::string driveLetters("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if ( (includeFile.size() >=2) && (includeFile[1] == ':' ) &&
       (driveLetters.find_first_of(includeFile[0]) != std::string::npos) )
  {
     retVal=true;
  }
  #endif

  return retVal;
}

//--------------------------------------------------------------------------
// Function      : handleIncludeLine
// Purpose       : Handle a netlist .include or .lib line, add the include file
//                 to includeFiles_ for later processing.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 01/10/2001
//--------------------------------------------------------------------------
void handleIncludeLine(
    const std::string& topLevelPath,
    const std::string& netlistFileName,
    TokenVector const& parsedLine, const ExtendedString & ES1,
    std::string& includeFile, std::string& libSelect, std::string& libInside)
{
  if ( parsedLine.size() < 2 )
  {
    Report::UserWarning0().at(netlistFileName, parsedLine[0].lineNumber_)
      << ES1 << " is missing argument(s), ignoring";
    return;
  }

  bool lib = ES1.substr(0,4) != ".INC";

  if (!lib || (lib && (parsedLine.size() >= 3)))
  {
    // Get the file name indicated by the .include line.
    std::string includeFileTmp = parsedLine[1].string_;
  
    // Strip off the enclosing double quotes if they are present.
    if ( (includeFileTmp[0] == '"' || includeFileTmp[0] == '{') &&
         (includeFileTmp[includeFileTmp.length()-1] == '"' || includeFileTmp[includeFileTmp.length()-1] == '}') )
    {
      includeFile = includeFileTmp.substr( 1, includeFileTmp.length()-2 );
    }
    else
    {
      includeFile = includeFileTmp;
    }

    // account for relative paths in include file name
    handleIncludeFilePath(topLevelPath, netlistFileName, includeFile);
  }
  else
  {
    includeFile = "";
  }

  if ( lib )
  {
    if ( parsedLine.size() > 3)
    {
      Report::UserWarning0().at(netlistFileName, parsedLine[0].lineNumber_)
        << "Extraneous data on .LIB ignored";
    }
    // This is a .lib line that is including selected entries of a library.
    if ( parsedLine.size() >= 3 )
    {
      libSelect = ExtendedString(parsedLine[2].string_).toUpper();
      libInside = "";
    }
    // This is a .lib line that is defining a block of netlist entries 
    // that can be selected with a '.lib <libSelect> <libInside>'
    else
    {
      libInside = ExtendedString(parsedLine[1].string_).toUpper();
      libSelect = "";
    }
  }
  else
  {
    if ( parsedLine.size() >= 3)
    {
      Report::UserWarning0().at(netlistFileName, parsedLine[0].lineNumber_)
        << "Extraneous data on .INCLUDE ignored";
    }
  }

  return;
}

//--------------------------------------------------------------------------
// Function      : handleEndlLine
// Purpose       : Handle a netlist .endl line
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 10/09/2009
//--------------------------------------------------------------------------
void handleEndlLine(
    const std::string& netlistFileName,
    TokenVector const& parsedLine,
    const std::string& libInside)
{
  if (libInside == "")
  {
    Report::UserError().at(netlistFileName, parsedLine[0].lineNumber_)
      << ".ENDL encountered without .LIB ";
    return;
  }

  if ( parsedLine.size() >= 2 )
  {
    ExtendedString libName ( parsedLine[1].string_ );
    libName.toUpper();
    if (libName != libInside)
    {
      Report::UserError().at(netlistFileName, parsedLine[0].lineNumber_)
        << ".ENDL encountered with name " << libName << ", which does not match .LIB name " << libInside;
    }
    if ( parsedLine.size() >2 )
    {
      Report::UserWarning().at(netlistFileName, parsedLine[0].lineNumber_)
        << "Extraneous field(s) following library name in .ENDL";
    }
  }
  else
  {
    Report::UserError().at(netlistFileName, parsedLine[0].lineNumber_)
      << ".ENDL encountered without library name, currently inside .LIB " << libInside;
  }
}

//--------------------------------------------------------------------------
// Function      : readLine
// Purpose       : Line-terminator-agnostic istream::getline() workalike for
//               : reading a single line from input stream
//
// Special Notes : file ptr is advanced
//               : line terminator is extracted
//--------------------------------------------------------------------------
void readLine( std::istream & in, std::string& line )
{
  line.clear();
  char c;

  // read all characters into line
  while( in.good() )
  {
    // read a char and append to line
    in.get( c );

    if ( in.eof() || c == '\n' )
    {
      return;
    }
    else if( c == '\r' )
    {
      if( in.peek() == '\n' )
      {
        // extract the CR LF pair and discard
        in.get();
      }

      return;
    }
    else
    {
      line.push_back( c );
    }
  }
}

//----------------------------------------------------------------------------
// Function       : splitLine
// Purpose        : Split an input character string into a tokenized vector.
// Special Notes  :
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 10/15/2014
//----------------------------------------------------------------------------
void splitLine( const std::string& charLine, 
                TokenVector& line )
{
  int lineLength = charLine.length();
  int currPtr = 0;
  char c=charLine[ currPtr ];
  const std::string nonid(" \t\n\r\0");
  line.clear();

  while ( currPtr < lineLength )
  {
    StringToken field;
    field.string_.reserve(16);

    if (nonid.find(c) == nonid.npos) //not a whitespace character
    {
      field.string_ += c;
      c = charLine[ ++currPtr ];
      while ( currPtr < lineLength )
      {
        if (nonid.find(c) == nonid.npos)
        {
          field.string_ += c;
          c = charLine[ ++currPtr ];
        }
        else
        {
          break;
        }
      }
    }
    else
    {
      if (c == '\n' || c == '\r' || c == '\0')
      {
        field.string_ += c;
      }
      else
      {
        field.string_ += c;
        c = charLine[ ++currPtr ];
      }
    }

    if (field.string_.length() > 0)
    {
      line.push_back(field);
    }
  }
}

//----------------------------------------------------------------------------
// Function       : splitLineNoWS
// Purpose        : Split an input character string into a tokenized vector.
// Special Notes  : No whitespace will be in the TokenVector
// Scope          : public
// Creator        : Heidi Thornquist
// Creation Date  : 10/15/2014
//----------------------------------------------------------------------------
void splitLineNoWS( const std::string& charLine, 
                    TokenVector& line )
{
  int lineLength = charLine.length();
  int currPtr = 0;
  char c=charLine[ currPtr ];
  const std::string nonid(" \t\n\r\0");
  line.clear();

  while ( currPtr < lineLength )
  {
    StringToken field;
    field.string_.reserve(16);

    if (nonid.find(c) == nonid.npos) //not a whitespace character
    {
      field.string_ += c;
      c = charLine[ ++currPtr ];
      while ( currPtr < lineLength )
      {
        if (nonid.find(c) == nonid.npos)
        {
          field.string_ += c;
          c = charLine[ ++currPtr ];
        }
        else
        {
          break;
        }
      }
    }
    else
    {
      if (c != '\n' && c != '\r' && c != '\0')
      {
        c = charLine[ ++currPtr ];
      }
    }

    if (field.string_.length() > 0)
    {
      line.push_back(field);
    }
  }
}

//--------------------------------------------------------------------------
// Function      : removeTwoTerminalDevice
// Purpose       : Given a two terminal device, this function checks to see
//                 if both nodes on the device are the same and, if so,
//                 decides whether or not the device should be removed
//                 from the circuit.  The decision to remove is based upon
//                 whether the specific device type was specified during the
//                 preprocessing phase as a device for which redundancies
//                 should be removed.  E.g., if the lines
//
//                 .PREPROCESS REMOVEUNUSED C
//                 C1 1 1 1
//                 R1 2 2 1
//
//                 appear in the netlist, the capacitor C1 will be removed
//                 from the netlist, whereas the resistor R1 will not.
//
// Special Notes : This function only determines *whether* a specific device
//                 should be removed; the process of removing the device from
//                 the netlist takes place in the handleLinePass1 and
//                 getLinePass2 functions.
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/08/2007
//--------------------------------------------------------------------------
bool removeTwoTerminalDevice(const std::vector<bool>& pFilter,
                             const char& linetype,
                             const ExtendedString & node1, 
                             const ExtendedString & node2)
{
  bool result=false;

  if ( ( pFilter[PreprocessType::REDUNDANT_C] && linetype == 'C') 
     || ( pFilter[PreprocessType::REDUNDANT_D] && linetype == 'D') 
     || ( pFilter[PreprocessType::REDUNDANT_I] && linetype == 'I') 
     || ( pFilter[PreprocessType::REDUNDANT_L] && linetype == 'L') 
     || ( pFilter[PreprocessType::REDUNDANT_R] && linetype == 'R')
     || ( pFilter[PreprocessType::REDUNDANT_V] && linetype == 'V') )
  {
    if (node1 == node2)
    {
      result=true;
    }
  }

  return result;
}

//--------------------------------------------------------------------------
// Function      : removeThreeTerminalDevice
// Purpose       : Given a three terminal device, this function checks to see
//                 if all three nodes on the device are the same and, if so,
//                 decides whether or not the device should be removed
//                 from the circuit.  The decision to remove is based upon
//                 whether the specific device type was specified during the
//                 preprocessing phase as a device for which redundancies
//                 should be removed.  E.g., if the lines
//
//                 .PREPROCESS REMOVEUNUSED M
//                 M1 1 1 1 4 NMOS
//                 Q1 2 2 2 5 NPN
//
//                 appear in the netlist, the MOSFET M1 will be removed
//                 from the netlist, whereas the BJT Q1 will not.
//
// Special Notes : This function only determines *whether* a specific device
//                 should be removed; the process of removing the device from
//                 the netlist takes place in the handleLinePass1 and
//                 getLinePass2 functions.
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool removeThreeTerminalDevice(const std::vector<bool>& pFilter,
                               const char& linetype,
                               const ExtendedString & node1,
                               const ExtendedString & node2,
                               const ExtendedString & node3)
{
  bool result=false;

  // Check three terminal devices.
  if ( (pFilter[PreprocessType::REDUNDANT_M] && linetype == 'M') 
        || (pFilter[PreprocessType::REDUNDANT_Q] && linetype == 'Q') )
  {
    if (node1 == node2 && node2 == node3)
    { 
      result=true;
    }
  }

  return result;
}


//--------------------------------------------------------------------------
// Function      : removeDevice
// Purpose       : Given a device block, this function checks to see
//                 if all the nodes on the device are the same and, if so,
//                 decides whether or not the device should be removed
//                 from the circuit.  The decision to remove is based upon
//                 whether the specific device type was specified during the
//                 preprocessing phase as a device for which redundancies
//                 should be removed.  E.g., if the lines
//
//                 .PREPROCESS REMOVEUNUSED M
//                 M1 1 1 1 4 NMOS
//                 Q1 2 2 2 5 NPN
//
//                 appear in the netlist, the MOSFET M1 will be removed
//                 from the netlist, whereas the BJT Q1 will not.
//
// Creator       : Heidi Thornquist
//
// Creation Date : 11/9/2015
//--------------------------------------------------------------------------
bool removeDevice(const std::vector<bool>& pFilter,
                  const DeviceBlock& device)
{ 
  bool result=false;
 
  // Get the vector of nodes from the device
  const std::vector<std::string> & nodes = device.getNodeValues();
  std::string linetype = device.getNetlistDeviceType();
  
  // Check two terminal devices.
  if ( ( pFilter[PreprocessType::REDUNDANT_C] && linetype == "C")
     || ( pFilter[PreprocessType::REDUNDANT_D] && linetype == "D")
     || ( pFilter[PreprocessType::REDUNDANT_I] && linetype == "I") 
     || ( pFilter[PreprocessType::REDUNDANT_L] && linetype == "L") 
     || ( pFilter[PreprocessType::REDUNDANT_R] && linetype == "R")
     || ( pFilter[PreprocessType::REDUNDANT_V] && linetype == "V") )
  { 
    if (nodes[0] == nodes[1])
    { 
      result=true;
    }
  }

  // Check three terminal devices.
  if ( (pFilter[PreprocessType::REDUNDANT_M] && linetype == "M") 
        || (pFilter[PreprocessType::REDUNDANT_Q] && linetype == "Q") )
  {
    if (nodes[0] == nodes[1] && nodes[1] == nodes[2])
    { 
      result=true;
    }
  }

  return result;
}

//---------------------------------------------------------------------------
// Function      : readExternalParamsFromFile
// Purpose       :
// Special Notes : Used to read parameters "tag" = "value" from an external
//                 file.  Any "tag"s found while reading the netlist during
//                 parsing will be replaced by "value".
// Scope         : public
// Creator       : Richard Schiek, 1437 Electrical Systems Modeling
// Creation Date : 07/24/2012
//---------------------------------------------------------------------------
void readExternalParamsFromFile( Parallel::Communicator& comm,
  std::string filename,
  std::vector< std::pair< std::string, std::string > > & paramList )
{
  // at this stage just support the Dakota params.in format of "value" = "tag".
  // we could support other formats as well.

  bool validFile = Util::checkIfValidFile(filename);

  // These parameters are read in on one processor and communicated to the others.
  if (comm.procID() == 0 && validFile)
  {
    const std::string allowedChars("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890+-_.$");
    const std::string whiteSpace(" \t=\n\r");    // note we will treat "=" as whitespace
    const std::string commentChars("*#;");  // if we find any of these then the rest of the line is a comment

    // open the params file and read its content
    std::ifstream paramFile(filename.c_str(), std::ios::in);
    if( paramFile.good() )
    {
      // loop over the file trying to gather names / values and response functions required.
      // general format:
      // white space tag/value white space or = tag/value
      //
      std::string aLine;
      getline(paramFile, aLine);
      bool word1IsNumeric = false;
      while ( paramFile.good() )
      {
        // getline worked, so try to parse the line
        // procedure (1) discard any comments
        //           (2) find tag/value boundaries by looking for allowedChars followed by whitespace
        //           (3) sort out order of tag/value pair
        //           (4) store pair in paramList
        //           (5) try to get another line

        // don't bother with lines that are blank or essentially blank (i.e. just \r\n)
        // shortest line we could have is x=1 or 3 chars.
        if( aLine.length() > 2 )
        {
          std::string::size_type commentLoc = aLine.find_first_of( commentChars, 0 );
          if( commentLoc != std::string::npos )
          {
            // chop off comment to end of line.
            aLine.erase( commentLoc, aLine.length()-commentLoc );
          }
          //check overall line length again.  This could have been just a comment line
          if( aLine.length() > 2 )
          {
            std::string::size_type word1Start = aLine.find_first_of( allowedChars, 0 );
            // check that we found a valid word1Start otherwise stop trying to parse this line.
            if( word1Start != std::string::npos )
            {
              std::string::size_type word1End   = aLine.find_first_of( whiteSpace, word1Start );
              // check that we found a valid word1End otherwise stop trying to parse this line.
              if( word1End != std::string::npos )
              {
                std::string::size_type word2Start = aLine.find_first_of( allowedChars, word1End );
                // check that we found a valid word2Start otherwise stop trying to parse this line.
                if( word2Start != std::string::npos )
                {
                  std::string::size_type word2End   = aLine.find_first_of( whiteSpace, word2Start );
                  // check that we found a valid word2End
                  if( word2End == std::string::npos )
                    word2End = aLine.length();
                  // if we get here then we have valid start,end indicies for word1 and word2
 
                  std::string word1=aLine.substr( word1Start, (word1End - word1Start) );
                  std::string word2=aLine.substr( word2Start, (word2End - word2Start) );

                  // sort out tag/value ordering.
                  // if word1=number assume format is value = tag
                  // otherwise assume format is tag = value
                  std::stringstream converter;
                  converter << word1;
                  double testvalue;
                  converter >> testvalue;
                  if( converter.fail() && !word1IsNumeric )
                  {
                    // couldn't convert tag1 to a double value so assume format is word1=tag word2=value
                    paramList.push_back( std::pair<std::string,std::string>(word1,word2) );
                  }
                  else
                  {
                    // tag1 was successfully converted to a value so assume format is word1=value word2=tag
                    paramList.push_back( std::pair<std::string,std::string>(word2,word1) );
                    word1IsNumeric=true;
                  }

                }  // if ( word2Start != std::string::npos )
              }    // if( word1End != std::string::npos )
            }      // if( word1Start != std::string::npos )
          }        // if( aLine.length() > 2 ) -- check after removing comments
        }          // if( aLine.length() > 2 ) -- outer check
        // try to get another line
        getline(paramFile, aLine);
      }
    }
    else
    {
      validFile = false;
    }
  }

  // Generate a warning if the user-specified params file does not exist, cannot
  // be opened, or is a directory name rather than a file name.  See SON Bugs 730
  // and 785 for more details.
  if ( comm.procID()==0 && !validFile )
  {
    Report::UserWarning() << "Could not find parameter file: " + filename + "  Attempting to continue.";
  }

  // Broadcast the parameters to all the other processors in parallel.
  if (comm.numProc() > 1)
  {
    int numParams = paramList.size();
    comm.bcast( &numParams, 1, 0 );

    for (int i=0; i<numParams; i++)
    {
      int* len = new int[2];

      // Processor zero packs the param.
      char *paramFirst = 0, *paramSecond = 0;
      if (comm.procID() == 0)
      {
        len[0] = paramList[i].first.length();
        len[1] = paramList[i].second.length();
        paramFirst = new char[len[0]+1];
        std::strcpy( paramFirst, paramList[i].first.c_str() );
        paramSecond = new char[len[1]+1];
        std::strcpy( paramSecond, paramList[i].second.c_str() );
      }

      // All other processors receive param lengths
      comm.bcast(len,2,0);
      if (paramFirst == 0)
      {
        paramFirst = new char[len[0]+1];
        paramSecond = new char[len[1]+1];
      }
      comm.bcast(paramFirst,len[0]+1,0);
      comm.bcast(paramSecond,len[1]+1,0);

      // Create pair.
      if (comm.procID() != 0)
      {
        std::string pF( paramFirst );
        std::string pS( paramSecond );
        paramList.push_back( std::pair<std::string,std::string>(pF,pS) );
      }

      delete [] len;
      delete [] paramFirst;
      delete [] paramSecond;
    }    
  }

  // for debug purposes.  output the params as read
  dout() << "Parameters read from \"" << filename << "\"" << std::endl;
  std::vector< std::pair< std::string, std::string > >::iterator listitr = paramList.begin();
  std::vector< std::pair< std::string, std::string > >::iterator enditr = paramList.end();
  while( listitr != enditr )
  {
    dout() << "  " << listitr->first << " , " << listitr->second << std::endl;
    listitr++;
  }
}

//-----------------------------------------------------------------------------
// Function      : combineParamValueFields
// Purpose       :
// Special Notes : 
//
// This function  sets up the combinedString, which is the concatenation of
// one or more fields/tokens from the parsed line.  Expressions that were 
// specified inside of curly braces or single quotes have already 
// been put into a single field string, but any expressions that don't 
// have these delimeters are potentially in a bunch of separate 
// fields.  To parse these properly, they need to be combined 
// into a single string.  
//
// This is to support issue 195.  It was prompted by this use 
// case, which happens a lot in industrial netlists:  
//
// .param PA=AGAUSS(1,1,1)
//
// Note that there is also an attempt to catch expressions without 
// curly braces earlier in the parsing process, in the 
// SpiceSeparatedFieldTool::getLine function.  That was added to address 
// bug 1692 on bugzilla.  However, that function is basically a lexer, 
// and doesn't have enough context to handle a lot of use cases.  
// As such, it does not attempt to operate on "." lines like .PARAM.
//
// When this function is called linePosition is at the first "value" 
// token for a parameter.  The value will ultimately be comprised of one 
// or more string tokens starting with this one.
//
// If there is more than one parameter on this line, the only 
// reliable delimeter between parameters is the equal sign.  
// Commas are not reliable, as they can be used inside of expressions.
//
// Creator       : Eric Keiter
// Creation Date : 04/18/2021
//-----------------------------------------------------------------------------
void combineParamValueFields(
    const TokenVector & parsed_line, 
    int & linePosition, 
    std::string & combinedString)
{
  const int numFields = parsed_line.size();

  combinedString = parsed_line[ linePosition ].string_ ;

  // search for the next equal signs.  If one is found, then
  // that means that there is more than one parameter defined
  // on this line.  The next equals sign is associated with 
  // the next parameter, rather than the current one being processed.
  bool foundEquals=false; bool foundComma=false; int equalsIndex=linePosition;
  for (;equalsIndex<numFields;equalsIndex++)
  {
    if ( parsed_line[ equalsIndex ].string_ == "=")
    { 
      // check if this is a comparison operator "==",">=", "<=" or "!="
      bool compOp=false;
      if (equalsIndex > 0 && !(parsed_line[ equalsIndex-1 ].string_.empty())  )
      {
        if (parsed_line[ equalsIndex-1 ].string_.back() == '<' || 
            parsed_line[ equalsIndex-1 ].string_.back() == '>' || 
            parsed_line[ equalsIndex-1 ].string_.back() == '!' || 
            parsed_line[ equalsIndex-1 ].string_.back() == '=') { compOp=true; }
      }
      if (!compOp)
      {
        if (equalsIndex < numFields-1 &&  !(parsed_line[ equalsIndex+1 ].string_.empty()))
        { 
          if (parsed_line[ equalsIndex+1 ].string_.front() == '=') {compOp=true; }
        }
      }
      if (!compOp) { foundEquals=true; break; }
    }
  }
  if (foundEquals) { if ( parsed_line[ equalsIndex-2 ].string_ == ",") { foundComma=true; } }

  if( (linePosition < numFields-1) )
  {
    int end = parsed_line.size();
    if(foundEquals)
    {
      if (foundComma) { end = equalsIndex-2; }
      else { end = equalsIndex-1; }
    }
    
    for (int lp1 = linePosition+1;lp1<end;lp1++) { combinedString += parsed_line[lp1].string_; }

    //if (modified) combinedString = std::string("{" + combinedString + "}");
  }
  combinedString = std::string("{" + combinedString + "}"); // always do this for now.  

  if (foundEquals) linePosition=equalsIndex-1;
  else linePosition = parsed_line.size();

}

//-----------------------------------------------------------------------------
// Function      : extractParamData
// Purpose       : Extract the parameters from a netlist .PARAM line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool extractParamData(
  CircuitBlock &            circuit_block,
  const std::string &       netlist_filename,
  const TokenVector &       parsed_line)
{
  const int numFields = parsed_line.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 4 )
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".param line has an unexpected number of fields";
    return false;
  }

  Util::OptionBlock option_block("PARAM", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);
 
  int linePosition = 1;   // Start of parameters on .param line.

  Util::Param parameter("", "");
  while ( linePosition < numFields )
  {
    parameter.setTag( parsed_line[linePosition].string_ );

    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") || (parameter.uTag() == "FREQ"))
    {
      Report::UserError0().at(netlist_filename, parsed_line[linePosition].lineNumber_)
        << "Parameter name " << parameter.uTag() << " is not permitted";
    }

    if ( linePosition + 2 >= numFields )
    {
      Report::UserError0().at(netlist_filename, parsed_line[linePosition].lineNumber_)
        << "Unexpectedly reached end of line while looking for parameters in .PARAM statement";
      return false;
    }

    if ( parsed_line[ linePosition+1].string_ != "=" )
    {
      Report::UserError0().at(netlist_filename, parsed_line[linePosition + 1].lineNumber_)
        << "Equal sign (=) required between parameter and value in .PARAM statement";
    }

    linePosition += 2;   // Advance to parameter value field.

    std::string combinedString;
    combineParamValueFields(parsed_line,linePosition,combinedString); // this will update the linePosition

    if (combinedString[0] == '{')
    {
      Util::Expression expPtr(circuit_block.getExpressionGroup(), combinedString);

      std::string msg;
      if (!(expPtr.getVoltageNodes().empty()) )
        msg = "Node Voltage";
      if (!(expPtr.getDeviceCurrents().empty()) )
        msg = "Device Current";
      if (!(expPtr.getLeadCurrents().empty()) )
        msg = "Lead Current";
      if (msg != "")
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << msg << " may not be used in parameter expression (" << parameter.tag() << ")";
      }
      parameter.setVal(combinedString);
    }
    else
    {
      // if we get here, that means that the 
      // (1) the user did not apply delimiters  (curly braces or single quotes)
      // (2) the above code did not add any delimeters.
      ExtendedString tmp ( combinedString );
      if (tmp.possibleParam())
        parameter.setVal(std::string("{" + combinedString + "}")); // if this might be a parameter, and we don't already have "{", add them
      else
        parameter.setVal(combinedString); // if we get here, this should be a pure number, and we don't want to involve the expr library.  ERK.  Note, check that this line ever happens now.  We may be at a stage where curly braces are always added.
    }

    option_block.addParam( parameter );
  }

  for (Util::ParamList::const_iterator paramIter = option_block.begin(), paramEnd = option_block.end(); paramIter != paramEnd; ++paramIter)
  {
    ExtendedString P ( (*paramIter).tag() );
    P.toUpper();
    if (!P.possibleParam())
    {
      Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
        << "Illegal parameter name: " << P;
    }
  }

  circuit_block.addParams(option_block);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : extractGlobalParamData
// Purpose       : Extract the parameters from a netlist .GLOBAL_PARAM line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/05/2001
//-----------------------------------------------------------------------------
bool extractGlobalParamData(
  CircuitBlock &                circuit_block,
  const std::string &           netlist_filename,
  const TokenVector &           parsed_line)
{
  const int numFields = parsed_line.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 4 )
  { 
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << ".param line has an unexpected number of fields";
    return false;
  }

  if (circuit_block.isSubcircuit())
  { 
    Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Attempt to assign global_param inside of subcircuit";
    return false;
  }

  Util::OptionBlock option_block("GLOBAL_PARAM", Util::OptionBlock::ALLOW_EXPRESSIONS,
                                  netlist_filename, parsed_line[0].lineNumber_);

  int linePosition = 1;   // Start of parameters on .param line.

  Util::Param parameter("", "");
  while ( linePosition < numFields )
  { 
    parameter.setTag( parsed_line[linePosition].string_ );
    
    if ( (parameter.uTag() == "TEMP") || (parameter.uTag() == "VT") ||
         (parameter.uTag() == "GMIN") || (parameter.uTag() == "TIME") || (parameter.uTag() == "FREQ"))
    { 
      Report::UserError0().at(netlist_filename, parsed_line[linePosition].lineNumber_)
        << "Parameter name " << parameter.uTag() << " is not permitted";
    }
    
    if ( linePosition + 2 >= numFields )
    { 
      Report::UserError0().at(netlist_filename, parsed_line[linePosition].lineNumber_) 
        << "Unexpectedly reached end of line while looking for parameters in .PARAM statement";
      return false;
    }
    
    if ( parsed_line[ linePosition+1].string_ != "=" )
    { 
      Report::UserError0().at(netlist_filename, parsed_line[linePosition + 1].lineNumber_)
        << "Equal sign (=) required between parameter and value in .PARAM statement";
    }
    
    linePosition += 2;   // Advance to parameter value field.

    std::string combinedString;
    combineParamValueFields(parsed_line,linePosition,combinedString); // this will update the linePosition

    if (combinedString[0] == '{')
    {
      Util::Expression expPtr(circuit_block.getExpressionGroup(), combinedString);

      std::string msg;
      if (!(expPtr.getVoltageNodes().empty()) )
        msg = "Node Voltage";
      if (!(expPtr.getDeviceCurrents().empty()) )
        msg = "Device Current";
      if (!(expPtr.getLeadCurrents().empty()) )
        msg = "Lead Current";
      if (msg != "")
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << msg << " may not be used in parameter expression (" << parameter.tag() << ")";
      }
      parameter.setVal(combinedString);
    }
    else
    {
      // if we get here, that means that the 
      // (1) the user did not apply delimiters  (curly braces or single quotes)
      // (2) the above code did not add any delimeters.
      ExtendedString tmp ( combinedString );
      if (tmp.possibleParam())
        parameter.setVal(std::string("{" + combinedString + "}")); // if this might be a parameter, and we don't already have "{", add them
      else
        parameter.setVal(combinedString); // if we get here, this should be a pure number, and we don't want to involve the expr library.  ERK.  Note, check that this line ever happens now.  We may be at a stage where curly braces are always added.
    }

    option_block.addParam( parameter );
  }

  Util::ParamList::const_iterator paramIter = option_block.begin(), paramEnd = option_block.end();
  for ( ; paramIter != paramEnd; ++paramIter)
  {
    ExtendedString P ( (*paramIter).tag() );
    P.toUpper();
    if (!P.possibleParam())
    {
      Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
        << "Illegal parameter name: " << P;
    }
  }

  circuit_block.addGlobalParams(option_block);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : extractOperatorData
// Purpose       : Extract data for operators like V(), I(), N(), P(), DNO(),
//                 S(), etc, and place them into the option_block as
//                 Util:Param objects
// Special Notes : This function return the position variable as the next
//                 index, in parsed_line, beyond the closing parenthesis of
//                 the operator.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 9/09/2021
//-----------------------------------------------------------------------------
bool extractOperatorData(const TokenVector &  parsed_line,
                         int&                 position,
                         Util::OptionBlock&   option_block,
                         std::ostringstream&  msg,
                         int&                 p_err)
{
  const int numFields = parsed_line.size();
  ExtendedString field("");

  // Note: the next if statement has to support I, IR, II, IM, IP and IDB.  That is why it
  // checks for strings of size 2 and 3.  A similar comment applies to the block for
  // V below.
  if (toupper(parsed_line[position].string_[0]) == 'I' && parsed_line[position].string_.size() <= 3)
  {
    if ((position+3 < numFields && parsed_line[position+3].string_ == ")") ||
        (position+4 < numFields && parsed_line[position+4].string_ == ")"))
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,1.0));

      if (parsed_line[position+3].string_ == ")")
      {
        field = parsed_line[position+2].string_;
        field.toUpper();
        option_block.addParam( Util::Param(field, 0.0) );

        position += 4;
      }
      else
      {
        // Note:  This block is here to handle the case where a user has
        // asked for I(YSOMETHING NAME) instead of I(YSOMETHING!NAME)
        field = parsed_line[position+2].string_ + " " + parsed_line[position+3].string_;
        field.toUpper();
        option_block.addParam( Util::Param(field,0.0) );

        position += 5;
      }
    }
    else
    {
      msg << "Unrecognized current specification";
      p_err = position;
    }
  }
  else if ( toupper(parsed_line[position].string_[0]) == 'V' && parsed_line[position].string_.size() <= 3 )
  {
    // position+3 < numFields test required to prevent a core dump on 
    // something invalid like .PRINT TRAN V(1
    if (position+3 < numFields && parsed_line[position+3].string_ == ")")
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,1.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      position += 4;
    }
    else if (position+5 < numFields && parsed_line[position+5].string_ == ")")
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,2.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      field = parsed_line[position+4].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

       position += 6;
    }
    else
    {
      msg << "Unrecognized voltage specification";
      p_err = position;
    }
  }
  else if ( toupper(parsed_line[position].string_[0]) == 'N'
            && parsed_line[position].string_.size() == 1 )
  {
    if( position+3 < numFields && parsed_line[position+3].string_ == ")" )
    {
      option_block.addParam( Util::Param("N",1.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      position += 4;
    }
    else
    {
      msg << "Unrecognized parenthetical specification";
      p_err = position;
    }
  }
  else if( (parsed_line[position].string_.size() == 3) &&
           (parsed_line[position].string_[0] == 'D' || parsed_line[position].string_[0] == 'd') )
  {
    if( position+3 < numFields && parsed_line[position+3].string_ == ")" )
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,1.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      position += 4;
    }
    else if (position+5 < numFields && parsed_line[position+5].string_ == ")")
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,2.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      field = parsed_line[position+4].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      position += 6;
    }
    else
    {
      msg << "Unrecognized noise specification";
      p_err = position;
    }
   }
  else if( parsed_line[position].string_ == "W" || parsed_line[position].string_ == "w" ||
          parsed_line[position].string_ == "P" || parsed_line[position].string_ == "p")
  {
    if( position+3 < numFields && parsed_line[position+3].string_ == ")" )
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,1.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      position += 4;
    }
    else
    {
      msg << "Unrecognized power specification";
      p_err = position;
    }
  }
  else if ( ((toupper(parsed_line[position].string_[0]) == 'S') ||
             (toupper(parsed_line[position].string_[0]) == 'Y') ||
             (toupper(parsed_line[position].string_[0]) == 'Z'))
             && parsed_line[position].string_.size() <= 3 )
  {
    if (position+5 < numFields && parsed_line[position+5].string_ == ")")
    {
      field = parsed_line[position].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,2.0) );

      field = parsed_line[position+2].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      field = parsed_line[position+4].string_;
      field.toUpper();
      option_block.addParam( Util::Param(field,0.0) );

      position += 6;
   }
   else
   {
      msg << "Unrecognized S-parameter specification";
      p_err = position;
   }
  }
  else
  {
    msg << "Unrecognized parenthetical specification";
    p_err = position;
  }

  return msg.str().empty();
}

} // namespace IO
} // namespace Xyce

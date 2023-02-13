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
// Purpose        : Provide some string and ParamList manipulating helper
//                  functions
//
// Special Notes  : These are intended for use with the external coupling
//                  output interface, but may have more general utility,
//                  so I'm putting them here.
//
// Creator        : Tom Russo
// Creation Date  : 02/07/2018
//-----------------------------------------------------------------------------


#include <Xyce_config.h>
#include <string>
#include <vector>
#include <cstddef>
#include <N_UTL_Param.h>
#include <N_UTL_ExtendedString.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_OutputAPIHelpers.h>

namespace Xyce {
namespace Util {
//-----------------------------------------------------------------------------
// Function      : Xyce::Util::stripWhiteSpace
// Purpose       : strip leading and trailing white space from a string
// Special Notes :
// Scope         : free function
// Creator       : Tom Russo
// Creation Date : 02/07/2018
//-----------------------------------------------------------------------------
/// Given a string, strip leading and trailng white space
///
/// @param s  string to strip.  Note that we modify this internally, and
///           are therefore counting on the fact that we are NOT passing by
///           reference.
/// @return modified string
std::string stripWhiteSpace(std::string s)
{
  std::size_t firstNonWhite=s.find_first_not_of(" \t\f\v\n\r");
  if (firstNonWhite != std::string::npos)
    s.erase(0,firstNonWhite);
  std::size_t lastNonWhite=s.find_last_not_of(" \t\f\v\n\r");
  if (lastNonWhite != std::string::npos)
    s.erase(lastNonWhite+1);

  return s;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::parseFunctionString
// Purpose       : try to recognize strings as function calls
// Special Notes :
// Scope         : free function
// Creator       : Tom Russo
// Creation Date : 02/07/2018
//-----------------------------------------------------------------------------
/// Given a output variable specifier, check if it looks like a function call
///
/// If the string passed in looks like <function>(arg[,arg]*) then we
/// return the function name and vector of arguments.
///
/// @param[in]  s       Reference to containing output specifier
/// @param[out] funName Reference to string into which to place recognized access function name
/// @param[out] funArgs Reference to string vector into which to place function arguments
/// @return true if recognized, false otherwise
bool parseFunctionString(const std::string &s, std::string &funName,
                         std::vector<std::string> &funArgs)
{
  bool success=false;
  
  funArgs.clear();
  funName="";
  
  std::size_t firstLeftParen=s.find_first_of("(");
  std::size_t lastLeftParen=s.find_last_of("(");
  // There Can Be Only One!
  if (firstLeftParen != std::string::npos && firstLeftParen == lastLeftParen)
  {
    // And we MUST have the final character be ')'
    if (s[s.length()-1] == ')')
    {
      std::string funString=s.substr(0,firstLeftParen);
      toUpper(funString);
      if (!funString.empty())
      {
        std::string argString=s.substr(firstLeftParen+1);
        argString.erase(argString.length()-1);
        argString = stripWhiteSpace(argString);
        toUpper(argString);
        if (!argString.empty())
        {
          funName=funString;
          success=true;
          while (!argString.empty())
          {
            std::size_t firstComma=argString.find_first_of(",");
            // If there's only one argument, we're done here
            if (firstComma == std::string::npos)
            {
              funArgs.push_back(argString);
              argString.clear();
            }
            else
            {
              std::string firstArg=argString.substr(0,firstComma);
              argString.erase(0,firstComma+1);
              funArgs.push_back(firstArg);
            }
          }
        }
      }
    }
  }
  return success;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Util::stringToParamList
// Purpose       : try to parse output variable specifiers into Params
// Special Notes :
// Scope         : free function
// Creator       : Tom Russo
// Creation Date : 02/07/2018
//-----------------------------------------------------------------------------
/// Given a output variable specifier, try to parse it into a ParamList
///
/// Tries to parse a output variable and turn it into a list of Param
/// objects.  If we can recognize it, push_back the parameters into the given
/// ParamList and return true.  Otherwise return false.
///
/// Recognizes:  V, I, N, P, W constructs (including the special ones
///               like VI, II, VR, VI, etc.) if and ONLY if they have the
///               right length names (1-3 for I and V, 1 for P, N, W) and
///               the right number of arguments (1-2 for V*, 1 for everything
///               else.
///
///              Brace-delimited expressions
///              Bare words (presumed to be parameter names) if they
///                contain no braces or parens.
///
/// @param[in] s          String reference containing output variable
/// @param[out] paramList References to ParamList object.
bool stringToParamList(const std::string & s, ParamList &paramList)
{
  bool success=false;
  std::string workingCopy=stripWhiteSpace(s);

  if (!workingCopy.empty())
  {
    
    // If the first character of the string is { and the last is },
    // we're done here --- this is a brace-delimited expression and will be
    // passed as-is to the Op builder
    if (workingCopy[0] == '{' && workingCopy[workingCopy.length()-1] == '}')
    {
      paramList.push_back(Param(workingCopy,(int)0));
      success=true;
    }
    else
    {
      // See if we're a recognized special operator and parse it
      // We must have ) as the last character in the string, and must
      // start with a recognized function followed by (.
      std::size_t firstLeftParen=workingCopy.find_first_of("(");
      if (firstLeftParen != std::string::npos &&
          workingCopy[workingCopy.length()-1] == ')')
      {
        std::string funName;
        std::vector<std::string> funArgs;
        if (parseFunctionString(workingCopy,funName,funArgs))
        {
          switch(funName[0])
          {
          case 'I':
            // I functions can have 1-3 char names and exactly one
            // argument
            if (funName.length() <= 3 && funArgs.size() == 1)
            {
              paramList.push_back(Param(funName,(int)1));
              paramList.push_back(Param(funArgs[0],(int)0));
              success=true;
            }
            else
              Report::UserError() << "Unrecognized " << funName[0]
                        << " specifier: " << workingCopy
                        << std::endl;
              
            break;
          case 'N':
          case 'W':
          case 'P':
            // N, P and W functions have only 1-char names, and one argument
            if (funName.length() == 1 && funArgs.size() == 1)
            {
              paramList.push_back(Param(funName,(int)1));
              paramList.push_back(Param(funArgs[0],(int)0));
              success=true;
            }
            else
              Report::UserError() << "Unrecognized " << funName[0]
                        << " specifier: " << workingCopy
                        << std::endl;
              
            break;
          case 'V':
            // V functions may have 1-3 character names, and one or two
            // arguments
            if (funName.length() <= 3 && funArgs.size() <= 2)
            {
              paramList.push_back(Param(funName,(int)funArgs.size()));
              for (int i=0;i<funArgs.size();i++)
                paramList.push_back(Param(funArgs[i],(int)0));
              
              success=true;
            }
            else
              Report::UserError() << "Unrecognized " << funName[0]
                        << " specifier: " << workingCopy
                        << std::endl;
            break;
          default:
            Report::UserError() << "Unrecognized function:" << workingCopy << std::endl;
          }
            
        }
      }
      else if (firstLeftParen == std::string::npos)
      {
        // there is no left paren, and we aren't a brace-wrapped expression
        // SO... we must be a bare parameter name, which means we
        // must have no left parens or braces present
        std::size_t firstInvalChar=workingCopy.find_first_of("){},");
        if (firstInvalChar == std::string::npos)
        {
          toUpper(workingCopy);
          paramList.push_back(Param(workingCopy,(int)0));
          success=true;
        }
      }
      else
      {
        Report::UserError() << " I don't know what to do with " << workingCopy
                  << std::endl;
      }
    }
  }
  return success;
}

/// Try to turn a list of output variable names into a ParamList
///
/// Given a vector of strings, try to process it into a ParamList
///
/// This is similar to---but simpler than---what is done in the output manager
/// function "extractPrintData".
///
/// @note This function does NOT clear the param list given to it.  It
///       always just pushes into the back of it.  The caller is responsible
///       for passing in an empty param list.
///
/// @param[in]  stringVec   Vector of output variable strings
/// @param[out] paramList   ParamList into which to push successfully parsed
///                         Params.
/// @param[out] stringValidities  Returned vector of bools representing
///                         whether or not strings could be parsed.  Will
///                         be the same length as stringVec.  An element will
///                         be true if the corresponding string could be parsed.
/// @return true if all strings were parseable, false if any string was rejected.
///
bool stringsToParamList(const std::vector<std::string> & stringVec,
                        ParamList & paramList,
                        std::vector<bool> & stringValidities)
{
  bool retval=true;
  
  stringValidities.clear();
  
  for (std::vector<std::string>::const_iterator it = stringVec.begin();
       it != stringVec.end();
       ++it)
  {
    if (!stringToParamList(*it,paramList))
    {
      stringValidities.push_back(false);
      retval=false;
    }
    else
    {
      stringValidities.push_back(true);
    }
  }
  return retval;
}

}
}

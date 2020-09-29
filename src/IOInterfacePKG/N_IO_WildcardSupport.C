//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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

#include <Xyce_config.h>

#include <N_IO_WildcardSupport.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce{
namespace IO {

//-----------------------------------------------------------------------------
// Function      : IO:makeRegexFromString
// Purpose       : Convert a string (e.g, for a node name or device name) into
//                 a C++11 std:regex object
// Special Notes : This allows us to control what types of regex's are supported
//                 in Xyce.  Only the * and ? characters are supported.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/07/2020
//-----------------------------------------------------------------------------
std::regex makeRegexFromString(const std::string &str)
{
  // These character are "special" to the std::regex class.  They need to be
  // escaped if they exist in str, as does the \ character.
  std::string specialChars="^$.+()[]{}|";

  std::string tmpStr(str);

  if ( (tmpStr.find("*") == std::string::npos) && (tmpStr.find("?") == std::string::npos) )
  {
    Report::DevelFatal().in("makeRegexFromString") <<
      "Xyce wildcard specification must contain at least one * or ? character";
  }

  // escape the \ characters first
  replaceAll(tmpStr, '\\', "\\\\");

  // escape the other special characters
  for (int i=0; i<specialChars.size(); i++)
  {
    std::string escString = "\\";
    escString.push_back(specialChars[i]);
    replaceAll(tmpStr, specialChars[i], escString);
  }

  // turn * and ? into the correct string for std::regex
  replaceAll(tmpStr, '*', "(.*)");
  replaceAll(tmpStr, '?', "(.{1})");
  
  // make regex object
  std::regex e;
  try
  {
    e.assign(tmpStr);
  }
  catch (std::regex_error& regexErr)
  {
    Report::DevelFatal().in("makeRegexFromString") <<
      "Error converting .PRINT wildcard specification " << str << " into std::regex object";
  }

  return e;
}

//-----------------------------------------------------------------------------
// Function      : IO:replaceAll
// Purpose       : Used to replace all instances of a character in a string
//                 within a different string.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/07/2020
//-----------------------------------------------------------------------------
void replaceAll(std::string& str, const char charToFind, const std::string& repStr)
{
  size_t pos = str.find(charToFind);
  while( pos != std::string::npos)
  {
    str.replace(pos, 1, repStr);
    pos = str.find(charToFind, pos+repStr.size());
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : IO:findWildCardMatch
// Purpose       : Is there at least one wild-card match for name in input_set
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/01/2019
//-----------------------------------------------------------------------------
bool findWildCardMatch(
  const std::string &name,
  const unordered_set<std::string> &input_set)
{
  bool ret = true;

  if ( ((name.size() == 1) && (name!="?")) || ((name.find_first_of("*") == std::string::npos) &&
			      (name.find_first_of("?") == std::string::npos)) )
  {
    // The single character ? is a wildcard, as is any multi-character string that has * or ?
    // in it.  So, return false if the string cannot be a valid wildcard in order to save the
    // cost of the regex_match() calls.
    ret = false;
  }
  else
  {
    std::regex e = makeRegexFromString(name);

    // search through the names in input_set to see if any of them match
    // the regex
    unordered_set<std::string>::const_iterator it;
    for (it=input_set.begin(); it!=input_set.end(); ++it)
    {
      if (std::regex_match((*it), e))
        return true;
    }

    if (it == input_set.end())
      ret = false;
  }

  return ret;
}

//-----------------------------------------------------------------------------
// Function      : IO:findAllWildCardMatch
// Purpose       : Find all wild-card matches for name in input_set
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/01/2019
//-----------------------------------------------------------------------------
bool findAllWildCardMatches(
  const std::string &name,
  const unordered_set<std::string> &input_set,
  std::vector<std::string> &matches)
{
  if ( !((name.size() == 1) || ((name.find_first_of("*") == std::string::npos) &&
                                (name.find_first_of("?") == std::string::npos))) )
  {
    std::regex e = makeRegexFromString(name);

    // search through the names in input_set to see if any of them match
    // the wild card pattern
    unordered_set<std::string>::const_iterator it;
    for (it=input_set.begin(); it!=input_set.end(); ++it)
    {
      if (std::regex_match((*it), e))
        matches.push_back(*it);
    }
  }

  return matches.size() != 0 ? true :false;
}

} // namespace IO
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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

namespace Xyce{
namespace IO {

//-----------------------------------------------------------------------------
// Function      : IO:splitWildCardString
// Purpose       : Pull out the non-* fragments of name into nameSubStrings
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/01/2019
//-----------------------------------------------------------------------------
bool splitWildCardString(const std::string& name, std::vector<std::string>& nameSubStrings)
{
  // default is success
  bool ret = true;

  // The stringstream ss recursively extracts each "word" from tempName, where the
  // "words" in tempName have been delimited by spaces.
  std::string tempName(name);
  std::replace(tempName.begin(), tempName.end(), '*', ' ');  // replace '*' by ' '
  std::stringstream ss(tempName);
  std::string temp;
  while (ss >> temp)
    nameSubStrings.push_back(temp);

  if (nameSubStrings.size() == 0)
    ret = false;

  return ret;
}

//-----------------------------------------------------------------------------
// Function      : IO:isWildCardMatch
// Purpose       : Is name a wild-card match for the substrings that
//                 were previouly extracted with splitWildCardString()
// Special Notes : This implements something akin to globbing (but just for
//                 the * character) rather than a regex.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/01/2019
//-----------------------------------------------------------------------------
bool isWildCardMatch(const std::string& name,
                     const std::vector<std::string>& nameSubStrings,
                     int firstStarPos,
                     bool trailingStar)
{
  // default is no match
  bool ret = false;

  // For the case where the wildcard ends with a specific string, rather
  // than a * character, there must be an exact match to that specific
  // string at the end of the name string.
  if (!trailingStar && (nameSubStrings.size() != 0))
  {
    std::string lastSubString = nameSubStrings[nameSubStrings.size()-1];
    if ( (name.size() < lastSubString.size()) ||
         (name.substr(name.size()-lastSubString.size(), lastSubString.size()) != lastSubString) )
      return false;
  }

  int numMatch=0;
  int i = 0;
  // position in node_name (*it) where the current match begins
  std::size_t matchStart=0;
  // position in node_name where the test for the current match should begin
  std::size_t pos=0;
  std::vector<std::string>::const_iterator nss_it;
  for (nss_it=nameSubStrings.begin(); nss_it!=nameSubStrings.end(); ++nss_it, ++i)
  {
    if ((i==0) && (firstStarPos!=0))
    {
      // If the first * is not the first character in name then the beginning of
      // the node_name must be an exact match for the first element in nameSubStrings.
      if (name.substr(0,(*nss_it).size()) == (*nss_it))
        ++numMatch;
    }
    else
    {
      // the "normal case" is that the element of nameSubStrings, being tested now,
      // has to find a match somewhere (starting at offset pos) in node_name,
      matchStart = name.substr(pos).find(*nss_it);
      if (matchStart != std::string::npos)
        ++numMatch;
    }
    if (numMatch != (i+1))
      break;

    // update pos to be the next character after where the matching
    // string (in node_name) ends
    pos = matchStart+(*nss_it).size();
    if (pos >= name.size())
      break;
  }

  if (numMatch == nameSubStrings.size())
    ret = true;

  return ret;
}

//-----------------------------------------------------------------------------
// Function      : IO:findWildCardMatch
// Purpose       : Is there at least one wild-card match for name in input_set
// Special Notes : This implements something akin to globbing (but just for
//                 the * character) rather than a regex.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 10/01/2019
//-----------------------------------------------------------------------------
bool findWildCardMatch(const std::string &name, const unordered_set<std::string> &input_set)
{
  bool ret = true;

  // determine where the first * is.
  std::size_t firstStarPos = name.find_first_of("*");
  bool trailingStar = (name.find_last_of("*") == name.size()-1);

  if ((name.size() == 1) || (firstStarPos == std::string::npos))
  {
    ret = false;
  }
  else
  {
    // Pull out the non-* fragments of name into nameSubStrings.
    std::vector<std::string> nameSubStrings;
    if (!splitWildCardString(name, nameSubStrings)) return false;

    // search through the names in input_set to see if any of them match
    // the wild card pattern
    unordered_set<std::string>::const_iterator it;
    for (it=input_set.begin(); it!=input_set.end(); ++it)
    {
      if (isWildCardMatch((*it), nameSubStrings, firstStarPos, trailingStar))
        return true;
    }

    if (it == input_set.end())
      ret = false;
  }

  return ret;
}

//-----------------------------------------------------------------------------
// Function      : IO:findAllWildCardMatches
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
  // determine where the first * is.
  std::size_t firstStarPos = name.find_first_of("*");
  bool trailingStar = (name.find_last_of("*") == name.size()-1);

  if ( !((name.size() == 1) || (firstStarPos == std::string::npos)) )
  {
    // Pull out the non-* fragments of name into nameSubStrings.
    std::vector<std::string> nameSubStrings;
    if (!splitWildCardString(name, nameSubStrings)) return false;

    // search through the node names to see if any of them match the
    // wild card pattern
    unordered_set<std::string>::const_iterator it;
    for (it=input_set.begin(); it!=input_set.end(); ++it)
    {
      if (isWildCardMatch((*it), nameSubStrings, firstStarPos, trailingStar))
        matches.push_back(*it);
    }
  }

  return matches.size() != 0 ? true :false;
}

} // namespace IO
} // namespace Xyce

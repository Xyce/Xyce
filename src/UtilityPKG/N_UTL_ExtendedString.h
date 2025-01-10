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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_ExtendedString_h
#define Xyce_N_UTL_ExtendedString_h

#include <string>
#include <algorithm>
#include <cstdlib>

#include <iostream>

#include <N_UTL_NoCase.h>

namespace Xyce {
namespace Util {

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to true then 1A=1e-18 rather than 1.
// The default is false.
extern bool useHspiceUnits;

//-----------------------------------------------------------------------------
// Function      : toUpper
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline void toUpper(std::string &s)
{
  for (std::string::iterator it = s.begin(); it != s.end(); ++it)
  {
    *it = toupper(*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : toLower
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline void toLower(std::string &s)
{
  for(std::string::iterator it = s.begin(); it != s.end(); ++it)
  {
    *it = tolower(*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : removeWhiteSpace
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline void removeWhiteSpace(std::string & tmp)
{
  std::string::iterator x = std::remove(tmp.begin(), tmp.end(), ' ');
  tmp.erase(x, tmp.end());
}

//-----------------------------------------------------------------------------
// Function      : removeBadChars
// Purpose       : Remove characters, from a string, that can cause trouble 
//                 during Xyce operator creation.
// Special Notes : The sequences \\ and \" in the string badChars put the
//                 characters backslash and double-quote into that string.
// Creator       : Pete Sholander, SNL
// Creation Date : 7/3/2018
//-----------------------------------------------------------------------------
inline void removeBadChars(std::string & tmp)
{
  std::string badChars("\\{}(),;*:$\"*");
  for (std::string::iterator it = tmp.begin(); it != tmp.end(); ) 
  {
    // erase() advances the iterator, when the "bad character" is found. 
    (badChars.find(*it) != std::string::npos) ? tmp.erase(it) : it++;
  }
}

bool isInt(const std::string & tmpStr);
bool possibleParam(const std::string & tmpStr);

double Value(const std::string & tmp);
bool isValue(const std::string & tmp);

bool isTableFileKeyword(const std::string & tmpStr);

//-----------------------------------------------------------------------------
// Function      : Bval
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline bool Bval(const std::string & tmpStr)
{
  if (isValue(tmpStr))
  {
    return Value(tmpStr) != 0;
  }

  return equal_nocase(tmpStr, "TRUE");
}

//-----------------------------------------------------------------------------
// Function      : isBool
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
inline bool isBool(const std::string & s)
{
  return  equal_nocase(s, "TRUE") || equal_nocase(s, "FALSE") || isValue(s);
}

//-----------------------------------------------------------------------------
// Function      : Ival
// Purpose       :
// Special Notes :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
inline int Ival(const std::string & tmpStr)
{
  if (isInt(tmpStr))
    return atoi(tmpStr.c_str());

  return 0;
}

//-----------------------------------------------------------------------------
// Function      : possibleParam
// Purpose       : Tests if a string is a valid parameter name
// Special Notes :
//
// Examples:
//
//   "RD" is a param
//
//   "{RD}" is NOT a param - anything that would otherwise be a param
//                           becomes NOT a param with curly braces
//
//   "V(1)" is not a param
//   "A+B" is not a param (basically any expression is not a parameter)
//
//   "" (empty string) is not a param (this function gets called with "" a lot)
//   "true"  is not a param
//   "false"  is not a param
//   "12.3" is not a param (isValue() is tested inside of "isBool" and would return true)
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/02/05
// Rewritten     : 04/16/2022, Eric Keiter
//-----------------------------------------------------------------------------
inline bool possibleParam (const std::string &tmpStr)
{
  if (tmpStr.empty()) { return false; }

  std::size_t found = tmpStr.find_first_of("(){}",0); // (redundant with the code below)
  if (found != std::string::npos) { return false; }

  // first character allows $, excludes numbers, excludes period.
  std::string tmp = tmpStr.substr(0,1);
  found = tmp.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_$", 0); 
  if (found != std::string::npos) { return false; }

  // check rest of the characters in tmpStr
  found = tmpStr.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789.", 1);
  if (found != std::string::npos) { return false; }

  //if (isBool(tmpStr)) { return false; } 
  // equivalent of isBool, above.  We should almost never get as far as isValue.
  if (equal_nocase(tmpStr, "TRUE")) { return false; }
  if (equal_nocase(tmpStr, "FALSE")) { return false; }
  if ( isValue(tmpStr) ) { return false; }

  return true;
}

} // namespace Util

//-----------------------------------------------------------------------------
// Class         : ExtendedString
// Purpose       :
//                  This  class  extends the C++ string class by adding
//                  functions toUpper, toLower, and removeWhiteSpace.
//                  Also includes spice style conversion to numeric Value.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class ExtendedString : public std::string
{
  public:
    ExtendedString(const char *X)
      : std::string(X)
    {}

    ExtendedString(const std::string & X)
      : std::string(X)
    {}

    ~ExtendedString()
    {}

    ExtendedString &toUpper() {
      Util::toUpper(*this);

      return *this;
    }

    ExtendedString &toLower() {
      Util::toLower(*this);

      return *this;
    }

    ExtendedString &removeWhiteSpace() {
      Util::removeWhiteSpace(*this);

      return *this;
    }

    ExtendedString &removeBadChars() {
      Util::removeBadChars(*this);

      return *this;
    };

    // This version is with the Xyce scaling factors
    double Value() const {
      return Util::Value(*this);
    }

    // This version is used with the scaling factors defined
    // in the IBIS standard.
    double IBISValue() const {
      return Util::Value(*this);
    }

    bool isValue() const {
      return Util::isValue(*this);
    }

    int Ival() const {
      return Util::Ival(*this);
    }

    bool isInt() const {
      return Util::isInt(*this);
    }

    bool Bval() const {
      return Util::Bval(*this);
    }

    bool isBool() const {
      return Util::isBool(*this);
    }

    bool possibleParam() {
      return Util::possibleParam(*this);
    }
};

} // namespace Xyce

#endif // Xyce_N_UTL_ExtendedString_h

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
// Purpose        : Case insensitive functions and functors
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NoCase_h
#define Xyce_N_UTL_NoCase_h

#include <algorithm>
#include <cctype>
#include <functional>
#include <string>

namespace Xyce {

//-----------------------------------------------------------------------------
// Function      : bit_tolower
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Jan 16 10:00:08 2015
//-----------------------------------------------------------------------------
///
/// Convert to lower case, the old way.
///
/// @invariant
///
/// @param c
///
/// @return
///
///
inline char bit_tolowercorrect(char c)
{
  return (c >= 'A' && c <= 'Z') ? c | 0x20 : c;
}

//-----------------------------------------------------------------------------
// Function      : bit_tolower
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Fri Jan 16 10:00:08 2015
//-----------------------------------------------------------------------------
///
/// Convert to lower case, ignoring original case.  This simply sets the 0x20 bit.
///
/// As a result,
///
/// NUL..US -> space..?
/// @ -> `
/// [ -> {
/// \ -> |
/// ] -> }
/// ^ -> ~
/// _ -> DEL
///
/// but it avoids two calls to tolower.
///
/// @invariant
///
/// Control characters, [, \, ], ^ and _ compare equally to space..?, {, |, |, ~ and DEL.
///
/// @param c
///
/// @return
///
///
inline char bit_tolowerhack(char c)
{
  return (c >= 'A' && c <= 'Z') ? c | 0x20 : c;
}

///
/// Compare two strings, case insensitively.
///
/// Return negative if the first differing character of s0 is less than
/// s1.  Return 0 is they are equal, otherize positive.
///
/// @param s0
/// @param s1
///
/// @return Negative if s0 is less than s1, positive if s0 is greater than s1, zero otherwise
///
int compare_nocase(const char *s0, const char *s1);

///
/// Test if s0 starts with s1, case insensitively.
///
/// Return true if s0 starts with s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return True if s0 starts with s1, case insensitively
///
bool startswith_nocase(const char *s0, const char *s1);

///
/// Test if s0 starts with s1, case insensitively.
///
/// Return true if s0 starts with s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return True if s0 starts with s1, case insensitively
///
inline bool startswith_nocase(const std::string &s0, const std::string &s1)
{
  return startswith_nocase(s0.c_str(), s1.c_str());
}

///
/// Test if fullString ends with ending, case insensitively.
///
/// Return true if fullString ends with ending, case insensitively.
///
/// @param fullString
/// @param ending
///
/// @return True if fullString ends with ending, case insensitively
///
bool endswith_nocase(const char *fullString, const char *ending);

///
/// Test if fullString ends with ending, case insensitively.
///
/// Return true if fullString ends with ending, case insensitively.
///
/// @param fullString
/// @param ending
///
/// @return True if fullString ends with ending, case insensitively
///
inline bool endswith_nocase (const std::string &fullString, const std::string &ending)
{
  return endswith_nocase(fullString.c_str(), ending.c_str());
}

///
/// Test if first string is less than second string, case insensitively.
///
/// Return true if s0 is less than s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return true if s0 is less than s1, case insensitively
///
inline bool less_nocase(const std::string &s0, const std::string &s1)
{
  return compare_nocase(s0.c_str(), s1.c_str()) < 0;
}

///
/// Functor to test if first string is less than second string, case insensitively.
///
/// Return true if s0 is less than s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return true if s0 is less than s1, case insensitively
///
struct LessNoCase 
{
  using result_type = bool;
  using first_argument_type = std::string;
  using second_argument_type = std::string;

  bool operator()(const std::string &s0, const std::string &s1 ) const
  {
    return less_nocase(s0, s1);
  }
};

///
/// Compare two strings for equality, case insensitively.
///
/// Return true if s0 is equal to s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return true if s0 is equal to s1, case insensitively
///
inline bool equal_nocase(const std::string &s0, const std::string &s1)
{
  return compare_nocase(s0.c_str(), s1.c_str()) == 0;
}

///
/// Compare two character for equality, case insensitively.
///
/// Return true if s0 is equal to s1, case insensitively.
///
/// @param c0
/// @param c1
///
/// @return true if s0 is equal to s1, case insensitively
///
inline bool equal_nocasechar(char c0, char c1)
{
  return bit_tolowercorrect(c0) == bit_tolowercorrect(c1);
}

///
/// Functor to hash a case insensitive string.
///
/// Return hash value of string.
///
/// The hashing function was lifted from boost source.
///
/// @param s0
/// @param s1
///
/// @return true if s0 is equal to s1, case insensitively
///
struct HashNoCase
{
  using result_type = size_t;
  using first_argument_type = std::string;

   size_t operator()(const std::string &s) const
  {
    size_t seed = 0;
    for (size_t i = 0, length = s.length(); i < length; ++i)
      seed ^= static_cast<size_t>(bit_tolowercorrect(s[i])) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
  }
};

///
/// Functor to compare two strings for equality, case insensitively.
///
/// Return true if s0 is equal to s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return true if s0 is equal to s1, case insensitively
///
struct EqualNoCase
{
  using result_type = bool;
  using first_argument_type = std::string;
  using second_argument_type = std::string;

  bool operator()(const std::string &s0, const std::string &s1 ) const
  {
    return equal_nocase(s0, s1);
  }
};

///
/// Predicate to compare string for equality, case insensitively.
///
/// Return true if predicate string is equal to s, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return true if predicate string is equal to s, case insensitively
///
struct EqualNoCasePred 
{
  using result_type = bool;
  using first_argument_type = std::string;
  using second_argument_type = std::string;

  EqualNoCasePred(const std::string &s)
    : s_(s)
  {}

  bool operator()(const std::string &s) const
  {
    return equal_nocase(s_, s);
  }

  private:
    const std::string s_;
};

///
/// Test is s0 start with s1, case insensitively.
///
/// Return true if the s0 start with s1, case insensitively.
///
/// @param s0
/// @param s1
///
/// @return True if s0 starts with s1, case insensitively
///
inline bool contains_nocase(const std::string &s0, const std::string &s1)
{
  return std::search(s0.begin(), s0.end(), s1.begin(), s1.end(), equal_nocasechar) != s0.end();
}

} // namepace Xyce

#endif // Xyce_N_UTL_NoCase_h

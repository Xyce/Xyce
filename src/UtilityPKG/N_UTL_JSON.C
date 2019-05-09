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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur, Raytheon
//
// Creation Date  : 8/17/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

#include <Xyce_config.h>

#include <N_UTL_JSON.h>

namespace Xyce {
namespace Util {

JSON::Sep JSON::sep;            ///< JSON separator (,)
JSON::Open JSON::open;          ///< JSON open brace ({)
JSON::Close JSON::close;        ///< JSON close brace (})

namespace {

//-----------------------------------------------------------------------------
// Function      : encode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:25:52 2014
//-----------------------------------------------------------------------------
///
/// Writes character t to output stream encoding quote and control characters.
///
/// @param os   output stream
/// @param c    character
///
///
void
encode(
  std::ostream &        os,
  const unsigned char & c)
{
  if (c == '"')
    os << "\\\"";
  else if (c == '\\')
    os << "\\\\";
  else if (c == '\b')
    os << "\\b";
  else if (c == '\f')
    os << "\\f";
  else if (c == '\n')
    os << "\\n";
  else if (c == '\r')
    os << "\\r";
  else if (c == '\t')
    os << "\\t";
  else if (::iscntrl(c))
    os << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<unsigned short>(c) << std::dec;
  else
    os << c;
}

};

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the character to the JSON stream
///
/// @param jout JSON stream
/// @param t    character
///
/// @return JSON stream
///
///
JSON &operator<<(JSON &jout, const signed char &t) {
  jout.os() << '"';
  encode(jout.os(), t);
  jout.os() << '"';
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the character to the JSON stream
///
/// @param jout JSON stream
/// @param t    character
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const unsigned char &t) {
  jout.os() << '"';
  encode(jout.os(), t);
  jout.os() << '"';
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the character to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const char &t) {
  jout.os() << '"';
  encode(jout.os(), t);
  jout.os() << '"';
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the string to the JSON stream
///
/// @param jout JSON stream
/// @param t    string
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const std::string &s)  {
  jout.os() << '"';
  for (std::string::const_iterator it = s.begin(), end = s.end(); it != end; ++it)
    encode(jout.os(), (*it));
  jout.os() << '"';
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const short &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const unsigned short &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const int &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const unsigned int &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const long &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const unsigned long &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const long long &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the integer to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const unsigned long long &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the float to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const float &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the double to the JSON stream
///
/// @param jout JSON stream
/// @param t    value to write
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const double &t) {
  jout.os() << t;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the separator to the JSON stream
///
/// @param jout JSON stream
/// @param t    unused
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const JSON::Sep &t) {
  jout.os() << ", ";
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the open brace to the JSON stream
///
/// @param jout JSON stream
/// @param t    unused
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const JSON::Open &t) {
  jout.os() << " { ";
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the close brace to the JSON stream
///
/// @param jout JSON stream
/// @param t    unused
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const JSON::Close &t) {
  jout.os() << " } ";
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:27:34 2014
//-----------------------------------------------------------------------------
///
/// Writes the encoded JSON streamopen brace to the JSON stream
///
/// @param jout JSON stream
/// @param t    input JSON stream
///
/// @return JSON stream
///
JSON &operator<<(JSON &jout, const JSON &jin) {
  jout.os() << jin.str();

  return jout;
}

} // namespace Util
} // namespace Xyce

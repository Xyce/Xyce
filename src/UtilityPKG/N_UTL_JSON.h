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

//-------------------------------------------------------------------------
//
// Purpose        :
//
//
//
// Special Notes  :
//
//
// Creator        : David Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_JSON_h
#define Xyce_N_UTL_JSON_h

#include <vector>
#include <list>
#include <map>
#include <set>
#include <sstream>

namespace Xyce {
namespace Util {


//-----------------------------------------------------------------------------
// Class         : JSON
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:33:48 2014
//-----------------------------------------------------------------------------
///
/// JSON is used to encode stuctures and POD into a JSON formatted
/// string.  The class basically serves as a tag to allow the put-to
/// operator (<<) to be implemented for JSON type output.
///
/// The class is no terribly clever and coupld likely use some
/// refactoring, but time and effort required a quick solution.  It's
/// failing is in encoding data that should be read as a composite
/// structure, but is not a structure for writing.  Accomplishing this
/// required the use of Sep, Open and Close to manually encode rather
/// than something clever maintain state.  I.E. can't fake a vector at
/// this point, only structures.  The binding of names to values is also
/// a little wonky.
///
/// The nameValuePair() function handles the binding ala introspection.
///
/// However, its really easy to write encodings for the STL containers
/// and for any regular structure.  The STL containers implemented
/// create JSON objects or vectors encoding.
///
class JSON
{
public:
  struct Sep
  {};

  struct Open
  {};

  struct Close
  {};

  static Sep sep;
  static Open open;
  static Close close;

  JSON()
    : os_(),
      jsonString_()
  {}

private:
  JSON(const JSON &json);
  JSON &operator=(const JSON &json);

public:
  std::ostream &os() {
    return os_;
  }

  const std::string &str() const
  {
    if (jsonString_.empty()) {
      jsonString_ = "{";
      jsonString_ += os_.str();
      jsonString_ += "}";
    }

    return jsonString_;
  }

private:
  std::ostringstream    os_;
  mutable std::string   jsonString_;
};

/// Plain old POD types and std::string
JSON &operator<<(JSON &jout, const signed char &t);
JSON &operator<<(JSON &jout, const unsigned char &t);
JSON &operator<<(JSON &jout, const char &t);
JSON &operator<<(JSON &jout, const short &t);
JSON &operator<<(JSON &jout, const unsigned short &t);
JSON &operator<<(JSON &jout, const int &t);
JSON &operator<<(JSON &jout, const unsigned int &t);
JSON &operator<<(JSON &jout, const long &t);
JSON &operator<<(JSON &jout, const unsigned long &t);
JSON &operator<<(JSON &jout, const long long &t);
JSON &operator<<(JSON &jout, const unsigned long long &t);
JSON &operator<<(JSON &jout, const float &t);
JSON &operator<<(JSON &jout, const double &t);
JSON &operator<<(JSON &jout, const std::string &s);
JSON &operator<<(JSON &jout, const JSON::Sep &t);
JSON &operator<<(JSON &jout, const JSON::Open &t);
JSON &operator<<(JSON &jout, const JSON::Close &t);
JSON &operator<<(JSON &jout, const JSON &jin);


//-----------------------------------------------------------------------------
// Class         : pair_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:44:15 2014
//-----------------------------------------------------------------------------
///
/// Helper class to pass name value pair to JSON output stream.
/// Essentially bind the name to the value.
///
/// @param name         name of object
/// @param value        object
///
template<class V>
struct pair_
{
  pair_(const std::string &name, const V &value)
    : name_(name),
      value_(value)
  {}

  const std::string &   name_;          ///< Object name
  const V &             value_;         ///< Object
};

//-----------------------------------------------------------------------------
// Function      : nameValuePair
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:45:27 2014
//-----------------------------------------------------------------------------
///
/// Creates bound name value pair for use in JSON stream via the
/// operator<<(JSON, pair<T>) operator.
///
/// @param name         Object name
/// @param value        Object
///
/// @return JSON stream
///
template<class V>
inline
const pair_<V>
nameValuePair(const std::string &name, const V &value) {
  return pair_<V>(name, value);
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:46:45 2014
//-----------------------------------------------------------------------------
///
/// Operator to write the pair<T> bound name and value pair to the JSON stream.
///
/// @param jout         JSON stream
/// @param p            bound name and value pair
///
/// @return JSON stream
///
template<class V>
inline
JSON &operator<<(JSON &jout, const pair_<V> &p) {
  jout << p.name_;
  jout.os() << " : ";
  jout << p.value_;
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:47:51 2014
//-----------------------------------------------------------------------------
///
/// Writes std::pair object to the JSON stream.
///
/// @tparam T           Type of std::pair::first
/// @tparam U           Type of std::pair::second
///
/// @param jout         JSON stream
/// @param p            pair
///
/// @return JSON stream
///
template<class T, class U>
JSON &operator<<(JSON &jout, const std::pair<T, U> &p) {
  jout.os() << " { first : ";
  jout << p.first;
  jout.os() << ", second : ";
  jout << p.second;
  jout.os() << " } ";
  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:47:51 2014
//-----------------------------------------------------------------------------
///
/// Writes std::vector object to the JSON streamas a vector.
///
/// @tparam T           Type of std::vector<T, A>
/// @tparam A           Allocator of std::vector<T, A>
///
/// @param jout         JSON stream
/// @param p            vector
///
/// @return JSON stream
///
template <class T, class A>
JSON &operator<<(JSON &jout, const std::vector<T, A> &v)
{
  jout.os() << " [ ";
  for (typename std::vector<T, A>::const_iterator it = v.begin(), end = v.end(); it != end; ++it)
  {
    if (it != v.begin())
      jout.os() << ", ";
    jout << (*it);
  }
  jout.os() << " ] ";

  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:47:51 2014
//-----------------------------------------------------------------------------
///
/// Writes std::list object to the JSON stream as a vector.
///
/// @tparam T           Type of std::list<T, A>
/// @tparam A           Allocator of std::list<T, A>
///
/// @param jout         JSON stream
/// @param p            list
///
/// @return JSON stream
///
template <class T, class A>
JSON &operator<<(JSON &jout, const std::list<T, A> &v)
{
  jout.os() << " [ ";
  for (typename std::list<T, A>::const_iterator it = v.begin(), end = v.end(); it != end; ++it)
  {
    if (it != v.begin())
      jout.os() << ", ";
    jout << (*it);
  }
  jout.os() << " ] ";

  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:47:51 2014
//-----------------------------------------------------------------------------
///
/// Writes std::set object to the JSON streamas a vector.
///
/// @tparam T           Type of std::set<T, C, A>
/// @tparam C           Comparator of std::set<T, C, A>
/// @tparam A           Allocator of std::set<T, C, A>
///
/// @param jout         JSON stream
/// @param p            set
///
/// @return JSON stream
///
template <class T, class C, class A>
JSON &operator<<(JSON &jout, const std::set<T, C, A> &s)
{
  jout.os() << " [ ";
  for (typename std::set<T, C, A>::const_iterator it = s.begin(), end = s.end(); it != end; ++it)
  {
    if (it != s.begin())
      jout.os() << ", ";
    jout << (*it);
  }
  jout.os() << " ] ";

  return jout;
}

//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 09:47:51 2014
//-----------------------------------------------------------------------------
///
/// Writes std::map object to the JSON streamas a vector.
///
/// @tparam K           Key of std::map<K, T, C, A>
/// @tparam T           Type of std::map<K, T, C, A>
/// @tparam C           Comparator of std::map<K, T, C, A>
/// @tparam A           Allocator of std::map<K, T, C, A>
///
/// @param jout         JSON stream
/// @param p            map
///
/// @return JSON stream
///
template <class K, class T, class C, class A>
JSON &operator<<(JSON &jout, const std::map<K, T, C, A> &m)
{
  jout.os() << " { ";
  for (typename std::map<K, T, C, A>::const_iterator it = m.begin(), end = m.end(); it != end; ++it)
  {
    if (it != m.begin())
      jout.os() << ", ";
    jout << nameValuePair((*it).first, (*it).second);
  }
  jout.os() << " } ";

  return jout;
}

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_JSON_h

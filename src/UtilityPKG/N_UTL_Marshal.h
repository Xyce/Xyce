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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 1/7/2014
//
//
//
//
//-------------------------------------------------------------------------

#ifndef N_UTL_MARSHAL_HPP
#define N_UTL_MARSHAL_HPP

#include <stdint.h>
#include <unordered_set>
using std::unordered_set;

#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <list>
#include <set>
#include <map>
#include <typeinfo>

#include <Xyce_config.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Util {

/**
 * @brief Struct <code>Marshal</code> is a data packer for sending and receiving parallel messages.
 * The data put-to (<<) is appended to the stream as a string of bytes, likewise data gotten-from
 * (>>) is extracted from the stream into the object as a string of bytes.
 *
 * The write() and read() functions perform the data movements to and from the packed stream.
 *
 * The common implementation is the create a << and >> operator for an object which properly appends
 * and extracts the object's members.
 *
 * The object can put-to and get-from it's typeid() to add type checking.  This operation ensures
 * that the data types being read was the data type written before the data is extracted.  This type
 * checking can be disabled since it may be desired to put-to an object of one type, but get-from
 * into an object of an extractable but different type.
 *
 * The TYPE_CHECK bit masks can be provided at put-to Marshal construction to activate the type
 * checking.  The Marshaller send the type check code as the first message to allow the get-from to
 * initialize properly.
 *
 * The put-to operator and get-from operators for plain old data, std::string, std::vector and
 * std::list have been implemented.  Additional ones could be added here, or left to the developer
 * using the marshaller.
 *
 * The stream and type_check members were left as public due to the extensive use.  If this proves
 * bothersome, getter/setter methods could be introduced.
 *
 */
class Marshal
{
public:
  /**
   * @brief Enumeration to activate type checking for std classes and plain old data.
   *
   */
  enum {
    TYPE_CHECK_NONE     = 0x00000000,
    TYPE_CHECK_POD      = 0x00000001,
    TYPE_CHECK_LIST     = 0x00000002,
    TYPE_CHECK_VECTOR   = 0x00000004,
    TYPE_CHECK_SET      = 0x00000008,
    TYPE_CHECK_MAP      = 0x00000010,
    TYPE_CHECK_UNORDERED_SET = 0x00000020,
    TYPE_CHECK_ALL      = 0xFFFFFFFF
  };

  /**
   * Creates a new <code>Marshal</code> instance for put-to operations.
   *
   */
  Marshal(unsigned type_check = TYPE_CHECK_NONE);

  /**
   * Creates a new <code>Marshal</code> instance for get-from operations.
   *
   * @param s			a <code>std::string</code> constant variable of packed bytes to
   *                            extract using the get-from operators.
   */
  explicit Marshal(const std::string &s);

  /**
   * @brief Member function <code>str</code> returns the string of packed bytes created by put-to
   * operations to the stream.
   *
   * @return			a <code>std::string</code> created from the packed byte stream.
   */
  std::string str() const;

  //-----------------------------------------------------------------------------
  // Function      : str
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Tue Feb  3 08:32:28 2015
  //-----------------------------------------------------------------------------
  ///
  /// 
  ///
  /// @invariant
  ///
  /// @param s 
  ///
  ///
  void str(const std::string &s);

  /**
   * @brief Member function <code>size</code> returns the byte count of the string of packed bytes
   * creates by put-to operations to the stream.
   *
   * @return			a <code>size_t</code> in bytes of the packed byte stream.
   */
  size_t size() const;

  /**
   * @brief Member function <code>write</code> writer bytes to the packed byte stream.
   *
   * @param byte_count		a <code>size_t</code> value of the number of packed bytes to write.
   *
   * @param address		a <code>char</code> constant pointer to get the bytes from.
   *
   */
  void write(const char *address, size_t byte_count);

  /**
   * @brief Member function <code>read</code> reads bytes from the packed byte stream.
   *
   * @param byte_count		a <code>size_t</code> value of the number of packed bytes to read.
   *
   * @param address		a <code>char</code> constant pointer to put the bytes to.
   *
   */
  void read(char *address, size_t byte_count);

  /**
   * @brief Member function <code>operator void *</code> returns the state of the packed byte stream.
   *
   * @return			a <code>void</code> const pointer which is non-zero if status is
   *                            good.
   */
  operator void * () const;

private:
  Marshal(const Marshal &marshal);                      ///< Not copyable
  Marshal &operator=(const Marshal &);                  ///< Not assignable

public:
  std::stringstream     stream;                         ///< Packed byte stream to put-to or get-from
  unsigned              m_typeCheck;                    ///< Type checking to activate
};


/**
 * @brief Function <code>operator<< </code> writes the object to the packed byte stream.  This is
 * the template class and has no implementation.  You must specialize this class to write an
 * object.
 *
 * @param mout  		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>T</code> const reference to the object to write.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template <typename T>
Marshal &operator<<(Marshal &mout, const T &t);

/**
 * @brief Function <code>operator>> </code> reads the object from the packed byte stream.  This is
 * the template class and has no implementation.  You must specialize this class to read an object.
 *
 * @param min     		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>T</code> const reference to the object to read.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template <typename T>
Marshal &operator>>(Marshal &min, T &t);

/**
 * @brief Function <code>operator<< </code> write the crc32 encoding of the name from the type
 * information to the packed byte stream.  When the bytes are read, the crc32 encoding of the type
 * being read is varified.
 *
 * @param mout  		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>std::type_info</code> const reference to the type
 *                              information to write for verification when read on extraction.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template<>
Marshal &operator<<(Marshal &mout, const std::type_info &t);

/**
 * @brief Function <code>operator<< </code> reads the crc32 encoding of the name from the type
 * information from the packed byte stream.  The read crc32 is compared to the crc32 encoding of the
 * name from the type information passed.  If the two are different and exception is thrown.
 *
 * @param min     		a <code>Marshal</code> reference to the marshaller.
 *
 * @param t       		a <code>std::type_info</code> const reference to the type
 *                              information to compare with the what was read from the packed byte
 *                              stream.
 *
 * @return			a <code>Marshal</code> reference to the marhsaller.
 */
template<>
Marshal &operator>>(Marshal &min, const std::type_info &t);

template<>
Marshal &operator<<(Marshal &mout, const bool &t);
template<>
Marshal &operator<<(Marshal &mout, const signed char &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned char &t);
template<>
Marshal &operator<<(Marshal &mout, const char &t);
template<>
Marshal &operator<<(Marshal &mout, const short &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned short &t);
template<>
Marshal &operator<<(Marshal &mout, const int &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned int &t);
template<>
Marshal &operator<<(Marshal &mout, const long &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned long &t);
template<>
Marshal &operator<<(Marshal &mout, const long long &t);
template<>
Marshal &operator<<(Marshal &mout, const unsigned long long &t);
template<>
Marshal &operator<<(Marshal &mout, const float &t);
template<>
Marshal &operator<<(Marshal &mout, const double &t);
template<>
Marshal &operator<<(Marshal &mout, const std::string &s);

template<>
Marshal &operator>>(Marshal &min, bool &t);
template<>
Marshal &operator>>(Marshal &min, signed char &t);
template<>
Marshal &operator>>(Marshal &min, unsigned char &t);
template<>
Marshal &operator>>(Marshal &min, char &t);
template<>
Marshal &operator>>(Marshal &min, short &t);
template<>
Marshal &operator>>(Marshal &min, unsigned short &t);
template<>
Marshal &operator>>(Marshal &min, int &t);
template<>
Marshal &operator>>(Marshal &min, unsigned int &t);
template<>
Marshal &operator>>(Marshal &min, long &t);
template<>
Marshal &operator>>(Marshal &min, unsigned long &t);
template<>
Marshal &operator>>(Marshal &min, long long &t);
template<>
Marshal &operator>>(Marshal &min, unsigned long long &t);
template<>
Marshal &operator>>(Marshal &min, float &t);
template<>
Marshal &operator>>(Marshal &min, double &t);
template<>
Marshal &operator>>(Marshal &min, std::string &s);

template<class T, class U>
Marshal &operator<<(Marshal &mout, const std::pair<T, U> &p) {
  mout << p.first << p.second;
  return mout;
}

template<class T, class U>
  Marshal &operator>>(Marshal &min, std::pair<T, U> &p) {
  min >> p.first >> p.second;
  return min;
}

template <class T, class A>
Marshal &operator<<(Marshal &mout, const std::vector<T, A> &v)
{
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_VECTOR)
    mout << typeid(v);

  size_t size = v.size();
  mout << size;
  for (typename std::vector<T, A>::const_iterator it = v.begin(); it != v.end(); ++it)
    mout << (*it);

  return mout;
}

template <class T, class A>
Marshal &operator>>(Marshal &min, std::vector<T, A> &v)
{
  if (min.m_typeCheck & Marshal::TYPE_CHECK_VECTOR)
    min >> typeid(v);

  size_t size = 0;
  min >> size;
  v.reserve(size);
  for (size_t i = 0; i < size; ++i)
  {
    T t;
    min >> t;
    v.push_back(t);
  }

  return min;
}

template <class T, class A>
Marshal &operator<<(Marshal &mout, const std::list<T, A> &l)
{
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_LIST)
    mout << typeid(l);

  size_t size = l.size();
  mout << size;
  for (typename std::list<T>::const_iterator it = l.begin(); it != l.end(); ++it)
    mout << (*it);

  return mout;
}

template <class T, class A>
Marshal &operator>>(Marshal &min, std::list<T, A> &l)
{
  if (min.m_typeCheck & Marshal::TYPE_CHECK_LIST)
    min >> typeid(l);

  size_t size;
  min >> size;
  for (size_t i = 0; i < size; ++i)
  {
    T t;
    min >> t;
    l.push_back(t);
  }

  return min;
}

template <class T, class C, class A>
Marshal &operator<<(Marshal &mout, const std::set<T, C, A> &s)
{
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_SET)
    mout << typeid(s);

  size_t size = s.size();
  mout << size;
  for (typename std::set<T, C, A>::const_iterator it = s.begin(); it != s.end(); ++it)
    mout << (*it);

  return mout;
}

template <class T, class C, class A>
Marshal &operator>>(Marshal &min, std::set<T, C, A> &s)
{
  if (min.m_typeCheck & Marshal::TYPE_CHECK_SET)
    min >> typeid(s);

  size_t size;
  min >> size;
  for (size_t i = 0; i < size; ++i)
  {
    T t;
    min >> t;
    s.insert(t);
  }

  return min;
}

template <class K, class T, class C, class A>
Marshal &operator<<(Marshal &mout, const std::map<K, T, C, A> &m)
{
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_MAP)
    mout << typeid(m);

  size_t size = m.size();
  mout << size;
  for (typename std::map<K, T, C, A>::const_iterator it = m.begin(); it != m.end(); ++it)
    mout << (*it).first << (*it).second;

  return mout;
}

template <class K, class T, class C, class A>
Marshal &operator>>(Marshal &min, std::map<K, T, C, A> &m)
{
  if (min.m_typeCheck & Marshal::TYPE_CHECK_MAP)
    min >> typeid(m);

  size_t size;
  min >> size;
  for (size_t i = 0; i < size; ++i)
  {
    K k;
    T t;
    min >> k >> t;
    m.insert(typename std::map<K, T, C, A>::value_type(k, t));
  }

  return min;
}

template <class T, class C, class A>
Marshal &operator<<(Marshal &mout, const unordered_set<T, C, A> &s)
{
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_UNORDERED_SET)
    mout << typeid(s);

  size_t size = s.size();
  mout << size;
  for (typename unordered_set<T, C, A>::const_iterator it = s.begin(); it != s.end(); ++it)
    mout << (*it);

  return mout;
}

template <class T, class C, class A>
Marshal &operator>>(Marshal &min, unordered_set<T, C, A> &s)
{
  if (min.m_typeCheck & Marshal::TYPE_CHECK_UNORDERED_SET)
    min >> typeid(s);

  size_t size;
  min >> size;
  for (size_t i = 0; i < size; ++i)
  {
    T t;
    min >> t;
    s.insert(t);
  }

  return min;
}

template <class T>
Marshal &write(Marshal &mout, const T &t)
{
  mout.write((const char *) &t, sizeof(T));

  return mout;
}

template <typename T>
Marshal &read(Marshal &min, T &t)
{
  t = T();

  min.read((char *) &t, sizeof(T));
  return min;
}

} // namespace Util
} // namespace Xyce

#endif // N_UTL_MARSHAL_HPP

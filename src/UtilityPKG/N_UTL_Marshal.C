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
// Creation Date  : 1/7/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <stdexcept>
#include <map>
#include <cstring>
#include <iostream>

#include <Xyce_config.h>

#include <N_UTL_Marshal.h>

namespace Xyce {
namespace Util {

namespace {

typedef std::map<uint32_t, std::string> TypeNameMap;

bool            s_crcInitialized = false;               ///< CRC table has been initialized
uint32_t        s_crcTable[256];                        ///< CRC lookup table
TypeNameMap     s_typeNameMap;                          ///< CRC'd type name map

uint32_t
crc32_reflect(
  uint32_t      seed,
  const char    c)
{
  uint32_t value = 0;

  // Swap bit 0 for bit 7, bit 1 For bit 6, etc....
  for(int i = 1; i < (c + 1); i++) {
    if (seed & 1)
      value |= (1 << (c - i));
    seed >>= 1;
  }

  return value;
}


void
crc32_initialize()
{
  s_crcInitialized = true;

  //0x04C11DB7 is the official polynomial used by PKZip, WinZip and Ethernet.
  uint32_t polynomial = 0x04C11DB7;

  std::fill(s_crcTable, s_crcTable + 256, 0);

  // 256 values representing ASCII character codes.
  for (int i = 0; i <= 0xFF; i++) {
    s_crcTable[i] = crc32_reflect(i, 8) << 24;

    for (int j = 0; j < 8; j++)
      s_crcTable[i] = (s_crcTable[i] << 1) ^ ((s_crcTable[i] & (1u << 31)) ? polynomial : 0);

    s_crcTable[i] = crc32_reflect(s_crcTable[i], 32);
  }
}


void
crc32_part(
  uint32_t &            crc,
  const unsigned char * s,
  unsigned              l)
{
  while (l--)
    crc = (crc >> 8)^s_crcTable[(crc & 0xFF)^*s++];
}


uint32_t
crc32(
  const unsigned char * s,
  unsigned              l)
{
  uint32_t crc = 0xffffffff;
  crc32_part(crc, s, l);
  return crc ^ 0xffffffff;
}


void
crc32_write(
  std::stringstream &           os,
  const std::type_info &        typeinfo)
{
  if (!s_crcInitialized)
    crc32_initialize();

  const char *name = typeinfo.name();
  unsigned length = std::strlen(name);

  uint32_t      crc = crc32((const unsigned char *) name, length);

  TypeNameMap::iterator it = s_typeNameMap.find(crc);
  if (it == s_typeNameMap.end())
    s_typeNameMap.insert(std::make_pair(crc, std::string(name)));

  os.write((const char *) &crc, sizeof(uint32_t));
}


void
crc32_check(
  std::stringstream &           is,
  const std::type_info &        typeinfo)
{
  if (!s_crcInitialized)
    crc32_initialize();

  const char *name = typeinfo.name();
  unsigned length = std::strlen(name);

  uint32_t      crc_check = crc32((const unsigned char *) name, length);
  uint32_t      crc = 0;

  {
    TypeNameMap::iterator it = s_typeNameMap.find(crc_check);
    if (it == s_typeNameMap.end())
      s_typeNameMap.insert(std::make_pair(crc_check, std::string(name)));
  }

  is.read((char *) &crc, sizeof(uint32_t));

  if (crc_check != crc) {
    std::ostringstream ss;
    ss << "Marshaller encountered type ";
    TypeNameMap::const_iterator it = s_typeNameMap.find(crc);
    if (it == s_typeNameMap.end())
      ss << "code " << std::hex << crc;
    else
      ss << (*it).second;

    ss << " when expecting type ";
    it = s_typeNameMap.find(crc_check);
    if (it == s_typeNameMap.end())
      ss << "code " << std::hex << crc_check;
    else
      ss << (*it).second;

    throw std::runtime_error(ss.str());
  }
}

} // namespace <empty>

Marshal::Marshal(
  unsigned              type_check)
  : stream(std::ios_base::out),
    m_typeCheck(TYPE_CHECK_NONE)
{
  m_typeCheck = type_check;
  write(reinterpret_cast<char *>(&m_typeCheck), sizeof(m_typeCheck));
}


Marshal::Marshal(
  const std::string &   s)
  : stream(s, std::ios_base::in),
    m_typeCheck(TYPE_CHECK_NONE)
{
  read(reinterpret_cast<char *>(&m_typeCheck), sizeof(m_typeCheck));
}

void
Marshal::str(
  const std::string &   s)
{
  stream.str(s);
  stream.seekg(0);
  if (!stream)
    stream.clear();

  read(reinterpret_cast<char *>(&m_typeCheck), sizeof(m_typeCheck));
}


std::string
Marshal::str() const
{
  return stream.str();
}


size_t
Marshal::size() const
{
  return stream.str().size();
}


Marshal::operator void * () const
{
  return static_cast<void *>(const_cast<std::stringstream*>(&stream));
}


void
Marshal::write(
  const char *          address,
  size_t                byte_count)
{
  stream.write(address, byte_count);
}


void
Marshal::read(
  char *                address,
  size_t                byte_count)
{
  stream.read(address, byte_count);
}


template<>
Marshal &operator<<(Marshal &mout, const std::type_info &t) {
  crc32_write(mout.stream, t);
  return mout;
}

template<>
Marshal &operator>>(Marshal &min, const std::type_info &t) {
  crc32_check(min.stream, t);
  return min;
}

template<>
Marshal &operator<<(Marshal &mout, const bool &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const signed char &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned char &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const char &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const short &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned short &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const int &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned int &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const long long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned long long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const float &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const double &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const std::string &s)  {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(s);

  size_t ul = s.size();
  mout << ul;
  mout.stream.write(s.data(), s.size());
  return mout;
}


template<>
Marshal &operator>>(Marshal &min, bool &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, signed char &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned char &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, char &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, short &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned short &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, int &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned int &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, long long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned long long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, float &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, double &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, std::string &s)  {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(s);

  size_t size = 0;
  min >> size;
  std::vector<char> c(size);

  min.stream.read(&c[0], size);
  s.assign(&c[0], size);

  return min;
}

} // namespace Util
} // namespace Xyce

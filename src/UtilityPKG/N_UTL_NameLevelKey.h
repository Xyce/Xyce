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

//-------------------------------------------------------------------------
//
// Purpose        : Map key for Device name and level
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NameLevelKey_h
#define Xyce_N_UTL_NameLevelKey_h

#include <utility>
#include <string>
#include <functional>

#include <N_UTL_NoCase.h>

namespace Xyce {

struct NameLevelKey : public std::pair<std::string, int>
{
  NameLevelKey()
    : std::pair<std::string, int>()
  {}

  NameLevelKey(const std::string &s, int i)
    : std::pair<std::string, int>(s, i)
  {}
};

struct NameLevelKeyLess : public std::binary_function<NameLevelKey, NameLevelKey, bool>
{
  bool operator()( const NameLevelKey &lhs , const NameLevelKey &rhs ) const
  {
    int i = compare_nocase(lhs.first.c_str(), rhs.first.c_str());

    if (i == 0)
      return lhs.second < rhs.second;
    else
      return i < 0;
  }
};

std::ostream &operator<<(std::ostream &os, const NameLevelKey &device_level_key);

} // namespace Xyce


namespace std {

// template<>
// struct less<Xyce::NameLevelKey> : public std::binary_function<Xyce::NameLevelKey, Xyce::NameLevelKey, bool>
// {
//   bool operator()(const Xyce::NameLevelKey &lhs, const Xyce::NameLevelKey &rhs) const
//   {
//     int i = Xyce::compare_nocase(lhs.first.c_str(), rhs.first.c_str());

//     if (i == 0)
//       return lhs.second < rhs.second;
//     else
//       return i < 0;
//   }
// };

template<>
struct equal_to<Xyce::NameLevelKey> : public std::binary_function<Xyce::NameLevelKey, Xyce::NameLevelKey, bool>
{
  bool operator()(const Xyce::NameLevelKey &lhs, const Xyce::NameLevelKey &rhs) const
  {
    Xyce::EqualNoCase x0;
    equal_to<int> x1;

    return x1(lhs.second, rhs.second) && x0(lhs.first, rhs.first);
  }
};

template<>
struct hash<Xyce::NameLevelKey> : public std::unary_function<Xyce::NameLevelKey, size_t>
{
  size_t operator()(const Xyce::NameLevelKey &node_id) const
  {
    Xyce::HashNoCase x0;
    hash<int> x1;

    return x0(node_id.first) ^ x1(node_id.second);
  }
};

} // std namespace

#endif // Xyce_N_UTL_NameLevelKey_h

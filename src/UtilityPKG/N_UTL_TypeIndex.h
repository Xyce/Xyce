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

#ifndef Xyce_N_UTL_TypeIndex_h
#define Xyce_N_UTL_TypeIndex_h

#include <typeinfo>

namespace Xyce {

/**
 * type_index wraps std::type_info for use as a key for maps and for comparison.
 *
 * When constructed, the 'undefined' value is assigned.  Do not try to use as a key if it has not been defined.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Feb  5 11:06:01 2014
 */
struct type_index
{
  /**
   * type_index constructs an undefined type_index.
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Wed Feb  5 11:08:56 2014
   */
  type_index()
    : type_(0)
  {}

  /**
   * type_index constructs a type_index for the type.
   *
   * @param type std::type_info reference
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Wed Feb  5 11:09:28 2014
   */
  type_index(const std::type_info &type)
    : type_(&type)
  {}

  /**
   * type_index copy the type_index
   *
   * @param t const reference to type_index to copy
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Wed Feb  5 11:14:09 2014
   */
  type_index(const type_index &t)
    : type_(t.type_)
  {}

  /**
   * operator= assigns the type_index
   *
   * @param t const reference to type_index to assign
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Wed Feb  5 11:14:40 2014
   */
  type_index &operator=(const type_index &t)
  {
    type_ = t.type_;

    return *this;
  }

  /**
   * type returns the typeid, undefined behavior if the type_index is undefined.
   *
   * @return const reference to the std::type_info
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Wed Feb  5 11:10:29 2014
   */
  const std::type_info &type() const 
  {
    return *type_;
  }

  /**
   * defined returns type if the type_index has been defined.
   *
   * @return type if the type_index has been defined
   *
   * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
   * @date   Wed Feb  5 11:11:58 2014
   */
  bool defined() const 
  {
    return type_ != 0;
  }

private:
  const std::type_info *    type_;    ///< pointer to the std::type_info or 0 if undefined
};

/**
 * operator== compare equal
 *
 * Undefined behavior if either type_index is undefined.
 *
 * @return true if the types contained in the type_indexes are equal.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Feb  5 11:15:34 2014
 */
inline
bool operator==(const type_index &type_left, const type_index &type_right) 
{
  return type_left.type() == type_right.type();
}

/**
 * operator!= compare not equal
 *
 * Undefined behavior if either type_index is undefined.
 *
 * @return true if the types contained in the type_indexes are not equal.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Feb  5 11:15:34 2014
 */
inline
bool operator!=(const type_index &type_left, const type_index &type_right) 
{
  return !(type_left.type() == type_right.type());
}

/**
 * operator< compare less
 *
 * Undefined behavior if either type_index is undefined.
 *
 * @return true if the type of the left type_index is less than the type of the right type_index.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Wed Feb  5 11:15:34 2014
 */
inline
bool operator<(const type_index &type_left, const type_index &type_right) 
{
  return type_left.type().before(type_right.type());
}

} // namespace Xyce

#endif // Xyce_N_UTL_TypeIndex_h

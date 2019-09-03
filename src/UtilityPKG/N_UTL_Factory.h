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

//-----------------------------------------------------------------------------
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Robert Hoekstra
//
// Creation Date : 03/01/2002
//
//
//
//
//-----------------------------------------------------------------------------

///
/// @file   N_UTL_Factory.h
/// @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
/// @date   Thu Jan 29 12:51:28 2015
///
/// @brief  Factory for creating analysis objects
///
///
#ifndef Xyce_N_UTL_Factory_h
#define Xyce_N_UTL_Factory_h

#include <typeinfo>

#include <N_LAS_fwd.h>

namespace Xyce {
namespace Util {

template<class B, class T>
class Factory;

//-----------------------------------------------------------------------------
// Template      : Factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:40:53 2015
//-----------------------------------------------------------------------------
///
/// The factory template defines an interface for type testing and object creation.
///
template<>
class Factory<void, void>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : Factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:41:54 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructor
  ///
  Factory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~Factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:43:11 2015
  //-----------------------------------------------------------------------------
  ///
  /// Destructor
  ///
  virtual ~Factory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : type
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:43:43 2015
  //-----------------------------------------------------------------------------
  ///
  /// Defines the interface to get the type info of the object created by the factory.
  ///
  /// @return type info of the object would be created.
  ///
  ///
  virtual const std::type_info &type() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : isType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:45:04 2015
  //-----------------------------------------------------------------------------
  ///
  /// Returns true if the type of the factory matches class U
  ///
  /// @param U object class type to test
  ///
  /// @return true if the object type of the factory matches class U
  ///
  ///
  template <class U>
  bool isType() const {
    return type() == typeid(U);
  }

private:
  Factory(const Factory &);                     ///< not copyable
  Factory &operator=(const Factory &);          ///< not assignable
};

//-----------------------------------------------------------------------------
// Template      : Factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:40:53 2015
//-----------------------------------------------------------------------------
///
/// The object factory template defines an interface for object type testing and
/// object creation.
///
template<class B>
class Factory<B, void> : public Factory<void, void>
{
public:
  typedef B Base;

  //-----------------------------------------------------------------------------
  // Function      : Factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:41:54 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructor
  ///
  Factory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~Factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:43:11 2015
  //-----------------------------------------------------------------------------
  ///
  /// Destructor
  ///
  virtual ~Factory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : type
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:43:43 2015
  //-----------------------------------------------------------------------------
  ///
  /// Defines the interface to get the type info of the object created by the factory.
  ///
  /// @return type info of the object that would be created.
  ///
  ///
  virtual const std::type_info &type() const = 0;

  //-----------------------------------------------------------------------------
  // Function      : isType
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:45:04 2015
  //-----------------------------------------------------------------------------
  ///
  /// Returns true if the object type of the factory matches class U
  ///
  /// @param U Object class type to test
  ///
  /// @return true if the object type of the factory matches class U
  ///
  ///
  template <class U>
  bool isType() const {
    return type() == typeid(U);
  }

private:
  Factory(const Factory &);                     ///< not copyable
  Factory &operator=(const Factory &);          ///< not assignable

public:
  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:47:13 2015
  //-----------------------------------------------------------------------------
  ///
  /// Creates the object object.
  ///
  /// @return the new object object
  ///
  ///
  virtual Base *create() const = 0;
};

//-----------------------------------------------------------------------------
// Template      : Factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:40:53 2015
//-----------------------------------------------------------------------------
///
/// The object factory template defines an interface for object type testing and
/// object creation.  This template implements the type checking.
///
/// @param T object type to be created
///
template<class B, class T>
class Factory : public Factory<B, void>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : Factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:49:44 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructor
  ///
  Factory()
    : Factory<B, void>()
  {}

  //-----------------------------------------------------------------------------
  // Function      : ~Factory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:50:03 2015
  //-----------------------------------------------------------------------------
  ///
  /// Destructor
  ///
  virtual ~Factory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : type
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:50:21 2015
  //-----------------------------------------------------------------------------
  ///
  /// Returns the type info the object type.
  ///
  /// @return the type info the object type.
  ///
  virtual const std::type_info &type() const {
    return typeid(T);
  }
};

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_Factory_h

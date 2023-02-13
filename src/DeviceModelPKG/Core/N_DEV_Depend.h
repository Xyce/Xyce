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

//-----------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 03/01/22
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Depend_h
#define Xyce_N_DEV_Depend_h

#include <string>
#include <vector>

namespace Xyce {
namespace Device {
//-----------------------------------------------------------------------------
// Class         : Depend
// Purpose       : Used to record information about dependent parameters
// Special Notes :
// Creator       : Dave Shirley
// Creation Date : 
//-----------------------------------------------------------------------------
///
///  The Depend struct is used to keep track of dependent parameters
///
struct Depend
{
  std::string                 name;      ///< parameter name
  Util::Expression *          expr;      ///< expression used comput value
  union resUnion
  {
    int *                  iresult;
    double *                result;
    std::vector<double> *   resVec;
  } resultU;                            ///< Holds a pointer to where the
                                        ///< parameter is stored.
  int                         vectorIndex; ///< Used if parameter is in a vector

  int                         n_vars, lo_var, n_global; 
  bool                     storeOriginal;    ///< true if original value stored
  int                      serialNumber;     ///< used if original value stored

  // Constructor
  Depend()
    : vectorIndex(-1), n_vars(0), lo_var(0), n_global(0)
  {};

};

// ERK.  this could be replaced by a lambda
struct MatchDependName
{
  MatchDependName(const std::string& name) : matchName_(name) {}
  bool operator()(const Depend & dep) const
  {
    return dep.name == matchName_;
  }
  private:
    const std::string& matchName_;
};

struct Depend_lesser
{
  bool operator ()(Xyce::Device::Depend const& a, Xyce::Device::Depend const& b) const 
  {
    return (a.name < b.name);
  }
};

struct Depend_greater
{
  bool operator ()(Xyce::Device::Depend const& a, Xyce::Device::Depend const& b) const 
  {
    return (a.name > b.name);
  }
};

struct Depend_equal
{
  bool operator ()(Xyce::Device::Depend const& a, Xyce::Device::Depend const& b) const 
  {
    return (a.name == b.name);
  }
};

}
}

#endif


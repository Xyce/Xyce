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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef  N_DEV_PARAM_H
#define  N_DEV_PARAM_H

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include <N_DEV_fwd.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_Param
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
class Param : public Util::Param
{
  friend class Pack<Param>;

public:
  Param()
    : Util::Param(),
      isGiven_(false),
      isDefault_(false)
  {}

  template <class T>
  Param(const std::string &tag, const T &value, bool is_given = false)
    : Util::Param(tag, value),
      isGiven_(is_given),
      isDefault_(false)
  {}

  Param(const Param &rhsParam)
    : Util::Param(rhsParam),
      isGiven_(rhsParam.isGiven_),
      isDefault_(rhsParam.isDefault_)
  {}

  Param &operator=(const Param &rhsParam) 
  {
    Util::Param::operator=(rhsParam);
    isGiven_ = rhsParam.isGiven_;
    isDefault_ = rhsParam.isDefault_;

    return *this;
  }

  virtual ~Param()
  {}

  void setGiven(bool is_given) 
  {
    isGiven_ = is_given;
  }

  void setDefault(bool is_default) 
  {
    isDefault_ = is_default;
  }

  bool given() const 
  {
    return isGiven_;
  }

  bool default_val() const 
  {
    return isDefault_;
  }

private:
  bool isGiven_;
  bool isDefault_;
};

inline void setParamValue(Param &param, const Param &from_param) 
{
  param.setVal(static_cast<const Util::Param &>(from_param));
}

inline void setParam(Param &param, const std::string &tag, const Param &from_param) 
{
  param.set(tag, static_cast<const Util::Param &>(from_param));
}

} // namespace Device
} // namespace Xyce

#endif

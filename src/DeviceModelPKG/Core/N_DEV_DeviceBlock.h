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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceBlock_h
#define Xyce_N_DEV_DeviceBlock_h

#include <string>
#include <vector>
#include <iosfwd>

#include <N_DEV_Param.h>
#include <N_DEV_InstanceName.h>
#include <N_UTL_NetlistLocation.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : ModelBlock
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
///
/// ModelBlock represents a .MODEL line from the netlist.
///
class ModelBlock
{ 
  friend class Pack<ModelBlock>;
  friend std::ostream& operator<<(std::ostream& os, const ModelBlock & mb);

public:
  ModelBlock(const std::string &name = "", const std::string &type = "", int level = 1);

  ModelBlock(const ModelBlock & right);
  ModelBlock &operator=(const ModelBlock & right);

  ~ModelBlock();

  const ModelName &getName() const
  {
    return name_;
  }

  void setName(const ModelName &name)
  {
    name_ = name;
  }

  const std::string &getType() const
  {
    return type_;
  }

  void setType(const std::string &type)
  {
    type_ = type;
  }

  int getLevel() const
  {
    return level_;
  }

  void setLevel(int level)
  {
    level_ = level;
  }

  const NetlistLocation &getNetlistLocation() const
  {
    return netlistLocation_;
  }

  void setNetlistLocation(const NetlistLocation &netlist_location)
  {
    netlistLocation_ = netlist_location;
  }

  bool operator==(const ModelBlock &right) const
  {
    return equal_nocase(name_, right.name_);
  }

  bool operator!=(const ModelBlock &right) const
  {
    return !equal_nocase(name_, right.name_);
  }

  void clear();

private:
  ModelName             name_;                  ///< Model name
  std::string           type_;                  ///< Model type
  int                   level_;                 ///< Device level
  NetlistLocation       netlistLocation_;       ///< Path and line number of .MODEL command

public:
  std::vector<Param>    params;                 ///< Parameters from the line
};

//-----------------------------------------------------------------------------
// Class         : InstanceBlock
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
///
/// InstanceBlock represent a device instance line from the netlist.
///
class InstanceBlock
{
  friend class Pack<InstanceBlock>;
  friend std::ostream& operator<<(std::ostream& os, const InstanceBlock & ib);

public:
  InstanceBlock(const std::string &name = std::string());

  InstanceBlock(const InstanceBlock &right);
  InstanceBlock &operator=(const InstanceBlock &right);

  ~InstanceBlock();

  const InstanceName &getInstanceName() const 
  {
    return name_;
  }

  void setInstanceName(const InstanceName &name) 
  {
    name_ = name;
  }

  const ModelName &getModelName() const 
  {
    return modelName_;
  }

  void setModelName(const ModelName &modelName) 
  {
    modelName_ = modelName;
  }

  const NetlistLocation &getNetlistLocation() const
  {
    return netlistLocation_;
  }

  void setNetlistLocation(const NetlistLocation &netlist_location)
  {
    netlistLocation_ = netlist_location;
  }

  bool operator==(const InstanceBlock &right) const
  {
    return equal_nocase(name_.getEncodedName(), right.name_.getEncodedName());
  }

  bool operator!=(const InstanceBlock &right) const
  {
    return !equal_nocase(name_.getEncodedName(), right.name_.getEncodedName());
  }

  void clear ();

private:
  InstanceName          name_;                  ///< Device instance name
  ModelName             modelName_;             ///< Model name if provided
  NetlistLocation       netlistLocation_;       ///< Path and line number of .MODEL command

public:
  std::vector<Param> params;

  int iNumNodes;
  int numIntVars;
  int numExtVars;
  int numStateVars;

  bool modelFlag;
  bool bsourceFlag;
  bool offFlag;
  bool off;
};

} // namespace Device
} // namespace Xyce

#endif

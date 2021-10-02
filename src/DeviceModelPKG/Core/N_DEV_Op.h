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

//-----------------------------------------------------------------------------
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Op_h
#define Xyce_N_DEV_Op_h

#include <iterator>

#include <N_DEV_fwd.h>
#include <N_UTL_Op.h>

namespace Xyce {
namespace Device {

class DeviceMgrGlobalParameterOp : public Util::Op::Op<DeviceMgrGlobalParameterOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  DeviceMgrGlobalParameterOp(const std::string &name, const DeviceMgr &device_manager, const double &global_parameter_value)
    : Base(name),
      globalParameterValue_(global_parameter_value)
  {}

  virtual ~DeviceMgrGlobalParameterOp()
  {}

  static complex get(const DeviceMgrGlobalParameterOp &op, const Util::Op::OpData &op_data);

  const double &                globalParameterValue_;
};

class DeviceEntityParameterOp : public Util::Op::Op<DeviceEntityParameterOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
  DeviceEntityParameterOp(const std::string &name, const DeviceEntity &device_entity, const std::string &device_parameter_name)
    : Base(name),
      deviceEntity_(device_entity),
      deviceParameterName_(device_parameter_name)
  {}

  virtual ~DeviceEntityParameterOp()
  {}

  static complex get(const DeviceEntityParameterOp &op, const Util::Op::OpData &op_data);

  const DeviceEntity &          deviceEntity_;
  const std::string             deviceParameterName_;
};

class ArtificialParameterOp : public Util::Op::Op<ArtificialParameterOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  ArtificialParameterOp(const std::string &name, const DeviceMgr &device_manager, const ArtificialParameters::ArtificialParameter &artificial_parameter, const std::string &artificial_parameter_name)
    : Base(name),
      deviceManager_(device_manager),
      artificialParameterName_(artificial_parameter_name),
      artificialParameter_(artificial_parameter)
  {}

  virtual ~ArtificialParameterOp()
  {}

  static complex get(const ArtificialParameterOp &op, const Util::Op::OpData &op_data);

  const DeviceMgr &                                     deviceManager_;
  const std::string                                     artificialParameterName_;
  const ArtificialParameters::ArtificialParameter &     artificialParameter_;
};

// class DeviceMgrParameterOp : public Util::Op::Op<DeviceMgrParameterOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
// {
// public:
//   DeviceMgrParameterOp(const std::string &name, const DeviceMgr &device_manager, const std::string &device_parameter_name)
//     : Base(name),
//       deviceManager_(device_manager),
//       deviceParameterName_(device_parameter_name)
//   {}


//   virtual ~DeviceMgrParameterOp()
//   {}

//   static complex get(const DeviceMgrParameterOp &op, const Util::Op::OpData &op_data);

//   const std::string     deviceParameterName_;
//   const DeviceMgr &     deviceManager_;
// };

class DeviceOptionsOp : public Util::Op::Op<DeviceOptionsOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  DeviceOptionsOp(const std::string &name, const DeviceOptions &device_options, const std::string &option_name)
    : Base(name),
      deviceOptions_(device_options),
      optionName_(option_name)
  {}


  virtual ~DeviceOptionsOp()
  {}

  static complex get(const DeviceOptionsOp &op, const Util::Op::OpData &op_data);

  const DeviceOptions &         deviceOptions_;
  const std::string             optionName_;
};

class MutualInductorInstancesOp : public Util::Op::Op<MutualInductorInstancesOp, Util::Op::ReduceSum, Util::Op::EvalNoop>
{
public:
    MutualInductorInstancesOp(const std::string &name,
       const std::string & inductor_name,
       const DeviceInstance &device_instance, int index)
    : Base(name),
      inductorName_(inductor_name),
      deviceInstance_(device_instance),
      inductorIndex_(index)
  {}

  virtual ~MutualInductorInstancesOp()
  {}

  static complex get(const MutualInductorInstancesOp &op, const Util::Op::OpData &op_data);

  const std::string               inductorName_;
  const DeviceInstance &          deviceInstance_;
  const int                       inductorIndex_;
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Op_h

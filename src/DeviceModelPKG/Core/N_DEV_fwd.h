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
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_fwd_h
#define Xyce_N_DEV_fwd_h

#if defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
using std::unordered_map;
#elif defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#else
#error neither unordered_map or tr1/unordered_map found
#endif

#if defined(HAVE_UNORDERED_SET)
#include <unordered_set>
using std::unordered_set;
#elif defined(HAVE_TR1_UNORDERED_SET)
#include <tr1/unordered_set>
using std::tr1::unordered_set;
#else
#error neither unordered_set or tr1/unordered_set found
#endif

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <N_UTL_NameLevelKey.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_TypeIndex.h>

namespace Xyce {
namespace Device {

class CompositeParam;
class Configuration;
struct Depend;
class Descriptor;
class Device;
class DeviceBuilder;
class DeviceEntity;
class DeviceInstance;
class DeviceMgr;
typedef DeviceMgr DeviceInterface; // DEPRECATED, use DeviceMgr

class DeviceModel;
class DeviceOptions;
class DeviceSensitivities;
class DeviceState;
class DeviceSupport;
class InstanceName;
struct ExternalSimulationData;
class ExternCodeInterface;
class ExternData;
struct FactoryBlock;
class InstanceBlock;
class MatrixLoadData;
class ModelBlock;
class NumericalJacobian;
class Param;
class ParameterBase;
class Region;
class RegionData;
class RxnRegion;
class RxnRegion2;
class RxnRegionData;

class SolverState;
struct Globals;
class SourceInstance;
class XyceInterface;

class ACData;
class ConstData;
class ExpData;
class PWLinData;
class PulseData;
class SFFMData;
class SinData;
class SourceData;

class DevicePDEInstance;
class DevicePDEModel;

class PDE_Electrode;
class PDE_1DElectrode;
class PDE_2DElectrode;

class XygraCoilData;

class SpecieSource;

class ScalingVars;

class DeviceMgrGlobalParameterOp;

typedef unordered_map<std::string, InstanceBlock, HashNoCase, EqualNoCase> DeviceNameInstanceBlockMap;

namespace GeneralExternal {
class Instance;
class Model;
}

 namespace Xygra {
class Instance;
class Model;
}

namespace ExternDevice {
class Instance;
class Model;
}

namespace Vsrc {
class Instance;
class Model;
}

namespace ArtificialParameters {
struct ArtificialParameter;
}

typedef type_index EntityTypeId;
typedef type_index ModelTypeId;
typedef type_index InstanceTypeId;

typedef std::string ModelName;

typedef unordered_map<std::string, Descriptor *, HashNoCase, EqualNoCase> ParameterMap;
typedef unordered_map<std::string, CompositeParam *, HashNoCase, EqualNoCase> CompositeMap;
typedef unordered_map<std::string, double, HashNoCase, EqualNoCase> GlobalParameterMap;
typedef unordered_map<std::string, ArtificialParameters::ArtificialParameter *, HashNoCase, EqualNoCase> ArtificialParameterMap;
typedef unordered_set<std::string, HashNoCase, EqualNoCase> PassthroughParameterSet;


typedef std::map<std::string, int, LessNoCase> DeviceCountMap;

typedef std::vector<CompositeParam *> CompositeVector;

typedef std::map<EntityTypeId, Device *> EntityTypeIdDeviceMap;

struct ModelGroupLeadTraits;

template <class M, class I, class G = ModelGroupLeadTraits>
class DeviceTraits;

template <class T>
class DeviceMaster;

template <class T>
class Config;

typedef std::vector<Device *> DeviceVector;
typedef std::vector<DeviceEntity *> EntityVector;
typedef std::vector<DeviceInstance *> InstanceVector;
typedef std::vector<DeviceModel *> ModelVector;

typedef std::vector<std::vector<int> > JacobianStamp;
typedef std::vector<int> IdVector;
typedef IdVector LocalIdVector;

typedef std::pair<int,int> IntPair;
typedef std::vector<IntPair> PairVector;
typedef std::map< IntPair, IntPair > PairMap;
typedef std::set< IntPair > PairSet;

// Loading can distinguish the type of devices to be loaded, for efficiency.
enum loadType { ALL, LINEAR, NONLINEAR, PDE, LINEAR_FREQ, NONLINEAR_FREQ };

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_fwd_h

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
// Purpose        : Anti-Windup Limiter: See Section C.2 of "Power System
//                  Modeling and Scripting" by F. Milano
//
// Special Notes  : Experimental new device for an LDRD.
//
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 11/12/15
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_AntiWindupLimiter_h
#define Xyce_N_DEV_AntiWindupLimiter_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace AntiWindupLimiter {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "AntiWindupLimiter";}
  static const char *deviceTypeName() {return "AntiWindupLimiter";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 0;} 
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
//
//	This  is  the  instance class  for an Anti-Windup Limiter  
//
// Special Notes :
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date  : 11/12/15
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
    
public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     Riter,
     const FactoryBlock &        factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  bool processParams ();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();
  bool loadDAEBVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

public:
  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> >  jacStamp;

  Model &       model_;         //< Owning model

  // user-specified parameters:
  double T_, lowerLimit_, upperLimit_;
  
  // LIDs
  int li_vin ,li_vout, li_BranCurr;

  // used in jac stamps
  int VI_I_Offset, VO_I_Offset; 
  int I_VI_Offset, I_VO_Offset;

  // parameters used to track internal variables
  double VIN_, VOUT_, BranCurr_, QV_;
  bool holdAtUpper_, holdAtLower_;
  bool upperLimitGiven_, lowerLimitGiven_;

};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
//
//
// Special Notes :
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date  : 11/12/15
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
    
public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &          MB,
     const FactoryBlock &        factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;
    
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  virtual bool processParams() 
  {
    return true;
  }

  virtual bool processInstanceParams() 
  {
    return true;
  }

public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

  void setupBaseInstanceContainer()
  {
    std::vector<Instance*>::iterator iter = instanceContainer.begin();
    std::vector<Instance*>::iterator end   = instanceContainer.end();
    for ( ; iter!=end; ++iter)
    {
      Xyce::Device::DeviceModel::baseInstanceContainer.push_back( static_cast<Xyce::Device::DeviceInstance *>(*iter) );
    }
  }

private:
  std::vector<Instance*> instanceContainer;

private:
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace AntiWindupLimiter
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_AntiWindupLimiter_h

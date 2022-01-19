//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : PowerGrid classes: provides a device that calculates the
//                  steady state power flow in a transmission grid bus shunt
//
// Special Notes  : Experimental new device for an LDRD.
//
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 12/4/14
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_PowerGridBusShunt_h
#define Xyce_N_DEV_PowerGridBusShunt_h

#include <complex>

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace PowerGridBusShunt {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "PowerGridBusShunt";}
  static const char *deviceTypeName() {return "PowerGridBusShunt level 1";}
  static int numNodes() {return 4;}
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
//	This  is  the  instance class  for a Power Grid Bus Shunt device.  
//
// Special Notes :
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date  : 12/4/14
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;

  // enum for the analysis type
  enum aType {IV, PQR, PQP};
    
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

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  //bool updateSecondaryState ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );
  bool processParams ();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

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
  std::string analysisTypeStr_;
  double shuntConductance_, shuntSusceptance_;

  // enum-based version of the analysisTypeStr_
  enum aType analysisType_;

  // definitions for power grid bus shunt device
  double g11, g12, g21, g22;
  double b11, b12, b21, b22;
  std::complex<double> y11, y12, y21, y22;

  // these are used only for the I=YV formulation
  double IR1, IR2, II1, II2;

  int li_VR1 ,li_VR2, li_VI1, li_VI2;

  int VR1_VR1_Offset;
  int VR1_VR2_Offset;
  int VR1_VI1_Offset;
  int VR1_VI2_Offset;

  int VR2_VR1_Offset;
  int VR2_VR2_Offset;
  int VR2_VI1_Offset;
  int VR2_VI2_Offset;

  int VI1_VR1_Offset;
  int VI1_VR2_Offset;
  int VI1_VI1_Offset;
  int VI1_VI2_Offset;

  int VI2_VR1_Offset;
  int VI2_VR2_Offset;
  int VI2_VI1_Offset;
  int VI2_VI2_Offset;

  // these are used only for the PQ formulation
  double P1, P2, Q1, Q2;

  int li_Theta1 ,li_Theta2, li_VM1, li_VM2;

  int Theta1_Theta1_Offset;
  int Theta1_Theta2_Offset;
  int Theta1_VM1_Offset;
  int Theta1_VM2_Offset;
 
  int Theta2_Theta1_Offset;
  int Theta2_Theta2_Offset;
  int Theta2_VM1_Offset;
  int Theta2_VM2_Offset;
    
  int VM1_Theta1_Offset;
  int VM1_Theta2_Offset;
  int VM1_VM1_Offset;
  int VM1_VM2_Offset;

  int VM2_Theta1_Offset;
  int VM2_Theta2_Offset;
  int VM2_VM1_Offset;
  int VM2_VM2_Offset;

  // only used for IV and PQ Real formulations
  double VR1, VR2, VI1, VI2;

  // only used for PQ Polar formulation
  double VM1, VM2, Theta1, Theta2;
  double dSin12, dSin21, dCos12, dCos21;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
//
//
// Special Notes :
// Creator        : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date  : 12/4/14
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

} // namespace PowerGridBusShunt
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_PowerGridBusShunt_h

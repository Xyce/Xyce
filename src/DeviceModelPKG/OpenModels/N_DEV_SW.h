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
// Purpose        : Switch base classes
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

#ifndef Xyce_N_DEV_SW_h
#define Xyce_N_DEV_SW_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace SW {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Controlled Switch";}
  static const char *deviceTypeName() {return "S level 1";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   SWiter,
     const FactoryBlock &      factory_block);


  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef );
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  bool processParams ();

  const std::vector<std::string> & getDepSolnVars();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  void setupPointers();

public:
  // iterator reference to the switch model which owns this instance
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  Util::Expression * Exp_ptr;
  int            expNumVars;
  int            expBaseVar;
  int            expNumDdt;
  std::list<std::string>   evnList;

  std::vector<double> expVarDerivs;
  std::vector<double> myVarVals;
  std::vector<double> ddtVals;
  double         expVal;


  // user specified parameters
  double R;     // resistance (ohms)
  double CONTROL;   // Value of control expression
  bool ON, OFF;      // whether switch is on or off initially

  // derived parameters
  double G;     // conductance (1.0/ohms)
  double dGdI,dGdV;

  // places to store node voltages
  double v_pos;
  double v_neg;

  double LeadCurrent;

  // double for current state of switch
  double SW_STATE;

  // and a state variable to save the SW_STATE
  double switch_state;

  std::vector<int>    li_ddt;

  int li_switch_state;

  // local indices (offsets)
  int li_Pos;
  int li_Neg;

  // store vector location for device lead current
  int li_branch_data;

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;

  // Offsets into the control nodes
  std::vector<int> APosEquControlNodeOffset;
  std::vector<int> ANegEquControlNodeOffset;

  // Ptr variables corresponding to the above declared indices.
  double * fPosEquPosNodePtr;
  double * fPosEquNegNodePtr;
  double * fNegEquPosNodePtr;
  double * fNegEquNegNodePtr;

  // Ptrs into the control nodes
  std::vector<double *> fPosEquControlNodePtr;
  std::vector<double *> fNegEquControlNodePtr;

  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &        MB,
     const FactoryBlock &      factory_block);

  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams ();
  bool processInstanceParams ();


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

  int dtype;      // device type: 1=SWITCH, 2=ISWITCH, 3=VSWITCH
  double VON;
  double VOFF;
  double ION;
  double IOFF;
  double RON;
  double ROFF;
  double ON;
  double OFF;
  double dInv;    // the inverse of (ON-OFF) or 1e-12, if too small.
  double Lm;      // log mean of resistances
  double Lr;      // log ratio of resistor values
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);
  virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

  // load functions, residual:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);

  // load functions, Jacobian:
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace SW
} // namespace Device
} // namespace Xyce

#endif

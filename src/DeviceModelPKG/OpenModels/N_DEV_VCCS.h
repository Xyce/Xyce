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
// Purpose        : Voltage  contolled current source classes.
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

#ifndef Xyce_N_DEV_VCCS_h
#define Xyce_N_DEV_VCCS_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace VCCS {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Linear Voltage Controlled Current Source";}
  static const char *deviceTypeName() {return "G level 1";}
  static int numNodes() {return 4;}
  static const char *primaryParameter() {return "T";}
  static bool isLinearDevice() {return true;}

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
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  bool isLinearDevice() const
  {
    if ( loadLeadCurrent )
    {
      return false;
    }
    return true;
  }

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef);
  void registerStateLIDs( const std::vector<int> & staLIDVecRef);
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef);
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;

  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool updateIntermediateVars () { return true; };
  bool updatePrimaryState ();

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  void setupPointers();

public:
  // iterator reference to the VCCS model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  double Transconductance;

  // Matrix equation index variables:
  // local indices (offsets)
  int li_Pos;
  int li_Neg;

  int li_ContPos;
  int li_ContNeg;

  int li_branch_data;   ///< Index for lead current and junction voltage (for power calculations)

  // Offset variables for all of the above index variables.
  int APosEquContPosVarOffset;
  int APosEquContNegVarOffset;
  int ANegEquContPosVarOffset;
  int ANegEquContNegVarOffset;

  double * f_PosEquContPosVarPtr;
  double * f_PosEquContNegVarPtr;
  double * f_NegEquContPosVarPtr;
  double * f_NegEquContNegVarPtr;
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
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1),
    separateInstances_(false)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec,
                               double * leadF, double * leadQ, double * junctionV)
  { return loadDAEVectors( solVec, fVec, qVec, bVec, leadF, leadQ, junctionV, ALL ); }
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec,
                               double * leadF, double * leadQ, double * junctionV, int loadType);

  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
  { return loadDAEMatrices( dFdx, dQdx, ALL ); }
  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType);

private:
  bool separateInstances_;
  InstanceVector      linearInstances_;            ///< List of owned linear voltage source instances
  InstanceVector      nonlinearInstances_;         ///< List of owned nonlinear voltage source instances
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace VCCS
} // namespace Device
} // namespace Xyce

#endif

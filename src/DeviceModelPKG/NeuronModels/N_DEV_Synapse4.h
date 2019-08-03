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
// Purpose        : Synapse4 classes
//
// Special Notes  :
//
// Creator        : Christy Warrender, SNL, Cognitive Modeling
//
// Creation Date  : 10/12/11
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Synapse4_h
#define Xyce_N_DEV_Synapse4_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Synapse.h>

#include <Sacado_No_Kokkos.hpp>

namespace Xyce {
namespace Device {
namespace Synapse4 {

// ---------- Forward Declarations -------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Synapse::Traits>
{
  static const char *name() {return "Synapse";}
  static const char *deviceTypeName() {return "YSYNAPSE level 4";}
  static int numNodes() {return 2;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : 4
// Purpose       :
//
//	This  is  the  instance class  for Synapse4s.  It
//	contains "unique" Synapse4  information - ie stuff that
//	will be true of only one  Synapse4 in the circuit, such
//	as the nodes to which it is connected.  A Synapse4 is
//	connected to only two circuit nodes.
//
//	This class  does not directly contain information about
//	its node indices. It contains indices into the 3 parts
//	(A, dx, and  b) of the matrix  problem A*dx = b, and
//	also x.  A is the Jacobian  matrix, dx is the update to
//	the solution vector x, and b is the right hand side
//	function vector.  These indices are global, and
//	determined by topology during  the initialization stage
//	of execution.
//
// Special Notes :
// Creator       : Christina Warrender, SNL,
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;

  // functions
public:
  static std::vector< std::vector<int> > jacStamp;

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Riter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStoreLIDs(const std::vector<int> & stoLIDVecRef );
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  bool processParams ();

  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  // enable the interface to produce plot files. - although we're not
  // actually using this for output
  bool plotfileFlag () {return true;}
  bool outputPlotFiles(bool force_final_output);

  void setupPointers();

  // is there currently a non-negligible synaptic current?
  // (used because there's no need to do calculations otherwise)
  bool active;

  // iterator reference to the Synapse4 model which owns this instance.
  Model &getModel() 
  {
    return model_;
  }

private:

  Model & model_;

  // user-specified parameters:
  double gMax;
  bool gMaxGiven;

  //Vector local index for Positive Node
  int li_Prev;
  //Vector local index for Negative Node
  int li_Post;

  // store vector quantities
  int li_A0_store;
  int li_B0_store;
  int li_t0_store;
  int li_branch_data;

#ifdef Xyce_FullSynapseJac
  // Offset variables corresponding to the above declared indices.
  int APostEquPostNodeOffset;

  // Pointers for Jacobian
  double *f_PostEquPostNodePtr;
#endif

  // vars used for load calculations
  double ipost;  // post Synapse4 current
  double didVpost;

  // flag to indicate random number generator was initialized
  bool randInitialized;

  bool ready;
  double respondTime;
};


//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
//
//
// Special Notes :
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
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

  // user-specified parameters
  double vThresh;
  double gMax;
  double delay;
  double eRev;
  double tau1;
  double tau2;
  double maxtau;

  // derived parameters
  double tp;		// time of EPSP peak, relative to start of postsynaptic response
  double factor;	// used to ensure peak conductance = 1S for weight (gMax) = 1
};

//-----------------------------------------------------------------------------
// Class         : Master4
// Purpose       :
// Special Notes :
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 07/16/12
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

  // load functions:
  virtual bool loadDAEVectors(double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Synapse4
} // namespace Device
} // namespace Xyce

#endif


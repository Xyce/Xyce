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
// Purpose        : Simple delay synapse
//
// Special Notes  :
//
// Creator        : Rich Schiek, SNL, Electrical Systems Modeling
//
// Creation Date  : 01/25/2011
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Synapse3_h
#define Xyce_N_DEV_Synapse3_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_UTL_RandomNumbers.h>

#include <N_DEV_Synapse.h>


#include <Sacado_No_Kokkos.hpp>

namespace Xyce {
namespace Device {
namespace Synapse3 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Synapse::Traits>
{
  static const char *name() {return "Synapse, Clopath";}
  static const char *deviceTypeName() {return "YSYNAPSE level 3";}
  static int numNodes() {return 2;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : 3
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

public:
  static std::vector< std::vector<int> > jacStamp;

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
  void registerStoreLIDs( const std::vector<int> & storeLIDVecRef );
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

  Model &getModel() 
  {
    return model_;
  }

  // make this static so that synapses share a random number generator and 
  // thus don't all seed the same generator with the same seed and then follow 
  // the same sequence.
  static Xyce::Util::RandomNumbers  * randomNumberGenerator_;
  
private:

  Model &       model_;         //< Owning model
  // make this public so that the master-class can potentially
  // access and output these values en-mass rather than
  // on a object by object basis.
  double synapticWeight;
  // user-specified parameters:
  double gMax;
  bool gMaxGiven;
  double transmissionProbability;
  bool transmissionProbabilityValueGiven;

  //Vector local index for Positive Node
  int li_Prev;
  //Vector local index for Negative Node
  int li_Post;

  // store variable quantities
  int li_A0_store;
  int li_B0_store;
  int li_t0_store;
  int li_weight_store;
  int li_VL1_store;
  int li_VL2_store;
  int li_VL3_store;
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

  bool ready;
  double respondTime;
  int transmissionFactor;  // 0 or 1 depending on random number < transmissionProbability

  // flag to indicate random number generator was initialized
  bool randInitialized;

  // initial weight can be specified on the instance line to override the model value
  double wInitialValue;   // initial value for weighting
  bool wInitialValueGiven; // flag set by parser if this value is given.

  double synapticWeightUpdate; // value calculated to update synaptic weight.
  double vl1Update;
  double vl2Update;
  double vl3Update;
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
  // parameters for weighting function
  double sParam;          // voltage threshold for a spike event
  double rParam;          // Resting voltage for resting event
  double wMin;            // minimum value for weighting
  double wMax;            // maximum value for weighting
  double wInitialValue;   // initial value for weighting
  double vL1tau1;         // rate for Longterm potentiation factor (LPF) based on post-synaptic voltage (rate 1)
  double vL2tau2;         // rate for Longterm potentiation factor (LPF) based on post-synaptic voltage (rate 2)
  double vL3tau3;         // rate for Longterm potentiation factor (LPF) based on pre-synaptic voltage (rate 3)
  double aLTD;            // long term depression coefficient
  double aLTP;            // long term potentiation coefficient
  double transmissionProbability;

  // derived parameters
  double tp;		// time of EPSP peak, relative to start of postsynaptic response
  double factor;	// used to ensure peak conductance = 1S for weight (gMax) = 1
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Synapse3
} // namespace Device
} // namespace Xyce

#endif

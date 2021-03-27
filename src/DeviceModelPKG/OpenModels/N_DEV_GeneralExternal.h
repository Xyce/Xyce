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

//-----------------------------------------------------------------------------
//
// Purpose        : General External Device classes.
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date  : 03/01/2017
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_GeneralExternal_h
#define Xyce_N_DEV_GeneralExternal_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_VectorComputeInterface.h>
#include <N_DEV_CompositeParam.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : GenExtDoubleData
// Purpose       : This is class is a CompositeParameter type for managing
//                 vector-composite parameters for doubles
// Special Notes : DPARAMS is the YGenExt vector-composite parameter
//                 allowing a netlist to specify name/value pairs that can
//                 be accessed by the coupled simulator.
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
class GenExtDoubleData : public CompositeParam
{
  friend class ParametricData<GenExtDoubleData>;

 public:
  static ParametricData<GenExtDoubleData> &getParametricData();
  GenExtDoubleData();
  void processParams();
  friend std::ostream & operator<<(std::ostream & os, const GenExtDoubleData & gEdd);

 private:
  std::string name_;
  double value_;

 public:
  std::string getName() const {return name_;};
  double getValue() const {return value_;};
};

//-----------------------------------------------------------------------------
// Class         : GenExtIntData
// Purpose       : This is class is a CompositeParameter type for managing
//                 vector-composite parameters for ints
// Special Notes : IPARAMS is the YGenExt vector-composite parameter
//                 allowing a netlist to specify name/value pairs that can
//                 be accessed by the coupled simulator.
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
class GenExtIntData : public CompositeParam
{
  friend class ParametricData<GenExtIntData>;

 public:
  static ParametricData<GenExtIntData> &getParametricData();
  GenExtIntData();
  void processParams();
  friend std::ostream & operator<<(std::ostream & os, const GenExtIntData & gEid);

 private:
  std::string name_;
  int value_;

 public:
  std::string getName() const {return name_;};
  int getValue() const {return value_;};
};

//-----------------------------------------------------------------------------
// Class         : GenExtBoolData
// Purpose       : This is class is a CompositeParameter type for managing
//                 vector-composite parameters for bools
// Special Notes : IPARAMS is the YGenExt vector-composite parameter
//                 allowing a netlist to specify name/value pairs that can
//                 be accessed by the coupled simulator.
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
class GenExtBoolData : public CompositeParam
{
  friend class ParametricData<GenExtBoolData>;

 public:
  static ParametricData<GenExtBoolData> &getParametricData();
  GenExtBoolData();
  void processParams();
  friend std::ostream & operator<<(std::ostream & os, const GenExtBoolData & gEid);

 private:
  std::string name_;
  bool value_;

 public:
  std::string getName() const {return name_;};
  bool getValue() const {return value_;};
};

//-----------------------------------------------------------------------------
// Class         : GenExtStringData
// Purpose       : This is class is a CompositeParameter type for managing
//                 vector-composite parameters for strings
// Special Notes : IPARAMS is the YGenExt vector-composite parameter
//                 allowing a netlist to specify name/value pairs that can
//                 be accessed by the coupled simulator.
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
class GenExtStringData : public CompositeParam
{
  friend class ParametricData<GenExtStringData>;

 public:
  static ParametricData<GenExtStringData> &getParametricData();
  GenExtStringData();
  void processParams();
  friend std::ostream & operator<<(std::ostream & os, const GenExtStringData & gEid);

 private:
  std::string name_;
  std::string value_;

 public:
  std::string getName() const {return name_;};
  std::string getValue() const {return value_;};
};


namespace GeneralExternal {
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "GeneralExternal";}
  static const char *deviceTypeName() {return "GeneralExternal device level 1";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 1000;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 GeneralExternal device.
// Special Notes :
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//-----------------------------------------------------------------------------
/// Instance class for General External device
///
/// The General External device is basically a stub device that depends entirely
/// on an interface object derived from VectorComputeInterface to do the
/// device computations.  An external code implements the VectorComputeInterface,
/// and associates it with a specific GeneralExteranal device using the
/// "setVectorLoader" method of the GenCouplingSimulator class.
/// The device information may also be modified using the "setNumInternalVars"
/// and "setJacStamp" methods of GenCouplingSimulator.
///
/// @see src/test/GenExtTestHarnesses/testGenCoup.C  for example usage.
///
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Master;
  friend struct Traits;

public:
  Instance(
     const Configuration &       configuration,
     const InstanceBlock &            IB,
     Model & Miter,
     const FactoryBlock &factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef );
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const;

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  int getNumVars();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();
  bool loadDAEBVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  bool setVectorLoader(Xyce::Device::VectorComputeInterface * vciPtr);
  void setNumInternalVars(int numInt);
  void setNumStateVars(int numState);
  void setNumBranchDataVarsIfAllocated(int numBranchDataIfAllocated);
  void setJacStamp(std::vector< std::vector<int> > & jS);
  void getSolution(std::vector<double> &sV);

  CompositeParam *constructComposite (const std::string &, const std::string &);

  void getDParams(std::vector<std::string> &names,std::vector<double> &values);
  void getIParams(std::vector<std::string> &names,std::vector<int> &values);
  void getBParams(std::vector<std::string> &names,std::vector<bool> &values);
  void getSParams(std::vector<std::string> &names,std::vector<std::string> &values);

protected:
private:
  void setupDenseJacStamp_();

public:
  // iterator reference to the general external model which owns this instance.
  // Getters and setters
  Model &getModel()
  {
    return model_;
  }


private:

  Model &       model_;         //< Owning model
  std::map<std::string, GenExtDoubleData *> doubleDataMap_;
  std::map<std::string, GenExtIntData *> intDataMap_;
  std::map<std::string, GenExtBoolData *> boolDataMap_;
  std::map<std::string, GenExtStringData *> stringDataMap_;


private:
  // parameter variables
  // For vector composite:
  std::vector<GenExtDoubleData*> doubleDataVec_;
  std::vector<GenExtIntData*> intDataVec_;
  std::vector<GenExtBoolData*> boolDataVec_;
  std::vector<GenExtStringData*> stringDataVec_;

  // state variables
  // This device has no state

  // local state indices (offsets)
  // This device has no state

  // local solution indices (offsets)
  // This device uses an array of li_ values instead of individually named
  // variables.
  std::vector<int> li_Nodes_;
  std::vector<int> li_States_;
  std::vector<int> li_Stores_;
  // these are the LIDs for lead currents
  std::vector<int> li_Branches_;

  // Matrix equation index variables:

  // Offset variables.  Again, this device uses an array instead of
  // discrete variables.
  // A_Equ_NodeOffests[equation][node] is the offset for node in
  // equation
  std::vector< std::vector<int> > A_Equ_NodeOffsets_;

  std::vector< std::vector<int> > jacStamp_;

  Xyce::Device::VectorComputeInterface * vciPtr_;

  /// Vectors and matrices for storing local representations of F,
  /// Q, and B vectors, used for passing to the computeXyceVectors
  /// function of the *vciPtr_ object.
  std::vector<double> solutionVars_;
  std::vector<double> flagSolutionVars_;
  std::vector<double> nextStoreVars_;
  std::vector<double> currStoreVars_;
  std::vector<double> FVec_;
  std::vector<double> QVec_;
  std::vector<double> BVec_;
  std::vector<double> dFdXdVpVec_;
  std::vector<double> dQdXdVpVec_;
  std::vector<std::vector<double> > dFdXMat_;
  std::vector<std::vector<double> > dQdXMat_;
  std::vector<std::vector<double> > storeVars;

  std::vector<std::complex<double> > solutionFDVars_;
  std::vector<std::complex<double> >fDFVec_;
  std::vector<std::complex<double> >fDBVec_;
  std::vector<std::vector<std::complex<double> > > dFdXFDMat_;

  // Storage for lead current components
  std::vector<double> leadCurrentF;
  std::vector<double> leadCurrentQ;


  // Flag to tell us whether the vector loader object implements frequency
  // domain loads
  bool haveFDLoads_;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/17/2017
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Master;
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

  // Additional Implementation Declarations
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       : This is class refers to the collected instances of the
//                 GeneralExternal device.
// Special Notes :
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 20 March 2018
//-----------------------------------------------------------------------------

class Master : public DeviceMaster<Traits>
{
  friend class Instance;                            ///< Don't force a lot of pointless getters
  friend class Model;                               ///< Don't force a lot of pointless getters

public:

  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
   : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec,
                            int loadType);

  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV, int loadType);

  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx,
                               int loadType);

  virtual bool updateFDIntermediateVars(double frequency,
                                        std::complex<double> * solVec);

  virtual bool loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec);
  virtual bool loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx);

};

//----------------------------------------------------------------------------
// Function      : Instance::setVectorLoader
// Purpose       : Associate a VectorCompute object with this device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 2/28/2017
//----------------------------------------------------------------------------
/// associate a vector loader object pointer with this instance
inline bool Instance::setVectorLoader(Xyce::Device::VectorComputeInterface * vciPtr)
{
  bool bsuccess=true;
  // Because we did not do this in the constructor as is typical for
  // most devices, make sure we get it done now.
  if (jacStamp_.empty())
      setupDenseJacStamp_();

  if (vciPtr)
  {
    vciPtr_=vciPtr;
    haveFDLoads_ = vciPtr_->haveFDLoads();
  }
  else
  {
    bsuccess=false;
    haveFDLoads_=false;
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function      : Instance::setNumInternalVars
// Purpose       : Set the number of internal variables for this device
// Special Notes : For God's sake, make sure to call this before
//                 topology needs to query the number!
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 3/1/2017
//----------------------------------------------------------------------------
inline void Instance::setNumInternalVars(int numInt)
{
  numIntVars=numInt;
}

//----------------------------------------------------------------------------
// Function      : Instance::setNumStateVars
// Purpose       : Set the number of internal variables for this device
// Special Notes : For God's sake, make sure to call this before
//                 topology needs to query the number!
// Scope         : public
// Creator       : Paul Kuberry, SNL
// Creation Date : 12/10/2020
//----------------------------------------------------------------------------
inline void Instance::setNumStateVars(int numState)
{
  numStateVars=numState;
}

//----------------------------------------------------------------------------
// Function      : Instance::setNumBranchDataIfAllocatedVars
// Purpose       : Set the number of internal variables for this device
// Special Notes : For God's sake, make sure to call this before
//                 topology needs to query the number!
// Scope         : public
// Creator       : Paul Kuberry, SNL
// Creation Date : 12/10/2020
//----------------------------------------------------------------------------
inline void Instance::setNumBranchDataVarsIfAllocated(int numBranchDataVarsIfAllocated)
{
  numBranchDataVarsIfAllocated=numBranchDataVarsIfAllocated;
}

//----------------------------------------------------------------------------
// Function      : Instance::setJacStamp
// Purpose       : Set the jacobian stamp
// Special Notes : For God's sake, make sure to call this before
//                 topology needs to query the jacStamp!
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 3/1/2017
//----------------------------------------------------------------------------
///
/// Set the jacobian stamp
///
/// The Jacobian Stamp is a sparse representation of the pattern of non-zero
/// elements of the jacobian associated with a device.
///
/// The jacobian stamp must have as many rows as the device has variables
/// (i.e. numIntVars+numExtVars).
///
/// Each row of the jacStamp is associated with an equation for the device.
/// The number of elements of a row is the number of non-zero Jacobian
/// entries for that equation.
///
/// The values of the elements in a row are the variable numbers for the
/// non-zero elements of the Jacobian.
///
/// Currently, Xyce uses the same jacStamp for both dFdX and dQdX, so
/// the jacobian stamp must be the union of nonzeros for these two
/// matrices.
///
/// A jacobian stamp row may be empty, indicating that the F and Q vector
/// elements for that row have no dependence on any variable, or the vector
/// elements are always zero (as would be the case for control nodes of a
/// controlled source, for example).
///
/// If this function is never called, a fully dense jacobian stamp
/// (in which every equation may depend on every variable) is constructed
/// at the time that setVectorLoader is called.
inline void Instance::setJacStamp(  std::vector< std::vector<int> > & jS)
{
  jacStamp_ = jS;
}

//----------------------------------------------------------------------------
// Function      : Instance::getNumVars
// Purpose       : Return the number of solution vars in a given instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical Models and Simulation
// Creation Date : 08/27/2008
//----------------------------------------------------------------------------
inline int Instance::getNumVars()
{
  return numExtVars+numIntVars;
}

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Resistor
} // namespace Device
} // namespace Xyce

#endif

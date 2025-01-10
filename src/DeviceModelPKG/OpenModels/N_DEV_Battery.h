//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Implementation of a non-ideal voltage source
//                  
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Models & Simulations
//
// Creation Date  : 10/28/2016
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Battery_h
#define Xyce_N_DEV_Battery_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>


namespace Xyce {
namespace Device {
namespace Battery {

class Model;
class Instance;

// sensitivity functor
// not yet implemented.
class BatterySensitivity :  public baseSensitivity
{
  public:
  BatterySensitivity() : 
    baseSensitivity() {};

  virtual ~BatterySensitivity() {};

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &name,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const ;
};

static BatterySensitivity memrSens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Battery";}
  static const char *deviceTypeName() {return "YBATTERY level 1";}
  static int numNodes() {return 3;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &p);
  static void loadInstanceParameters(ParametricData<Instance> &p);
};

//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Battery::Instance
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Battery device instance class.
//
// An instance is created for each occurance of the device in the netlist.
//
// It contains "unique" device information - ie stuff that will be
// true of only one battery in the circuit, such as the nodes to
// which it is connected.  A battery is connected to only two
// circuit nodes with a third node for the temperature
//
// This class does not directly contain information about its node
// indices. It contains indices into the 5 parts (dFdx, dQdx, dx, F,
// and Q) of the matrix problem A*dx = b, and also the solution
// vector x.  A is the Jacobian matrix that will be formed from dFdx
// and d(dQ/dt)dx, dx is the update to x, and b is right hand side
// function vector that will be formed from F and dQ/dt.  These
// indices are global, and determined by topology during the
// initialization stage of execution.
//
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;              
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class BatterySensitivity;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Battery::Instance::~Instance
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Destroys this instance
  //
  // @author Eric Keiter, SNL
  // @date   3/16/00
  ~Instance() {}

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Battery::Instance::getModel
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Gets the resistor model that owns this instance.
  //
  // @return reference to the owning Battery::Model
  //
  // @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  // @date   Mon Aug 12 08:36:37 2013
  Model &getModel()   { return model_;  }

  virtual void registerLIDs(const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef) /* override */;
  virtual void registerStateLIDs(const std::vector<int> & staLIDVecRef) /* override */;
  virtual void registerStoreLIDs(const std::vector<int> & stoLIDVecRef) /* override */;
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef) /* override */;
  virtual void registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec) /* override */;

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  virtual bool processParams() /* override */;
  virtual bool updateTemperature(const double & temp_tmp) /* override */;
  virtual bool updateIntermediateVars() /* override */;
  virtual bool updatePrimaryState() /* override */;

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Battery::Instance::jacobianStamp
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Return Jacobian stamp that informs topology of the layout of the
  // resistor jacobian.
  //
  // The Jacobian stamp describes the shape of the Jacobian to the
  // Topology subsystem.  The Topology subsystem, in turn, returns
  // the offsets into the matrix and solution vectors where this
  // instance data is located.
  //
  // @return const reference to a std::vector of std::vector of
  // integers describing Jacobian stamp shape
  //
  // @author Robert Hoekstra
  // @date 8/20/2001
  virtual const std::vector< std::vector<int> > &jacobianStamp() const  /* override */ {
    return jacStamp;
  }

  virtual bool loadDAEQVector() /* override */;
  virtual bool loadDAEFVector() /* override */;
  virtual bool loadDAEdQdx() /* override */;
  virtual bool loadDAEdFdx() /* override */;


  virtual void setupPointers() /* override */;

private:
  static std::vector< std::vector<int> >  jacStamp; //< All Battery instances have a common Jacobian Stamp
  static void initializeJacobianStamp();

  Model &     model_;                 //< Owning model

  // User-specified parameters:

  // Derived parameters:
  double      G;                      //< Conductance(1.0/ohms)
  double      Rcell;                  //< Cell resistance 
  double      i0;                     //< branch current

  // F/Q load values and jacboian load values
  double      v1FEqu;
  double      dv1dv1FEqu, dv2dv1FEqu, dCellTempdv1FEqu, dIdv1Equ, dCapUseddv1Equ;
  double      v2Equ;
  double      dv1dv2FEqu, dv2dv2FEqu, dCellTempdv2FEqu, dIdv2Equ, dCapUseddv2Equ;
  double      cellTempEqu, dCellTempdCellTempEqu;
  double      currentEqu, dV1dCurrentEqu, dV2dCurrentEqu;
  double      capUsedEqu, dCurrentdCapUsedEqu;
  
  // indices for array loads
  int         li_Pos;                 //< Index for Positive Node
  int         li_Neg;                 //< Index for Negative Node
  int         li_CellTemp;            //< Index for cell temperature 
  int         li_current;             //< Index for branch current
  int         li_capUsed;             //< Index for internal x, thickness, variable
  int         li_branch_data;         //< Index for Lead Current and junction voltage (for power calculations)

  // Offset variables corresponding to the above declared indices.
  int         APosEquPosNodeOffset;
  int         APosEquNegNodeOffset;     
  int         APosEquCellTempOffset; 
  int         APosEquCurrentOffset; 
  int         APosEquCapUsedNodeOffset; 
  
  int         ANegEquPosNodeOffset;     
  int         ANegEquNegNodeOffset;     
  int         ANegEquCellTempOffset;
  int         ANegEquCurrentOffset;     
  int         ANegEquCapUsedNodeOffset; 
  
  int         CellTempEquCellTempOffset; 
  
  int         CurrentEquPosNodeOffset;
  int         CurrentEquNegNodeOffset;
  int         CurrentEquCurrentOffset;
  
  int         CapUsedEquCurrentOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double *    f_PosEquPosNodePtr;
  double *    f_PosEquNegNodePtr;
  double *    f_PosEquCellTempPtr;
  double *    f_PosEquCurrentPtr;
  double *    f_PosEquCapUsedNodePtr;
  
  double *    f_NegEquPosNodePtr;
  double *    f_NegEquNegNodePtr;
  double *    f_NegEquCellTempPtr;
  double *    f_NegEquCurrentPtr;
  double *    f_NegEquCapUsedNodePtr;
  
  double *    f_CellTempEquCellTempPtr;
  
  double *    f_CurrentEquPosNodePtr;
  double *    f_CurrentEquNegNodePtr;
  double *    f_CurrentEquCurrentPtr;
    
  double *    f_CapUsedEquCurrentPtr;
  double *    q_CapUsedEquCurrentPtr;

#endif

};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Battery::Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// Battery model class
//
class Model : public DeviceModel
{
  friend class ParametricData<Model>;               //< Allow ParametricData to changes member values
  friend class Instance;                            //< Don't force a lot of pointless getters
  friend struct Traits;
  friend class Master;                              //< Don't force a lot of pointless getters

public:
  typedef std::vector<Instance *> InstanceVector;

  Model(
     const Configuration &     configuration,
     const ModelBlock &        model_block,
     const FactoryBlock &      factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Battery::Model::addInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014 
  //---------------------------------------------------------------------------
  //
  // Add an instance to the list of instances associated with this model
  //
  // @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
  // @date   8/12/2013
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

  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

  virtual bool processParams() /* override */;
  virtual bool processInstanceParams() /* override */;

private:
  InstanceVector      instanceContainer;            //< List of owned resistor instances

  // model parameters for Yakopcic model
  double      VCapOffset_;
  double      VCap0_;
  double      VCap1_;
  double      VCap2_;
  double      VCap3_;
  
  double      IC0_;
  double      IComp0_;
  double      IComp1_;
  
  double      Tnom_;
  double      TComp0_;
  double      TComp1_;
  
  double      IF0_;
  double      IFact0_;
  double      IFact1_;
  
  double      TFactNom_;
  double      TFact0_;
  double      TFact1_;
};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Battery::Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014 
//-----------------------------------------------------------------------------
//
// Battery master
//
// The "master" class is the one that contains the updateState, loadDAEVectors
// and loadDAEMatrices methods that are actually called when it is time to
// compute and load device contributions.
//
// The default implementations of these methods in the DeviceMaster
// template class simply loops over all instances and calls their
// updatePrimaryState, loadDAEFVector/loadDAEQVector, and
// loadDAEdFdx/loadDAEdQdx methods, respectively.
//
// For efficiency, the Battery class reimplements these methods to do the
// work directly, instead of calling instance-level functions.
//
class Master : public DeviceMaster<Traits>
{
  friend class Instance;                           
  friend class Model;                             

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Battery::Master::Master
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Construct a Battery Device.
  //
  // @param configuration
  // @param factory_block
  // @param solver_state
  // @param device_options
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options)
  {}

  virtual bool updateState(double * solVec, double * staVec, double * stoVec);
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx);

};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Battery
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Battery_h

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
// Purpose        : Implementation of the Yakopcic memristor model.  See.
//                  
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Models & Simulations
//
// Creation Date  : 10/23/2014
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MemristorYakopcic_h
#define Xyce_N_DEV_MemristorYakopcic_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_MemristorTEAM.h>

#include <random>

namespace Xyce {
namespace Device {
namespace MemristorYakopcic {

class Model;
class Instance;

// sensitivity functor
// not yet implemented.
class memristorYakopcicSensitivity :  public baseSensitivity
{
  public:
  memristorYakopcicSensitivity() : 
    baseSensitivity() {};

  virtual ~memristorYakopcicSensitivity() {};

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

static memristorYakopcicSensitivity memrSens;

struct Traits : public DeviceTraits<Model, Instance, MemristorTEAM::Traits>
{
  static const char *name() {return "MemristorYakopcic";}
  static const char *deviceTypeName() {return "YMEMRISTOR level 3";}
  static int numNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &p);
  static void loadInstanceParameters(ParametricData<Instance> &p);
};

//-----------------------------------------------------------------------------
// Class         : Xyce::Device::MemristorYakopcic::Instance
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// MemristorYakopcic device instance class.
//
// An instance is created for each occurance of the device in the netlist.
//
// It contains "unique" device information - ie stuff that will be
// true of only one memristor in the circuit, such as the nodes to
// which it is connected.  A memristor is connected to only two
// circuit nodes.
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
  friend class memristorYakopcicSensitivity;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::MemristorYakopcic::Instance::~Instance
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
  // Function      : Xyce::Device::MemristorYakopcic::Instance::getModel
  // Purpose       : destructor
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Gets the resistor model that owns this instance.
  //
  // @return reference to the owning MemristorYakopcic::Model
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
  // Function      : Xyce::Device::MemristorYakopcic::Instance::jacobianStamp
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
  static std::vector< std::vector<int> >  jacStamp; //< All MemristorYakopcic instances have a common Jacobian Stamp
  static void initializeJacobianStamp();

  Model &     model_;                 //< Owning model

  // User-specified parameters:
  double      XO_; 

  // Derived parameters:
  double      G;                      //< Conductance(1.0/ohms)
  double      Reff;                   //< Effective resistance 
  double      dReffdvpos;             //< derivative of Reff with respect to Vpos
  double      dReffdvneg;             //< derivative of Reff with respect to Vneg
  double      dReffdx;                //< derivative of Reff with respect to x  
  double      dIdx;                //< derivative of I with respect to x  
  double      i0;                     //< Current(ohms)
  double      xVarFContribution;      //< x, internal variable for thickness of conductive layer, F-vector contribution
  double      dxFEqdVpos;             //< derivative of X F equation with respect to Vpos
  double      dxFEqdVneg;             //< derivative of X F equation with respect to Vneg
  double      dxFEqdx;                //< derivative of X F equation with respect to X
  int         resNoiseLastUpdateStep; //< The last step on which the random noise was udpated.  
  double      resNoiseLastUpdateTime; //< The last time when the noise term in the resistence was updated.
  double      resNoiseNextUpdateTime; //< The next time to update the resistance noise.
  int         resNoiseHiLoState;      //< Flag indicating if we are currently in the low or high RTN state
  double      resNoiseRfactor;        //< factor that changes the effective Ron/Roff until the next random update
  int         xNoiseLastUpdateStep; //< The last step on which the random noise was udpated.  
  double      xNoiseLastUpdateTime; //< The last time when the noise term in the Xistence was updated.
  double      xNoiseNextUpdateTime; //< The next time to update the Xistance noise.
  int         xNoiseHiLoState;      //< Flag indicating if we are currently in the low or high RTN state
  double      xNoiseFactor;        //< factor that changes the effective growth rate until the next random update

  int         li_Pos;                 //< Index for Positive Node
  int         li_Neg;                 //< Index for Negative Node
  int         li_x;                   //< Index for internal x, thickness, variable
  int         li_store_R;             //< Index to store resistence value
  int         li_store_tdt;           //< Index to store for next RTN time time delta t 
  int         li_branch_data;         //< Index for Lead Current and junction voltage (for power calculations)

  // Offset variables corresponding to the above declared indices.
  int         APosEquPosNodeOffset;   //< Column index into matrix of Pos/Pos conductance
  int         APosEquNegNodeOffset;   //< Column index into matrix of Pos/Neg conductance
  int         APosEquXNodeOffset;     //< Column index into matrix for internal varaible x, layer thickness
  int         ANegEquPosNodeOffset;   //< Column index into matrix of Neg/Pos conductance
  int         ANegEquNegNodeOffset;   //< Column index into matrix of Neg/Neg conductance
  int         ANegEquXNodeOffset;     //< Column index into matrix for internal varaible x, layer thickness
  int         XEquVPosOffset;         //< Thickness governing equation, VPos dependence
  int         XEquVNegOffset;         //< Thickness governing equation, VNeg dependence
  int         XEquXOffset;            //< Thickness variable, in thickness governing equation equation

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double *    f_PosEquPosNodePtr;
  double *    f_PosEquNegNodePtr;
  double *    f_PosEquXNodePtr;
  double *    f_NegEquPosNodePtr;
  double *    f_NegEquNegNodePtr;
  double *    f_NegEquXNodePtr;
  double *    f_XEquPosNodePtr;
  double *    f_XEquNegNodePtr;
  double *    f_XEquXNodePtr;
  double *    q_XEquXNodePtr;
#
#endif

};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::MemristorYakopcic::Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014
//-----------------------------------------------------------------------------
//
// MemristorYakopcic model class
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
  // Function      : Xyce::Device::MemristorYakopcic::Model::addInstance
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
  double      Eta_;
  double      Vp_;
  double      Vn_;
  double      Ap_;
  double      An_;
  double      A1_;
  double      A2_;
  double      B_;
  double      AlphaP_;
  double      AlphaN_;
  double      XP_;
  double      XN_;
  
  
  // old model parameters for TEAM model

  double      xScaling_;
  bool        randomResNoiseOn_;
  int         randomResNoiseSeed_;
  double      randomResNoiseLambda_;
  double      randomResNoiseMean_;
  double      randomResNoiseSD_;
  double      randomResUpdateTime_;
  double      randomResEpsilonUpdateTime_;
  double      randomResDelta_;
  double      randomResDeltaGrad_;
  bool        randomXNoiseOn_;
  int         randomXNoiseSeed_;
  double      randomXNoiseLambda_;
  double      randomXNoiseMean_;
  double      randomXNoiseSD_;
  double      randomXUpdateTime_;
  double      randomXEpsilonUpdateTime_;
  double      randomXDelta_;
  double      randomXDeltaGrad_;
  std::mt19937 * randomNumberGenPtr_;
  std::uniform_real_distribution<double> * uniformRandomPtr_;
  std::normal_distribution<double> * gaussianRandomPtr_;
  
};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::MemristorYakopcic::Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical Models & Simulations
// Creation Date : 10/23/2014 
//-----------------------------------------------------------------------------
//
// MemristorYakopcic master
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
// For efficiency, the MemristorYakopcic class reimplements these methods to do the
// work directly, instead of calling instance-level functions.
//
class Master : public DeviceMaster<Traits>
{
  friend class Instance;                           
  friend class Model;                             

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::MemristorYakopcic::Master::Master
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Richard Schiek, Electrical Models & Simulations
  // Creation Date : 10/23/2014
  //---------------------------------------------------------------------------
  //
  // Construct a MemristorYakopcic Device.
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

} // namespace MemristorYakopcic
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_MemristorYakopcic_h

//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        :
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

#ifndef Xyce_N_DEV_Resistor_h
#define Xyce_N_DEV_Resistor_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceMaster.h>

namespace Xyce {
namespace Device {
namespace Resistor {

class Model;
class Instance;

/// sensitivity functor
class resistorSensitivity :  public baseSensitivity
{
  public:
  resistorSensitivity() :
    baseSensitivity() {};

  virtual ~resistorSensitivity() {};

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

static resistorSensitivity resSens;

/// matrix sensitivity functor
class resistorMatrixSensitivity :  public baseMatrixSensitivity
{
  public:
  resistorMatrixSensitivity() :
    baseMatrixSensitivity() {};

  virtual ~resistorMatrixSensitivity() {};

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLids,
    std::vector< std::vector<int> > & Q_jacLids
    ) const ;
};

static resistorMatrixSensitivity resMatrixSens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Resistor";}
  static const char *deviceTypeName() {return "R level 1";}
  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "R";}
  static const char *instanceDefaultParameter() {return "R";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &p);
  static void loadInstanceParameters(ParametricData<Instance> &p);
};

//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Resistor::Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Resistor device instance class.
///
/// An instance is created for each occurance of the device in the netlist.
///
/// It contains "unique" resistor information - ie stuff that will be
/// true of only one resistor in the circuit, such as the nodes to
/// which it is connected.  A resistor is connected to only two
/// circuit nodes.
///
/// This class does not directly contain information about its node
/// indices. It contains indices into the 5 parts (dFdx, dQdx, dx, F,
/// and Q) of the matrix problem A*dx = b, and also the solution
/// vector x.  A is the Jacobian matrix that will be formed from dFdx
/// and d(dQ/dt)dx, dx is the update to x, and b is right hand side
/// function vector that will be formed from F and dQ/dt.  These
/// indices are global, and determined by topology during the
/// initialization stage of execution.
///
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;              ///< Allow ParametricData to changes member values
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class resistorSensitivity;
  friend class resistorMatrixSensitivity;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  /// Destroys this instance
  ~Instance() {}

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:

  bool isLinearDevice() const;

  /// Gets the resistor model that owns this instance.
  Model &getModel()   { return model_;  }

  virtual void registerLIDs(const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef) /* override */;
  virtual void registerStateLIDs(const std::vector<int> & staLIDVecRef) /* override */;
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef) /* override */;
  virtual void registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec) /* override */;

  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  virtual bool processParams() /* override */;
  virtual bool updateTemperature(const double & temp_tmp) /* override */;
  virtual bool updateIntermediateVars() { return true; }
  virtual bool updatePrimaryState() { return true; }

  int getNumNoiseSources () const { return 1; }
  void setupNoiseSources (Xyce::Analysis::NoiseData & noiseData);
  void getNoiseSources (Xyce::Analysis::NoiseData & noiseData);

  /// Return Jacobian stamp that informs topology of the layout of the
  /// resistor jacobian.
  virtual const std::vector< std::vector<int> > &jacobianStamp() const  /* override */ 
  {
    if (solVarDep) { return jacStamp_solVarDep; }
    else { return jacStamp; }
  }

  virtual bool loadDAEFVector() /* override */;
  virtual bool loadDAEdFdx() /* override */;

  /// Load Q vector
  /// Since the Resistor does no charge storage, this is a no-op.
  virtual bool loadDAEQVector()
  {
    return true;
  }

  /// Load derivative of Q vector with respect to solution vector
  /// Since the Resistor does no charge storage, this is a no-op.
  virtual bool loadDAEdQdx()
  {
    return true;
  }

  virtual void setupPointers() /* override */;

private:
  static std::vector< std::vector<int> >  jacStamp; ///< All Resistor instances have a common Jacobian Stamp, except when solution dependent
  static void initializeJacobianStamp();
  
  std::vector< std::vector<int> >  jacStamp_solVarDep; 

  Model &     model_;                 ///< Owning model

  // Stuff for handling solution-variable-dependent capacitance
  Util::Expression * expPtr;
  int                expNumVars;
  std::vector<double> expVarDerivs;
  bool solVarDep;

  // User-specified parameters:
  double      R;                      ///< Resistance (ohms)
  double      multiplicityFactor;     ///< multiplicity factor (M)
  double      factor;                 ///< computed in  updateTemperature

  // These are for the semiconductor resistor
  double      length;                 ///< Resistor length.
  double      width;                  ///< Resistor width.
  double      temp;                   ///< Temperature of this instance

  // Temperature dependence parameters, these can override values specified in the model
  double      tempCoeff1;             ///< First order temperature coeff.
  double      tempCoeff2;             ///< Second order temperature coeff.
  double      tempCoeffExp;           ///< Exponential temperature coeff.
  double      dtemp;                  ///< Externally specified device temperature.
  ///<   NOT used, only here for compatibility in parsing
  ///<   netlist from simulators that support it

  // Flags used to tell if the user has specified one of these values on the command line.
  bool        tempCoeff1Given;        ///< First order temperature coeff was given in netlist
  bool        tempCoeff2Given;        ///< Second order temperature coeff was given in netlist
  bool        tempCoeffExpGiven;      ///< Exponential temperature coeff was given in netlist
  bool        dtempGiven;             ///< Externally specified device temperature was given in netlist

  // Derived parameters:
  double      G;                      ///< Conductance(1.0/ohms)
  double      i0;                     ///< Current(ohms)

  int         li_Pos;                 ///< Index for Positive Node
  int         li_Neg;                 ///< Index for Negative Node
  int         li_branch_data;         ///< Index for Lead Current and junction voltage (for power calculations)

  // Offset variables corresponding to the above declared indices.
  int         APosEquPosNodeOffset;   ///< Column index into force matrix of Pos/Pos conductance
  int         APosEquNegNodeOffset;   ///< Column index into force matrix of Pos/Neg conductance
  int         ANegEquPosNodeOffset;   ///< Column index into force matrix of Neg/Pos conductance
  int         ANegEquNegNodeOffset;   ///< Column index into force matrix of Neg/Neg conductance

  // Offsets into the control nodes
  std::vector<int> APosEquControlNodeOffset;
  std::vector<int> ANegEquControlNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double *    f_PosEquPosNodePtr;
  double *    f_PosEquNegNodePtr;
  double *    f_NegEquPosNodePtr;
  double *    f_NegEquNegNodePtr;

  // Ptrs into the control nodes
  std::vector<double *> fPosEquControlNodePtr;
  std::vector<double *> fNegEquControlNodePtr;
#endif

};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Resistor::Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Resistor model class
///
class Model : public DeviceModel
{
  friend class ParametricData<Model>;               ///< Allow ParametricData to changes member values
  friend class Instance;                            ///< Don't force a lot of pointless getters
  friend struct Traits;
  friend class Master;                              ///< Don't force a lot of pointless getters

public:
  typedef std::vector<Instance *> InstanceVector;

  Model(
     const Configuration &       configuration,
     const ModelBlock &        model_block,
     const FactoryBlock &      factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Model::addInstance
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David Baur
  // Creation Date : 8/12/2013
  //---------------------------------------------------------------------------
  ///
  /// Add an instance to the list of instances associated with this model
  ///
  /// @author David G. Baur  Raytheon  Sandia National Laboratories 1355
  /// @date   8/12/2013
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
  InstanceVector      instanceContainer;            ///< List of owned resistor instances

  // Semiconductor resistor parameters
  double      tempCoeff1;     ///< First order temperature coefficient
  double      tempCoeff2;     ///< Second order temperature coefficient
  double      tempCoeffExp;   ///< Exponential temperature coefficient
  bool        tempCoeffExpModelGiven; ///< Exponential temperature coeff given in model
  double      sheetRes;       ///< Sheet resistance
  double      resistanceMultiplier;       ///<  resistance multiplier
  double      defWidth;       ///< Default width
  double      narrow;         ///< Narrowing due to side etching
  double      tnom;           ///< Parameter measurement temperature
};


//-----------------------------------------------------------------------------
// Class         : Xyce::Device::Resistor::Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Resistor master
///
/// The "master" class is the one that contains the updateState, loadDAEVectors
/// and loadDAEMatrices methods that are actually called when it is time to
/// compute and load device contributions.
///
/// The default implementations of these methods in the DeviceMaster
/// template class simply loops over all instances and calls their
/// updatePrimaryState, loadDAEFVector/loadDAEQVector, and
/// loadDAEdFdx/loadDAEdQdx methods, respectively.
///
/// For efficiency, the Resistor class reimplements these methods to do the
/// work directly, instead of calling instance-level functions.
///
class Master : public DeviceMaster<Traits>
{
  friend class Instance;                            ///< Don't force a lot of pointless getters
  friend class Model;                               ///< Don't force a lot of pointless getters

public:

  //---------------------------------------------------------------------------
  // Function      : Xyce::Device::Resistor::Master::Master
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter
  // Creation Date : 11/26/08
  //---------------------------------------------------------------------------
  ///
  /// Construct a Resistor Device.
  ///
  /// @param configuration
  /// @param factory_block
  /// @param solver_state
  /// @param device_options
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options),
      separateInstances_(false)
  {}

  virtual bool updateState(double * solVec, double * staVec, double * stoVec)
  {
    return true;
  }

  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec)
  {
    return true;
  }

  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV)
  {
    return loadDAEVectors( solVec, fVec, qVec, bVec,
                           leadF, leadQ, junctionV, ALL );
  }

  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV, int loadType);

  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
  {
    return loadDAEMatrices( dFdx, dQdx, ALL );
  }

  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx,
                               int loadType);

  virtual bool loadFreqDAEVectors(double frequency,
                                  std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec);

  virtual bool loadFreqDAEMatrices(double frequency,
                                   std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx);

  private:
    InstanceVector      linearInstances_;            ///< List of owned linear resistor instances
    InstanceVector      nonlinearInstances_;         ///< List of owned nonlinear resistor instances
    bool separateInstances_;
};

void registerDevice(const DeviceCountMap& deviceMap = DeviceCountMap(),
                    const std::set<int>& levelSet = std::set<int>());

} // namespace Resistor
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Resistor_h

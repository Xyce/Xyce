//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_DEV_Capacitor_h
#define Xyce_N_DEV_Capacitor_h

#include <N_DEV_fwd.h>
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Capacitor {

class Model;
class Instance;

/// sensitivity functor
class capSensitivity :  public baseSensitivity
{
  public:
  capSensitivity() : 
    baseSensitivity() {};

  virtual ~capSensitivity() {};

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

static capSensitivity capSens;

/// sensitivity functor
class capMatrixSensitivity :  public baseMatrixSensitivity
{
  public:
  capMatrixSensitivity() : 
    baseMatrixSensitivity() {};

  virtual ~capMatrixSensitivity() {};

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

static capMatrixSensitivity capMatrixSens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Capacitor";}
  static const char *deviceTypeName() {return "C level 1";}

  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "C";}
  static const char *instanceDefaultParameter() {return "C";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Capacitor::Instance
// Special Notes : A capacitor  will have two circuit nodes.
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
/// Capacitor instance
///
/// This class refers to a single instance of the capacitor device.  It
/// contains indicies into the matrix equation.  See the comments for the
/// Resistor::Instance class for more details.
///
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class capSensitivity;
  friend class capMatrixSensitivity;

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
  bool isLinearDevice() const;
  void registerLIDs( const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef );
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars () { return true; };
  bool updatePrimaryState ();

  bool setIC ();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  Model &       model_;         //< Owning model

  // Stuff for handling solution-variable-dependent capacitance
  Util::Expression * expPtr;
  int                expNumVars;

  std::vector<double> expVarDerivs;

  // user-specified parameters:
  double C;    // User specified capacitance. (Farads)
  double Q;    // User specified charge. (Coulomb)
  double multiplicityFactor;    // multiplicity factor (M)
  double IC;   // Optional initial value capacitor voltage (V).

  // These are for the semiconductor capacitor
  double length;    // capacitor length
  double width;     // capacitor width
  double temp;      // temperature of this instance

  // Genie 121412. temperature dependence parameters
  // these can override values specified in the model
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.

  // flags used to tell if the user has specified one of these values
  // on the command line.
  bool tempCoeff1Given;
  bool tempCoeff2Given;

  // These are for the age-aware capacitor
  double age;                 ///< age in hours
  double ageCoef;             ///< degradation coeficient.
  double baseCap;             ///< the baseline capacitance before aging

  bool tempGiven;
  bool ICGiven;
  bool solVarDepC;
  bool solVarDepQ;
  bool UIC_;                  ///< Set only on first iteration of first time step of NOOP or UIC transients

  // state variables:
  double q0;                  ///< charge in the capacitor
  double q0_Jdxp;             ///< Jdxp term for Q (Only used when doing NOOP transients

  // now held in the store vector at li_store_dev_i
  double vcap; // voltage drop across capacitor

  //local id's (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;                 ///< for the "voltage source" when IC is specified
  int li_branch_data;   ///< Index for lead current and junction voltage (for power calculations)

  int li_QState;

  std::vector<int> li_dQdXState;
  std::vector<int> li_dCdXState;
  int li_vcapState;
  int li_capState;

  // Offsets for Jacobian
  int APosEquPosNodeOffset;
  int ANegEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquNegNodeOffset;

  // offsets for when C is solution-variable dependent
  std::vector<int> APosEquDepVarOffsets;
  std::vector<int> ANegEquDepVarOffsets;

  int ABraEquPosNodeOffset;
  int ABraEquNegNodeOffset;
  int ABraEquBraNodeOffset;
  int APosEquBraNodeOffset;
  int ANegEquBraNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointers for Jacobian
  double * qPosEquPosNodePtr;
  double * qNegEquPosNodePtr;
  double * qPosEquNegNodePtr;
  double * qNegEquNegNodePtr;

  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fBraEquBraNodePtr;
  double * fPosEquBraNodePtr;
  double * fNegEquBraNodePtr;

  std::vector<double *> qPosEquDepVarsPtrs;
  std::vector<double *> qNegEquDepVarsPtrs;
#endif

  std::vector< std::vector<int> > jacStamp;
  std::vector< std::vector<int> > jacStamp_IC;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
/// Capacitor Model class
///
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;
  friend class capSensitivity;
  friend class capMatrixSensitivity;

public:
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
  bool processParams ();
  bool processInstanceParams ();
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

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

  // for the semiconductor capacitor (cap with a model card)

  double capacitanceMultiplier; // capacitance multiplier
  double cj;     // junction bottom capacitance
  double cjsw;   // junction sidewall capacitance
  double defWidth; // default width
  double narrow;   // narrowing due to side etching
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.
  double baseCap;
  double tnom;

  bool tnomGiven;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Capacitor Master class
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
/// For efficiency, the Capacitor class reimplements these methods to do the
/// work directly, instead of calling instance-level functions.
///
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       ss1,
     const DeviceOptions &     do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1),
      separateInstances_(false)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec)
  { return updateState( solVec, staVec, stoVec, ALL ); }
  virtual bool updateState (double * solVec, double * staVec, double * stoVec, int loadType);
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, 
                               double * leadF, double * leadQ, double * junctionV)
  { return loadDAEVectors( solVec, fVec, qVec, bVec, leadF, leadQ, junctionV, ALL ); }
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, 
                               double * leadF, double * leadQ, double * junctionV, int loadType);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
  { return loadDAEMatrices( dFdx, dQdx, ALL ); }
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType);

  private:
  bool separateInstances_;
  InstanceVector      linearInstances_;            ///< List of owned linear capacitor instances
  InstanceVector      nonlinearInstances_;         ///< List of owned nonlinear capacitor instances
};

void registerDevice(const DeviceCountMap& deviceMap = DeviceCountMap(),
                    const std::set<int>& levelSet = std::set<int>());

} // namespace Capacitor
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Capacitor_h


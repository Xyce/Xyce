//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Inductor classes.
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

#ifndef Xyce_N_DEV_Inductor_h
#define Xyce_N_DEV_Inductor_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Inductor {

class Model;
class Instance;

/// sensitivity functor
class indSensitivity :  public baseSensitivity
{
  public:
  indSensitivity() : 
    baseSensitivity() {};

  virtual ~indSensitivity() {};

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

static indSensitivity indSens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Inductor";}
  static const char *deviceTypeName() {return "L level 1";}
  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "L";}
  static const char *instanceDefaultParameter() {return "L";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 inductor device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the InductorInstance
//                 class for a more detailed explanation.
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
  friend class indSensitivity;

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

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars () { return true; }
  bool updatePrimaryState () { return true; }

  bool setIC ();

  void varTypes( std::vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

public:
  // iterator reference to the inductor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp_BASE;


  Model &       model_;         //< Owning model

  // parameter variables
  double L;  // User specified inductance.
  double multiplicityFactor;     // multiplicity factor (M)
  double IC; // Initial condition: initial, time-zero inductor current(A)
  bool ICGiven;
  double baseL;     // the baseline inductance before temperature effects
  double temp;      // temperature of this instance
  bool tempGiven;
  double dtemp;      // delta temperature of this instance
  bool dtempGiven;   // delta temperature given parameter

  // Genie 121412. temperature dependence parameters
  // these can override values specified in the model
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.

  // flags used to tell if the user has specified one of these values
  // on the command line.
  bool tempCoeff1Given;
  bool tempCoeff2Given;

  // state variables
  double f0; // most recent value for the  flux through the inductor.

  // local state indices (offsets)
  int li_fstate;

  // local solution indices (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;

  int li_branch_data;   ///< Index for lead current and junction voltage (for power calculations)

  // Matrix equation index variables:
  std::vector<int> xLBraVar_J;
  std::vector<int> li_LBra;

  int ABraEquLBraVar_I; // Row index for the branch current
  // contribution of inductors this instance
  // is coupled to.

  // Offset variables for all of the above index variables.
  int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
  // branch current equ.

  int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
  // branch current equ.

  int ABraEquBraVarOffset;  // Offset, branch current variable
  // contribution, branch current equation.

  int APosEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the pos node

  int ANegEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the neg node

  int AEPosEquEBraVarOffset;

  int AENegEquEBraVarOffset;

  int AEBraEquEPosNodeOffset;

  int AEBraEquENegNodeOffset;

  int AEBraEquLNegNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Pointer variables for the Jacobian matrices
  double * fPosEquBraVarPtr;
  double * fNegEquBraVarPtr;
  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fBraEquBraVarPtr;
  double * qBraEquBraVarPtr;
#endif

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
  friend class indSensitivity;

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
  virtual bool processParams ();
  virtual bool processInstanceParams ();


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

  double inductanceMultiplier;  // User specified inductance multiplier.
  double IC; // Initial condition: initial, time-zero inductor current(A)
  double tempCoeff1;     // first order temperature coeff.
  double tempCoeff2;     // second order temperature coeff.
  double baseL;
  double tnom;
  bool tc1Given;
  bool tc2Given;
  bool tnomGiven;
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
public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1),
      separateInstances_(false)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec) { return true; }
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

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
  InstanceVector      linearInstances_;            ///< List of owned linear inductor instances
  InstanceVector      nonlinearInstances_;         ///< List of owned nonlinear inductor instances
};

void registerDevice();

} // namespace Inductor
} // namespace Device
} // namespace Xyce

#endif

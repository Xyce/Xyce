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
// Purpose        : Non-Linear Mutual Inductor classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/21/05
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MutIndNonLin_h
#define Xyce_N_DEV_MutIndNonLin_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

// defines simple container InductorInstanceData
#include <N_DEV_MutIndLin.h>

#include <Teuchos_RCP.hpp>
#include <fstream>

namespace Xyce {
namespace Device {
namespace MutIndNonLin {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Nonlinear Mutual Inductor";}
  static const char *deviceTypeName() {return "K level 1";}
  static int numNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the nonlinear
//                 mutual inductor device.
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
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
     const InstanceBlock &       IB,
     Model &                     Iiter,
     const FactoryBlock &        factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & staLIDVecRef );

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);
  void updateInductanceMatrix();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool setIC ();

  bool plotfileFlag () {return true;}

  void varTypes( std::vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // DAE load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void auxDAECalculations ();

  bool outputPlotFiles(bool force_final_output);

  // iterator reference to the inductor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

  // function is used in deviceMgr to help determine if the instance's
  // loadLeadCurrent variable needs to be set for a given mutual inductor instance.
  std::vector< std::string > getInductorNames() const
  {
    return inductorNames;
  }

  std::vector< double > getInductorInductances() const  {return inductorInductances;}

  void setInductorInductances(std::vector< double > & set) 
  {
    if (set.size() == inductorInductances.size()) { for (int ii=0;ii<set.size();ii++) { inductorInductances[ii] = set[ii]; } }
  }

private:

  Model &       model_;         //< Owning model

  // This container bundles up the physical data for each inductor
  // involved in this mutual inductor
  int numInductors;
  double L;
  std::vector< InductorInstanceData* > instanceData;

  // These vectors let the new param options load and set inductor data
  // the parser passes all of these to us
  std::vector< std::string > inductorNames;
  std::vector< double > inductorInductances;
  std::vector< std::string > inductorsNode1;
  std::vector< std::string > inductorsNode2;
  // and here's the list of ones we are coupling
  std::vector< std::string > couplingInductor;
  std::vector< double > couplingCoefficient;
  std::vector< double > initialCondition;
  std::vector< bool > initialConditionGiven;
  //std::vector< std::vector< double > > mutualCouplingCoef;

  // local indices for extra equations
  int li_MagVar;
  int li_RVar;

  // offsets in the jacobian
  int mEquVPosOffset, mEquVNegOffset;
  std::vector< int > mEquInductorOffsets;
  int mEquMOffset, mEquROffset;

  int rEquROffset;
  std::vector< int > rEquInductorOffsets;

  // state variable for mag, h and r
  int li_MagVarState;
  int li_MagVarDerivState;
  int li_RVarStore;
  int li_BVarStore;
  int li_HVarStore;

  double nonlinFlag;   // flag created by parser.  Don't need it but must read it in
  bool nonlinFlagGiven;
  double mutualCup;    // mutaul coupling value
  bool mutualCupGiven;

  std::vector< double > inductanceVals;      // the inductances of the inductors
  std::vector< std::vector< double > > LO;        // L' * L (matrix)
  std::vector< double > inductorCurrents;    // currents through inductors (col vec.)
  std::vector< double > LOI;                 // LO * I (col vector).

  double temp;         // temperature of this instance
  bool tempGiven;      // flag if temp was given

  // intermediate values needed for nonlinear model
  double qV;
  double delM0;        // modeling constant
  double Happ;
  double He;
  double Heo;
  double delM;
  double Mirrp;
  double Manp;
  double P;
  double dP_dM;
  double dP_dVp;
  double dP_dVn;

  // variables for limiting of non-linear, internal vars
  double Mag_orig;
  double Rvar_orig;

  // these vectors are used repeadily in loadDAEdFdx
  // so rather than create and destroy them on each call
  // we will allocate them in the constructor.
  std::vector< double > dHe_dI;
  std::vector< double > dManp_dI;
  std::vector< double > ddelM_dI;
  std::vector< double > dMirrp_dI;
  std::vector< double > dP_dI;

  // scaling crontrol for new equations
  double scalingRHS;   // scaling for loading DAE or RHS components

  // used in scaling the tanh() approximation to sgn()
  double maxVoltageDrop;

  std::vector< std::vector<int> > jacStamp;

  // output stream for output of internal state if requested by user
  Teuchos::RCP< std::ofstream > outputFileStreamPtr;
  bool outputStateVarsFlag;

};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
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

  // Data Members for Associations

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

  double A;                  // Thermal energy parameter (amp/m)
  double Alpha;              // domain coupling parameter (dimensionless)
  double AreaInCm2;          // mean magnetic cross-sectional area (cm^2)
  double Area;               // mean magnetic cross-sectional area (m^2)
  double BetaH;              // modeling constant (dimensionless)
  double BetaM;              // modeling constant (dimensionless)
  double C;                  // domain flesing parameter (dimensionless)
  double CLim;               // Value below which domain flesing parameter (dimensionless) will be treated as zero.
  double DeltaVScaling;      // smoothing coefficient for V_1 in tanh
  double GapInCm;            // effective air gap (cm)
  double Gap;                // effective air gap (m)
  double Kirr;               // domain anisotropy parameter (amp/m)
  double Ms;                 // saturation magnetization (amp/m)
  double LevelIgnored;       // for pspice compatibility -- ignored
  double PackIgnored;        // for pspice compatibility -- ignored
  double PathInCm;           // total mean magnetic path (cm)
  double Path;               // total mean magnetic path (m)
  double tempCoeff1;         // first order temperature coeff.
  double tempCoeff2;         // second order temperature coeff.
  double tnom;               // reference temperature
  double pZeroTol;           // absolute value below which to consider P=0
  double HCgsFactor;         // conversion factor to change H from SI units to CGS units (H/M to Oersted)
  double BCgsFactor;         // conversion factor to change B form SI units (Tesla) to CGS units (Gauss)
  double mVarScaling;        // scaling for M variable
  double rVarScaling;        // scaling for R variable
  double mEqScaling;         // scaling for M equation
  double rEqScaling;         // scaling for r equation
  bool mVarScalingGiven;     // flag to indicate that mVarScaling was given in the model block
  bool rVarScalingGiven;     // flag to indicate that rVarScaling was given in the model block
  bool mEqScalingGiven;      // flag to indicate that mEqScaling was given in the model block
  bool rEqScalingGiven;      // flag to indicate that rEqScaling was given in the model block
  int outputStateVars;       // flag indicating if user wants M,H and R output
  int factorMS;              // flag to factor Ms out of M
  bool factorMSGiven;
  int BHSiUnits;             // flag to indicate that B and H should be output in SI units. Default is CGS
                             // units for output while SI units are used for calculations.

  // flags indicating if temperature parameters were given
  bool tc1Given;
  bool tc2Given;
  bool tnomGiven;
  bool UseConstantDeltaVScaling;
  bool includeMEquation;
  bool includeMEquationGiven;
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

void registerDevice();

} // namespace MutIndNonLin
} // namespace Device
} // namespace Xyce

#endif

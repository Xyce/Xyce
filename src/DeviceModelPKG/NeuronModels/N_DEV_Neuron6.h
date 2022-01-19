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
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 06/10/09
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron6_h
#define Xyce_N_DEV_Neuron6_h

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_MembraneModel.h>

#include <N_DEV_Neuron.h>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;

namespace Xyce {
namespace Device {
namespace Neuron6 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Neuron::Traits>
{
  static const char *name() {return "Neuron";}
  static const char *deviceTypeName() {return "YNEURON level 6";}
  static int numNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Neuron device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the NeuronInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
    
public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       IB,
     Model &                     Miter,
     const FactoryBlock &        factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool setIC ();

  void varTypes( std::vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

public:
  // iterator reference to the Neuron model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // model level parameters that can be overridden at the instance level
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  double segArea;  // segment area (derrived from radius, length and nSeg)
  int    nSeg;     // number of segments
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool nSegGiven;

  // conductance between segments -- calculated from radius, rInt and length and number of segments
  double gSeg;

  int numIntVarsPerSegment;
  int numStateVarsPerSegment;

  // storage for local ID's of internal vars and jacobian offsets
  std::vector< int > li_internalVars;
  std::vector< std::vector< int > > jacobianOffsets;

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions (no q terms on the external nodes)
  double kcl1Fvalue;
  double kcl2Fvalue;
  // internal segments
  std::vector<double> segFvalue;
  std::vector<double> segQvalue;
  std::vector<double> segNEquFvalue, segNEquQvalue;
  std::vector<double> segMEquFvalue, segMEquQvalue;
  std::vector<double> segHEquFvalue, segHEquQvalue;
  std::vector<double> segAEquFvalue, segAEquQvalue;
  std::vector<double> segBEquFvalue, segBEquQvalue;
  std::vector<double> segM_EquFvalue, segM_EquQvalue;
  std::vector<double> segH_EquFvalue, segH_EquQvalue;
  std::vector<double> segCEquFvalue, segCEquQvalue;
  std::vector<double> segCaEquFvalue, segCaEquQvalue;

  // jacobian terms
  double dkcl1F_dVin, dkcl1F_dVs0;
  double dkcl2F_dVout, dkcl2F_dVsn;
  // internal equations
  std::vector<double> segF_dVp, segF_dV, segF_dVn, segF_dn, segF_dm, segF_dh, segF_da, segF_db, segF_dM, segF_dH, segF_dc;
  std::vector<double> segQ_dV;
  std::vector<double> dnF_dV, dnF_dn, dnQ_dn;
  std::vector<double> dmF_dV, dmF_dm, dmQ_dm;
  std::vector<double> dhF_dV, dhF_dh, dhQ_dh;
  std::vector<double> daF_dV, daF_da, daQ_da;
  std::vector<double> dbF_dV, dbF_db, dbQ_db;
  std::vector<double> dMF_dV, dMF_dM, dMQ_dM;
  std::vector<double> dHF_dV, dHF_dH, dHQ_dH;
  std::vector<double> dcF_dV, dcF_dc, dcF_dCa, dcQ_dc;
  std::vector<double> dCaF_dV, dCaF_dM, dCaF_dH, dCaF_dCa, dCaQ_dCa;

  // state variables
  std::vector<double> potassiumCurrent;
  std::vector<double> sodiumCurrent;

  // local state indices (offsets)
  std::vector<int> li_KCurrentState;
  std::vector<int> li_NaCurrentState;

  // local solution indices (offsets)
  int li_Pos;      // local index to positive node on this device
  int li_Neg;      // local index to negative node on this device
  // local solution indices for internal vars (variable number of these)
  std::vector<int> li_Vol;      // local index to segment voltage
  std::vector<int> li_nPro;     // local index to n promoter value (Na current)
  std::vector<int> li_mPro;     // local index to m promoter value (K current)
  std::vector<int> li_hPro;     // local index to h promoter value (K current)
  std::vector<int> li_aPro;     // local index to a promoter value
  std::vector<int> li_bPro;     // local index to a promoter value
  std::vector<int> li_MPro;     // local index to a promoter value
  std::vector<int> li_HPro;     // local index to a promoter value
  std::vector<int> li_cPro;     // local index to a promoter value
  std::vector<int> li_CaPro;     // local index to a promoter value

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset, APosEquNextNodeOffset;
  int ANegEquNegNodeOffset, ANegEquLastNodeOffset;
  std::vector<int> SegVEqnVpreOffset;
  std::vector<int> SegVEqnVsegOffset;
  std::vector<int> SegVEqnVnexOffset;
  std::vector<int> SegVEqnNOffset;
  std::vector<int> SegVEqnMOffset;
  std::vector<int> SegVEqnHOffset;
  std::vector<int> SegVEqnAOffset;
  std::vector<int> SegVEqnBOffset;
  std::vector<int> SegVEqnM_Offset;
  std::vector<int> SegVEqnH_Offset;
  std::vector<int> SegVEqnCOffset;
  std::vector<int> NEquVNodeOffset;
  std::vector<int> NEquNNodeOffset;
  std::vector<int> MEquVNodeOffset;
  std::vector<int> MEquMNodeOffset;
  std::vector<int> HEquVNodeOffset;
  std::vector<int> HEquHNodeOffset;
  std::vector<int> AEquVNodeOffset;
  std::vector<int> AEquANodeOffset;
  std::vector<int> BEquVNodeOffset;
  std::vector<int> BEquBNodeOffset;
  std::vector<int> M_EquVNodeOffset;
  std::vector<int> M_EquM_NodeOffset;
  std::vector<int> H_EquVNodeOffset;
  std::vector<int> H_EquH_NodeOffset;
  std::vector<int> CEquVNodeOffset;
  std::vector<int> CEquCNodeOffset;
  std::vector<int> CEquCaNodeOffset;
  std::vector<int> CaEquVNodeOffset;
  std::vector<int> CaEquM_NodeOffset;
  std::vector<int> CaEquH_NodeOffset;
  std::vector<int> CaEquCaNodeOffset;

  // maps to track the appropriate jacobian offsets for each segment's previous, current, and next segment
  std::map <int, int> prevMap;
  std::map<int, int> segMap;
  std::map<int, int> nextMap;

  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  friend class ParametricData<Model>;
  typedef std::vector<Instance *> InstanceVector;


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

private:

  // parameter variables
  double cMem;     // membrane capacitance
  double gMem;     // membrane conductance
  double vRest;    // resting potential
  double eNa;      // sodium rest potential
  double gNa;      // sodium base conductance
  double eK;       // potassium rest potential
  double gK;       // potassium base conductance
  double eA;       // a-current rest potential
  double gA;       // a-current base conductance
  double eCa;      // Calcium rest potential
  double gCa;      // Calcium base conductance
  double eKCa;     // potassium-calcium rest potential
  double gKCa;     // potassium-calcium base conductance
  double CaInit;  // initial intra-cellular calcium concentration
  double CaGamma;  // calcium current to concentration multiplier
  double CaTau;    // calcium removal time constant
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  std::string ionChannelModel; // what model will be used for the ion channels
  int    nSeg;     // number of segments

  // Value of current expression for user-defined membranem model
  double I;

  // these are vectors of strings to allow the user to specify independant vars and
  // equations for a given membrane
  std::vector<std::string> membraneCurrentEqus;
  std::vector<std::string> membraneIndpVars;
  std::vector<std::string> membraneIndpFEqus;
  std::vector<std::string> membraneIndpQEqus;
  std::vector<std::string> membraneFunctions;
  std::vector<std::string> membraneParameters;

  // flags that parameters were given
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool ionChannelModelGiven;
  bool nSegGiven;

  // flags that parameters were given
  bool cMemGiven;
  bool gMemGiven;
  bool vRestGiven;
  bool eNaGiven;
  bool gNaGiven;
  bool eKGiven;
  bool gKGiven;
  bool eAGiven;
  bool gAGiven;
  bool eCaGiven;
  bool gCaGiven;
  bool eKCaGiven;
  bool gKCaGiven;
  bool CaInitGiven;
  bool CaGammaGiven;
  bool CaTauGiven;

  bool membraneCurrentEqusGiven;
  bool membraneIndpVarsGiven;
  bool membraneIndpFEqusGiven;
  bool membraneIndpQEqusGiven;
  bool membraneFunctionsGiven;
  bool membraneParametersGiven;

  // these are meta flags.  If one component of the a sodium current is given
  // then all of the required equaitons will be used.  Otherwise they are off
  // Or, if hodgenHuxleyOn_ is true, then all of the H odgenHuxley equations are loaded.
  // by default all of these are off.
  bool hodgenHuxleyOn_;
  bool ConnorStevensOn_;
  bool sodiumOn_;
  bool potassiumOn_;
  bool aCurrentOn_;
  bool calciumOn_;
  RCP< MembraneModel > membraneModel_;


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
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Neuron6
} // namespace Device
} // namespace Xyce

#endif

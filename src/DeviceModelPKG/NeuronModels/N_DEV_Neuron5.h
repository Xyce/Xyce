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
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/02/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron5_h
#define Xyce_N_DEV_Neuron5_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Neuron.h>

namespace Xyce {
namespace Device {
namespace Neuron5 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Neuron::Traits>
{
  static const char *name() {return "Neuron";}
  static const char *deviceTypeName() {return "YNEURON level 5";}
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
// Creation Date : 01/02/08
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

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions
  double kcl1Fvalue, kcl1Qvalue;
  double kcl2Fvalue, kcl2Qvalue;
  double nEquFvalue, nEquQvalue;
  double mEquFvalue, mEquQvalue;
  double hEquFvalue, hEquQvalue;
  double aEquFvalue, aEquQvalue;
  double bEquFvalue, bEquQvalue;
  double M_EquFvalue, M_EquQvalue;
  double H_EquFvalue, H_EquQvalue;
  double cEquFvalue, cEquQvalue;
  double CaEquFvalue, CaEquQvalue;

  double dkcl1F_dV1, dkcl1F_dV2, dkcl1F_dn, dkcl1F_dm, dkcl1F_dh,
    dkcl1F_da, dkcl1F_db, dkcl1F_dM, dkcl1F_dH, dkcl1F_dc, dkcl1Q_dV1, dkcl1Q_dV2;
  double dkcl2F_dV1, dkcl2F_dV2, dkcl2F_dn, dkcl2F_dm, dkcl2F_dh,
    dkcl2F_da, dkcl2F_db, dkcl2F_dM, dkcl2F_dH, dkcl2F_dc, dkcl2Q_dV1, dkcl2Q_dV2;
  double dnF_dV1, dnF_dn, dnQ_dn;
  double dmF_dV1, dmF_dm, dmQ_dm;
  double dhF_dV1, dhF_dh, dhQ_dh;
  double daF_dV1, daF_da, daQ_da;
  double dbF_dV1, dbF_db, dbQ_db;
  double dMF_dV1, dMF_dM, dMQ_dM;
  double dHF_dV1, dHF_dH, dHQ_dH;
  double dcF_dV1, dcF_dc, dcF_dCa, dcQ_dc;
  double dCaF_dV1, dCaF_dV2, dCaF_dM, dCaF_dH, dCaF_dCa, dCaQ_dCa;

  // state variables
  double potassiumCurrent;
  double sodiumCurrent;

  // local state indices (offsets)
  int li_KCurrentState;
  int li_NaCurrentState;

  // local solution indices (offsets)
  int li_Pos;      // local index to positive node on this device
  int li_Neg;      // local index to negative node on this device
  int li_nPro;     // local index to n promoter value (Na current)
  int li_mPro;     // local index to m promoter value (K current)
  int li_hPro;     // local index to h promoter value (K current)
  int li_aPro;     // local index to
  int li_bPro;     // local index
  int li_M_Pro;    // local index
  int li_H_Pro;    // local index
  int li_cPro;     // local index
  int li_CaPro;    // local index

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int APosEquNNodeOffset;
  int APosEquMNodeOffset;
  int APosEquHNodeOffset;
  int APosEquANodeOffset;
  int APosEquBNodeOffset;
  int APosEquM_NodeOffset;
  int APosEquH_NodeOffset;
  int APosEquCNodeOffset;

  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;
  int ANegEquNNodeOffset;
  int ANegEquMNodeOffset;
  int ANegEquHNodeOffset;
  int ANegEquANodeOffset;
  int ANegEquBNodeOffset;
  int ANegEquM_NodeOffset;
  int ANegEquH_NodeOffset;
  int ANegEquCNodeOffset;

  int ANEquPosNodeOffset;
  int ANEquNNodeOffset;

  int AMEquPosNodeOffset;
  int AMEquMNodeOffset;

  int AHEquPosNodeOffset;
  int AHEquHNodeOffset;

  int AAEquPosNodeOffset;
  int AAEquANodeOffset;

  int ABEquPosNodeOffset;
  int ABEquBNodeOffset;

  int AM_EquPosNodeOffset;
  int AM_EquM_NodeOffset;

  int AH_EquPosNodeOffset;
  int AH_EquH_NodeOffset;

  int ACEquPosNodeOffset;
  int ACEquCNodeOffset;
  int ACEquCaNodeOffset;

  int ACaEquPosNodeOffset;
  int ACaEquNegNodeOffset;
  int ACaEquM_NodeOffset;
  int ACaEquH_NodeOffset;
  int ACaEquCaNodeOffset;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/02/08
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

} // namespace Neuron5
} // namespace Device
} // namespace Xyce

#endif

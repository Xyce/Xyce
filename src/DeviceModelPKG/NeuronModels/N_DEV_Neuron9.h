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
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Christy Warrender, SNL, Cognitive Modeling
//
// Creation Date  : 06/22/12
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron9_h
#define Xyce_N_DEV_Neuron9_h

#include <N_DEV_Configuration.h>
#include <Sacado_No_Kokkos.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Neuron.h>

namespace Xyce {
namespace Device {
namespace Neuron9 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Neuron::Traits>
{
  static const char *name() {return "Neuron";}
  static const char *deviceTypeName() {return "YNEURON level 9";}
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
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;friend class Master;

public:
  static std::vector< std::vector<int> > jacStamp;

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Miter,
     const FactoryBlock &      factory_block);


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

private:
  // These functions represent the equations that need to be solved
  // for this device.  Since Xyce loads an F and Q contribution, the
  // equations are broken up into their F and Q components.  Thus there
  // is a kcl1EquF() and kcl1EquQ().  Automatic differentiation will
  // be used to generate all the derivatives of these equations for the
  // dF/dX and dQ/dX loads

  // first we list some utility functions for calculating coefficients.
  // alpha and beta equations are taken from Brette et al 07.
  // They're generally functions of membrane voltage; here, the membrane voltage is
  // the difference between Vn1 and Vn2.
  // These functions expect V to be in milli-volts and then return values that
  // are in 1/ms.  Thus the extra factor's of 1000 here and there

  // potassium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaN( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2);  // convert voltage to milli-volts
    ScalarT r;
    r = 0.032*(15.0 - vDiff + VT)/( std::exp( (15.0 - vDiff + VT)/5.0) - 1.0);
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaN( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2);  // convert voltage to milli-volts
    ScalarT r;
    r = 0.5*( std::exp( (10.0 - vDiff + VT)/40.0) );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  // sodium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaM( const ScalarT & Vn1, const ScalarT  & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2);  // convert voltage to milli-volts
    ScalarT r;
    r = 0.32*(13.0 - vDiff + VT )/( std::exp((13.0 - vDiff + VT )/4.0) - 1.0);
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaM( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2);  // convert voltage to milli-volts
    ScalarT r;
    r = 0.28*(vDiff - VT - 40.0)/( std::exp((vDiff - VT - 40.0)/5.0) - 1.0);
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT alphaH( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2);  // convert voltage to milli-volts
    ScalarT r;
    r = 0.128 * std::exp( (17.0 - vDiff + VT)/18.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaH( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2);  // convert voltage to milli-volts
    ScalarT r;
    r = 4.0/( 1.0 + std::exp( (40.0 - vDiff + VT)/5.0) );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  // now the device equations
  // KCL equation 1
  template <typename ScalarT>
  static ScalarT kcl1EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& memG, const ScalarT& leakE, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE )
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    ScalarT r = memG * (Vn1 - Vn2 - leakE) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE );
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl1EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
  {
    ScalarT r = memC * (Vn1 - Vn2);
    return r;
  }

  // KCL equation 2 -- -1 * equation 1 because of device symmetry
  template <typename ScalarT>
  static ScalarT kcl2EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& memG, const ScalarT& leakE, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE )
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    ScalarT r = -1.0*(memG * (Vn1 - Vn2 - leakE) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE ));
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl2EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
  {
    ScalarT r = -1.0 * memC * (Vn1 - Vn2);
    return r;
  }

  // n conservation equation
  template <typename ScalarT>
  static ScalarT nEquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& Vrest )
  {
    ScalarT r = alphaN<ScalarT>( Vn1, Vn2, Vrest ) * (1.0 - n ) - betaN<ScalarT>( Vn1, Vn2, Vrest ) * n;
    return r;
  }

  template <typename ScalarT>
  static ScalarT nEquQ( const ScalarT& n )
  {
    ScalarT r = -n;
    return r;
  }

  // m conservation equation
  template <typename ScalarT>
  static ScalarT mEquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& m, const ScalarT& Vrest )
  {
    ScalarT r = alphaM<ScalarT>( Vn1, Vn2, Vrest ) * (1.0 - m ) - betaM<ScalarT>( Vn1, Vn2, Vrest ) * m;
    return r;
  }

  template <typename ScalarT>
  static ScalarT mEquQ( const ScalarT& m )
  {
    ScalarT r = -m;
    return r;
  }

  // h conservation equation
  template <typename ScalarT>
  static ScalarT hEquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& h, const ScalarT& Vrest )
  {
    ScalarT r = alphaH<ScalarT>( Vn1, Vn2, Vrest ) * (1.0 - h ) - betaH<ScalarT>( Vn1, Vn2, Vrest ) * h;
    return r;
  }

  template <typename ScalarT>
  static ScalarT hEquQ( const ScalarT& h )
  {
    ScalarT r = -h;
    return r;
  }

public:
  // iterator reference to the Neuron model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // constant threshold adjustment
  static const double VT;

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions
  double kcl1Fvalue, kcl1Qvalue;
  double kcl2Fvalue, kcl2Qvalue;
  double nEquFvalue, nEquQvalue;
  double mEquFvalue, mEquQvalue;
  double hEquFvalue, hEquQvalue;
  double dkcl1F_dV1, dkcl1F_dV2, dkcl1F_dn, dkcl1F_dm, dkcl1F_dh, dkcl1Q_dV1, dkcl1Q_dV2;
  double dkcl2F_dV1, dkcl2F_dV2, dkcl2F_dn, dkcl2F_dm, dkcl2F_dh, dkcl2Q_dV1, dkcl2Q_dV2;
  double dnF_dV1, dnF_dV2, dnF_dn, dnQ_dn;
  double dmF_dV1, dmF_dV2, dmF_dm, dmQ_dm;
  double dhF_dV1, dhF_dV2, dhF_dh, dhQ_dh;

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

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int APosEquNNodeOffset;
  int APosEquMNodeOffset;
  int APosEquHNodeOffset;

  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;
  int ANegEquNNodeOffset;
  int ANegEquMNodeOffset;
  int ANegEquHNodeOffset;

  int ANEquPosNodeOffset;
  int ANEquNegNodeOffset;
  int ANEquNNodeOffset;

  int AMEquPosNodeOffset;
  int AMEquNegNodeOffset;
  int AMEquMNodeOffset;

  int AHEquPosNodeOffset;
  int AHEquNegNodeOffset;
  int AHEquHNodeOffset;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/22/12
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

  bool processParams ();
  bool processInstanceParams ();

private:

  // parameter variables
  double cMem;     // membrane capacitance
  double gMem;     // membrane conductance of leak current
  double eLeak;    // reversal potential of leak current
  double eNa;      // sodium reversal potential
  double gNa;      // sodium base conductance
  double eK;       // potassium reversal potential
  double gK;       // potassium base conductance
  double vRest;    // resting potential

  // flags that parameters were given
  bool cMemGiven;
  bool gMemGiven;
  bool eLeakGiven;
  bool eNaGiven;
  bool gNaGiven;
  bool eKGiven;
  bool gKGiven;
  bool vRestGiven;


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


//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 06/01/12
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Neuron9
} // namespace Device
} // namespace Xyce

#endif

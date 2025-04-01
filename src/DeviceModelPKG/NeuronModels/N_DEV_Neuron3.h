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

#ifndef Xyce_N_DEV_Neuron3_h
#define Xyce_N_DEV_Neuron3_h

#include <N_DEV_Configuration.h>
#include <Sacado_No_Kokkos.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Neuron.h>

namespace Xyce {
namespace Device {
namespace Neuron3 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Neuron::Traits>
{
  static const char *name() {return "Neuron";}
  static const char *deviceTypeName() {return "YNEURON level 3";}
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

protected:
private:
  // These functions represents the equations that need to be solved
  // for this device.  Since Xyce loads an F and Q contribution, the
  // equations are broken up into their F and Q components.  Thus there
  // is a kcl1EquF() and kcl1EquQ().  Automatic differentiation will
  // be used to generate all the derivatives of these equations for the
  // dF/dX and dQ/dX loads

  // first we list some utility functions for calculating coefficients
  //
  // These functions expect V to be in milli-volts and then return values that
  // are in 1/ms.  Thus the extra factor's of 1000 here and there

  // potassium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaN( const ScalarT & Vn1, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vrest);  // convert voltage to milli-volts
    ScalarT r;
    //       if ((vDiff > 9.99) && (vDiff < 10.01) )
    //       {
    //         r = 1.0/(10.0 * ( std::exp( (10.0 - vDiff)/10.0 )));
    //       }
    //       else
    //       {
    r = (10.0 - vDiff) /
      (100.0 * ( std::exp( (10.0 - vDiff)/10.0 ) - 1.0 ));
    //      }
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaN( const ScalarT & Vn1, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vrest);
    ScalarT r = 0.125 * std::exp( -vDiff/80.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  // sodium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaM( const ScalarT & Vn1, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vrest);
    ScalarT r;
    //       if ((vDiff > 24.99) && (vDiff < 25.01) )
    //       {
    //         r = (1.0) /
    //                   (( std::exp( (25.0 - vDiff)/10.0 )));
    //       }
    //       else
    //       {
    r = (25.0 - vDiff) /
      (10.0 * ( std::exp( (25.0 - vDiff)/10.0 ) - 1.0 ));
    //      }
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaM( const ScalarT & Vn1, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vrest);
    ScalarT r = 4.0 * std::exp( -vDiff/18.0 );
    r *= 1000.0; // change from 1/ms to 1/s

    return r;
  }

  template <typename ScalarT>
  static ScalarT alphaH( const ScalarT & Vn1, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vrest);
    ScalarT r = 0.07 * std::exp( -vDiff/20.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaH( const ScalarT & Vn1, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vrest);
    ScalarT r = 1.0 / ( std::exp( (30.0 - vDiff)/10.0 ) + 1.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  // now the device equations
  // KCL equation 1
  // VSeg -- segment voltage
  // VSegP -- previous segment voltage
  // VSegN -- next segment voltage
  // n, m, h -- vars which determine current through themembrane
  // gPrev, gNext -- conductivity to previous and next segment (can be zero)
  //
  // Full equation is:
  // i_trans_membrane  - i_externally_injected - gNext (VSegN - VSeg) - gPrev (VSegP - VSeg ) + Cmem dVseg/dt = 0
  //
  // i_tarns_membrane is determined by the ion channel equations -- hodgkin-huxley in this case.
  // i_externally_injeted is not supported here but could be an externally injected current.
  //
  template <typename ScalarT>
  static ScalarT kcl1EquF( const ScalarT& VSeg, const ScalarT& VSegP, const ScalarT & VSegN, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& gPrev, const ScalarT& gNext, const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE )
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    //ScalarT r = memG * (VSeg - restV) + Kg * powN * (VSeg - Ke ) + NaG * powM * h * (VSeg - NaE ) - gNext * (VSegN - VSeg) - gPrev * (VSegP - VSeg);
    ScalarT r = memG * (VSeg - restV) + Kg * powN * (VSeg - Ke ) + NaG * powM * h * (VSeg - NaE )
      - gNext * (VSegN - VSeg) - gPrev * (VSegP - VSeg);
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl1EquQ( const ScalarT& VSeg, const ScalarT& memC )
  {
    ScalarT r = memC * VSeg;
    return r;
  }

  // n conservation equation
  template <typename ScalarT>
  static ScalarT nEquF( const ScalarT& Vn1, const ScalarT& n, const ScalarT& Vrest )
  {
    ScalarT r = alphaN<ScalarT>( Vn1, Vrest ) * (1.0 - n ) - betaN<ScalarT>( Vn1, Vrest ) * n;
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
  static ScalarT mEquF( const ScalarT& Vn1, const ScalarT& m, const ScalarT& Vrest )
  {
    ScalarT r = alphaM<ScalarT>( Vn1, Vrest ) * (1.0 - m ) - betaM<ScalarT>( Vn1, Vrest ) * m;
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
  static ScalarT hEquF( const ScalarT& Vn1, const ScalarT& h, const ScalarT& Vrest )
  {
    ScalarT r = alphaH<ScalarT>( Vn1, Vrest ) * (1.0 - h ) - betaH<ScalarT>( Vn1, Vrest ) * h;
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

  std::vector< std::vector<int> > jacStamp;

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

  // the cable equation is dependent on parameters from the previous and next segment of other
  // cables attached to the input/ouput nodes (i.e. when we branch a cable).  This lets
  // one set the parameters for the previous/next cable segment.
  double rIntPrevious;
  double radiusPrevious;
  double lengthPrevious;
  double rIntNext;
  double radiusNext;
  double lengthNext;
  bool   rIntPreviousGiven;
  bool   radiusPreviousGiven;
  bool   lengthPreviousGiven;
  bool   rIntNextGiven;
  bool   radiusNextGiven;
  bool   lengthNextGiven;

  // conductance between segments -- calculated from radius, rInt and length and number of segments
  // gBackward connects to previous node. gForward connects to the next node
  std::vector<double> gBackward;
  std::vector<double> gForward;

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
  // jacobian terms
  double dkcl1F_dVin, dkcl1F_dVs0;
  double dkcl2F_dVout, dkcl2F_dVsn;
  // internal equations
  std::vector<double> segF_dVp, segF_dV, segF_dVn, segF_dn, segF_dm, segF_dh;
  std::vector<double> segQ_dV;
  std::vector<double> dnF_dV, dnF_dn, dnQ_dn;
  std::vector<double> dmF_dV, dmF_dm, dmQ_dm;
  std::vector<double> dhF_dV, dhF_dh, dhQ_dh;

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
  std::vector<int> NEquVNodeOffset;
  std::vector<int> NEquNNodeOffset;
  std::vector<int> MEquVNodeOffset;
  std::vector<int> MEquMNodeOffset;
  std::vector<int> HEquVNodeOffset;
  std::vector<int> HEquHNodeOffset;
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
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  int    nSeg;     // number of segments

  // the cable equation is dependent on parameters from the previous and next segment of other
  // cables attached to the input/ouput nodes (i.e. when we branch a cable).  This lets
  // one set the parameters for the previous/next cable segment.
  double rIntPrevious;
  double radiusPrevious;
  double lengthPrevious;
  double rIntNext;
  double radiusNext;
  double lengthNext;
  bool   rIntPreviousGiven;
  bool   radiusPreviousGiven;
  bool   lengthPreviousGiven;
  bool   rIntNextGiven;
  bool   radiusNextGiven;
  bool   lengthNextGiven;

  // flags that parameters were given
  bool cMemGiven;
  bool gMemGiven;
  bool vRestGiven;
  bool eNaGiven;
  bool gNaGiven;
  bool eKGiven;
  bool gKGiven;
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool nSegGiven;


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

} // namespace Neuron3
} // namespace Device
} // namespace Xyce

#endif

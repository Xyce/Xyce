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
// Purpose        : Resistor with zero resistance.
//
// Special Notes  : This is a special case of the resistor device when it is
//                  specified with zero resistance.  It will then behave like
//                  a voltage source with zero voltage difference so that a
//                  branch current is calculated.
//
// Creator        : Richard Schiek, Electrical and Microsystems Modeling.
//
// Creation Date  : 02/01/10
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Resistor3_h
#define Xyce_N_DEV_Resistor3_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Resistor.h>

namespace Xyce {
namespace Device {
namespace Resistor3 {

// ---------- Forward Declarations ----------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Resistor::Traits>
{
  static const char *name() {return "Resistor";}
  static const char *deviceTypeName() {return "R level 3";}
  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "R";}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
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
     const InstanceBlock &     IB,
     Model &                   Viter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  bool isLinearDevice() const { return true; }

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateIntermediateVars () { return true; }
  bool updatePrimaryState () { return true; }

  // load functions, residual:
  bool loadDAEQVector () { return true; }
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () { return true; }
  bool loadDAEdFdx ();

  void setupPointers ();

  void varTypes( std::vector<char> & varTypeVec );

  void getLIDs(int & lpos, int & lneg,int & lbra)
  {lpos = li_Pos; lneg = li_Neg; lbra = li_Bra;}

public:
  // iterator reference to the vsrc model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;
  static std::vector< std::vector<int> > jacStampPDE;
  static ParametricData<Instance>       parMap_;

  Model &       model_;         //< Owning model

  // Parameters
  // user-specified paramters:
  double R;  // resistance  (ohms)
  double multiplicityFactor;     // multiplicityFactor (M)
  // these are for the semiconductor resistor
  double length;      // resistor length.
  double width;      // resistor width.
  double temp;   // temperature of this instance
  // temperature dependence parameters
  // these can override values specified in the model
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.
  double tempCoeffExp; // Exponential temperature coeff.
  double dtemp;        // externally specified device temperature.  This parameter
  // is NOT used and is only here for compatibility in parsing
  // netlist from simulators that do support it.
  // flags used to tell if the user has specified one of these values
  // on the command line.
  bool tempCoeff1Given;   // First order temperature coeff was given in netlist
  bool tempCoeff2Given;   // Second order temperature coeff was given in netlist
  bool tempCoeffExpGiven; // Exponential temperature coeff was given in netlist
  bool dtempGiven;

  // indices into state vector:
  int istate_I;  // index for i0;

  // Matrix equation index variables:

  //local indices (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;

  // Jacobian matrix indices:
  //Locally indexed offsets for jacobian
  int ABraEquPosNodeOffset; // Offset, pos. node voltage contribution,
  // branch current equ.

  int ABraEquNegNodeOffset; // Offset, neg. node voltage contribution,
  // branch current equ.

  int APosEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the pos node

  int ANegEquBraVarOffset;  // Offset, branch current variable
  // contribution, KCL equation of the neg node

  //  The following jacobian offsets are only neccessary
  // for 2-level newton.
  int APosEquPosNodeOffset;  // Offset, positive node variable
  // contribution, positive node KCL.

  int ANegEquNegNodeOffset;  // Offset, negative node variable
  // contribution, negative node KCL.

  int ABraEquBraVarOffset;  // Offset, branch current variable
  // contribution, branch current equation.

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Jacobian matrix pointers:
  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fPosEquBraVarPtr;
  double * fNegEquBraVarPtr;

  //  The following jacobian pointers are only neccessary for 2-level newton.
  double * fPosEquPosNodePtr;
  double * fNegEquNegNodePtr;
  double * fBraEquBraVarPtr;
#endif
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;

  typedef std::vector<Instance *> InstanceVector;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &        MB,
     const FactoryBlock &      factory_block);
  ~Model ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;
  virtual bool processParams() 
  {
    return true;
  }

  virtual bool processInstanceParams() 
  {
    return true;
  }

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

  // This is the dc and transient analysis value of the source.
  double DC_TRAN;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
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

  virtual bool updateState (double * solVec, double * staVec, double * stoVec) { return true; }
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Resistor3
} // namespace Device
} // namespace Xyce

#endif

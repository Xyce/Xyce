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
// Purpose        : ThermalResistor classes
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

#ifndef Xyce_N_DEV_ThermalResistor_h
#define Xyce_N_DEV_ThermalResistor_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceInstance.h>

#include <N_DEV_Resistor.h>

namespace Xyce {
namespace Device {
namespace ThermalResistor {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, Resistor::Traits>
{
  static const char *name() {return "Resistor";}
  static const char *deviceTypeName() {return "R level 2";}
  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "R";}
  static const char *instanceDefaultParameter() {return "R";}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
//
//	This  is  the  instance class  for resistors.  It
//	contains "unique" resistor  information - ie stuff that
//	will be true of only one  resistor in the circuit, such
//	as the nodes to which it is connected.  A resistor is
//	connected to only two circuit nodes.
//
//	This class  does not directly contain information about
//	its node indices. It contains indices into the 3 parts
//	(A, dx, and  b) of the matrix  problem A*dx = b, and
//	also x.  A is the Jacobian  matrix, dx is the update to
//	the solution vector x, and b is the right hand side
//	function vector.  These indices are global, and
//	determined by topology during  the initialization stage
//	of execution.
//
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

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Riter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef );
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef) /* override */;

  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  bool processParams ();

  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool plotfileFlag () {return true;}

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  void setupPointers();

  bool outputPlotFiles(bool force_final_output);

  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  // user-specified paramters:
  double R;  // resistance  (ohms)
  double multiplicityFactor;     // multiplicity factor (M)
  // these are for the semiconductor resistor
  double length;      // resistor length.
  double width;      // resistor width.

  double area;        // resistor width.
  double thermalLength; // Length of material thermally coupled to resistor.
  double thermalArea;   // Width of material thermally coupled to resistor.

  // these four params are copied from the model:
  double resistivity;    // material resistivity
  double density;        // material density
  double heatCapacity;   // conductor volumetric heat capacity
  double thermalHeatCapacity;   // volumetric heat capacity of material thermally coupled to resistor

  double temp;   // temperature of this instance

  // derived parameters:
  double G;  // conductance (1.0/ohms)
  double i0; // current (ohms)

  //Vector local index for Positive Node
  int li_Pos;
  //Vector local index for Negative Node
  int li_Neg;

  bool tempModelEnabled;
  bool outputInternalVarsFlag;
  int li_TempState;
  int li_branch_data;         ///< Index for Lead Current and junction voltage (for power calculations)

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;

  // Pointers for Jacobian
  double *f_PosEquPosNodePtr;
  double *f_PosEquNegNodePtr;
  double *f_NegEquPosNodePtr;
  double *f_NegEquNegNodePtr;
};


//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
//
//
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

  // Semiconductor resistor parameters
  double tempCoeff1;   // first order temperature coeff.
  double tempCoeff2;   // second order temperature coeff.
  double tempCoeffExp;  // exponential temperature coeff.
  bool   tempCoeffExpModelGiven;
  double sheetRes;   // sheet resistance
  double resistanceMultiplier;   //  resistance multiplier

  double resistivity;    // material resistivity
  double density;        // material density
  double heatCapacity;   // conductor volumetric heat capacity
  double thermalHeatCapacity;   // volumetric heat capacity of material thermally coupled to resistor
  double defArea;        // default area
  double defLength;      // default length

  double defWidth;   // default width
  double narrow;     // narrowing due to side etching
  double tnom;       // parameter measurement temperature
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
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace ThermalResistor
} // namespace Device
} // namespace Xyce

#endif


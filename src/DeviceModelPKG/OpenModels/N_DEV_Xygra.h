//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Xygra classes.
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

#ifndef Xyce_N_DEV_Xygra_h
#define Xyce_N_DEV_Xygra_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <Sacado.hpp>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_CompositeParam.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : XygraCoilData
// Purpose       : This is class is a CompositeParameter type for managing
//                 coil vector-composite data
// Special Notes :
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
class XygraCoilData : public CompositeParam
{
  friend class ParametricData<XygraCoilData>;

public:
  static ParametricData<XygraCoilData> &getParametricData();

  XygraCoilData();

  void processParams();
  friend std::ostream & operator<<(std::ostream & os, const XygraCoilData & xcd);

private:
  std::string name;
  int numWindings;

public:
  std::string getName() const { return name;};
  int getNumWindings() const { return numWindings;};
};


namespace Xygra {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Xygra *DEPRECATED*";}
  static const char *deviceTypeName() {return "Xygra level 1";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 1000;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Xygra device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the ResistorInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;

  typedef Sacado::Fad::DFad<double> XygraFadType;

public:
  Instance(
     const Configuration &       configuration,
     const InstanceBlock &            IB,
     Model & Miter,
     const FactoryBlock &factory_block);

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

  bool getVoltages(std::vector<double> &voltageValues);
  bool setConductances(const std::vector< std::vector<double> > &conductanceMatrix);
  bool setK(const std::vector< std::vector<double> > &kMatrix, const double t=0);
  bool setSources(const std::vector<double> &sourceVector,const double t=0);
  int getNumNodes();
  int getNumWindings();
  void getCoilWindings(std::vector<int> &coilWindings);
  void getCoilNames(std::vector<std::string> &coilNames);

  void varTypes( std::vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  CompositeParam *constructComposite (const std::string &, const std::string &);

protected:
private:
  void setupJacStamp_();
  void interpolateSandK_();

public:
  // iterator reference to the Xygra model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  std::map<std::string, XygraCoilData *> coilDataMap;

private:
  // parameter variables

  // state variables
  // This device has no state

  // local state indices (offsets)
  // This device has no state

  // local solution indices (offsets)
  // This device uses an array of li_ values instead of individually named
  // variables.
  std::vector<int> li_Nodes_;

  // Matrix equation index variables:

  // Offset variables.  Again, this device uses an array instead of
  // discrete variables.
  // A_Equ_NodeOffests[equation][node] is the offset for node in
  // equation
  std::vector< std::vector<int> > A_Equ_NodeOffsets_;

  std::vector< std::vector<int> > jacStamp_;

  // These guys hold the Alegra input
  std::vector< std::vector<double> > theConductanceMatrix_;
  std::vector< std::vector<double> > theKMatrix_;
  std::vector< std::vector<double> > k0_;
  std::vector< std::vector<double> > k1_;
  std::vector<double> theSourceVector_;
  std::vector<double> s0_;
  std::vector<double> s1_;
  // times that (s0,k0) and (s1,k1) apply to.
  double t0_;
  double t1_;

  // For vector composite:
  std::vector<XygraCoilData*> coilDataVec;
  // total number of coils
  int nCoils;
  // number of windngs in each coil
  std::vector<int> nWindings;
  // names of each coil
  std::vector<std::string> coilNames;
  // sum over coils of number of windings per coil
  int totalNumWindings;
  // offsets into global node array of start of each coil's external vars
  std::vector<int> coilExtStart;
  // offsets into global node array of start of each coil's Internal vars
  std::vector<int> coilIntStart;
  // vector of pairs of nodes (pos,neg) for every winding
  std::vector<std::pair<int,int> > windingNodes;

  // For computation of RHS/F vector and jacobian/dFdX
  // We copy solution vars here so we can differentiate w.r.t them.
  std::vector<XygraFadType> solutionVars;
  // This is the vector of winding dv's
  std::vector<XygraFadType> dV;
  // This is the vector of winding currents
  std::vector<XygraFadType> windingCurrents;
  // and finally the vector of contributions into F:
  std::vector<XygraFadType> fContributions;
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

  static int numOrig;
  static int numSer;

  // Additional Implementation Declarations
};

//----------------------------------------------------------------------------
// Function      : Instance::getNumNodes
// Purpose       : Return the number of nodes in a given instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/2008
//----------------------------------------------------------------------------
inline int Instance::getNumNodes()
{
  return numExtVars+numIntVars;
}
//----------------------------------------------------------------------------
// Function      : Instance::getNumWindings()
// Purpose       : Return the number of windings in a given instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/2008
//----------------------------------------------------------------------------
inline int Instance::getNumWindings()
{
  return totalNumWindings;
}

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Resistor
} // namespace Device
} // namespace Xyce

#endif

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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 12/11/09
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_TransLine_h
#define Xyce_N_DEV_TransLine_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>


#define TRANS_MOD_RLC 1
#define TRANS_MOD_LC  2

namespace Xyce {
namespace Device {
namespace TransLine {

// ---------- Forward Declarations ----------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Lumped Transmission Line";}
  static const char *deviceTypeName() {return "YTRANSLINE level 1";}
  static int numNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

// note that V1 and V3 are shared by neighbor nodes.
struct lumpData
{
public:
  int indexV1;
  int indexV2;
  int indexI;
  int indexV3;

  int li_V1;
  int li_V2;
  int li_I;
  int li_V3;

  int offset_v1_v2m1; // from previous lump
  int offset_v1_iim1; // from previous lump
  int offset_v1_v1;
  int offset_v1_ii;

  int offset_v2_v2;
  int offset_v2_ii;
  int offset_v2_v3;

  int offset_ii_v1;
  int offset_ii_v2;
  int offset_ii_ii;

  int offset_v3_v2; // only used if second external node
  int offset_v3_v3;
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
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
     Model &                   Citer,
     const FactoryBlock &      factory_block);


  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  bool isLinearDevice() const { return true; }

  // Additional Public Declarations
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp) { return true; }

  bool updateIntermediateVars () { return true; }
  bool updatePrimaryState () { return true; }

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

public:
  // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // Data Members for Class Attributes

  // User-specified parameters:
  int numLumps;
  int numTransLineVars;

  double length;

  bool numLumpsGiven;
  bool lengthGiven;

  double L,C,G,R;

  std::vector<lumpData>lumpVec;

  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
class Model  : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;

public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &      MB,
     const FactoryBlock &    factory_block);
  ~Model   ();

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
  // input parameters:
  int elevNumber;
  double resist;
  double induct;
  double conduct;
  double capac;

  bool elevNumberGiven;
  bool resistGiven;
  bool inductGiven;
  bool conductGiven;
  bool capacGiven;

  int specialCase;     // what kind of model (RC, RLC, RL, ...)

  std::vector<Instance*> instanceContainer;

private:

};

//-----------------------------------------------------------------------------
// Class         : Xyce::Device::TransLine::Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// TransLine master
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
/// For efficiency, the TransLine class reimplements these methods to do the
/// work directly, instead of calling instance-level functions.
///
class Master : public DeviceMaster<Traits>
{
  friend class Instance;                            ///< Don't force a lot of pointless getters
  friend class Model;                               ///< Don't force a lot of pointless getters

public:
  ///
  /// Construct a TransLine Device.
  ///
  /// @param configuration
  /// @param factory_block
  /// @param solver_state
  /// @param device_options
  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
    : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options),
    separateInstances_(false)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec) { return true; }
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec,
                               double * leadF, double * leadQ, double * junctionV)
  { return loadDAEVectors( solVec, fVec, qVec, bVec, leadF, leadQ, junctionV, ALL ); }
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec,
                               double * leadF, double * leadQ, double * junctionV, int loadType);
  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx)
  { return loadDAEMatrices( dFdx, dQdx, ALL ); }
  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType);

private:
  bool separateInstances_;
  InstanceVector      linearInstances_;            ///< List of owned linear transline instances
  InstanceVector      nonlinearInstances_;         ///< List of owned nonlinear transline instances
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace TransLine
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_TransLine_h


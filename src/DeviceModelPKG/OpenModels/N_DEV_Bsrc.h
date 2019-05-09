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
// Purpose        : General expression dependent source
//
// Special Notes  :
//
// Creator        : Robert J Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/05/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Bsrc_h
#define Xyce_N_DEV_Bsrc_h

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {
namespace Bsrc {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Expression Based Voltage or Current Source";}
  static const char *deviceTypeName() {return "B level 1";}
  static int numNodes() {return 2;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
     Model &                   BMiter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs(const std::vector<int> & intLIDVecRef,
                    const std::vector<int> & extLIDVecRef );
  void registerStateLIDs(const std::vector<int> & staLIDVecRef);
  void registerStoreLIDs(const std::vector<int> & stoLIDVecRef );
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector<std::string> & getDepSolnVars();


  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateTemperature(const double & temp);
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

public:
  // iterator reference to the bsrc model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  Util::Expression * Exp_ptr;

  int            expNumVars;
  int            expBaseVar;
  int            expNumDdt;
  std::list<std::string>   evnList;

  std::vector<double> expVarDerivs;
  std::vector<double> myVarVals;
  std::vector<double> ddtVals;
  double         expVal;

  InstanceBlock IB;

  // flag for voltage src, needs current variable
  bool isVSRC;

  // Value of voltage or current expression
  double V;
  double I;

  double temp;
  
  // scale factor
  double scale;
  int nlstep;

  // indices into state vector:
  std::vector<int>    li_ddt;

  // solution vector indices:
  // rhs vector indices:
  int li_Pos;
  int li_Neg;
  int li_Bra;
  int li_branch_data;   ///< Index for lead current and junction voltage (for power calculations)
  // if it is not part of the solution vector

  // Local offset variables for all of the above index variables.
  int ABraEquPosNodeOffset;
  int ABraEquNegNodeOffset;
  int APosEquBraVarOffset;
  int ANegEquBraVarOffset;


 
  double rc_const;
  std::vector<int> APosEquExpVarOffsets;
  std::vector<int> ANegEquExpVarOffsets;
  std::vector<int> ABraEquExpVarOffsets;

   

  int APosEquPosNodeOffset;
  int ANegEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquNegNodeOffset;

  bool newABM;
 // Local offset variables for all of the above index variables.
  double * fBraEquPosNodePtr;
  double * fBraEquNegNodePtr;
  double * fPosEquBraVarPtr;
  double * fNegEquBraVarPtr;

  std::vector<double *> fPosEquExpVarPtrs;
  std::vector<double *> fNegEquExpVarPtrs;
  std::vector<double *> fBraEquExpVarPtrs;

  std::vector< std::vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
  ~Model ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;
  virtual bool processParams();
  virtual bool processInstanceParams();

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


  // This is the dc and transient analysis value of the
  // source.
  double DC_TRAN;

  // This is the AC magnitude
  double ACMAG;

  // This is the AC phase.
  double ACPHASE;

  // This parameter is part of the specification that the
  // source has distortion inputs at a frequency of this
  // magnitude.  It is triggered by the DISTOF1 keyword in
  // the netlist.
  double F1MAG;

  // This parameter is part of the specification that the
  // source has distortion inputs at a frequency of this
  // magnitude.  It is triggered by the DISTOF2 keyword in
  // the netlist.
  double F2MAG;

  // This parameter is associated with the specification that
  // the source has a distortion input.  This is the phase
  // associated with DISTOF1.
  double F1PHASE;

  // This parameter is associated with the specification that
  // the source has a distortion input.  This is the phase
  // associated with DISTOF2.
  double F2PHASE;
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

void registerDevice(const DeviceCountMap& deviceMap = DeviceCountMap(),
                    const std::set<int>& levelSet = std::set<int>());

} // namespace Bsrc
} // namespace Device
} // namespace Xyce

#endif


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
// Purpose        : Independent current source
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

#ifndef Xyce_N_DEV_ISRC_h
#define Xyce_N_DEV_ISRC_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Param.h>
#include <N_DEV_Source.h>

#include <Teuchos_RCP.hpp>
#include <N_UTL_FFTInterface.hpp>

namespace Xyce {
namespace Device {
namespace ISRC {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Independent Current Source";}
  static const char *deviceTypeName() {return "I level 1";}
  static int numNodes() {return 2;}
  static const char *primaryParameter() {return "DCV0";}
  static const char *instanceDefaultParameter() {return "DCV0";}
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
class Instance : public SourceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;

public:
  Instance(
     const Configuration &     configuration,
     const InstanceBlock &     instance_block,
     Model &                   model,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
 
  bool isLinearDevice() const { return true; }

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef);
  void registerStateLIDs( const std::vector<int> & staLIDVecRef);

  void registerStoreLIDs(const std::vector<int> & stoLIDVecRef );
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;

  bool processParams ();

  bool updateIntermediateVars () { return true; }
  bool updatePrimaryState () { return true; }

  // load functions, residual:
  bool loadDAEQVector () { return true; }
  bool loadDAEFVector () { return true; }
  bool loadDAEBVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () { return true; }
  bool loadDAEdFdx () { return true; }

  bool loadBVectorsforAC (double * bVecReal, double * bVecImag);

  bool loadFreqBVector(double frequency,
                       std::vector<Util::FreqVecEntry>& bVec);
 
  bool calculateFDVars ();

protected:
private:

public:
  // iterator reference to the ISRC model that owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;
  static std::vector< std::vector<int> > jacStampPDE;

  Model &       model_;         //< Owning model

  // index variables:
  int li_Pos;
  int li_Neg;

  // Store variables
  int li_branch_data;   ///< Index for lead current and junction voltage (for power calculations)

  bool          HBSpecified_;
  bool          ACSpecified_;

  // Parameters
  double DCV0;
  double par0;
  double par1;
  double par2;
  double par3;
  double par4;
  double par5;
  double par6;
  double REPEATTIME;
  double T;
  double V;
  std::string DATA;
  int RB;
  int NUM;
  bool REPEAT;
  int TRANSIENTSOURCETYPE;
  bool TRANSIENTSOURCETYPEgiven;
  int ACSOURCETYPE;
  bool ACSOURCETYPEgiven;
  int DCSOURCETYPE;
  bool DCSOURCETYPEgiven;
  bool gotParams;

  double ACMAG;
  double ACPHASE;

  double mag;
  double freq, v0;
  double phase;

                   
  bool freqVarsLoaded;
  int size_;
  Teuchos::RCP<N_UTL_FFTInterface<std::vector<double> > > ftInterface_;
  std::vector<double> ftInData_, ftOutData_, iftInData_, iftOutData_;
                   
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

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
    const FactoryBlock &        factory_block,
    const SolverState &         solver_state,
    const DeviceOptions &       device_options);

  virtual bool updateState (double * solVec, double * staVec, double * stoVec) { return true; }
  virtual bool updateSecondaryState (double * staDerivVec, double * stoVec) { return true; }

  // new DAE stuff:
  // new DAE load functions, residual:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec,
                               double * leadF, double * leadQ, double * junctionV);

  // new DAE load functions, Jacobian:
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);

private:

  bool          HBSpecified_;
  bool          ACSpecified_;
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace ISRC
} // namespace Device
} // namespace Xyce

#endif


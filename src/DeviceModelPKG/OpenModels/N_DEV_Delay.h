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
// Purpose        : Delay element classes
//
// Special Notes  :
//
// Creator        : Tom Russo
//
// Creation Date  : 17 Nov 2020
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Delay_h
#define Xyce_N_DEV_Delay_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_Source.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Delay {

class Model;
class Instance;
class History;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Delay element";}
  static const char *deviceTypeName() {return "Delay level 1";}
  static int numNodes() {return 4;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IBref,
     Model &                   Viter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  bool isLinearDevice() const
  {
    if (loadLeadCurrent)
    {
      return false;
    }
    return true;
  }

  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef);
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();
  bool loadDAEBVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  bool getInstanceBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);
  void acceptStep();

  void setupPointers();

  double getMaxTimeStepSize();
  virtual bool maxTimeStepSupported () {return false;};

  void varTypes( std::vector<char> & varTypeVec );

  DeviceState * getInternalState();
  bool setInternalState( const DeviceState & state );

  bool isConverged();

private:
  bool interpVoltageFromHistory_(double t, double & voltage,
                                const double currentT, const double currentV);

public:
  // iterator reference to the vcvs model which owns this instance.
  // Getters and setters
  Model &getModel()
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  double LeadCurrent;

  double TD_;   //< time delay

  // Matrix equation index variables:

  // local indices (offsets)
  int li_Pos;
  int li_Neg;
  int li_Bra;

  int li_ContPos;
  int li_ContNeg;

  int li_branch_data;   ///< Index for lead current and junction voltage (for power calculations)

  // Offset variables for all of the above index variables.
  int ABraEquPosNodeOffset;
  int ABraEquNegNodeOffset;
  int ABraEquContPosNodeOffset;
  int ABraEquContNegNodeOffset;
  int APosEquBraVarOffset;
  int ANegEquBraVarOffset;

  // Ptr variables for all of the above index variables.
  double * f_BraEquPosNodePtr;
  double * f_BraEquNegNodePtr;
  double * f_BraEquContPosNodePtr;
  double * f_BraEquContNegNodePtr;
  double * f_PosEquBraVarPtr;
  double * f_NegEquBraVarPtr;

  std::vector<History> history_;
  double timeOld_;

  bool newBreakPoint_;
  double newBreakPointTime_;

  bool canSetBreakPoints_;

  bool useExtrapolation_;
  bool useOnlyLinearInterpolation_;
  bool lastInterpolationConverged_;  //< True if interpolation obtained solely from converged history

  double v_drop_;
};

//-----------------------------------------------------------------------------
// Class         : History
// Purpose       : Provide a structure to save internal state history
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 17 Nov 2020
//-----------------------------------------------------------------------------
class History
{
  friend class Instance;
  friend struct Traits;

public:
  History();
  History(const History &right);
  History(double t, double v);
  ~History();
  inline bool operator<(const double &test_t) const { return (t_ < test_t);};

private:
  double t_;
  double v_;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Tom Russo
// Creation Date : 14 Apr 2020
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
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Delay
} // namespace Device
} // namespace Xyce

#endif

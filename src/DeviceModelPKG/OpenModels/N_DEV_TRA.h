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
// Purpose        : Transmission line.
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

#ifndef Xyce_N_DEV_TRA_h
#define Xyce_N_DEV_TRA_h

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace TRA {

class Model;
class Instance;
class History;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Ideal Transmission Line";}
  static const char *deviceTypeName() {return "T level 1";}
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Master;
  friend struct Traits;

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
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & st0LIDVecRef );
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  bool getInstanceBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);
  void acceptStep();

  double getMaxTimeStepSize();
  virtual bool maxTimeStepSupported () {return false;};

  DeviceState * getInternalState();
  bool setInternalState( const DeviceState & state );

  bool isConverged();

private:
  void pruneHistory(double t);
  void InterpV1V2FromHistory(double t, double * v1p , double * v2p);

public:
  // Getters and setters
  Model &getModel()
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  double Z0;   // Characteristic impedence
  double ZO;
  double G0;   // Conductance
  double td;   // Time delay  (= NL/freq if given that way)
  double freq; // frequency
  double NL;   // Normalized length
  // Flags
  // not supporting initial condition just yet... this is An Issue, I think
  // double IC_V12, IC_V34, IC_I1, IC_I2;

  bool DCMODE;
  // local indices (offsets)
  int li_Pos1;
  int li_Neg1;
  int li_Int1;
  int li_Ibr1;
  int li_Pos2;
  int li_Neg2;
  int li_Int2;
  int li_Ibr2;

  // indices into store vec for lead currents if needed
  int li_branch_data_1;
  int li_branch_data_2;


  // Matrix elements
  // Matrix equation offset variables
  int APos1EquPos1NodeOffset;
  int APos1EquInt1NodeOffset;
  int AInt1EquPos1NodeOffset;
  int AInt1EquInt1NodeOffset;
  int AInt1EquIbr1NodeOffset;
  int ANeg1EquIbr1NodeOffset;
  int AIbr1EquInt1NodeOffset;
  int AIbr1EquNeg1NodeOffset;
  // for DC simulations these 6 pairs get filled because v1 and v2 reduce
  // to V1=Ibr2*Z0 and V2=Ibr1*Z0 and no time delay
  int APos2EquPos2NodeOffset;
  int APos2EquInt2NodeOffset;
  int AInt2EquPos2NodeOffset;
  int AInt2EquInt2NodeOffset;
  int AInt2EquIbr2NodeOffset;
  int ANeg2EquIbr2NodeOffset;
  int AIbr2EquInt2NodeOffset;
  int AIbr2EquNeg2NodeOffset;
  int AIbr1EquPos2NodeOffset;
  int AIbr1EquNeg2NodeOffset;
  int AIbr1EquIbr2NodeOffset;
  int AIbr2EquPos1NodeOffset;
  int AIbr2EquNeg1NodeOffset;
  int AIbr2EquIbr1NodeOffset;

  double Vpos1,Vpos2,Vneg1,Vneg2,Vint1,Vint2,Ibr1,Ibr2; // solution vars
  double last_t;     // "primary" state variables
  double v1;    //
  double v2;    //
  bool first_BP_call_done;

  std::vector<History> history;
  double timeOld;

  bool newBreakPoint;
  double newBreakPointTime;
};

//-----------------------------------------------------------------------------
// Class         : History
// Purpose       : Provide a structure to save internal state history
// Special Notes :
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 6/14/2001
//-----------------------------------------------------------------------------
class History
{
  friend class Instance;
  friend struct Traits;

public:
  History();
  History(const History &right);
  History(double t, double v1, double v2);
  ~History();
  inline bool operator<(const double &test_t) const;

private:
  double t;
  double v1;
  double v2;
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
  friend class Master;
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
// Purpose       : This is class refers to the collected instances of the
//                 Transmission line device.
// Special Notes :
// Creator       : Tom Russo, ANL
// Creation Date : 26 April 2018
//-----------------------------------------------------------------------------

class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:

  Master(
     const Configuration &     configuration,
     const FactoryBlock &      factory_block,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options)
   : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec,
                            int loadType);

  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV, int loadType);

  virtual bool loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx,
                               int loadType);

  virtual bool loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec);
  virtual bool loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx);

};

//-----------------------------------------------------------------------------
// Function      : History::operator<
// Purpose       : compare used in the lower_bound operation.  Returns
//                 true if the time in the history object is less than
//                 the given time (double)
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 7/27/2005
//-----------------------------------------------------------------------------
inline bool History::operator<(const double &test_t) const
{
  return (t < test_t);
}

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace TRA
} // namespace Device
} // namespace Xyce

#endif

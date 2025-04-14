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

//----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 07/26/2002
//
//----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ADC_h
#define Xyce_N_DEV_ADC_h

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {
namespace ADC {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "ADC";}
  static const char *deviceTypeName() {return "YADC level 1 (Analog to Digital Interface)";};
  static int numNodes() {return 2;}
  static bool isLinearDevice() {return true;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//----------------------------------------------------------------------------
// Class          : Instance
// Purpose        : This class refers to a single instance of the ADC device.
//                  It contains indices into the matrix equation. See comments
//                  for the ResistorInstance class for more details.
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
    friend class ParametricData<Instance>;
    friend class Model;
    friend struct Traits;
    friend class Master;
    friend class DeviceInstanceOp;
    
public:
  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       instance_block,
     Model &                     model,
     const FactoryBlock &        factory_block);

  Instance(const Instance &right);

  ~Instance();

  // Additional Public Declarations
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );
  void registerStoreLIDs(const std::vector<int> & stoLIDVecRef );

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  virtual void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  bool processParams ();

  bool updateIntermediateVars () { return true; }
  bool updatePrimaryState () { return true; }

  void getTVVEC(std::vector< std::pair<double, double> > & TVVEC_Out);
  void getAndDontClearTVVEC(std::vector< std::pair<double, double> > & TVVEC_Out);
  void trimTVVEC(double earliestTime);
  bool getInstanceBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);
  int deltaVToStateVal(double deltaV);
  double deltaVTovFrac(double deltaV);
  void acceptStep();

  bool getInstanceParamsMap(std::map<std::string,double>& paramsMap);
  bool setNumberQuantLevels();
  int getNumberQuantLevels();
  bool setBitVectorWidth(int width);
  int getBitVectorWidth() {return outputBitVectorWidth_;} 

  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  Model &       model_;         //< Owning model

public:

  bool loadDAEQVector () {return true;};
  bool loadDAEFVector ();

  bool loadDAEdQdx () {return true;};
  bool loadDAEdFdx ();

  int getLIPos() const {
    return li_Pos;
  }

  int getLINeg() const {
    return li_Neg;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  // parameters:
  double R;

  // derived parameters:
  double G;  // conductance (1.0/ohms)
  double i0; // current (amps)

  // other variables:
  std::vector< std::pair<double, double> > TVVEC;
  int outputBitVectorWidth_; // number of bits on digital side
  bool outputBitVectorWidthGiven_; // width was specified on the instance line
  int nQuantLevels_;         // 2^(outputBitVectorWidth_)
  int lastOutputLevel_;      // save last value, so we know when we're going to cause a change

  // Indices into the state vector:

  //local indices (offsets)
  int li_Pos;
  int li_Neg;

  //Locally indexed offsets for jacobian
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;
  
  // store vector storage for device output state
  int li_store_output_state;
};

//----------------------------------------------------------------------------
// Function       : Model
// Purpose        :
// Special Notes  :
// Creator        : Lon Waters
// Creation Date  : 07/26/2002
//----------------------------------------------------------------------------
class Model : public DeviceModel
{
    friend class ParametricData<Model>;
    friend class Instance;
    friend struct Traits;
    friend class Master;

    typedef std::vector<Instance *> InstanceVector;

public:
  Model(
     const Configuration &      configuration,
     const ModelBlock &         model_block,
     const FactoryBlock &       factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  bool processParams ();
  bool processInstanceParams ();
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

private:
  // Model Parameters
  double lowerVoltageLimit_;
  double upperVoltageLimit_;
  double settlingTime_;
  int outputBitVectorWidth_; // number of bits on digital side

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
// Creator       : Richard Schiek, Electrical and Microsystems modeling
// Creation Date : 02/25/2009
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

  bool getBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);
};

//-----------------------------------------------------------------------------
// Function      : Instance::getNumberQuantLevels
// Purpose       : return number of possible quantization levels of an ADC
// Special Notes : Assumes that width has already been set!
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/07/2004
//-----------------------------------------------------------------------------

inline int Instance::getNumberQuantLevels()
{
  return nQuantLevels_;
}

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace ADC
} // namespace Device
} // namespace Xyce

#endif

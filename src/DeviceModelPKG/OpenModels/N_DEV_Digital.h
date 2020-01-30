//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 01/05/06
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Digital_h
#define Xyce_N_DEV_Digital_h

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace Digital {

class Model;
class Instance;
class GateData;
class DeviceData;
class InvData;
class AndData;
class NandData;
class OrData;
class NorData;
class AddData;
class XorData;
class NxorData;
class DffData;
class JkffData;
class TffData;
class DltchData;
class BufData;

// this enum is in the digital namespace so it can be used by all of
// Gate classes. NOT is deprecated now, and replaced by INV.
enum gType {INV, NOT, AND, NAND, OR, NOR, ADD, XOR, NXOR, DFF, JKFF, TFF, DLTCH, BUF};

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Behavioral Digital";}
  static const char *deviceTypeName() {return "Digital level 1";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 20;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This class refers to a single instance of a digital
//                 device.  It contains indicies into the matrix equation.
//                 See the comments for the ResistorInstance class for
//                 more details.
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class GateData;
  friend class AndData;
  friend class NandData;
  friend class OrData;
  friend class NorData;
  friend class AddData;
  friend class XorData;
  friend class NxorData;
  friend class DffData;
  friend class JkffData;
  friend class TffData;
  friend class DltchData;
  friend class BufData;

public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &       instance_block,
     Model &                     model,
     const FactoryBlock &        factory_block);


  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  // Additional Public Declarations
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateIntermediateVars () { return true; };
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool getInstanceBreakPoints(std::vector<Util::BreakPoint> &);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  // added for debug.  Allows tracking of digital device transitions since
  // acceptStep() is only executed after a time-step is accepted by time integration
  void acceptStep();

public:
  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // Data Members for Class Attributes
  // For more explanation of these attributes
  //   see the resistor classes.

  // state variables:

  std::vector<double> qlo;      // charge in the capacitor
  std::vector<double> ilo;      // current throught the capacitor
  std::vector<double> vcaplo;   // voltage drop across capacitor
  std::vector<double> qhi;
  std::vector<double> ihi;
  std::vector<double> vcaphi;
  std::vector<double> qref;
  std::vector<double> iref;
  std::vector<double> vcapref;

  std::vector<double> rilo;
  std::vector<double> rihi;
  std::vector<double> riref;
  std::vector<double> currentOut;
  std::vector<double> currentIn;

  std::vector<double> glo;
  std::vector<double> ghi;

  std::vector<double> qInp;      // charge in the capacitor
  std::vector<double> iInp;      // current throught the capacitor
  std::vector<double> vcapInp;   // voltage drop across capacitor

  std::vector<double> currentInp;

  // input params:

  bool ic1;
  bool ic2;
  bool ic3;

  int numInput;       // Number of input leads
  int numOutput;      // Number of output leads
  enum gType gate;

  //local id's (offsets)
  int li_Lo;
  int li_Hi;
  int li_Ref;
  std::vector<int> li_Inp;
  std::vector<int> li_Out;

  // Input state vars
  std::vector<int> li_currentStateInp;
  std::vector<int> li_transitionTimeInp;
  std::vector<int> li_QinpState;
  std::vector<int> li_IinpState;

  // Output state vars
  std::vector<int> li_currentStateOut;
  std::vector<int> li_transitionTimeOut;
  std::vector<int> li_QloState;
  std::vector<int> li_IloState;
  std::vector<int> li_QhiState;
  std::vector<int> li_IhiState;

  std::vector<bool> inpL;
  std::vector<double> iTime;
  std::vector<bool> outL;
  std::vector<double> oTime;

  double breakTime;

  // Offsets for Jacobian
  int row_Lo;
  int row_Hi;
  int row_Ref;
  std::vector< std::vector<int> > li_jac_Ref;
  std::vector< std::vector<int> > li_jac_Lo;
  std::vector< std::vector<int> > li_jac_Hi;

  std::vector< std::vector<int> > jacStamp;

  // added for debug purposes.  Allows for tracking of input state changes
  double prevInputStateChangeTime_;
  bool inputStateChange_;

  //Genie 110812. change state var
  //std::vector<bool> changeState;

  // added for compatibility with PSpice DIGINITSTATE
  int digInitState_;
  bool supportsXState_;
  std::vector<bool> icGiven_;

  // added for separating each gate type into its own class
  GateData * gateInfo_;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
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
     const ModelBlock &          MB,
     const FactoryBlock &        factory_block);
  ~Model   ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  bool processParams ();
  bool processInstanceParams ();
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;

  virtual std::ostream &printOutInstances(std::ostream &os) const;

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

  // Input Parameters
  double vlo;
  double vhi;
  double vref;
  double clo;
  double chi;
  double cload;
  double rload;
  double s0rlo;
  double s0rhi;
  double s0tsw;
  double s0vlo;
  double s0vhi;
  double s1rlo;
  double s1rhi;
  double s1tsw;
  double s1vlo;
  double s1vhi;
  double delay;

  // Dependent Parameters

  double gload;
};


//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Genie Hsieh, SNL, Parallel Computational Sciences
// Creation Date : 10/23/12
//-----------------------------------------------------------------------------
class Master : public DeviceMaster<Traits>
{
  friend class Instance;
  friend class Model;

public:
  Master(
     std::vector< std::pair<std::string,double> > & parNames,
     const Configuration &       configuration,
     const FactoryBlock &        factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

};

//-----------------------------------------------------------------------------
// Class         : DeviceData
// Purpose       : This class implements the functions that vary based on 
//                 the device type (U or Y).  It will be deprecated once the
//                 Y Digital Devices are removed from Xyce 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/07/15
//-----------------------------------------------------------------------------
class DeviceData
{
public:
  DeviceData(const char devLetter_);
  virtual ~DeviceData();
  char getDeviceLetter() {return devLetter_;}
  virtual void updatePowerPinLi(int&, int&, int&, int&, const bool, const bool, const bool);
  
private:
  DeviceData(const DeviceData &right);
  DeviceData &operator=(const DeviceData &right);

protected:
  const char devLetter_; // device type (U or Y)
};

//-----------------------------------------------------------------------------
// Class         : UDeviceData
// Purpose       : This class implements functions specific to the U Device 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/07/15
//-----------------------------------------------------------------------------
class UDeviceData : public DeviceData
{
public:
  UDeviceData(const char devLetter_);
  ~UDeviceData();
  void updatePowerPinLi(int&, int&, int&, int&, const bool, const bool, const bool);

private:
  UDeviceData(const UDeviceData &right);
  UDeviceData &operator=(const UDeviceData &right);
};

//-----------------------------------------------------------------------------
// Class         : YDeviceData
// Purpose       : This class implements functions specific to the Y Device 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/07/15
//-----------------------------------------------------------------------------
class YDeviceData : public DeviceData
{
public:
  YDeviceData(const char devLetter_);
  ~YDeviceData();
  void updatePowerPinLi(int&, int&, int&, int&, const bool, const bool, const bool);

private:
  YDeviceData(const YDeviceData &right);
  YDeviceData &operator=(const YDeviceData &right);
};

//-----------------------------------------------------------------------------
// Class         : GateData
// Purpose       : This class implements the functions that vary based on 
//                 the gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class GateData
{
  //friend class Instance;

public:
  GateData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  virtual ~GateData();
  int getNumInput() {return numInput_;}
  int getNumOutput() {return numOutput_;}
  gType getType() {return type_;}
  bool getSupportsXState() {return supportsXState_;};
  // getNumIO may be deprecated
  void getNumIO(int&, int&);

  void virtual checkErrors(const Instance&, const InstanceBlock&, const int&, const int&);
  bool virtual isClockLine(const int);
  void virtual setIC(Instance&, const int);
  void virtual evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
			      std::vector<double>&, const double, const double,
                              const bool, const bool, const std::vector<bool>&);
private:
  GateData(const GateData &right);
  GateData &operator=(const GateData &right);

protected:
  const std::string gateType_; // gate type (NOT, AND, ...)
  const char devLetter_; // device type (U or Y)
  const int ilNumInput_; // number of inputs found on the instance line
  int numInput_;
  int numOutput_;
  gType type_;
  bool supportsXState_;
};

//-----------------------------------------------------------------------------
// Class         : InvData
// Purpose       : This class implements functions specific to the INV and NOT
//                 gate types 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class InvData : public GateData
{
  friend class Instance;
public:
  InvData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~InvData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
                      const bool, const bool, const std::vector<bool>&);
  void checkErrors(const Instance&, const InstanceBlock&, const int&, const int&); 

private:
  InvData(const InvData &right);
  InvData &operator=(const InvData &right);
};

//-----------------------------------------------------------------------------
// Class         : AndData
// Purpose       : This class implements functions specific to the AND gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class AndData : public GateData
{
  friend class DeviceInstance;

public:
  AndData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~AndData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
                      const bool, const bool, const std::vector<bool>&);
  void checkErrors(const Instance&, const InstanceBlock&, const int&, const int&);

private:
  AndData(const AndData &right);
  AndData &operator=(const AndData &right);
};

//-----------------------------------------------------------------------------
// Class         : NandData
// Purpose       : This class implements functions specific to the NAND gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class NandData : public GateData
{
public:
  NandData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~NandData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
		      const bool, const bool, const std::vector<bool>&);
  void checkErrors(const Instance&, const InstanceBlock&, const int&, const int&);

private:
  NandData(const NandData &right);
  NandData &operator=(const NandData &right);
};

//-----------------------------------------------------------------------------
// Class         : OrData
// Purpose       : This class implements functions specific to the OR gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class OrData : public GateData
{
public:
  OrData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~OrData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
                      const bool, const bool, const std::vector<bool>&);
  void checkErrors(const Instance&, const InstanceBlock&, const int&, const int&);

private:
  OrData(const OrData &right);
  OrData &operator=(const OrData &right);
};

//-----------------------------------------------------------------------------
// Class         : NorData
// Purpose       : This class implements functions specific to the NOR gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class NorData : public GateData
{
public:
  NorData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~NorData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
                      const bool, const bool, const std::vector<bool>&);
  void checkErrors(const Instance&, const InstanceBlock&, const int&, const int&);

private:
  NorData(const NorData &right);
  NorData &operator=(const NorData &right);
};

//-----------------------------------------------------------------------------
// Class         : AddData
// Purpose       : This class implements functions specific to the ADD gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class AddData : public GateData
{
public:
  AddData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~AddData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
		      const bool, const bool, const std::vector<bool>&);

private:
  AddData(const AddData &right);
  AddData &operator=(const AddData &right);
};

//-----------------------------------------------------------------------------
// Class         : XorData
// Purpose       : This class implements functions specific to the XOR gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class XorData : public GateData
{
public:
  XorData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~XorData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
                      const bool, const bool, const std::vector<bool>&);

private:
  XorData(const XorData &right);
  XorData &operator=(const XorData &right);
};

//-----------------------------------------------------------------------------
// Class         : NxorData
// Purpose       : This class implements functions specific to the NXOR gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class NxorData : public GateData
{
public:
  NxorData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~NxorData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
                      const bool, const bool, const std::vector<bool>&);

private:
  NxorData(const NxorData &right);
  NxorData &operator=(const NxorData &right);
};

//-----------------------------------------------------------------------------
// Class         : DffData
// Purpose       : This class implements functions specific to the DFF gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class DffData : public GateData
{
public:
  DffData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~DffData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double, 
                      const bool, const bool, const std::vector<bool>&);
  bool isClockLine(const int); 
  void setIC(Instance&, const int);

private:
  DffData(const DffData &right);
  DffData &operator=(const DffData &right);
  int clockPin_; // used within the isClockLine() function
};

//-----------------------------------------------------------------------------
// Class         : JkffData
// Purpose       : This class implements functions specific to the JKFF gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 07/21/16
//-----------------------------------------------------------------------------
class JkffData : public GateData
{
public:
  JkffData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~JkffData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double, 
                      const bool, const bool, const std::vector<bool>&);
  bool isClockLine(const int); 
  void setIC(Instance&, const int);

private:
  JkffData(const JkffData &right);
  JkffData &operator=(const JkffData &right);
  int clockPin_; // used within the isClockLine() function
};

//-----------------------------------------------------------------------------
// Class         : TffData
// Purpose       : This class implements functions specific to the TFF gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 07/21/16
//-----------------------------------------------------------------------------
class TffData : public GateData
{
public:
  TffData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~TffData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double, 
                      const bool, const bool, const std::vector<bool>&);
  bool isClockLine(const int); 
  void setIC(Instance&, const int);

private:
  TffData(const TffData &right);
  TffData &operator=(const TffData &right);
  int clockPin_; // used within the isClockLine() function
};


//-----------------------------------------------------------------------------
// Class         : DltchData
// Purpose       : This class implements functions specific to the DLTCH gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class DltchData : public GateData
{
public:
  DltchData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~DltchData();

  void setIC(Instance&, const int);
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double, 
                      const bool, const bool, const std::vector<bool>&);

private:
  DltchData(const DltchData &right);
  DltchData &operator=(const DltchData &right);
};

//-----------------------------------------------------------------------------
// Class         : BufData
// Purpose       : This class implements functions specific to the Buf gate type 
//
// Special Notes :
// Creator       : Pete Sholander
// Creation Date : 08/06/15
//-----------------------------------------------------------------------------
class BufData : public GateData
{
public:
  BufData(const std::string gateType_, const char devLetter_, const int ilNumInput_);
  ~BufData();
  void evalTruthTable(const std::vector<bool>, std::vector<bool>&, 
                      std::vector<double>&, const double, const double,
		      const bool, const bool, const std::vector<bool>&);

private:
  BufData(const BufData &right);
  BufData &operator=(const BufData &right);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace Digital
} // namespace Device
} // namespace Xyce

#endif

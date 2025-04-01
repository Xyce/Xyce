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

//-----------------------------------------------------------------------------
//
// Purpose        : IBIS (I/O Buffer Information Specification) device model
//
// Special Notes  : This is the YIBIS device
//
// Creator        : Peter Sholander, SNL, Electrical Models & Simulation
//
// Creation Date  : 06/01/18
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_IBIS_h
#define Xyce_N_DEV_IBIS_h

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
namespace IBIS {

enum specType {IBIS_TYP, IBIS_MIN, IBIS_MAX, IBIS_SPEC_INVALID};

enum IBISModelType{IBIS_INPUT, IBIS_OUTPUT, IBIS_IO, IBIS_3STATE, IBIS_OPEN_DRAIN,
                   IBIS_IO_OPEN_DRAIN, IBIS_OPEN_SINK, IBIS_IO_OPEN_SINK,
                   IBIS_OPEN_SOURCE, IBIS_IO_OPEN_SOURCE,
                   IBIS_INPUT_ECL, IBIS_OUTPUT_ECL, IBIS_IO_ECL,
                   IBIS_3STATE_ECL, IBIS_TERMINATOR, IBIS_SERIES,
                   IBIS_SERIES_SWITCH, IBIS_INPUT_DIFF, IBIS_OUTPUT_DIFF,
                   IBIS_IO_DIFF, IBIS_3STATE_DIFF,IBIS_MODEL_INVALID};

enum IBISModelPolarity{IBIS_POLARITY_INVERTING, IBIS_POLARITY_NONINVERTING, IBIS_POLARITY_INVALID };

//-----------------------------------------------------------------------------
// Class         : tmmParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/027/18
//-----------------------------------------------------------------------------
struct tmmParam
{
  tmmParam()
    : type_(IBIS_SPEC_INVALID),
      val_(0),
      given_(false)
  {}

  tmmParam(const tmmParam &tp)
    : type_(tp.type_),
      val_(tp.val_),
      given_(tp.given_)
  {}

  tmmParam &operator=(const tmmParam & tp)
  {
    if (this != &tp)
    {
      type_ = tp.type_;
      val_ = tp.val_;
      given_ = tp.given_;
    }
    return *this;
  }

public:
  virtual ~tmmParam()
  {}

  IBIS::specType type_;
  double val_;
  bool given_;
};

//-----------------------------------------------------------------------------
// Class         : pkgRLC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/06/18
//-----------------------------------------------------------------------------
struct pkgRLC
{
  pkgRLC()
    : type_(IBIS_SPEC_INVALID),
      R_pkg(0),
      L_pkg(0),
      C_pkg(0)
  {}

  pkgRLC(const pkgRLC &pkg_RLC)
    : type_(pkg_RLC.type_),
      R_pkg(pkg_RLC.R_pkg),
      L_pkg(pkg_RLC.L_pkg),
      C_pkg(pkg_RLC.C_pkg)
  {}

  pkgRLC &operator=(const pkgRLC & pkg_RLC)
  {
    if (this != &pkg_RLC)
    {
      type_ = pkg_RLC.type_;
      R_pkg = pkg_RLC.R_pkg;
      L_pkg = pkg_RLC.L_pkg;
      C_pkg = pkg_RLC.C_pkg;
    }
    return *this;
  }

public:
  virtual ~pkgRLC()
  {}

  IBIS::specType type_;
  double R_pkg;
  double L_pkg;
  double C_pkg;
};

//-----------------------------------------------------------------------------
// Class         : Pin
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/06/18
//-----------------------------------------------------------------------------
struct Pin
{
  Pin()
  : pinNum(-1),
    signal_name(""),
    model_name(""),
    R_pin(0),
    L_pin(0),
    C_pin(0),
    R_pin_given(false),
    L_pin_given(false),
    C_pin_given(false)
  {}

  Pin(const Pin &p)
  : pinNum(p.pinNum),
    signal_name(p.signal_name),
    model_name(p.model_name),
    R_pin(p.R_pin),
    L_pin(p.L_pin),
    C_pin(p.C_pin),
    R_pin_given(p.R_pin_given),
    L_pin_given(p.L_pin_given),
    C_pin_given(p.C_pin_given)
  {}

  Pin &operator=(const Pin & p)
  {
    if (this != &p)
    {
      pinNum = p.pinNum;
      signal_name = p.signal_name;
      model_name = p.model_name;
      R_pin = p.R_pin;
      L_pin = p.L_pin;
      C_pin = p.C_pin;
      R_pin_given = p.R_pin_given;
      L_pin_given = p.L_pin_given;
      C_pin_given = p.C_pin_given;
    }
    return *this;
  }

public:
  virtual ~Pin()
  {}

  int pinNum;
  std::string signal_name;
  std::string model_name;
  double R_pin;
  double L_pin;
  double C_pin;
  bool R_pin_given;
  bool L_pin_given;
  bool C_pin_given;
};

//-----------------------------------------------------------------------------
// Class         : ibisBufferModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/06/18
//-----------------------------------------------------------------------------
struct ibisBufferModel
{
  ibisBufferModel()
  : Model_type(IBIS_MODEL_INVALID),
    Polarity(IBIS_POLARITY_INVALID),
    Vinl(0),
    Vinl_given(false),
    Vinh(0),
    Vinh_given(false),
    Vmeas(0),
    Vmeas_given(false),
    Rref(0),
    Rref_given(false),
    Vref(0),
    Vref_given (false),
    C_comp(),
    Voltage_Range(),
    Temperature_Range(),
    ramp_dV_dt_r_num(),
    ramp_dV_dt_r_den(),
    ramp_dV_dt_f_num(),
    ramp_dV_dt_f_den()
  {}

  ibisBufferModel(const ibisBufferModel &imp)
  : Model_type(imp.Model_type),
    Polarity(imp.Polarity),
    Vinl(imp.Vinl),
    Vinl_given(imp.Vinl_given),
    Vinh(imp.Vinh),
    Vinh_given(imp.Vinh_given),
    Vmeas(imp.Vmeas),
    Vmeas_given(imp.Vmeas_given),
    Rref(imp.Rref),
    Rref_given(imp.Rref_given),
    Vref(imp.Vref),
    Vref_given(imp.Vref_given),
    C_comp(imp.C_comp),
    Voltage_Range(imp.Voltage_Range),
    Temperature_Range(imp.Temperature_Range),
    ramp_dV_dt_r_num(imp.ramp_dV_dt_r_num),
    ramp_dV_dt_r_den(imp.ramp_dV_dt_r_den),
    ramp_dV_dt_f_num(imp.ramp_dV_dt_f_num),
    ramp_dV_dt_f_den(imp.ramp_dV_dt_f_den)
  {}

  ibisBufferModel &operator=(const ibisBufferModel & imp)
  {
    if (this != &imp)
    {
      Model_type = imp.Model_type,
        Polarity = imp.Polarity,
        Vinl = imp.Vinl,
        Vinl_given = imp.Vinl_given,
        Vinh = imp.Vinh,
        Vinh_given= imp.Vinh_given,
        Vmeas = imp.Vmeas,
        Vmeas_given = imp.Vmeas_given,
        Rref = imp.Rref,
        Rref_given = imp.Rref_given,
        Vref= imp.Vref,
        Vref_given = imp.Vref_given,
        C_comp = imp.C_comp;
      Voltage_Range = imp.Voltage_Range;
      Temperature_Range = imp.Temperature_Range;
      ramp_dV_dt_r_num = imp.ramp_dV_dt_r_num;
      ramp_dV_dt_r_den = imp.ramp_dV_dt_r_den;
      ramp_dV_dt_f_num = imp.ramp_dV_dt_f_num;
      ramp_dV_dt_f_den = imp.ramp_dV_dt_f_den;
    }
    return *this;
  }

public:
  virtual ~ibisBufferModel()
  {}

  IBISModelType Model_type;
  IBISModelPolarity Polarity;
  double Vinl;
  bool Vinl_given;
  double Vinh;
  bool Vinh_given;
  double Vmeas;
  bool Vmeas_given;
  double Rref;
  bool Rref_given;
  double Vref;
  bool Vref_given;

  std::vector<tmmParam> C_comp;
  std::vector<tmmParam> Voltage_Range;
  std::vector<tmmParam> Temperature_Range;
  std::vector<tmmParam> GND_Clamp_Reference;
  std::vector<tmmParam> ramp_dV_dt_r_num;
  std::vector<tmmParam> ramp_dV_dt_r_den;
  std::vector<tmmParam> ramp_dV_dt_f_num;
  std::vector<tmmParam> ramp_dV_dt_f_den;
};


class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "IBIS";}
  static const char *deviceTypeName() {return "IBIS level 1";}
  static int numNodes() {return 2;}
  static int numOptionalNodes() {return 20;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Peter Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/01/18
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
  bool updateIntermediateVars() { return true; }
  bool updatePrimaryState() { return true; }
  bool updateSecondaryState();

  // load functions, residual:
  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  void setupPointers();

  void varTypes( std::vector<char> & varTypeVec );

  // functions specific to the IBIS device model
  bool readIbsFile();
  bool ibisStrToVal(const std::string& valStr,
                    double& val,
                    bool& valFound,
                    const int lineNum);
  bool makeTmmVec(const IO::TokenVector& parsedLine,
                  std::vector<tmmParam>& tp,
                  const int offset,
                  const int lineNum);
  IBIS::IBISModelType setIBISModelType(std::string modelTypeStr);
  IBIS::IBISModelPolarity setIBISModelPolarity(std::string modelTypeStr);
  void splitIBISFileLine(const std::string& aLine, IO::TokenVector & parsedLine);
  void readIBISFileLine(std::istream & in, std::string& line, int& lineNum);

public:
  // iterator reference to the YIBIS model which owns this instance.
  // Getters and setters
  Model &getModel()
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // pointers for table expressions
  Util::Expression * Exp_ptr_pc;  // pointer to expression for power clamp table
  Util::Expression * Exp_ptr_gc;  // pointer to expression for ground clamp table
  Util::Expression * Exp_ptr_pu;  // pointer to expression for pullup table
  Util::Expression * Exp_ptr_pd;  // pointer to expression for pulldown table

  int            expNumVars;
  int            expNumVars_pc;
  int            expNumVars_gc;
  int            expNumVars_pu;
  int            expNumVars_pd;

  int            expNumDdt;
  int            expNumDdt_pc;
  int            expNumDdt_gc;
  int            expNumDdt_pu;
  int            expNumDdt_pd;

  std::vector<double> expVarDerivs;
  std::vector<double> expVarDerivs_pc;
  std::vector<double> expVarDerivs_gc;
  std::vector<double> expVarDerivs_pu;
  std::vector<double> expVarDerivs_pd;

  double         expVal;
  double         expVal_pc;
  double         expVal_gc;
  double         expVal_pu;
  double         expVal_pd;

  InstanceBlock IB;

  // variables specific to the IBIS device
  std::string fileName_;
  bool fileName_given;
  std::string modelName_;
  std::vector<std::string> nodeList_;
  bool modelName_given;
  char ibisCommentChar_;
  std::string ibisVer_;

  //variables common to the component
  std::string ibisComponent_;
  std::vector<pkgRLC> pkgRLCVec_;
  std::vector<Pin> pinVec_;

  // variables specific to a buffer model
  IBISModelType Model_type;
  IBISModelPolarity Polarity;
  double Vinl;
  bool Vinl_given;
  double Vinh;
  bool Vinh_given;
  double Vmeas;
  bool Vmeas_given;
  double Rref;
  bool Rref_given;
  double Vref;
  bool Vref_given;
  double C_comp;
  bool C_comp_given;
  double GND_Clamp_Reference;
  bool GND_Clamp_Reference_given;
  double ramp_dV_dt_r_num;
  bool ramp_dV_dt_r_num_given;
  double ramp_dV_dt_r_den;
  bool ramp_dV_dt_r_den_given;
  double ramp_dV_dt_f_num;
  bool ramp_dV_dt_f_num_given;
  double ramp_dV_dt_f_den;
  bool ramp_dV_dt_f_den_given;
  ibisBufferModel bufferModel_;

  // current variable for VI tables
  double gndClampI_;
  double powerClampI_;
  double pulldownI_;
  double pullupI_;

  // Value of current expressions
  double I;

  double temp;

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

  std::vector<int> APosEquExpVarOffsets;
  std::vector<int> ANegEquExpVarOffsets;
  std::vector<int> ABraEquExpVarOffsets;

  int APosEquPosNodeOffset;
  int ANegEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int ANegEquNegNodeOffset;

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
// Creator       : Peter Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/01/18
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

private:
  std::vector<Instance*> instanceContainer;

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

void registerDevice();

} // namespace IBIS
} // namespace Device
} // namespace Xyce

#endif


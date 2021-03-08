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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Diode_h
#define Xyce_N_DEV_Diode_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {
namespace Diode {

class Model;
class Instance;

template <typename ScalarT> 
inline ScalarT Xycemax ( ScalarT f1, ScalarT f2) { return f1 > f2 ? f1 : f2; }

template <typename ScalarT> 
inline ScalarT Xycemin ( ScalarT f1, ScalarT f2) { return f1 < f2 ? f1 : f2; }

/// general sensitivity functor for all model params.
class diodeSensitivity :  public baseSensitivity
{
  public:
  diodeSensitivity() : 
    baseSensitivity() {};

  virtual ~diodeSensitivity() {};

  virtual void operator()(
    const ParameterBase &entity,
    const std::string &param,
    std::vector<double> & dfdp, 
    std::vector<double> & dqdp, 
    std::vector<double> & dbdp, 
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const ;
};

static diodeSensitivity diodeSens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Diode";}
  static const char *deviceTypeName() {return "D level 1,2";}
  static int numNodes() {return 2;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};


template <typename ScalarT> 
bool processParams (
    ScalarT & M,
    ScalarT & EG,
    ScalarT & FC,
    const ScalarT & RS,
    ScalarT & COND,
    ScalarT & F2,
    ScalarT & F3
    );

template <typename ScalarT> 
bool updateTemperature 
(
 //const double & temp, Instance & instance, Model & model_,
  // instance params:
   ScalarT & Temp,
   ScalarT & tJctCap,
   ScalarT & tJctPot,
   ScalarT & tDepCap,
   ScalarT & tF1,
   ScalarT & tSatCur,
   ScalarT & tSatCurR,
   ScalarT & tVcrit,
   ScalarT & tRS,
   ScalarT & tCOND,
   ScalarT & tIRF,
   ScalarT & tIKF,
   ScalarT & tBrkdwnV,

   // model variables/params:
   const ScalarT & TNOM,
   const ScalarT & VJ,
   const ScalarT & CJO,
   const ScalarT & M,
   const ScalarT & N,
   const ScalarT & IS,
   const ScalarT & EG,
   const ScalarT & XTI,
   const ScalarT & RS,
   const ScalarT & COND,
   const ScalarT & IRF,
   const ScalarT & NR,
   const ScalarT & IKF,
   const ScalarT & TIKF,
   const ScalarT & ISR,
   const ScalarT & IBV,
   const ScalarT & BV,
   const bool & BVGiven,
   const bool & IRFGiven,
   const ScalarT & TBV1,
   const ScalarT & TBV2,
   const ScalarT & TRS1,
   const ScalarT & TRS2,
   const ScalarT & FC,
   const int  level

 );

template <typename ScalarT> 
bool updateIntermediateVars 
  (
   // inputs:
   const ScalarT & Vp,
   const ScalarT & Vpp,
   const ScalarT & Vn,
   const ScalarT & Vd,

  // instance params:
   const ScalarT & Temp,
   const ScalarT & tJctCap,
   const ScalarT & tJctPot,
   const ScalarT & tDepCap,
   const ScalarT & tF1,
   const ScalarT & tSatCur,
   const ScalarT & tSatCurR,
   const ScalarT & tVcrit,
   const ScalarT & tRS,
   const ScalarT & tCOND,
   const ScalarT & tIRF,
   const ScalarT & tIKF,
   const ScalarT & tBrkdwnV,
  
   // instance variables:
   const ScalarT & Area,
   const int & lambertWFlag,
   const double & gmin,

  // model params:
  const ScalarT M   , // grading parameter
  const ScalarT BV  , // breakdown voltage
  const ScalarT IBV , // reverse breakdown current
  const ScalarT NBV , // reverse breakdown ideality factor
  const ScalarT IBVL, // low-level reverse breakdown current
  const ScalarT NBVL, // low-level reverse breakdown ideality factor
  const ScalarT N   , // non-ideality factor.
  const ScalarT NR  , // emission coeff. for ISR.
  const ScalarT TT  , // transit time.
  const ScalarT F2  , // capacitive polynomial factor
  const ScalarT F3  , // capacitive polynomial factor

  const int  level,

  // outputs:
  ScalarT & Id,
  ScalarT & Gd,
  ScalarT & Qd,
  ScalarT & Cd,
  ScalarT & Gspr
  ); 

template <typename ScalarT> 
bool applyLimiters
( 
   DeviceSupport & devSupport,

   // inputs:
   const ScalarT & Vp,
   const ScalarT & Vpp,
   const ScalarT & Vn,

   // parameters:
   const ScalarT & tVcrit,

   // output
   ScalarT & Vd,
   ScalarT & Vd_orig,
   ScalarT & Vd_old,

   const ScalarT & currVd_old,
   const ScalarT & nextVd_old,

   const double InitCond,
   const bool InitCondGiven,
   const bool BVGiven,

   // instance vars:
   const int  off,
   bool & origFlag,

   const bool dotICapplies, // check all the "flagSol". if any == 1, then true

  // solver state variables:
   const int & newtonIter,
   const bool initJctFlag,
   const bool voltageLimiterFlag,
   const bool dcopFlag,
   const bool locaEnabledFlag
  ); 

//-----------------------------------------------------------------------------
// Class         : Instance
//
// Purpose       : This class represents a single instance  of the  diode
//                 device.  It mainly contains indices and pointers into
//                 the matrix equation (see the resistor instance class for
//                 more details).
//
//                 The diode will have 1 internal node, in addition to the
//                 nodes which are connected to it.  This is so that  it
//                 can be placed in series with a resistor that represents
//                 the resistance of intrinsic Si.
//
// Special Notes :
// Creator       : Eric Keiter, Parallel Computational Sciences
// Creation Date : 02/28/00
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class diodeSensitivity;

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
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef);
  void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  void loadNodeSymbols(Util::SymbolTable &symbol_table) const;

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature ( const double & temp = -999.0 );
  bool lambertWCurrent (double Isat, double Vte, double RS);
  bool lambertWBreakdownCurrent (double Isat, double Vte, double RS);
  bool lambertWLinearReverseBias (double Isat, double Vte, double RS);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  int getNumNoiseSources () const;
  void setupNoiseSources (Xyce::Analysis::NoiseData & noiseData);
  void getNoiseSources (Xyce::Analysis::NoiseData & noiseData);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

public:
  // iterator reference to the diode model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp_RS;
  static std::vector< std::vector<int> > jacStamp;

  static std::vector<int> jacMap_RS;
  static std::vector<int> jacMap;

  static std::vector< std::vector<int> > jacMap2_RS;
  static std::vector< std::vector<int> > jacMap2;


  Model &       model_;         //< Owning model

  int  off;
  double Area;
  double InitCond;
  double Temp;
  int lambertWFlag;
  bool InitCondGiven;

  double tJctPot;
  double tJctCap;
  double tDepCap;
  double tSatCur;
  double tVcrit;
  double tF1;
  double tBrkdwnV;
  double tSatCurR;
  double tIKF;
  double tRS;
  double tCOND;
  double tIRF;

  double Id;     //diode current
  double Gd;     //diode conductivity
  double Cd;     //depletion capacitance
  double Gcd;    //dep cap conductivity
  double Qd;     //capacitor charge
  double Icd;    //capacitor current
  double Gspr;
  //double LeadCurrent;

  double Vpp;
  double Vp;
  double Vn;
  double Vc;

  double Vd;
  double Vd_old;
  double Vd_orig;

  // end of intermediate variables

  // state variables:
  double q0;  // charge in the capacitor
  double i0;  // current throught the capacitor

  //local indices (offsets)
  // int li_QState;

  // for voltage limiting
  int li_storevd;

  // for lead current
  int li_branch_data;         ///< Index for Lead Current and junction voltage (for power calculations)

  int li_Pos;
  int li_Neg;
  int li_Pri;

  // Matrix equation local offset variables
  int APosEquPosNodeOffset;
  int APosEquPriNodeOffset;
  int ANegEquNegNodeOffset;
  int ANegEquPriNodeOffset;
  int APriEquPosNodeOffset;
  int APriEquNegNodeOffset;
  int APriEquPriNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Matrix equation local pointer variables
  double * fPosEquPosNodePtr;
  double * fPosEquPriNodePtr;
  double * fNegEquNegNodePtr;
  double * fNegEquPriNodePtr;
  double * fPriEquPosNodePtr;
  double * fPriEquNegNodePtr;
  double * fPriEquPriNodePtr;

  double * qPosEquPosNodePtr;
  double * qPosEquPriNodePtr;
  double * qNegEquNegNodePtr;
  double * qNegEquPriNodePtr;
  double * qPriEquPosNodePtr;
  double * qPriEquNegNodePtr;
  double * qPriEquPriNodePtr;
#endif


  // Flags
  bool TEMP_GIVEN;
  bool AREA_GIVEN;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;
  friend class diodeSensitivity;

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

  double IS;   // saturation current (A)
  double RS;   // ohmic resistance (ohms)
  double COND; // corresponding conductance
  double N;    // emission coefficient
  double ISR;  // recombination saturation current (A)
  double NR;   // emission coefficient for ISR
  double IKF;  // high-injection knee current (A)
  double TT;   // transit time (sec)
  double CJO;  // zero-bias junction capacitance (F)
  double VJ;   // built-in junction potential (V)
  double M;    // grading coefficient
  double EG;   // activation  energy (eV).
  //    For Si, EG = 1.11
  //        Ge, EG = 0.67
  //        Sbd, EG = 0.69
  double XTI;  // isaturation-current temp. exp
  double TIKF; // IKF temperature coeff.
  double TBV1; // BV linear temperature coeff.
  double TBV2; // BV quadratic temperature coeff.
  double TRS1; // RS linear temperature coeff.
  double TRS2; // RS quadratic temperature coeff.
  double FC;   // coefficient for forward-bias depletion capacitance
  // formula
  double BV;   // reverse breakdown voltage
  double IBV;  // current at  breakdown voltage (A)
  double IRF;  // adjustment for linear portion of reverse current
  double NBV;  // reverse breakdown ideality factor
  double IBVL; // low-level current at  breakdown voltage (A)
  double NBVL; // low-level reverse breakdown ideality factor
  double F2;
  double F3;
  double TNOM; // parameter measurement temperature (C)
  double KF;   // flicker noise coefficient
  double AF;   // flicker noise exponent

  bool BVGiven;
  bool IRFGiven;
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

  // new DAE stuff:
  // new DAE load functions, residual:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);

  // new DAE load functions, Jacobian:
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
};

void registerDevice(const DeviceCountMap& deviceMap = DeviceCountMap(),
                    const std::set<int>& levelSet = std::set<int>());

} // namespace Diode
} // namespace Device
} // namespace Xyce

#endif


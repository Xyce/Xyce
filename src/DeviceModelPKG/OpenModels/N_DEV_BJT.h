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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_BJT_H
#define Xyce_N_DEV_BJT_H

#include <Sacado_No_Kokkos.hpp>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {
namespace BJT {

class Model;
class Instance;

typedef Sacado::Fad::SFad<double, 1> fadType;

template <typename ScalarT> 
inline ScalarT Xycemax ( ScalarT f1, ScalarT f2) { return f1 > f2 ? f1 : f2; }

template <typename ScalarT> 
inline ScalarT Xycemin ( ScalarT f1, ScalarT f2) { return f1 < f2 ? f1 : f2; }

/// general sensitivity functor for all instance params.
class bjtInstanceSensitivity :  public baseSensitivity
{
  public:
  bjtInstanceSensitivity() : 
    baseSensitivity() {};

  virtual ~bjtInstanceSensitivity() {};

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

/// general sensitivity functor for all model params.
class bjtModelSensitivity :  public baseSensitivity
{
  public:
  bjtModelSensitivity() : 
    baseSensitivity() {};

  virtual ~bjtModelSensitivity() {};

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

static bjtInstanceSensitivity bjtInstanceSens;
static bjtModelSensitivity bjtModelSens;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Bipolar Junction Transistor";}
  static const char *deviceTypeName() {return "Q level 1";}
  static int numNodes() {return 3;}
  static int numOptionalNodes() {return 4;}
  static int numFillNodes() {return 1;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

template <typename ScalarT> 
bool updateTemperature( 
    // inputs
    // instance
    const ScalarT & TEMP,
  
    // model
    const ScalarT & TNOM,
    const ScalarT & energyGap,
    const ScalarT & tempExpIS,
    const ScalarT & betaExp,
    const ScalarT & potBE,
    const ScalarT & depCapBE,
    const ScalarT & juncExpBE,
    const ScalarT & depCapBC,
    const ScalarT & juncExpBC,
    const ScalarT & depCapCoeff,
    const ScalarT & satCur,
    const ScalarT & betaF,
    const ScalarT & betaR,

    const ScalarT & c2,
    const ScalarT & c4,
    bool c2Given,
    bool c4Given,

    bool leakBECurrentGiven, 
    bool leakBCCurrentGiven,

    const ScalarT & leakBEEmissionCoeff,
    const ScalarT & leakBCEmissionCoeff,

    const ScalarT & rollOffExp,
    const ScalarT & baseResist,
    const ScalarT & collectorResist,
    const ScalarT & emitterResist,
    const ScalarT & potBC,
    const ScalarT & rollOffF,
    const ScalarT & rollOffR,
    const ScalarT & earlyVoltF,
    const ScalarT & earlyVoltR,

    // outputs
    ScalarT & vt,
    ScalarT & leakBECurrent,
    ScalarT & leakBCCurrent,
    ScalarT & tBELeakCur,
    ScalarT & tBCLeakCur,
    ScalarT & tleakBEEmissionCoeff,
    ScalarT & tleakBCEmissionCoeff,

    ScalarT & tRollOffExp,
    ScalarT & tInvRollOffF,
    ScalarT & tInvRollOffR,
    ScalarT & tBaseResist,
    ScalarT & tCollectorResist,
    ScalarT & tEmitterResist,

    ScalarT & tBECap,
    ScalarT & tBEPot,
    ScalarT & tBCCap,
    ScalarT & tBCPot,
    ScalarT & tDepCap,
    ScalarT & tF1,
    ScalarT & tF4,
    ScalarT & tF5,
    ScalarT & tVCrit,
    ScalarT & tSatCur,
    ScalarT & tBetaF,
    ScalarT & tBetaR,
    ScalarT & tInvEarlyVoltF,
    ScalarT & tInvEarlyVoltR
    );

template <typename ScalarT> 
bool processParams(
   // inputs 
    bool leakBECurrentGiven,
    bool leakBCCurrentGiven,
    bool c2Given,
    bool c4Given,
    bool minBaseResistGiven,
    bool VAFgiven,
    bool IKFgiven,
    bool VARgiven,
    bool IKRgiven,
    bool VTFgiven,
    bool FCgiven,

    const ScalarT & c2,
    const ScalarT & c4,
    const ScalarT & satCur,
    const ScalarT & baseResist,
    const ScalarT & earlyVoltF,
    const ScalarT & rollOffF,
    const ScalarT & earlyVoltR,
    const ScalarT & rollOffR,

    const ScalarT & collectorResist,
    const ScalarT & emitterResist,
    const ScalarT & transTimeFVBC,
    const ScalarT & excessPhase,
    const ScalarT & transTimeF,
    const ScalarT & juncExpBE,
    const ScalarT & juncExpBC,

   //outputs
    ScalarT & leakBECurrent,
    ScalarT & leakBCCurrent,
    ScalarT & minBaseResist,
    ScalarT & invEarlyVoltF,
    ScalarT & invRollOffF,
    ScalarT & invEarlyVoltR,
    ScalarT & invRollOffR,

    ScalarT & collectorConduct,
    ScalarT & emitterConduct,
    ScalarT & transTimeVBCFac,
    ScalarT & excessPhaseFac,
    ScalarT & depCapCoeff,
 
    ScalarT & f2,
    ScalarT & f3,
    ScalarT & f6,
    ScalarT & f7 
    );

template <typename ScalarT>
void oldDAEExcessPhaseCalculation1 (
    const ScalarT & td,
    const ScalarT & qB,
    const ScalarT & iBE,
    bool dcopFlag,
    bool beginIntegrationFlag,
    double * currStaVec,
    double * lastStaVec,
    const int li_istateCEXBC 
    )
{};

template <>
void oldDAEExcessPhaseCalculation1 (
    const fadType & td,
    const fadType & qB,
    const fadType & iBE,
    bool dcopFlag,
    bool beginIntegrationFlag,
    double * currStaVec,
    double * lastStaVec,
    const int li_istateCEXBC 
    );


template <>
void oldDAEExcessPhaseCalculation1 (
    const double & td,
    const double & qB,
    const double & iBE,
    bool dcopFlag,
    bool beginIntegrationFlag,
    double * currStaVec,
    double * lastStaVec,
    const int li_istateCEXBC 
    );

template <typename ScalarT>
void oldDAEExcessPhaseCalculation2 
   (const ScalarT & td,
    const ScalarT & qB,
    const ScalarT & iBE,
    const ScalarT & gBE,

    const double dt0, //dt0 = getSolverState().currTimeStep;
    const double dt1, //dt1 = getSolverState().lastTimeStep;

    bool dcopFlag,
    bool beginIntegrationFlag,

    double * nextStaVec, // raw pointers fine here
    const double * currStaVec, // raw pointers fine here
    const double * lastStaVec, // raw pointers fine here

    const int li_istateCEXBC,

    ScalarT & iEX, 
    ScalarT & gEX, 
    ScalarT & iC_local)
{};


template <>
void oldDAEExcessPhaseCalculation2 
   (const fadType & td,
    const fadType & qB,
    const fadType & iBE,
    const fadType & gBE,

    const double dt0, //dt0 = getSolverState().currTimeStep;
    const double dt1, //dt1 = getSolverState().lastTimeStep;

    bool dcopFlag,
    bool beginIntegrationFlag,

    double * nextStaVec, // raw pointers fine here
    const double * currStaVec, // raw pointers fine here
    const double * lastStaVec, // raw pointers fine here

    const int li_istateCEXBC,

    fadType & iEX, 
    fadType & gEX, 
    fadType & iC_local);

template <>
void oldDAEExcessPhaseCalculation2 
   (const double & td,
    const double & qB,
    const double & iBE,
    const double & gBE,

    const double dt0, //dt0 = getSolverState().currTimeStep;
    const double dt1, //dt1 = getSolverState().lastTimeStep;

    bool dcopFlag,
    bool beginIntegrationFlag,

    double * nextStaVec, // raw pointers fine here
    const double * currStaVec, // raw pointers fine here
    const double * lastStaVec, // raw pointers fine here

    const int li_istateCEXBC,

    double & iEX, 
    double & gEX, 
    double & iC_local);

template <typename ScalarT> 
void auxDAECalculations (
    const ScalarT & i_fx, 
    const ScalarT & td,
    const ScalarT & iBE,
    const ScalarT & iBEleak,
    const ScalarT & iBC,
    const ScalarT & iBCleak,
    const ScalarT & qB,
    const ScalarT & invqB,
    const ScalarT & tBetaF,
    const ScalarT & tBetaR,

    const ScalarT & gBC,
    const ScalarT & gBE,

    const ScalarT & dqBdvBp,
    const ScalarT & dqBdvCp,
    const ScalarT & dqBdvEp,

    bool dcopFlag,

    ScalarT & iCE,
    ScalarT & iC,
    ScalarT & iB,
    ScalarT & iE,

    ScalarT & diCEdvBp,
    ScalarT & diCEdvCp,
    ScalarT & diCEdvEp,

    ScalarT & diBEdvBp,
    ScalarT & diBEdvCp,
    ScalarT & diBEdvEp
    );

template <typename ScalarT> 
bool updateIntermediateVars 
  (
  // inputs:
   const ScalarT & vBE,
   const ScalarT & vBC,
   const ScalarT & vBX,
   const ScalarT & vCS,

   const ScalarT & i_fx,

  // instance params:
   const ScalarT & AREA,

  // instance variables:
   const ScalarT & tSatCur,
   const ScalarT & vt,
   const ScalarT & tleakBEEmissionCoeff,
   const ScalarT & tleakBCEmissionCoeff,
   const ScalarT & tBELeakCur,
   const ScalarT & tBCLeakCur,
   const ScalarT & tInvRollOffF,
   const ScalarT & tInvRollOffR,
   const ScalarT & tRollOffExp,
   const ScalarT & tInvEarlyVoltF,
   const ScalarT & tInvEarlyVoltR,
   const ScalarT & tBECap,
   const ScalarT & tBCCap,
   const ScalarT & tDepCap,
   const ScalarT & tBEPot,
   const ScalarT & tBCPot,
   const ScalarT & tF1,
   const ScalarT & tF4,
   const ScalarT & tF5,

   // from device options:
   const double & gmin,
   const bool newExcessPhase,

   // from solver state:
   const bool dcopFlag,
   const bool tranopFlag,
   const bool acopFlag,
   const bool initTranFlag,
   const bool beginIntegrationFlag,
   const int newtonIter,
   const double pdt,

  // model params:
  const ScalarT & emissionCoeffF,
  const ScalarT & emissionCoeffR,
  const ScalarT & baseFracBCCap,
  const ScalarT & CJS,
  const ScalarT & depCapCoeff,
  const ScalarT & potBC,
  const ScalarT & transTimeF,
  const ScalarT & transTimeR,
  const ScalarT & transTimeBiasCoeffF,
  const ScalarT & transTimeVBCFac,
  const ScalarT & transTimeHighCurrF,
  const ScalarT & juncExpBE,
  const ScalarT & juncExpBC,
  const ScalarT & potSubst,
  const ScalarT & expSubst,

  const ScalarT & emitterConduct,
  const ScalarT & collectorConduct,
  const ScalarT & minBaseResist,
  const ScalarT & baseResist,
  const ScalarT & baseCurrHalfResist,

  const ScalarT & excessPhaseFac,

  // model params from processParams:
  const ScalarT & f2,
  const ScalarT & f3,
  const ScalarT & f6,
  const ScalarT & f7,

  const int  level,

  // these are here b/c of the excess phase old-DAE form
  double * nextStaVec,
  double * currStaVec,
  double * lastStaVec,
  const int li_istateCEXBC,
  const double dt0, //  getSolverState().currTimeStep, 
  const double dt1, //  getSolverState().lastTimeStep, 

  // outputs:
   ScalarT & iB,
   ScalarT & iC,
   ScalarT & iE,

   ScalarT & iBE,
   ScalarT & gBE,
   ScalarT & iBC,
   ScalarT & gBC,
   ScalarT & iCE,
   ScalarT & iBEleak,
   ScalarT & gBEleak,
   ScalarT & iBCleak,
   ScalarT & gBCleak,
   ScalarT & qB,
   ScalarT & invqB,
   ScalarT & iBEhighCurr,
   ScalarT & gBEhighCurr,

   ScalarT & capeqCB,
   ScalarT & geqCB,

   ScalarT & qBEdep,
   ScalarT & capBEdep,
   ScalarT & qBEdiff,
   ScalarT & capBEdiff,

   ScalarT & qBCdep,
   ScalarT & capBCdep,
   ScalarT & qBCdiff,
   ScalarT & capBCdiff,

   ScalarT & qBX,
   ScalarT & capBX,
   ScalarT & qCS,
   ScalarT & capCS,

   ScalarT & gEpr,
   ScalarT & gCpr,
   ScalarT & gX,

   ScalarT & diBrdvB,
   ScalarT & diBrdvCp,
   ScalarT & diBrdvEp,
   ScalarT & diBrdvBp,

   ScalarT & diCEdvEp,
   ScalarT & diCEdvCp,
   ScalarT & diCEdvBp,

   ScalarT & diBEdvBp,
   ScalarT & diBEdvCp,
   ScalarT & diBEdvEp,

   ScalarT & gBEtot,
   ScalarT & gBCtot,
   ScalarT & tBetaF,
   ScalarT & tBetaR 
  );


//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This class refers to a single instance of the BJT
//                 device.  It contains indices into the matrix equation.
//                 See the comments for the ResistorInstance class for
//                 more details.
//
//                 The bjt will have 4 external nodes: collector, base,
//                 emitter, and substrate, and 3 internal nodes:
//                 collectorPrime, basePrime, and emitterPrime.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class bjtInstanceSensitivity;
  friend class bjtModelSensitivity;

  // functions
public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   it_MB,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & stateLIDVecRef );
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef);
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  
  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature (const double & temp = -999.0 );
  bool lambertWCurrent (double &Id, double &Gd, double Vd, double Vte, double Isat);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  int getNumNoiseSources () const;
  void setupNoiseSources (Xyce::Analysis::NoiseData & noiseData);
  void getNoiseSources (Xyce::Analysis::NoiseData & noiseData);

  void loadErrorWeightMask ();

#ifdef Xyce_DEBUG_EXCESS_PHASE
  bool plotfileFlag () {return true;}
#else
  bool plotfileFlag () {return false;}
#endif

  void oldDAEExcessPhaseCalculation1 ();
  void oldDAEExcessPhaseCalculation2
  (double & iEX, double & gEX, double & iC_local);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  // Debugging Excess Phase function:
  bool outputPlotFiles(bool force_final_output);

  void setupPointers();

protected:
private:

  //attributes
public:
  //iterator reference to the BJT model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  //external instance params
  double AREA;  // The normalized emitter area (AREA)
  double icVBE; // the inital base-emitter voltage (ICVBE)
  double icVCE; // the inital collector-emitter voltage (ICVCE)
  double TEMP;  // instance temperature (TEMP)
  bool   OFF;   // initial mode of operation (OFF)
  bool lambertWFlag;
  bool IC_GIVEN;
  bool externalNodeMode;
  bool offFlag;

  //generated instance params
  double vt;           //thermal junc. volt.
  double tSatCur;      //saturation current (temp. adj.)
  double tBetaF;       //forward beta       (temp. adj.)
  double tBetaR;       //reverse beta       (temp. adj.)
  double tBELeakCur;   //BE leakage current (temp. adj.)
  double tBCLeakCur;   //BC leakage current (temp. adj.)
  double tBECap;       //BE capacitance     (temp. adj.)
  double tBCCap;       //BC capacitance     (temp. adj.)
  double tBEPot;       //BE potential       (temp. adj.)
  double tBCPot;       //BC potential       (temp. adj.)
  double tDepCap;      //join pt in diode curve (temp. adj.)
  double tF1;          //polynomial coeff.  (temp. adj.)
  double tF4;          //polynomial coeff.  (temp. adj.)
  double tF5;          //polynomial coeff.  (temp. adj.)
  double tVCrit;       //critical voltage   (temp. adj.)

  // additional temperature-adjusted parameters
  double tleakBEEmissionCoeff;
  double tleakBCEmissionCoeff;
  double tRollOffExp;
  double tInvRollOffF;
  double tInvRollOffR;
  double tInvEarlyVoltF;
  double tInvEarlyVoltR;
  double tBaseResist;
  double tCollectorResist;
  double tEmitterResist;

  //generated intermediate variables
  double vEEp;      // e-e' voltage
  double vBBp;      // b-b' voltage
  double vCCp;      // c-c' voltage

  double vBE;       // b'-e' voltage
  double vBC;       // b'-c' voltage
  double vBX;       // b-c'  voltage
  double vCS;       // c'-s  voltage

  double vBE_old;   // b'-e' voltage, from previous newton step.
  double vBC_old;   // b'-c' voltage, from previous newton step.

  double vBE_orig;  // b'-e' voltage, before pinning.
  double vBC_orig;  // b'-c' voltage, before pinning.

  double qB;        // Base charge factor
  double invqB;     // inverse of qB
  double dqBdvEp;   // d(qB)/d(vE')
  double dqBdvBp;   // d(qB)/d(vB')
  double dqBdvCp;   // d(qB)/d(vC')

  double iBE;       // b-e current (not including capacitors)
  double iBC;       // b-c current (not including capacitors)
  double iBEleak;
  double iBCleak;
  double iCE;       // c-e current

  double iB;        // total current to base
  double iC;        // total current to collector
  double iE;        // total current to emitter

  // high current versions of iBE, gBE.
  double iBEhighCurr;
  double gBEhighCurr;

  double gBE;
  double gBC;
  double gBEleak;
  double gBCleak;

  double gEpr;      // conductance for e-e' resistance
  double gCpr;      // conductance for c-c' resistance

  double gX;        // conductance for b-b' resistance

  double geqCB;     // high current forward transit effect (conductance)
  double capeqCB;   // high current forward transit effect (capacitance)

  // Partial derivatives that originally were local to load function
  double diBrdvB;
  double diBrdvEp;
  double diBrdvCp;
  double diBrdvBp;
  double diCEdvEp;
  double diCEdvCp;
  double diCEdvBp;
  double diBEdvEp;
  double diBEdvCp;
  double diBEdvBp;

  double gBEtot;
  double gBCtot;

  //state variables
  double qBEdiff;   // charge in the b-e diffusion capacitor
  double iBEdiff;   // current through the b-e diffusion capacitor
  double capBEdiff; // capacitance for b-e diffusion capacitor
  double qBEdep;    // charge in the b-e depletion capacitor
  double iBEdep;    // current through the b-e diffusion capacitor
  double capBEdep; // capacitance for b-e depletion capacitor
  double qCS;       // charge in the c-s capacitor
  double iCS;       // current through the c-s capacitor
  double capCS;     // capacitance for c-s capacitor
  double qBCdiff;   // charge in the b-c diffusion capacitor
  double iBCdiff;   // current through the b-c diffusion capacitor
  double capBCdiff; // capacitance for b-c diffusion capacitor
  double qBCdep;    // charge in the b-c depletion capacitor
  double iBCdep;    // current through the b-c depletion capacitor
  double capBCdep;  // capacitance for b-c depletion capacitor
  double qBX;       // charge in the b-cP capacitor
  double iBX;       // current through the b-cP capacitor
  double capBX;     // capacitance for b-cP capacitor

  //local indexing of solution and state variables
  int li_Coll;
  int li_CollP;
  int li_Base;
  int li_BaseP;
  int li_Emit;
  int li_EmitP;
  int li_Subst;

  //new variables for full new DAE integration of excess phase term
  int li_Ifx;
  int li_dIfx;

  int li_qstateBEdiff;
  int li_qstateBEdep;
  int li_qstateCS;
  int li_qstateBCdiff;
  int li_qstateBCdep;
  int li_qstateBX;

  // for the "old" excess phase calculation.
  int li_istateCEXBC;

  // stored data for limiting 
  int li_storevBE;
  int li_storevBC;
  int li_store_capeqCB;
  
  // branch data vector indexes for lead currents and power
  /* int li_branchvBE; */
  /* int li_branchvBC; */
  /* int li_branch_capeqCB; */
  int li_branch_dev_ib;
  int li_branch_dev_ie;
  int li_branch_dev_ic;
  int li_branch_dev_is;

  //conductance values
  double gcpr;
  double gepr;
  double gx;
  double gm;
  double go;
  double gmu;
  double gpi;
  double gccs;
  double geqbx;
  double geqbc;

  // excess phase variables
  double nextCexbc;
  double currCexbc;
  double lastCexbc;
  double phaseScalar;
  double dt0,dt1;

  // Offset variables corresponding to the above declared indices.
  int AEmitEquEmitPNodeOffset;
  int AEmitPEquEmitNodeOffset;
  int ABaseEquBasePNodeOffset;
  int ABasePEquBaseNodeOffset;
  int ACollEquCollPNodeOffset;
  int ACollPEquCollNodeOffset;
  int AEmitEquEmitNodeOffset;
  int AEmitPEquEmitPNodeOffset;
  int ABaseEquBaseNodeOffset;
  int ABasePEquBasePNodeOffset;
  int ACollEquCollNodeOffset;
  int ACollPEquCollPNodeOffset;
  int AEmitPEquBasePNodeOffset;
  int ABasePEquEmitPNodeOffset;
  int AEmitPEquCollPNodeOffset;
  int ACollPEquEmitPNodeOffset;
  int ABasePEquCollPNodeOffset;
  int ACollPEquBasePNodeOffset;
  int ABaseEquCollPNodeOffset;
  int ACollPEquBaseNodeOffset;
  int ASubstEquSubstNodeOffset;
  int ASubstEquCollPNodeOffset;
  int ACollPEquSubstNodeOffset;
  int ABaseEquEmitPNodeOffset;

  //new offsets for full new DAE integration of excess phase term
  int ACollPEquIfxNodeOffset;
  int AEmitPEquIfxNodeOffset;

  // ERK.  These 3 are only needed for dcop.
  int AIfxEquCollPNodeOffset;
  int AIfxEquBasePNodeOffset;
  int AIfxEquEmitPNodeOffset;

  int AIfxEquIfxNodeOffset;
  int AIfxEqudIfxNodeOffset;

  int AdIfxEquCollPNodeOffset;
  int AdIfxEquBasePNodeOffset;
  int AdIfxEquEmitPNodeOffset;
  int AdIfxEquIfxNodeOffset;
  int AdIfxEqudIfxNodeOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // f-matrix pointers:
  double * f_EmitEquEmitPNodePtr;
  double * f_EmitPEquEmitNodePtr;
  double * f_BaseEquBasePNodePtr;
  double * f_BasePEquBaseNodePtr;
  double * f_CollEquCollPNodePtr;
  double * f_CollPEquCollNodePtr;
  double * f_EmitEquEmitNodePtr;
  double * f_EmitPEquEmitPNodePtr;
  double * f_BaseEquBaseNodePtr;
  double * f_BasePEquBasePNodePtr;
  double * f_CollEquCollNodePtr;
  double * f_CollPEquCollPNodePtr;
  double * f_EmitPEquBasePNodePtr;
  double * f_BasePEquEmitPNodePtr;
  double * f_EmitPEquCollPNodePtr;
  double * f_CollPEquEmitPNodePtr;
  double * f_BasePEquCollPNodePtr;
  double * f_CollPEquBasePNodePtr;
  double * f_BaseEquCollPNodePtr;
  double * f_CollPEquBaseNodePtr;
  double * f_SubstEquSubstNodePtr;
  double * f_SubstEquCollPNodePtr;
  double * f_CollPEquSubstNodePtr;
  double * f_BaseEquEmitPNodePtr;

  //new offsets for full new DAE integration of excess phase term
  double * f_CollPEquIfxNodePtr;
  double * f_EmitPEquIfxNodePtr;

  // ERK.  These 3 are only needed for dcop.
  double * f_IfxEquCollPNodePtr;
  double * f_IfxEquBasePNodePtr;
  double * f_IfxEquEmitPNodePtr;

  double * f_IfxEquIfxNodePtr;
  double * f_IfxEqudIfxNodePtr;

  double * f_dIfxEquCollPNodePtr;
  double * f_dIfxEquBasePNodePtr;
  double * f_dIfxEquEmitPNodePtr;
  double * f_dIfxEquIfxNodePtr;
  double * f_dIfxEqudIfxNodePtr;


  // q-matrix pointers:
  double * q_EmitEquEmitPNodePtr;
  double * q_EmitPEquEmitNodePtr;
  double * q_BaseEquBasePNodePtr;
  double * q_BasePEquBaseNodePtr;
  double * q_CollEquCollPNodePtr;
  double * q_CollPEquCollNodePtr;
  double * q_EmitEquEmitNodePtr;
  double * q_EmitPEquEmitPNodePtr;
  double * q_BaseEquBaseNodePtr;
  double * q_BasePEquBasePNodePtr;
  double * q_CollEquCollNodePtr;
  double * q_CollPEquCollPNodePtr;
  double * q_EmitPEquBasePNodePtr;
  double * q_BasePEquEmitPNodePtr;
  double * q_EmitPEquCollPNodePtr;
  double * q_CollPEquEmitPNodePtr;
  double * q_BasePEquCollPNodePtr;
  double * q_CollPEquBasePNodePtr;
  double * q_BaseEquCollPNodePtr;
  double * q_CollPEquBaseNodePtr;
  double * q_SubstEquSubstNodePtr;
  double * q_SubstEquCollPNodePtr;
  double * q_CollPEquSubstNodePtr;
  double * q_BaseEquEmitPNodePtr;

  //new offsets for full new DAE integration of excess phase term
  double * q_CollPEquIfxNodePtr;
  double * q_EmitPEquIfxNodePtr;

  // ERK.  These 3 are only needed for dcop.
  double * q_IfxEquCollPNodePtr;
  double * q_IfxEquBasePNodePtr;
  double * q_IfxEquEmitPNodePtr;

  double * q_IfxEquIfxNodePtr;
  double * q_IfxEqudIfxNodePtr;

  double * q_dIfxEquCollPNodePtr;
  double * q_dIfxEquBasePNodePtr;
  double * q_dIfxEquEmitPNodePtr;
  double * q_dIfxEquIfxNodePtr;
  double * q_dIfxEqudIfxNodePtr;
#endif


  static std::vector< std::vector<int> > jacStamp_RB_RC_RE_;
  static std::vector< std::vector<int> > jacStamp_RB_RC_;
  static std::vector< std::vector<int> > jacStamp_RB_RE_;
  static std::vector< std::vector<int> > jacStamp_RC_RE_;
  static std::vector< std::vector<int> > jacStamp_RB_;
  static std::vector< std::vector<int> > jacStamp_RC_;
  static std::vector< std::vector<int> > jacStamp_RE_;
  static std::vector< std::vector<int> > jacStamp_;

  static std::vector<int> jacMap_RB_RC_RE_;
  static std::vector<int> jacMap_RB_RC_;
  static std::vector<int> jacMap_RB_RE_;
  static std::vector<int> jacMap_RC_RE_;
  static std::vector<int> jacMap_RB_;
  static std::vector<int> jacMap_RC_;
  static std::vector<int> jacMap_RE_;
  static std::vector<int> jacMap_;

  static std::vector< std::vector<int> > jacMap2_RB_RC_RE_;
  static std::vector< std::vector<int> > jacMap2_RB_RC_;
  static std::vector< std::vector<int> > jacMap2_RB_RE_;
  static std::vector< std::vector<int> > jacMap2_RC_RE_;
  static std::vector< std::vector<int> > jacMap2_RB_;
  static std::vector< std::vector<int> > jacMap2_RC_;
  static std::vector< std::vector<int> > jacMap2_RE_;
  static std::vector< std::vector<int> > jacMap2_;

  int callsOutputPlot;
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
  friend class Master;
  friend class bjtInstanceSensitivity;
  friend class bjtModelSensitivity;

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
  void updateIntermediateParams();

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

  //external model params
  int TYPE;                  //+1 = NPN, -1 = PNP
  double TNOM;               //nominal temperature
  double satCur;             //Saturation Current (IS)

  double betaF;              //forward beta (BF)
  bool   BFgiven;            // given flag for forward beta (spice/Pspice)
  bool   BFMgiven;           // given flag for forward beta (hspice)
  double emissionCoeffF;     //forward current emmission coeff (NF)
  double earlyVoltF;         //forward early voltage (VAF)
  bool   VAgiven;            // given flag for fwd early voltage (Hspice/Pspice)
  bool   VAFgiven;           // given flag for fwd early voltage (Spice)
  bool   VBFgiven;           // given flag for fwd early voltage (Hspice)

  double rollOffF;           //forward high current roll-off (IKF)
  bool IKFgiven;             // given flag for high current roll-off (spice)
  bool IKgiven;              // given flag for high current roll-off (pspice/hspice)
  bool JBFgiven;             // given flag for high current roll-off (hspice)

  double leakBECurrent;      //BE leakage saturation current (ISE)
  double leakBEEmissionCoeff;//BE leakage emission coeff. (NE)
  bool NEgiven;              // given flag for leakage emission coef.(NE) (spice/pspice)
  bool NLEgiven;             // given flag for leakage emission coef.(NLE) (hspice)

  double betaR;              // reverse beta (BR)
  bool BRgiven;              // given flag for reverse beta, BR (spice/pspice)
  bool BRMgiven;             // given flag for reverse beta, BR (hspice)

  double emissionCoeffR;     // reverse current emmission coeff (NR)
  double earlyVoltR;         // reverse early voltage (VAR/VB/VRB/BV)
  bool VARgiven;             // given flag for reverse early voltage(VAR)  (spice)
  bool VBgiven;              // given flag for reverse early voltage(VB)  (pspice/hspice)
  bool VRBgiven;             // given flag for reverse early voltage(VRB)  (hspice)
  bool BVgiven;              // given flag for reverse early voltage(BV)  (hspice)


  double rollOffR;           //reverse high current roll-off (IKR)
  bool IKRgiven;             // given flag  for reverse high current roll-off (IKR) (spice)
  bool JBRgiven;             // given flag  for reverse high current roll-off (JBR) (Hspice)
  double leakBCCurrent;      //BC leakage saturation current (ISC)
  double leakBCEmissionCoeff;//BC leakage emission coeff. (NC)

  double baseResist;         //zero bias base resistance (RB)
  double baseCurrHalfResist; //current for 1/2 base resistance (IRB)
  bool IRBgiven;             // given flag for 1/2 base resist. (IRB, spice)
  bool JRBgiven;             // given flag for 1/2 base resist. (JRB, spice)
  bool IOBgiven;             // given flag for 1/2 base resist. (IOB, spice)

  double minBaseResist;      //min base resistance for high current (RBM)
  double emitterResist;      //emitter resistance (RE)
  double collectorResist;    //collector resistance (RC)

  double depCapBE;           //BE zero bias depletion capacitance (CJE)
  double potBE;              //BE built-in potential (VJE)
  bool VJEgiven;             // given flag for BE built-in potential (VJE,spice)
  bool PEgiven;              // given flag for BE built-in potential (PE,pspice, hspice)

  double juncExpBE;          //BE junction exponential factor (MJE)
  bool MJEgiven;             // given flag for BE exponential factor (MJE,spice)
  bool MEgiven;              // given flag for BE exponential factor (ME,pspice, hspice)

  double transTimeF;         //ideal forward transit time (TF)
  double transTimeBiasCoeffF;//bias dependent coefficient for TF (XTF)
  double transTimeFVBC;      //VBC dependence for TF (VTF)
  double transTimeHighCurrF; //high current parameter for TF (ITF)
  bool ITFgiven;             // given flag for ITF (spice)
  bool JTFgiven;             // given flag for JTF (hspice)
  bool VTFgiven;             // given flag for VTF
  double excessPhase;        //excess phase at freq=1.0/(TF*2PI) Hz (PTF)

  double depCapBC;           //BC zero bias depletion capacitance (CJC)
  double potBC;              //BC built-in potential (VJC)
  bool VJCgiven;             // given flag for BC built-in potential (VJE,spice)
  bool PCgiven;              // given flag for BC built-in potential (PE,pspice, hspice)

  double juncExpBC;          //BC junction exponential factor (MJC)
  bool MJCgiven;             // given flag for BC exponential factor (MJC,spice)
  bool MCgiven;              // given flag for BC exponential factor (MC,pspice, hspice)

  double baseFracBCCap;      //fraction of BC cap. to int. base node (XCJC)
  bool XCJCgiven;            // given flag for spice/pspice.
  bool CDISgiven;            // given flag for hspice

  double transTimeR;         //ideal reverse transit time (TR)

  double CJS;                //zero-bias coll-subst capacitance (CJS)
  bool CJSgiven;             //zero-bias coll-subst capacitance (CJS) given flag (spice)
  bool CCSgiven;             //zero-bias coll-subst capacitance (CCS) given flag (Pspice/Hspice)
  bool CSUBgiven;            //zero-bias coll-subst capacitance (CSUB) given flag (Hspice)

  double potSubst;           //substrate junction built-in potential (VJS)
  bool VJSgiven;             // given flag for spice format (VJS)
  bool PSgiven;              // given flag for spice format (PS)
  bool PSUBgiven;            // given flag for spice format (PSUB)

  double expSubst;           //subst. junction exponential factor (MJS)
  bool MJSgiven;             // given flag for spice format (MJS)
  bool MSgiven;              // given flag for spice format (MS)
  bool ESUBgiven;            // given flag for spice format (ESUB)

  double betaExp;            //beta temperature exponent (XTB)
  bool XTBgiven;             // given for spice/pspice
  bool TBgiven;              // given for hspice
  bool TCBgiven;             // given for hspice

  double energyGap;          //energy gap for temp. effect on IS (EG)
  double tempExpIS;          //temp. exponent for IS (XTI)
  bool XTIgiven;             // given for spice/hspice
  bool PTgiven;              // given for pspice

  double depCapCoeff;        //coeff. for fwd bias depletion cap. (FC)
  double fNCoeff;            //flicker-noise coeff. (KF)
  double fNExp;              //flicker-noise exponent (AF)

  double rollOffExp;         //Pspice high-current rolloff parameter (NK)
  
  bool FCgiven;              // given flag for FC
  bool NKgiven;              // spice/pspice
  bool NKFgiven;             // hspice

  double c2;
  double c4;

  bool leakBECurrentGiven; // specified as ISE (spice)
  bool JLEgiven; // hspice version of leakBECurrent
  bool leakBCCurrentGiven; // specified as ISC (spice)
  bool JLCgiven; // hspice version of leakBCCurrent

  bool c2Given;
  bool c4Given;
  bool minBaseResistGiven;

  //generated model params
  double invEarlyVoltF;      //inverse of fwd early voltage
  double invEarlyVoltR;      //inverse of rvs early voltage
  double invRollOffF;        //inverse of fwd high curr. roll-off
  double invRollOffR;        //inverse of rvs high curr. roll-off
  double collectorConduct;   //collector conductance
  double emitterConduct;     //emitter conductance
  double transTimeVBCFac;    //VBC transit time factor
  double excessPhaseFac;     //excess phase factor

  double f2;
  double f3;
  double f6;
  double f7;
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

  friend class Instance;
  friend struct Traits;
  friend class Model;
};

void registerDevice(const DeviceCountMap& deviceMap = DeviceCountMap(),
                    const std::set<int>& levelSet = std::set<int>());

} // namespace BJT
} // namespace Device
} // namespace Xyce

#endif

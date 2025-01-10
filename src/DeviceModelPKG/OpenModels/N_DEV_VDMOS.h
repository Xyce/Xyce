//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Vertical, double-dffused power MOSFET
//                  This device model is based on the Uniform Charge Control
//                  Model of Fjeldly, et. al.
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

#ifndef Xyce_N_DEV_VDMOS_h
#define Xyce_N_DEV_VDMOS_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>

#include <N_DEV_MOSFET1.h>

namespace Xyce {
namespace Device {
namespace VDMOS {

// ---------- Forward Declarations -------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, MOSFET1::Traits>
{
  static const char *name() {return  "Power MOSFET";}
  static const char *deviceTypeName() {return "M level 18";}
  static int numNodes() {return 4;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : pmc
// Creation Date : 1/16/2004
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
     Model &                   Miter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  void registerStateLIDs( const std::vector<int> & staLIDVecRef);
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef);
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers();

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool applyScale ();
  bool processParams ();
  bool updateTemperature(const double & temp_tmp);

  bool UCCMqmeyer(double vgps, double vgpdd, double vgpb, double von_local,
                  double vddsat_local, double & capgs_local, double & capgdd_local,
                  double & capgb_local, double phi, double cox );

  bool UCCMMeyercap(double vgps, double vgpdd, double vgpb,
                    double & cgs, double & cgd, double & cgb);
  bool UCCMCharges(double vgps, double vgpdd, double vgpb,
                   double & qD, double & qS, double & qB);
  bool UCCMcvon(double vbs_local, double & von_local, double & dvonvbs_local);
  bool UCCMmosa1(double vgps, double vdds, double dvonvbs,
                 double & cdraindrift_loc, double & vsate);
  bool UCCMmosa2(double vgps, double vdds, double dvonvbs,
                 double & cdraindrift_loc, double & vsate);

  inline bool isConverged();

public:
  // iterator reference to the model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  int states;     // index into state table for this device
  int dNode;      // number of the gate node of the mosfet
  int gNode;      // number of the gate node of the mosfet
  int sNode;      // number of the source node of the mosfet
  int bNode;      // number of the bulk node of the mosfet
  int dNodePrime; // number of the internal drain node of the mosfet
  int gNodePrime; // number of the internal gate node of the mosfet
  int sNodePrime; // number of the internal source node of the mosfet
  int dDriftNode; // number of the internal draindrift node of the mosfet

  double l;                   // the length of the channel region
  double w;                   // the width of the channel region
  double drainArea;           // the area of the drain diffusion
  double sourceArea;          // the area of the source diffusion
  double drainSquares;        // the length of the drain in squares
  double sourceSquares;       // the length of the source in squares
  double drainPerimeter;
  double sourcePerimeter;
  double sourceCond;          //conductance of source(or 0):set in setup
  double gateCond;            //conductance of gate(or 0):set in setup
  double drainCond;           //conductance of drain(or 0):set in setup
  double draindriftCond;      //conductance of draindrift channel
  double numberParallel;      // number simulated parallel mosfets
  double vt;                  // CONSTKoverQ*temp: set in updateTemperature

  double temp;                // operating temperature of this instance
  double dtemp;               // delta temperature of this instance
  bool dtempGiven;            // delta temperature given
  double tSurfMob;            // temperature corrected surface mobility
  double tPhi;                // temperature corrected Phi
  double tVto;                // temperature corrected Vto
  double tSatCur;             // temperature corrected saturation Cur.
  double tSatCurDens;         // temperature corrected sat. cur. density
  double tCbd;                // temperature corrected B-D Capacitance
  double tCbs;                // temperature corrected B-S Capacitance
  double tCj;                 // temperature corrected Bulk bottom Cap
  double tCjsw;               // temperature corrected Bulk side Cap
  double tBulkPot;            // temperature corrected Bulk potential
  double tDepCap;             // temperature adjusted transition point in
  // the cureve matching Fc * Vj
  double tVbi;                // temperature adjusted Vbi

  double von;
  double vdsat;
  double vddsat;
  double sourceVcrit;         // vcrit for pos. vds
  double drainVcrit;          // vcrit for neg. vds
  double draindriftVcrit;     // vcrit for neg. vdd
  double cdd;                 // draindrift current
  double cd;                  // drain current
  double gmbs;                // bulk-source transconductance
  double gm;                  // transconductance
  double gddd;                // draindrift conductance
  double dIdd_dVd;            // draindrift derivative
  double gds;                 // drain-source conductance
  double gdds;                // drain-source conductance
  double gdsshunt;            // drain-source shunt conductance
  double gbs;                 // bulk-source conductance
  double gbd;                 // bulk-drain conductance
  double cbd;
  double Cbd;
  double Cbdsw;
  double cbs;
  double Cbs;
  double Cbssw;
  double f2d;
  double f3d;
  double f4d;
  double f2s;
  double f3s;
  double f4s;
  double n0;
  double vp;
  double gammas;
  double gammal;
  double gchi0;
  double vtoo;
  double vthLimit;

  int mode;                   // device mode : 1 = normal, -1 = inverse
  double mode_low;
  double mode_high;
  bool off ;                  // non-zero to indicate device is off for
  // dc analysis
  bool dNodePrimeSet;
  bool sNodePrimeSet;

  bool limitedFlag;           // for convergence testing


  ////////////////////////////////////////////////////////////////////
  // these are intermediate variables added to the instance class instead
  // of leaving them to be calculated repeatedly in the load function

  // some caluclated quantities
  double EffectiveLength;
  double DrainSatCur;
  double SourceSatCur;
  double GateSourceOverlapCap;
  double GateDrainOverlapCap;
  double GateBulkOverlapCap;
  double OxideCap;

  // Solution variables and intermediate quantities
  // drain,source,gate, bulk, drainprime and sourceprime voltages
  double Vd;
  double Vs;
  double Vg;
  double Vb;
  double Vdp;
  double Vgp;
  double Vsp;
  double Vdd;
  double Vd1p;
  // voltage drops between pairs of nodes
  double Vddp;   // drain-drain'
  double Vddd;   // drain-draindrift
  double Vdddp;  // draindrift-drain'
  double Vssp;   // source-source'
  double Vbsp;   // bulk-source'
  double Vbdp;   // bulk-drain'
  double Vggp;   // gate-gate'
  double Vgpsp;  // gate'-source'
  double Vgpdp;  // gate'-drain'
  double Vgpb;   // gate'-bulk
  double Vdpsp;  // drop across channel
  double Vbdd;   // bulk-draindrift
  double D1vd;   // D1prime-drain voltage drop

  // the gate-drain voltage drop isn't actually a state variable, but it
  // is calculated at the same time and in the same manner as the state
  // vars.  So here we go, sticking it in the instance class.
  double vgpd;

  // the variables capgs, capgd and capgb are the raw output of
  // qmeyer.  They get massaged into total capacitances in
  // updateIntermediateVars, and get used in updatePrimaryState to get
  // charges on the capacitors.

  double Capgs;   // total gate-source capacitance
  double Capgdd;  // total gate-drain capacitance
  double Capgb;   // total gate-bulk capacitance

  // current through source and drain resistors
  double Isource;
  double Igate;
  double Idrain;
  double Idraindrift;
  double Irdsshunt;
  double Ird1rs;

  // ISUBMOD calculation
  double mm1;
  double dmm1vgs;
  double dmm1vds;
  double dmm1vbs;
  double ISUB;
  double GMSUB;
  double GDDSSUB;
  double GBSSUB;

  double cdrain;       // the channel current
  double cdraindrift;  // the draindrift current

  // these are calculated in loadRHS and used in the jacobian load
  double Gm,Gmbs;  // we do this so we don't really need the xnrm/xrev vars
  double revsum;   // described in comments at the end of
  double nrmsum;   // updateIntermediateVars (uIVB in remaining comments)
  double cdreq;

  // variables related to diode 1, which connects the source to the drain
  int    D1DIOoff1;           // 'off' flag for diode
  double D1DIOarea;           // area factor for the diode
  double D1DIOinitCond;       // initial condition
  double D1DIOtemp;           // temperature of the instance
  double D1DIOtJctPot;        // temperature adjusted junction potential
  double D1DIOtJctCap;        // temperature adjusted junction capacitance
  double D1DIOtDepCap;        // temperature adjusted transition point in
  //              the curve matching (Fc * Vj )
  double D1DIOtSatCur;        // temperature adjusted saturation current
  double D1DIOtSatRCur;
  double D1DIOtVcrit;         // temperature adjusted V crit
  double D1DIOtF1;            // temperature adjusted f1
  double D1DIOtBrkdwnV;
  double D1gspr;              // area-scaled conductance
  double D1gd;
  double D1cdeq;
  double D1vt;      // K t / Q
  double D1vte;

  // end of intermediate variables that aren't state variables
  ////////////////////////////
  int li_Drain;
  int li_DrainPrime;
  int li_Source;
  int li_SourcePrime;
  int li_Gate;
  int li_GatePrime;
  int li_Bulk;
  int li_DrainDrift;
  int li_D1Prime;

  //  The VDMOS has a 9x9 matrix block;  35 entries are nonzero:
  // --------------------------------------------------------------------------+
  // | #NZ     |       |                                                       |
  // | entries |       |  V_d   V_g   V_s   V_b   V_d'  V_g'  V_s'  V_dd  V_d1'|
  // --------------------------------------------------------------------------+
  // |    4    | KCL_d |  aa          -gs                           -ab    A   |
  // |    2    | KCL_g |         a                      -a                     |
  // |    4    | KCL_s |  -gs          g                       h           B   |
  // |    4    | KCL_b |                     j     k     i     l               |
  // |    5    | KCL_d'|                     o     p     n     q     b         |
  // |    5    | KCL_g'|        -a           d     e     c     f               |
  // |    5    | KCL_s'|               s     t     u     r     v               |
  // |    3    | KCL_dd|  -ab                      b                 w         |
  // |    3    |KCL_D1'|   C           D                                   E   |
  // --------------------------------------------------------------------------+
  //     35 total
  ////////////////////////////////////////////////////////////////////
  // Offset variables corresponding to the above declared indices.

  // Jacobian Matrix Offset:

  // V_d Row:
  int ADrainEquDrainNodeOffset;             // aa
  int ADrainEquSourceNodeOffset;            // -gs
  int ADrainEquDrainDriftNodeOffset;        // -ab
  int ADrainEquD1PrimeNodeOffset;           // A

  // V_g Row:
  int AGateEquGateNodeOffset;               // a
  int AGateEquGatePrimeNodeOffset;          // -a

  // V_s Row:
  int ASourceEquDrainNodeOffset;            // -gs
  int ASourceEquSourceNodeOffset;           // g
  int ASourceEquSourcePrimeNodeOffset;      // h
  int ASourceEquD1PrimeNodeOffset;          // B

  // V_b Row:
  int ABulkEquBulkNodeOffset;               // j
  int ABulkEquDrainPrimeNodeOffset;         // k
  int ABulkEquGatePrimeNodeOffset;          // i
  int ABulkEquSourcePrimeNodeOffset;        // l

  // V_d' Row:
  int ADrainPrimeEquBulkNodeOffset;         // o
  int ADrainPrimeEquDrainPrimeNodeOffset;   // p
  int ADrainPrimeEquGatePrimeNodeOffset;    // n
  int ADrainPrimeEquSourcePrimeNodeOffset;  // q
  int ADrainPrimeEquDrainDriftNodeOffset;   // b

  // V_g' Row:
  int AGatePrimeEquGateNodeOffset;          // -a
  int AGatePrimeEquBulkNodeOffset;          // d
  int AGatePrimeEquDrainPrimeNodeOffset;    // e
  int AGatePrimeEquGatePrimeNodeOffset;     // c
  int AGatePrimeEquSourcePrimeNodeOffset;   // f

  // V_s' Row:
  int ASourcePrimeEquSourceNodeOffset;      // s
  int ASourcePrimeEquBulkNodeOffset;        // t
  int ASourcePrimeEquDrainPrimeNodeOffset;  // u
  int ASourcePrimeEquGatePrimeNodeOffset;   // r
  int ASourcePrimeEquSourcePrimeNodeOffset; // v

  // V_dd Row:
  int ADrainDriftEquDrainNodeOffset;        // -ab
  int ADrainDriftEquDrainPrimeNodeOffset;   // b
  int ADrainDriftEquDrainDriftNodeOffset;   // w


  // V_d1' row:
  int AD1PrimeEquDrainNodeOffset;          // C
  int AD1PrimeEquSourceNodeOffset;         // D
  int AD1PrimeEquD1PrimeNodeOffset;        // E

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ////////////////////////////////////////////////////////////////////

  // F-matrix pointers:
  double * f_DrainEquDrainNodePtr;
  double * f_DrainEquSourceNodePtr;
  double * f_DrainEquDrainDriftNodePtr;
  double * f_DrainEquD1PrimeNodePtr;

  double * f_GateEquGateNodePtr;
  double * f_GateEquGatePrimeNodePtr;

  double * f_SourceEquDrainNodePtr;
  double * f_SourceEquSourceNodePtr;
  double * f_SourceEquSourcePrimeNodePtr;
  double * f_SourceEquD1PrimeNodePtr;

  double * f_BulkEquBulkNodePtr;
  double * f_BulkEquDrainPrimeNodePtr;
  double * f_BulkEquGatePrimeNodePtr;
  double * f_BulkEquSourcePrimeNodePtr;

  double * f_DrainPrimeEquBulkNodePtr;
  double * f_DrainPrimeEquDrainPrimeNodePtr;
  double * f_DrainPrimeEquGatePrimeNodePtr;
  double * f_DrainPrimeEquSourcePrimeNodePtr;
  double * f_DrainPrimeEquDrainDriftNodePtr;

  double * f_GatePrimeEquGateNodePtr;
  double * f_GatePrimeEquBulkNodePtr;
  double * f_GatePrimeEquDrainPrimeNodePtr;
  double * f_GatePrimeEquGatePrimeNodePtr;
  double * f_GatePrimeEquSourcePrimeNodePtr;

  double * f_SourcePrimeEquSourceNodePtr;
  double * f_SourcePrimeEquBulkNodePtr;
  double * f_SourcePrimeEquDrainPrimeNodePtr;
  double * f_SourcePrimeEquGatePrimeNodePtr;
  double * f_SourcePrimeEquSourcePrimeNodePtr;

  double * f_DrainDriftEquDrainNodePtr;
  double * f_DrainDriftEquDrainPrimeNodePtr;
  double * f_DrainDriftEquDrainDriftNodePtr;

  double * f_D1PrimeEquDrainNodePtr;
  double * f_D1PrimeEquSourceNodePtr;
  double * f_D1PrimeEquD1PrimeNodePtr;

  // Q-matrix pointers:
  double * q_DrainEquDrainNodePtr;
  double * q_DrainEquSourceNodePtr;
  double * q_DrainEquDrainDriftNodePtr;
  double * q_DrainEquD1PrimeNodePtr;

  double * q_GateEquGateNodePtr;
  double * q_GateEquGatePrimeNodePtr;

  double * q_SourceEquDrainNodePtr;
  double * q_SourceEquSourceNodePtr;
  double * q_SourceEquSourcePrimeNodePtr;
  double * q_SourceEquD1PrimeNodePtr;

  double * q_BulkEquBulkNodePtr;
  double * q_BulkEquDrainPrimeNodePtr;
  double * q_BulkEquGatePrimeNodePtr;
  double * q_BulkEquSourcePrimeNodePtr;

  double * q_DrainPrimeEquBulkNodePtr;
  double * q_DrainPrimeEquDrainPrimeNodePtr;
  double * q_DrainPrimeEquGatePrimeNodePtr;
  double * q_DrainPrimeEquSourcePrimeNodePtr;
  double * q_DrainPrimeEquDrainDriftNodePtr;

  double * q_GatePrimeEquGateNodePtr;
  double * q_GatePrimeEquBulkNodePtr;
  double * q_GatePrimeEquDrainPrimeNodePtr;
  double * q_GatePrimeEquGatePrimeNodePtr;
  double * q_GatePrimeEquSourcePrimeNodePtr;

  double * q_SourcePrimeEquSourceNodePtr;
  double * q_SourcePrimeEquBulkNodePtr;
  double * q_SourcePrimeEquDrainPrimeNodePtr;
  double * q_SourcePrimeEquGatePrimeNodePtr;
  double * q_SourcePrimeEquSourcePrimeNodePtr;

  double * q_DrainDriftEquDrainNodePtr;
  double * q_DrainDriftEquDrainPrimeNodePtr;
  double * q_DrainDriftEquDrainDriftNodePtr;

  double * q_D1PrimeEquDrainNodePtr;
  double * q_D1PrimeEquSourceNodePtr;
  double * q_D1PrimeEquD1PrimeNodePtr;
#endif

  ////////////////////////////////////////////////////////////////////
  // 3f5 State Variables & related quantities:
  // voltage drops
  double vbdd;
  double vbs;
  double vgpdd;
  double vgps;
  double vdds;

  // "original" versions of various voltage drop variables:
  double vbdd_orig;
  double vbs_orig;
  double vgpdd_orig;
  double vgps_orig;
  double vdds_orig;
  double D1vd_orig;

  // "old" versions of various voltage drop variables:
  double vbdd_old;
  double vbs_old;
  double vgpdd_old;
  double vgps_old;
  double vdds_old;
  double D1vd_old;

  // meyer capacitances
  //gate-source capacitor
  double capgs;   //value
  double qgs;     // charge
  double cqgs;    // current
  // gate-drain capacitor
  double capgdd;  //value
  double qgdd;    //charge
  double cqgdd;   // current
  //gate-bulk capacitor
  double capgb;   //value
  double qgb;     //charge
  double cqgb;    //current

  // diode capacitances
  double capbd;
  double qbd;     // bulk-drain capacitor charge
  double cqbd;    // bulk-drain capacitor current

  double capbs;
  double qbs;     // bulk-source capacitor charge
  double cqbs;    // bulk-source capacitor current

  // diode #1
  double D1DIOcapCharge;
  double D1DIOcapCurrent;
  double D1capd;

  // indices into the state vector.
  int li_state_vbdd;
  int li_state_vbs;
  int li_state_vgps;
  int li_state_vdds;
  int li_state_D1vd;

  int li_state_capgs;
  int li_state_capgdd;
  int li_state_capgb;

  int li_state_qgs;
  int li_state_qgdd;
  int li_state_qgb;

  int li_state_qbd;
  int li_state_qbs;

  int li_state_D1DIOcapCharge;

  int li_state_von;

  // store vector space for lead currents
  int li_branch_data_d;
  int li_branch_data_g;
  int li_branch_data_s;
  int li_branch_data_b;

  static std::vector< std::vector<int> > jacStamp_DC_SC_GC;
  static std::vector< std::vector<int> > jacStamp_DC_GC;
  static std::vector< std::vector<int> > jacStamp_SC_GC;
  static std::vector< std::vector<int> > jacStamp_DC_SC;
  static std::vector< std::vector<int> > jacStamp_GC;
  static std::vector< std::vector<int> > jacStamp_SC;
  static std::vector< std::vector<int> > jacStamp_DC;
  static std::vector< std::vector<int> > jacStamp;

  static std::vector<int> jacMap_DC_SC_GC;
  static std::vector<int> jacMap_DC_GC;
  static std::vector<int> jacMap_SC_GC;
  static std::vector<int> jacMap_DC_SC;
  static std::vector<int> jacMap_GC;
  static std::vector<int> jacMap_SC;
  static std::vector<int> jacMap_DC;
  static std::vector<int> jacMap;

  static std::vector< std::vector<int> > jacMap2_DC_SC_GC;
  static std::vector< std::vector<int> > jacMap2_DC_GC;
  static std::vector< std::vector<int> > jacMap2_SC_GC;
  static std::vector< std::vector<int> > jacMap2_DC_SC;
  static std::vector< std::vector<int> > jacMap2_GC;
  static std::vector< std::vector<int> > jacMap2_SC;
  static std::vector< std::vector<int> > jacMap2_DC;
  static std::vector< std::vector<int> > jacMap2;

  // duplicating all of the above, but with the RD1RS resistor nonzero
  static std::vector< std::vector<int> > jacStamp_D1C_DC_SC_GC;
  static std::vector< std::vector<int> > jacStamp_D1C_DC_GC;
  static std::vector< std::vector<int> > jacStamp_D1C_SC_GC;
  static std::vector< std::vector<int> > jacStamp_D1C_DC_SC;
  static std::vector< std::vector<int> > jacStamp_D1C_GC;
  static std::vector< std::vector<int> > jacStamp_D1C_SC;
  static std::vector< std::vector<int> > jacStamp_D1C_DC;
  static std::vector< std::vector<int> > jacStamp_D1C;

  static std::vector<int> jacMap_D1C_DC_SC_GC;
  static std::vector<int> jacMap_D1C_DC_GC;
  static std::vector<int> jacMap_D1C_SC_GC;
  static std::vector<int> jacMap_D1C_DC_SC;
  static std::vector<int> jacMap_D1C_GC;
  static std::vector<int> jacMap_D1C_SC;
  static std::vector<int> jacMap_D1C_DC;
  static std::vector<int> jacMap_D1C;

  static std::vector< std::vector<int> > jacMap2_D1C_DC_SC_GC;
  static std::vector< std::vector<int> > jacMap2_D1C_DC_GC;
  static std::vector< std::vector<int> > jacMap2_D1C_SC_GC;
  static std::vector< std::vector<int> > jacMap2_D1C_DC_SC;
  static std::vector< std::vector<int> > jacMap2_D1C_GC;
  static std::vector< std::vector<int> > jacMap2_D1C_SC;
  static std::vector< std::vector<int> > jacMap2_D1C_DC;
  static std::vector< std::vector<int> > jacMap2_D1C;
};

//--------------------Class Model-----------------------------------

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : pmc
// Creation Date : 1/16/2004
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

  int dtype;                  // device type : 1 = nmos,  -1 = pmos
  int gateType;

  double l0;                  // the length of the channel region
  double w0;                  // the width of the channel region
  double tnom;                // temperature at which parameters measured
  double latDiff;
  double jctSatCurDensity;    // input - use tSatCurDens
  double jctSatCur;           // input - use tSatCur instead
  double drainResistance;
  double gateResistance;
  double sourceResistance;
  double sheetResistance;
  double gateSourceOverlapCapFactor;
  double gateDrainOverlapCapFactor;
  double gateBulkOverlapCapFactor;
  double oxideCapFactor;
  double vt0;                 // input - use tVto
  double capBD;               // input - use tCbs
  double capBS;               // input - use tCbd
  double timeScale;
  double bulkCapFactor;       // input - use tCj
  double sideWallCapFactor;   // input - use tCjsw
  double bulkJctPotential;    // input - use tBulkPot
  double bulkJctBotGradingCoeff;
  double bulkJctSideGradingCoeff;
  double fwdCapDepCoeff;
  double phi;                 // input - use tPhi
  double gamma;
  double lambda;
  double substrateDoping;
  double surfaceStateDensity;
  double oxideThickness;
  double surfaceMobility;    // input - use tSurfMob
  double surfaceMobility0;

  bool capBDGiven;
  bool capBSGiven;
  bool bulkCapFactorGiven;
  bool sideWallCapFactorGiven;
  bool vpGiven;

  // extra parameters needed for VDMOS
  double maxDriftVel;         // max carrier drift vel.
  double junctionDepth;       // xj
  double rdi;                 // intrinsic ohmic res. per unit width
  double rsi;                 // intrinsic ohmic res. per unit width

  double delta;
  double eta;                 // subthreshold ideality factor
  double m;
  double mc;
  double sigma0;              // threshold voltage coeff.
  double vsigmat;             // DIBL ouble rsub;
  double vsigma;              // DIBL parameter
  double theta;
  double gammas0;
  double gammal0;
  double lgammas;
  double wgammas;
  double lgammal_;
  double wgammal;
  double kacc;
  double gb;
  double knit;
  double nitd;
  double nits;
  double mm;
  double k;
  double deltaSqr;
  double kvt;
  double mdtemp;
  double kvs;
  double tvs;
  double mth;
  double artd;
  double brtd;
  double crtd;
  double drtd;
  double nrtd;
  double n2;

  // WARD-LIKE MODEL PARAMETERS
  double xqc;                 // CV_mod_2 partitioning
  double mcv;                 // CV_mod_2 fit parameter
  double vfb;                 // Flatband voltage
  int    fpe;                 // select partitioning
  double alpha;
  // END WARD-LIKE MODEL PARAMETERS
  int    cv;                  // 1= Meyer, 2= Meyer-like capacitances
  int    cve;                 // 1= Meyer, 2= Ward capacitances
  double ls;
  double rsub;
  double vp;
  double ai;
  double bi;
  double delmax;
  double md;
  int    isubmod;

  double kaccd;
  double kaccs;
  double invm;
  double invmc;
  double invmd;

  double fact1;
  double vtnom;
  double egfet1;
  double pbfact1;

  // drift region parameters
  double driftParamA;
  double driftParamB;
  double rdsshunt;

  // variables related to diode 1 connecting source to drain
  double D1DIOsatCur;               // saturation current
  double D1DIOresist;               // ohmic series resistance
  double D1DIOconductance;          // conductance corresponding to ohmic R
  double D1DIOemissionCoeff;        // emission coefficient (N)
  double D1DIOtransitTime;          // transit time (TT)
  double D1DIOjunctionCap;          // Junction Capacitance (Cj0)
  double D1DIOjunctionPot;          // Junction Potecomp
  double D1DIOgradingCoeff;         // grading coefficient (m)
  double D1DIOactivationEnergy;     // activation energy (EG)
  double D1DIOsaturationCurrentExp; // Saturation current exponential (XTI)
  double D1DIOdepletionCapCoeff;    // Depletion Cap fraction coeff (FC)
  double D1DIObreakdownVoltage;     // Voltage at reverse breakdown
  double D1DIObreakdownCurrent;     // Current at above voltage
  double D1DIOf2;                   // coeff for capacitance equation precomp
  double D1DIOf3;                   // coeff for capacitance equation precomp
  double D1DIOnomTemp;              // nominal temperature
  double D1DIOfNcoef;
  double D1DIOfNexp;
  double D1DIOikf;
  double D1DIOisr;
  double D1DIOnr;
  bool   D1DIObreakdownVoltageGiven;
};

//-----------------------------------------------------------------------------
// Function      : Instance:isConverged ()
// Purpose       : Return whether a BSIM3 device has done something that should
//                  be interpreted as invalidating other convergence tests
//                  In case of bsim3, just do it if the limiter function
//                  pnjlim.  This actually agrees with how the Check flag
//                  is used in Spice3F5 b3ld.c
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
inline bool Instance::isConverged()
{
  return (!limitedFlag);
}
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

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);

};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace VDMOS
} // namespace Device
} // namespace Xyce

#endif

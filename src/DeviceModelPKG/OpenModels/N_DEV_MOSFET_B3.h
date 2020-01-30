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

#ifndef Xyce_N_DEV_MOSFET_B3_h
#define Xyce_N_DEV_MOSFET_B3_h

#include <Sacado_No_Kokkos.hpp>

// ----------   Standard Includes   ----------
#include <N_DEV_Configuration.h>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>

#include <N_DEV_MOSFET1.h>

namespace Xyce {
namespace Device {
namespace MOSFET_B3 {

class Model;
class Instance;

typedef Sacado::Fad::SFad<double, 1> fadType;

/// general sensitivity functor for all instance params.
class bsim3InstanceSensitivity :  public baseSensitivity
{
  public:
  bsim3InstanceSensitivity() : 
    baseSensitivity() {};

  virtual ~bsim3InstanceSensitivity() {};

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
class bsim3ModelSensitivity :  public baseSensitivity
{
  public:
  bsim3ModelSensitivity() : 
    baseSensitivity() {};

  virtual ~bsim3ModelSensitivity() {};

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

static bsim3InstanceSensitivity bsim3InstanceSens;
static bsim3ModelSensitivity bsim3ModelSens;


struct Traits : public DeviceTraits<Model, Instance, MOSFET1::Traits>
{
  static const char *name() {return  "BSIM3";}
  static const char *deviceTypeName() {return "M level 9";}
  static int numNodes() {return 4;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

#if 1
//-----------------------------------------------------------------------------
// Struct        : bsim3SizeDependParam
// Purpose       : copied over from the 3f5 code.  Almost no changes.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/2018
//-----------------------------------------------------------------------------
template <typename T = double>
class SizeDependParam
{
  friend class Model;
  friend class Instance;
  friend class Master;

//private:
public:
  T Width;
  T Length;

  T cdsc;
  T cdscb;
  T cdscd;
  T cit;
  T nfactor;
  T xj;
  T vsat;
  T at;
  T a0;
  T ags;
  T a1;
  T a2;
  T keta;
  T nsub;
  T npeak;
  T ngate;
  T gamma1;
  T gamma2;
  T vbx;
  T vbi;
  T vbm;
  T vbsc;
  T xt;
  T phi;
  T litl;
  T k1;
  T kt1;
  T kt1l;
  T kt2;
  T k2;
  T k3;
  T k3b;
  T w0;
  T nlx;
  T dvt0;
  T dvt1;
  T dvt2;
  T dvt0w;
  T dvt1w;
  T dvt2w;
  T drout;
  T dsub;
  T vth0;
  T ua;
  T ua1;
  T ub;
  T ub1;
  T uc;
  T uc1;
  T u0;
  T ute;
  T voff;
  T vfb;
  T delta;
  T rdsw;
  T rds0;
  T prwg;
  T prwb;
  T prt;
  T eta0;
  T etab;
  T pclm;
  T pdibl1;
  T pdibl2;
  T pdiblb;
  T pscbe1;
  T pscbe2;
  T pvag;
  T wr;
  T dwg;
  T dwb;
  T b0;
  T b1;
  T alpha0;
  T alpha1;
  T beta0;


  // CV model
  T elm;
  T cgsl;
  T cgdl;
  T ckappa;
  T cf;
  T clc;
  T cle;
  T vfbcv;
  T noff;
  T voffcv;
  T acde;
  T moin;


  // Pre-calculated constants
  T dw;
  T dl;
  T leff;
  T weff;

  T dwc;
  T dlc;
  T leffCV;
  T weffCV;
  T abulkCVfactor;
  T cgso;
  T cgdo;
  T cgbo;
  T tconst;

  T u0temp;
  T vsattemp;
  T sqrtPhi;
  T phis3;
  T Xdep0;
  T sqrtXdep0;
  T theta0vb0;
  T thetaRout;

  T cof1;
  T cof2;
  T cof3;
  T cof4;
  T cdep0;
  T vfbzb;
  T ldeb;
  T k1ox;
  T k2ox;

  // ERK.  I added this to make temperature sweeps work.
  T referenceTemperature;
};
#else
//-----------------------------------------------------------------------------
// Struct        : bsim3SizeDependParam
// Purpose       : copied over from the 3f5 code.  Almost no changes.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/14/01
//-----------------------------------------------------------------------------
class SizeDependParam
{
  friend class Model;
  friend class Instance;
  friend class Master;

private:
  double Width;
  double Length;

  double cdsc;
  double cdscb;
  double cdscd;
  double cit;
  double nfactor;
  double xj;
  double vsat;
  double at;
  double a0;
  double ags;
  double a1;
  double a2;
  double keta;
  double nsub;
  double npeak;
  double ngate;
  double gamma1;
  double gamma2;
  double vbx;
  double vbi;
  double vbm;
  double vbsc;
  double xt;
  double phi;
  double litl;
  double k1;
  double kt1;
  double kt1l;
  double kt2;
  double k2;
  double k3;
  double k3b;
  double w0;
  double nlx;
  double dvt0;
  double dvt1;
  double dvt2;
  double dvt0w;
  double dvt1w;
  double dvt2w;
  double drout;
  double dsub;
  double vth0;
  double ua;
  double ua1;
  double ub;
  double ub1;
  double uc;
  double uc1;
  double u0;
  double ute;
  double voff;
  double vfb;
  double delta;
  double rdsw;
  double rds0;
  double prwg;
  double prwb;
  double prt;
  double eta0;
  double etab;
  double pclm;
  double pdibl1;
  double pdibl2;
  double pdiblb;
  double pscbe1;
  double pscbe2;
  double pvag;
  double wr;
  double dwg;
  double dwb;
  double b0;
  double b1;
  double alpha0;
  double alpha1;
  double beta0;


  // CV model
  double elm;
  double cgsl;
  double cgdl;
  double ckappa;
  double cf;
  double clc;
  double cle;
  double vfbcv;
  double noff;
  double voffcv;
  double acde;
  double moin;


  // Pre-calculated constants
  double dw;
  double dl;
  double leff;
  double weff;

  double dwc;
  double dlc;
  double leffCV;
  double weffCV;
  double abulkCVfactor;
  double cgso;
  double cgdo;
  double cgbo;
  double tconst;

  double u0temp;
  double vsattemp;
  double sqrtPhi;
  double phis3;
  double Xdep0;
  double sqrtXdep0;
  double theta0vb0;
  double thetaRout;

  double cof1;
  double cof2;
  double cof3;
  double cof4;
  double cdep0;
  double vfbzb;
  double ldeb;
  double k1ox;
  double k2ox;

  // ERK.  I added this to make temperature sweeps work.
  double referenceTemperature;
};
#endif

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class SizeDependParam<double>;
  //friend class SizeDependParam;
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;
  friend class bsim3InstanceSensitivity;
  friend class bsim3ModelSensitivity;

  // functions
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

  bool processParams ();

  bool updateTemperature(const double & temp_tmp);
  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  double StrongInversionNoiseEval(double Vds, double freq, double temp);
  int getNumNoiseSources () const;
  void setupNoiseSources (Xyce::Analysis::NoiseData & noiseData);
  void getNoiseSources (Xyce::Analysis::NoiseData & noiseData);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  bool auxChargeCalculations ();
  bool setupCapacitors_newDAE ();

  bool setupCapacitors_oldDAE ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupPointers ();

  bool setIC ();

  inline bool isConverged();

public:
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         ///< Owning model

  int dNode;
  int gNode;
  int sNode;
  int bNode;
  int dNodePrime;
  int sNodePrime;
  int qNode;

  double ueff;
  double thetavth;
  double von;
  double vdsat;
  double cgdo;
  double cgso;
  double vjsm;
  double IsEvjsm;
  double vjdm;
  double IsEvjdm;

  double l;
  double w;
  double numberParallel;   // PSpiceish "parallel copies" kludge
  double drainArea;
  double sourceArea;
  double drainSquares;
  double sourceSquares;
  double drainPerimeter;
  double sourcePerimeter;
  double sourceConductance;
  double drainConductance;

  double icVBS;
  double icVDS;
  double icVGS;
  bool OFF;
  int mode;
  int nqsMod;

  // OP point
  double qinv;
  double cd;
  double cbs;
  double cbd;
  double csub;
  double cdrain;
  double gm;
  double gds;
  double gmbs;
  double gbd;
  double gbs;

  double gbbs;
  double gbgs;
  double gbds;

  double cggb;
  double cgdb;
  double cgsb;
  double cbgb;
  double cbdb;
  double cbsb;
  double cdgb;
  double cddb;
  double cdsb;
  double capbd;
  double capbs;

  double cqgb;
  double cqdb;
  double cqsb;
  double cqbb;

  double qgate;
  double qbulk;
  double qdrn;

  double gtau;
  double gtg;
  double gtd;
  double gts;
  double gtb;

  double rds;  // Noise bugfix 
  double Vgsteff;
  double Vdseff;
  double Abulk;
  double AbovVgst2Vtm;

  bool limitedFlag;

  SizeDependParam<double>  * paramPtr;
  //SizeDependParam  * paramPtr;

  bool icVBSGiven;
  bool icVDSGiven;
  bool icVGSGiven;

  bool dNodePrimeSet;
  bool sNodePrimeSet;

  // end of original 3f5 stuff

  // Variables from the 3f5 b3ld function, which  were local to that
  // function but are more appropriate as instance variables:
  bool ChargeComputationNeeded;

  double gcbdb;
  double gcbgb;
  double gcbsb;
  double gcddb;
  double gcdgb;
  double gcdsb;
  double gcgdb;
  double gcggb;
  double gcgsb;
  double gcsdb;
  double gcsgb;
  double gcssb;

  double qgd, qgs, qgb;
  double qgdo, qgso;

  double qsrc, CoxWL;
  double Cgg, Cgd;
  double Cgb, Cdg, Cdd, Cds;
  double Csg, Csd, Css, Csb, Cbg, Cbd;
  double Cbb;

  // new DAE stuff:
  // Capacitance variables, which should correspond to the
  // conductance variables used by the old-DAE.
  //
  //  In general  gcggb -->  CAPcggb.
  //              gcgdb -->  CAPcgdb.  etc.
  double CAPcggb;
  double CAPcgdb;
  double CAPcgsb;
  double CAPcbgb;
  double CAPcbdb;
  double CAPcbsb;
  double CAPcdgb;
  double CAPcddb;
  double CAPcdsb;
  double CAPcsgb;
  double CAPcsdb;
  double CAPcssb;

  double Qeqqd_Jdxp;
  double Qeqqb_Jdxp;
  double Qeqqg_Jdxp;
  // end of new-DAE stuff.

  double dxpart;
  double sxpart;
  double ggtg;
  double ggtd;
  double ggts;
  double ggtb;
  double ddxpart_dVd;
  double ddxpart_dVg;
  double ddxpart_dVb;
  double ddxpart_dVs;
  double dsxpart_dVd;
  double dsxpart_dVg;
  double dsxpart_dVb;
  double dsxpart_dVs;

  double gbspsp;
  double gbbdp;
  double gbbsp;
  double gbspg;
  double gbspb;
  double gbspdp;
  double gbdpdp;
  double gbdpg;
  double gbdpb;
  double gbdpsp;

  double cdreq;
  double ceqbd;
  double ceqbs;
  double cdreq_Jdxp;
  double ceqbd_Jdxp;
  double ceqbs_Jdxp;

  double Gm;
  double Gmbs;
  double FwdSum;
  double RevSum;
  double T1global;

  double dVgst_dVg;
  double dVgst_dVb;
  double dVgs_eff_dVg;

  double dDeltaPhi_dVg;
  double dDeltaPhi_dVd;
  double dDeltaPhi_dVb;

  double gqdef;
  double gcqdb, gcqsb, gcqgb;
  double gcqbb;
  double ScalingFactor;
  double cqgate, cqbulk, cqdrn;

  // Variables which were model variables in the 3f5 version of the BSIM3,
  // but are temperature dependent and should be instance variables (as
  // temperature is an instance-level quantity.
  double vtm;
  double jctTempSatCurDensity;
  double jctSidewallTempSatCurDensity;

  double unitAreaJctCapTemp;
  double unitLengthSidewallJctCapTemp;
  double unitLengthGateSidewallJctCapTemp;
  double PhiBTemp;
  double PhiBSWTemp;
  double PhiBSWGTemp;
  // end of 3f5 stuff

  double temp;

  // solution variables, and intermediate quantities.
  //double mode;         // mode=1:normal mode.  mode=-1:inverse mode.
  double Vd;           // drain node voltage
  double Vs;           // source node voltage
  double Vg;           // gate node voltage
  double Vb;           // bulk node voltage
  double Vsp;          // source prime voltage
  double Vdp;          // drain prime voltage

  double Qtotal;       // total charge variable.

  double Vddp;         // voltage drop between drain and drain'
  double Vssp;         // voltage drop between source and source'

  double Vbsp;         // voltage drop, bulk-source prime
  double Vbdp;         // voltage drop, bulk-drain  prime

  double Vgsp;         // voltage drop, gate-Source prime
  double Vgdp;         // voltage drop, gate-Drain prime
  double Vgb;          // voltage drop, gate-Bulk

  double Vdpsp;        // voltage drop accross the channel
  double Vgt;          // Vgs-Vthresh

  // resistor currents:
  double Idrain;       // current through drain resistor
  double Isource;      // current through source resistor

  // channel current stuff:
  double df1dVdp;
  double df2dVdp;

  double df1dVsp;
  double df2dVsp;

  double df1dVg;
  double df2dVg;

  double df1dVb;
  double df2dVb;

  double vgb, vgd;
  double cqdef;
  double ceqqd;
  double ceqqb;
  double ceqqg;

  double cqdef_Jdxp;
  double ceqqd_Jdxp;
  double ceqqb_Jdxp;
  double ceqqg_Jdxp;

  // state variables, voltage drops
  double vbd;   // Volage drop, bulk-drain'
  double vbs;   // Volage drop, bulk-source'
  double vgs;   // Volage drop, gate-source'
  double vds;   // Volage drop, drain-source'

  // old versions of state variables, voltage drops
  // here "old" refers to the previous newton iteration.
  double vbd_old;   // Volage drop, bulk-drain'
  double vbs_old;   // Volage drop, bulk-source'
  double vgs_old;   // Volage drop, gate-source'
  double vds_old;   // Volage drop, drain-source'

  // "original" versions of various voltage drop variables:
  // original refers to the beginning of the newton iterations,
  // before any limits are imposed on the change in voltage drop.
  double vgs_orig;
  double vds_orig;
  double vbs_orig;
  double vbd_orig;
  double vgd_orig;

  // state variables, intrinsic capacitors
  double qb;
  double qg;
  double qd;

  // state variables, parasitic capacitors
  double qbs;
  double qbd;

  // state variables, cheq
  double qcheq;

  // state variables, cdump
  double qcdump;

  // state variable, qdef
  double qdef;

  // Indices into the state vector:

  // state variables, voltage drops
  int li_store_vbd;
  int li_store_vbs;
  int li_store_vgs;
  int li_store_vds;
  int li_store_von;

  // internal output vars, including transconductance:
  int li_store_gm;
  int li_store_Vds;
  int li_store_Vgs;
  int li_store_Vbs;
  int li_store_Vdsat;
  int li_store_Vth;
  int li_store_Gds;
  int li_store_Cgs;
  int li_store_Cgd;

  double Vds;
  double Vgs; 
  double Vbs;
  double Vdsat;
  double Vth;
  int li_branch_dev_id;
  int li_branch_dev_ig;
  int li_branch_dev_is;
  int li_branch_dev_ib;

  // state variables, intrinsic capacitors
  int li_state_qb;
  int li_state_qg;
  int li_state_qd;

  // state variables, parasitic capacitors
  int li_state_qbs;
  int li_state_qbd;

  // state variables, cheq
  int li_state_qcheq;

  // state variables, cdump
  int li_state_qcdump;

  // state variable, qdef
  int li_state_qdef;

  ////////////////////////////////////////////////////////////////////
  //  Local variable indices
  int li_Drain;
  int li_Gate;
  int li_Source;
  int li_Bulk;
  int li_DrainPrime;
  int li_SourcePrime;
  int li_Charge;

  // local indies
  int li_Ibs;
  int li_Ids;
  int li_Igs;

  ////////////////////////////////////////////////////////////////////
  //  Jacobian matrix indices:
  //  This is a 7x7 to 10x10 matrix block with 37 entries
  //
  // ---------------------------------------------------------------------------
  // | #NZ     |       |                                                        |
  // | entries |       |  V_d   V_g   V_s   V_b   V_d'  V_s'  Q  Ibs Ids Igs    |
  // ---------------------------------------------------------------------------|
  // |    3/2  | KCL_d |   a                       b                 i1         |
  // |    6/5  | KCL_g |         c           d     e     f    1          i2     |
  // |    5/2  | KCL_s |               g                 h       i3  i4  i5     |
  // |    6/5  | KCL_b |         i           j     k     l    2  i6             |
  // |    6    | KCL_d'|   m     n           o     p     q    3                 |
  // |    6    | KCL_s'|         r     s     t     u     v    4                 |
  // |    5    | Q_equ |         7           9     6     8    5                 |
  // |    2/1  | icVBS |               i7    i8                  (i9)           |
  // |    2/1  | icVDS |   i10         i11                          (i12)       |
  // |    2/1  | icVGS |         i13   i14                              (i15)   |
  // ---------------------------------------------------------------------------
  //     43 at operating point, 34 normally
  //
  //     Terms beginning with "i", not in parenthesis are for
  //     initial conditions and apply only at the operating point.
  //     These terms reduce to the terms in parenthesis during time
  //     integration -- i.e. the other "i" terms go to zero during
  //     time integration.
  ////////////////////////////////////////////////////////////////////

  // Offset variables corresponding to the above declared indices.

  // Jacobian Matrix Offsets:

  // V_d Row:
  int ADrainEquDrainNodeOffset;             // a
  int ADrainEquDrainPrimeNodeOffset;        // b
  int ADrainEquIdsOffset;                   // i1

  // V_g Row:
  int AGateEquGateNodeOffset;               // c
  int AGateEquBulkNodeOffset;               // d
  int AGateEquDrainPrimeNodeOffset;         // e
  int AGateEquSourcePrimeNodeOffset;        // f
  int AGateEquChargeVarOffset;              // 1
  int AGateEquIgsOffset;                    // i2

  // V_s Row:
  int ASourceEquSourceNodeOffset;           // g
  int ASourceEquSourcePrimeNodeOffset;      // h
  int ASourceEquIbsOffset;                  // i3
  int ASourceEquIdsOffset;                  // i4
  int ASourceEquIgsOffset;                  // i5

  // V_b Row:
  int ABulkEquGateNodeOffset;               // i
  int ABulkEquBulkNodeOffset;               // j
  int ABulkEquDrainPrimeNodeOffset;         // k
  int ABulkEquSourcePrimeNodeOffset;        // l
  int ABulkEquChargeVarOffset;              // 2
  int ABulkEquIbsOffset;                    // i6

  // V_d' Row:
  int ADrainPrimeEquDrainNodeOffset;        // m
  int ADrainPrimeEquGateNodeOffset;         // n
  int ADrainPrimeEquBulkNodeOffset;         // o
  int ADrainPrimeEquDrainPrimeNodeOffset;   // p
  int ADrainPrimeEquSourcePrimeNodeOffset;  // q
  int ADrainPrimeEquChargeVarOffset;        // 3

  // V_s' Row:
  int ASourcePrimeEquGateNodeOffset;        // r
  int ASourcePrimeEquSourceNodeOffset;      // s
  int ASourcePrimeEquBulkNodeOffset;        // t
  int ASourcePrimeEquDrainPrimeNodeOffset;  // u
  int ASourcePrimeEquSourcePrimeNodeOffset; // v
  int ASourcePrimeEquChargeVarOffset;       // 4

  // MOSFET charge (Q) Row:
  int AChargeEquChargeVarOffset;            // 5
  int AChargeEquDrainPrimeNodeOffset;       // 6
  int AChargeEquGateNodeOffset;             // 7
  int AChargeEquSourcePrimeNodeOffset;      // 8
  int AChargeEquBulkNodeOffset;             // 9

  // icVBS
  int icVBSEquVsOffset;                     // i7
  int icVBSEquVbOffset;                     // i8
  int icVBSEquIbsOffset;                    // i9

  // icVDS
  int icVDSEquVdOffset;                     // i10
  int icVDSEquVsOffset;                     // i11
  int icVDSEquIdsOffset;                    // i12

  // icVGS
  int icVGSEquVgOffset;                     // i13
  int icVGSEquVsOffset;                     // i14
  int icVGSEquIgsOffset;                    // i15

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ///////////////////////////////////////////////////
  // Jacobian Matrix Pointers:

  // V_d Row:
  double * f_DrainEquDrainNodePtr;             // a
  double * f_DrainEquDrainPrimeNodePtr;        // b
  double * f_DrainEquIdsPtr;                   // i1

  // V_g Row:
  double * f_GateEquGateNodePtr;               // c
  double * f_GateEquBulkNodePtr;               // d
  double * f_GateEquDrainPrimeNodePtr;         // e
  double * f_GateEquSourcePrimeNodePtr;        // f
  double * f_GateEquChargeVarPtr;              // 1
  double * f_GateEquIgsPtr;                    // i2

  // V_s Row:
  double * f_SourceEquSourceNodePtr;           // g
  double * f_SourceEquSourcePrimeNodePtr;      // h
  double * f_SourceEquIbsPtr;                  // i3
  double * f_SourceEquIdsPtr;                  // i4
  double * f_SourceEquIgsPtr;                  // i5

  // V_b Row:
  double * f_BulkEquGateNodePtr;               // i
  double * f_BulkEquBulkNodePtr;               // j
  double * f_BulkEquDrainPrimeNodePtr;         // k
  double * f_BulkEquSourcePrimeNodePtr;        // l
  double * f_BulkEquChargeVarPtr;              // 2
  double * f_BulkEquIbsPtr;                    // i6

  // V_d' Row:
  double * f_DrainPrimeEquDrainNodePtr;        // m
  double * f_DrainPrimeEquGateNodePtr;         // n
  double * f_DrainPrimeEquBulkNodePtr;         // o
  double * f_DrainPrimeEquDrainPrimeNodePtr;   // p
  double * f_DrainPrimeEquSourcePrimeNodePtr;  // q
  double * f_DrainPrimeEquChargeVarPtr;        // 3

  // V_s' Row:
  double * f_SourcePrimeEquGateNodePtr;        // r
  double * f_SourcePrimeEquSourceNodePtr;      // s
  double * f_SourcePrimeEquBulkNodePtr;        // t
  double * f_SourcePrimeEquDrainPrimeNodePtr;  // u
  double * f_SourcePrimeEquSourcePrimeNodePtr; // v
  double * f_SourcePrimeEquChargeVarPtr;       // 4

  // MOSFET charge (Q) Row:
  double * f_ChargeEquChargeVarPtr;            // 5
  double * f_ChargeEquDrainPrimeNodePtr;       // 6
  double * f_ChargeEquGateNodePtr;             // 7
  double * f_ChargeEquSourcePrimeNodePtr;      // 8
  double * f_ChargeEquBulkNodePtr;             // 9

  // icVBS
  double * f_icVBSEquVsPtr;                   // i7
  double * f_icVBSEquVbPtr;                   // i8
  double * f_icVBSEquIbsPtr;                  // i9

  // icVDS
  double * f_icVDSEquVdPtr;                   // i10
  double * f_icVDSEquVsPtr;                   // i11
  double * f_icVDSEquIdsPtr;                  // i12

  // icVGS
  double * f_icVGSEquVgPtr;                   // i13
  double * f_icVGSEquVsPtr;                   // i14
  double * f_icVGSEquIgsPtr;                  // i15

  // V_d Row:
  double * q_DrainEquDrainNodePtr;             // a
  double * q_DrainEquDrainPrimeNodePtr;        // b
  double * q_DrainEquIdsPtr;                   // i1

  // V_g Row:
  double * q_GateEquGateNodePtr;               // c
  double * q_GateEquBulkNodePtr;               // d
  double * q_GateEquDrainPrimeNodePtr;         // e
  double * q_GateEquSourcePrimeNodePtr;        // f
  double * q_GateEquChargeVarPtr;              // 1
  double * q_GateEquIgsPtr;                    // i2

  // V_s Row:
  double * q_SourceEquSourceNodePtr;           // g
  double * q_SourceEquSourcePrimeNodePtr;      // h
  double * q_SourceEquIbsPtr;                  // i3
  double * q_SourceEquIdsPtr;                  // i4
  double * q_SourceEquIgsPtr;                  // i5

  // V_b Row:
  double * q_BulkEquGateNodePtr;               // i
  double * q_BulkEquBulkNodePtr;               // j
  double * q_BulkEquDrainPrimeNodePtr;         // k
  double * q_BulkEquSourcePrimeNodePtr;        // l
  double * q_BulkEquChargeVarPtr;              // 2
  double * q_BulkEquIbsPtr;                    // i6

  // V_d' Row:
  double * q_DrainPrimeEquDrainNodePtr;        // m
  double * q_DrainPrimeEquGateNodePtr;         // n
  double * q_DrainPrimeEquBulkNodePtr;         // o
  double * q_DrainPrimeEquDrainPrimeNodePtr;   // p
  double * q_DrainPrimeEquSourcePrimeNodePtr;  // q
  double * q_DrainPrimeEquChargeVarPtr;        // 3

  // V_s' Row:
  double * q_SourcePrimeEquGateNodePtr;        // r
  double * q_SourcePrimeEquSourceNodePtr;      // s
  double * q_SourcePrimeEquBulkNodePtr;        // t
  double * q_SourcePrimeEquDrainPrimeNodePtr;  // u
  double * q_SourcePrimeEquSourcePrimeNodePtr; // v
  double * q_SourcePrimeEquChargeVarPtr;       // 4

  // MOSFET charge (Q) Row:
  double * q_ChargeEquChargeVarPtr;            // 5
  double * q_ChargeEquDrainPrimeNodePtr;       // 6
  double * q_ChargeEquGateNodePtr;             // 7
  double * q_ChargeEquSourcePrimeNodePtr;      // 8
  double * q_ChargeEquBulkNodePtr;             // 9

  // icVBS
  double * q_icVBSEquVsPtr;                   // i7
  double * q_icVBSEquVbPtr;                   // i8
  double * q_icVBSEquIbsPtr;                  // i9

  // icVDS
  double * q_icVDSEquVdPtr;                   // i10
  double * q_icVDSEquVsPtr;                   // i11
  double * q_icVDSEquIdsPtr;                  // i12

  // icVGS
  double * q_icVGSEquVgPtr;                   // i13
  double * q_icVGSEquVsPtr;                   // i14
  double * q_icVGSEquIgsPtr;                  // i15
#endif

  // flag for updateTemperature call.  Needed for .STEP temperature sweeps.
  bool updateTemperatureCalled_;

  static std::vector< std::vector<int> > jacStamp_DC_SC;
  static std::vector< std::vector<int> > jacStamp_DC;
  static std::vector< std::vector<int> > jacStamp_SC;
  static std::vector< std::vector<int> > jacStamp;

  static std::vector<int> jacMap_DC_SC;
  static std::vector<int> jacMap_DC;
  static std::vector<int> jacMap_SC;
  static std::vector<int> jacMap;

  static std::vector< std::vector<int> > jacMap2_DC_SC;
  static std::vector< std::vector<int> > jacMap2_DC;
  static std::vector< std::vector<int> > jacMap2_SC;
  static std::vector< std::vector<int> > jacMap2;

  int blockHomotopyID; // For homotopy
  double randomPerturb; // For homotopy

  // Most versions of this device will use one of the four
  // jacobian's defined as static above.  However, if nqsMod is true
  // or or any of the three initial conditions are specified (icVBS,
  // icVDS or icVGS), then the number of jacobian permutations goes
  // up to 112 (4 from static versions above * 2 for nqsMod on/off *
  // 7 for IC's = 4 * 2 * 7 = 56) Rather than enumerate all of these
  // as static functions, we'll make non-static, member variables of
  // the jacStamp for those few devices that specify nqsMod or
  // initial conditoins.

  std::vector< std::vector< int > > jacStampSpecial;    // for simulation

  // used when we need to merge the d' and or s' nodes
  std::vector< std::vector<int> > jacStampSpecialMerged;
  std::vector<int>           jacSpecialMap;
  std::vector< std::vector<int> > jacSpecialMap2;
  std::vector<int>           jacSpecialMergedMap;
  std::vector< std::vector<int> > jacSpecialMergedMap2;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/14/00
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class SizeDependParam<double>;
  //friend class SizeDependParam;
  friend class ParametricData<Model>;
  friend class Instance;
  friend struct Traits;
  friend class Master;
  friend class bsim3InstanceSensitivity;
  friend class bsim3ModelSensitivity;

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

  bool clearTemperatureData ();

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

  // 3f5 stuff:
  int    modType;
  int    dtype;

  int    mobMod;
  int    capMod;
  int    noiMod;
  int    binUnit;
  int    paramChk;

  std::string version;

  bool binPrefixFlag;
  // Bogus PSpice-y thing where you can set L and W on the model *or* the
  // instance line, with instance overriding model
  double model_l;
  double model_w;
  //

  double tox;
  double toxm;
  double cdsc;
  double cdscb;
  double cdscd;
  double cit;
  double nfactor;
  double xj;
  double vsat;
  double at;
  double a0;
  double ags;
  double a1;
  double a2;
  double keta;
  double nsub;
  double npeak;
  double ngate;
  double gamma1;
  double gamma2;
  double vbx;
  double vbm;
  double xt;
  double k1;
  double kt1;
  double kt1l;
  double kt2;
  double k2;
  double k3;
  double k3b;
  double w0;
  double nlx;
  double dvt0;
  double dvt1;
  double dvt2;
  double dvt0w;
  double dvt1w;
  double dvt2w;
  double drout;
  double dsub;
  double vth0;
  double ua;
  double ua1;
  double ub;
  double ub1;
  double uc;
  double uc1;
  double u0;
  double ute;
  double voff;
  double delta;
  double rdsw;
  double prwg;
  double prwb;
  double prt;
  double eta0;
  double etab;
  double pclm;
  double pdibl1;
  double pdibl2;
  double pdiblb;
  double pscbe1;
  double pscbe2;
  double pvag;
  double wr;
  double dwg;
  double dwb;
  double b0;
  double b1;
  double alpha0;
  double alpha1;
  double beta0;
  double ijth;
  double vfb;

  // CV model
  double elm;
  double cgsl;
  double cgdl;
  double ckappa;
  double cf;
  double vfbcv;
  double clc;
  double cle;
  double dwc;
  double dlc;
  double noff;
  double voffcv;
  double acde;
  double moin;
  double tcj;
  double tcjsw;
  double tcjswg;
  double tpb;
  double tpbsw;
  double tpbswg;

  // Length Dependence
  double lcdsc;
  double lcdscb;
  double lcdscd;
  double lcit;
  double lnfactor;
  double lxj;
  double lvsat;
  double lat;
  double la0;
  double lags;
  double la1;
  double la2;
  double lketa;
  double lnsub;
  double lnpeak;
  double lngate;
  double lgamma1;
  double lgamma2;
  double lvbx;
  double lvbm;
  double lxt;
  double lk1;
  double lkt1;
  double lkt1l;
  double lkt2;
  double lk2;
  double lk3;
  double lk3b;
  double lw0;
  double lnlx;
  double ldvt0;
  double ldvt1;
  double ldvt2;
  double ldvt0w;
  double ldvt1w;
  double ldvt2w;
  double ldrout;
  double ldsub;
  double lvth0;
  double lua;
  double lua1;
  double lub;
  double lub1;
  double luc;
  double luc1;
  double lu0;
  double lute;
  double lvoff;
  double ldelta;
  double lrdsw;
  double lprwg;
  double lprwb;
  double lprt;
  double leta0;
  double letab;
  double lpclm;
  double lpdibl1;
  double lpdibl2;
  double lpdiblb;
  double lpscbe1;
  double lpscbe2;
  double lpvag;
  double lwr;
  double ldwg;
  double ldwb;
  double lb0;
  double lb1;
  double lalpha0;
  double lalpha1;
  double lbeta0;
  double lvfb;

  // CV model
  double lelm;
  double lcgsl;
  double lcgdl;
  double lckappa;
  double lcf;
  double lclc;
  double lcle;
  double lvfbcv;
  double lnoff;
  double lvoffcv;
  double lacde;
  double lmoin;

  // Width Dependence
  double wcdsc;
  double wcdscb;
  double wcdscd;
  double wcit;
  double wnfactor;
  double wxj;
  double wvsat;
  double wat;
  double wa0;
  double wags;
  double wa1;
  double wa2;
  double wketa;
  double wnsub;
  double wnpeak;
  double wngate;
  double wgamma1;
  double wgamma2;
  double wvbx;
  double wvbm;
  double wxt;
  double wk1;
  double wkt1;
  double wkt1l;
  double wkt2;
  double wk2;
  double wk3;
  double wk3b;
  double ww0;
  double wnlx;
  double wdvt0;
  double wdvt1;
  double wdvt2;
  double wdvt0w;
  double wdvt1w;
  double wdvt2w;
  double wdrout;
  double wdsub;
  double wvth0;
  double wua;
  double wua1;
  double wub;
  double wub1;
  double wuc;
  double wuc1;
  double wu0;
  double wute;
  double wvoff;
  double wdelta;
  double wrdsw;
  double wprwg;
  double wprwb;
  double wprt;
  double weta0;
  double wetab;
  double wpclm;
  double wpdibl1;
  double wpdibl2;
  double wpdiblb;
  double wpscbe1;
  double wpscbe2;
  double wpvag;
  double wwr;
  double wdwg;
  double wdwb;
  double wb0;
  double wb1;
  double walpha0;
  double walpha1;
  double wbeta0;
  double wvfb;

  // CV model
  double welm;
  double wcgsl;
  double wcgdl;
  double wckappa;
  double wcf;
  double wclc;
  double wcle;
  double wvfbcv;
  double wnoff;
  double wvoffcv;
  double wacde;
  double wmoin;

  // Cross-term Dependence
  double pcdsc;
  double pcdscb;
  double pcdscd;
  double pcit;
  double pnfactor;
  double pxj;
  double pvsat;
  double pat;
  double pa0;
  double pags;
  double pa1;
  double pa2;
  double pketa;
  double pnsub;
  double pnpeak;
  double pngate;
  double pgamma1;
  double pgamma2;
  double pvbx;
  double pvbm;
  double pxt;
  double pk1;
  double pkt1;
  double pkt1l;
  double pkt2;
  double pk2;
  double pk3;
  double pk3b;
  double pw0;
  double pnlx;
  double pdvt0;
  double pdvt1;
  double pdvt2;
  double pdvt0w;
  double pdvt1w;
  double pdvt2w;
  double pdrout;
  double pdsub;
  double pvth0;
  double pua;
  double pua1;
  double pub;
  double pub1;
  double puc;
  double puc1;
  double pu0;
  double pute;
  double pvoff;
  double pdelta;
  double prdsw;
  double pprwg;
  double pprwb;
  double pprt;
  double peta0;
  double petab;
  double ppclm;
  double ppdibl1;
  double ppdibl2;
  double ppdiblb;
  double ppscbe1;
  double ppscbe2;
  double ppvag;
  double pwr;
  double pdwg;
  double pdwb;
  double pb0;
  double pb1;
  double palpha0;
  double palpha1;
  double pbeta0;
  double pvfb;

  // CV model
  double pelm;
  double pcgsl;
  double pcgdl;
  double pckappa;
  double pcf;
  double pclc;
  double pcle;
  double pvfbcv;
  double pnoff;
  double pvoffcv;
  double pacde;
  double pmoin;

  double tnom;
  double cgso;
  double cgdo;
  double cgbo;
  double xpart;
  double cFringOut;
  double cFringMax;

  double sheetResistance;
  double jctSatCurDensity;
  double jctSidewallSatCurDensity;
  double bulkJctPotential;
  double bulkJctBotGradingCoeff;
  double bulkJctSideGradingCoeff;
  double bulkJctGateSideGradingCoeff;
  double sidewallJctPotential;
  double GatesidewallJctPotential;
  double unitAreaJctCap;
  double unitLengthSidewallJctCap;
  double unitLengthGateSidewallJctCap;
  double jctEmissionCoeff;
  double jctTempExponent;

  double Lint;
  double Ll;
  double Llc;
  double Lln;
  double Lw;
  double Lwc;
  double Lwn;
  double Lwl;
  double Lwlc;
  double Lmin;
  double Lmax;

  double Wint;
  double Wl;
  double Wlc;
  double Wln;
  double Ww;
  double Wwc;
  double Wwn;
  double Wwl;
  double Wwlc;
  double Wmin;
  double Wmax;


  // Pre-calculated constants
  // MCJ: move to size-dependent param.
  double vtm;
  double cox;
  double cof1;
  double cof2;
  double cof3;
  double cof4;
  double vcrit;
  double factor1;
  double PhiB;
  double PhiBSW;
  double PhiBSWG;

  double oxideTrapDensityA;
  double oxideTrapDensityB;
  double oxideTrapDensityC;
  double em;
  double ef;
  double af;
  double kf;
  double lintnoi;

  std::list<SizeDependParam<double> *> sizeDependParamList;
  //std::list<SizeDependParam*> sizeDependParamList;

  bool  npeakGiven;
  bool  gamma1Given;
  bool  gamma2Given;
  bool  k1Given;
  bool  k2Given;
  bool  nsubGiven;
  bool  xtGiven;
  bool  vbxGiven;
  bool  vbmGiven;
  bool  vfbGiven;
  bool  vth0Given;

  // Variables from the 3f5 BSIM3 function b3temp, but are  more
  // appropriate as model variables.
  double Vtm0;
  double Eg0;
  double ni;
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
  friend class SizeDependParam<double>;
  //friend class SizeDependParam;
public:
  Master(
     const Configuration &       configuration,
     const FactoryBlock &      factory_block,
     const SolverState & ss1,
     const DeviceOptions & do1)
    : DeviceMaster<Traits>(configuration, factory_block, ss1, do1)
  {}

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);

  // load functions, residual:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);

  // load functions, Jacobian:
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);

};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);


} // namespace MOSFET_B3
} // namespace Device
} // namespace Xyce

#endif


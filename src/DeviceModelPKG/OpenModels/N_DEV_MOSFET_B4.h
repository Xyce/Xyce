//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Creation Date  : 11/25/06
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MOSFET_B4_h
#define Xyce_N_DEV_MOSFET_B4_h

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
namespace MOSFET_B4 {

enum noiseType {RDNOIZ, RSNOIZ, RGNOIZ, RBPSNOIZ, RBPDNOIZ, RBPBNOIZ, RBSBNOIZ, RBDBNOIZ, IDNOIZ, FLNOIZ, IGSNOIZ, IGDNOIZ, IGBNOIZ, NUMNOIZ };

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, MOSFET1::Traits>
{
  static const char *name() {return  "BSIM4";}
  static const char *deviceTypeName() {return "M level 14";}
  static int numNodes() {return 4;}
  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Struct        : SizeDependParam
// Purpose       : copied over from the 3f5 code.  Almost no changes.
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
class SizeDependParam
{
  friend class Model;
  friend class Instance;
  friend class Master;

private:
  double Width;
  double Length;
  double NFinger;

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
  double ndep;
  double nsd;
  double phin;
  double ngate;
  double gamma1;
  double gamma2;
  double vbx;
  double vbi;
  double vbm;
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
  double dvtp0;
  double dvtp1;
  double lpe0;
  double lpeb;
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
  double ud;
  double ud1;
  double up;
  double lp;
  double u0;
  double eu;
  double ute;
  double voff;
  double tvoff;
  double minv;
  double minvcv;
  double vfb;
  double delta;
  double rdsw;
  double rds0;
  double rs0;
  double rd0;
  double rsw;
  double rdw;
  double prwg;
  double prwb;
  double prt;
  double eta0;
  double etab;
  double pclm;
  double pdibl1;
  double pdibl2;
  double pdiblb;
  double fprout;
  double pdits;
  double pditsd;
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
  double agidl;
  double bgidl;
  double cgidl;
  double egidl;
  double agisl;
  double bgisl;
  double cgisl;
  double egisl;
  double aigc;
  double bigc;
  double cigc;
  double aigs;
  double bigs;
  double cigs;
  double aigd;
  double bigd;
  double cigd;
  double aigbacc;
  double bigbacc;
  double cigbacc;
  double aigbinv;
  double bigbinv;
  double cigbinv;
  double nigc;
  double nigbacc;
  double nigbinv;
  double ntox;
  double eigbinv;
  double pigcd;
  double poxedge;
  double xrcrg1;
  double xrcrg2;
  double lambda; // overshoot
  double vtl; // thermal velocity limit
  double xn; // back scattering parameter
  double lc; // back scattering parameter
  double tfactor;  // ballistic transportation factor
  double vfbsdoff;  // S/D flatband offset voltage
  double tvfbsdoff;

  // added for stress effect
  double ku0;
  double kvth0;
  double ku0temp;
  double rho_ref;
  double inv_od_ref;
  // added for well proximity effect
  double kvth0we;
  double k2we;
  double ku0we;

  // CV model
  double cgsl;
  double cgdl;
  double ckappas;
  double ckappad;
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
  double dwj;
  double leffCV;
  double weffCV;
  double weffCJ;
  double abulkCVfactor;
  double cgso;
  double cgdo;
  double cgbo;

  double u0temp;
  double vsattemp;
  double sqrtPhi;
  double phis3;
  double Xdep0;
  double sqrtXdep0;
  double theta0vb0;
  double thetaRout;
  double mstar;
  double mstarcv;
  double voffcbn;
  double voffcbncv;
  double rdswmin;
  double rdwmin;
  double rswmin;
  double vfbsd;

  double cof1;
  double cof2;
  double cof3;
  double cof4;
  double cdep0;
  double ToxRatio;
  double Aechvb;
  double Bechvb;
  double ToxRatioEdge;
  double AechvbEdgeS;
  double AechvbEdgeD;
  double BechvbEdge;
  double ldeb;
  double k1ox;
  double k2ox;
  double vfbzbfactor;

  // ERK.  I added this to make temperature sweeps work.
  double referenceTemperature;
};


//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class SizeDependParam;
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;

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
  void registerStoreLIDs(const std::vector<int> & stoLIDVecRef);
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  void setupPointers ();

  bool processParams ();

  void debugJacStampOutput ();

  bool updateTemperature(const double & temp_tmp);
  bool updateIntermediateVars ();
  bool updatePrimaryState ();


  double Eval1ovFNoise(double Vds, double freq, double temp);
  int getNumNoiseSources () const;
  void setupNoiseSources (Xyce::Analysis::NoiseData & noiseData);
  void getNoiseSources (Xyce::Analysis::NoiseData & noiseData);

  int polyDepletion(
     double phi, double ngate, double epsgate,
     double coxe, double Vgs_arg,
     double & Vgs_eff, double & dVgs_eff_dVg);

  int NumFingerDiff(
     double nf_arg, int minSD,
     double & nuIntD, double & nuEndD, double & nuIntS, double & nuEndS);

  int PAeffGeo(
     double nf_arg, int geo, int minSD,
     double Weffcj, double DMCG, double DMCI, double DMDG,
     double & Ps, double & Pd, double & As, double & Ad);

  int RdseffGeo(
     double nf_arg,
     int geo, int rgeo, int minSD,
     double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
     int Type, double & Rtot);

  int RdsEndIso(
     double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
     double nuEnd, int rgeo, int Type, double & Rend);

  int RdsEndSha(
     double Weffcj, double Rsh, double DMCG, double DMCI, double DMDG,
     int rgeo, int Type, double nuEnd, double & Rend);

  int DioIjthVjmEval (
     double Nvtm, double Ijth,
     double Isb, double XExpBV, double & Vjm);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void setupFVectorVars ();

  bool auxChargeCalculations ();
  bool setupCapacitors_newDAE ();

  bool setupCapacitors_oldDAE ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  bool setIC ();

  inline bool isConverged();

  // Beginning of 3f5 stuff:
  // attributes:
public:
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  double ueff;
  double thetavth;
  double von;
  double vdsat;
  double cgdo;
  double qgdo;
  double cgso;
  double qgso;
  double grbsb;
  double grbdb;
  double grbpb;
  double grbps;
  double grbpd;

  double vjsmFwd;
  double vjsmRev;
  double vjdmFwd;
  double vjdmRev;
  double XExpBVS;
  double XExpBVD;
  double SslpFwd;
  double SslpRev;
  double DslpFwd;
  double DslpRev;
  double IVjsmFwd;
  double IVjsmRev;
  double IVjdmFwd;
  double IVjdmRev;

  double grgeltd;
  double Pseff;
  double Pdeff;
  double Aseff;
  double Adeff;

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
  // stress effect instance param
  double sa;
  double sb;
  double sd;
  bool SDgiven;
  double sca;
  double scb;
  double scc;
  double sc;

  double rbdb;
  double rbsb;
  double rbpb;
  double rbps;
  double rbpd;

  bool RBDBgiven;
  bool RBSBgiven;
  bool RBPBgiven;
  bool RBPSgiven;
  bool RBPDgiven;

  double delvto;
  double xgw;
  double ngcon;

  bool XGWgiven;
  bool NGCONgiven;

  // added here to account stress effect instance dependence
  double u0temp;
  double vsattemp;
  double vth0;
  double vfb;
  double vfbzb;
  double vtfbphi1;
  double vtfbphi2;
  double k2;
  double vbsc;
  double k2ox;
  double eta0;

  double icVDS;
  double icVGS;
  double icVBS;
  double nf;
  bool OFF;
  int mode;
  int trnqsMod;
  int acnqsMod;
  int rbodyMod;
  int rgateMod;
  int geoMod;
  int rgeoMod;
  int min;

  bool RBODYMODgiven;
  bool RGATEMODgiven;
  bool GEOMODgiven;
  bool RGEOMODgiven;
  bool TRNQSMODgiven;
  bool ACNQSMODgiven;

  // OP point
  double Vgsteff;
  double Vgsteff_forNoise;  // ERK.  see similar comment about Abulk.
  double vgs_eff;
  double vgd_eff;
  double dvgs_eff_dvg;
  double dvgd_eff_dvg;
  double Vdseff;
  double Vdseff_forNoise; // ERK.  see similar comment about Abulk.
  double nstar;
  double Abulk;
  double Abulk_forNoise; // ERK. this is a bit screwy.  But in spice3/ngspice, Abulk is saved as an instance variable and used later in noise calculations.  But, within the load function the "local" copy of Abulk gets modified further afterwards.  So it is then "wrong" for noise at that point.  As we've made Abulk a class variable, there needs to be an extra copy for noise, saved at the right time.
  double EsatL;
  double AbovVgst2Vtm;
  double qinv;
  double cd;
  double cbs;
  double cbd;
  double csub;
  double Igidl;
  double Igisl;
  double gm;
  double gds;
  double gmbs;
  double gbd;
  double gbs;

  double gbbs;
  double gbgs;
  double gbds;
  double ggidld;
  double ggidlg;
  double ggidls;
  double ggidlb;
  double ggisld;
  double ggislg;
  double ggisls;
  double ggislb;

  double Igcs;
  double gIgcsg;
  double gIgcsd;
  double gIgcss;
  double gIgcsb;
  double Igcd;
  double gIgcdg;
  double gIgcdd;
  double gIgcds;
  double gIgcdb;

  double Igs;
  double gIgsg;
  double gIgss;
  double Igd;
  double gIgdg;
  double gIgdd;

  double Igb;
  double gIgbg;
  double gIgbd;
  double gIgbs;
  double gIgbb;

  double grdsw;
  double IdovVds;
  double gcrg;
  double gcrgd;
  double gcrgg;
  double gcrgs;
  double gcrgb;

  double gstot;
  double gstotd;
  double gstotg;
  double gstots;
  double gstotb;

  double gdtot;
  double gdtotd;
  double gdtotg;
  double gdtots;
  double gdtotb;

  double cggb;
  double cgdb;
  double cgsb;
  double cbgb;
  double cbdb;
  double cbsb;
  double cdgb;
  double cddb;
  double cdsb;
  double csgb;
  double csdb;
  double cssb;
  double cgbb;
  double cdbb;
  double csbb;
  double cbbb;
  double capbd;
  double capbs;

  double cqgb;
  double cqdb;
  double cqsb;
  double cqbb;

  double qgate;
  double qbulk;
  double qdrn;
  double qsrc;

  double qchqs;
  double taunet;
  double gtau;
  double gtg;
  double gtd;
  double gts;
  double gtb;
  double SjctTempRevSatCur;
  double DjctTempRevSatCur;
  double SswTempRevSatCur;
  double DswTempRevSatCur;
  double SswgTempRevSatCur;
  double DswgTempRevSatCur;

  bool limitedFlag;

  SizeDependParam  * paramPtr;

  bool icVBSGiven;
  bool icVDSGiven;
  bool icVGSGiven;

  bool scaGiven;
  bool scbGiven;
  bool sccGiven;
  bool scGiven;
  bool sourcePerimeterGiven;
  bool drainPerimeterGiven;
  bool sourceAreaGiven;
  bool drainAreaGiven;
  bool sourceSquaresGiven;
  bool drainSquaresGiven;

  bool drainMOSFET_B4Exists;
  bool sourceMOSFET_B4Exists;

  // Variables that were local to the b4ld function, but are
  // more appropriate as class variables.
  bool ChargeComputationNeeded;

  double temp;
  bool TEMPgiven;

  // solution variables, and intermediate quantities.
  double Vd;           // drain node voltage
  double Vs;           // source node voltage
  double Vb;           // bulk node voltage

  double Vdp;          // drain prime voltage
  double Vsp;          // source prime voltage
  double Vgp;          // gate prime voltage
  double Vbp;          // bulk prime voltage

  double Vge;
  double Vgm;
  double Vdb;
  double Vsb;
  double Vds;
  double Vgs;
  double Vbs;

  double Qtotal;       // total charge variable.

  double Vddp;         // voltage drop between drain and drain'
  double Vssp;         // voltage drop between source and source'

  double Vbsp;         // voltage drop, bulk-source prime
  double Vbdp;         // voltage drop, bulk-drain  prime

  double Vgsp;         // voltage drop, gate-Source prime
  double Vgdp;         // voltage drop, gate-Drain prime
  double Vgb;          // voltage drop, gate-Bulk

  double Vdpsp;        // voltage drop accross the channel

  // Voltage drops in the substrate network:
  double Vdbb;         // drain-body/body
  double Vdbbp;        // drain-body/body-prime
  double Vbpb;         // body-prime/body
  double Vsbb;         // source-body/body
  double Vsbbp;        // source-body/body-prime

  // resistor currents:
  double Idrain;       // current through drain resistor
  double Isource;      // current through source resistor
  double Idbb;         // current through drain-body/body resistor
  double Idbbp;        // current through drain-body/body-prime resistor
  double Ibpb;         // current through body-prime/body resistor
  double Isbb;         // current through source-body/body resistor
  double Isbbp;        // current through source-body/body-prime resistor

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
  double vbd;
  double vbs;
  double vgs;
  double vds;
  double vges;
  double vgms;
  double vdes;
  double vses;
  double vdbs;
  double vsbs;
  double vdbd;
  double vged;
  double vgmd;

  // old versions of state variables, voltage drops
  // here "old" refers to the previous newton iteration.
  double vbd_old;
  double vbs_old;
  double vgs_old;
  double vds_old;
  double vges_old;
  double vgms_old;
  double vdes_old;
  double vses_old;
  double vdbs_old;
  double vsbs_old;
  double vdbd_old;
  double vged_old;
  double vgmd_old;

  // "original" versions of various voltage drop variables:
  // original refers to the beginning of the newton iterations,
  // before any limits are imposed on the change in voltage drop.
  double vbd_orig;
  double vbs_orig;
  double vgs_orig;
  double vds_orig;
  double vgd_orig;
  double vges_orig;
  double vgms_orig;
  double vdes_orig;
  double vses_orig;
  double vdbs_orig;
  double vsbs_orig;
  double vdbd_orig;
  double vbs_jct_orig;
  double vbd_jct_orig;
  double vgmb_orig;
  double vgb_orig;
  double vged_orig;
  double vgmd_orig;

  // b4ld function variables
  double Gm, Gmbs, FwdSum, RevSum;
  double ceqdrn, cdrain, ceqbd, ceqbs;
  double cqgate, cqbody, cqdrn;

  double gbbdp, gbbsp, gbdpg, gbdpdp, gbdpb;
  double gbdpsp, gbspg, gbspdp, gbspb, gbspsp;

  double Istoteq, gIstotg, gIstotd, gIstots, gIstotb;
  double Idtoteq, gIdtotg, gIdtotd, gIdtots, gIdtotb;
  double Ibtoteq, gIbtotg, gIbtotd, gIbtots, gIbtotb;
  double Igtoteq, gIgtotg, gIgtotd, gIgtots, gIgtotb;

  double ceqgcrg, ceqgstot, ceqgdtot;
  double ceqjs, ceqjd;
  double vbs_jct, vbd_jct;
  double ceqqjs, ceqqjd;
  double ceqqgmid;
  double gjbd, gjbs, gdpr, gspr, geltd, gcggb;
  double ggtg, gcgdb, ggtd, gcgsb;
  double ggts, gcgbb, ggtb;
  double dxpart;
  double sxpart;
  double gqdef;

  double ddxpart_dVd, ddxpart_dVg, ddxpart_dVb, ddxpart_dVs;
  double dsxpart_dVd, dsxpart_dVg, dsxpart_dVb, dsxpart_dVs;
  double vgmb;

  double CoxWL, ScalingFactor;
  double DMCGeff, DMCIeff, DMDGeff;

  double ceqqgmid_Jdxp;
  double ceqqjs_Jdxp;
  double ceqqjd_Jdxp;

  double qgmb;
  double qgb;
  double Cgg, Cgd;
  double Cgb, Cdg, Cdd, Cds;
  double Csg, Csd, Css;
  double Csb, Cbg, Cbd;
  double Cbb;

  // Capacitance variables
  //  In general  gcggb -->  CAPcggb.
  //              gcgdb -->  CAPcgdb.  etc.
  double CAPcggb;
  double CAPcgdb;
  double CAPcgsb;
  double CAPcbgb;
  double CAPcbdb;
  double CAPcbsb;
  double CAPcdgb;
  double CAPcdgmb;
  double CAPcddb;
  double CAPcdbdb;
  double CAPcdsb;
  double CAPcsgb;
  double CAPcsdb;
  double CAPcssb;
  double CAPcgmdb;
  double CAPcgmsb;
  double CAPcgmgmb;
  double CAPcbgmb;
  double CAPcsbsb;
  double CAPcqgb;
  double CAPcqdb;
  double CAPcqsb;

  double CAPcgmbb;
  double CAPcsgmb;
  double CAPcgbb;
  double CAPcdbb;
  double CAPcsbb;
  double CAPcbbb;
  double CAPcqbb;

  double Qeqqd_Jdxp;
  double Qeqqb_Jdxp;
  double Qeqqg_Jdxp;

  double Qeqqgmid_Jdxp;
  double Qeqqjs_Jdxp;
  double Qeqqjd_Jdxp;
  double Qqcheq_Jdxp;
  // end of new-DAE stuff.

  // gate resistor model currents (not needed in spice version)
  //
  // rgateMod==0   no gate resistor.
  // rgateMod==1   linear gate resistor
  // rgateMod==2   nonlinear gate resistor
  // rgateMod==3   2 gate resistors, in series.
  //
  double Igate;     // used by rgateMod= 1,2, and 3
  double IgateMid;  // used by rgateMod= 3 only

  // gate model voltage drops:
  double Vgegp;
  double Vgegm;
  double Vgmgp;

  // state variables, intrinsic capacitors
  double qb;
  double qg;
  double qd;
  double qgmid;

  // state variables, parasitic capacitors
  double qbs;
  double qbd;

  // state variables, cheq
  double qcheq;
  double cqcheq;
  double cqcheq_Jdxp;

  // state variables, cdump  There is no "qcdump" variable, just cqcdump
  double cqcdump;

  // this is a state variable in SPICE, but we use it only as an
  // instance variable --- it is always type*Qtotal
  double qdef;

  // Indices into the state vector:

  // state variables, voltage drops
  int li_store_vbd;
  int li_store_vbs;
  int li_store_vgs;
  int li_store_vds;
  int li_store_vges;
  int li_store_vgms;
  int li_store_vdes;
  int li_store_vses;
  int li_store_vdbs;
  int li_store_vsbs;
  int li_store_vdbd;
  int li_store_vged;
  int li_store_vgmd;
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

  double Vdsat;
  double Vth;

  // state variables, intrinsic capacitors
  int li_state_qb;
  int li_state_qg;
  int li_state_qd;
  int li_state_qgmid;

  // state variables, parasitic capacitors
  int li_state_qbs;
  int li_state_qbd;

  // state variables, cheq
  int li_state_qcheq;

  // state variables, cdump
  int li_state_qcdump;

  // state variable, qdef
  int li_state_qdef;
  
  // branch variables for lead current and power
  int li_branch_dev_id;
  int li_branch_dev_ig;
  int li_branch_dev_is;
  int li_branch_dev_ib;

  ////////////////////////////////////////////////////////////////////
  //  Local variable indices
  int li_Drain;      // dNode;
  int li_GateExt;    // gNodeExt;
  int li_Source;     // sNode;
  int li_Body;       // bNode;
  int li_DrainPrime; // dNodePrime;
  int li_GatePrime;  // gNodePrime;
  int li_GateMid;    // gNodeMid;
  int li_SourcePrime;// sNodePrime;
  int li_BodyPrime;  // bNodePrime;
  int li_DrainBody;  // dbNode;
  int li_SourceBody; // sbNode;
  int li_Charge;     // qNode;

  // local indies
  int li_Ibs;
  int li_Ids;
  int li_Igs;

  ////////////////////////////////////////////////////////////////////
  //  Jacobian matrix offsets:
  int GEge, GEgp, GEdp, GEsp, GEbp, GEgm, GEigs;
  int GPge, GPgp, GPdp, GPsp, GPbp, GPq, GPgm;
  int GMge, GMgp, GMdp, GMsp, GMbp, GMgm;
  int DPgm, DPdp,  DPd, DPgp, DPsp, DPbp, DPdb, DPq;
  int   Dd,  Dgp,  Ddp,  Dsp,  Dbp, Dids;
  int SPgm, SPdp, SPgp, SPsp, SPs, SPbp, SPsb, SPq;
  int   Ss,  Sdp,  Sgp,  Ssp, Sbp, Sibs, Sids, Sigs;
  int BPgm, BPdp, BPgp, BPsp, BPb, BPbp, BPdb, BPsb;
  int DBdp, DBdb, DBbp, DBb;
  int SBsp, SBbp,  SBb, SBsb;
  int  Bdb,  Bbp,  Bsb, Bb, Bibs;
  int   Qq,  Qgp,  Qdp, Qsp, Qbp;

  // icVBS
  int IBSb;
  int IBSs;
  int IBSibs;

  // icVDS
  int IDSd;
  int IDSs;
  int IDSids;

  // icVGS
  int IGSg;
  int IGSs;
  int IGSigs;


#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Jacobian Matrix Pointers

  double * f_DdPtr  ;	  double * q_DdPtr  ;
  double * f_DdpPtr ;	  double * q_DdpPtr ;
  double * f_DspPtr ;	  double * q_DspPtr ;
  double * f_DgpPtr ;	  double * q_DgpPtr ;
  double * f_DbpPtr ;	  double * q_DbpPtr ;
  double * f_DidsPtr;

  double * f_GEgePtr ;	  double * q_GEgePtr ;
  double * f_GEdpPtr ;	  double * q_GEdpPtr ;
  double * f_GEspPtr ;	  double * q_GEspPtr ;
  double * f_GEgpPtr ;	  double * q_GEgpPtr ;
  double * f_GEgmPtr ;	  double * q_GEgmPtr ;
  double * f_GEbpPtr ;	  double * q_GEbpPtr ;
  double * f_GEigsPtr;

  double * f_SsPtr  ;	  double * q_SsPtr  ;
  double * f_SdpPtr ;	  double * q_SdpPtr ;
  double * f_SspPtr ;	  double * q_SspPtr ;
  double * f_SgpPtr ;	  double * q_SgpPtr ;
  double * f_SbpPtr ;	  double * q_SbpPtr ;
  double * f_SibsPtr;
  double * f_SidsPtr;
  double * f_SigsPtr;

  double * f_BbPtr  ;	  double * q_BbPtr  ;
  double * f_BbpPtr ;	  double * q_BbpPtr ;
  double * f_BsbPtr ;	  double * q_BsbPtr ;
  double * f_BdbPtr ;	  double * q_BdbPtr ;
  double * f_BibsPtr;

  double * f_DPdPtr  ;	  double * q_DPdPtr  ;
  double * f_DPdpPtr ;	  double * q_DPdpPtr ;
  double * f_DPspPtr ;	  double * q_DPspPtr ;
  double * f_DPgpPtr ;	  double * q_DPgpPtr ;
  double * f_DPgmPtr ;	  double * q_DPgmPtr ;
  double * f_DPbpPtr ;	  double * q_DPbpPtr ;
  double * f_DPdbPtr ;	  double * q_DPdbPtr ;


  double * f_DPqPtr ;	    double * q_DPqPtr ;


  double * f_SPsPtr  ;	  double * q_SPsPtr  ;
  double * f_SPdpPtr ;	  double * q_SPdpPtr ;
  double * f_SPspPtr ;	  double * q_SPspPtr ;
  double * f_SPgpPtr ;	  double * q_SPgpPtr ;
  double * f_SPgmPtr ;	  double * q_SPgmPtr ;
  double * f_SPbpPtr ;	  double * q_SPbpPtr ;
  double * f_SPsbPtr ;	  double * q_SPsbPtr ;

  double * f_SPqPtr ;	    double * q_SPqPtr ;

  double * f_GPgePtr;	  double * q_GPgePtr;
  double * f_GPdpPtr;	  double * q_GPdpPtr;
  double * f_GPspPtr;	  double * q_GPspPtr;
  double * f_GPgpPtr;	  double * q_GPgpPtr;
  double * f_GPgmPtr;	  double * q_GPgmPtr;
  double * f_GPbpPtr;	  double * q_GPbpPtr;

  double * f_GPqPtr;	    double * q_GPqPtr;

  double * f_GMgePtr;	  double * q_GMgePtr;
  double * f_GMdpPtr;	  double * q_GMdpPtr;
  double * f_GMspPtr;	  double * q_GMspPtr;
  double * f_GMgpPtr;	  double * q_GMgpPtr;
  double * f_GMgmPtr;	  double * q_GMgmPtr;
  double * f_GMbpPtr;	  double * q_GMbpPtr;

  double * f_BPbPtr ;	  double * q_BPbPtr ;
  double * f_BPdpPtr;	  double * q_BPdpPtr;
  double * f_BPspPtr;	  double * q_BPspPtr;
  double * f_BPgpPtr;	  double * q_BPgpPtr;
  double * f_BPgmPtr;	  double * q_BPgmPtr;
  double * f_BPbpPtr;	  double * q_BPbpPtr;
  double * f_BPsbPtr;	  double * q_BPsbPtr;
  double * f_BPdbPtr;	  double * q_BPdbPtr;

  double * f_SBbPtr ;	  double * q_SBbPtr ;
  double * f_SBspPtr;	  double * q_SBspPtr;
  double * f_SBbpPtr;	  double * q_SBbpPtr;
  double * f_SBsbPtr;	  double * q_SBsbPtr;

  double * f_DBbPtr ;	  double * q_DBbPtr ;
  double * f_DBdpPtr;	  double * q_DBdpPtr;
  double * f_DBbpPtr;	  double * q_DBbpPtr;
  double * f_DBdbPtr;	  double * q_DBdbPtr;

  double * f_QdpPtr;	    double * q_QdpPtr;
  double * f_QspPtr;	    double * q_QspPtr;
  double * f_QgpPtr;	    double * q_QgpPtr;
  double * f_QbpPtr;	    double * q_QbpPtr;
  double * f_QqPtr ;	    double * q_QqPtr ;

  double * f_IBSbPtr;
  double * f_IBSsPtr;
  double * f_IBSibsPtr;
  double * f_IDSdPtr;
  double * f_IDSsPtr;
  double * f_IDSidsPtr;
  double * f_IGSgPtr;
  double * f_IGSsPtr;
  double * f_IGSigsPtr;
#endif

  // flag for updateTemperature call.  Needed for .STEP temperature sweeps.
  bool updateTemperatureCalled_;

  // Jacobian stamp related structures:
  std::vector< std::vector<int> > jacStamp;
  std::vector<int> jacMap;
  std::vector< std::vector<int> > jacMap2;

  int blockHomotopyID; // For homotopy
  double randomPerturb; // For homotopy

  double ceqdrn_Jdxp, ceqbd_Jdxp, ceqbs_Jdxp;
  double Istoteq_Jdxp, Idtoteq_Jdxp;
  double Ibtoteq_Jdxp, Igtoteq_Jdxp;
  double ceqgcrg_Jdxp, ceqgstot_Jdxp, ceqgdtot_Jdxp;
  double ceqjs_Jdxp, ceqjd_Jdxp;
  double T0;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/06
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class SizeDependParam;
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
  int modType;
  int dtype;

  int    mobMod;
  int    cvchargeMod;
  int    capMod;
  int    dioMod;
  int    trnqsMod;
  int    acnqsMod;
  int    fnoiMod;
  int    tnoiMod;
  int    rdsMod;
  int    rbodyMod;
  int    rgateMod;
  int    perMod;
  int    geoMod;
  int    mtrlMod;
  int    igcMod;
  int    igbMod;
  int    tempMod;
  int    binUnit;
  int    paramChk;
  std::string version;
  double eot;
  double vddeot;
  double ados;
  double bdos;
  double toxe;
  double toxp;
  double toxm;
  double dtox;
  double epsrox;
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
  double phig;
  double epsrgate;
  double easub;
  double epsrsub;
  double ni0sub;
  double bg0sub;
  double tbgasub;
  double tbgbsub;
  double ndep;
  double nsd;
  double phin;
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
  double dvtp0;
  double dvtp1;
  double lpe0;
  double lpeb;
  double dvt0;
  double dvt1;
  double dvt2;
  double dvt0w;
  double dvt1w;
  double dvt2w;
  double drout;
  double dsub;
  double vth0;
  double eu;
  double ua;
  double ua1;
  double ub;
  double ub1;
  double uc;
  double uc1;
  double ud;
  double ud1;
  double up;
  double lp;
  double u0;
  double ute;
  double voff;
  double tvoff;
  double minv;
  double minvcv;
  double voffl;
  double voffcvl;
  double delta;
  double rdsw;
  double rdswmin;
  double rdwmin;
  double rswmin;
  double rsw;
  double rdw;
  double prwg;
  double prwb;
  double prt;
  double eta0;
  double etab;
  double pclm;
  double pdibl1;
  double pdibl2;
  double pdiblb;
  double fprout;
  double pdits;
  double pditsd;
  double pditsl;
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
  double agidl;
  double bgidl;
  double cgidl;
  double egidl;
  double agisl;
  double bgisl;
  double cgisl;
  double egisl;
  double aigc;
  double bigc;
  double cigc;
  double aigsd;
  double bigsd;
  double cigsd;
  double aigs;
  double bigs;
  double cigs;
  double aigd;
  double bigd;
  double cigd;
  double aigbacc;
  double bigbacc;
  double cigbacc;
  double aigbinv;
  double bigbinv;
  double cigbinv;
  double nigc;
  double nigbacc;
  double nigbinv;
  double ntox;
  double eigbinv;
  double pigcd;
  double poxedge;
  double toxref;
  double ijthdfwd;
  double ijthsfwd;
  double ijthdrev;
  double ijthsrev;
  double xjbvd;
  double xjbvs;
  double bvd;
  double bvs;

  double jtss;
  double jtsd;
  double jtssws;
  double jtsswd;
  double jtsswgs;
  double jtsswgd;
  double njts;
  double njtssw;
  double njtsswg;
  double njtsd;
  double njtsswd;
  double njtsswgd;
  double xtss;
  double xtsd;
  double xtssws;
  double xtsswd;
  double xtsswgs;
  double xtsswgd;
  double tnjts;
  double tnjtssw;
  double tnjtsswg;
  double tnjtsd;
  double tnjtsswd;
  double tnjtsswgd;
  double vtss;
  double vtsd;
  double vtssws;
  double vtsswd;
  double vtsswgs;
  double vtsswgd;

  double xrcrg1;
  double xrcrg2;
  double lambda;
  double vtl;
  double lc;
  double xn;
  double vfbsdoff;  // S/D flatband offset voltage
  double lintnoi;  // lint offset for noise calculation
  double tvfbsdoff;

  double vfb;
  double gbmin;
  double rbdb;
  double rbsb;
  double rbpb;
  double rbps;
  double rbpd;

  double rbps0;
  double rbpsl;
  double rbpsw;
  double rbpsnf;

  double rbpd0;
  double rbpdl;
  double rbpdw;
  double rbpdnf;

  double rbpbx0;
  double rbpbxl;
  double rbpbxw;
  double rbpbxnf;
  double rbpby0;
  double rbpbyl;
  double rbpbyw;
  double rbpbynf;

  double rbsbx0;
  double rbsby0;
  double rbdbx0;
  double rbdby0;

  double rbsdbxl;
  double rbsdbxw;
  double rbsdbxnf;
  double rbsdbyl;
  double rbsdbyw;
  double rbsdbynf;

  double tnoia;
  double tnoib;
  double rnoia;
  double rnoib;
  double ntnoi;

  // CV model and Parasitics
  double cgsl;
  double cgdl;
  double ckappas;
  double ckappad;
  double cf;
  double vfbcv;
  double clc;
  double cle;
  double dwc;
  double dlc;
  double xw;
  double xl;
  double dlcig;
  double dlcigd;
  double dwj;
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
  double dmcg;
  double dmci;
  double dmdg;
  double dmcgt;
  double xgw;
  double xgl;
  double rshg;
  double ngcon;

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
  double lndep;
  double lnsd;
  double lphin;
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
  double ldvtp0;
  double ldvtp1;
  double llpe0;
  double llpeb;
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
  double lud;
  double lud1;
  double lup;
  double llp;
  double lu0;
  double leu;
  double lute;
  double lvoff;
  double ltvoff;
  double lminv;
  double lminvcv;
  double ldelta;
  double lrdsw;
  double lrsw;
  double lrdw;
  double lprwg;
  double lprwb;
  double lprt;
  double leta0;
  double letab;
  double lpclm;
  double lpdibl1;
  double lpdibl2;
  double lpdiblb;
  double lfprout;
  double lpdits;
  double lpditsd;
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
  double lagidl;
  double lbgidl;
  double lcgidl;
  double legidl;
  double lagisl;
  double lbgisl;
  double lcgisl;
  double legisl;
  double laigc;
  double lbigc;
  double lcigc;
  double laigsd;
  double lbigsd;
  double lcigsd;
  double laigs;
  double lbigs;
  double lcigs;
  double laigd;
  double lbigd;
  double lcigd;
  double laigbacc;
  double lbigbacc;
  double lcigbacc;
  double laigbinv;
  double lbigbinv;
  double lcigbinv;
  double lnigc;
  double lnigbacc;
  double lnigbinv;
  double lntox;
  double leigbinv;
  double lpigcd;
  double lpoxedge;
  double lxrcrg1;
  double lxrcrg2;
  double llambda;
  double lvtl;
  double lxn;
  double lvfbsdoff;
  double ltvfbsdoff;

  // CV model
  double lcgsl;
  double lcgdl;
  double lckappas;
  double lckappad;
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
  double wndep;
  double wnsd;
  double wphin;
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
  double wdvtp0;
  double wdvtp1;
  double wlpe0;
  double wlpeb;
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
  double wud;
  double wud1;
  double wup;
  double wlp;
  double wu0;
  double weu;
  double wute;
  double wvoff;
  double wtvoff;
  double wminv;
  double wminvcv;
  double wdelta;
  double wrdsw;
  double wrsw;
  double wrdw;
  double wprwg;
  double wprwb;
  double wprt;
  double weta0;
  double wetab;
  double wpclm;
  double wpdibl1;
  double wpdibl2;
  double wpdiblb;
  double wfprout;
  double wpdits;
  double wpditsd;
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
  double wagidl;
  double wbgidl;
  double wcgidl;
  double wegidl;
  double wagisl;
  double wbgisl;
  double wcgisl;
  double wegisl;
  double waigc;
  double wbigc;
  double wcigc;
  double waigsd;
  double wbigsd;
  double wcigsd;
  double waigs;
  double wbigs;
  double wcigs;
  double waigd;
  double wbigd;
  double wcigd;
  double waigbacc;
  double wbigbacc;
  double wcigbacc;
  double waigbinv;
  double wbigbinv;
  double wcigbinv;
  double wnigc;
  double wnigbacc;
  double wnigbinv;
  double wntox;
  double weigbinv;
  double wpigcd;
  double wpoxedge;
  double wxrcrg1;
  double wxrcrg2;
  double wlambda;
  double wvtl;
  double wxn;
  double wvfbsdoff;
  double wtvfbsdoff;

  // CV model
  double wcgsl;
  double wcgdl;
  double wckappas;
  double wckappad;
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
  double pndep;
  double pnsd;
  double pphin;
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
  double pdvtp0;
  double pdvtp1;
  double plpe0;
  double plpeb;
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
  double pud;
  double pud1;
  double pup;
  double plp;
  double pu0;
  double peu;
  double pute;
  double pvoff;
  double ptvoff;
  double pminv;
  double pminvcv;
  double pdelta;
  double prdsw;
  double prsw;
  double prdw;
  double pprwg;
  double pprwb;
  double pprt;
  double peta0;
  double petab;
  double ppclm;
  double ppdibl1;
  double ppdibl2;
  double ppdiblb;
  double pfprout;
  double ppdits;
  double ppditsd;
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
  double pagidl;
  double pbgidl;
  double pcgidl;
  double pegidl;
  double pagisl;
  double pbgisl;
  double pcgisl;
  double pegisl;
  double paigc;
  double pbigc;
  double pcigc;
  double paigsd;
  double pbigsd;
  double pcigsd;
  double paigs;
  double pbigs;
  double pcigs;
  double paigd;
  double pbigd;
  double pcigd;
  double paigbacc;
  double pbigbacc;
  double pcigbacc;
  double paigbinv;
  double pbigbinv;
  double pcigbinv;
  double pnigc;
  double pnigbacc;
  double pnigbinv;
  double pntox;
  double peigbinv;
  double ppigcd;
  double ppoxedge;
  double pxrcrg1;
  double pxrcrg2;
  double plambda;
  double pvtl;
  double pxn;
  double pvfbsdoff;
  double ptvfbsdoff;

  // CV model
  double pcgsl;
  double pcgdl;
  double pckappas;
  double pckappad;
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
  double SjctSatCurDensity;
  double DjctSatCurDensity;
  double SjctSidewallSatCurDensity;
  double DjctSidewallSatCurDensity;
  double SjctGateSidewallSatCurDensity;
  double DjctGateSidewallSatCurDensity;
  double SbulkJctPotential;
  double DbulkJctPotential;
  double SbulkJctBotGradingCoeff;
  double DbulkJctBotGradingCoeff;
  double SbulkJctSideGradingCoeff;
  double DbulkJctSideGradingCoeff;
  double SbulkJctGateSideGradingCoeff;
  double DbulkJctGateSideGradingCoeff;
  double SsidewallJctPotential;
  double DsidewallJctPotential;
  double SGatesidewallJctPotential;
  double DGatesidewallJctPotential;
  double SunitAreaJctCap;
  double DunitAreaJctCap;
  double SunitLengthSidewallJctCap;
  double DunitLengthSidewallJctCap;
  double SunitLengthGateSidewallJctCap;
  double DunitLengthGateSidewallJctCap;
  double SjctEmissionCoeff;
  double DjctEmissionCoeff;
  double SjctTempExponent;
  double DjctTempExponent;
  double njtsstemp;
  double njtsswstemp;
  double njtsswgstemp;
  double njtsdtemp;
  double njtsswdtemp;
  double njtsswgdtemp;

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

  // added for stress effect
  double saref;
  double sbref;
  double wlod;
  double ku0;
  double kvsat;
  double kvth0;
  double tku0;
  double llodku0;
  double wlodku0;
  double llodvth;
  double wlodvth;
  double lku0;
  double wku0;
  double pku0;
  double lkvth0;
  double wkvth0;
  double pkvth0;
  double stk2;
  double lodk2;
  double steta0;
  double lodeta0;

  double web;
  double wec;
  double kvth0we;
  double k2we;
  double ku0we;
  double scref;
  double wpemod;
  double lkvth0we;
  double lk2we;
  double lku0we;
  double wkvth0we;
  double wk2we;
  double wku0we;
  double pkvth0we;
  double pk2we;
  double pku0we;

  // Pre-calculated constants
  // move to size-dependent param
  double Eg0;
  double vtm;
  double vtm0;
  double coxe;
  double coxp;
  double cof1;
  double cof2;
  double cof3;
  double cof4;
  double vcrit;
  double factor1;
  double PhiBS;
  double PhiBSWS;
  double PhiBSWGS;
  double SjctTempSatCurDensity;
  double SjctSidewallTempSatCurDensity;
  double SjctGateSidewallTempSatCurDensity;
  double PhiBD;
  double PhiBSWD;
  double PhiBSWGD;
  double DjctTempSatCurDensity;
  double DjctSidewallTempSatCurDensity;
  double DjctGateSidewallTempSatCurDensity;
  double SunitAreaTempJctCap;
  double DunitAreaTempJctCap;
  double SunitLengthSidewallTempJctCap;
  double DunitLengthSidewallTempJctCap;
  double SunitLengthGateSidewallTempJctCap;
  double DunitLengthGateSidewallTempJctCap;

  double oxideTrapDensityA;
  double oxideTrapDensityB;
  double oxideTrapDensityC;
  double em;
  double ef;
  double af;
  double kf;

  double ni;
  double Vtm0;

  // given variables:
  bool vtlGiven;
  bool ndepGiven;
  bool gamma1Given;
  bool k1Given;
  bool k2Given;
  bool nsubGiven;
  bool phigGiven;
  bool xtGiven;
  bool vbxGiven;
  bool gamma2Given;
  bool vfbGiven;
  bool vth0Given;
  bool rbps0Given;
  bool rbpd0Given;
  bool rbsbx0Given;
  bool rbsby0Given;
  bool rbdbx0Given;
  bool rbdby0Given;
  bool lambdaGiven;
  bool pigcdGiven;

  bool toxeGiven;
  bool toxpGiven;
  bool dtoxGiven;
  bool cgdoGiven;
  bool dlcGiven;
  bool cgsoGiven;
  bool cgboGiven;

  std::list<SizeDependParam*> sizeDependParamList;
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
  friend class SizeDependParam;
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

  // load functions, residual:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);

  // load functions, Jacobian:
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace MOSFET_B4
} // namespace Device
} // namespace Xyce

#endif


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
// Creator        : Dave Shirley
//
// Creation Date  : 05/20/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MOSFET_B3SOI_h
#define Xyce_N_DEV_MOSFET_B3SOI_h

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
namespace MOSFET_B3SOI {

// ---------- Forward Declarations -------
class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, MOSFET1::Traits>
{
  static const char *name() {return "BSIM3 SOI";}
  static const char *deviceTypeName() {return "M level 10";}
  static int numNodes() {return 4;}
  static int numOptionalNodes() {return 3;}
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
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class SizeDependParam
{
  friend class Model;
  friend class Instance;
  friend class Master;

private:
  double Width;
  double Length;
  double Rth0;
  double Cth0;

  double ntrecf;
  double ndif;
  double betaGB1;
  double betaGB2;
  double alphaGB1;
  double alphaGB2;
  double ahli;
  double aely;
  double vabjt;
  double lbjt0;
  double nbjt;
  double vtun0;
  double vrec0;
  double istun;
  double isrec;
  double isdif;
  double isbjt;
  double nrecr0;
  double nrecf0;
  double ndiode;
  double ntun;
  double ngidl;
  double bgidl;
  double agidl;
  double siid;
  double sii2;
  double sii1;
  double sii0;
  double esatii;
  double lii;
  double vdsatii0;
  double beta2;
  double beta1;
  double fbjtii;
  double ketas;
  double kb1;
  double k1w1;
  double k1w2;
  double vsdfb;
  double vfbb;
  double jtun;
  double xtun;
  double jrec;
  double jdif;
  double jbjt;
  double ahli0;
  double xrec;
  double xdif;
  double xbjt;
  double cgeo;
  double wdiosCV;
  double wdiodCV;
  double oxideRatio;
  double rbody;
  double cth;
  double rth;
  double rds0denom;
  double uctemp;
  double ubtemp;
  double uatemp;
  double delvt;
  double vsdth;
  double xrcrg2;
  double xrcrg1;
  double poxedge;
  double pigcd;
  double cigsd;
  double bigsd;
  double aigsd;
  double cigc;
  double bigc;
  double aigc;
  double nigc;
  double dt3;
  double dt2;
  double st3;
  double st2;
  double sdt1;
  double k1eff;
  double Bechvb;
  double ToxRatio;
  double Aechvb;
  double BechvbEdge;
  double ToxRatioEdge;
  double dlcig;
  double AechvbEdge;
  double qsi;
  double vearly;
  double lratiodif;
  double lratio;
  double arfabjt;
  double leffCVbg;
  double ntrecr;
  double leffCVb;
  double wdios;
  double wdiod;
  double vfbsd;

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
  double pvag;
  double wr;
  double dwg;
  double dwb;
  double b0;
  double b1;
  double alpha0;
  double beta0;

  // CV model
  double cgsl;
  double cgdl;
  double ckappa;
  double cf;
  double clc;
  double cle;
  double noff;
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

  double u0temp;
  double vsattemp;
  double sqrtPhi;
  double phis3;
  double Xdep0;
  double sqrtXdep0;
  double theta0vb0;
  double thetaRout;

  double cdep0;
  double vfbzb;
  double ldeb;

  // ERK.  I added this to make temperature sweeps work.
  double referenceTemperature;
};


//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Dave Shirley
// Creation Date : 05/20/04
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
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
  void registerStoreLIDs( const std::vector<int> & stoLIDVecRef);
  virtual void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);
  
  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();

  bool updateTemperature(const double & temp_tmp);
  bool updateIntermediateVars ();
  bool updatePrimaryState ();

  double Eval1ovFNoise (double vds, double freq);

  int getNumNoiseSources () const;
  void setupNoiseSources (Xyce::Analysis::NoiseData & noiseData);
  void getNoiseSources (Xyce::Analysis::NoiseData & noiseData);

  bool loadMatrix (Linear::Matrix & JMat);
  bool checkModel ();
  double B3SOIlimit(double vnew, double vold, double limit, int *check);

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void setupJacStamp ();
  void debugOutputModelParams();

  bool auxChargeCalculations ();
  bool setupCapacitors_newDAE ();

  void setupPointers ();

  bool setupCapacitors_oldDAE ();

  bool setIC ();

  // Beginning of 3f5 stuff:
  // attributes:
public:
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:
  static std::vector< std::vector< std::vector<int> > > jacStamp_v;
  static std::vector< std::vector<int> > jacMap_v;
  static std::vector< std::vector< std::vector<int> > > jacMap2_v;


  Model &       model_;         //< Owning model

  int dNode;
  int gNode;
  int sNode;
  int eNode;
  int pNode;
  bool pNodeMappedToB;
  int bNode;
  int tNode;
  int dNodePrime;
  int sNodePrime;
  int gNodePrime;
  int gNodeMid;
  int P_index;
  int B_index;
  int T_index;
  int jacID;

  double ueff;
  double thetavth;
  double von;
  double vdsat;
  double cgdo;
  double cgso;

  double l;
  double w;
  double numberParallel;
  double drainArea;
  double sourceArea;
  double drainSquares;
  double sourceSquares;
  double drainPerimeter;
  double sourcePerimeter;

  //3f5 bsimsoi varaibles used in creating jacobian:
  double drainConductance;
  double FwdSum;
  double gbbb;
  double gbbdp;
  double gbbe;
  double gbbg;
  double gbbsp;
  double gbbT;
  double gcbdb;
  double gcbeb;
  double gcbgb;
  double gcbsb;
  double gcbT;
  double gcddb;
  double gcdeb;
  double gcdgb;
  double gcdgmb;
  double gcdsb;
  double gcdT;
  double gcedb;
  double gceeb;
  double gcegb;
  double gcegmb;
  double gcesb;
  double gceT;
  double gcgbb;
  double gcgdb;
  double gcgeb;
  double gcggb;
  double gcgmdb;
  double gcgmeb;
  double gcgmgmb;
  double gcgmsb;
  double gcgsb;
  double gcgT;
  double gcrg;
  double gcrgb;
  double gcrgd;
  double gcrgg;
  double gcrgs;
  double gcrg_jac;
  double gcrgb_jac;
  double gcrgd_jac;
  double gcrgg_jac;
  double gcrgs_jac;
  double gcsdb;
  double gcseb;
  double gcsgb;
  double gcsgmb;
  double gcssb;
  double gcsT;
  double gcTt;
  double gddpb;
  double gddpdp;
  double gddpe;
  double gddpg;
  double gddpsp;
  double gddpT;
  double gds;
  double geltd;
  double gigb;
  double gigd;
  double gige;
  double gigg;
  double gigs;
  double gigT;
  double gigb_jac;
  double gigd_jac;
  double gige_jac;
  double gigg_jac;
  double gigs_jac;
  double gigT_jac;
  double gIdtotb;
  double gIgtotb;
  double gIgtotd;
  double gIgtotg;
  double gIgtots;
  double gIstotb;
  double gIstotd;
  double gIstotg;
  double gIstots;
  double gppb;
  double gppp;
  double gsspb;
  double gsspdp;
  double gsspe;
  double gsspg;
  double gsspsp;
  double gsspT;
  double gTtb;
  double gTtdp;
  double gTte;
  double gTtg;
  double gTtsp;
  double gTtt;
  double Gm;
  double Gmbs;
  double Gme;
  double Gmin;
  double GmT;
  double RevSum;
  double sourceConductance;
  double gIdtotg;
  double gIdtotd;
  double gIdtots;

  // new variable from updateTempreature
  double rbodyext;
  double csesw;
  double dt4;
  double st4;
  double cdmin;
  double cdbox;
  double csbox;
  double csmin;
  double grgeltd;
  double phi;
  double cdesw;

  int mode;
  int bjtoff;
  int debugMod;
  bool OFF;
  double rth0;
  double cth0;
  double bodySquares;
  double frbody;
  int soiMod;
  double nbc;
  double nseg;
  double pdbcp;
  double psbcp;
  double agbcp;
  double aebcp;
  double vbsusr;
  int tnodeout;
  int rgateMod;

  double cdrain;
  double gIgsg;
  double gIgss;
  double gIgcdg;
  double gIgcds;
  double gIgcdd;
  double gIgcdb;
  double Igs;
  double Igcd;
  double gIgdg;
  double gIgcsg;
  double gIgdd;
  double gIgcss;
  double gIgcsd;
  double gIgcsb;
  double Igd;
  double Igcs;

  // OP point
  double qinv;
  double cb;
  double cd;
  double cbd;
  double gm;
  double gmbs;

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

  double qbulk;

  // initial condition variables
  double icVDS;
  double icVGS;
  double icVBS;
  double icVES;
  double icVPS;

  bool icVDSGiven;
  bool icVGSGiven;
  bool icVBSGiven;
  bool icVESGiven;
  bool icVPSGiven;

  SizeDependParam  * paramPtr;

  // end of original 3f5 stuff

  // Variables from the 3f5 b3ld function, which  were local to that
  // function but are more appropriate as instance variables:
  bool ChargeComputationNeeded;

  bool selfheat;
  int bodyMod;
  int floating;

  double dxpart;
  double sxpart;

  double cdreq;
  double ceqbd;
  double ceqbs;

  double qgdo;
  double qgso;
  double qgd;
  double qgs;
  double qge;
  double qgme;
  double qgate;
  double qbody;
  double qdrn;
  double qsub;
  double qsrc;

  // added for soi:
  double ceqbody;
  double ceqgate;
  double ceqgcrg;
  double ceqqe;
  double ceqqgmid;
  double ceqbodcon;
  double ceqth;
  double ceqqth;
  double Igtoteq;
  double Idtoteq;
  double Istoteq;


  double dVgst_dVg;
  double dVgst_dVb;
  double dVgs_eff_dVg;

  double dDeltaPhi_dVg;
  double dDeltaPhi_dVd;
  double dDeltaPhi_dVb;

  double cqdrn;
  double cqgate;
  double cqsub;
  double cqbody;
  double cqtemp;

  // Variables which were model variables in the 3f5 version of the B3SOI,
  // but are temperature dependent and should be instance variables (as
  // temperature is an instance-level quantity.
  double vtm;

  // end of 3f5 stuff

  double temp;

  // solution variables, and intermediate quantities.
  //double mode;
  double Vd;
  double Vg;
  double Vs;
  double Ve;
  double Vp;
  double Vb;
  double Vsp;
  double Vdp;
  double Vgp;
  double Vgm;

  double Qtotal;

  // resistor currents:
  double Idrain;
  double Isource;

  double Igate;
  double IgateMid;
  double Ith;

  double vgb;
  double vgd;
  double ceqqd;
  double ceqqb;
  double ceqqg;

  double ceqqg_Jdxp;
  double ceqqb_Jdxp;
  double ceqqd_Jdxp;
  double ceqqe_Jdxp;
  double ceqqth_Jdxp;
  double ceqgcrg_Jdxp;
  double ceqqgmid_Jdxp;
  double cth_Jdxp;
  double cjs_Jdxp;
  double cjd_Jdxp;
  double cdreq_Jdxp;
  double cbody_Jdxp;
  double cgate_Jdxp;
  double Idtoteq_Jdxp;
  double Istoteq_Jdxp;
  double Igtoteq_Jdxp;

  double cbodcon_Jdxp;
  double ceqbody_Jdxp;
  double ceqgate_Jdxp;
  double ceqbs_Jdxp;
  double ceqbd_Jdxp;
  double ceqbodcon_Jdxp;

  double Idrain_Jdxp;
  double Isource_Jdxp;
  double Igate_Jdxp;
  double IgateMid_Jdxp;

  double ceqth_Jdxp;

  // state variables, voltage drops
  double vbd;
  double vbs;
  double vps;
  double vpd;
  double ved;
  double veb;
  double ves;
  double vgs;
  double vge;
  double vds;
  double vged;
  double vgmd;
  double vgme;
  double vgmb;

  double vg;
  double vd;
  double vs;
  double vp;
  double ve;
  double deltemp;
  double delTemp;
  double TempRatioMinus1;

  double vges;
  double vgms;
  double vgge;
  double vggm;

  // "original" versions of various voltage drop variables:
  // original refers to the beginning of the newton iterations,
  // before any limits are imposed on the change in voltage drop.
  double vbd_orig;
  double vbs_orig;
  double vps_orig;
  double vpd_orig;
  double ves_orig;
  double ved_orig;
  double vgs_orig;
  double vds_orig;

  // and these are "mode-aware" versions of voltage drops.  If vds is
  // negative then drain and source are reversed and one needs to
  // use them differently.
  double Vds, Vgs, Vbs;
  double Vbd, Ves, Vps;
  double Vds_orig, Vgs_orig, Vbs_orig;
  double Vbd_orig, Ves_orig, Vps_orig;

  double delTemp_orig;

  double vges_orig;
  double vgms_orig;

  // "orig" variables that are node-based, rather than junction.
  double Vd_orig;
  double Vg_orig;
  double Vs_orig;
  double Ve_orig;
  double Vb_orig;
  double Vp_orig;
  double Vsp_orig;
  double Vdp_orig;
  double Vgp_orig;
  double Vgm_orig;

  // These two are added to allow loadRHS to compile
  double vgd_orig;

  //PMC: global variables needed in intermediateVars
  double Vgsteff;
  double Vdseff;
  double ni;
  double Abulk;
  double vbseff;
  double nstar;
  double rds;
  double AbovVgst2Vtm;
  double ids;
  double igidl;
  double ic;
  double ig;
  double itun;
  double ibs;
  double ibd;
  double iii;
  double ibp;
  double gbpbs;
  double gbpps;
  double gbpT;
  double cbodcon;

  double gme;
  double gmT;
  double gtempg;
  double gtempb;
  double gtempe;
  double gtempT;
  double gtempd;
  double cth;
  double cjs;
  double cjd;
  double cbody;
  double cgate;
  double gjdb;
  double gjdd;
  double gjdg;
  double gjde;
  double gjdT;
  double gjsb;
  double gjsd;
  double gjsg;
  double gjsT;
  double gbes;
  double gbps;
  double gbT;

  double cgT;
  double cbT;
  double ceT;
  double cdT;
  double cbeb;
  double ceeb;
  double cdeb;
  double qse;
  double qde;
  double qbf;
  double qjs;
  double qjd;
  double cbb;
  double cbg;
  double gcse;
  double gcde;

  // state variables, intrinsic capacitors
  double qb;
  double qg;
  double qd;
  double qe;
  double qgmid;
  double qth;

  double wdiosCV_NoSwap;
  double wdiodCV_NoSwap;

  double CAPcgmgmb;
  double CAPcgmdb;
  double CAPcgmsb;
  double CAPcgmeb;
  double CAPcdgmb;
  double CAPcsgmb;
  double CAPcegmb;
  double CAPcggb;
  double CAPcgdb;
  double CAPcgsb;
  double CAPcgeb;
  double CAPcgbb;
  double CAPcdgb;
  double CAPcegb;
  double CAPcsgb;
  double CAPcbgb;
  double CAPcddb;
  double CAPcdsb;
  double CAPcdeb;
  double CAPcdT;
  double CAPcsdb;
  double CAPcssb;
  double CAPcseb;
  double CAPcsT;
  double CAPcgT;
  double CAPcbdb;
  double CAPcbsb;
  double CAPcbeb;
  double CAPcbT;
  double CAPcedb;
  double CAPcesb;
  double CAPceeb;
  double CAPceT;
  double CAPcTt;

  double Qeqqg;
  double Qeqqb;
  double Qeqqd;
  double Qeqqe;
  double Qeqqth;
  double Qeqqgmid;

  double Qeqqg_Jdxp;
  double Qeqqb_Jdxp;
  double Qeqqd_Jdxp;
  double Qeqqe_Jdxp;
  double Qeqqth_Jdxp;
  double Qeqqgmid_Jdxp;

  // Indices into the state vector:

  // state variables, voltage drops
  int li_store_vbd;
  int li_store_vbs;
  int li_store_vgs;
  int li_store_vds;
  int li_store_ves;
  int li_store_vps;

  int li_store_vg;
  int li_store_vd;
  int li_store_vs;
  int li_store_vp;
  int li_store_ve;
  int li_store_vgp;
  int li_store_vgm;
  int li_store_deltemp;

  int li_store_vges;
  int li_store_vgms;

  // lead current vars if needed
  int li_branch_dev_id;
  int li_branch_dev_ig;
  int li_branch_dev_is;
  int li_branch_dev_ie;
  int li_branch_dev_ib;

  // state variables, intrinsic capacitors
  int li_state_qb;
  int li_state_qg;
  int li_state_qd;
  int li_state_qe;
  int li_state_qgmid;
  int li_state_qth;

  ////////////////////////////////////////////////////////////////////
  //  Local variable indices
  int li_Drain;
  int li_Gate;
  int li_Source;
  int li_Substrate;
  int li_ExtBody;
  int li_Body;
  int li_Temperature;
  int li_DrainPrime;
  int li_SourcePrime;
  int li_GatePrime;
  int li_GateMid;

  // local indies for currents on voltages sources used for initial conditions
  int li_Ids;
  int li_Igs;
  int li_Ibs;
  int li_Ies;
  int li_Ips;


  ////////////////////////////////////////////////////////////////////
  // Offset variables corresponding to the above declared indices.

  // Jacobian Matrix Offsets:

  //  drain row:
  int ADrainEquDrainNodeOffset;
  int ADrainEquDrainPrimeNodeOffset;
  int ADrainEquIdsOffset;

  //  gate row:
  int AGateEquGateNodeOffset;
  int AGateEquBodyNodeOffset;
  int AGateEquDrainPrimeNodeOffset;
  int AGateEquSourcePrimeNodeOffset;
  int AGateEquGatePrimeNodeOffset;
  int AGateEquGateMidNodeOffset;
  int AGateEquIgsOffset;

  //  source row:
  int ASourceEquSourceNodeOffset;
  int ASourceEquSourcePrimeNodeOffset;
  int ASourceEquIdsOffset;
  int ASourceEquIgsOffset;
  int ASourceEquIbsOffset;
  int ASourceEquIesOffset;
  int ASourceEquIpsOffset;

  //  substrate row:
  int ASubstrateEquSubstrateNodeOffset;
  int ASubstrateEquBodyNodeOffset;
  int ASubstrateEquTemperatureNodeOffset;
  int ASubstrateEquDrainPrimeNodeOffset;
  int ASubstrateEquSourcePrimeNodeOffset;
  int ASubstrateEquGatePrimeNodeOffset;
  int ASubstrateEquGateMidNodeOffset;
  int ASubstrateEquIesOffset;

  // external body row:
  int AExtBodyEquExtBodyNodeOffset;
  int AExtBodyEquBodyNodeOffset;
  int AExtBodyEquIpsOffset;

  // body row:
  int ABodyEquSubstrateNodeOffset;
  int ABodyEquExtBodyNodeOffset;
  int ABodyEquBodyNodeOffset;
  int ABodyEquTemperatureNodeOffset;
  int ABodyEquDrainPrimeNodeOffset;
  int ABodyEquSourcePrimeNodeOffset;
  int ABodyEquGatePrimeNodeOffset;
  int ABodyEquIbsOffset;

  // temperature row:
  int ATemperatureEquSubstrateNodeOffset;
  int ATemperatureEquBodyNodeOffset;
  int ATemperatureEquTemperatureNodeOffset;
  int ATemperatureEquDrainPrimeNodeOffset;
  int ATemperatureEquSourcePrimeNodeOffset;
  int ATemperatureEquGatePrimeNodeOffset;

  // drain' row:
  int ADrainPrimeEquDrainNodeOffset;
  int ADrainPrimeEquSubstrateNodeOffset;
  int ADrainPrimeEquBodyNodeOffset;
  int ADrainPrimeEquTemperatureNodeOffset;
  int ADrainPrimeEquDrainPrimeNodeOffset;
  int ADrainPrimeEquSourcePrimeNodeOffset;
  int ADrainPrimeEquGatePrimeNodeOffset;
  int ADrainPrimeEquGateMidNodeOffset;

  // source' row:
  int ASourcePrimeEquSourceNodeOffset;
  int ASourcePrimeEquSubstrateNodeOffset;
  int ASourcePrimeEquBodyNodeOffset;
  int ASourcePrimeEquTemperatureNodeOffset;
  int ASourcePrimeEquDrainPrimeNodeOffset;
  int ASourcePrimeEquSourcePrimeNodeOffset;
  int ASourcePrimeEquGatePrimeNodeOffset;
  int ASourcePrimeEquGateMidNodeOffset;

  // gate' row:
  int AGatePrimeEquGateNodeOffset;
  int AGatePrimeEquSubstrateNodeOffset;
  int AGatePrimeEquBodyNodeOffset;
  int AGatePrimeEquTemperatureNodeOffset;
  int AGatePrimeEquDrainPrimeNodeOffset;
  int AGatePrimeEquSourcePrimeNodeOffset;
  int AGatePrimeEquGatePrimeNodeOffset;
  int AGatePrimeEquGateMidNodeOffset;

  // gate mid row:
  int AGateMidEquGateNodeOffset;
  int AGateMidEquSubstrateNodeOffset;
  int AGateMidEquBodyNodeOffset;
  int AGateMidEquDrainPrimeNodeOffset;
  int AGateMidEquSourcePrimeNodeOffset;
  int AGateMidEquGatePrimeNodeOffset;
  int AGateMidEquGateMidNodeOffset;

  // These offset are for the voltage sources that represent initial
  // conditions on Vds, Vgs, Vbs, Ves and Vps

  // icVDS
  int icVDSEquVsOffset;
  int icVDSEquVdOffset;
  int icVDSEquIdsOffset;

  // icVGS
  int icVGSEquVsOffset;
  int icVGSEquVgOffset;
  int icVGSEquIgsOffset;

  // icVBS
  int icVBSEquVsOffset;
  int icVBSEquVbOffset;
  int icVBSEquIbsOffset;

  // icVES
  int icVESEquVsOffset;
  int icVESEquVeOffset;
  int icVESEquIesOffset;

  // icVPS
  int icVPSEquVsOffset;
  int icVPSEquVpOffset;
  int icVPSEquIpsOffset;

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // Jacobian Matrix f-Ptrs:

  //  drain row:
  double * f_DrainEquDrainNodePtr;
  double * f_DrainEquDrainPrimeNodePtr;
  double * f_DrainEquIdsPtr;

  //  gate row:
  double * f_GateEquGateNodePtr;
  double * f_GateEquBodyNodePtr;
  double * f_GateEquDrainPrimeNodePtr;
  double * f_GateEquSourcePrimeNodePtr;
  double * f_GateEquGatePrimeNodePtr;
  double * f_GateEquGateMidNodePtr;
  double * f_GateEquIgsPtr;

  //  source row:
  double * f_SourceEquSourceNodePtr;
  double * f_SourceEquSourcePrimeNodePtr;
  double * f_SourceEquIdsPtr;
  double * f_SourceEquIgsPtr;
  double * f_SourceEquIbsPtr;
  double * f_SourceEquIesPtr;
  double * f_SourceEquIpsPtr;

  //  substrate row:
  double * f_SubstrateEquSubstrateNodePtr;
  double * f_SubstrateEquBodyNodePtr;
  double * f_SubstrateEquTemperatureNodePtr;
  double * f_SubstrateEquDrainPrimeNodePtr;
  double * f_SubstrateEquSourcePrimeNodePtr;
  double * f_SubstrateEquGatePrimeNodePtr;
  double * f_SubstrateEquGateMidNodePtr;
  double * f_SubstrateEquIesPtr;

  // external body row:
  double * f_ExtBodyEquExtBodyNodePtr;
  double * f_ExtBodyEquBodyNodePtr;
  double * f_ExtBodyEquIpsPtr;

  // body row:
  double * f_BodyEquSubstrateNodePtr;
  double * f_BodyEquExtBodyNodePtr;
  double * f_BodyEquBodyNodePtr;
  double * f_BodyEquTemperatureNodePtr;
  double * f_BodyEquDrainPrimeNodePtr;
  double * f_BodyEquSourcePrimeNodePtr;
  double * f_BodyEquGatePrimeNodePtr;
  double * f_BodyEquIbsPtr;

  // temperature row:
  double * f_TemperatureEquSubstrateNodePtr;
  double * f_TemperatureEquBodyNodePtr;
  double * f_TemperatureEquTemperatureNodePtr;
  double * f_TemperatureEquDrainPrimeNodePtr;
  double * f_TemperatureEquSourcePrimeNodePtr;
  double * f_TemperatureEquGatePrimeNodePtr;

  // drain' row:
  double * f_DrainPrimeEquDrainNodePtr;
  double * f_DrainPrimeEquSubstrateNodePtr;
  double * f_DrainPrimeEquBodyNodePtr;
  double * f_DrainPrimeEquTemperatureNodePtr;
  double * f_DrainPrimeEquDrainPrimeNodePtr;
  double * f_DrainPrimeEquSourcePrimeNodePtr;
  double * f_DrainPrimeEquGatePrimeNodePtr;
  double * f_DrainPrimeEquGateMidNodePtr;

  // source' row:
  double * f_SourcePrimeEquSourceNodePtr;
  double * f_SourcePrimeEquSubstrateNodePtr;
  double * f_SourcePrimeEquBodyNodePtr;
  double * f_SourcePrimeEquTemperatureNodePtr;
  double * f_SourcePrimeEquDrainPrimeNodePtr;
  double * f_SourcePrimeEquSourcePrimeNodePtr;
  double * f_SourcePrimeEquGatePrimeNodePtr;
  double * f_SourcePrimeEquGateMidNodePtr;

  // gate' row:
  double * f_GatePrimeEquGateNodePtr;
  double * f_GatePrimeEquSubstrateNodePtr;
  double * f_GatePrimeEquBodyNodePtr;
  double * f_GatePrimeEquTemperatureNodePtr;
  double * f_GatePrimeEquDrainPrimeNodePtr;
  double * f_GatePrimeEquSourcePrimeNodePtr;
  double * f_GatePrimeEquGatePrimeNodePtr;
  double * f_GatePrimeEquGateMidNodePtr;

  // gate mid row:
  double * f_GateMidEquGateNodePtr;
  double * f_GateMidEquSubstrateNodePtr;
  double * f_GateMidEquBodyNodePtr;
  double * f_GateMidEquDrainPrimeNodePtr;
  double * f_GateMidEquSourcePrimeNodePtr;
  double * f_GateMidEquGatePrimeNodePtr;
  double * f_GateMidEquGateMidNodePtr;

  // These offset are for the voltage sources that represent initial
  // conditions on Vds, Vgs, Vbs, Ves and Vps

  // f_icVDS
  double * f_icVDSEquVsPtr;
  double * f_icVDSEquVdPtr;
  double * f_icVDSEquIdsPtr;

  // f_icVGS
  double * f_icVGSEquVsPtr;
  double * f_icVGSEquVgPtr;
  double * f_icVGSEquIgsPtr;

  // f_icVBS
  double * f_icVBSEquVsPtr;
  double * f_icVBSEquVbPtr;
  double * f_icVBSEquIbsPtr;

  // f_icVES
  double * f_icVESEquVsPtr;
  double * f_icVESEquVePtr;
  double * f_icVESEquIesPtr;

  // f_icVPS
  double * f_icVPSEquVsPtr;
  double * f_icVPSEquVpPtr;
  double * f_icVPSEquIpsPtr;

  // Jacobian Matrix q-Ptrs:

  //  drain row:
  double * q_DrainEquDrainNodePtr;
  double * q_DrainEquDrainPrimeNodePtr;
  double * q_DrainEquIdsPtr;

  //  gate row:
  double * q_GateEquGateNodePtr;
  double * q_GateEquBodyNodePtr;
  double * q_GateEquDrainPrimeNodePtr;
  double * q_GateEquSourcePrimeNodePtr;
  double * q_GateEquGatePrimeNodePtr;
  double * q_GateEquGateMidNodePtr;
  double * q_GateEquIgsPtr;

  //  source row:
  double * q_SourceEquSourceNodePtr;
  double * q_SourceEquSourcePrimeNodePtr;
  double * q_SourceEquIdsPtr;
  double * q_SourceEquIgsPtr;
  double * q_SourceEquIbsPtr;
  double * q_SourceEquIesPtr;
  double * q_SourceEquIpsPtr;

  //  substrate row:
  double * q_SubstrateEquSubstrateNodePtr;
  double * q_SubstrateEquBodyNodePtr;
  double * q_SubstrateEquTemperatureNodePtr;
  double * q_SubstrateEquDrainPrimeNodePtr;
  double * q_SubstrateEquSourcePrimeNodePtr;
  double * q_SubstrateEquGatePrimeNodePtr;
  double * q_SubstrateEquGateMidNodePtr;
  double * q_SubstrateEquIesPtr;

  // external body row:
  double * q_ExtBodyEquExtBodyNodePtr;
  double * q_ExtBodyEquBodyNodePtr;
  double * q_ExtBodyEquIpsPtr;

  // body row:
  double * q_BodyEquSubstrateNodePtr;
  double * q_BodyEquExtBodyNodePtr;
  double * q_BodyEquBodyNodePtr;
  double * q_BodyEquTemperatureNodePtr;
  double * q_BodyEquDrainPrimeNodePtr;
  double * q_BodyEquSourcePrimeNodePtr;
  double * q_BodyEquGatePrimeNodePtr;
  double * q_BodyEquIbsPtr;

  // temperature row:
  double * q_TemperatureEquSubstrateNodePtr;
  double * q_TemperatureEquBodyNodePtr;
  double * q_TemperatureEquTemperatureNodePtr;
  double * q_TemperatureEquDrainPrimeNodePtr;
  double * q_TemperatureEquSourcePrimeNodePtr;
  double * q_TemperatureEquGatePrimeNodePtr;

  // drain' row:
  double * q_DrainPrimeEquDrainNodePtr;
  double * q_DrainPrimeEquSubstrateNodePtr;
  double * q_DrainPrimeEquBodyNodePtr;
  double * q_DrainPrimeEquTemperatureNodePtr;
  double * q_DrainPrimeEquDrainPrimeNodePtr;
  double * q_DrainPrimeEquSourcePrimeNodePtr;
  double * q_DrainPrimeEquGatePrimeNodePtr;
  double * q_DrainPrimeEquGateMidNodePtr;

  // source' row:
  double * q_SourcePrimeEquSourceNodePtr;
  double * q_SourcePrimeEquSubstrateNodePtr;
  double * q_SourcePrimeEquBodyNodePtr;
  double * q_SourcePrimeEquTemperatureNodePtr;
  double * q_SourcePrimeEquDrainPrimeNodePtr;
  double * q_SourcePrimeEquSourcePrimeNodePtr;
  double * q_SourcePrimeEquGatePrimeNodePtr;
  double * q_SourcePrimeEquGateMidNodePtr;

  // gate' row:
  double * q_GatePrimeEquGateNodePtr;
  double * q_GatePrimeEquSubstrateNodePtr;
  double * q_GatePrimeEquBodyNodePtr;
  double * q_GatePrimeEquTemperatureNodePtr;
  double * q_GatePrimeEquDrainPrimeNodePtr;
  double * q_GatePrimeEquSourcePrimeNodePtr;
  double * q_GatePrimeEquGatePrimeNodePtr;
  double * q_GatePrimeEquGateMidNodePtr;

  // gate mid row:
  double * q_GateMidEquGateNodePtr;
  double * q_GateMidEquSubstrateNodePtr;
  double * q_GateMidEquBodyNodePtr;
  double * q_GateMidEquDrainPrimeNodePtr;
  double * q_GateMidEquSourcePrimeNodePtr;
  double * q_GateMidEquGatePrimeNodePtr;
  double * q_GateMidEquGateMidNodePtr;

  // These offset are for the voltage sources that represent initial
  // conditions on Vds, Vgs, Vbs, Ves and Vps

  // q_icVDS
  double * q_icVDSEquVsPtr;
  double * q_icVDSEquVdPtr;
  double * q_icVDSEquIdsPtr;

  // q_icVGS
  double * q_icVGSEquVsPtr;
  double * q_icVGSEquVgPtr;
  double * q_icVGSEquIgsPtr;

  // q_icVBS
  double * q_icVBSEquVsPtr;
  double * q_icVBSEquVbPtr;
  double * q_icVBSEquIbsPtr;

  // q_icVES
  double * q_icVESEquVsPtr;
  double * q_icVESEquVePtr;
  double * q_icVESEquIesPtr;

  // q_icVPS
  double * q_icVPSEquVsPtr;
  double * q_icVPSEquVpPtr;
  double * q_icVPSEquIpsPtr;

  // end of matrix pointer section
#endif

  // Most versions of this device will use one of the jacobian's defined
  // as static above.  However, if any of the five initial conditions are
  // specified (icVDS, icVGS, icVBS, icVES or icVPS), then the number of
  // jacobian permutations goes way up. Rather than enumerate all of
  // these as static vectors, we'll make non-static, member variables of
  // the jacStamp for those few devices that specify initial conditoins.

  std::vector< std::vector< int > > jacStampIC;
  std::vector<int>  jacMapIC;
  std::vector< std::vector<int> > jacMapIC2;

  // boolean for new voltlim debugging:
  bool vlDebug;

  int blockHomotopyID; // For homotopy
  double randomPerturb; // For homotopy
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Dave Shirley
// Creation Date : 05/20/04
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
private:
  std::vector<Instance*> instanceContainer;

private:

  int    dtype;

  double cbox;
  double csi;

  int    mobMod;
  int    capMod;
  int    binUnit;
  int    paramChk;
  std::string version;
  // Bogus PSpice-y thing where you can set L and W on the model *or* the
  // instance line, with instance overriding model
  double model_l;
  double model_w;
  // range of channel dimensions for which model is valid
  double Lmax;
  double Lmin;
  double Wmax;
  double Wmin;
  //
  //soi specific variables:
  double dtoxcv;
  double npeak;
  double pdibl1;
  double pdibl2;
  double pdiblb;
  int shMod;
  double tbox;
  double tsi;
  double rth0;
  double cth0;
  double ngidl;
  double agidl;
  double bgidl;
  double ndiode;
  double xbjt;
  double xdif;
  double xrec;
  double xtun;
  double bodyJctGateSideGradingCoeff;
  int fnoiMod;
  int tnoiMod;
  double tnoia;
  double tnoib;
  double rnoia;
  double rnoib;
  double ntnoi;
  double noif;
  double k1w1;
  double k1w2;
  double ketas;
  double dwbc;
  double beta1;
  double beta2;
  double vdsatii0;
  double tii;
  double lii;
  double sii0;
  double sii1;
  double sii2;
  double siid;
  double fbjtii;
  double esatii;
  double ntun;
  double nrecf0;
  double nrecr0;
  double isbjt;
  double isdif;
  double isrec;
  double istun;
  double ln;
  double vrec0;
  double vtun0;
  double nbjt;
  double lbjt0;
  double ldif0;
  double vabjt;
  double aely;
  double ahli;
  double rbody;
  double rbsh;
  double cgeo;
  double tt;
  double ndif;
  double vsdfb;
  double vsdth;
  double csdmin;
  double asd;
  double csdesw;
  double ntrecf;
  double ntrecr;
  double dlcb;
  double fbody;
  double delvt;
  double kb1;
  double dlbg;
  int igbMod;
  int igcMod;
  double toxqm;
  double wth0;
  double rhalo;
  double ntox;
  double toxref;
  double ebg;
  double vevb;
  double alphaGB1;
  double betaGB1;
  double vgb1;
  double vecb;
  double alphaGB2;
  double betaGB2;
  double vgb2;
  double voxh;
  double deltavox;
  double aigc;
  double bigc;
  double cigc;
  double aigsd;
  double bigsd;
  double cigsd;
  double nigc;
  double pigcd;
  double poxedge;
  double dlcig;
  int soiMod;
  double vbs0pd;
  double vbs0fd;
  double vbsa;
  double nofffd;
  double vofffd;
  double k1b;
  double k2b;
  double dk2b;
  double dvbd0;
  double dvbd1;
  double moinFD;
  int rgateMod;
  int bug1830fix;

  double xrcrg1;
  double xrcrg2;
  double rshg;
  double ngcon;
  double xgw;
  double xgl;
  double lalphaGB1;
  double lbetaGB1;
  double lalphaGB2;
  double lbetaGB2;
  double lndif;
  double lntrecf;
  double lntrecr;
  double lxbjt;
  double lxdif;
  double lxrec;
  double lxtun;
  double laigc;
  double lbigc;
  double lcigc;
  double laigsd;
  double lbigsd;
  double lcigsd;
  double lnigc;
  double lpigcd;
  double lpoxedge;
  double lnpeak;
  double lk1w1;
  double lk1w2;
  double lketas;
  double lpdibl1;
  double lpdibl2;
  double lpdiblb;
  double lfbjtii;
  double lbeta1;
  double lbeta2;
  double lvdsatii0;
  double llii;
  double lesatii;
  double lsii0;
  double lsii1;
  double lsii2;
  double lsiid;
  double lkb1;
  double lagidl;
  double lbgidl;
  double lngidl;
  double lntun;
  double lndiode;
  double lnrecf0;
  double lnrecr0;
  double lisbjt;
  double lisdif;
  double lisrec;
  double listun;
  double lvrec0;
  double lvtun0;
  double lnbjt;
  double llbjt0;
  double lvabjt;
  double laely;
  double lahli;
  double lvsdfb;
  double lvsdth;
  double ldelvt;
  double lxrcrg1;
  double lxrcrg2;
  double walphaGB1;
  double wbetaGB1;
  double walphaGB2;
  double wbetaGB2;
  double wndif;
  double wntrecf;
  double wntrecr;
  double wxbjt;
  double wxdif;
  double wxrec;
  double wxtun;
  double waigc;
  double wbigc;
  double wcigc;
  double waigsd;
  double wbigsd;
  double wcigsd;
  double wnigc;
  double wpigcd;
  double wpoxedge;
  double wnpeak;
  double wk1w1;
  double wk1w2;
  double wkb1;
  double wketas;
  double wpdibl1;
  double wpdibl2;
  double wpdiblb;
  double wfbjtii;
  double wbeta1;
  double wbeta2;
  double wvdsatii0;
  double wlii;
  double wesatii;
  double wsii0;
  double wsii1;
  double wsii2;
  double wsiid;
  double wagidl;
  double wbgidl;
  double wngidl;
  double wntun;
  double wndiode;
  double wnrecf0;
  double wnrecr0;
  double wisbjt;
  double wisdif;
  double wisrec;
  double wistun;
  double wvrec0;
  double wvtun0;
  double wnbjt;
  double wlbjt0;
  double wvabjt;
  double waely;
  double wahli;
  double wvsdfb;
  double wvsdth;
  double wdelvt;
  double wxrcrg1;
  double wxrcrg2;
  double palphaGB1;
  double pbetaGB1;
  double palphaGB2;
  double pbetaGB2;
  double pndif;
  double pntrecf;
  double pntrecr;
  double pxbjt;
  double pxdif;
  double pxrec;
  double pxtun;
  double paigc;
  double pbigc;
  double pcigc;
  double paigsd;
  double pbigsd;
  double pcigsd;
  double pnigc;
  double ppigcd;
  double ppoxedge;
  double pnpeak;
  double pk1w1;
  double pk1w2;
  double pkb1;
  double pketas;
  double ppdibl1;
  double ppdibl2;
  double ppdiblb;
  double pfbjtii;
  double pbeta1;
  double pbeta2;
  double pvdsatii0;
  double plii;
  double pesatii;
  double psii0;
  double psii1;
  double psii2;
  double psiid;
  double pagidl;
  double pbgidl;
  double pngidl;
  double pntun;
  double pndiode;
  double pnrecf0;
  double pnrecr0;
  double pisbjt;
  double pisdif;
  double pisrec;
  double pistun;
  double pvrec0;
  double pvtun0;
  double pnbjt;
  double plbjt0;
  double pvabjt;
  double paely;
  double pahli;
  double pvsdfb;
  double pvsdth;
  double pdelvt;
  double pxrcrg1;
  double pxrcrg2;

  bool npeakGiven;
  bool vsdthGiven;
  bool vsdfbGiven;
  bool csdminGiven;
  bool  gamma1Given;
  bool  gamma2Given;
  bool  vbxGiven;
  bool  vbmGiven;
  bool  xtGiven;
  bool  k1Given;
  bool  k2Given;

  double vcrit;
  double vtm;

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
  double pvag;
  double wr;
  double dwg;
  double dwb;
  double b0;
  double b1;
  double alpha0;
  double beta0;

  // CV model
  double cgsl;
  double cgdl;
  double ckappa;
  double cf;
  double clc;
  double cle;
  double dwc;
  double dlc;
  double noff;
  double acde;
  double moin;
  double tcjswg;
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
  double lngate;
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
  double lpvag;
  double lwr;
  double ldwg;
  double ldwb;
  double lb0;
  double lb1;
  double lalpha0;
  double lbeta0;

  // CV model
  double lcgsl;
  double lcgdl;
  double lckappa;
  double lnoff;
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
  double wngate;
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
  double wpvag;
  double wwr;
  double wdwg;
  double wdwb;
  double wb0;
  double wb1;
  double walpha0;
  double wbeta0;

  // CV model
  double wcgsl;
  double wcgdl;
  double wckappa;
  double wnoff;
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
  double pngate;
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
  double ppvag;
  double pwr;
  double pdwg;
  double pdwb;
  double pb0;
  double pb1;
  double palpha0;
  double pbeta0;

  // CV model
  double pcgsl;
  double pcgdl;
  double pckappa;
  double pnoff;
  double pacde;
  double pmoin;

  double tnom;
  double cgso;
  double cgdo;
  double xpart;

  double sheetResistance;
  double GatesidewallJctPotential;
  double unitLengthGateSidewallJctCap;

  double Lint;
  double Ll;
  double Llc;
  double Lln;
  double Lw;
  double Lwc;
  double Lwn;
  double Lwl;
  double Lwlc;

  double Wint;
  double Wl;
  double Wlc;
  double Wln;
  double Ww;
  double Wwc;
  double Wwn;
  double Wwl;
  double Wwlc;

  // Pre-calculated constants
  // MCJ: move to size-dependent param.
  double cox;
  double factor1;

  double oxideTrapDensityA;
  double oxideTrapDensityB;
  double oxideTrapDensityC;
  double em;
  double ef;
  double af;
  double kf;

  std::list<SizeDependParam*> sizeDependParamList;

  bool  vth0Given;
  bool igbModGiven;
  // end of original 3f5 stuff.

  // Variables from the 3f5 B3SOI function b3soitemp, but are  more
  // appropriate as model variables.
  double Vtm0;
  double Eg0;
  double ni;
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

  // load functions:
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * bVec, double * leadF, double * leadQ, double * junctionV);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);

};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace MOSFET_B3SOI
} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_MOSFET_B3SOI_h


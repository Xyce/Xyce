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
// Purpose        : Transmission line.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 06/16/10
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_LTRA_h
#define Xyce_N_DEV_LTRA_h

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

// device parameters
#define LTRA_MOD_LTRA 0
#define LTRA_MOD_R 1
#define LTRA_MOD_L 2
#define LTRA_MOD_G 3
#define LTRA_MOD_C 4
#define LTRA_MOD_LEN 5
#define LTRA_V1 6
#define LTRA_I1 7
#define LTRA_V2 8
#define LTRA_I2 9
#define LTRA_IC 10
#define LTRA_MOD_RELTOL 11
#define LTRA_MOD_ABSTOL 12
#define LTRA_POS_NODE1 13
#define LTRA_NEG_NODE1 14
#define LTRA_POS_NODE2 15
#define LTRA_NEG_NODE2 16
#define LTRA_INPUT1 17
#define LTRA_INPUT2 18
#define LTRA_DELAY 19
#define LTRA_BR_EQ1 20
#define LTRA_BR_EQ2 21
#define LTRA_MOD_NL 22
#define LTRA_MOD_FREQ 23

#define LTRA_MOD_FULLCONTROL 26
#define LTRA_MOD_HALFCONTROL 27
#define LTRA_MOD_NOCONTROL 28
#define LTRA_MOD_PRINT 29
#define LTRA_MOD_NOPRINT 30

#define LTRA_MOD_STEPLIMIT  32
#define LTRA_MOD_NOSTEPLIMIT 33
#define LTRA_MOD_LININTERP 34
#define LTRA_MOD_QUADINTERP 35
#define LTRA_MOD_MIXEDINTERP 36
#define LTRA_MOD_RLC  37
#define LTRA_MOD_RC 38
#define LTRA_MOD_RG 39
#define LTRA_MOD_LC 40
#define LTRA_MOD_RL 41
#define LTRA_MOD_STLINEREL 42
#define LTRA_MOD_STLINEABS 43
#define LTRA_MOD_CHOPREL 44
#define LTRA_MOD_CHOPABS 45
#define LTRA_MOD_TRUNCNR 46
#define LTRA_MOD_TRUNCDONTCUT 47

namespace Xyce {
namespace Device {
namespace LTRA {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance>
{
  static const char *name() {return "Lossy Transmission Line";}
  static const char *deviceTypeName() {return "O level 1";}
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend struct Traits;
  friend class Master;

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

  void loadNodeSymbols(Util::SymbolTable &symbol_table) const; // override

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();

  bool loadDAEQVector () {return true;}
  bool loadDAEFVector ();

  bool loadDAEdQdx () {return true;}
  bool loadDAEdFdx ();

  void setupPointers();

  bool getInstanceBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes);
  void acceptStep();

  double getMaxTimeStepSize();
  virtual bool maxTimeStepSupported () {return true;};

  DeviceState * getInternalState();
  bool setInternalState( const DeviceState & state );
  bool setIC();

public:
  Model &getModel()
  {
    return model_;
  }

private:
  static std::vector< std::vector<int> > jacStamp;

  Model &       model_;         //< Owning model

  void calculateMaxTimeStep_(void);   // Calculate a maximum time step size to minimize errors

  double input1;	// accumulated excitation for port 1
  double input2;	// accumulated excitation for port 2

  double currp1;	// Current at port 1
  double currp2;	// Current at port 2

  double vpos1;	// Voltage at port 1 positive node
  double vneg1;	// Voltage at port 1 negative node

  double vpos2;	// Voltage at port 2 positive node
  double vneg2;	// Voltage at port 2 negative node

  double initVolt1;	// initial condition:  voltage on port 1
  double initVolt2;	// initial condition:  voltage on port 2

  double initCur1;	// initial condition:  current at port 1
  double initCur2;	// initial condition:  current at port 2

  std::vector<double> v1;	// past values of voltage at port 1
  std::vector<double> v2;	// past values of voltage at port 2

  std::vector<double> i1;	// past values of current at port 1
  std::vector<double> i2;	// past values of current at port 2

  int listSize;	// Size of variables vectors above

  bool initVolt1Given;
  bool initVolt2Given;

  bool initCur1Given;
  bool initCur2Given;

  // local indices (offsets)
  int li_Pos1;
  int li_Neg1;

  int li_Pos2;
  int li_Neg2;

  int li_Ibr1;
  int li_Ibr2;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Matrix elements
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //    Matrix equation offset variables
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Positive Node1
  int APos1EquPos1NodeOffset;
  int APos1EquIbr1NodeOffset;
  // Negative Node1
  int ANeg1EquNeg1NodeOffset;
  int ANeg1EquIbr1NodeOffset;
  // Positive Node2
  int APos2EquPos2NodeOffset;
  int APos2EquIbr2NodeOffset;
  // Negative Node1
  int ANeg2EquNeg2NodeOffset;
  int ANeg2EquIbr2NodeOffset;
  // Branch1 equation due to Jajeet LTRA
  int AIbr1EquPos1NodeOffset;
  int AIbr1EquNeg1NodeOffset;
  int AIbr1EquPos2NodeOffset;
  int AIbr1EquNeg2NodeOffset;
  int AIbr1EquIbr1NodeOffset;
  int AIbr1EquIbr2NodeOffset;
  // Branch2 equation due to Jajeet LTRA
  int AIbr2EquPos1NodeOffset;
  int AIbr2EquNeg1NodeOffset;
  int AIbr2EquPos2NodeOffset;
  int AIbr2EquNeg2NodeOffset;
  int AIbr2EquIbr1NodeOffset;
  int AIbr2EquIbr2NodeOffset;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //    Entries in matrix corresponding to offsets above
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  double *pos1Pos1Ptr;
  double *pos1Ibr1Ptr;

  double *neg1Neg1Ptr;
  double *neg1Ibr1Ptr;

  double *pos2Pos2Ptr;
  double *pos2Ibr2Ptr;

  double *neg2Neg2Ptr;
  double *neg2Ibr2Ptr;

  double *ibr1Pos1Ptr;
  double *ibr1Neg1Ptr;
  double *ibr1Pos2Ptr;
  double *ibr1Neg2Ptr;
  double *ibr1Ibr1Ptr;
  double *ibr1Ibr2Ptr;

  double *ibr2Pos1Ptr;
  double *ibr2Neg1Ptr;
  double *ibr2Pos2Ptr;
  double *ibr2Neg2Ptr;
  double *ibr2Ibr1Ptr;
  double *ibr2Ibr2Ptr;

  bool first_BP_call_done;

  bool newBreakPoint;
  double newBreakPointTime;
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

private:

  int quadInterp_ (double t, double t1, double t2, double t3, double& c1, double& c2, double& c3);
  int linInterp_  (double t, double t1, double t2, double& c1, double& c2);
  double intlinfunc_ (double lolimit, double hilimit,
                      double lovalue, double hivalue,
                      double t1, double t2);
  double twiceintlinfunc_(double lolimit, double hilimit,
                          double otherlolimit, double lovalue,
                          double hivalue, double t1, double t2);
  double thriceintlinfunc_(double lolimit, double hilimit,
                           double secondlolimit, double thirdlolimit,
                           double lovalue, double  hivalue,
                           double t1, double t2);

  bool modelCalculations_(int& isaved, double& qf1, double& qf2, double& qf3,
                          double& lf2, double& lf3);

  double rlcH1dashFunc_(double time, double T, double alpha, double beta);
  double rlcH2Func_(double time, double T, double alpha, double beta);
  double rlcH3dashFunc_(double time, double T, double alpha, double beta);
  double rlcH1dashTwiceIntFunc_(double time, double beta);
  double rlcH3dashIntFunc_(double time, double T, double beta);
  double rcH1dashTwiceIntFunc_(double time, double cbyr);
  double rcH2TwiceIntFunc_(double time, double rclsqr);

  double rcH3dashTwiceIntFunc_(double time, double cbyr, double rclsqr);

  // coefficient setups:
  void rcCoeffsSetup_(
     double & h1dashfirstcoeff,
     double & h2firstcoeff,
     double & h3dashfirstcoeff,
     std::vector<double> & h1dashcoeffs,
     std::vector<double> & h2coeffs,
     std::vector<double> & h3dashcoeffs,
     int listsize, double cbyr, double rclsqr, double curtime,
     const std::vector<double> & timelist, int timeindex, double reltol);

  void rlcCoeffsSetup_(
     double & h1dashfirstcoeff, double & h2firstcoeff, double & h3dashfirstcoeff,
     std::vector<double> & h1dashcoeffs, std::vector<double> & h2coeffs, std::vector<double> & h3dashcoeffs,
     int listsize, double T, double alpha, double beta, double curtime,
     const std::vector<double> & timelist, int timeindex, double reltol, int *auxindexptr);

  bool straightLineCheck_(double x1, double y1,
                          double x2, double y2,
                          double x3, double y3,
                          double reltol, double abstol);



  double lteCalculate_ ( Instance & instance,
                         double curtime );

  double SECONDDERIV_(int i, double a, double b, double c);


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


  double h1dashFirstVal; // first needed value of h1dasg at current timepoint
  double h2FirstVal;	   // first needed value of h2 at current timepoint
  double h3dashFirstVal; // first needed value of h3dash at current timepoint

  double h1dashFirstCoeff; // first needed coeff of h1dash for the current timepoint
  double h2FirstCoeff;     // first needed coeff of h2 for the current timepoint
  double h3dashFirstCoeff; // first needed coeff of h3dash for the current timepoint

  std::vector<double> h1dashCoeffs; // list of other coefficients for h1dash
  std::vector<double> h2Coeffs;     // list of other coefficients for h2
  std::vector<double> h3dashCoeffs; // list of other coefficients for h3dash

  int listSize;       // size of above lists

  // input parameters:
  double resist;
  double induct;
  double conduct;
  double capac;
  double length;
  double reltol;
  double abstol;

  bool noStepLimit;	// Don't limit step size to less than line delay
  bool stepLimit;	// Do limit step size to less than line delay
  int stepLimitType;	//
                        // (These could be combined into a single
                        //  bool, but preserving spice3 compatibility)

  bool linInterp;
  bool quadInterp;
  bool mixedInterp;

  double stLineReltol;
  double stLineAbstol;

  bool lteTimeStepControl;      // indicates whether full time step
  // control using local-truncation
  // error estimation is performed

  bool truncNR;
  bool truncDontCut;

  bool resistGiven;
  bool inductGiven;
  bool conductGiven;
  bool capacGiven;
  bool lengthGiven;
  bool reltolGiven;
  bool abstolGiven;
  bool noStepLimitGiven;
  bool stepLimitGiven;
  bool linInterpGiven;
  bool quadInterpGiven;
  bool mixedInterpGiven;
  bool stLineReltolGiven;
  bool stLineAbstolGiven;
  bool lteTimeStepControlGiven;
  bool truncNRGiven;
  bool truncDontCutGiven;

  // calculated parameters
  double td;           // propagation delay T - calculated
  double imped;        // impedance Z - calculated
  double admit;        // admittance Y - calculated
  double alpha;        // alpha - calculated
  double beta;         // beta - calculated
  double attenuation;  // e^(-beta T) - calculated
  double cByR;         // C/R - for the RC line - calculated
  double rclsqr;       // RCl^2 - for the RC line - calculated
  double intH1dash;    // \int_0^\inf h'_1(\tau) d \tau - calculated
  double intH2;        // \int_0^\inf h_2(\tau) d \tau - calculated
  double intH3dash;    // \int_0^\inf h'_3(\tau) d \tau - calculated

  double coshlrootGR;  // cosh(l*sqrt(G*R)), used for DC analysis
  double rRsLrGRorG;   // sqrt(R)*sinh(l*sqrt(G*R))/sqrt(G)
  double rGsLrGRorR;   // sqrt(G)*sinh(l*sqrt(G*R))/sqrt(R)

  int auxIndex;        // auxiliary index for h2 and h3dash
  double chopReltol;   // separate reltol for truncation of impulse responses
  double chopAbstol;   // separate abstol for truncation of impulse responses

  double maxSafeStep;  // maximum safe step for impulse response calculations
  double maxTimeStep;  // maximum time step to be used
  int howToInterp;     // indicates how to interpolate for delayed timepoint
  bool printFlag;      // flag to indicate whether debugging output should be printed

  int specialCase;     // what kind of model (RC, RLC, RL, ...)

  bool tdover;

  bool restartStoredFlag;   // flag to indicate if this particular model has
  // been saved to the restart file already (restart
  // functions are instance functions, so the potential
  // exists for storing same model multiple times.)
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
  virtual bool updateState (double * solVec, double * staVec, double * stoVec,
                            int loadType);
  virtual bool updateSecondaryState (double * staDeriv, double * stoVec) {return true;}

  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV);
  virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec,
                               double * bVec, double * leadF, double * leadQ,
                               double * junctionV, int loadType);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx);
  virtual bool loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx,
                                int loadType);
  virtual bool loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                  std::vector<Util::FreqVecEntry>& fVec,
                                  std::vector<Util::FreqVecEntry>& bVec);
  virtual bool loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                   std::vector<Util::FreqMatEntry>& dFdx);

};

void registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet);

} // namespace LTRA
} // namespace Device
} // namespace Xyce

#endif

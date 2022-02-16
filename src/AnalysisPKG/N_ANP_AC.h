//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : AC analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Ting Mei
//
// Creation Date  : 01/11
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_AC_h
#define Xyce_N_ANP_AC_h

#include <iosfwd>
#include <Teuchos_SerialDenseMatrix.hpp>  
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_RCP.hpp>
#include <N_IO_OutputMOR.h>

#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_LAS_fwd.h>
#include <N_TOP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_ANP_StepEvent.h>
#include <N_TIA_TIAParams.h>
#include <N_UTL_FixedQueue.h>
#include <N_UTL_Listener.h>
#include <N_UTL_Op.h>
#include <N_UTL_OptionBlock.h>

#include <expressionGroup.h>
#include <mainXyceExpressionGroup.h>

namespace Xyce {
namespace Nonlinear {

template <typename ScalarT>
class objectiveFunctionData;
}
}

namespace Xyce {
namespace Analysis {

std::ostream& sensStdOutput (
       const std::string idString,
       const std::vector<double> & paramVals,
       const std::vector<double> & sensitivities,
       const std::vector<double> & scaled_sensitivities,
       const std::vector<std::string> & paramNameVec,
       const std::vector<double> & param_dP,
       const std::vector<int> & numericalDiff,
       const std::vector<Xyce::Nonlinear::objectiveFunctionData<std::complex<double> > *> & objFuncDataVec,
       const std::vector<double> & objectiveVec,
       OutputMgrAdapter & outputManagerAdapter,
       std::ostream& os
       );


//-------------------------------------------------------------------------------
// This is a test class for objective function expressions.  
//
// A real class will need to be developed.
//
// This class is necessary for now b/c the "main" group doesn't know about the 
// AC solution vector and the function for getting voltage node values won't work.
//
// Possibly the "outputs" group would work as well.  Try that later.
//-------------------------------------------------------------------------------
class ACExpressionGroup : public Xyce::Util::mainXyceExpressionGroup
{
  public:
    ACExpressionGroup (
      const Teuchos::RCP<Xyce::Util::mainXyceExpressionGroup> & mainGroup,
      const Linear::BlockVector & X) 
      : 
      mainXyceExpressionGroup( 
        mainGroup->comm_, mainGroup->top_,
        mainGroup->analysisManager_,
        mainGroup->deviceManager_,
        mainGroup->outputManager_),
      X_(X)
    {};

    ~ACExpressionGroup () {};

  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval);

  private:
    const Linear::BlockVector & X_;

};

//-------------------------------------------------------------------------
// Class         : AC
// Purpose       : AC analysis class
// Special Notes :
// Creator       : Ting Mei
// Creation Date : 05/24/11
//-------------------------------------------------------------------------
class AC: public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
public:
  AC(
    AnalysisManager &                   analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Device::DeviceMgr &                   device_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager);

  virtual ~AC();

  void notify(const StepEvent &event);

  void setTIAParams(const TimeIntg::TIAParams &tia_params)
  {
    tiaParams_ = tia_params;
  }

  const TimeIntg::TIAParams &getTIAParams() const
  {
    return tiaParams_;
  }

  TimeIntg::TIAParams &getTIAParams()
  {
    return tiaParams_;
  }

  double getCurrentFreq()
  {
    return currentFreq_;
  }

  void setRFParamsRequested(const std::string & type);

  bool setACLinOptions(const Util::OptionBlock & OB);
  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setTimeIntegratorOptions(const Util::OptionBlock &option_block);
  bool setACLinSolOptions(const Util::OptionBlock &option_block);
  bool setDCLinSolOptions(const Util::OptionBlock &option_block);
  bool setDataStatements(const Util::OptionBlock & paramsBlock);
  bool convertDataToSweepParams();
  bool setSensitivityOptions(const Util::OptionBlock &option_block);
  bool setSensAnalysisParams(const Util::OptionBlock &option_block);

  bool getDCOPFlag() const;

  void finalExpressionBasedSetup();

protected:
  bool doRun();
  bool doInit();
  bool doLoopProcess();
  bool doProcessSuccessfulStep();
  bool doProcessFailedStep();
  bool doFinish();
  bool doHandlePredictor();

  void printStepHeader(std::ostream &os)
  {}

  void printProgress(std::ostream &os)
  {}

private:
  int setupSweepParam_();

  bool updateDataParams_(int stepNumber);
  bool updateCurrentFreq_(int stepNumber);
  bool createLinearSystem_();


  bool updateLinearSystem_C_and_G_();
  bool updateLinearSystemFreq_();
  bool updateLinearSystemMagAndPhase_();

  bool solveLinearSystem_();

  // sensitivity functions
  bool precomputeDCsensitivities_ ();
  bool solveSensitivity_();
  bool solveDirectSensitivity_();
  bool solveAdjointSensitivity_();
  bool loadSensitivityRHS_(int ipar);
  bool applyOmega_dJdp(int ipar);
  bool unapplyOmega_dJdp(int ipar);
  bool outputSensitivity_();

  void solve_mag_phase_Sensitivities_(
      const double dxdpReal,
      const double dxdpImag,
      const double xr,
      const double xi,
      double & dxdp_mag,
      double & dxdp_phase
      );

private:
  AnalysisManager &     analysisManager_;
  Loader::Loader &      loader_;
  Linear::System &      linearSystem_;
  Nonlinear::Manager &  nonlinearManager_;
  Topo::Topology &      topology_;





  Device::DeviceMgr &                   deviceManager_;

  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::OutputMOR                         outputMOR_;
  OutputMgrAdapter &    outputManagerAdapter_;
  TimeIntg::TIAParams   tiaParams_;
              
  Linear::Vector *      bVecRealPtr;
  Linear::Vector *      bVecImagPtr;

  int                   acLoopSize_;

  std::vector<int>      acSweepFailures_;

  bool                          stepFlag_;
  std::string                   type_;
  double                        np_;
  double                        fStart_;
  double                        fStop_;
  double                        stepMult_;
  double                        fstep_;

  double                        currentFreq_;

  
         
  std::map<int, std::pair<int, double> > portMap_;

  std::vector<int>  portPID_;
  bool sparcalc_;
  int  numPorts_;
  bool hParamsRequested_;
  bool sParamsRequested_;
  bool zParamsRequested_;

  // Y and S parameters 
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Yparams_;     
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Sparams_;
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Zparams_;
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Hparams_;

  Util::Op::RFparamsData RFparams_;

  std::vector<int> bMatEntriesVec_, portNumVec_;

  std::vector<double> Z0sVec_;

  Linear::Matrix *                C_;
  Linear::Matrix *                G_;
  Linear::BlockMatrix *           ACMatrix_;
  Linear::BlockVector *           B_;
  Linear::BlockVector *           X_;

  // sensitivity vectors
  Linear::Vector *      dbdpVecRealPtr;
  Linear::Vector *      dbdpVecImagPtr;
  Linear::Vector *      dOdxVecRealPtr;
  Linear::Vector *      dOdxVecImagPtr;

  std::vector<Linear::BlockMatrix *>  dJdpVector_;

  Linear::Matrix *                dCdp_;
  Linear::Matrix *                dGdp_;

  Linear::Matrix *                origC_;
  Linear::Matrix *                origG_;

  Linear::BlockVector *           dBdp_;
  Linear::BlockVector *           dXdp_;
  Linear::BlockVector *           lambda_;
  Linear::BlockVector *           dOdXreal_;
  Linear::BlockVector *           dOdXimag_;
  Linear::BlockVector *           savedX_;

  Linear::BlockVector *           sensRhs_;

  Linear::Solver *              blockSolver_;
  Linear::Problem *             blockProblem_;
  Util::OptionBlock             acLinSolOptionBlock_;

  Linear::Solver *              directSensSolver_;
  Linear::Problem *             directSensProblem_;
  Util::OptionBlock             dcLinSolOptionBlock_;

  SweepVector                   acSweepVector_;
  std::map< std::string, std::vector<std::string> > dataNamesMap_;
  std::map< std::string, std::vector< std::vector<double> > > dataTablesMap_;

  // parameter sensitivities 
  bool sensFlag_;
  bool solveAdjointSensitivityFlag_;
  bool solveDirectSensitivityFlag_;
  bool outputScaledFlag_; // include scaled sensitivities in IO 
  bool outputUnscaledFlag_; // include unscaled sensitivities in IO
  int maxParamStringSize_;
  bool stdOutputFlag_;

  int numSensParams_;
  std::vector<double>   objectiveVec_;
  std::vector<double>   dOdpVec_;
  std::vector<double>   dOdpAdjVec_;
  std::vector<double>   scaled_dOdpVec_;
  std::vector<double>   scaled_dOdpAdjVec_;

  // objective function objects
  //std::vector<std::string> objFuncStrings_;
  // these are for the new-expression method of specifying objective functions
  bool objFuncGiven_;
  bool objFuncGIDsetup_;
  std::vector<Xyce::Nonlinear::objectiveFunctionData<std::complex<double> > *> objFuncDataVec_;
  std::vector<std::string> objFuncStrings_; // needed for output

  std::vector<std::string> paramNameVec_;
  std::vector<double> param_dP_;
  std::vector<int> numericalDiff_;

  // finite difference variables
  //int difference_;
  bool forceFD_;
  bool forceDeviceFD_;
  bool forceAnalytic_;
  bool newLowMem_;
  bool sparseAdjointStorage_;
  bool reuseFactors_;

};

bool registerACFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_AC_h

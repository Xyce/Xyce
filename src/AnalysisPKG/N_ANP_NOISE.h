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
// Purpose        : NOISE analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Eric Keiter
//
// Creation Date  : 01/11
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_NOISE_h
#define Xyce_N_ANP_NOISE_h

#include <Teuchos_SerialDenseMatrix.hpp>  

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
#include <N_UTL_OptionBlock.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Analysis {

// used for the output for the -noise_names_file command line option
bool outputNoiseNameFile(Parallel::Machine comm, const std::string & path,
                         Loader::Loader & circuitLoader);

//-------------------------------------------------------------------------
// Class         : NOISE
// Purpose       : NOISE analysis class
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 12/8/2014
//-------------------------------------------------------------------------
class NOISE: public AnalysisBase, public Util::ListenerAutoSubscribe<StepEvent>
{
  public:
  NOISE(
    AnalysisManager &                   analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager);

  ~NOISE();

  void notify(const StepEvent &event);

  const TimeIntg::TIAParams &getTIAParams() const
  {
    return tiaParams_;
  }

  TimeIntg::TIAParams &getTIAParams()
  {
    return tiaParams_;
  }

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setTimeIntegratorOptions(const Util::OptionBlock &option_block);
  bool setACLinSolOptions(const Util::OptionBlock &option_block);
  bool setDataStatements(const Util::OptionBlock & paramsBlock);
  bool convertDataToSweepParams();

  bool getDCOPFlag() const;

protected:
  bool doRun();
  bool doInit();
  bool doLoopProcess();
  bool doProcessSuccessfulStep();
  bool doProcessFailedStep();
  bool doFinish();
  bool doHandlePredictor();

public:
  void printStepHeader(std::ostream &os)
  {}

  void printProgress(std::ostream &os)
  {}

private:
  int setupSweepParam_();

  bool updateDataParams_(int stepNumber);
  bool updateCurrentFreq_(int stepNumber);

  bool createACLinearSystem_();

  bool updateACLinearSystem_C_and_G_();
  bool updateACLinearSystemFreq_();
  bool updateACLinearSystemMagAndPhase_();
  bool solveACLinearSystem_();

  void resetAdjointNOISELinearSystem_();
  bool solveAdjointNOISE_();
  void setupAdjointRHS_();

  double noiseIntegral( double noizDens, double lnNdens, double lnNlstDens, 
      double delLnFreq, double delFreq, double lnFreq, double lnLastFreq);

  void processOutputNodes ();

  void hackTecplotOutput();

  std::ostream& noiseOutputToScreen_(std::ostream& os);

private:
  AnalysisManager &                     analysisManager_;
  Loader::Loader &                      loader_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  TimeIntg::TIAParams                   tiaParams_;

  // AC B-vectors
  Linear::Vector * bVecRealPtr;
  Linear::Vector * bVecImagPtr;

  // NOISE B-vectors
  Linear::Vector * bNoiseVecRealPtr;
  Linear::Vector * bNoiseVecImagPtr;

  bool hackOutputCalledBefore_;
  bool outputNodeSingle_;
  std::string outputNode1_;
  std::string outputNode2_;
  std::string specifiedSource_;
  int noiseLoopSize_;

  std::vector<int>      acSweepFailures_;
  std::vector<int>              noiseSweepFailures_;

  bool                          stepFlag_;
  std::string                   type_;
  double                        np_;
  double                        fStart_;
  double                        fStop_;

  double                        stepMult_;
  double                        fstep_;
  int                           pts_per_summary_;
  bool                          pts_per_summary_Given;

  // noise integrals are not calculated for DATA=<name> case if the
  // specified frequencies are not monotonically increasing
  bool                          calcNoiseIntegrals_;

  double                        delFreq_;
  double                        lastFreq_;
  double                        currentFreq_;
  double                        lnFreq_;
  double                        lnLastFreq_;
  double                        delLnFreq_;

  double                        GainSqInv_;
  double                        lnGainInv_;

  double                        totalOutputNoise_;
  double                        totalInputNoise_;

  double                        totalOutputNoiseDens_;
  double                        totalInputNoiseDens_;

  Linear::Matrix *                C_;
  Linear::Matrix *                G_;
  Linear::BlockMatrix *           ACMatrix_;
  Linear::BlockVector *           B_;
  Linear::BlockVector *           X_;
  Linear::BlockVector *           saved_AC_X_;

  Linear::Solver *              blockSolver_;
  Linear::Problem *             blockProblem_;
  Util::OptionBlock             acLinSolOptionBlock_;

  SweepVector                   noiseSweepVector_;
  std::map< std::string, std::vector<std::string> > dataNamesMap_;
  std::map< std::string, std::vector< std::vector<double> > > dataTablesMap_;

  std::vector<std::string> outputVarNames_;
  std::vector<int>    outputVarGIDs_;

  std::vector<Xyce::Analysis::NoiseData*> noiseDataVec_;
  Teuchos::SerialDenseMatrix<int, std::complex<double> > Sparams_;

  std::vector<double> objectiveVec_; 
  std::vector<double> dOdpVec_; 
  std::vector<double> dOdpAdjVec_;
  std::vector<double> scaled_dOdpVec_;
  std::vector<double> scaled_dOdpAdjVec_;
};

bool registerNOISEFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_NOISE_h

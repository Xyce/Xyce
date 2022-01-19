//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electronic Simulator.
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
// Purpose        : EmbeddedSampling analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 8/14/2017
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_EmbeddedSampling_h
#define Xyce_N_ANP_EmbeddedSampling_h

#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_IO_fwd.h>

#include <N_ANP_UQ_fwd.h>
#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>

#ifdef Xyce_STOKHOS_ENABLE
// make sure linking against the correct trilinos!
#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"
#endif

#include <N_IO_OptionBlock.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : EmbeddedSampling
// Purpose       : EmbeddedSampling analysis class
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-------------------------------------------------------------------------
class EmbeddedSampling : public AnalysisBase
{
public:
  EmbeddedSampling(
      AnalysisManager &analysis_manager, 
      Linear::System & linear_system,
      Nonlinear::Manager & nonlinear_manager,
      Device::DeviceMgr & device_manager,
      Linear::Builder & builder,
      Loader::Loader &loader, 
      Topo::Topology & topology, 
      IO::InitialConditionsManager & initial_conditions_manager,
      AnalysisBase &child_analysis);

  virtual ~EmbeddedSampling();

  void notify(const StepEvent &event);

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setEmbeddedSamplingOptions(const Util::OptionBlock & option_block);
  bool setDCOptions (const Util::OptionBlock & paramsBlock);

  const TimeIntg::TIAParams &getTIAParams() const; // override
  TimeIntg::TIAParams &getTIAParams(); // override

  virtual bool getDCOPFlag() const;

  bool setLinSol(const Util::OptionBlock & OB) { saved_lsOB_ = OB; return true; }
  bool setESLinSol(const Util::OptionBlock & OB) { saved_lsESOB_ = OB; return true; }

  void stepCallBack();

  void setBeginningIntegrationFlag(bool bif)
  {
    childAnalysis_.setBeginningIntegrationFlag(bif);
  }

  bool getBeginningIntegrationFlag()
  {
    return childAnalysis_.getBeginningIntegrationFlag();
  }

  void setIntegrationMethod(int im)
  {
    childAnalysis_.setIntegrationMethod(im);
  }

  int getIntegrationMethod ()
  {
    return childAnalysis_.getIntegrationMethod();
  }

  void setStepNumber(int step)
  {
    childAnalysis_.setStepNumber(step);
  }

  unsigned int getStepNumber()
  {
    return childAnalysis_.getStepNumber();
  }

  void setTranStepNumber(int step)
  {
    childAnalysis_.setTranStepNumber(step);
  }

  int getTranStepNumber()
  {
    return childAnalysis_.getTranStepNumber();
  }

  virtual void finalExpressionBasedSetup() {};

protected:
  virtual bool doRun();

  virtual bool doInit();
  void setupBlockSystemObjects (); // called from doInit
  void setupStokhosObjects ();

  virtual bool doLoopProcess();
  virtual bool doProcessSuccessfulStep();
  virtual bool doProcessFailedStep();
  virtual bool doFinish();
  virtual bool doHandlePredictor() { return true; }


  void computeEnsembleOutputs();
  void hackEnsembleOutput();

private:
  AnalysisManager &     analysisManager_;
  Loader::Loader &      loader_;
  Linear::System &      linearSystem_;
  Nonlinear::Manager &  nonlinearManager_;
  Device::DeviceMgr &   deviceManager_;
  Linear::Builder &     builder_;
  Topo::Topology &      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  OutputMgrAdapter &    outputManagerAdapter_;
  IO::Measure::Manager &    measureManager_;
  AnalysisBase &        childAnalysis_;


  Loader::ESLoader *                    esLoaderPtr_; 
  Teuchos::RCP<Linear::ESBuilder>       esBuilderPtr_;
  Linear::System *                      esLinearSystem_;
  Linear::ESSolverFactory *             solverFactory_;

  // Linear solver options
  Util::OptionBlock                     saved_lsOB_;
  Util::OptionBlock                     saved_lsESOB_;

  SweepVector           samplingVector_;
  SweepVector           dcSweepVector_;
  int                   maxParamStringSize_;

  Parallel::Manager *   pdsMgrPtr_;

  std::vector<std::string> paramNameVec_;
  std::vector<std::string> typeVec_;

  // normal
  std::vector<double> meanVec_;
  std::vector<double> stdDevVec_;

  // uniform (and all the others, actually as bounds are applied to all)
  std::vector<double> lower_bounds_Vec_;
  std::vector<double> upper_bounds_Vec_;
  bool lower_bounds_Given_;
  bool upper_bounds_Given_;

  // gamma
  std::vector<double> alphaVec_;
  std::vector<double> betaVec_;

  // sample points arrays
  std::vector<double> X_; // uncorrelated
  std::vector<double> Y_; // correlated
  std::vector< double > covMatrix_;

  bool covMatrixGiven_;

  int numSamples_;
  bool numSamplesGiven_;
  UQ::SampleType  sampleType_;

  int userSeed_;
  bool userSeedGiven_;

  std::string hackOutputFormat_;
  bool hackOutputCalledBefore_;
  bool hackOutputAllSamples_;
  bool outputtersCalledBefore_;
  bool outputSampleStats_;

  bool paramsOuterLoop_;

#ifdef Xyce_STOKHOS_ENABLE
  bool regressionPCEenable_;
  bool projectionPCEenable_;
  int PCEorder_;

  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > regrBases; 
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > regrBasis;

  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > quadBases; 
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > quadBasis;

  // Quadrature method
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quadMethod;

  // Triple product tensor
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > regrCijk;
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > quadCijk;

  // Expansion method
  Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > regrExpn;
  Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > quadExpn;

  bool resamplePCE_;
  bool outputPCECoeffs_;
  std::vector<std::string> regressionPCEcoeffs_;
  std::vector<std::string> projectionPCEcoeffs_;

  int numResamples_;
  bool numResamplesGiven_;

  bool useSparseGrid_;
#endif

  bool stdOutputFlag_;
  int debugLevel_;

  bool outputsGiven_;
  bool outputsSetup_;
  std::vector<UQ::outputFunctionData*> outFuncDataVec_;

  bool measuresGiven_;
  std::vector<UQ::outputFunctionData*> measFuncDataVec_;

  bool outFuncGIDsetup_;

  int outputIndex_;

  // non-block objects:
  Linear::Vector * nextSolutionPtr_;
  Linear::Vector * nextStatePtr_;
  Linear::Vector * nextStorePtr_;
 
  bool resetForStepCalledBefore_;
  bool useExpressionSamples_;
};

bool registerEmbeddedSamplingFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_EmbeddedSampling_h

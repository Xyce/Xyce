//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2017 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Purpose        : PCE analysis class
//
// Special Notes  : 
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 6/3/2019
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_PCE_h
#define Xyce_N_ANP_PCE_h

#ifdef Xyce_STOKHOS_ENABLE
#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>

#include <N_ANP_UQ_fwd.h>
#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>

#ifdef Xyce_STOKHOS_ENABLE
// make sure linking against the correct trilinos!
#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"
#include <Stokhos_Sparse3TensorUtilities.hpp>
#endif

#include <N_IO_OptionBlock.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : PCE
// Purpose       : PCE analysis class
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-------------------------------------------------------------------------
class PCE : public AnalysisBase
{
public:
  PCE(
      AnalysisManager &analysis_manager, 
      Linear::System & linear_system,
      Nonlinear::Manager & nonlinear_manager,
      Device::DeviceMgr & device_manager,
      Linear::Builder & builder,
      Loader::Loader &loader, 
      Topo::Topology & topology, 
      IO::InitialConditionsManager & initial_conditions_manager,
      AnalysisBase &child_analysis);

  virtual ~PCE();

  void notify(const StepEvent &event);

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setPCEOptions(const Util::OptionBlock & option_block);
  bool setDCOptions (const Util::OptionBlock & paramsBlock);

  const TimeIntg::TIAParams &getTIAParams() const; // override
  TimeIntg::TIAParams &getTIAParams(); // override

  virtual bool getDCOPFlag() const;

  bool setLinSol(const Util::OptionBlock & OB) { saved_lsOB_ = OB; return true; }
  bool setPCELinSol(const Util::OptionBlock & OB) { saved_lsPCEOB_ = OB; return true; }

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


  void convertPointToPCE (int index, Stokhos::OrthogPolyApprox<int,double> & pceApprox);
  void evaluateVector (Teuchos::RCP<Linear::BlockVector> & bX_quad_ptr_);
  void outputXvectors(); // hack output, X vectors to stdout
  void hackPCEOutput2(); // new one, outputs *everything* to a plot file

  void computePCEOutputs();
  void hackPCEOutput(); // derived from the ES hack output

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


  Loader::PCELoader *                    pceLoaderPtr_; 
  Teuchos::RCP<Linear::PCEBuilder>       pceBuilderPtr_;
  Linear::System *                      pceLinearSystem_;
  Linear::PCESolverFactory *            solverFactory_;

  // Linear solver and nonlinear solver options
  Util::OptionBlock                     saved_lsOB_;
  Util::OptionBlock                     saved_lsPCEOB_;

  SweepVector           samplingVector_;
  SweepVector           dcSweepVector_;
  int                   maxParamStringSize_;

  Parallel::Manager *       pdsMgrPtr_;

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

  int numBlockRows_; // size of the expansion (number of PCE coefs)
  int numQuadPoints_;
  UQ::SampleType  sampleType_;

  int userSeed_;
  bool userSeedGiven_;

  std::string hackOutputFormat_;
  bool hackOutputCalledBefore_;
  bool hackOutputCalledBefore2_;
  bool hackOutputAllSamples_;
  bool outputtersCalledBefore_;
  bool outputSampleStats_;

  bool coefsOuterLoop_;

#ifdef Xyce_STOKHOS_ENABLE
  int PCEorder_;

  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases; 
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis;

  // Quadrature method
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quadMethod;

  // Triple product tensor
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;

  // Expansion method
  Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > expnMethod;

  bool resamplePCE_;
  bool outputPCECoeffs_;

  int numResamples_;
  bool numResamplesGiven_;

  bool useSparseGrid_;

  Teuchos::RCP<Linear::Graph> pceGraph;
#endif

  bool stdOutputFlag_;
  int debugLevel_;
  bool outputStochasticMatrix_;
  int voltLimAlgorithm_;
  bool stdOutputFlag_xvec_;

  bool outputsGiven_;
  bool outputsSetup_;
  std::vector<UQ::outputFunctionData*> outFuncDataVec_;

  bool measuresGiven_;
  std::vector<UQ::outputFunctionData*> measFuncDataVec_;

  bool outFuncGIDsetup_;

  int outputIndex_;
  int outputIndex2_;

  // non-block objects:
  Linear::Vector * nextSolutionPtr_;
  Linear::Vector * nextStatePtr_;
  Linear::Vector * nextStorePtr_;
 
  bool resetForStepCalledBefore_;
  bool useExpressionSamples_;
};

bool registerPCEFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // stokhos

#endif // Xyce_N_ANP_PCE_h

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
// Purpose        : Sampling analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 8/14/2017
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Sampling_h
#define Xyce_N_ANP_Sampling_h

#include <N_ANP_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>

#include <N_ANP_UQ_fwd.h>
#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>

#include <N_UTL_ExpressionData.h>
#include <N_UTL_MachDepParams.h>

#ifdef Xyce_STOKHOS_ENABLE
// make sure linking against the correct trilinos!
#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"
#endif

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : Sampling
// Purpose       : Sampling analysis class
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-------------------------------------------------------------------------
class Sampling : public AnalysisBase
{
public:
  Sampling(AnalysisManager &analysis_manager, Loader::Loader &loader, 
  Topo::Topology & topology, AnalysisBase &child_analysis);

  virtual ~Sampling();

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setSamplingOptions(const Util::OptionBlock & option_block);

  const TimeIntg::TIAParams &getTIAParams() const; // override
  TimeIntg::TIAParams &getTIAParams(); // override

  virtual bool getDCOPFlag() const;

protected:
  virtual bool doRun();
  virtual bool doInit();
  virtual bool doLoopProcess();
  virtual bool doProcessSuccessfulStep();
  virtual bool doProcessFailedStep();
  virtual bool doFinish();
  virtual bool doHandlePredictor() { return true; }

  void updateEnsembleOutputs();
  void completeEnsembleOutputs();
  void hackEnsembleOutput();

private:
  AnalysisManager &     analysisManager_;
  Loader::Loader &      loader_;
  Topo::Topology &      topology_;
  OutputMgrAdapter &    outputManagerAdapter_;
  IO::Measure::Manager &    measureManager_;
  AnalysisBase &        childAnalysis_;
  SweepVector           samplingVector_;
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

  int numSamples_;
  bool numSamplesGiven_;
  UQ::SampleType  sampleType_;

  int userSeed_;
  bool userSeedGiven_;

  std::string hackOutputFormat_;

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

  bool outputSampleStats_;

  bool useExpressionSamples_;
};

bool registerSamplingFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_Sampling_h

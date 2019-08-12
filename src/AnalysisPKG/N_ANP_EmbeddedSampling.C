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
// Purpose       : Embedded Sampling class analysis functions.
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#if __cplusplus>=201103L
#else
#include <sstream>
#endif

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_EmbeddedSampling.h>
#include <N_ANP_StepEvent.h>
#include <N_ANP_UQSupport.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_Transient.h>

#include <N_DEV_DeviceMgr.h>

#include <N_ERH_Message.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_MeasureManager.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>
#include <N_UTL_MachDepParams.h>
#include <N_LOA_CktLoader.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>

// new Embedded Sampler loader class
#include <N_LOA_ESLoader.h>

#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_System.h>
// new Embedded sampling headers: (patterned, kind of, off the HB headers)  ES = Embedded Sampling
#include <N_LAS_ESBuilder.h>
#include <N_LAS_ESSolverFactory.h>
#include <N_LAS_TranSolverFactory.h>

#include <N_PDS_fwd.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_LAPACK.hpp>


#if Xyce_STOKHOS_ENABLE
#include <Sacado_No_Kokkos.hpp>
#include <Stokhos_Sacado.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#endif

// ---------- Standard Includes ----------
#if( defined HAVE__ISNAN_AND__FINITE_SUPPORT )
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#else
#define isnan(x) std::isnan(x)
#define isinf(x) std::isinf(x)
#endif

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::EmbeddedSampling
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
EmbeddedSampling::EmbeddedSampling(
    AnalysisManager &analysis_manager, 
    Linear::System & linear_system,
    Nonlinear::Manager & nonlinear_manager,
    Device::DeviceMgr & device_manager,
    Linear::Builder & builder,
    Loader::Loader &loader, 
    Topo::Topology & topology, 
    IO::InitialConditionsManager & initial_conditions_manager,
    AnalysisBase &child_analysis)
    : AnalysisBase(analysis_manager, "EmbeddedSampling"),
      analysisManager_(analysis_manager),
      loader_(loader),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      deviceManager_(device_manager),
      builder_(builder),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
      measureManager_(outputManagerAdapter_.getMeasureManager()),
      childAnalysis_(child_analysis),

      esLoaderPtr_(0),
      esBuilderPtr_(0),
      esLinearSystem_(0),
      solverFactory_(0),

      samplingVector_(),
      maxParamStringSize_(0),
      lower_bounds_Given_(false),
      upper_bounds_Given_(false),
      covMatrixGiven_(false),
      numSamples_(1),
      numSamplesGiven_(false),
      sampleType_(UQ::MC),
      userSeed_(0),
      userSeedGiven_(false),
      hackOutputFormat_("TECPLOT"),
      hackOutputCalledBefore_(false),
      hackOutputAllSamples_(false),
      outputtersCalledBefore_(false),
      outputSampleStats_(true),
#if Xyce_STOKHOS_ENABLE
      regressionPCEenable_(false),
      projectionPCEenable_(false),
      PCEorder_(4),
      resamplePCE_(false),
      outputPCECoeffs_(false),
      regressionPCEcoeffs_(),
      projectionPCEcoeffs_(),
      numResamples_(1000),
      numResamplesGiven_(false),
      useSparseGrid_(false),
#endif
      stdOutputFlag_(false),
      debugLevel_(0),
      outputsGiven_(false),
      outputsSetup_(false),
      measuresGiven_(false),
      outFuncGIDsetup_(false),
      outputIndex_(0),
      resetForStepCalledBefore_(false)
{
  pdsMgrPtr_ = analysisManager_.getPDSManager();

  nextSolutionPtr_ = builder_.createVector();
  nextStatePtr_ = builder_.createStateVector();
  nextStorePtr_ = builder_.createStoreVector();
}


//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::~EmbeddedSampling
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
EmbeddedSampling::~EmbeddedSampling()
{
  if (nextSolutionPtr_)
  {
    delete nextSolutionPtr_;
  }
  if (nextStatePtr_)
  {
    delete nextStatePtr_;
  }
  if (nextStorePtr_)
  {
    delete nextStorePtr_;
  }

  if (esLoaderPtr_)
  {
    delete esLoaderPtr_;
  }

  if (esLinearSystem_)
  {
    delete esLinearSystem_;
  }

  if (solverFactory_)
  {
    delete solverFactory_;
  }

  covMatrix_.clear();

  {
  std::vector<UQ::outputFunctionData*>::iterator iter = outFuncDataVec_.begin();
  std::vector<UQ::outputFunctionData*>::iterator end = outFuncDataVec_.end();
  for (; iter!=end; ++iter)
  {
    if ( (*iter) )
    {
      delete (*iter);
    }
  }
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void EmbeddedSampling::notify(const StepEvent &event) 
{
  if (event.state_ == StepEvent::STEP_STARTED)
  {
    AnalysisBase::resetForStepAnalysis();

    if (resetForStepCalledBefore_)
    {
      analysisManager_.setNextOutputTime(0.0);

      nonlinearManager_.resetAll(Nonlinear::DC_OP);
      nonlinearManager_.setMatrixFreeFlag( false );
      nonlinearManager_.setLinSolOptions( saved_lsOB_ );
      nonlinearManager_.registerSolverFactory( NULL );

      analysisManager_.initializeSolverSystem(
          getTIAParams(), loader_, linearSystem_, nonlinearManager_, deviceManager_);

      deviceManager_.initializeAll(linearSystem_);

      nonlinearManager_.initializeAll(
        analysisManager_,
        analysisManager_.getNonlinearEquationLoader(),
        linearSystem_,
        *analysisManager_.getDataStore(),
        *analysisManager_.getPDSManager(),
        initialConditionsManager_,
        analysisManager_.getOutputManagerAdapter().getOutputManager(),
        topology_);

      TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
      esLoaderPtr_->loadDeviceErrorWeightMask(dsPtr->deviceErrorWeightMask_);

      analysisManager_.getXyceTranTimer().resetStartTime();
    }

    resetForStepCalledBefore_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  bool bsuccess = true;

  bool paramGiven=false;
  bool typeGiven=false;

  bool meanGiven=false;
  bool stdDevGiven=false;

  lower_bounds_Given_=false; // class variables
  upper_bounds_Given_=false; // class variables 

  bool alphaGiven=false;
  bool betaGiven=false;

  bool normalSpecified=false;
  bool uniformSpecified=false;
  bool gammaSpecified=false;

  Util::ParamList::const_iterator iter = paramsBlock.begin();
  Util::ParamList::const_iterator end   = paramsBlock.end();

  for ( ; iter != end; ++ iter)
  {
    if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      paramGiven=true;
      ExtendedString tag = iter->stringValue();
      tag.toUpper();
      paramNameVec_.push_back(tag);
      if (tag.size() > maxParamStringSize_) { maxParamStringSize_ = tag.size(); }
    }
    else if (std::string( iter->uTag() ,0,4) == "TYPE")
    {
      typeGiven=true;
      ExtendedString tag = iter->stringValue();
      tag.toUpper();
      typeVec_.push_back(tag);


      if(tag=="NORMAL")  normalSpecified=true;
      if(tag=="UNIFORM") uniformSpecified=true;
      if(tag=="GAMMA")   gammaSpecified=true;
    }
    else if (std::string( iter->uTag() ,0,5) == "MEANS")
    {
      meanGiven=true;
      double mean = iter->getImmutableValue<double>();
      meanVec_.push_back(mean);
    }
    else if (std::string( iter->uTag() ,0,14) == "STD_DEVIATIONS")
    {
      stdDevGiven=true;
      double stdDev = iter->getImmutableValue<double>();
      stdDevVec_.push_back(stdDev);
    }
    else if (std::string( iter->uTag() ,0,12) == "LOWER_BOUNDS")
    {
      lower_bounds_Given_=true;
      double start = iter->getImmutableValue<double>();
      lower_bounds_Vec_.push_back(start);
    }
    else if (std::string( iter->uTag() ,0,12) == "UPPER_BOUNDS")
    {
      upper_bounds_Given_=true;
      double stop = iter->getImmutableValue<double>();
      upper_bounds_Vec_.push_back(stop);
    }
    else if (std::string( iter->uTag() ,0,5) == "ALPHA")
    {
      alphaGiven=true;
      double alpha = iter->getImmutableValue<double>();
      alphaVec_.push_back(alpha);
    }
    else if (std::string( iter->uTag() ,0,4) == "BETA")
    {
      betaGiven=true;
      double beta = iter->getImmutableValue<double>();
      betaVec_.push_back(beta);
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() 
        << " is not a recognized sampling option.\n" << std::endl;
    }
  }


  int paramSize = paramNameVec_.size();
  int typeSize = typeVec_.size();

  int meanSize = meanVec_.size();
  int stdDevSize = stdDevVec_.size();

  int lower_bounds_Size = lower_bounds_Vec_.size();
  int upper_bounds_Size = upper_bounds_Vec_.size();

  int alphaSize = alphaVec_.size();
  int betaSize = betaVec_.size();

  if ( paramSize != typeSize ) { Report::DevelFatal() <<  "parameter and type arrays must be equal sizes."; }

  if (meanGiven || normalSpecified) { if (paramSize != meanSize) { Report::DevelFatal() <<  "parameter and mean arrays must be equal sizes. parameter size = " << paramSize << " means size = " << meanSize ; } }
  if (stdDevGiven || normalSpecified) { if (paramSize != stdDevSize) { Report::DevelFatal() <<  "parameter and stdDev arrays must be equal sizes."; } }
  if (lower_bounds_Given_ || uniformSpecified) { if (paramSize != lower_bounds_Size) { Report::DevelFatal() <<  "parameter and lower bounds arrays must be equal sizes."; } }
  if (upper_bounds_Given_ || uniformSpecified) { if (paramSize != upper_bounds_Size) { Report::DevelFatal() <<  "parameter and upper bounds arrays must be equal sizes."; } }
  if (alphaGiven || gammaSpecified) { if (paramSize != alphaSize) { Report::DevelFatal() <<  "parameter and alpha arrays must be equal sizes."; } }
  if (betaGiven || gammaSpecified) { if (paramSize != betaSize) { Report::DevelFatal() <<  "parameter and beta arrays must be equal sizes."; } }


  // check if the lower bounds are always < upper bounds, assuming both given
  if (lower_bounds_Given_ && upper_bounds_Given_)
  {
    for (int ibound=0;ibound<lower_bounds_Size;++ibound)
    {
      if(lower_bounds_Vec_[ibound] >= upper_bounds_Vec_[ibound])
      {
        Report::DevelFatal() <<  paramNameVec_[ibound] << " lower_bounds must be smaller than upper_bounds.";
      }
    }
  }
  else // if they are NOT given, then this is a problem for uniform distributions.
  { 
    // how to check this?  Currently this gets muddled if each param has a different type of distribution.  This probably needs a refactor.
  }

  // now put all this information into the sampling vector.
  for (int ip=0;ip<paramSize;++ip)
  {
    SweepParam sampling_param;
    sampling_param.type = typeVec_[ip]; // type = normal, uniform, etc
    sampling_param.name = paramNameVec_[ip]; // param name

    if (sampling_param.type == "UNIFORM")
    {
      sampling_param.startVal = lower_bounds_Vec_[ip];
      sampling_param.stopVal  = upper_bounds_Vec_[ip];
    }
    else if (sampling_param.type == "NORMAL") 
    {
      sampling_param.mean     = meanVec_[ip];
      sampling_param.stdDev   = stdDevVec_[ip];

      if ( !(lower_bounds_Vec_.empty()) )
      {
        sampling_param.lower_bound = lower_bounds_Vec_[ip];
        sampling_param.lower_boundGiven = true;

      }

      if ( !(upper_bounds_Vec_.empty()) )
      {
        sampling_param.upper_bound = upper_bounds_Vec_[ip];
        sampling_param.upper_boundGiven = true;
      }
    }
#if __cplusplus>=201103L
    else if (sampling_param.type == "GAMMA") 
    {
      sampling_param.alpha    = alphaVec_[ip];
      sampling_param.beta     = betaVec_[ip];

      if ( !(lower_bounds_Vec_.empty()) )
      {
        sampling_param.lower_bound = lower_bounds_Vec_[ip];
        sampling_param.lower_boundGiven = true;

      }

      if ( !(upper_bounds_Vec_.empty()) )
      {
        sampling_param.upper_bound = upper_bounds_Vec_[ip];
        sampling_param.upper_boundGiven = true;
      }
    }
#endif    
    else
    {
      Report::DevelFatal().in("parseEmbeddedSamplingParam") << "Unsupported SAMPLING type";
    }
    samplingVector_.push_back(sampling_param);

  }

  outputManagerAdapter_.setStepSweepVector(samplingVector_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::setEmbeddedSamplingOptions
// Purpose       :
// Special Notes : These are from '.options EMBEDDEDSAMPLES'
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::setEmbeddedSamplingOptions(const Util::OptionBlock & option_block)
{
  bool bsuccess = true;
  numSamplesGiven_ = false;
  Util::ParamList::const_iterator it  = option_block.begin();
  Util::ParamList::const_iterator end = option_block.end();
  for ( ; it != end; ++ it)
  {
    if ((*it).uTag() == "NUMSAMPLES")
    {
      numSamples_ = (*it).getImmutableValue<int>();
      numSamplesGiven_ = true;
      if (numSamples_ <= 0)
        Report::UserError() << "NUMSAMPLES parameter on .EMBEDDEDSAMPLES line must > 0";
    }
    else if (std::string((*it).uTag() ,0,9) == "COVMATRIX" ) // this is a vector
    {
      covMatrixGiven_ = true;
      covMatrix_.push_back( (*it).getImmutableValue<double>() );
    }
    else if ((*it).uTag() == "OUTPUTFORMAT" ) 
    {
      ExtendedString tag = (*it).stringValue();
      hackOutputFormat_ = tag.toUpper();
    }
    else if ((*it).uTag() == "OUTPUTALLSAMPLES")
    {
      hackOutputAllSamples_=static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTSAMPLESTATS")
    {
      outputSampleStats_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
#if Xyce_STOKHOS_ENABLE
    else if ((*it).uTag() == "REGRESSION_PCE")
    {
      regressionPCEenable_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "PROJECTION_PCE")
    {
      projectionPCEenable_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "ORDER")
    {
      PCEorder_ = (*it).getImmutableValue<int>();
      if (PCEorder_ < 0)
       Report::UserError() << "ORDER parameter on .EMBEDDEDSAMPLES line must >= 0";
    }
#endif
    else if ((*it).uTag() == "SAMPLE_TYPE")
    {
      if ((*it).isNumeric())
      {
        int tmp = (*it).getImmutableValue<int>();
        if (tmp==0)
        {
          sampleType_ = UQ::MC;
        }
        else if (tmp==1)
        {
          sampleType_ = UQ::LHS;
        }
        else
        {
          Xyce::Report::UserWarning() << (*it).uTag() 
            << " = " << tmp << " is not a recognized sampling option.  Setting " << (*it).uTag() << " = RANDOM.\n" << std::endl;
          sampleType_ = UQ::MC;
        }
      }
      else
      {
        ExtendedString p((*it).stringValue());
        p.toUpper();
        if (p.substr(0,2) == "RANDOM")
        {
          sampleType_ = UQ::MC;
        }
        else if (p.substr(0,3) == "LHS")
        {
          sampleType_ = UQ::LHS;
        }
        else
        {
          Xyce::Report::UserWarning() << (*it).uTag() 
            << " = " << p << " is not a recognized sampling option.  Setting " << (*it).uTag() << " = RANDOM.\n" << std::endl;
          sampleType_ = UQ::MC;
        }
      }
    }
    else if ((*it).uTag() == "SEED")
    {
      userSeed_ = (*it).getImmutableValue<int>();
      userSeedGiven_ = true;
    }
#if Xyce_STOKHOS_ENABLE
    else if ((*it).uTag() == "RESAMPLE")
    {
      resamplePCE_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUT_PCE_COEFFS")
    {
      outputPCECoeffs_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "SPARSE_GRID")
    {
      useSparseGrid_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
#endif
    else if ((*it).uTag() == "STDOUTPUT")
    {
      stdOutputFlag_ =
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DEBUGLEVEL")
    {
      if (DEBUG_ES)
      {
        debugLevel_ = (*it).getImmutableValue<int>();
      }
    }
    else if (std::string((*it).uTag() ,0,7) == "OUTPUTS" )// this is a vector of expressions/solution variables
    {
      outputsGiven_ = true;
      UQ::outputFunctionData * ofDataPtr = new UQ::outputFunctionData();
      ofDataPtr->outFuncString = (*it).stringValue();
      outFuncDataVec_.push_back(ofDataPtr);
    }
    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() 
        << " is not a recognized sampling option, or may only be supported for PCE\n" << std::endl;
    }
  }

  // parse the expression now, so if there are any errors, they will come
  // up early in the simulation.
  if (outputsGiven_)
  {
    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      outFuncDataVec_[iout]->expDataPtr = new Util::ExpressionData(outFuncDataVec_[iout]->outFuncString);
    }
  }
  else if (measuresGiven_)
  {
    for (int iout=0;iout<measFuncDataVec_.size();++iout)
    {
      measFuncDataVec_[iout]->measureResponseFound 
        = measureManager_.find(measFuncDataVec_[iout]->outFuncString);

      Report::UserWarning0() << "Measure response " << measFuncDataVec_[iout]->outFuncString << " was not found.";
    }
    Report::UserWarning0() << "Output function was not specified";
  }
  else // give a warning.
  {
    Report::UserWarning0() << "Output function was not specified";
  }

  // signal the OutputMgr that Embedded Sampling is enabled
  outputManagerAdapter_.setEnableEmbeddedSamplingFlag(true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::getTIAParams() const
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
const TimeIntg::TIAParams & EmbeddedSampling::getTIAParams() const
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::getTIAParams()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
TimeIntg::TIAParams & EmbeddedSampling::getTIAParams()
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::getDCOPFlag()
// Purpose       :
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::getDCOPFlag() const
{
  return childAnalysis_.getDCOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::stepCallBack
// Purpose       : 
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void EmbeddedSampling::stepCallBack ()
{
  computeEnsembleOutputs();

#if Xyce_STOKHOS_ENABLE
  if (regressionPCEenable_) 
  {
    const int numParams = paramNameVec_.size();
    std::vector< std::vector<double> > x;
    UQ::unScaleSampleValues(numSamples_, samplingVector_, covMatrix_, meanVec_, Y_, x);

    std::vector< Stokhos::OrthogPolyApprox<int,double> > pceVec;

    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

      Stokhos::OrthogPolyApprox<int,double> & regressionPCE = outFunc.regressionPCE;
      regressionPCE.reset(regrBasis);

      std::vector<double> & f = outFunc.sampleOutputs;
      UQ::solveRegressionPCE( paramNameVec_.size(), PCEorder_, x, f, regressionPCE);

      if (stdOutputFlag_)
      {
        if (outputPCECoeffs_)
        {
          Xyce::lout() << "PCE coefs for " << outFunc.outFuncString << " = " << regressionPCE << std::endl;
          regressionPCE.print(Xyce::lout());
        }

        Xyce::lout() << std::endl;
        Xyce::lout() << "(embedded sampling) regression PCE mean of " << outFunc.outFuncString << " = " << regressionPCE.mean() << std::endl;
        Xyce::lout() << "(embedded sampling) regression PCE stddev of " << outFunc.outFuncString << " = " << regressionPCE.standard_deviation() << std::endl;
      }

      if (resamplePCE_)
      {
        pceVec.push_back(regressionPCE);
      }
    }

    if (resamplePCE_)
    {
      std::vector < std::vector<double> > fvec (outFuncDataVec_.size());
      std::vector <UQ::statisticalMoments> statVec (outFuncDataVec_.size());

      long theSeed = UQ::getTheSeed( analysisManager_.getComm(), analysisManager_.getCommandLine(), userSeed_, userSeedGiven_);
      UQ::sampleApproximationPCE(theSeed, sampleType_, samplingVector_, covMatrix_, meanVec_, numResamples_, numParams, pceVec, fvec, statVec);

      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
        Xyce::lout() << std::endl;
        Xyce::lout() << "Statistics from re-sampling regression PCE approximation of " << outFunc.outFuncString << ":" << std::endl;;
        Xyce::lout() << statVec[iout];
        Xyce::lout() << std::endl;
      }
    }
  }

//
  if (projectionPCEenable_)
  {
    std::vector< Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > > pceVec;

    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

      Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = outFunc.projectionPCE;
      projectionPCE.init(0.0);
      projectionPCE.reset(quadExpn);

      std::vector<double> & f = outFunc.sampleOutputs;

      UQ::solveProjectionPCE(quadBasis, quadMethod, f, projectionPCE);

      if (stdOutputFlag_)
      {
        if (outputPCECoeffs_)
        {
          Xyce::lout() << "PCE coefs for " << outFunc.outFuncString << " = " << projectionPCE << std::endl;
          projectionPCE.print(Xyce::lout());
        }

        Xyce::lout() << std::endl;
        Xyce::lout() << "(embedded sampling) projection PCE mean of " << outFunc.outFuncString << " = " << projectionPCE.mean() << std::endl;
        Xyce::lout() << "(embedded sampling) projection PCE stddev of " << outFunc.outFuncString << " = " << projectionPCE.standard_deviation() << std::endl;
      }

      if (resamplePCE_)
      {
        pceVec.push_back(projectionPCE);
      }
    }

    if (resamplePCE_)
    {
      std::vector < std::vector<double> > fvec (outFuncDataVec_.size());
      std::vector <UQ::statisticalMoments> statVec (outFuncDataVec_.size());

      long theSeed = UQ::getTheSeed( analysisManager_.getComm(), analysisManager_.getCommandLine(), userSeed_, userSeedGiven_);
      const int numParams = paramNameVec_.size();
      UQ::sampleApproximationPCE(theSeed, sampleType_, samplingVector_, covMatrix_, meanVec_, numResamples_, numParams, pceVec, fvec, statVec);

      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
        Xyce::lout() << std::endl;
        Xyce::lout() << "Statistics from re-sampling projection PCE approximation of " << outFunc.outFuncString << ":" << std::endl;;
        Xyce::lout() << statVec[iout];
        Xyce::lout() << std::endl;
      }
    }
  }
#endif
//
  // only call outputter functions if OUTPUTS= was used on the
  // .OPTIONS EMBEDDEDSAMPLES line.
  if (outputsGiven_)
  {
#if Xyce_STOKHOS_ENABLE
    if (!outputtersCalledBefore_)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[0]);
      if (regressionPCEenable_)
      {
        Stokhos::OrthogPolyApprox<int,double> & regressionPCE = outFunc.regressionPCE;
        int NN=regressionPCE.size();
        for (int ii=0;ii<NN;ii++)
        {
          const Stokhos::MultiIndex<int>& trm = regrBasis->term(ii);
          std::string coefString = "_coef(";
          for (int jj=0; jj< trm.size()-1; jj++)
            coefString += std::to_string(trm[jj]) + ",";
          coefString += std::to_string(trm[trm.size()-1]) + ")";

          regressionPCEcoeffs_.push_back(coefString);
        }
      }

      if (projectionPCEenable_)
      {
        Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = outFunc.projectionPCE;
        int NN=projectionPCE.size();
        std::vector<std::string> coeffs;
        for (int ii=0;ii<NN;ii++)
        {
          const Stokhos::MultiIndex<int>& trm = quadBasis->term(ii);
          std::string coefString = "_coef(";
          for (int jj=0; jj< trm.size()-1; jj++)
            coefString += std::to_string(trm[jj]) + ",";
          coefString += std::to_string(trm[trm.size()-1]) + ")";

          projectionPCEcoeffs_.push_back(coefString);
        }
      }
    }
#endif

#if Xyce_STOKHOS_ENABLE
    outputManagerAdapter_.outputEmbeddedSampling(regressionPCEenable_, projectionPCEenable_,
                 numSamples_, regressionPCEcoeffs_, projectionPCEcoeffs_, outFuncDataVec_);
#else
    // regressionPCEcoeffs_ and projectionPCEcoeffs_ vectors aren't defined for
    // non-Stokhos builds.  So, make two empty vectors.
    std::vector<std::string> emptyVec1, emptyVec2;
    outputManagerAdapter_.outputEmbeddedSampling(0, 0, numSamples_, emptyVec1,
                 emptyVec2, outFuncDataVec_);
#endif

    hackEnsembleOutput ();

    outputtersCalledBefore_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::doRun()
// Purpose       : This is the main controlling loop for sampling analysis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::doInit()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::doInit()
{
  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "EmbeddedSampling::init" << std::endl;
  }

  UQ::checkParameterList(
      analysisManager_.getComm(), 
      loader_, 
      samplingVector_.begin(), 
      samplingVector_.end());

#if Xyce_STOKHOS_ENABLE
  if( !numSamplesGiven_  && !projectionPCEenable_)
  {
    Report::UserFatal0() << "Number of samples not specified, and quadrature PCE isn't being used (quadrature PCE is only method that doesn't need it)";
  }
#else
  if( !numSamplesGiven_ )
  {
    Report::UserFatal0() << "Number of samples not specified";
  }
#endif
  else
  {
    // Deal with the random number seed, and set up random samples.
    // Don't bother with this is projection PCE has been specified.
    long theSeed = UQ::getTheSeed(
        analysisManager_.getComm(), 
        analysisManager_.getCommandLine(), userSeed_, userSeedGiven_);

    UQ::setupSampleValues(theSeed, sampleType_,
        numSamples_, samplingVector_, covMatrix_, meanVec_, X_, Y_);
  }

#if Xyce_STOKHOS_ENABLE
  // determine the samples. 
  //
  // non-intrusive spectral projection(NISP) samples are determined 
  // by the quadrature points
  if (projectionPCEenable_) 
  {
    const int d = paramNameVec_.size();
    const int p = PCEorder_;
    quadBases.resize(d); 
    for (int i=0; i<d; i++)
    {
      SweepParam & sp = samplingVector_[i];

      if (sp.type == "UNIFORM")
      {
        quadBases[i] = rcp(new Stokhos::LegendreBasis<int,double>(p));
      }
      else if (sp.type == "NORMAL") 
      {
        quadBases[i] = rcp(new Stokhos::HermiteBasis<int,double>(p,true));
      }
      else
      {
        Report::UserFatal0() << "Polynomial Chaos only works for normal and uniform distributions.";
      }
    }

    if (useSparseGrid_)
    {
      // ERK.  Understand these better
      double drop = 1.0e-12;
      const Stokhos::TotalOrderIndexSet<int> index_set(d, p);
      typedef Stokhos::TotalOrderLess< Stokhos::MultiIndex<int> > total_less;
      typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > lexo_less;

      quadBasis = rcp(new Stokhos::SmolyakBasis<int,double,total_less>( quadBases, index_set, drop));
      quadMethod = rcp(new Stokhos::SmolyakSparseGridQuadrature<int,double>(quadBasis, index_set));
    }
    else
    {
      quadBasis = rcp(new Stokhos::CompletePolynomialBasis<int,double>(quadBases));
      quadMethod = rcp(new Stokhos::TensorProductQuadrature<int,double>(quadBasis));
    }

    quadCijk = quadBasis->computeTripleProductTensor();
    quadExpn = rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(quadBasis, quadCijk, quadMethod));

    UQ::setupPCEQuadPoints ( quadBasis, quadMethod, quadExpn, samplingVector_, covMatrix_, meanVec_, X_, Y_);
    numSamples_ = quadMethod->size();
  }

  // if regression PCE specified, then create the basis, 
  // determine the basis size, and then force # of samples to be larger
  // than the basis size so that the least squares solve will be
  // successful.
  if (regressionPCEenable_) 
  {
    const int d = paramNameVec_.size();
    const int p = PCEorder_;
    regrBases.resize(d); 
    for (int i=0; i<d; i++)
    {
      SweepParam & sp = samplingVector_[i];

      if (sp.type == "UNIFORM")
      {
        regrBases[i] = rcp(new Stokhos::LegendreBasis<int,double>(p));
      }
      else if (sp.type == "NORMAL") 
      {
        regrBases[i] = rcp(new Stokhos::HermiteBasis<int,double>(p,true));
      }
      else
      {
        Report::UserFatal0() << "Polynomial Chaos only works for normal and uniform distributions.";
      }
    }

    regrBasis = rcp(new Stokhos::CompletePolynomialBasis<int,double>(regrBases));
    int basisSize = regrBasis->size();

    if (numSamples_ < basisSize)
    {
      Report::UserFatal0()
            << "Number of samples = " << numSamples_ << ", which is smaller than the basis size = " 
            << basisSize << ".  Increase the number of samples to use this basis";
    }
  }
#endif

  setupBlockSystemObjects ();
  setupEnsembles ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::setupBlockSystemObjects 
// Purpose       :
// 
// set things up related to the linear system and solver.  The goal here
// is to replace the linear solver objects used by the underlying child
// analysis with block versions, where the # of blocks equals the number
// of samples.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
void  EmbeddedSampling::setupBlockSystemObjects ()
{
  analysisManager_.resetSolverSystem();

  esBuilderPtr_ = rcp(new Linear::ESBuilder(numSamples_));

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << "EmbeddedSampling::setupBlockSystemObjects():  Generate Maps,etc\n";
  }
  {
    Stats::StatTop _setupStepStat("Setup Maps/Graphs");
    Stats::TimeBlock _setupStepTimer(_setupStepStat);

    esBuilderPtr_->generateMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::SOLUTION ), false),
                             rcp(pdsMgrPtr_->getParallelMap( Parallel::SOLUTION_OVERLAP_GND ), false) );
    esBuilderPtr_->generateStateMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::STATE ),false) );
    esBuilderPtr_->generateStoreMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::STORE ),false) );
    esBuilderPtr_->generateLeadCurrentMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::LEADCURRENT ),false) );
    esBuilderPtr_->generateGraphs( *pdsMgrPtr_->getMatrixGraph( Parallel::JACOBIAN ));
  }

  // Create ES Loader.
  delete esLoaderPtr_;
  esLoaderPtr_ = new Loader::ESLoader(deviceManager_, builder_, numSamples_, samplingVector_, Y_);
  esLoaderPtr_->registerESBuilder(esBuilderPtr_);
  esLoaderPtr_->registerAppLoader( rcp(&loader_, false) );
  //-----------------------------------------

  //Finish setup of ES Loader
  //-----------------------------------------
  {
    //-----------------------------------------
    //Construct Solvers, etc.
    //-----------------------------------------
    delete esLinearSystem_;
    esLinearSystem_ = new Linear::System();
    //-----------------------------------------

    //hack needed by TIA initialization currently
    esBuilderPtr_->registerPDSManager( pdsMgrPtr_ );
    esLinearSystem_->registerPDSManager( pdsMgrPtr_ );
    esLinearSystem_->registerBuilder( &*esBuilderPtr_ );

    // the linear system class will only create a "linear problem" object if the 
    // builder object (in this case ES builder) will create a non-NULL matrix pointer.
    esLinearSystem_->initializeSystem();

    nonlinearManager_.setLinSolOptions( saved_lsESOB_ );
    nonlinearManager_.setMatrixFreeFlag( false ); // NOT TRUE for ES!

    if (!solverFactory_)
    {
      // Generate the ES solver factory.
      //solverFactory_ = new Linear::ESSolverFactory( builder_ );
      solverFactory_ = new Linear::TranSolverFactory();
    }

#if 0
    // Register application loader with solver factory
    solverFactory_->registerESLoader( rcp(esLoaderPtr_, false) );
    solverFactory_->registerESBuilder( esBuilderPtr_ );
#endif

#if 0
    // ERK:  I have not created a preconditioner factory
    if (!precFactory_)
    {
      // Generate the ES preconditioner factory.
      precFactory_ = new Linear::ESPrecondFactory(saved_lsESOB_, builder_);
    }

    // Register application loader with preconditioner factory
    precFactory_->registerESLoader( rcp(esLoaderPtr_, false) );
    precFactory_->registerESBuilder( esBuilderPtr_ );
    precFactory_->setFastTimes( goodTimePoints_ );
    precFactory_->setTimeSteps( timeSteps_ );
    precFactory_->setESFreqs( freqPoints_ );
#endif

    nonlinearManager_.registerSolverFactory( solverFactory_ );
    //nonlinearManager_.registerPrecondFactory( precFactory_ );

  //Initialization of Solvers, etc. 
    analysisManager_.initializeSolverSystem(
        getTIAParams(),
        *esLoaderPtr_,
        *esLinearSystem_,
        nonlinearManager_,
        deviceManager_);

    nonlinearManager_.initializeAll(
        analysisManager_, 
        analysisManager_.getNonlinearEquationLoader(), 
        *esLinearSystem_,
        *analysisManager_.getDataStore(),
        *analysisManager_.getPDSManager(),
        initialConditionsManager_,
        analysisManager_.getOutputManagerAdapter().getOutputManager(),
        topology_);

    childAnalysis_.registerLinearSystem(esLinearSystem_);
 
    // don't have a mode yet.  Do I need one??  The real mode should be 
    // the underlying analysis (tran, dc, etc), so probably not
    //nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_ES)); 
  }
  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  esLoaderPtr_->loadDeviceErrorWeightMask(dsPtr->deviceErrorWeightMask_);
  
  childAnalysis_.registerParentAnalysis(this);
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::setupEnsembles
// Purpose       : This function sets up the stokhos ensemble objects, using
//                 the computed sample values.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void EmbeddedSampling::setupEnsembles ()
{

}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::doLoopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::doLoopProcess()
{

  // The idea of this type of analysis is to perform "simultaneous propagation", 
  // in which the sampling is performed, but is essentially the inner loop rather 
  // than the outer loop of the simulation.
  //
  // This is being set up by considering the original analysis (DC or TRAN or whatever) 
  // to be the "child" analysis, but to then replace the linear objects with block 
  // versions and replace the loader with a block loader.  Otherwise, all the control 
  // is handled in the child process.
  Xyce::lout() << "***** Beginning Embedded Sampling (simultaneous propagation) simulation....\n" << std::endl;

#if Xyce_STOKHOS_ENABLE
  if (projectionPCEenable_)
  {
    Xyce::lout() << "***** Projection PCE enabled.  Number of quadrature points = " << numSamples_ << "\n" << std::endl;
  }
  else
#endif
  {
    Xyce::lout() << "***** Number of sample points = " << numSamples_ << "\n" << std::endl;
  }

  // test:
  analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

  // solve the loop.
  bool status = childAnalysis_.run();

  return status;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::doProcessSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::doProcessSuccessfulStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::doProcessFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::doProcessFailedStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::doFinish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::doFinish()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::computeEnsembleOutputs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
void EmbeddedSampling::computeEnsembleOutputs()
{
  Parallel::Machine comm = analysisManager_.getComm();

  if (outputsGiven_)
  {
    // first zero out the val, mean and stddev for this time/dc point.
    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      outFuncDataVec_[iout]->sampleOutputs.clear();
    }

    IO::OutputMgr & output_manager = outputManagerAdapter_.getOutputManager();
    double circuit_time = analysisManager_.getStepErrorControl().nextTime;
    double circuit_dt = analysisManager_.getStepErrorControl().currentTimeStep;

    // loop over the params
    TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
    Xyce::Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextSolutionPtr) );
    Xyce::Linear::BlockVector & bSta = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextStatePtr) );
    Xyce::Linear::BlockVector & bSto = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextStorePtr) );

    nextSolutionPtr_->putScalar(0.0);
    nextStatePtr_->putScalar(0.0);
    nextStorePtr_->putScalar(0.0);

    // setup
    if (!outputsSetup_)
    {
      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

        outFunc.expDataPtr->setup( comm,
                 output_manager.getOpBuilderManager(),
                 output_manager.getMainContextFunctionMap(),
                 output_manager.getMainContextParamMap(),
                 output_manager.getMainContextGlobalParamMap());
      }

      outputsSetup_ = true;
    }

    int BlockCount = bX.blockCount();
    for( int i = 0; i < BlockCount; ++i )
    {
      *nextSolutionPtr_ = bX.block(i);
      *nextStatePtr_ = bSta.block(i);
      *nextStorePtr_ = bSto.block(i);

      //get expression value and compute means, etc
      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

        double val = outFunc.expDataPtr->evaluate(comm, circuit_time, circuit_dt, 
              nextSolutionPtr_, nextStatePtr_, nextStorePtr_ );

        outFunc.sampleOutputs.push_back(val);
      }
    }

    // finish mean and std dev
    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
      //outFunc.completeStatistics(BlockCount);
      outFunc.completeStatistics();

      if (stdOutputFlag_ && outputSampleStats_)
      {
        // histrogram is a hack that doesn't work yet
        //UQ::histrogram(std::cout, outFunc.outFuncString, outFunc.sampleOutputs);

        std::string sampleTypeStr = "Embedded MC";
        if ( sampleType_ == UQ::MC)
        {
          sampleTypeStr = "Embedded MC";
        }
        else
        {
          sampleTypeStr = "Embedded LHS";
        }
        outFunc.output(Xyce::lout(), sampleTypeStr);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::hackEnsembleOutput ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
void EmbeddedSampling::hackEnsembleOutput ()
{
    std::string fileName;

    if (hackOutputFormat_=="TECPLOT")
    {
      fileName = analysisManager_.getNetlistFilename() + "_ensemble.dat";
    }
    else if (hackOutputFormat_=="STD")
    {
      fileName = analysisManager_.getNetlistFilename() + "_ensemble.prn";
    }
    else
    {
      Xyce::Report::UserWarning() << hackOutputFormat_ 
        << " is not a recognized ensemble output option.\n" << std::endl;
    }

  //
    std::ofstream output_stream(fileName.c_str(), 
            !hackOutputCalledBefore_ ? std::ios_base::out : std::ios_base::app);

    if (!hackOutputCalledBefore_)
    {
      // header output
      output_stream << "TITLE = \"embedded sampling output\"\tVARIABLES= "<<std::endl;

      //if ( analysisManager_.getTransientFlag() || analysisManager_.getTranOPFlag() ) // this doesn't work!
      //if ( childAnalysis_.isAnalysis(ANP_MODE_TRANSIENT) )
      if ( !(childAnalysis_.getDCOPFlag()) )
      {
        output_stream << "\t\" TIME \""<<std::endl;
      }
      else
      {
        output_stream << "\t\" INDEX \""<<std::endl;
      }

      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

        if (outputSampleStats_)
        {
          std::string meanString = outFunc.outFuncString + "_mean";
          std::string meanStringPlus = outFunc.outFuncString + "_meanPlus";
          std::string meanStringMinus = outFunc.outFuncString + "_meanMinus";

          std::string stddevString = outFunc.outFuncString + "_stddev";
          std::string varianceString = outFunc.outFuncString + "_variance";

          output_stream << "\t\" " << meanString << "\""<<std::endl;
          output_stream << "\t\" " << meanStringPlus << "\""<<std::endl;
          output_stream << "\t\" " << meanStringMinus << "\""<<std::endl;

          output_stream << "\t\" " << stddevString << "\""<<std::endl;
          output_stream << "\t\" " << varianceString << "\""<<std::endl;
        }

#if Xyce_STOKHOS_ENABLE
        if (regressionPCEenable_)
        {
          std::string meanString = outFunc.outFuncString + "_regr_pce_mean";
          std::string meanStringPlus = outFunc.outFuncString + "_regr_pce_meanPlus";
          std::string meanStringMinus = outFunc.outFuncString + "_regr_pce_meanMinus";

          std::string stddevString = outFunc.outFuncString + "_regr_pce_stddev";
          std::string varianceString = outFunc.outFuncString + "_regr_pce_variance";

          output_stream << "\t\" " << meanString << "\""<<std::endl;
          output_stream << "\t\" " << meanStringPlus << "\""<<std::endl;
          output_stream << "\t\" " << meanStringMinus << "\""<<std::endl;

          output_stream << "\t\" " << stddevString << "\""<<std::endl;
          output_stream << "\t\" " << varianceString << "\""<<std::endl;

          if (outputPCECoeffs_)
          {
            std::vector<std::string>::const_iterator it;
            for (it=regressionPCEcoeffs_.begin();it!=regressionPCEcoeffs_.end();++it)
            {
              output_stream << "\t\" " << outFunc.outFuncString + *it << "\""<<std::endl;
            }
          }
        }

        if (projectionPCEenable_)
        {
          std::string meanString = outFunc.outFuncString + "_quad_pce_mean";
          std::string meanStringPlus = outFunc.outFuncString + "_quad_pce_meanPlus";
          std::string meanStringMinus = outFunc.outFuncString + "_quad_pce_meanMinus";

          std::string stddevString = outFunc.outFuncString + "_quad_pce_stddev";
          std::string varianceString = outFunc.outFuncString + "_quad_pce_variance";

          output_stream << "\t\" " << meanString << "\""<<std::endl;
          output_stream << "\t\" " << meanStringPlus << "\""<<std::endl;
          output_stream << "\t\" " << meanStringMinus << "\""<<std::endl;

          output_stream << "\t\" " << stddevString << "\""<<std::endl;
          output_stream << "\t\" " << varianceString << "\""<<std::endl;

          if (outputPCECoeffs_)
          {
            std::vector<std::string>::const_iterator it;
            for (it=projectionPCEcoeffs_.begin();it!=projectionPCEcoeffs_.end();++it)
            {
              output_stream << "\t\" " << outFunc.outFuncString + *it << "\""<<std::endl;
            }
          }
        }
#endif
        if (hackOutputAllSamples_)
        {
          for(int i=0;i<numSamples_; ++i)
          {
#if __cplusplus>=201103L
            std::string sampleString = outFunc.outFuncString + "_"+std::to_string(i);
#else
            std::stringstream ss;
            ss << i;
            std::string sampleString = outFunc.outFuncString + "_"+ss.str();
#endif
            output_stream << "\t\" " << sampleString << "\""<<std::endl;
          }
        }
      }
   
      output_stream << "ZONE F=POINT  T=\"Xyce data\""<<std::endl;
      hackOutputCalledBefore_ = true;
    }

    // data output
    output_stream.setf(std::ios::scientific);

    //if ( analysisManager_.getTransientFlag() || analysisManager_.getTranOPFlag() ) // this doesn't work!
    //if ( childAnalysis_.isAnalysis(ANP_MODE_TRANSIENT) )
    if ( !(childAnalysis_.getDCOPFlag()) )
    {
      output_stream << analysisManager_.getStepErrorControl().currentTime;
    }
    else
    {
      output_stream << outputIndex_;
      outputIndex_++;
    }

    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

      if (outputSampleStats_)
      {
        output_stream << "\t" << outFunc.sm.mean;

        output_stream << "\t" << (outFunc.sm.mean+outFunc.sm.stddev);
        output_stream << "\t" << (outFunc.sm.mean-outFunc.sm.stddev);

        output_stream << "\t" << outFunc.sm.stddev;
        output_stream << "\t" << outFunc.sm.variance;
      }

#if Xyce_STOKHOS_ENABLE
      if (regressionPCEenable_)
      {
        Stokhos::OrthogPolyApprox<int,double> & regressionPCE = outFunc.regressionPCE;

        double pce_mean = regressionPCE.mean();
        double pce_stddev = regressionPCE.standard_deviation();
        double pce_variance = pce_stddev*pce_stddev;

        if ( isinf(pce_mean) || isnan(pce_mean) )
        {
          pce_mean = 0.0;
        }

        if ( isinf(pce_stddev) || isnan(pce_stddev) )
        {
          pce_stddev = 0.0;
        }

        if ( isinf(pce_variance) || isnan(pce_variance) )
        {
          pce_variance = 0.0;
        }

        output_stream << "\t" << pce_mean;

        output_stream << "\t" << (pce_mean+pce_stddev);
        output_stream << "\t" << (pce_mean-pce_stddev);

        output_stream << "\t" << pce_stddev;
        output_stream << "\t" << pce_variance;

        if (outputPCECoeffs_)
        {
          int NN=regressionPCE.size();
          for (int ii=0;ii<NN;ii++)
          {
            output_stream << "\t" << regressionPCE[ii];
          }
        }
      }

      if (projectionPCEenable_)
      {
        Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = outFunc.projectionPCE;

        double pce_mean = projectionPCE.mean();
        double pce_stddev = projectionPCE.standard_deviation();
        double pce_variance = pce_stddev*pce_stddev;

        if ( isinf(pce_mean) || isnan(pce_mean) )
        {
          pce_mean = 0.0;
        }

        if ( isinf(pce_stddev) || isnan(pce_stddev) )
        {
          pce_stddev = 0.0;
        }

        if ( isinf(pce_variance) || isnan(pce_variance) )
        {
          pce_variance = 0.0;
        }

        output_stream << "\t" << pce_mean;

        output_stream << "\t" << (pce_mean+pce_stddev);
        output_stream << "\t" << (pce_mean-pce_stddev);

        output_stream << "\t" << pce_stddev;
        output_stream << "\t" << pce_variance;

        if (outputPCECoeffs_)
        {
          int NN=projectionPCE.size();
          for (int ii=0;ii<NN;ii++)
          {
            //output_stream << "\t" << projectionPCE[ii];
            output_stream << "\t" << projectionPCE.fastAccessCoeff(ii);
          }
        }
      }
#endif

      // output individual samples
      if (hackOutputAllSamples_)
      {
        for(int i=0;i<numSamples_; ++i)
        {
          output_stream << "\t" << outFunc.sampleOutputs[i];
        }
      }
    }

    output_stream << std::endl;

  return;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSampling::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool EmbeddedSampling::setDCOptions(const Util::OptionBlock & paramsBlock)
{
  dcSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));

  if (DEBUG_ES)
  {
    if ( (paramsBlock.size() > 0) )
    {
      Xyce::Util::ParamList::const_iterator iter = paramsBlock.begin();
      int i=0;
      for ( ; iter != paramsBlock.end(); ++iter,++i)
      {
        const Xyce::Util::Param & par = *(iter);
        Xyce::lout() << "EmbeddedSampling::setDCOptions.  paramsBlock["<<i<<"] = " << par;
      }
    }
  }

  return true;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : EmbeddedSamplingFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
///
/// Factory for parsing Step parameters from the netlist and creating Step analysis.
///
class EmbeddedSamplingFactory : public Util::Factory<AnalysisBase, EmbeddedSampling>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : EmbeddedSamplingFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 5/26/2018
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the EmbeddedSampling analysis factory
  ///
  /// @invariant Stores the results of parsing.  Multiple EmbeddedSampling analysis options may be
  /// applied and each generates and additional step.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  ///
  EmbeddedSamplingFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Device::DeviceMgr &   device_manager,
    Linear::Builder &     builder,
    Loader::Loader &            loader,
    Topo::Topology &            topology,
    IO::InitialConditionsManager & initial_conditions_manager)
    : Util::Factory<AnalysisBase, EmbeddedSampling>(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      deviceManager_(device_manager),
      builder_(builder),
      loader_(loader),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager)
  {}

  virtual ~EmbeddedSamplingFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 5/26/2018
  //-----------------------------------------------------------------------------
  ///
  /// Create a new EmbeddedSampling analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new EmbeddedSampling analysis object
  ///
  EmbeddedSampling *create() const
  {
    // don't have a mode yet
    //analysisManager_.setAnalysisMode(ANP_MODE_ES);

    EmbeddedSampling *es = new EmbeddedSampling(
        analysisManager_, 
    linearSystem_, nonlinearManager_, deviceManager_, builder_, 
        loader_, 
        topology_, 
        initialConditionsManager_,
        analysisManager_.getAnalysisObject());

    for (std::vector<Util::OptionBlock>::const_iterator it = samplingSweepAnalysisOptionBlock_.begin(), end = samplingSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      es->setAnalysisParams(*it);
    }
    es->setEmbeddedSamplingOptions(samplingOptionBlock_);
    es->setDCOptions(dcOptionBlock_);
    es->setESLinSol(esLinSolOptionBlock_);
    es->setLinSol(linSolOptionBlock_);

    return es;
  }

  //-----------------------------------------------------------------------------
  // Function      : setEmbeddedSamplingAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 5/26/2018
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Appends to any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setEmbeddedSamplingAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    for (std::vector<Util::OptionBlock>::iterator it = samplingSweepAnalysisOptionBlock_.begin(), end = samplingSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      if (Util::compareParamLists(option_block, *it))
      {
        (*it) = option_block;
        return;
      }
    }

    // save the new one.
    samplingSweepAnalysisOptionBlock_.push_back(option_block); // save a copy for later.
  }

  //-----------------------------------------------------------------------------
  // Function      : setEmbeddedSamplingOptionBlock
  // Purpose       : 
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter
  // Creation Date : 5/26/2018
  //-----------------------------------------------------------------------------
  ///
  /// Saves the sampling parsed option block.
  ///
  /// @invariant Overwrites any previously specified sampling option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setEmbeddedSamplingOptionBlock(const Util::OptionBlock &option_block)
  {
    samplingOptionBlock_ = option_block;
    return true;
  }

  bool setDCOptionBlock(const Util::OptionBlock &option_block)
  {
    dcOptionBlock_ = option_block;
    return true;
  }

  bool setESLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    esLinSolOptionBlock_ = option_block;
    return true;
  }

  bool setLinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    linSolOptionBlock_ = option_block;
    return true;
  }

public:
  AnalysisManager &                 analysisManager_;
  Linear::System &                  linearSystem_;
  Nonlinear::Manager &              nonlinearManager_;
  Device::DeviceMgr &               deviceManager_;
  Linear::Builder &                 builder_;
  Loader::Loader &                  loader_;
  Topo::Topology &                  topology_;
  IO::InitialConditionsManager &    initialConditionsManager_;

private:
  std::vector<Util::OptionBlock>    samplingSweepAnalysisOptionBlock_;
  Util::OptionBlock                 samplingOptionBlock_;
  Util::OptionBlock                 dcOptionBlock_;
  Util::OptionBlock                 esLinSolOptionBlock_;
  Util::OptionBlock                 linSolOptionBlock_;
};

//-----------------------------------------------------------------------------
// .SAMPLING
struct EmbeddedSamplingAnalysisReg : public IO::PkgOptionsReg
{
  EmbeddedSamplingAnalysisReg(
    EmbeddedSamplingFactory &             factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setEmbeddedSamplingAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  EmbeddedSamplingFactory &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractEmbeddedSamplingData
// Purpose       : Extract the parameters from a netlist .SAMPLING line held in
//                 parsed_line.
//
// Special Notes : used by the parser
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool extractEmbeddedSamplingData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("EMBEDDEDSAMPLING", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  int parameterStartPos = 1;

  // Create an option block to temporarily store the default options.
  Util::OptionBlock defaultOptions;

  // Get the default options from metadata.
  addDefaultOptionsParameters(options_manager, defaultOptions, "EMBEDDEDSAMPLING" );

  // Extract the parameters from parsed_line.
  int parameterEndPos = numFields - 1;
  Util::ParamList inputParameters;
  Util::Param parameter("", "");
  int intervalParameterStart = -1;
  int i = parameterStartPos;
  std::string paramBaseName;
  while (i <= parameterEndPos-1)
  {
    // Check for equal sign.
    if ( parsed_line[i+1].string_ != "=" )
    {
      // Stop after the tagged parameters have been extracted
      // from a .OPTIONS RESTART or .OPTIONS OUTPUT line, they
      // will be handled later.
      intervalParameterStart = i;
      break;
    }

    // Extract parameter name and value from parsed_line and add to
    // parameter list. Check to see if the parameter is "VECTOR"
    // valued and treat accordingly.
    parameter.set( parsed_line[i].string_, "" );
    Util::Param *parameterPtr = Util::findParameter(defaultOptions.begin(), defaultOptions.end(), parameter.tag());
    if (parameterPtr == NULL)
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "No options parameter " << parameter.tag() << " found, parameter will be ignored.";
      i+= 3;
    }
    else if (parameterPtr->stringValue() != "VECTOR")
    {
      parameter.setVal( parsed_line[i+2].string_ );
      inputParameters.push_back( parameter );
      i+= 3;

      if (i < parameterEndPos-1 && parsed_line[i].string_ == ",")
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << "Options parameter " << parameter.tag() << " is flagged as not VECTOR, but has comma in value.";
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      std::ostringstream paramName;
      paramBaseName = ExtendedString(parsed_line[i].string_).toUpper();
      int j = 1;

      paramName << paramBaseName << j;
      i += 2;
      parameter.set(paramName.str(), parsed_line[i].string_);
      option_block.addParam(parameter);

      // This while loop is dangerous if we are near the end of the
      // parsed line because it still assumed that the format is
      // option=value,value and not option=value,value,
      // that is no trailing comma.  It does work if
      // option=value,value is at the end of a line now (RLS 5/07)
      int testSize = parsed_line.size()-1;
      while ((i < testSize) && (parsed_line[i+1].string_ == ",") )
      {
        paramName.str("");
        ++j;
        paramName << paramBaseName << j;
        i += 2;
        parameter.set(paramName.str(), parsed_line[i].string_);
        option_block.addParam(parameter);
      }

      ++i;
    }
  }

  // For each input parameter, check that it is in the default
  // set and if so, set its value in "parameters" to the input
  // value, otherwise flag it as an unknown parameter.
  for (Util::ParamList::const_iterator it = inputParameters.begin(), end = inputParameters.end(); it != end; ++it)
  {
    Util::Param *parameterPtr = Util::findParameter(defaultOptions.begin(), defaultOptions.end(), (*it).tag());
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal(*it);
      option_block.addParam( *parameterPtr );
    }
    else
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "No options parameter " << (*it).tag() << " found, parameter will be ignored.";
    }
  }

  circuit_block.addOptions(option_block);

  return true; 
}

//-----------------------------------------------------------------------------
// Function      : populateMetadata
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
void populateMetadata(IO::PkgOptionsMgr & options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("EMBEDDEDSAMPLING");

    parameters.insert(Util::ParamMap::value_type("PARAM", Util::Param("PARAM", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("TYPE", Util::Param("TYPE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MEANS", Util::Param("MEANS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("STD_DEVIATIONS", Util::Param("STD_DEVIATIONS", "VECTOR")));

    parameters.insert(Util::ParamMap::value_type("LOWER_BOUNDS", Util::Param("LOWER_BOUNDS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("UPPER_BOUNDS", Util::Param("UPPER_BOUNDS", "VECTOR")));

    parameters.insert(Util::ParamMap::value_type("ALPHA", Util::Param("ALPHA", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("BETA", Util::Param("BETA", "VECTOR")));

    //parameters.insert(Util::ParamMap::value_type("COVMATRIX", Util::Param("COVMATRIX", "VECTOR")));
  }

  {
    // Must use "EMBEDDEDSAMPLES" instead of "SAMPLING!"
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("EMBEDDEDSAMPLES");

    parameters.insert(Util::ParamMap::value_type("NUMSAMPLES", Util::Param("NUMSAMPLES", 1)));
    parameters.insert(Util::ParamMap::value_type("COVMATRIX", Util::Param("COVMATRIX", "VECTOR")));

    parameters.insert(Util::ParamMap::value_type("OUTPUTFORMAT", Util::Param("OUTPUTFORMAT", "STD")));
    parameters.insert(Util::ParamMap::value_type("OUTPUTS", Util::Param("OUTPUTS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("OUTPUTALLSAMPLES", Util::Param("OUTPUTALLSAMPLES", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTSAMPLESTATS", Util::Param("OUTPUTSAMPLESTATS", true)));
#if Xyce_STOKHOS_ENABLE
    parameters.insert(Util::ParamMap::value_type("REGRESSION_PCE", Util::Param("REGRESSION_PCE", false)));
    parameters.insert(Util::ParamMap::value_type("PROJECTION_PCE", Util::Param("PROJECTION_PCE", false)));
    parameters.insert(Util::ParamMap::value_type("ORDER", Util::Param("ORDER", 4)));
#endif
    parameters.insert(Util::ParamMap::value_type("SAMPLE_TYPE", Util::Param("SAMPLE_TYPE", 0)));
    parameters.insert(Util::ParamMap::value_type("SEED", Util::Param("SEED", 0)));

    parameters.insert(Util::ParamMap::value_type("RESAMPLE", Util::Param("RESAMPLE", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_PCE_COEFFS", Util::Param("OUTPUT_PCE_COEFFS", false)));
    parameters.insert(Util::ParamMap::value_type("SPARSE_GRID", Util::Param("SPARSE_GRID", false)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("STDOUTPUT", Util::Param("STDOUTPUT", false)));
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
bool registerEmbeddedSamplingFactory(FactoryBlock & factory_block)
{
  EmbeddedSamplingFactory *factory = new EmbeddedSamplingFactory(
      factory_block.analysisManager_, 
      factory_block.linearSystem_, 
      factory_block.nonlinearManager_, 
      factory_block.deviceManager_, 
      factory_block.builder_, 
      factory_block.loader_, 
      factory_block.topology_, 
      factory_block.initialConditionsManager_
      );

  addAnalysisFactory(factory_block, factory);

  populateMetadata(factory_block.optionsManager_);

  factory_block.optionsManager_.addCommandParser(".EMBEDDEDSAMPLING", extractEmbeddedSamplingData);
  factory_block.optionsManager_.addCommandProcessor("EMBEDDEDSAMPLING", new EmbeddedSamplingAnalysisReg(*factory));

  // WARNING:  If you use "EMBEDDEDSAMPLING" here for .options, it will 
  // confuse the "processor" code, 
  // mixing and matching .SAMPLING with .options SAMPLING
  // So, using "EMBEDDEDSAMPLES" instead as identifier.
  factory_block.optionsManager_.addOptionsProcessor("EMBEDDEDSAMPLES", IO::createRegistrationOptions<EmbeddedSamplingFactory>(*factory, &EmbeddedSamplingFactory::setEmbeddedSamplingOptionBlock)); 

  // DC options might not be needed
  factory_block.optionsManager_.addOptionsProcessor("DC", IO::createRegistrationOptions <EmbeddedSamplingFactory> (*factory, &EmbeddedSamplingFactory::setDCOptionBlock)); 

  factory_block.optionsManager_.addOptionsProcessor("LINSOL-ES", IO::createRegistrationOptions <EmbeddedSamplingFactory> (*factory, &EmbeddedSamplingFactory::setESLinSolOptionBlock));
  factory_block.optionsManager_.addOptionsProcessor("LINSOL", IO::createRegistrationOptions <EmbeddedSamplingFactory> (*factory, &EmbeddedSamplingFactory::setLinSolOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce

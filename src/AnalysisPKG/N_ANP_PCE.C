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
// Purpose       : PCE class analysis functions.
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#if Xyce_STOKHOS_ENABLE

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ANP_PCE.h>
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
#include <N_PDS_EpetraHelpers.h>

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

// new PCE loader class
#include <N_LOA_PCELoader.h>

#include <N_LAS_Graph.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_System.h>
// new PCE headers: (patterned, kind of, off the HB headers)  PCE = Polynomial Chaos Expansion
#include <N_LAS_PCEBuilder.h>
#include <N_LAS_PCESolverFactory.h>
#include <N_LAS_TranSolverFactory.h>

#include <N_PDS_fwd.h>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_Utils.hpp>
#include <Teuchos_LAPACK.hpp>


#if Xyce_STOKHOS_ENABLE
#include <Sacado_No_Kokkos.hpp>
#include <Stokhos_Sacado.hpp>
#endif



namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : PCE::PCE
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
PCE::PCE(
    AnalysisManager &analysis_manager, 
    Linear::System & linear_system,
    Nonlinear::Manager & nonlinear_manager,
    Device::DeviceMgr & device_manager,
    Linear::Builder & builder,
    Loader::Loader &loader, 
    Topo::Topology & topology, 
    IO::InitialConditionsManager & initial_conditions_manager,
    AnalysisBase &child_analysis)
    : AnalysisBase(analysis_manager, "PCE"),
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

      pceLoaderPtr_(0),
      pceBuilderPtr_(0),
      pceLinearSystem_(0),
      solverFactory_(0),

      samplingVector_(),
      maxParamStringSize_(0),
      lower_bounds_Given_(false),
      upper_bounds_Given_(false),
      covMatrixGiven_(false),
      numBlockRows_(1),
      numQuadPoints_(1),
      sampleType_(UQ::LHS), 
      userSeed_(0),
      userSeedGiven_(false),
      hackOutputFormat_("TECPLOT"),
      hackOutputCalledBefore_(false),
      hackOutputCalledBefore2_(false),
      hackOutputAllSamples_(false),
      outputSampleStats_(false),
#if Xyce_STOKHOS_ENABLE
      PCEorder_(4),
      resamplePCE_(false),
      outputPCECoeffs_(false),
      numResamples_(1000),
      numResamplesGiven_(false),
      useSparseGrid_(false),
#endif
      stdOutputFlag_(false),
      debugLevel_(0),
      outputStochasticMatrix_(false),
      voltLimAlgorithm_(2),
      stdOutputFlag_xvec_(false),
      outputsGiven_(false),
      outputsSetup_(false),
      measuresGiven_(false),
      outFuncGIDsetup_(false),
      outputIndex_(0),
      outputIndex2_(0),
      resetForStepCalledBefore_(false),
      useExpressionSamples_(false)
{
  pdsMgrPtr_ = analysisManager_.getPDSManager();

  nextSolutionPtr_ = builder_.createVector();
  nextStatePtr_ = builder_.createStateVector();
  nextStorePtr_ = builder_.createStoreVector();
}

//-----------------------------------------------------------------------------
// Function      : PCE::~PCE
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
PCE::~PCE()
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

  if (pceLoaderPtr_)
  {
    delete pceLoaderPtr_;
  }

  if (pceLinearSystem_)
  {
    delete pceLinearSystem_;
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
// Function      : PCE::notify
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void PCE::notify(const StepEvent &event) 
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
      pceLoaderPtr_->loadDeviceErrorWeightMask(dsPtr->deviceErrorWeightMask_);

      analysisManager_.getXyceTranTimer().resetStartTime();
    }

    resetForStepCalledBefore_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : PCE::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::setAnalysisParams(const Util::OptionBlock & paramsBlock)
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
      if (stdDev < 0) { Report::DevelFatal() << "STD_DEVIATIONS values for .PCE must be >= 0";}
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
      if (alpha <= 0) { Report::DevelFatal() << "ALPHA values for .PCE must be > 0";}
      alphaVec_.push_back(alpha);
    }
    else if (std::string( iter->uTag() ,0,4) == "BETA")
    {
      betaGiven=true;
      double beta = iter->getImmutableValue<double>();
      if (beta <= 0) { Report::DevelFatal() << "BETA values for .PCE must be > 0";}
      betaVec_.push_back(beta);
    }
    else if (iter->uTag() == "USEEXPR")
    {
      useExpressionSamples_ = static_cast<bool>(iter->getImmutableValue<bool>());
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() 
        << " is not a recognized sampling option.\n" << std::endl;
    }
  }

  if (useExpressionSamples_)
  {
    N_PDS_Manager &pds_manager = *analysisManager_.getPDSManager();
    N_PDS_Comm & pdsComm = *(pds_manager.getPDSComm());

    SweepVector exprSamplingVector_;
    loader_.getRandomParams(exprSamplingVector_, pdsComm);
    samplingVector_.insert
      (samplingVector_.end(), exprSamplingVector_.begin(), exprSamplingVector_.end());

    paramNameVec_.resize(samplingVector_.size());
    for (int ii=0;ii<samplingVector_.size();ii++) { paramNameVec_[ii] = samplingVector_[ii].name; }
  }
  else
  {
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
      else
      {
        Report::DevelFatal().in("parsePCEParam") << "Unsupported SAMPLING type";
      }
      samplingVector_.push_back(sampling_param);
    }
  }

  outputManagerAdapter_.setStepSweepVector(samplingVector_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCE::setPCEOptions
// Purpose       :
// Special Notes : These are from '.options PCES'
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::setPCEOptions(const Util::OptionBlock & option_block)
{
  bool bsuccess = true;

  Util::ParamList::const_iterator it  = option_block.begin();
  Util::ParamList::const_iterator end = option_block.end();
  for ( ; it != end; ++ it)
  {
    if (std::string((*it).uTag() ,0,9) == "COVMATRIX" ) // this is a vector
    {
      covMatrixGiven_ = true;
      covMatrix_.push_back( (*it).getImmutableValue<double>() );
    }
    else if ((*it).uTag() == "OUTPUTFORMAT" ) 
    {
      ExtendedString tag = (*it).stringValue();
      hackOutputFormat_ = tag.toUpper();
    }
#if 0
    else if ((*it).uTag() == "OUTPUT_ALL_SAMPLES")
    {
      hackOutputAllSamples_=static_cast<bool>((*it).getImmutableValue<bool>());
    }
#endif
    else if ((*it).uTag() == "OUTPUT_SAMPLE_STATS")
    {
      outputSampleStats_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
#if Xyce_STOKHOS_ENABLE
    else if ((*it).uTag() == "ORDER")
    {
      PCEorder_ = (*it).getImmutableValue<int>();
      if (PCEorder_ < 0)
        Report::UserError() << "ORDER parameter on .OPTIONS PCES line must >= 0";
    }
#endif
    else if ((*it).uTag() == "SAMPLE_TYPE")
    {
      ExtendedString p((*it).stringValue()); p.toUpper();
      if (p == "MC")
      {
        sampleType_ = UQ::MC;
      }
      else if (p == "LHS")
      {
        sampleType_ = UQ::LHS;
      }
      else
      {
        Xyce::Report::UserWarning() << (*it).uTag() 
          << " = " << p << " is not a recognized sampling option.  Setting " << (*it).uTag() << " = MC.\n" << std::endl;
        sampleType_ = UQ::MC;
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
      if (DEBUG_PCE)
      {
        debugLevel_ = (*it).getImmutableValue<int>();
      }
    }
    else if ((*it).uTag() == "OUTPUT_STOCHASTIC_MATRIX")
    {
      outputStochasticMatrix_ = (*it).getImmutableValue<bool>();
    }
    else if ((*it).uTag() == "VOLTLIM_ALG")
    {
      voltLimAlgorithm_ = (*it).getImmutableValue<int>();
    }
    else if ((*it).uTag() == "STDOUTPUT_XVEC")
    {
      stdOutputFlag_xvec_ = (*it).getImmutableValue<bool>();
    }

    else if (std::string((*it).uTag() ,0,7) == "OUTPUTS" )// this is a vector of expressions/solution variables
    {
      outputsGiven_ = true;
      UQ::outputFunctionData * ofDataPtr = new UQ::outputFunctionData();
      ExtendedString funcName = (*it).stringValue();
      ofDataPtr->outFuncString = funcName.toUpper();
      outFuncDataVec_.push_back(ofDataPtr);
    }
    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() 
        << " is not a recognized sampling option.\n" << std::endl;
    }
  }

  // parse the expression now, so if there are any errors, they will come
  // up early in the simulation.
  if (outputsGiven_)
  {
    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      outFuncDataVec_[iout]->expDataPtr = new Util::ExpressionData(
          analysisManager_.getExpressionGroup(),
          outFuncDataVec_[iout]->outFuncString);
    }
  }
  else if (measuresGiven_)
  {
    for (int iout=0;iout<measFuncDataVec_.size();++iout)
    {
      measFuncDataVec_[iout]->measureResponseFound = measureManager_.find(measFuncDataVec_[iout]->outFuncString);

      if (!(measFuncDataVec_[iout]->measureResponseFound))
      {
        Report::UserWarning0() << "Measure response " << measFuncDataVec_[iout]->outFuncString << " was not found.";
      }
    }
  }
  else // give a warning.
  {
    Report::UserWarning0() << "Output function was not specified";
  }

  // signal the OutputMgr that PCE is enabled
  outputManagerAdapter_.setEnablePCEFlag(true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCE::getTIAParams() const
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
const TimeIntg::TIAParams & PCE::getTIAParams() const
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : PCE::getTIAParams()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
TimeIntg::TIAParams & PCE::getTIAParams()
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : PCE::getDCOPFlag()
// Purpose       :
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::getDCOPFlag() const
{
  return childAnalysis_.getDCOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : PCE::stepCallBack
// Purpose       : 
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void PCE::stepCallBack ()
{
  outputXvectors(); // temporary
  computePCEOutputs();

#if Xyce_STOKHOS_ENABLE
  std::vector< Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > > pceVec;

  for (int iout=0;iout<outFuncDataVec_.size();++iout)
  {
    UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);

    Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = outFunc.projectionPCE;
    projectionPCE.init(0.0);
    projectionPCE.reset(expnMethod);

    std::vector<double> & f = outFunc.sampleOutputs;

    UQ::solveProjectionPCE(basis, quadMethod, f, projectionPCE);

    if (outputPCECoeffs_)
    {
      Xyce::lout() << "PCE coefs for " << outFunc.outFuncString << " = " << projectionPCE << std::endl;
      projectionPCE.print(Xyce::lout());
    }

    if (stdOutputFlag_)
    {
      Xyce::lout() << std::endl;
      Xyce::lout() << "Intrusive PCE mean of " << outFunc.outFuncString << " = " << projectionPCE.mean() << std::endl;
      Xyce::lout() << "Intrusive PCE stddev of " << outFunc.outFuncString << " = " << projectionPCE.standard_deviation() << std::endl;
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
      Xyce::lout() << "Statistics from re-sampling PCE approximation of " << outFunc.outFuncString << ":" << std::endl;;
      Xyce::lout() << statVec[iout];
      Xyce::lout() << std::endl;
    }
  }
#endif

  // only call outputter functions if OUTPUTS= was used on the
  // .OPTIONS PCES line.
  if (outputsGiven_)
  {
    // output for .PRINT PCE
#if Xyce_STOKHOS_ENABLE
    outputManagerAdapter_.outputPCE(numQuadPoints_, outFuncDataVec_);
#endif

    hackPCEOutput ();
  }
  hackPCEOutput2 ();
}

//-----------------------------------------------------------------------------
// Function      : PCE::doRun()
// Purpose       : This is the main controlling loop for sampling analysis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : PCE::doInit()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::doInit()
{
  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "PCE::init" << std::endl;
  }

  UQ::checkParameterList(
      analysisManager_.getComm(), 
      loader_, 
      samplingVector_.begin(), 
      samplingVector_.end());

  setupStokhosObjects ();
  setupBlockSystemObjects ();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCE::setupBlockSystemObjects 
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
void  PCE::setupBlockSystemObjects ()
{
  analysisManager_.resetSolverSystem();

  pceBuilderPtr_ = rcp(new Linear::PCEBuilder(numBlockRows_,numQuadPoints_));

  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << "PCE::setupBlockSystemObjects():  Generate Maps,etc\n";
  }
  {
    Stats::StatTop _setupStepStat("Setup Maps/Graphs");
    Stats::TimeBlock _setupStepTimer(_setupStepStat);

    pceBuilderPtr_->generateMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::SOLUTION ), false),
                                  rcp(pdsMgrPtr_->getParallelMap( Parallel::SOLUTION_OVERLAP_GND ), false) );

    pceBuilderPtr_->generateStateMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::STATE ),false) );
    pceBuilderPtr_->generateStoreMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::STORE ),false) );
    pceBuilderPtr_->generateLeadCurrentMaps( rcp(pdsMgrPtr_->getParallelMap( Parallel::LEADCURRENT ),false) );

    pceBuilderPtr_->generateGraphs( *pceGraph, *pdsMgrPtr_->getMatrixGraph( Parallel::JACOBIAN ));
  }

  // Create PCE Loader.
  delete pceLoaderPtr_;
  pceLoaderPtr_ = new Loader::PCELoader(deviceManager_, builder_, numQuadPoints_, numBlockRows_, samplingVector_, Y_, analysisManager_.getCommandLine(), voltLimAlgorithm_, useExpressionSamples_);
  pceLoaderPtr_->registerPCEBuilder(pceBuilderPtr_);
  pceLoaderPtr_->registerAppLoader( rcp(&loader_, false) );
  pceLoaderPtr_->registerPCEbasis (basis) ;
  pceLoaderPtr_->registerPCEquadMethod (quadMethod) ;
  pceLoaderPtr_->registerPCEexpnMethod (expnMethod) ;
  pceLoaderPtr_->registerPCEtripleProductTensor (Cijk) ;

  //-----------------------------------------

  //Finish setup of PCE Loader
  //-----------------------------------------
  {
    //-----------------------------------------
    //Construct Solvers, etc.
    //-----------------------------------------
    delete pceLinearSystem_;
    pceLinearSystem_ = new Linear::System();
    //-----------------------------------------

    //hack needed by TIA initialization currently
    pceBuilderPtr_->registerPDSManager( pdsMgrPtr_ );
    pceLinearSystem_->registerPDSManager( pdsMgrPtr_ );
    pceLinearSystem_->registerBuilder( &*pceBuilderPtr_ );

    // the linear system class will only create a "linear problem" object if the 
    // builder object (in this case PCE builder) will create a non-NULL matrix pointer.
    pceLinearSystem_->initializeSystem();

    nonlinearManager_.setLinSolOptions( saved_lsPCEOB_ );
    nonlinearManager_.setMatrixFreeFlag( false ); // NOT TRUE for ES!

    if (!solverFactory_)
    {
      solverFactory_ = new Linear::PCESolverFactory( *pceBuilderPtr_ );
      solverFactory_->registerPCELoader( rcp(pceLoaderPtr_, false) );
      solverFactory_->registerPCEBuilder( pceBuilderPtr_ );
      solverFactory_->setNumCoefs( numBlockRows_ );
      solverFactory_->setNumQuadPoints( numQuadPoints_ );
      solverFactory_->setCoefsOuterLoop (coefsOuterLoop_);
    }

    nonlinearManager_.registerSolverFactory( solverFactory_ );

  //Initialization of Solvers, etc. 
    analysisManager_.initializeSolverSystem(
        getTIAParams(),
        *pceLoaderPtr_,
        *pceLinearSystem_,
        nonlinearManager_,
        deviceManager_);

    nonlinearManager_.initializeAll(
        analysisManager_, 
        analysisManager_.getNonlinearEquationLoader(), 
        *pceLinearSystem_,
        *analysisManager_.getDataStore(),
        *analysisManager_.getPDSManager(),
        initialConditionsManager_,
        analysisManager_.getOutputManagerAdapter().getOutputManager(),
        topology_);

    childAnalysis_.registerLinearSystem(pceLinearSystem_);
 
    // don't have a mode yet.  Do I need one??  The real mode should be 
    // the underlying analysis (tran, dc, etc), so probably not
    //nonlinearManager_.setAnalysisMode(nonlinearAnalysisMode(ANP_MODE_PCE)); 
  }
  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  pceLoaderPtr_->loadDeviceErrorWeightMask(dsPtr->deviceErrorWeightMask_);
  pceLoaderPtr_->registerSolverFactory ( solverFactory_ );
  //pceLoaderPtr_->registerLinearSystem( &linearSystem_ );
  pceLoaderPtr_->setLinSolOptions( saved_lsOB_ );

  childAnalysis_.registerParentAnalysis(this);
}

//-----------------------------------------------------------------------------
// Function      : PCE::setupStokhosObjects 
// Purpose       : This function sets up the stokhos ensemble objects, using
//                 the computed sample values.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void PCE::setupStokhosObjects ()
{
#if Xyce_STOKHOS_ENABLE
  // determine the samples. 
  //
  // spectral projection samples are determined 
  // by the quadrature points
  //
  const int d = paramNameVec_.size();
  const int p = PCEorder_;
  bases.resize(d); 
  for (int i=0; i<d; i++)
  {
    SweepParam & sp = samplingVector_[i];

    if (sp.type == "UNIFORM")
    {
      bases[i] = rcp(new Stokhos::LegendreBasis<int,double>(p));
    }
    else if (sp.type == "NORMAL") 
    {
      bases[i] = rcp(new Stokhos::HermiteBasis<int,double>(p,true));
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

    basis = rcp(new Stokhos::SmolyakBasis<int,double,total_less>( bases, index_set, drop));
    quadMethod = rcp(new Stokhos::SmolyakSparseGridQuadrature<int,double>(basis, index_set));
  }
  else
  {
    basis = rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    quadMethod = rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
  }

  Cijk = basis->computeTripleProductTensor();
  expnMethod = rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, Cijk, quadMethod));

  UQ::setupPCEQuadPoints ( basis, quadMethod, expnMethod, samplingVector_, covMatrix_, meanVec_, X_, Y_);
  numQuadPoints_ = quadMethod->size();

  Epetra_Comm* petraComm = Parallel::getEpetraComm( analysisManager_.getPDSManager()->getPDSComm() );

  if (outputStochasticMatrix_)
  {
    // this is to give a nice SPY plot of the stochastic block matrix structure
    std::string file="A.mm";
    Stokhos::sparse3Tensor2MatrixMarket(
        *basis, 
        *Cijk, 
        *petraComm, file);
  }

  pceGraph = rcp( new Linear::Graph( Stokhos::sparse3Tensor2CrsGraph(*basis, *Cijk, *petraComm ) ) );

  numBlockRows_ = pceGraph->numLocalEntities();

#if 0
  std::cout << "Cijk:" <<std::endl;
  Cijk->print(std::cout);

  // this is the matrix that gets output to A.mm
  Epetra_CrsMatrix mat(Copy, *pceGraph->epetraObj());
  pceGraph->epetraObj()->PrintGraphData(std::cout);
  mat.Print(std::cout);
#endif
#endif
}

//-----------------------------------------------------------------------------
// Function      : PCE::doLoopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::doLoopProcess()
{
  // This type of analysis solves a fully intrusive PCE system.  It builds upon
  // the child analysis type, but instead of solving directly for circuit unknowns, 
  // it solves for coefficients of the PCE expansion

  // This is being set up by considering the original analysis (DC or TRAN or whatever) 
  // to be the "child" analysis, but to then replace the linear objects with block 
  // versions and replace the loader with a block loader.  Otherwise, all the control 
  // is handled in the child process.
  Xyce::lout() << "***** Beginning Intrusive PCE simulation....\n" << std::endl;
  Xyce::lout() << "***** Number of quadrature points = " << numQuadPoints_ << "\n" << std::endl;
  Xyce::lout() << "***** Number of linear system block rows = " << numBlockRows_ << "\n" << std::endl;

  // test:
  analysisManager_.setAnalysisMode(Xyce::Analysis::ANP_MODE_TRANSIENT);

  // solve the loop.
  bool status = childAnalysis_.run();

  return status;
}

//-----------------------------------------------------------------------------
// Function      : PCE::doProcessSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::doProcessSuccessfulStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCE::doProcessFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::doProcessFailedStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PCE::doFinish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::doFinish()
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : PCE::convertPointToPCE
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/28/2019
//-----------------------------------------------------------------------------
void PCE::convertPointToPCE (int index, Stokhos::OrthogPolyApprox<int,double> & pceApprox)
{
  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextSolutionPtr) );

  pceApprox.reset(basis);
  int basisSize = basis->size();
  for (int icoef=0;icoef<basisSize;++icoef)
  {
    *nextSolutionPtr_ = bX.block(icoef);
    pceApprox[icoef] = (*nextSolutionPtr_)[index];
  }
}

//-----------------------------------------------------------------------------
// Function      : PCE::evaluateVector
// Purpose       : 
// Special Notes : maybe rewrite this as stand-alone function
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/28/2019
//-----------------------------------------------------------------------------
void PCE::evaluateVector ( Teuchos::RCP<Linear::BlockVector> & bX_quad_ptr_ )
{
  std::vector< Stokhos::OrthogPolyApprox<int,double> > pceVec(1);

  int solutionSize = bX_quad_ptr_->block(0).localLength(); // serial only

  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextSolutionPtr) );

  for (int isol=0;isol<solutionSize;++isol)
  {
    convertPointToPCE(isol, pceVec[0]);

    std::vector < std::vector<double> > xSamples(pceVec.size(), std::vector<double>(numQuadPoints_,0.0) );
    Xyce::Analysis::UQ::evaluateApproximationPCE(samplingVector_, Y_, numQuadPoints_, pceVec, xSamples);

    for (int iquad=0;iquad<numQuadPoints_;++iquad)
    {
      (bX_quad_ptr_->block(iquad))[isol] = xSamples[0][iquad];
    }
  }
  bX_quad_ptr_->assembleGlobalVector();
}

//-----------------------------------------------------------------------------
// Function      : PCE::outputXvectors
// Purpose       : Outputs both the coefficient vector and also a converted quad vector to std output
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/28/2019
//-----------------------------------------------------------------------------
void PCE::outputXvectors()
{
  if (stdOutputFlag_xvec_)
  {
    TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
    Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextSolutionPtr) );

    Xyce::lout() << "--------------------------------------------------------------" <<std::endl;
    Xyce::lout() << "X coef vector:" <<std::endl;
    bX.print(Xyce::lout());
    Xyce::lout() << "--------------------------------------------------------------" <<std::endl;

    Teuchos::RCP<Linear::BlockVector> bXNext_quad_ptr_ = Teuchos::rcp( pceBuilderPtr_->createQuadVector() ); 
    evaluateVector(bXNext_quad_ptr_);

    Xyce::lout() << "--------------------------------------------------------------" <<std::endl;
    Xyce::lout() << "X quad vector:" <<std::endl;
    bXNext_quad_ptr_->print(Xyce::lout());
    Xyce::lout() << "--------------------------------------------------------------" <<std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : PCE::hackPCEOutput2
// Purpose       : Outputs both the coefficient vector and also a converted quad vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/28/2019
//-----------------------------------------------------------------------------
void PCE::hackPCEOutput2()
{
  TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
  Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextSolutionPtr) );

  int solutionSize = bX.block(0).localLength();

  std::string fileName; 

  if (hackOutputFormat_=="TECPLOT")
  {
    fileName = analysisManager_.getNetlistFilename() + "_pce2.dat";
  }
  else if (hackOutputFormat_=="STD")
  {
    fileName = analysisManager_.getNetlistFilename() + "_pce2.prn";
  }
  else
  {
    Xyce::Report::UserWarning() << hackOutputFormat_ 
      << " is not a recognized intrusive PCE output option.\n" << std::endl;
  }

  std::ofstream output_stream(fileName.c_str(), 
          !hackOutputCalledBefore2_ ? std::ios_base::out : std::ios_base::app);

  if (!hackOutputCalledBefore2_)
  {
    output_stream << "TITLE = \"intrusive PCE output, solution vector coefficients \"\tVARIABLES= "<<std::endl;

    if ( !(childAnalysis_.getDCOPFlag()) )
    {
      output_stream << "\t\" TIME \""<<std::endl;
    }
    else
    {
      output_stream << "\t\" INDEX \""<<std::endl;
    }

    for (int isol=0;isol<solutionSize;++isol)
    {
      std::string varName = "var_" + std::to_string(isol);

      std::string meanString = varName + "_quad_pce_mean";
      std::string meanStringPlus = varName + "_quad_pce_meanPlus";
      std::string meanStringMinus = varName + "_quad_pce_meanMinus";
      std::string meanStringPlusTwoSigma = varName + "_quad_pce_meanPlusTwoSigma";
      std::string meanStringMinusTwoSigma = varName + "_quad_pce_meanMinusTwoSigma";

      std::string stddevString = varName + "_quad_pce_stddev";
      std::string varianceString = varName + "_quad_pce_variance";

      output_stream << "\t\" " << meanString << "\""<<std::endl;
      output_stream << "\t\" " << meanStringPlus << "\""<<std::endl;
      output_stream << "\t\" " << meanStringMinus << "\""<<std::endl;
      output_stream << "\t\" " << meanStringPlusTwoSigma << "\""<<std::endl;
      output_stream << "\t\" " << meanStringMinusTwoSigma << "\""<<std::endl;

      output_stream << "\t\" " << stddevString << "\""<<std::endl;
      output_stream << "\t\" " << varianceString << "\""<<std::endl;

      int basisSize = basis->size();
      for (int icoef=0;icoef<basisSize;++icoef)
      {
          std::string sampleString = varName + "_coef_"+std::to_string(icoef);
          output_stream << "\t\" " << sampleString << "\""<<std::endl;
      }
    }
 
    output_stream << "ZONE F=POINT  T=\"Xyce data\""<<std::endl;
    hackOutputCalledBefore2_ = true;
  }

  // data output
  output_stream.setf(std::ios::scientific);

  if ( !(childAnalysis_.getDCOPFlag()) )
  {
    output_stream << analysisManager_.getStepErrorControl().currentTime;
  }
  else
  {
    output_stream << outputIndex2_;
    outputIndex2_++;
  }

  for (int isol=0;isol<solutionSize;++isol)
  {
    Stokhos::OrthogPolyApprox<int,double> pceApprox;
    convertPointToPCE(isol,pceApprox);

    double pce_mean = pceApprox.mean();
    double pce_stddev = pceApprox.standard_deviation();
    double pce_variance = pce_stddev*pce_stddev;

    if ( std::isinf(pce_mean) || std::isnan(pce_mean) )
    {
      pce_mean = 0.0;
    }

    if ( std::isinf(pce_stddev) || std::isnan(pce_stddev) )
    {
      pce_stddev = 0.0;
    }

    if ( std::isinf(pce_variance) || std::isnan(pce_variance) )
    {
      pce_variance = 0.0;
    }

    output_stream << "\t" << pce_mean;

    output_stream << "\t" << (pce_mean+pce_stddev);
    output_stream << "\t" << (pce_mean-pce_stddev);
    output_stream << "\t" << (pce_mean+2*pce_stddev);
    output_stream << "\t" << (pce_mean-2*pce_stddev);

    output_stream << "\t" << pce_stddev;
    output_stream << "\t" << pce_variance;

    int basisSize = basis->size();
    for (int icoef=0;icoef<basisSize;++icoef)
    {
      output_stream << "\t" << pceApprox[icoef];
    }
  }

  output_stream << std::endl;

  return;
}


//-----------------------------------------------------------------------------
// Function      : PCE::computePCEOutputs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
void PCE::computePCEOutputs()
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
    //Linear::BlockVector & bX = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextSolutionPtr) );
    Linear::BlockVector & bSta = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextStatePtr) );
    Linear::BlockVector & bSto = *dynamic_cast<Xyce::Linear::BlockVector*>( (dsPtr->nextStorePtr) );

    Teuchos::RCP<Linear::BlockVector> bXNext_quad_ptr_ = Teuchos::rcp( pceBuilderPtr_->createQuadVector() ); 
    evaluateVector(bXNext_quad_ptr_);

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

    int BlockCount = bXNext_quad_ptr_->blockCount();
    for( int i = 0; i < BlockCount; ++i )
    {
      *nextSolutionPtr_ = bXNext_quad_ptr_->block(i);
      *nextStatePtr_ = bSta.block(i);
      *nextStorePtr_ = bSto.block(i);

      //get expression value and compute means, etc
      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
        Util::Op::OpData opDataTmp(0, nextSolutionPtr_, 0, nextStatePtr_, nextStorePtr_, 0);
        double val;
        outFunc.expDataPtr->evaluate(comm, circuit_time, circuit_dt, opDataTmp, val);
        outFunc.sampleOutputs.push_back(val);
      }
    }

    // finish mean and std dev
    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
      outFunc.completeStatistics();

      if (stdOutputFlag_ && outputSampleStats_)
      {
        // histrogram is a hack that doesn't work yet
        //UQ::histrogram(std::cout, outFunc.outFuncString, outFunc.sampleOutputs);

        std::string sampleTypeStr = "PCE Quad Points";
        outFunc.output(Xyce::lout(), sampleTypeStr);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PCE::hackPCEOutput ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
void PCE::hackPCEOutput ()
{
    std::string fileName; 

    if (hackOutputFormat_=="TECPLOT")
    {
      fileName = analysisManager_.getNetlistFilename() + "_pce.dat";
    }
    else if (hackOutputFormat_=="STD")
    {
      fileName = analysisManager_.getNetlistFilename() + "_pce.prn";
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
      output_stream << "TITLE = \"intrusive PCE output\"\tVARIABLES= "<<std::endl;

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
          std::string meanStringPlusTwoSigma = outFunc.outFuncString + "_meanPlusTwoSigma";
          std::string meanStringMinusTwoSigma = outFunc.outFuncString + "_meanMinusTwoSigma";

          std::string stddevString = outFunc.outFuncString + "_stddev";
          std::string varianceString = outFunc.outFuncString + "_variance";

          output_stream << "\t\" " << meanString << "\""<<std::endl;
          output_stream << "\t\" " << meanStringPlus << "\""<<std::endl;
          output_stream << "\t\" " << meanStringMinus << "\""<<std::endl;
          output_stream << "\t\" " << meanStringPlusTwoSigma << "\""<<std::endl;
          output_stream << "\t\" " << meanStringMinusTwoSigma << "\""<<std::endl;

          output_stream << "\t\" " << stddevString << "\""<<std::endl;
          output_stream << "\t\" " << varianceString << "\""<<std::endl;
        }

#if Xyce_STOKHOS_ENABLE
        std::string meanString = outFunc.outFuncString + "_quad_pce_mean";
        std::string meanStringPlus = outFunc.outFuncString + "_quad_pce_meanPlus";
        std::string meanStringMinus = outFunc.outFuncString + "_quad_pce_meanMinus";
        std::string meanStringPlusTwoSigma = outFunc.outFuncString + "_quad_pce_meanPlusTwoSigma";
        std::string meanStringMinusTwoSigma = outFunc.outFuncString + "_quad_pce_meanMinusTwoSigma";

        std::string stddevString = outFunc.outFuncString + "_quad_pce_stddev";
        std::string varianceString = outFunc.outFuncString + "_quad_pce_variance";

        output_stream << "\t\" " << meanString << "\""<<std::endl;
        output_stream << "\t\" " << meanStringPlus << "\""<<std::endl;
        output_stream << "\t\" " << meanStringMinus << "\""<<std::endl;
        output_stream << "\t\" " << meanStringPlusTwoSigma << "\""<<std::endl;
        output_stream << "\t\" " << meanStringMinusTwoSigma << "\""<<std::endl;

        output_stream << "\t\" " << stddevString << "\""<<std::endl;
        output_stream << "\t\" " << varianceString << "\""<<std::endl;
#endif
        if (hackOutputAllSamples_)
        {
          for(int i=0;i<numQuadPoints_; ++i)
          {
            std::string sampleString = outFunc.outFuncString + "_"+std::to_string(i);
            output_stream << "\t\" " << sampleString << "\""<<std::endl;
          }
        }
      }
   
      output_stream << "ZONE F=POINT  T=\"Xyce data\""<<std::endl;
      hackOutputCalledBefore_ = true;
    }

    // data output
    output_stream.setf(std::ios::scientific);

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
        output_stream << "\t" << (outFunc.sm.mean+2*outFunc.sm.stddev);
        output_stream << "\t" << (outFunc.sm.mean-2*outFunc.sm.stddev);

        output_stream << "\t" << outFunc.sm.stddev;
        output_stream << "\t" << outFunc.sm.variance;
      }

#if Xyce_STOKHOS_ENABLE
      Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = outFunc.projectionPCE;

      double pce_mean = projectionPCE.mean();
      double pce_stddev = projectionPCE.standard_deviation();
      double pce_variance = pce_stddev*pce_stddev;

      if ( std::isinf(pce_mean) || std::isnan(pce_mean) )
      {
        pce_mean = 0.0;
      }

      if ( std::isinf(pce_stddev) || std::isnan(pce_stddev) )
      {
        pce_stddev = 0.0;
      }

      if ( std::isinf(pce_variance) || std::isnan(pce_variance) )
      {
        pce_variance = 0.0;
      }

      output_stream << "\t" << pce_mean;

      output_stream << "\t" << (pce_mean+pce_stddev);
      output_stream << "\t" << (pce_mean-pce_stddev);
      output_stream << "\t" << (pce_mean+2*pce_stddev);
      output_stream << "\t" << (pce_mean-2*pce_stddev);

      output_stream << "\t" << pce_stddev;
      output_stream << "\t" << pce_variance;
#endif

      // output individual samples
      if (hackOutputAllSamples_)
      {
        for(int i=0;i<numQuadPoints_; ++i)
        {
          output_stream << "\t" << outFunc.sampleOutputs[i];
        }
      }
    }

    output_stream << std::endl;

    return;
}

//-----------------------------------------------------------------------------
// Function      : PCE::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool PCE::setDCOptions(const Util::OptionBlock & paramsBlock)
{
  dcSweepVector_.push_back(parseSweepParams(paramsBlock.begin(), paramsBlock.end()));

  if (DEBUG_PCE)
  {
    if ( (paramsBlock.size() > 0) )
    {
      Xyce::Util::ParamList::const_iterator iter = paramsBlock.begin();
      int i=0;
      for ( ; iter != paramsBlock.end(); ++iter,++i)
      {
        const Xyce::Util::Param & par = *(iter);
        Xyce::lout() << "PCE::setDCOptions.  paramsBlock["<<i<<"] = " << par;
      }
    }
  }

  return true;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : PCEFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
///
/// Factory for parsing Step parameters from the netlist and creating Step analysis.
///
class PCEFactory : public Util::Factory<AnalysisBase, PCE>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : PCEFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 5/26/2018
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the PCE analysis factory
  ///
  /// @invariant Stores the results of parsing.  Multiple PCE analysis options may be
  /// applied and each generates and additional step.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  ///
  PCEFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Device::DeviceMgr &   device_manager,
    Linear::Builder &     builder,
    Loader::Loader &            loader,
    Topo::Topology &            topology,
    IO::InitialConditionsManager & initial_conditions_manager)
    : Util::Factory<AnalysisBase, PCE>(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      deviceManager_(device_manager),
      builder_(builder),
      loader_(loader),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager)
  {}

  virtual ~PCEFactory()
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
  /// Create a new PCE analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new PCE analysis object
  ///
  PCE *create() const
  {
    // don't have a mode yet
    //analysisManager_.setAnalysisMode(ANP_MODE_PCE);

    PCE *pce = new PCE(
        analysisManager_, 
    linearSystem_, nonlinearManager_, deviceManager_, builder_, 
        loader_, 
        topology_, 
        initialConditionsManager_,
        analysisManager_.getAnalysisObject());

    for (std::vector<Util::OptionBlock>::const_iterator it = samplingSweepAnalysisOptionBlock_.begin(), end = samplingSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      pce->setAnalysisParams(*it);
    }
    pce->setPCEOptions(samplingOptionBlock_);
    pce->setDCOptions(dcOptionBlock_);
    pce->setPCELinSol(pceLinSolOptionBlock_);
    pce->setLinSol(linSolOptionBlock_);

    return pce;
  }

  //-----------------------------------------------------------------------------
  // Function      : setPCEAnalysisOptionBlock
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
  void setPCEAnalysisOptionBlock(const Util::OptionBlock &option_block)
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
  // Function      : setPCEOptionBlock
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
  bool setPCEOptionBlock(const Util::OptionBlock &option_block)
  {
    samplingOptionBlock_ = option_block;
    return true;
  }

  bool setDCOptionBlock(const Util::OptionBlock &option_block)
  {
    dcOptionBlock_ = option_block;
    return true;
  }

  bool setPCELinSolOptionBlock(const Util::OptionBlock &option_block)
  {
    pceLinSolOptionBlock_ = option_block;
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
  Util::OptionBlock                 pceLinSolOptionBlock_;
  Util::OptionBlock                 linSolOptionBlock_;
};

//-----------------------------------------------------------------------------
// .SAMPLING
struct PCEAnalysisReg : public IO::PkgOptionsReg
{
  PCEAnalysisReg(
    PCEFactory &             factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setPCEAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  PCEFactory &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractPCEData
// Purpose       : Extract the parameters from a netlist .SAMPLING line held in
//                 parsed_line.
//
// Special Notes : used by the parser
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/26/2018
//-----------------------------------------------------------------------------
bool extractPCEData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("PCE", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  int parameterStartPos = 1;

  // Create an option block to temporarily store the default options.
  Util::OptionBlock defaultOptions;

  // Get the default options from metadata.
  addDefaultOptionsParameters(options_manager, defaultOptions, "PCE" );

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
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("PCE");

    parameters.insert(Util::ParamMap::value_type("USEEXPR", Util::Param("USEEXPR", true)));
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
    // Must use "PCES" instead of "PCE!"
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("PCES");

    parameters.insert(Util::ParamMap::value_type("NUMSAMPLES", Util::Param("NUMSAMPLES", 1)));
    parameters.insert(Util::ParamMap::value_type("COVMATRIX", Util::Param("COVMATRIX", "VECTOR")));

    parameters.insert(Util::ParamMap::value_type("OUTPUTFORMAT", Util::Param("OUTPUTFORMAT", "STD")));
    parameters.insert(Util::ParamMap::value_type("OUTPUTS", Util::Param("OUTPUTS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_ALL_SAMPLES", Util::Param("OUTPUT_ALL_SAMPLES", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_SAMPLE_STATS", Util::Param("OUTPUT_SAMPLE_STATS", true)));
#if Xyce_STOKHOS_ENABLE
    parameters.insert(Util::ParamMap::value_type("ORDER", Util::Param("ORDER", 4)));
#endif
    //parameters.insert(Util::ParamMap::value_type("SAMPLE_TYPE", Util::Param("SAMPLE_TYPE", 1))); // default=LHS
    parameters.insert(Util::ParamMap::value_type("SAMPLE_TYPE", Util::Param("SAMPLE_TYPE", "LHS"))); 
    parameters.insert(Util::ParamMap::value_type("SEED", Util::Param("SEED", 0)));

    parameters.insert(Util::ParamMap::value_type("RESAMPLE", Util::Param("RESAMPLE", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_PCE_COEFFS", Util::Param("OUTPUT_PCE_COEFFS", false)));
    parameters.insert(Util::ParamMap::value_type("SPARSE_GRID", Util::Param("SPARSE_GRID", false)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("STDOUTPUT", Util::Param("STDOUTPUT", false)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_STOCHASTIC_MATRIX", Util::Param("OUTPUT_STOCHASTIC_MATRIX", false)));
    parameters.insert(Util::ParamMap::value_type("VOLTLIM_ALG", Util::Param("VOLTLIM_ALG", 2)));
    parameters.insert(Util::ParamMap::value_type("STDOUTPUT_XVEC", Util::Param("STDOUTPUT_XVEC", false)));
  }


  // -------------------------------------------------------------------------------------------
  // copied over from LINSOL, with a few things added
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("LINSOL-PCE");

    parameters.insert(Util::ParamMap::value_type("AZ_max_iter", Util::Param("AZ_max_iter", 200)));
    parameters.insert(Util::ParamMap::value_type("AZ_precond", Util::Param("AZ_precond", 14)));
    parameters.insert(Util::ParamMap::value_type("AZ_solver", Util::Param("AZ_solver", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_conv", Util::Param("AZ_conv", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_pre_calc", Util::Param("AZ_pre_calc", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_keep_info", Util::Param("AZ_keep_info", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_orthog", Util::Param("AZ_orthog", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_subdomain_solve", Util::Param("AZ_subdomain_solve", 9)));
    parameters.insert(Util::ParamMap::value_type("AZ_ilut_fill", Util::Param("AZ_ilut_fill", 3.0)));
    parameters.insert(Util::ParamMap::value_type("AZ_drop", Util::Param("AZ_drop", 1.0E-3)));
    parameters.insert(Util::ParamMap::value_type("AZ_reorder", Util::Param("AZ_reorder", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_scaling", Util::Param("AZ_scaling", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_kspace", Util::Param("AZ_kspace", 50)));
    parameters.insert(Util::ParamMap::value_type("AZ_tol", Util::Param("AZ_tol", 1.0E-9)));
    parameters.insert(Util::ParamMap::value_type("AZ_output", Util::Param("AZ_output", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_diagnostics", Util::Param("AZ_diagnostics", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_overlap", Util::Param("AZ_overlap", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_rthresh", Util::Param("AZ_rthresh", 1.0001)));
    parameters.insert(Util::ParamMap::value_type("AZ_athresh", Util::Param("AZ_athresh", 1.0E-4)));
    parameters.insert(Util::ParamMap::value_type("AZ_filter", Util::Param("AZ_filter", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TR_partition", Util::Param("TR_partition", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_partition_type", Util::Param("TR_partition_type", "HYPERGRAPH")));
#ifdef Xyce_SHYLU
    parameters.insert(Util::ParamMap::value_type("ShyLU_rthresh", Util::Param("ShyLU_rthresh", 1.0E-3)));
#endif
    parameters.insert(Util::ParamMap::value_type("TR_reindex", Util::Param("TR_reindex", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_solvermap", Util::Param("TR_solvermap", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_amd", Util::Param("TR_amd", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_btf", Util::Param("TR_btf", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf", Util::Param("TR_global_btf", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf_droptol", Util::Param("TR_global_btf_droptol", 1.0E-16)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf_verbose", Util::Param("TR_global_btf_verbose", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_amd", Util::Param("TR_global_amd", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_amd_verbose", Util::Param("TR_global_amd_verbose", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_singleton_filter", Util::Param("TR_singleton_filter", 0)));
    parameters.insert(Util::ParamMap::value_type("SLU_EQUILIBRATE", Util::Param("SLU_EQUILIBRATE", 1)));
    parameters.insert(Util::ParamMap::value_type("SLU_REFACTOR", Util::Param("SLU_REFACTOR", 1)));
    parameters.insert(Util::ParamMap::value_type("SLU_PERMUTE", Util::Param("SLU_PERMUTE", 2)));
    parameters.insert(Util::ParamMap::value_type("SLU_PIVOT_THRESH", Util::Param("SLU_PIVOT_THRESH", -1.0)));
    parameters.insert(Util::ParamMap::value_type("SLU_FILL_FAC", Util::Param("SLU_FILL_FAC", -1)));
    parameters.insert(Util::ParamMap::value_type("BTF", Util::Param("BTF", 0)));
    parameters.insert(Util::ParamMap::value_type("BTF_VERBOSE", Util::Param("BTF_VERBOSE", 0)));
    parameters.insert(Util::ParamMap::value_type("BTF_ATHRESH", Util::Param("BTF_ATHRESH", 0.0)));
    parameters.insert(Util::ParamMap::value_type("BTF_RTHRESH", Util::Param("BTF_RTHRESH", 0.0)));
    parameters.insert(Util::ParamMap::value_type("BTF_RTHRESH_INIT", Util::Param("BTF_RTHRESH_INIT", 0.0)));
    parameters.insert(Util::ParamMap::value_type("BTF_INIT", Util::Param("BTF_INIT", 0)));
    parameters.insert(Util::ParamMap::value_type("BTF_THRESHOLDING", Util::Param("BTF_THRESHOLDING", 0)));
    parameters.insert(Util::ParamMap::value_type("BTF_RNTHRESHFAC", Util::Param("BTF_RNTHRESHFAC", 1.0e-3)));
    parameters.insert(Util::ParamMap::value_type("adaptive_solve", Util::Param("adaptive_solve", 0)));
    parameters.insert(Util::ParamMap::value_type("use_aztec_precond", Util::Param("use_aztec_precond", 1)));
    parameters.insert(Util::ParamMap::value_type("use_ifpack_factory", Util::Param("use_ifpack_factory", 0)));
    parameters.insert(Util::ParamMap::value_type("ifpack_type", Util::Param("ifpack_type", "Amesos")));
    parameters.insert(Util::ParamMap::value_type("diag_perturb", Util::Param("diag_perturb", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TR_rcm", Util::Param("TR_rcm", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale", Util::Param("TR_scale", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_left", Util::Param("TR_scale_left", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_right", Util::Param("TR_scale_right", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_exp", Util::Param("TR_scale_exp", 1.0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_iter", Util::Param("TR_scale_iter", 0)));
    parameters.insert(Util::ParamMap::value_type("TYPE", Util::Param("TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("PREC_TYPE", Util::Param("PREC_TYPE", "DEFAULT")));
#ifdef Xyce_BELOS
    parameters.insert(Util::ParamMap::value_type("BELOS_SOLVER_TYPE", Util::Param("BELOS_SOLVER_TYPE", "Block GMRES")));
#endif
    parameters.insert(Util::ParamMap::value_type("KLU_REPIVOT", Util::Param("KLU_REPIVOT", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_LS", Util::Param("OUTPUT_LS", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_BASE_LS", Util::Param("OUTPUT_BASE_LS", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_FAILED_LS", Util::Param("OUTPUT_FAILED_LS", 1)));

    parameters.insert(Util::ParamMap::value_type("DIRECT_SOLVER", Util::Param("DIRECT_SOLVER", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_LS", Util::Param("OUTPUT_LS", 1)));
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
bool registerPCEFactory(FactoryBlock & factory_block)
{
  PCEFactory *factory = new PCEFactory(
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

  factory_block.optionsManager_.addCommandParser(".PCE", extractPCEData);
  factory_block.optionsManager_.addCommandProcessor("PCE", new PCEAnalysisReg(*factory));

  // WARNING:  If you use "PCE" here for .options, it will 
  // confuse the "processor" code so using "PCES" instead as identifier.
  factory_block.optionsManager_.addOptionsProcessor("PCES", IO::createRegistrationOptions<PCEFactory>(*factory, &PCEFactory::setPCEOptionBlock)); 

  // DC options might not be needed
  factory_block.optionsManager_.addOptionsProcessor("DC", IO::createRegistrationOptions <PCEFactory> (*factory, &PCEFactory::setDCOptionBlock)); 

  factory_block.optionsManager_.addOptionsProcessor("LINSOL-PCE", IO::createRegistrationOptions <PCEFactory> (*factory, &PCEFactory::setPCELinSolOptionBlock));
  factory_block.optionsManager_.addOptionsProcessor("LINSOL", IO::createRegistrationOptions <PCEFactory> (*factory, &PCEFactory::setLinSolOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce

#endif // stokhos

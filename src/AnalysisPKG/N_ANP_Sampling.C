//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose       : Sampling Sweep class analysis functions.
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ANP_Sampling.h>
#include <N_ANP_StepEvent.h>

#include <N_ANP_UQSupport.h>

#include <N_ERH_Message.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_MeasureManager.h>

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
// Function      : Sampling::Sampling
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/1/2017
//-----------------------------------------------------------------------------
Sampling::Sampling(AnalysisManager &analysis_manager, Loader::Loader &loader, 
  Topo::Topology & topology, AnalysisBase &child_analysis)
    : AnalysisBase(analysis_manager, "Sampling"),
      analysisManager_(analysis_manager),
      loader_(loader),
      topology_(topology),
      outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
      measureManager_(outputManagerAdapter_.getMeasureManager()),
      childAnalysis_(child_analysis),
      samplingVector_(),
      maxParamStringSize_(0),
      lower_bounds_Given_(false),
      upper_bounds_Given_(false),
      covMatrixGiven_(false),
      numSamples_(1),
      numSamplesGiven_(false),
      sampleType_(UQ::LHS),
      userSeed_(0),
      userSeedGiven_(false),
      hackOutputFormat_("STD"),
#if Xyce_STOKHOS_ENABLE
      regressionPCEenable_(false),
      projectionPCEenable_(false),
      PCEorder_(4),
      resamplePCE_(false),
      outputPCECoeffs_(false),
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
      outputSampleStats_(true),
      useExpressionSamples_(false)
{
  pdsMgrPtr_ = analysisManager_.getPDSManager();
}


//-----------------------------------------------------------------------------
// Function      : Sampling::~Sampling
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/1/2017
//-----------------------------------------------------------------------------
Sampling::~Sampling()
{
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

  {
  std::vector<UQ::outputFunctionData*>::iterator iter = measFuncDataVec_.begin();
  std::vector<UQ::outputFunctionData*>::iterator end = measFuncDataVec_.end();
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
// Function      : Sampling::setAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::setAnalysisParams(const Util::OptionBlock & paramsBlock)
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
      if (stdDev < 0) { Report::DevelFatal() << "STD_DEVIATIONS values for .SAMPLING must be >= 0";}
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
      if (alpha <= 0) { Report::DevelFatal() << "ALPHA values for .SAMPLING must be > 0";}
      alphaVec_.push_back(alpha);
    }
    else if (std::string( iter->uTag() ,0,4) == "BETA")
    {
      betaGiven=true;
      double beta = iter->getImmutableValue<double>();
      if (beta <= 0) { Report::DevelFatal() << "BETA values for .SAMPLING must be > 0";}
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
    Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
    Parallel::Communicator &pdsComm = *(pds_manager.getPDSComm());

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
        Report::DevelFatal().in("parseSamplingParam") << "Unsupported SAMPLING type";
      }
      samplingVector_.push_back(sampling_param);
    }
  }

  outputManagerAdapter_.setStepSweepVector(samplingVector_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sampling::setSamplingOptions
// Purpose       :
// Special Notes : These are from '.options sampling'
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Sampling::setSamplingOptions(const Util::OptionBlock & option_block)
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
    else if (std::string((*it).uTag() ,0,7) == "OUTPUTS" ) 
    {
      outputsGiven_ = true;
      UQ::outputFunctionData * ofDataPtr = new UQ::outputFunctionData();
      ExtendedString funcName = (*it).stringValue();
      ofDataPtr->outFuncString = funcName.toUpper();
      outFuncDataVec_.push_back(ofDataPtr);
    }
    else if (std::string((*it).uTag() ,0,8) == "MEASURES" ) 
    {
      measuresGiven_ = true;
      UQ::outputFunctionData * ofDataPtr = new UQ::outputFunctionData();
      ExtendedString measName = (*it).stringValue();
      ofDataPtr->outFuncString = measName.toUpper();
      measFuncDataVec_.push_back(ofDataPtr);
    }
    else if ((*it).uTag() == "OUTPUT_SAMPLE_STATS")
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
      if (DEBUG_SAMPLING)
      {
        debugLevel_ = (*it).getImmutableValue<int>();
      }
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
    Report::UserWarning0() << "Neither output functions nor measures functions were specified";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sampling::getTIAParams() const
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
const TimeIntg::TIAParams & Sampling::getTIAParams() const
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : Sampling::getTIAParams()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
TimeIntg::TIAParams & Sampling::getTIAParams()
{
  return childAnalysis_.getTIAParams();
}

//-----------------------------------------------------------------------------
// Function      : Sampling::getDCOPFlag()
// Purpose       :
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::getDCOPFlag() const
{
  return childAnalysis_.getDCOPFlag();
}

//-----------------------------------------------------------------------------
// Function      : Sampling::doRun()
// Purpose       : This is the main controlling loop for sampling analysis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::doRun()
{
  return doInit() && doLoopProcess() && doFinish();
}

//-----------------------------------------------------------------------------
// Function      : Sampling::doInit()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::doInit()
{
  if (DEBUG_SAMPLING)
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "Sampling::init" << std::endl;
  }

  // check that all the specified params exist
#if 0
  std::cout << "Sampling::doInit.  Size of samplingVector = " << samplingVector_.size() <<std::endl;
  for (int ii=0;ii<samplingVector_.size();ii++)
  {
    std::cout << "samplingVector_["<<ii<<"].name = " << samplingVector_[ii].name << std::endl;
  }
#endif
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
    // Don't bother with this if projection PCE has been specified.
    long theSeed = UQ::getTheSeed(
        analysisManager_.getComm(), 
        analysisManager_.getCommandLine(), userSeed_, userSeedGiven_);
#if 0
Parallel::Machine comm = analysisManager_.getComm();
N_ERH_ErrorMgr::safeBarrier(comm);

Parallel::Manager &pds_manager = *analysisManager_.getPDSManager();
Parallel::Communicator &pdsComm = *(pds_manager.getPDSComm());
int myPID = pdsComm.procID();
int numProc = pdsComm.numProc();

for (int jj=0;jj<numProc;jj++)
{
  if (jj == myPID)
  {
    std::cout << "proc ID = " << myPID << " theSeed = " << theSeed << std::endl;
  }
}
#endif

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

  Util::publish<StepEvent>(analysisManager_,
      StepEvent(StepEvent::INITIALIZE, samplingVector_, numSamples_));

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Sampling::doLoopProcess()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::doLoopProcess()
{
  bool integration_status = true;

  if ( sampleType_ == UQ::MC)
  {
    Xyce::lout() << "***** Beginning Monte Carlo Sampling simulation....\n" << std::endl;
  }
  else if ( sampleType_ == UQ::LHS)
  {
    Xyce::lout() << "***** Beginning Latin Hypercube Sampling simulation....\n" << std::endl;
  }

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

  for (int i = 0; i < numSamples_; ++i)
  {
    // Tell the manager if any of our sweeps are being reset in this loop iteration.
    // ERK:  This reset boolean always is set to "false" - holdover from sweeps.
    // It should probably be always true in the sampling case.  Check this.
    bool reset = false;
    if (useExpressionSamples_)
    {
      reset = UQ::updateExpressionSamplingTerms(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numSamples_, false);
    }
    else
    {
      reset = UQ::updateSamplingParams(loader_, i, samplingVector_.begin(), samplingVector_.end(), Y_, numSamples_, false);
    }


    analysisManager_.setSweepSourceResetFlag(reset);

    outputManagerAdapter_.setStepSweepVector(samplingVector_);

    StepEvent step_event(StepEvent::STEP_STARTED, samplingVector_, i);
    Util::publish<StepEvent>(analysisManager_, step_event);

    // solve the loop.
    integration_status = childAnalysis_.run();

    step_event.state_ = StepEvent::STEP_COMPLETED;
    step_event.finalSimTime_ = getTIAParams().finalTime;
    Util::publish<StepEvent>(analysisManager_, step_event);

    // update the ensemble output functions 
    updateEnsembleOutputs ();
  }

  return integration_status;
}

//-----------------------------------------------------------------------------
// Function      : Sampling::doProcessSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::doProcessSuccessfulStep()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sampling::doProcessFailedStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::doProcessFailedStep()
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Sampling::doFinish()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool Sampling::doFinish()
{
  Util::publish<StepEvent>(analysisManager_, StepEvent(StepEvent::FINISH, samplingVector_, numSamples_));

  completeEnsembleOutputs();
  hackEnsembleOutput();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sampling::updateEnsembleOutputs
// Purpose       : This is essentially acting like a .RESULT, except that it 
//                 computes means, variances and std dev 
//                 of the quantities of interest.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Sampling::updateEnsembleOutputs()
{
  Parallel::Machine comm = analysisManager_.getComm();

  if (outputsGiven_)
  {
    IO::OutputMgr & output_manager = outputManagerAdapter_.getOutputManager();
    double circuit_time = analysisManager_.getStepErrorControl().nextTime;
    double circuit_dt = analysisManager_.getStepErrorControl().currentTimeStep;

    TimeIntg::DataStore * dsPtr = analysisManager_.getDataStore();
    Linear::Vector & solution_vector = *(dsPtr->nextSolutionPtr);
    Linear::Vector & state_vector = *(dsPtr->nextStatePtr);
    Linear::Vector & store_vector = *(dsPtr->nextStorePtr);

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

    //get expression value and derivatives, and compute means, etc
    for (int iout=0;iout<outFuncDataVec_.size();++iout)
    {
      UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
      Util::Op::OpData opDataTmp(0, &solution_vector, 0, &state_vector, &store_vector, 0);
      double val = 0.0;
      outFunc.expDataPtr->evaluate(comm, circuit_time, circuit_dt, opDataTmp, val);
      outFunc.sampleOutputs.push_back(val);
    }
  }

  if (measuresGiven_)
  {
    if (Parallel::rank(comm) == 0) // might not be necessary. For expData objects (above) evaluating on only proc=0 caused MPI errors.
    {
      // get measure value
      for (int iout=0;iout<measFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & measFunc = *(measFuncDataVec_[iout]);

        double val = 0.0;
        measureManager_.getMeasureValue(measFunc.outFuncString, val);
        measFunc.sampleOutputs.push_back(val);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Sampling::completeEnsembleOutputs
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Sampling::completeEnsembleOutputs()
{
  if (outputsGiven_)
  {
#if Xyce_STOKHOS_ENABLE
    // the seed is needed for resampling a PCE approximation
    long theSeed=0;
    if (resamplePCE_)
    {
      theSeed = UQ::getTheSeed( analysisManager_.getComm(), analysisManager_.getCommandLine(), userSeed_, userSeedGiven_);
    }
#endif

    Parallel::Machine comm = analysisManager_.getComm();
    if (Parallel::rank(comm) == 0)
    {
      for (int iout=0;iout<outFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & outFunc = *(outFuncDataVec_[iout]);
        outFunc.completeStatistics();

        if (stdOutputFlag_ && outputSampleStats_)
        {
          std::string sampleTypeStr = "MC";

          if ( sampleType_ == UQ::MC)
          {
            sampleTypeStr = "MC";
          }
          else
          {
            sampleTypeStr = "LHS";
          }

          outFunc.output(Xyce::lout(), sampleTypeStr);

          // histrogram is a hack that doesn't work yet
          //UQ::histrogram(Xyce::lout(), outFunc.outFuncString, outFunc.sampleOutputs);
        }
      }

#if Xyce_STOKHOS_ENABLE
      if (regressionPCEenable_)
      {
        // set up the lower-case "x" vector.  Capital "X" (or "Y") contains the 
        // actual parameter values used by the device models.
        //
        // Lower case "x" contains the unscaled values used by PCE.
        // So, for Hermite polynomials (normal dist), they must be converted to standard normals.
        // and for Legendre polynomials (uniform dist), they must be converted to domain and range over (-1,+1)
        //
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
          UQ::solveRegressionPCE( paramNameVec_.size(), x, f, regressionPCE);

          if (outputPCECoeffs_)
          {
            Xyce::lout() << "PCE coefs for " << outFunc.outFuncString << " = " << regressionPCE << std::endl;
            regressionPCE.print(Xyce::lout());
          }

          if (stdOutputFlag_)
          {
            Xyce::lout() << std::endl;
            Xyce::lout() << "(traditional sampling) regression PCE mean of " << outFunc.outFuncString << " = " << regressionPCE.mean() << std::endl;
            Xyce::lout() << "(traditional sampling) regression PCE stddev of " << outFunc.outFuncString << " = " << regressionPCE.standard_deviation() << std::endl;
          }

          if (resamplePCE_)
          {
            pceVec.push_back(regressionPCE);
          }
        } // for loop over outFuncDataVec_
     
        if (resamplePCE_)
        {
          std::vector < std::vector<double> > fvec (outFuncDataVec_.size());
          std::vector <UQ::statisticalMoments> statVec (outFuncDataVec_.size());

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
      } // regressionPCEenable_

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

          if (outputPCECoeffs_)
          {
            Xyce::lout() << "PCE coefs for " << outFunc.outFuncString << " = " << projectionPCE << std::endl;
            projectionPCE.print(Xyce::lout());
          }

          if (stdOutputFlag_)
          {
            Xyce::lout() << std::endl;
            Xyce::lout() << "(traditional sampling) projection PCE mean of " << outFunc.outFuncString << " = " << projectionPCE.mean() << std::endl;
            Xyce::lout() << "(traditional sampling) projection PCE stddev of " << outFunc.outFuncString << " = " << projectionPCE.standard_deviation() << std::endl;
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
    }
  }

  if (measuresGiven_)
  {
#if Xyce_STOKHOS_ENABLE
    // the seed is needed for resampling a PCE approximation
    long theSeed=0;
    if (resamplePCE_)
    {
      theSeed = UQ::getTheSeed( analysisManager_.getComm(), analysisManager_.getCommandLine(), userSeed_, userSeedGiven_);
    }
#endif

    Parallel::Machine comm = analysisManager_.getComm();
    if (Parallel::rank(comm) == 0)
    {
      for (int iout=0;iout<measFuncDataVec_.size();++iout)
      {
        UQ::outputFunctionData & measFunc = *(measFuncDataVec_[iout]);
        measFunc.completeStatistics();

        if (stdOutputFlag_ && outputSampleStats_)
        {
          std::string sampleTypeStr = "MC";

          if ( sampleType_ == UQ::MC)
          {
            sampleTypeStr = "MC";
          }
          else
          {
            sampleTypeStr = "LHS";
          }

          measFunc.output(Xyce::lout(), sampleTypeStr);

          // histrogram is a hack that doesn't work yet
          //UQ::histrogram(Xyce::lout(), measFunc.measFuncString, measFunc.measFuncEvalVec);
        }
      }

#if Xyce_STOKHOS_ENABLE
      if (regressionPCEenable_)
      {
        // set up the lower-case "x" vector.  Capital "X" (or "Y") contains the 
        // actual parameter values used by the device models.
        //
        // Lower case "x" contains the unscaled values used by PCE.
        // So, for Hermite polynomials (normal dist), they must be converted to standard normals.
        // and for Legendre polynomials (uniform dist), they must be converted to domain and range over (-1,+1)
        //
        const int numParams = paramNameVec_.size();
        std::vector< std::vector<double> > x;
        UQ::unScaleSampleValues(numSamples_, samplingVector_, covMatrix_, meanVec_, Y_, x);

        std::vector< Stokhos::OrthogPolyApprox<int,double> > pceVec;

        for (int iout=0;iout<measFuncDataVec_.size();++iout)
        {
          UQ::outputFunctionData & measFunc = *(measFuncDataVec_[iout]);

          Stokhos::OrthogPolyApprox<int,double> & regressionPCE = measFunc.regressionPCE;
          regressionPCE.reset(regrBasis);

          std::vector<double> & f = measFunc.sampleOutputs;
          UQ::solveRegressionPCE( paramNameVec_.size(), x, f, regressionPCE);

          if (outputPCECoeffs_)
          {
            Xyce::lout() << "PCE coefs for " << measFunc.outFuncString << " = " << regressionPCE << std::endl;
            regressionPCE.print(Xyce::lout());
          }

          if (stdOutputFlag_)
          {
            Xyce::lout() << std::endl;
            Xyce::lout() << "(traditional sampling) regression PCE mean of " << measFunc.outFuncString << " = " << regressionPCE.mean() << std::endl;
            Xyce::lout() << "(traditional sampling) regression PCE stddev of " << measFunc.outFuncString << " = " << regressionPCE.standard_deviation() << std::endl;
          }

          if (resamplePCE_)
          {
            pceVec.push_back(regressionPCE);
          }
        } // for loop over measFuncDataVec_
     
        if (resamplePCE_)
        {
          std::vector < std::vector<double> > fvec (measFuncDataVec_.size());
          std::vector <UQ::statisticalMoments> statVec (measFuncDataVec_.size());

          UQ::sampleApproximationPCE(theSeed, sampleType_, samplingVector_, covMatrix_, meanVec_, numResamples_, numParams, pceVec, fvec, statVec);

          for (int iout=0;iout<measFuncDataVec_.size();++iout)
          {
            UQ::outputFunctionData & measFunc = *(measFuncDataVec_[iout]);
            Xyce::lout() << std::endl;
            Xyce::lout() << "Statistics from re-sampling regression PCE approximation of " << measFunc.outFuncString << ":" << std::endl;;
            Xyce::lout() << statVec[iout];
            Xyce::lout() << std::endl;
          }
        }
      } // regressionPCEenable_

      if (projectionPCEenable_)
      {
        std::vector< Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > > pceVec;

        for (int iout=0;iout<measFuncDataVec_.size();++iout)
        {
          UQ::outputFunctionData & measFunc = *(measFuncDataVec_[iout]);

          Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = measFunc.projectionPCE;
          projectionPCE.init(0.0);
          projectionPCE.reset(quadExpn);

          std::vector<double> & f = measFunc.sampleOutputs;

          UQ::solveProjectionPCE(quadBasis, quadMethod, f, projectionPCE);

          if (outputPCECoeffs_)
          {
            Xyce::lout() << "PCE coefs for " << measFunc.outFuncString << " = " << projectionPCE << std::endl;
            projectionPCE.print(Xyce::lout());
          }

          if (stdOutputFlag_)
          {
            Xyce::lout() << std::endl;
            Xyce::lout() << "(traditional sampling) projection PCE mean of " << measFunc.outFuncString << " = " << projectionPCE.mean() << std::endl;
            Xyce::lout() << "(traditional sampling) projection PCE stddev of " << measFunc.outFuncString << " = " << projectionPCE.standard_deviation() << std::endl;
          }

          if (resamplePCE_)
          {
            pceVec.push_back(projectionPCE);
          }
        }

        if (resamplePCE_)
        {
          std::vector < std::vector<double> > fvec (measFuncDataVec_.size());
          std::vector <UQ::statisticalMoments> statVec (measFuncDataVec_.size());

          const int numParams = paramNameVec_.size();
          UQ::sampleApproximationPCE(theSeed, sampleType_, samplingVector_, covMatrix_, meanVec_, numResamples_, numParams, pceVec, fvec, statVec);

          for (int iout=0;iout<measFuncDataVec_.size();++iout)
          {
            UQ::outputFunctionData & measFunc = *(measFuncDataVec_[iout]);
            Xyce::lout() << std::endl;
            Xyce::lout() << "Statistics from re-sampling projection PCE approximation of " << measFunc.outFuncString << ":" << std::endl;;
            Xyce::lout() << statVec[iout];
            Xyce::lout() << std::endl;
          }
        }
      }
#endif

    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Sampling::hackEnsembleOutput ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void Sampling::hackEnsembleOutput ()
{
  if (outputsGiven_)
  {
    std::string fileName; 

    if (hackOutputFormat_=="TECPLOT")
    {
      fileName = analysisManager_.getNetlistFilename() + "_sampling.dat";
    }
    else if (hackOutputFormat_=="STD")
    {
      fileName = analysisManager_.getNetlistFilename() + "_sampling.prn";
    }
    else
    {
      Xyce::Report::UserWarning() << hackOutputFormat_ 
        << " is not a recognized sampling output option.\n" << std::endl;
    }
  }

  return;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : SamplingFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
///
/// Factory for parsing Step parameters from the netlist and creating Step analysis.
///
class SamplingFactory : public Util::Factory<AnalysisBase, Sampling>
{
public:
  //-----------------------------------------------------------------------------
  // Function      : SamplingFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 8/14/2017
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the Sampling analysis factory
  ///
  /// @invariant Stores the results of parsing.  Multiple Sampling analysis options may be
  /// applied and each generates and additional step.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  ///
  SamplingFactory(
    Analysis::AnalysisManager & analysis_manager,
    Linear::System &            linear_system,
    Nonlinear::Manager &        nonlinear_manager,
    Loader::Loader &            loader,
    Topo::Topology &            topology)
    : Util::Factory<AnalysisBase, Sampling>(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      topology_(topology)
  {}

  virtual ~SamplingFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 8/14/2017
  //-----------------------------------------------------------------------------
  ///
  /// Create a new Sampling analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new Sampling analysis object
  ///
  Sampling *create() const
  {
    Sampling *sampling = new Sampling(analysisManager_, loader_, topology_, analysisManager_.getAnalysisObject());
    for (std::vector<Util::OptionBlock>::const_iterator it = samplingSweepAnalysisOptionBlock_.begin(), end = samplingSweepAnalysisOptionBlock_.end(); it != end; ++it)
    {
      sampling->setAnalysisParams(*it);
    }
    sampling->setSamplingOptions(samplingOptionBlock_);

    return sampling;
  }

  //-----------------------------------------------------------------------------
  // Function      : setSamplingAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter, SNL
  // Creation Date : 8/14/2017
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Appends to any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setSamplingAnalysisOptionBlock(const Util::OptionBlock &option_block)
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
  // Function      : setSamplingOptionBlock
  // Purpose       : 
  // Special Notes :
  // Scope         : public
  // Creator       : Eric Keiter
  // Creation Date : 
  //-----------------------------------------------------------------------------
  ///
  /// Saves the sampling parsed option block.
  ///
  /// @invariant Overwrites any previously specified sampling option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setSamplingOptionBlock(const Util::OptionBlock &option_block)
  {
    samplingOptionBlock_ = option_block;
    return true;
  }

public:
  AnalysisManager &             analysisManager_;
  Linear::System &              linearSystem_;
  Nonlinear::Manager &          nonlinearManager_;
  Loader::Loader &              loader_;
  Topo::Topology &              topology_;

private:
  std::vector<Util::OptionBlock>        samplingSweepAnalysisOptionBlock_;
  Util::OptionBlock                     samplingOptionBlock_;
};

//-----------------------------------------------------------------------------
// .SAMPLING
struct SamplingAnalysisReg : public IO::PkgOptionsReg
{
  SamplingAnalysisReg(
    SamplingFactory &             factory)
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setSamplingAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  SamplingFactory &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractSamplingData
// Purpose       : Extract the parameters from a netlist .SAMPLING line held in
//                 parsed_line.
//
// Special Notes : used by the parser
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 8/14/2017
//-----------------------------------------------------------------------------
bool extractSamplingData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("SAMPLING", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  int parameterStartPos = 1;

  // Create an option block to temporarily store the default options.
  Util::OptionBlock defaultOptions;

  // Get the default options from metadata.
  addDefaultOptionsParameters(options_manager, defaultOptions, "SAMPLING" );

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
// Creation Date :
//-----------------------------------------------------------------------------
void populateMetadata(IO::PkgOptionsMgr & options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("SAMPLING");

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
    // Must use "SAMPLES" instead of "SAMPLING!"
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("SAMPLES");

    parameters.insert(Util::ParamMap::value_type("NUMSAMPLES", Util::Param("NUMSAMPLES", 1)));
    parameters.insert(Util::ParamMap::value_type("COVMATRIX", Util::Param("COVMATRIX", "VECTOR")));

    parameters.insert(Util::ParamMap::value_type("OUTPUTFORMAT", Util::Param("OUTPUTFORMAT", "STD")));
    parameters.insert(Util::ParamMap::value_type("OUTPUTS", Util::Param("OUTPUTS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MEASURES", Util::Param("MEASURES", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_SAMPLE_STATS", Util::Param("OUTPUT_SAMPLE_STATS", true)));
#if Xyce_STOKHOS_ENABLE
    parameters.insert(Util::ParamMap::value_type("REGRESSION_PCE", Util::Param("REGRESSION_PCE", true)));
    parameters.insert(Util::ParamMap::value_type("PROJECTION_PCE", Util::Param("PROJECTION_PCE", true)));
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
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
bool registerSamplingFactory(FactoryBlock & factory_block)
{
  SamplingFactory *factory = new SamplingFactory(
      factory_block.analysisManager_, 
      factory_block.linearSystem_, 
      factory_block.nonlinearManager_, 
      factory_block.loader_, 
      factory_block.topology_
      );

  addAnalysisFactory(factory_block, factory);

  populateMetadata(factory_block.optionsManager_);

  factory_block.optionsManager_.addCommandParser(".SAMPLING", extractSamplingData);
  factory_block.optionsManager_.addCommandProcessor("SAMPLING", new SamplingAnalysisReg(*factory));

  // WARNING:  If you use "SAMPLING" here for .options, it will 
  // confuse the "processor" code, 
  // mixing and matching .SAMPLING with .options SAMPLING
  // So, using "SAMPLES" instead as identifier.
  factory_block.optionsManager_.addOptionsProcessor("SAMPLES", IO::createRegistrationOptions<SamplingFactory>(*factory, &SamplingFactory::setSamplingOptionBlock));  

  return true;
}

} // namespace Analysis
} // namespace Xyce

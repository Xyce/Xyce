//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose       : This file contains the functions to manage measure objects
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <utility>
#include <sstream>

#include <N_ANP_fwd.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_MeasureExtrema.h>
#include <N_IO_MeasureFFT.h>
#include <N_IO_MeasureStats.h>
#include <N_IO_MeasureTranStats.h>
#include <N_IO_MeasureAverage.h>
#include <N_IO_MeasureDerivativeEvaluation.h>
#include <N_IO_MeasureDuty.h>
#include <N_IO_MeasureEquationEvaluation.h>
#include <N_IO_MeasureErrorFunctions.h>
#include <N_IO_MeasureFindWhen.h>
#include <N_IO_MeasureFourier.h>
#include <N_IO_MeasureFrequency.h>
#include <N_IO_MeasureIntegralEvaluation.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_MeasureMax.h>
#include <N_IO_MeasureMin.h>
#include <N_IO_MeasureOffTime.h>
#include <N_IO_MeasureOnTime.h>
#include <N_IO_MeasurePeakToPeak.h>
#include <N_IO_MeasureRMS.h>
#include <N_IO_MeasureError.h>
#include <N_IO_MeasureRiseFallDelay.h>
#include <N_IO_MeasureTrigTarg.h>
#include <N_IO_Remeasure.h>
#include <N_IO_NetlistImportTool.h>
#include <N_IO_OpBuilders.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_OutputPrn.h>
#include <N_IO_OutputCsd.h>
#include <N_IO_ParsingHelpers.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_PDS_ParHelpers.h>
#include <N_LAS_SystemHelpers.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Expression.h>

#include <expressionGroup.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : constructor
// Special Notes : The output of all measure values to both the  mt (or ms or
//                 ma) files and stdout is enabled by default in the constructor
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Manager::Manager(const CmdParse &cp)
  : commandLine_(cp),
    measureOutputFileSuffix_(".mt"),
    use_cont_files_(true),
    enableMeasGlobalPrint_(true),
    enableMeasGlobalVerbosePrint_(true),
    measDgt_(6),
    measDgtGiven_(false),
    measFail_(true),
    measOut_(true),
    measOutGiven_(false),
    useLTTM_(false),
    measGlobalDefaultVal_(0),
    measGlobalDefaultValGiven_(false),
    firstSweepValueFound_(false),
    startSweepValue_(0.0),
    endSweepValue_(0.0)
{}

//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Manager::~Manager()
{
  for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    delete (*it);
}

//-----------------------------------------------------------------------------
// Function      : Manager::makeMeasureOps
// Purpose       : Make the operators required for each measure
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : Unknown
//-----------------------------------------------------------------------------
void
Manager::makeMeasureOps(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager) 
{
  for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
    (*it)->makeMeasureOps(comm, op_builder_manager);
}

//-----------------------------------------------------------------------------
// Function      : Manager::notify
// Purpose       : Reset each measure at the start of a STEP, and output the measure
//                 value at the end of a STEP.  Output for the non-step case is
//                 currently handled by outputMacroResults().
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : Unknown
//-----------------------------------------------------------------------------
void
Manager::notify(
  const Analysis::StepEvent &   step_event)
{
  switch (step_event.state_) {
    case Analysis::StepEvent::INITIALIZE:
      break;

    case Analysis::StepEvent::STEP_STARTED:
      for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it) 
      {
        (*it)->reset();
      }
      // move any inactive measures back to the active list
      activeMeasuresList_ = allMeasuresList_;
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      {
        if( isMeasureActive() )
        {
          // do output to mt (or ms or ma) file and stdout, based on MEASPRINT option
          if (enableMeasGlobalPrint_)
	  {
            outputResultsToMTFile(step_event.count_);
          }
          // do stdout also
          if (enableMeasGlobalVerbosePrint_)
	  {
            outputVerboseResults(Xyce::lout(), step_event.finalSimTime_  );
          }
        }
      }

      break;

    case Analysis::StepEvent::FINISH:
      break;
  }
}


//-----------------------------------------------------------------------------
// Function      : Manager::addMeasure
// Purpose       : Entry point when .measure lines are passed in the netlist
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Manager::addMeasure(const Manager &measureMgr, const Util::OptionBlock & measureBlock )
{
  // based on what's in the option block passed in and any .OPTIONS MEASURE lines,
  // we create the needed measure instance
  if (DEBUG_IO)
    Xyce::dout() << "In Manager::addMeasure" << std::endl
                 << "measureLine passed it was: " << std::endl << measureBlock << std::endl;

  // here's the base data we should pull from the option block while
  std::string type;
  std::string mode;
  bool contMode = false;

  Util::ParamList::const_iterator it = std::find_if(measureBlock.begin(), measureBlock.end(), Util::EqualParam("TYPE"));
  if (it != measureBlock.end())
  {
    type = (*it).stringValue();
  }
  else
  {
    // this shouldn't happen, but catch it if does
    Report::UserError0() << "Missing TYPE in Measure Manager";
  }

  it = std::find_if(measureBlock.begin(), measureBlock.end(), Util::EqualParam("MODE"));
  if (it != measureBlock.end())
  {
    mode = (*it).stringValue();
    if (mode == "TRAN_CONT" || mode=="AC_CONT" || mode == "NOISE_CONT" || mode=="DC_CONT")
      contMode = true;
  }
  else
  {
    // this shouldn't happen, but catch it if does
    Report::UserError0() << "Missing MODE in Measure Manager";
  }

  Base  * theMeasureObject = 0;
  // have enough info to make the correct measure class
  if( type=="TRIG" || type=="TARG" )
  {
    bool fracMaxFound = false;
    it = std::find_if(measureBlock.begin(), measureBlock.end(), Util::EqualParam("FRAC_MAX"));
    if (it != measureBlock.end())
      fracMaxFound = true;

    if (contMode)
      theMeasureObject = new Measure::TrigTargCont(measureMgr, measureBlock);
    else if ( (mode == "TRAN") && (useLTTM_ || fracMaxFound) )
      theMeasureObject = new Measure::RiseFallDelay(measureMgr, measureBlock);
    else
      theMeasureObject = new Measure::TrigTarg(measureMgr, measureBlock);
  }
  else if( type=="AVG" )
  {
    theMeasureObject = new Measure::Average(measureMgr, measureBlock);
  }
  else if( type=="MAX" )
  {
    theMeasureObject = new Measure::Max(measureMgr, measureBlock);
  }
  else if( type=="MIN" )
  {
    theMeasureObject = new Measure::Min(measureMgr, measureBlock);
  }
  else if( type=="PP" )
  {
    theMeasureObject = new Measure::PeakToPeak(measureMgr, measureBlock);
  }
  else if( type=="RMS" )
  {
    theMeasureObject = new Measure::RMS(measureMgr, measureBlock);
  }
  else if( type=="FREQ" )
  {
    theMeasureObject = new Measure::Frequency(measureMgr, measureBlock);
  }
  else if( type=="DUTY" )
  {
    theMeasureObject = new Measure::Duty(measureMgr, measureBlock);
  }
  else if( type=="ON_TIME" )
  {
    theMeasureObject = new Measure::OnTime(measureMgr, measureBlock);
  }
  else if( type=="OFF_TIME" )
  {
    theMeasureObject = new Measure::OffTime(measureMgr, measureBlock);
  }
  else if( type=="FIND" || type=="WHEN" )
  {
    if (contMode)
      theMeasureObject = new Measure::FindWhenCont(measureMgr, measureBlock);
    else if (mode == "FFT")
      theMeasureObject = new Measure::FFTFind(measureMgr, measureBlock);
    else
      theMeasureObject = new Measure::FindWhen(measureMgr, measureBlock);
  }
  else if( type=="PARAM" || type=="EQN"  )
  {
    theMeasureObject = new Measure::EquationEvaluation(measureMgr, measureBlock);
  }
  else if( type=="DERIVATIVE" || type=="DERIV" )
  {
    if (contMode)
      theMeasureObject = new Measure::DerivativeEvaluationCont(measureMgr, measureBlock);
    else
      theMeasureObject = new Measure::DerivativeEvaluation(measureMgr, measureBlock);
  }
  else if( type=="INTEGRAL" || type=="INTEG" )
  {
    theMeasureObject = new Measure::IntegralEvaluation(measureMgr, measureBlock);
  }
  else if( type=="ERROR" )
  {
    theMeasureObject = new Measure::Error(measureMgr, measureBlock);
  }
  else if( type=="ERR1" || type=="ERR" )
  {
    theMeasureObject = new Measure::Err1(measureMgr, measureBlock);
  }
  else if( type=="ERR2" )
  {
    theMeasureObject = new Measure::Err2(measureMgr, measureBlock);
  }
  else if( type=="FOUR" )
  {
    theMeasureObject = new Measure::Fourier(measureMgr, measureBlock);
  }
  else if( type=="ENOB")
  {
    theMeasureObject = new Measure::ENOB(measureMgr, measureBlock);
  }
  else if( type=="SFDR")
  {
    theMeasureObject = new Measure::SFDR(measureMgr, measureBlock);
  }
  else if( type=="SNDR")
  {
    theMeasureObject = new Measure::SNDR(measureMgr, measureBlock);
  }
  else if( type=="SNR")
  {
    theMeasureObject = new Measure::SNR(measureMgr, measureBlock);
  }
  else if( type=="THD")
  {
    theMeasureObject = new Measure::THD(measureMgr, measureBlock);
  }
  else
  {
    // unknown type issue warning.
    Report::UserWarning0() << "Unknown MEASURE type requested, \"" << type << "\".  This measure will be ignored";
  }

  // if the measure object is supported, then add it to the active and all lists
  if (theMeasureObject && theMeasureObject->getTypeSupported() )
  {
    // Check for previous measure definition with this object's name.  If found then
    // remove the previous definitions, from all of the lists, and issue a warning message.
    int offset=0;
    for (MeasurementVector::iterator it = allMeasuresList_.begin(); it!=allMeasuresList_.end(); ++it, ++offset)
    {
      if (theMeasureObject->getMeasureName() == (*it)->getMeasureName())
      {
        // check both output lists
	MeasurementVector::iterator itOL;
        int offsetOL=0;
        for (itOL = measureOutputList_.begin(); itOL!=measureOutputList_.end(); ++itOL, ++offsetOL)
        {
          if (theMeasureObject->getMeasureName() == (*itOL)->getMeasureName())
	  {
            measureOutputList_.erase(measureOutputList_.begin()+offsetOL);
	    break;
          }
        }

        offsetOL=0;
        for (itOL = contMeasureOutputList_.begin(); itOL!=contMeasureOutputList_.end(); ++itOL, ++offsetOL)
        {
          if (theMeasureObject->getMeasureName() == (*itOL)->getMeasureName())
	  {
            contMeasureOutputList_.erase(contMeasureOutputList_.begin()+offsetOL);
            break;
          }
        }

        delete (*it);
        allMeasuresList_.erase(allMeasuresList_.begin()+offset);
        activeMeasuresList_.erase(activeMeasuresList_.begin()+offset);
        Report::UserWarning0() << "Measure \"" << theMeasureObject->getMeasureName() << "\" redefined, ignoring any previous definitions";

        break;
      }
    }

    allMeasuresList_.push_back( theMeasureObject );
    activeMeasuresList_.push_back( theMeasureObject );

    // Make list of continuous mode measures, if they will use separate output files.
    // The measureOutputList_ will use the <netlistName>.mtX files.
    if (contMode && use_cont_files_)
      contMeasureOutputList_.push_back( theMeasureObject );
    else
      measureOutputList_.push_back( theMeasureObject );

    // Used to help register lead current requests with device manager.
    // Measure manager keeps the combined list, based on parsing of
    // dependent solution variable vector for each measure.
    getLeadCurrentDevices(theMeasureObject->getDepSolVarIterVector(), devicesNeedingLeadCurrents_);
  }
  else
  {
    delete theMeasureObject;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::checkMeasureModes
// Purpose       : Used to check agreement between the analyis type and each 
//               : measure's mode.  Xyce will error out if they do not agree. 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 05/16/2017
//-----------------------------------------------------------------------------
bool Manager::checkMeasureModes(const Analysis::Mode analysisMode) const
{
  bool bsuccess = true;
  if ( isMeasureActive() )
  {
    // loop over all measure objects
    for (MeasurementVector::const_iterator it = allMeasuresList_.begin(), end = allMeasuresList_.end(); it != end; ++it)
    {
      // Check agreement between each measure's specified mode (from the .MEASURE line) and
      // the analysis type.  This could be condensed into one if statement but it seems more
      // readable with the conditions for analysis type broken out.
      if ( (((*it)->getMeasureMode() == "TRAN") || ((*it)->getMeasureMode() == "TRAN_CONT") || ((*it)->getMeasureMode() == "FFT")) &&
           (analysisMode != Xyce::Analysis::ANP_MODE_TRANSIENT) )
      {
        Report::UserError0() << "Netlist analysis statement and measure mode (" << (*it)->getMeasureMode()
                             << ") for measure " << (*it)->getMeasureName() << " do not agree";
        bsuccess = false;
      }
      else if ( (((*it)->getMeasureMode() == "AC") || ((*it)->getMeasureMode() == "AC_CONT")) && (analysisMode != Xyce::Analysis::ANP_MODE_AC) )
      {
        Report::UserError0() << "Netlist analysis statement and measure mode (" << (*it)->getMeasureMode()
                             << ") for measure " << (*it)->getMeasureName() << " do not agree";
        bsuccess = false;
      }
      else if ( (((*it)->getMeasureMode() == "DC") || ((*it)->getMeasureMode() == "DC_CONT")) && !(analysisMode == Xyce::Analysis::ANP_MODE_DC_SWEEP) )
      {
        Report::UserError0() << "Netlist analysis statement and measure mode (" << (*it)->getMeasureMode()
                             << ") for measure " << (*it)->getMeasureName() << " do not agree";
        bsuccess = false;
      }
      else if ( (((*it)->getMeasureMode() == "NOISE") || ((*it)->getMeasureMode() == "NOISE_CONT")) && !(analysisMode == Xyce::Analysis::ANP_MODE_NOISE) )
      {
        Report::UserError0() << "Netlist analysis statement and measure mode (" << (*it)->getMeasureMode()
                             << ") for measure " << (*it)->getMeasureName() << " do not agree";
        bsuccess = false;
      }
    }
  }

  return bsuccess ;
}

//-----------------------------------------------------------------------------
// Function      : Manager::fixupFFTMeasures
// Purpose       : Used to associate .MEASURE FFT lines with their .FFT line
// Special Notes : This function is not used for EQN/PARAM measures, which are
//                 treated like a TRAN mode measure even if entered on a
//                 .MEASURE FFT line.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 01/18/2021
//-----------------------------------------------------------------------------
void Manager::fixupFFTMeasures(Parallel::Machine comm, const FFTMgr& FFTMgr)
{
  if (!allMeasuresList_.empty())
  {
    const IO::FFTAnalysisVector FFTAnalysisList = FFTMgr.getFFTAnalysisList();
    IO::FFTAnalysisVector::const_reverse_iterator fftal_it;

    // loop over all measure objects
    for (MeasurementVector::const_iterator it = allMeasuresList_.begin(), end = allMeasuresList_.end(); it != end; ++it)
    {
      if ( ((*it)->getMeasureMode() == "FFT")  && ((*it)->getMeasureType() != "EQN") )
      {
        if (FFTAnalysisList.size() == 0)
        {
          Report::UserError0() << "No .FFT statement in netlist for measure " << (*it)->getMeasureName();
        }
        else
	{
          Util::Op::OpList::const_iterator ov_it = (*it)->getOutputVars()->begin();
	  std::string measureVarName = (*ov_it)->getName();
          // Need to adjust measure name, if it is an operator like VR(1), to remove the 'R'.
          // This is needed to make FIND measures work.
          size_t parenIdx = measureVarName.find_first_of('(');
          if ((measureVarName[0] != '{') && (parenIdx != 1) && isComplexCurrentOp(measureVarName,parenIdx))
            measureVarName = measureVarName[0] + measureVarName.substr(parenIdx);

          for (fftal_it = FFTAnalysisList.rbegin(); fftal_it != FFTAnalysisList.rend(); fftal_it++)
	  {
            if ( (*fftal_it)->getOutputVarName() == measureVarName )
	    {
              (*it)->fixupFFTMeasure(*fftal_it);
              break;
	    }
	  }
          if (fftal_it == FFTAnalysisList.rend())
	  {
            Report::UserError0() << "No matching .FFT statement with output variable "
                                 << measureVarName << " for measure " << (*it)->getMeasureName();
          }
        }
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setMeasureOutputFileSuffix
// Purpose       : The output files for different measure modes use different
//               : suffixes.  It's .mt for TRAN, .ma for AC and NOISE, and .ms for DC.
//               : The default of .mt was set in the constructor for the Measure 
//               : Manager.
// Special Notes : This function is used for both measure and remeasure
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 05/16/2017
//-----------------------------------------------------------------------------
void Manager::setMeasureOutputFileSuffix(const Analysis::Mode analysisMode)
{
  if ( analysisMode == Xyce::Analysis::ANP_MODE_TRANSIENT ) 
  {
    measureOutputFileSuffix_=".mt";
  }
  else if ( (analysisMode == Xyce::Analysis::ANP_MODE_AC) || (analysisMode == Xyce::Analysis::ANP_MODE_NOISE) )
  {
    measureOutputFileSuffix_=".ma";
  }  
  else if ( analysisMode == Xyce::Analysis::ANP_MODE_DC_SWEEP )   
  {
    measureOutputFileSuffix_=".ms";
  }
 
  return;
}

//-----------------------------------------------------------------------------
// Function      : Manager::updateTranMeasures
// Purpose       : Called during the simulation to update the measure objects
//               : for .TRAN analyses
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::updateTranMeasures(
  Parallel::Machine comm,
  double circuitTime,
  double endSimTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  // loop over active masure objects and get them to update themselves.
  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateTran(comm, circuitTime, endSimTime, solnVec, stateVec, storeVec,
                      lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), [](Measure::Base* meas){ return meas->finishedCalculation(); }),
                            activeMeasuresList_.end());
}


//-----------------------------------------------------------------------------
// Function      : Manager::updateDCMeasures
// Purpose       : Called during the simulation to update the measure objects
//               : for .DC analyses
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::updateDCMeasures(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  // Used in descriptive output to stdout. Store first/last values of first variable
  // found in the DC sweep vector, or the first/last table row indexes for a data sweep.
  if ( dcParamsVec.size() > 0 )
    recordStartEndSweepVals(getDCSweepVal(dcParamsVec));

  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateDC(comm, dcParamsVec, solnVec, stateVec, storeVec, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), [](Measure::Base* meas){ return meas->finishedCalculation(); }),
                            activeMeasuresList_.end());
}

//-----------------------------------------------------------------------------
// Function      : Manager::updateACMeasures
// Purpose       : Called during the simulation to update the measure objects
//               : for .AC analyses
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/21/2014
//-----------------------------------------------------------------------------
void Manager::updateACMeasures(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
  const Linear::Vector *real_solution_vector,
  const Linear::Vector *imaginary_solution_vector,
  const Util::Op::RFparamsData *RFparams)
{
  // Used in descriptive output to stdout. Store first/last frequency values
  recordStartEndSweepVals(frequency);

  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it) 
  {
    (*it)->updateAC(comm, frequency, fStart, fStop, real_solution_vector, imaginary_solution_vector, RFparams);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), [](Measure::Base* meas){ return meas->finishedCalculation(); }),
                            activeMeasuresList_.end()); }

//-----------------------------------------------------------------------------
// Function      : Manager::updateNoiseMeasures
// Purpose       : Called during the simulation to update the measure objects
//               : for .NOISE analyses
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 05/12/2020
//-----------------------------------------------------------------------------
void Manager::updateNoiseMeasures(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
  const Linear::Vector *real_solution_vector,
  const Linear::Vector *imaginary_solution_vector,
  double totalOutputNoiseDens,
  double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec
  )
{
  // Used in descriptive output to stdout. Store first/last frequency values
  recordStartEndSweepVals(frequency);

  for (MeasurementVector::iterator it = activeMeasuresList_.begin(); it != activeMeasuresList_.end(); ++it)
  {
    (*it)->updateNoise(comm, frequency, fStart, fStop, real_solution_vector, imaginary_solution_vector,
		       totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec);
  }
  activeMeasuresList_.erase(std::remove_if(activeMeasuresList_.begin(), activeMeasuresList_.end(), [](Measure::Base* meas){ return meas->finishedCalculation(); }),
                            activeMeasuresList_.end()); }

//-----------------------------------------------------------------------------
// Function      : Manager::outputResultsToMTFile
// Purpose       : Opens the .mt (or .ms or .ma) file, at end of simulation
//                 or at the end of each step.  It also makes the separate
//                 .mt (or .ms or .ma) files for the continuous measures, if
//                 requested.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/28/2018
//-----------------------------------------------------------------------------
void Manager::outputResultsToMTFile(int stepNumber) const
{
  // output non-continuous mode measures
  if ( !measureOutputList_.empty() )
  {
    // Make file name. The file suffix is mt for TRAN, ma for AC and ms for DC.
    std::string filename = IO::makeOutputFileNameWithStepNum(commandLine_, measureOutputFileSuffix_, stepNumber);

    // open file
    std::ofstream outputFileStream;
    outputFileStream.open( filename.c_str() );

    // output the measure values
    // loop over measure objects and get the results
    for (MeasurementVector::const_iterator it = measureOutputList_.begin(), end = measureOutputList_.end(); it != end; ++it)
    {
      // only output results to .mt file if measurePrintOption_ for that
      // measure is set to "ALL"
      if ( (*it)->getMeasurePrintOption() == "ALL")
      {
        (*it)->printMeasureResult( outputFileStream );
      }
    }

    // close file
    outputFileStream.close();
  }

  // output continuous mode measures
  if ( !contMeasureOutputList_.empty() )
  {
    for (MeasurementVector::const_iterator it = contMeasureOutputList_.begin(), end = contMeasureOutputList_.end(); it != end; ++it)
    {
      if ( (*it)->getMeasurePrintOption() == "ALL")
      {
        // Make file name. The file suffix is mt for TRAN, ma for AC and ms for DC.
        ExtendedString measNameLowerCase((*it)->getMeasureName());
        measNameLowerCase.toLower();
        std::string filename = IO::makeOutputFileNameWithStepNum(commandLine_,
                                     "_" + measNameLowerCase + measureOutputFileSuffix_,
                                     stepNumber);

        // open file
        std::ofstream outputFileStream;
        outputFileStream.open( filename.c_str() );

        // output the measure values
        (*it)->printMeasureResult( outputFileStream );

        // close file
        outputFileStream.close();
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Manager::outputAllResultsToLogFile
// Purpose       : Output measure results.  This function is intended to output
//                 all of the measure values to the log file.  It is most useful
//                 for parallel debugging with the -per-processor command line
//                 option.
// Special Notes : outputVerboseResults() is used to output measure info to stdout.
//                 outputResultsToMTFile also handles opening/closing the .mt file.
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void Manager::outputAllResultsToLogFile() const
{
  if ( isMeasureActive() )
  {
    // loop over measure objects and get the results
    for (MeasurementVector::const_iterator it = allMeasuresList_.begin(), end = allMeasuresList_.end(); it != end; ++it)
    {
      // only output results to .mt file if measurePrintOption_ for that
      // measure is set to "ALL"
      if ( (*it)->getMeasurePrintOption() == "ALL")
      {
        (*it)->printMeasureResult(Xyce::lout());
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : Manager::outputVerboseResults
// Purpose       : Output measure results, to std output, at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystem Modeling
// Creation Date : 02/05/2015
//-----------------------------------------------------------------------------
std::ostream &Manager::outputVerboseResults( std::ostream& outputStream, double endSimTime ) const
{
  if ( isMeasureActive() )
  {
    outputStream << std::endl
                 << " ***** Measure Functions ***** " << std::endl
                 << std::endl;

    // loop over active measure objects and get the results
    for (MeasurementVector::const_iterator it = allMeasuresList_.begin(), end = allMeasuresList_.end(); it != end; ++it)
    {
      // only output results to stdout file if measurePrintOption for that measure is 
      // set to "ALL" or "STDOUT"
      if ( ((*it)->getMeasurePrintOption() == "ALL") || ((*it)->getMeasurePrintOption() == "STDOUT") )
      { 
        // these two functions print out diagnostic information if the measure failed.  An example
        // is a nonsensical time window.  AT qualifier can either a time or a frequency.  So, it
	// gets its own function.
        (*it)->printMeasureWarnings( endSimTime, startSweepValue_, endSweepValue_ );
        (*it)->printMeasureWarningsForAT(endSimTime);
        // this function prints the measure value, and measure times if appropriate for a given measure such as
        // MIN, MAX, PP, TRIG/TARG
        (*it)->printVerboseMeasureResult(outputStream);
        // this function print out information about the measurement window used.  The sweep values
        // refer to the AC, DC or NOISE sweep values.
        (*it)->printMeasureWindow( outputStream, endSimTime, startSweepValue_, endSweepValue_ );
        // this function prints out information about the time window used for Rise/Fall/Cross
        // if one of those keywords was specified
        (*it)->printRFCWindow( outputStream );
        outputStream << std::endl; // separate verbose output by one line
      }
    }
  }
  
  return outputStream;
}

//-----------------------------------------------------------------------------
// Function      : Manager::recordStartEndACDCNoiseSweepVals
// Purpose       : Used to record start/end sweep values for AC, DC and NOISE measures.
// Special Notes : For AC and NOISE measures, sweepVal is frequency. For DC measures, it
//                 is the value of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 06/30/2020
//-----------------------------------------------------------------------------
void Manager::recordStartEndSweepVals(double sweepVal)
{
  if (!firstSweepValueFound_)
  {
    startSweepValue_ = sweepVal;
    firstSweepValueFound_ = true;
  }
  endSweepValue_ = sweepVal;

  return;
}

//-----------------------------------------------------------------------------
// Function      : Manager::getMeasureValue
// Purpose       : Get the value of a .measure test
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/29/2010
//-----------------------------------------------------------------------------
bool Manager::getMeasureValue(const std::string &name, double &value) const
{
  for (MeasurementVector::const_iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
  {
    if (equal_nocase((*it)->getMeasureName(), name))
    {
      value = (*it)->getMeasureResult();
      return true;
    }
  }
  return false;
}

const Base *Manager::find(const std::string &name) const
{
  for (MeasurementVector::const_iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it)
  {
    if (equal_nocase((*it)->getMeasureName(), name))
    {
      return *it;
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Manager::remeasure
// Purpose       : Recompute measure functions based on existing Xyce output
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 4/8/13
//-----------------------------------------------------------------------------
void Manager::remeasure(
  Parallel::Communicator &pds_comm,
  const std::string &netlist_filename,
  const std::string &remeasure_path,
  const char& analysisName,
  Util::Op::BuilderManager &op_builder_manager,
  OutputMgr &output_manager,
  FFTMgr &fft_manager,
  Analysis::AnalysisManager &analysis_manager,
  Analysis::AnalysisCreatorRegistry &analysis_registry,
  Util::SymbolTable &symbol_table)
{
  // check that remeasure is supported for the specified analysis mode.  Note that
  // the analysis mode might be 'O' if .OP was specified and it comes before the
  // .TRAN (or .AC or .DC) statement in the netlist.
  RemeasureBase* remeasureObj(0);
  if ( analysisName == 'O' ) 
  {
    Report::UserFatal0() << "Please comment out the .OP statement to get remeasure to work";
    return;
  }
  else if (analysisName == 'A')
  {
    remeasureObj = new RemeasureAC(pds_comm,*this,output_manager,analysis_manager,
                                   analysis_registry);
  }
  else if (analysisName == 'D')
  {
    remeasureObj = new RemeasureDC(pds_comm,*this,output_manager,analysis_manager,
                                   analysis_registry);
  }
  else if (analysisName == 'T')
  {
    remeasureObj = new RemeasureTRAN(pds_comm,*this,output_manager,analysis_manager,
                                     analysis_registry);
  }
  else
  {
    // unsupported analysis type.  Report error and exit remeasure routine
    Report::UserFatal0() << "Remeasure only supports .TRAN, .AC and .DC analysis modes";
    return;
  } 

  // output "re-measure header text" to stdout
  Xyce::lout() << "In OutputMgr::remeasure " << std::endl
               << "file to reprocess through measure and/or fft functions: " << remeasure_path << std::endl;

  // Some error checking on the file name specified by -remeasure command-line option.  
  // Will check for a valid file extension (e.g., .PRN) below.
  if ( !(Util::checkIfValidFile(remeasure_path)) )
  {
    // Error out if the user-specified file does not exist, cannot be opened,
    // or is a directory name rather than a file name.  See SON Bugs 730 
    // and 785 for more details.
    Report::UserFatal0() << "Could not find remeasure file: " << remeasure_path;
    return;
  }
  else if ( remeasure_path.length() < 4 )
  {
    // the ExtendedString fileExt() line below gives a "cryptic error message" without this check
    Report::UserFatal0() << "Remeasure filename must end in .PRN, .CSV or .CSD";
    return;
  }

  // Open file for reading.
  // Support PRN, CSV and CSD format now.  The function logic is to first test for a valid extension.
  // If the extension is invalid then the function errors out.  If the extension is valid, then
  // make the correct object type and try to open the file.
  OutputFileBase * fileToReMeasure;
  ExtendedString fileExt(remeasure_path.substr(remeasure_path.length()-4));
  // The Xyce format for a .CSV file is basically the same as a .PRN file made with the
  // -delim COMMA command line option.  The only difference is that the .PRN file has something
  // like "End of Xyce(TM) Simulation" as its last line.  That difference has no effect 
  // on re-measure.  So, re-measure of .PRN and .CSV files currently use the same function.
  if ( (fileExt.toUpper() == ".PRN")  || (fileExt.toUpper() == ".CSV") )
  {
    fileToReMeasure = new OutputPrn();
  }
  else if (fileExt.toUpper() == ".CSD")
  {
    fileToReMeasure = new OutputCsd();
  }
  else
  {
    // unsupported file format.  Report error and exit remeasure routine
    Report::UserFatal0() << "Remeasure only supports .PRN, .CSV or .CSD formats";
    return;
  }

  if (!fileToReMeasure->openFileForRead(remeasure_path))
  {
    // open failed.  Report error and exit remeasure routine
    Report::UserFatal0() << "Could not open remeasure file: " << remeasure_path;
    return;
  }

  // load data-names & allocate space for a line
  std::vector<std::string> fileVarNames;
  if (!(fileToReMeasure->getOutputVarNames(fileVarNames)))
  {
    // reading var names failed.  report error and exit remeasure routine.
    Report::UserFatal0() << "Problem reading variable names in remeasure file";
    return;
  }

  // Get column names for building up a symbol table.  Use the solSymIdx and brVarSymIdx 
  // vectors to track which column names will be treated as solution variables and which
  // columns will be treated as lead currents.  This approach handles V() and I() operators.
  std::vector<int> solSymIdx, brVarSymIdx;
  fileToReMeasure->convertOutputNamesToSolVarNames(fileVarNames, solSymIdx, brVarSymIdx);
  
  int numVars = fileVarNames.size();
  int numLocalVars = 0;

  // create an Linear::Vector to hold the data from the file.  This is
  // needed to support getPrintValue_ for evaluating the measure functions.
  std::vector<int> lbMap;
  if (Parallel::rank(pds_comm.comm()) == 0)
  {
    numLocalVars = numVars;
  }
  else
  {
    numLocalVars = 0;
  }

  // set up lbMap to map variable names to indices in the solution vector
  lbMap.resize(numLocalVars);
  for (int i=0;i<numLocalVars; i++)
  {
    lbMap[i]=i;
  }

  // Columns with I() operators are entered as "lead currents".  This is true even for
  // V and L devices, where those values are actually solution variables.
  for (std::vector<int>::iterator it = brVarSymIdx.begin(); it!= brVarSymIdx.end(); it++)
  {
    symbol_table[Util::BRANCH_SYMBOL][fileVarNames[*it]] = (Parallel::rank(pds_comm.comm()) == 0) ? *it : -1;
    remeasureObj->setIndepVarCol(Parallel::rank(pds_comm.comm()), *it, fileVarNames[*it]);
  }

  // All other columns (e.g., with V() operators or expressions) are entered as solution variables
  for (std::vector<int>::iterator it = solSymIdx.begin(); it!= solSymIdx.end(); it++)
  {
    symbol_table[Util::SOLUTION_SYMBOL][fileVarNames[*it]] = (Parallel::rank(pds_comm.comm()) == 0) ? *it : -1;
  
    // while scanning fileVarNames look for "TIME" or "FREQ" as a name in the 0th or 1st
    // column.  The column will be used later to "sense" when a step (caused by a .STEP line) has 
    // occurred in the data, for remeasure for TRAN and AC.  For remeasure for DC, the Index
    // column is used.
    remeasureObj->setIndepVarCol(Parallel::rank(pds_comm.comm()), *it, fileVarNames[*it]);
  }

  int IndepVarCol = remeasureObj->getIndepVarCol();

  // This function call will throw a FatalError if the required column (FREQ, Index or TIME)
  // was not found in the remeasured file.
  remeasureObj->checkIndepVarCol(Parallel::rank(pds_comm.comm()), IndepVarCol);

  // call make_ops to get the measure functions ready to get data
  makeMeasureOps(pds_comm.comm(), op_builder_manager);
  if ( !checkMeasureModes(remeasureObj->getAnalysisMode()) )
  {
     Report::UserFatal0() << "Error making measures during remeasure";
     return;
  }
  setMeasureOutputFileSuffix(remeasureObj->getAnalysisMode());

  // support remeasure of both FFT lines and FFT measures
  fft_manager.enableFFTAnalysis(remeasureObj->getAnalysisMode());
  if (fft_manager.isFFTActive())
  {
    // StepErrorControl is not made by the AnalysisManager during -remeasure
    TimeIntg::StepErrorControl* sec(NULL);
    fft_manager.fixupFFTParametersForRemeasure(pds_comm.comm(), op_builder_manager,
				   analysis_manager.getFinalTimeForRemeasure(), *sec);
    fixupFFTMeasures(pds_comm.comm(), fft_manager);
  }

  // this safeBarrier will cause remeasure to error out if any of the MeasureOps are invalid.
  Report::safeBarrier(pds_comm.comm());

  Teuchos::RCP<Parallel::ParMap> aParMap = Teuchos::rcp( Parallel::createPDSParMap(numVars, numLocalVars, lbMap, 0, pds_comm) );
  Teuchos::RCP<Linear::Vector> varValuesVec = Teuchos::rcp( Linear::createVector(*aParMap) );
  varValuesVec->putScalar(0);
  
  // Variables used to handle a STEP in the re-measured data.  For .DC, we need to 
  // handle "steps" from both the .DC line and the .STEP line.  That is done with
  // both the dcStepCount variable defined in the RemeasureDC object and the 
  // remeasureStepCount variable defined below.
  double prevIndepVar=0;  // used to sense steps from .STEP line
  int remeasureStepCount=0; // current step number, for steps caused by a .STEP line

  // For .DC, the steps caused by the .DC and .STEP lines are sensed as follows.
  // Consider this combo of .DC, .STEP and .PRINT lines (where the voltages for
  // V1, V2 and V3 do NOT have to appear in the output .PRN file):
  //
  //   .DC V1 1 2 1 V2 3 4 1
  //   .STEP V3 5 6 1
  //   .PRINT DC V1:DCV0 V2:DCV0 V3:DCV0
  //
  // The .PRN file will then have:
  // 
  // Index    V1:DCV0           V2:DCV0           V3:DCV0
  // 0        1.000e+00         3.000e+00         5.000e+00
  // 1        2.000e+00         3.000e+00         5.000e+00
  // 2        1.000e+00         4.000e+00         5.000e+00
  // 3        2.000e+00         4.000e+00         5.000e+00
  // 0        1.000e+00         3.000e+00         6.000e+00
  // 1        2.000e+00         3.000e+00         6.000e+00
  // 2        1.000e+00         4.000e+00         6.000e+00
  // 3        2.000e+00         4.000e+00         6.000e+00
  //
  // So, the steps caused by the .STEP line can be sensed from the INDEX column.
  // Since VSRC1:DCV0 does not have to appear in the output file, the steps caused
  // by the first variable on the .DC line will be sensed with the combo of the  
  // dcStepCountvariable (which counts position in the current step loop for V1, which 
  // is 0 or 1 in this case) and the dcParamsVec[0].maxStep variable which is the 
  // total number of steps (e.g., 2 in this case) for V1.
  
  // run though lines in the re-measured output file calling update measure as we go.
  int reading[1], summedReading[1];
  reading[0] = 0;
  summedReading[0] = 1;

  while (summedReading[0]==1)
  {
    if (Parallel::rank(pds_comm.comm()) == 0)
    {  
      reading[0] = fileToReMeasure->getOutputNextVarValuesParallel(varValuesVec.get());
      summedReading[0] = reading[0];
    }
    
    // reduce all procs on reading[] so all procs know when we are done.
    pds_comm.sumAll( reading, summedReading,  1);
    
    if( summedReading[0] == 1 )
    {
        double lt[1], sumlt[1];
        lt[0] = 0.0;
        sumlt[0] = 0.0;
        if( numLocalVars > 0 )
        {
          lt[0] = (*varValuesVec)[IndepVarCol];
        }
        pds_comm.sumAll( lt, sumlt,  1);

        remeasureObj->setIndepVar(sumlt[0]);

        // Xyce's "arrow of time" is always increasing.  So, sense a .STEP in the
        // re-measured data if the previous time (or frequency) variable is greater than the current
        // time variable.  The processing sequence is then:
        //    1) output the measure data, since it hasn't been updated or reset yet.
        //    2) reset all of the measures
        //    3) move any inactive measures back to the "active list"
        //    4) increment the step count (remeasureStepCount)
        //
        // For .DC output, the Index column will reset to 0, each time a "step",
        // based on the .STEP line, occurs.
        //
        // Note that this if statement is not processed for the last step
        // (or if no .STEP statement is used).  In that case, the data is
        // output near the end of this function.
	if ( prevIndepVar > remeasureObj->getIndepVar() )
        {
          if (fft_manager.isFFTActive() || !allMeasuresList_.empty())
             Xyce::lout() << "***** Processing Data for Step " << remeasureStepCount << std::endl;

          // This also re-calculates all of the .FFT analyses for the current step
          fft_manager.outputResultsToFFTfile(remeasureStepCount);
          fft_manager.outputVerboseResults(Xyce::lout());

          // Output the measure info to both mt file and stdout.
          // At present, .OPTIONS MEASURE MEASPRINT does not apply to -remeasure
          if( !allMeasuresList_.empty() )
          {
            outputResultsToMTFile(remeasureStepCount);
            outputVerboseResults(Xyce::lout(), prevIndepVar);
          }

          // reset all of the measures and FFT analyses (if any)
          fft_manager.resetFFTAnalyses();
          for (MeasurementVector::iterator it = allMeasuresList_.begin(); it != allMeasuresList_.end(); ++it) 
          {
            (*it)->reset();
          }

          // move any inactive measures back to the active list
          activeMeasuresList_ = allMeasuresList_;

          // Increment the step count, since we are not getting this info
          // from the analysis manager during re-measure.  Do this after outputting the
          // .mt, .ma or .ms file, since the current measured data is for the previous step.
          remeasureStepCount++;

          // Reset current value of first sweep variable on .DC line, at the start of each step.
          // This is currently only used by DC remeasure.
          remeasureObj->resetSweepVars();
        }
        prevIndepVar = remeasureObj->getIndepVar(); 

        // update each FFT analysis object, and then each measure
        fft_manager.updateFFTData(pds_comm.comm(), remeasureObj->getIndepVar(),
                                  varValuesVec.get(), 0, 0, varValuesVec.get(), 0, 0);
        remeasureObj->updateMeasures(*varValuesVec);

        // Update the Sweep Vector.  This is currently only used by DC remeasure.
        remeasureObj->updateSweepVars();  
      }

      varValuesVec->putScalar(0);
  }

  // Output the Measure results to Xyce::dout().
  // Make sure the function gets called on all processors, but only one outputs it.
  // output step count in stdout if .STEP was used.
  if (remeasureStepCount > 0)
  {
    Xyce::lout() << "***** Processing Data for Step " << remeasureStepCount << std::endl;
  }

  // This also re-calculates the FFT analyses
  fft_manager.outputResultsToFFTfile(remeasureStepCount);
  fft_manager.outputVerboseResults(Xyce::lout());

  // Output the Measure results to file.
  if (Parallel::rank(pds_comm.comm()) == 0)
  {
    // If .STEP is not used, then remeasureStepCount will be "0", and the output
    // file extension will be (for TRAN) .mt0.  Otherwise this IF statement will output the data
    // for the last step, with the appropriate number (X) in .mtX (or .maX or .mSX) for the file name.
    outputResultsToMTFile(remeasureStepCount);
  }
  else 
  {
    outputAllResultsToLogFile();
  }

  // At present, .OPTIONS MEASURE MEASPRINT does not affect -remeasure
  if (Parallel::rank(pds_comm.comm()) == 0)
  {
    outputVerboseResults( Xyce::lout(), output_manager.getCircuitTime() );
  }
  else
  {
    outputVerboseResults( Xyce::lout(), output_manager.getCircuitTime() );
  }

  fileToReMeasure->closeFileForRead();

  // clean-up pointers
  delete remeasureObj;
  remeasureObj=0;
  delete fileToReMeasure;
  fileToReMeasure = 0;
}

//-----------------------------------------------------------------------------
// Function      : Manager::registerMeasureOptions
// Purpose       : registers set of variables to set for .OPTIONS MEASURE
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/24/2018
//-----------------------------------------------------------------------------
bool Manager::registerMeasureOptions(const Util::OptionBlock &option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; )
  {
    if ( (*it).tag() == "USE_CONT_FILES")
    {
      use_cont_files_ = (*it).getImmutableValue<int>();
      ++it;
    }
    else if ((*it).tag() == "MEASPRINT")
    {
      // Note both of the variables were set to true, by default, in the
      // constructor for the measure manager.  This option is specific to Xyce.
      std::string measPrintStr = (*it).getImmutableValue<std::string>();
      if ( measPrintStr == "ALL" )
      {
        enableMeasGlobalPrint_ = true;
        enableMeasGlobalVerbosePrint_ = true;
      }
      else if (measPrintStr == "STDOUT")
      {
        enableMeasGlobalPrint_ = false;
        enableMeasGlobalVerbosePrint_ = true;
      }
      else if (measPrintStr == "NONE")
      {
        enableMeasGlobalPrint_ = false;
        enableMeasGlobalVerbosePrint_ = false;
      }
      else
      {
        // default to ALL case, for an unrecognized value
        Report::UserWarning0() << "Unknown option value " << measPrintStr <<
	                           " ignored for .OPTIONS MEASURE MEASPRINT";
        enableMeasGlobalPrint_ = true;
        enableMeasGlobalVerbosePrint_ = true;
      }

      // need to point at next parameter
      ++it;
    }
    else if ((*it).tag() == "MEASOUT")
    {
      // This controls whether the .mt0 files is made.  It takes precedence
      // over MEASPRINT. 
      int measOutIntVal = (*it).getImmutableValue<int>();
      measOutGiven_ = true;
      if ( measOutIntVal == 0 )
      {
	measOut_ = false;
      }
      else if ( measOutIntVal == 1 )
      {
	measOut_ = true;
      }
      else
      {
	// default to true, for any other value
        Report::UserWarning0() << ".OPTIONS MEASURE MEASOUT value must be 0 or 1. Setting value to 1";
        measOut_ = true;
      }
      // need to point at next parameter
      ++it;
    }
    else if ((*it).tag() == "MEASDGT")
    {
      // this control the precision of the .MEASURE output to stdout 
      // and the .mt0 (or .ma0 or .ms0) file
      measDgt_ = (*it).getImmutableValue<int>();
      measDgtGiven_ = true;
      // need to point at next parameter
      ++it;
    }
    else if ((*it).tag() == "MEASFAIL")
    {
      // this control the value output to the the .mt0 (or .ma0 or .ms0) file
      // for a failed measure
      int measFailIntVal = (*it).getImmutableValue<int>();
      if ( measFailIntVal == 0 )
      {
	measFail_ = false;
      }
      else if ( measFailIntVal == 1 )
      {
	measFail_ = true;
      }
      else
      {
        // default to 1, for any other value
        Report::UserWarning0() << ".OPTIONS MEASURE MEASFAIL value must be 0 or 1. Setting value to 1";
        measFail_ = true;
      }
      // need to point at next parameter
      ++it;
    }
    else if ((*it).tag() == "USE_LTTM")
    {
      // this control the value output to the the .mt0 (or .ma0 or .ms0) file
      // for a failed measure
      int lttmVal = (*it).getImmutableValue<int>();
      if ( lttmVal == 0 )
      {
	useLTTM_ = false;
      }
      else if ( lttmVal == 1 )
      {
	useLTTM_ = true;
      }
      else
      {
        // default to true, for any other value
        Report::UserWarning0() << ".OPTIONS MEASURE USE_LTTM value must be 0 or 1. Setting value to 1";
        useLTTM_ = true;
      }
      // need to point at next parameter
      ++it;
    }
    else if ((*it).tag() == "DEFAULT_VAL")
    {
      // This sets the default value for all failed measures.  This
      // option is specific to Xyce.
      measGlobalDefaultVal_ = (*it).getImmutableValue<double>();
      measGlobalDefaultValGiven_ = true;
      // need to point at next parameter
      ++it;
    }
    else
    {
      // silently ignore?
      ++it;
    }
  }

  // consolidate info from MEASPRINT and MEASOUT now that all of the option
  // blocks have been processed.  The HSPICE MEASOUT value takes precedence.
  if (measOutGiven_)
  {
    enableMeasGlobalPrint_ = measOut_;
  }

  return true;
}

namespace {

//-----------------------------------------------------------------------------
// Class         : MeasureOptionsReg
// Purpose       : functor for registering Measure options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct MeasureOptionsReg : public PkgOptionsReg
{
  MeasureOptionsReg( Measure::Manager &mgr )
    : measureManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { 
    return measureManager_.addMeasure( measureManager_ , options );
  }

  Measure::Manager &measureManager_;
};

//-----------------------------------------------------------------------------
// Function      : Manager::handleVALQualifier
// Purpose       : HSPICE allows TRIG and TARG syntax of the form:
//               :       TRIG {exp} VAL=0.1 
//               :  
//               : This function skips past the VAL= part if it exists, 
//               : and parsing of the .MEASURE line is in a TRIG or TARG block.
//               : It also catches the syntax error of VAL, without the = sign,
//               : for that case. 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 02/13/2019
//-----------------------------------------------------------------------------
bool handleValQualifier(const std::string& netlist_filename,
  int numFields,
  const std::string & parsedType,
  const IO::TokenVector & parsed_line,
  int& position)
{

  bool retval = true;
  if ( ((position+1) < numFields) && (parsedType=="TRIG" || parsedType=="TARG") &&
       (parsed_line[position+1].string_ == "VAL") )
  {
    if ( ((position+2) < numFields) && (parsed_line[position+2].string_ == "="))
    {
      position +=2;
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) <<
          "Invalid VAL= syntax in TRIG or TARG block on .MEASURE line";
      retval = false;
    }
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : extractMEASUREData
// Purpose       : Convert a .measure line to an options block
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool
extractMEASUREData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  // .measure has a lot of optonal and complex syntax.
  //
  // .measure mode name type output_var... options...
  //
  // where mode = AC, DC, FFT, NOISE or TRAN
  //       name = any text string for a name (can't be AC, DC or TRAN)
  //       type = keywords that let us figure out what type of measurement is needed
  //              TRIG and/or TARG = RiseFallDelay
  //              AVG, MAX, MIN, PP, RMS, INTEG, ON_TIME, OFF_TIME = Statistics
  //              FIND/WHEN = FindWhen
  //              EQN or PARAM = Equation evaluation
  //              DERIVATIVE or DERIV = Derivative
  //              INTEGRAL = Integral
  //              ERROR = Error
  //              ERR = Err function
  //              ERR1 = Err1 function
  //              ERR2 = Err2 function
  //              FOUR = Fourier analysis (similar to .FOUR)
  //       output_var = simulation variale to be measured.  This can be any of the following
  //              v(a), v(a,b), v(a)=number v(a)=v(b), i(a) ix(a) or an expression.
  //              So it really could be one or more
  //              variables and or a real number.  The comparison is always equity.
  //       options = these are keywords=value pairs that set limits on when or how
  //              the measurement is done. Value is usually a number, but can be a string in
  //              at least one case.
  //

  // we will use these sets to test for measure types and keywords.
  // These are separate sets for the measure TYPES for each allowed MODE.
  // typeSet is then the union of those sets.
  std::set<std::string> typeSet;
  std::set<std::string> typeSetTran;
  std::set<std::string> typeSetAc;
  std::set<std::string> typeSetDc;
  std::set<std::string> typeSetFft;
  std::set<std::string> typeSetNoise;
  std::set<std::string> typeSetCont;

  // allowed types for the TRAN mode
  typeSetTran.insert( std::string("TRIG") );
  typeSetTran.insert( std::string("TARG") );
  typeSetTran.insert( std::string("AVG") );
  typeSetTran.insert( std::string("MAX") );
  typeSetTran.insert( std::string("MIN") );
  typeSetTran.insert( std::string("PP") );
  typeSetTran.insert( std::string("RMS") );
  typeSetTran.insert( std::string("FREQ") );
  typeSetTran.insert( std::string("FIND") );
  typeSetTran.insert( std::string("WHEN") );
  typeSetTran.insert( std::string("PARAM") );
  typeSetTran.insert( std::string("EQN") );
  typeSetTran.insert( std::string("DERIVATIVE") );
  typeSetTran.insert( std::string("DERIV") );
  typeSetTran.insert( std::string("DUTY") );
  typeSetTran.insert( std::string("INTEGRAL") );
  typeSetTran.insert( std::string("INTEG") );
  typeSetTran.insert( std::string("ERROR") );
  typeSetTran.insert( std::string("ON_TIME") );
  typeSetTran.insert( std::string("OFF_TIME") );
  typeSetTran.insert( std::string("FOUR") );
  typeSetTran.insert( std::string("ERR") );
  typeSetTran.insert( std::string("ERR1") );
  typeSetTran.insert( std::string("ERR2") );

  // allowed types for the AC mode
  typeSetAc.insert( std::string("AVG") );
  typeSetAc.insert( std::string("DERIVATIVE") );
  typeSetAc.insert( std::string("DERIV") );
  typeSetAc.insert( std::string("ERROR") );
  typeSetAc.insert( std::string("EQN") );
  typeSetAc.insert( std::string("ERR") );
  typeSetAc.insert( std::string("ERR1") );
  typeSetAc.insert( std::string("ERR2") );
  typeSetAc.insert( std::string("FIND") );
  typeSetAc.insert( std::string("INTEGRAL") );
  typeSetAc.insert( std::string("INTEG") );
  typeSetAc.insert( std::string("MAX") );
  typeSetAc.insert( std::string("MIN") );
  typeSetAc.insert( std::string("PARAM") );
  typeSetAc.insert( std::string("PP") );
  typeSetAc.insert( std::string("RMS") );
  typeSetAc.insert( std::string("TRIG") );
  typeSetAc.insert( std::string("TARG") );
  typeSetAc.insert( std::string("WHEN") );

  // allowed types for the DC mode
  typeSetDc.insert( std::string("AVG") );
  typeSetDc.insert( std::string("DERIVATIVE") );
  typeSetDc.insert( std::string("DERIV") );
  typeSetDc.insert( std::string("ERROR") );
  typeSetDc.insert( std::string("EQN") );
  typeSetDc.insert( std::string("ERR") );
  typeSetDc.insert( std::string("ERR1") );
  typeSetDc.insert( std::string("ERR2") );
  typeSetDc.insert( std::string("FIND") );
  typeSetDc.insert( std::string("INTEGRAL") );
  typeSetDc.insert( std::string("INTEG") );
  typeSetDc.insert( std::string("MAX") );
  typeSetDc.insert( std::string("MIN") );
  typeSetDc.insert( std::string("PARAM") );
  typeSetDc.insert( std::string("PP") );
  typeSetDc.insert( std::string("RMS") );
  typeSetDc.insert( std::string("TRIG") );
  typeSetDc.insert( std::string("TARG") );
  typeSetDc.insert( std::string("WHEN") );

  // allowed types for the FFT mode
  typeSetFft.insert( std::string("ENOB") );
  typeSetFft.insert( std::string("EQN") );
  typeSetFft.insert( std::string("FIND") );
  typeSetFft.insert( std::string("PARAM") );
  typeSetFft.insert( std::string("SFDR") );
  typeSetFft.insert( std::string("SNDR") );
  typeSetFft.insert( std::string("SNR") );
  typeSetFft.insert( std::string("THD") );

  // allowed types for the NOISE mode
  typeSetNoise.insert( std::string("AVG") );
  typeSetNoise.insert( std::string("DERIVATIVE") );
  typeSetNoise.insert( std::string("DERIV") );
  typeSetNoise.insert( std::string("ERROR") );
  typeSetNoise.insert( std::string("EQN") );
  typeSetNoise.insert( std::string("ERR") );
  typeSetNoise.insert( std::string("ERR1") );
  typeSetNoise.insert( std::string("ERR2") );
  typeSetNoise.insert( std::string("FIND") );
  typeSetNoise.insert( std::string("INTEGRAL") );
  typeSetNoise.insert( std::string("INTEG") );
  typeSetNoise.insert( std::string("MAX") );
  typeSetNoise.insert( std::string("MIN") );
  typeSetNoise.insert( std::string("PARAM") );
  typeSetNoise.insert( std::string("PP") );
  typeSetNoise.insert( std::string("RMS") );
  typeSetNoise.insert( std::string("TRIG") );
  typeSetNoise.insert( std::string("TARG") );
  typeSetNoise.insert( std::string("WHEN") );

  // allowed types for the "continuous" measures
  typeSetCont.insert( std::string("DERIVATIVE") );
  typeSetCont.insert( std::string("DERIV") );
  typeSetCont.insert( std::string("FIND") );
  typeSetCont.insert( std::string("TRIG") );
  typeSetCont.insert( std::string("TARG") );
  typeSetCont.insert( std::string("WHEN") );

  // Make a union for the TYPE sets.  This is useful, once we know that the TYPE is valid
  // for one of the allowed modes (e.g, TRAN, AC, DC, FFT or NOSIE).  This happens after we've parsed
  // the MODE and TYPE fields on the measure line.  In those cases, we just care that the
  // TYPE is in the combined union set.
  set_union( typeSetTran.begin(), typeSetTran.end(), typeSetAc.begin(),
             typeSetAc.end(), std::inserter<std::set<std::string> >(typeSet, typeSet.begin()) );
  set_union( typeSetDc.begin(), typeSetDc.end(), typeSet.begin(),
             typeSet.end(), std::inserter<std::set<std::string> >(typeSet, typeSet.begin()) );
  set_union( typeSetFft.begin(), typeSetFft.end(), typeSet.begin(),
             typeSet.end(), std::inserter<std::set<std::string> >(typeSet, typeSet.begin()) );
  set_union( typeSetNoise.begin(), typeSetNoise.end(), typeSet.begin(),
             typeSet.end(), std::inserter<std::set<std::string> >(typeSet, typeSet.begin()) );
  set_union( typeSetCont.begin(), typeSetCont.end(), typeSet.begin(),
             typeSet.end(), std::inserter<std::set<std::string> >(typeSet, typeSet.begin()) );

  // Sets of allowed keywords.  simpleKeywords must have a numeric value.
  // numOrTextKeywords can have a numeric value or a string value (like LAST).
  // keywords is the union of those two sets.
  std::set<std::string> keywords;
  std::set<std::string> simpleKeywords;
  std::set<std::string> numOrTextKeywords;

  // these keywords will contain a number
  simpleKeywords.insert( std::string("TD") );
  simpleKeywords.insert( std::string("GOAL") );
  simpleKeywords.insert( std::string("WEIGHT") );
  simpleKeywords.insert( std::string("MINVAL") );
  simpleKeywords.insert( std::string("AT") );
  simpleKeywords.insert( std::string("FROM") );
  simpleKeywords.insert( std::string("TO") );
  simpleKeywords.insert( std::string("IGNORE") );
  simpleKeywords.insert( std::string("IGNOR") );
  simpleKeywords.insert( std::string("YMIN") );
  simpleKeywords.insert( std::string("YMAX") );
  simpleKeywords.insert( std::string("ON") );
  simpleKeywords.insert( std::string("OFF") );
  simpleKeywords.insert( std::string("FRAC_MAX") );
  simpleKeywords.insert( std::string("MIN_THRESH") );
  simpleKeywords.insert( std::string("MAX_THRESH") );
  simpleKeywords.insert( std::string("NUMFREQ") );
  simpleKeywords.insert( std::string("GRIDSIZE") );
  simpleKeywords.insert( std::string("DEFAULT_VAL") );
  simpleKeywords.insert( std::string("INDEPVARCOL") );
  simpleKeywords.insert( std::string("INDEPVAR2COL") );
  simpleKeywords.insert( std::string("DEPVARCOL") );
  simpleKeywords.insert( std::string("RFC_LEVEL") );
  simpleKeywords.insert( std::string("BINSIZ") );
  simpleKeywords.insert( std::string("MINFREQ") );
  simpleKeywords.insert( std::string("MAXFREQ") );
  simpleKeywords.insert( std::string("NBHARM") );

  // these keywords may contain string, or in some cases (like 
  // RISE, FALL and CROSS) may contain a string or a number
  numOrTextKeywords.insert( std::string("RISE") );
  numOrTextKeywords.insert( std::string("FALL") );
  numOrTextKeywords.insert( std::string("CROSS") );
  numOrTextKeywords.insert( std::string("FILE") );
  numOrTextKeywords.insert( std::string("COMP_FUNCTION") );
  numOrTextKeywords.insert( std::string("INDEPVARCOL") );
  numOrTextKeywords.insert( std::string("DEPVARCOL") );
  numOrTextKeywords.insert( std::string("PRECISION") );
  numOrTextKeywords.insert( std::string("PRINT") );
  numOrTextKeywords.insert( std::string("OUTPUT") );

  // make a union for the keywords set
  set_union( simpleKeywords.begin(), simpleKeywords.end(), numOrTextKeywords.begin(),
             numOrTextKeywords.end(), std::inserter<std::set<std::string> >(keywords, keywords.begin()) );

  // Note: option blocks with the name DOT_MEASURE_LINE are generated by .MEASURE lines.
  // option blocks with the name MEASURE come from .OPTIONS MEASURE lines.
  Util::OptionBlock option_block("DOT_MEASURE_LINE", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  const int numFields = parsed_line.size();
  std::string parsedType;  // used later to enable VAL= syntax for TRIG and TARG

  if( numFields < 4 )
  {
    Report::UserError().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Too few items on .MEASURE line.  Need at least .MEASURE <mode> <name> <type>";
    return false;
  }

  Util::Param parameter;

  // Look for the fixed items.
  // First word should be the MODE (e.g., TRAN, DC, AC or NOISE)
  ExtendedString currentWord(parsed_line[1].string_);
  currentWord.toUpper();
  if (currentWord == "TR")
  {
    currentWord = "TRAN"; // TR is a synonym for TRAN
  }
  // parsedMeasureMode is a local variable will be used to allow DC, AC and NOISE mode to only
  // be enabled for selected measure types.
  std::string parsedMode = currentWord;

  if( (currentWord == "TRAN" ) || (currentWord == "AC" ) || (currentWord == "DC") || (currentWord == "FFT") ||
      (currentWord == "NOISE") || (currentWord == "TRAN_CONT" ) || (currentWord == "AC_CONT" ) ||
      (currentWord == "DC_CONT") || (currentWord == "NOISE_CONT"))
  {
    parameter.set("MODE", std::string(currentWord));
    option_block.addParam(parameter);
  }
  else
  {
    Report::UserError0().at(netlist_filename, parsed_line[1].lineNumber_) << "Unknown mode in .MEASURE line.  "
	<< "Should be TRAN/TR, DC, AC, FFT, NOISE, TRAN_CONT, DC_CONT, AC_CONT or NOISE_CONT";
    return false;
  }

  // Second word should be the NAME of the measure.  The NAME can not be
  // one of the allowed modes. 
  currentWord = parsed_line[2].string_;
  currentWord.toUpper();
  if (currentWord == "TR")
  {
    currentWord = "TRAN"; // TR is a synonym for TRAN
  }
  if( currentWord != "DC" && currentWord != "TRAN" && currentWord != "AC" &&
      currentWord != "NOISE" && currentWord != "FFT" &&
      currentWord != "DC_CONT" && currentWord != "TRAN_CONT" &&
      currentWord != "AC_CONT" && currentWord != "NOISE_CONT")
  {
    parameter.set("NAME", std::string(currentWord));
    option_block.addParam(parameter);
  }
  else
  {
    Report::UserError0().at(netlist_filename, parsed_line[2].lineNumber_) << "Illegal name in .MEASURE line.  "
       << "Cannot be AC, DC, FFT, NOISE, TRAN/TR, AC_CONT, DC_CONT, NOISE_CONT or TRAN_CONT";
    return false;
  }

  // Third word should be the measure TYPE.  Different types may be allowed for
  // different modes.  So, these IF-ELSE blocks must use the mode-specific typeSet's.
  // Also, the error message strings may need to be updated as support for new/different
  // types are added for each mode.
  currentWord = parsed_line[3].string_;
  currentWord.toUpper();

  if (parsedMode == "TRAN")
  {
    if( typeSetTran.find( currentWord ) != typeSetTran.end() )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
      parsedType=currentWord;  // used later to enable the VAL= syntax for TRIG and TARG
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[3].lineNumber_) << "Illegal type in .MEASURE line for TRAN mode.  "
	   << "Must be one of: AVG, DERIV/DERIVATIVE, DUTY, EQN/PARAM, ERR, ERR1, ERR2, ERROR, FIND, WHEN, FOUR, "
           << "INTEG/INTEGRAL, MIN, MAX, OFF_TIME, ON_TIME, PP, RMS, TRIG, TARG";
      return false;
    }
  } 
  else if (parsedMode == "AC")
  {
    if(  typeSetAc.find( currentWord ) != typeSetAc.end()  )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
      parsedType=currentWord;  // used later to enable the VAL= syntax for TRIG and TARG
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[3].lineNumber_) << "Only AVG, DERIV, EQN/PARAM, ERR, ERR1, ERR2, "
	 << "ERROR, FIND, INTEG, MIN, MAX, PP, RMS, TRIG, TARG and WHEN measure types are supported for AC measure mode";
      return false;
    }
  }
  else if (parsedMode == "DC")
  {
    if(  typeSetDc.find( currentWord ) != typeSetDc.end()  )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
      parsedType=currentWord;  // used later to enable the VAL= syntax for TRIG and TARG
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[3].lineNumber_) << "Only AVG, DERIV, EQN/PARAM, ERR, ERR1, ERR2, "
	 << "ERROR, FIND, INTEG, MIN, MAX, PP, RMS, TRIG, TARG and WHEN measure types are supported for DC measure mode";
      return false;
    }
  }
  else if (parsedMode == "FFT")
  {
    if (typeSetFft.find( currentWord ) != typeSetFft.end()  )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
      parsedType=currentWord;  // used later to enable the VAL= syntax for TRIG and TARG
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[3].lineNumber_) << "Only ENOB, EQN/PARAM, FIND, SDFR, SNDR, SNR and THD "
	 << "measure types are supported for FFT measure mode";
      return false;
    }
  }
  else if (parsedMode == "NOISE")
  {
    if(  typeSetNoise.find( currentWord ) != typeSetNoise.end()  )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
      parsedType=currentWord;  // used later to enable the VAL= syntax for TRIG and TARG
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[3].lineNumber_) << "Only AVG, DERIV, EQN/PARAM, ERR, ERR1, ERR2, "
	 << "ERROR, FIND, INTEG, MIN, MAX, PP, RMS, TRIG, TARG and WHEN measure types are supported for NOISE measure mode";
      return false;
    }
  }
  else if ( (parsedMode == "TRAN_CONT") || (parsedMode == "AC_CONT") || (parsedMode == "DC_CONT") || (parsedMode == "NOISE_CONT"))
  {
    if( typeSetCont.find( currentWord ) != typeSetCont.end() )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
      parsedType=currentWord;  // used later to enable the VAL= syntax for TRIG and TARG
    }
    else
    {
      Report::UserError0().at(netlist_filename, parsed_line[3].lineNumber_) << "Only DERIV, FIND, TRIG, TARG and "
	 << "WHEN measure types are supported for continuous (CONT) measure modes";
      return false;
    }
  }

  // already got MEASURE <DC|AC|NOISE|TRAN> name TYPE
  // now try and parse of the rest of the line catching keywords when they're found
  int position = 4;

  // but first handle alternate PARAM syntaxes of PARAM=<exp> and PARAM <exp>
  // by skipping over the = token
  if ((parsedType =="PARAM") && (position < numFields) && (parsed_line[position].string_ ==  "=")) position++;

  while ( position < numFields )
  {
    currentWord = parsed_line[position].string_;
    currentWord.toUpper();

    // The type tag, TARG can get repeated, so check for it with other params.
    if( typeSet.find( currentWord ) != typeSet.end() )
    {
      parameter.set("TYPE", std::string(currentWord));
      option_block.addParam(parameter);
    }
    else if( numOrTextKeywords.find( currentWord ) != numOrTextKeywords.end() )
    {
      // these can be in the form of TAG=value or TAG=keyword where keyword = LAST
      // and potentially the "=" is optional.
      // IF and first ELSE-IF block needed to protect against core dumps from invalid .MEASURE lines.
      if ( (position+1) >= numFields )
      {
        // will core dump later if this is true, so error out now
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Invalid Measure Line";
        return false;
      } 
      else if ( (parsed_line[(position+1)].string_ == "=") && ( (position+2) >= numFields ) )
      {
        // will core dump later if this is true, so error out now
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Invalid Measure Line";
        return false;
      }
      else if( ((position+1) < numFields) || ((parsed_line[(position+1)].string_ == "=") && ((position+2) < numFields)) )
      {
        // value is in next or next + 1 position
        int valPosition = ++position;
        if( parsed_line[valPosition].string_ == "=" )
        {
          valPosition = ++position;
        }
        const std::string & value = parsed_line[valPosition].string_;
        if( Util::isInt(value) )
        {
          parameter.set(currentWord, Util::Ival(value) );
        }
        else if( Util::isValue(value) )
        {
          int valAsInt = static_cast<int>(Util::Value(value));
          parameter.set(currentWord, valAsInt);
        }
        else
        {
          parameter.set(currentWord, value);
        }
        option_block.addParam(parameter);
      }
      else
      {
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Incomplete .MEASURE line.  " 
	  << "RISE, FALL, CROSS must be followed by a value or LAST";
        return false;
      }

    }
    else if( simpleKeywords.find( currentWord ) != simpleKeywords.end() )
    {
      // these are in the form of TAG=value and potentially the "=" is optional.
      // IF and first ELSE-IF block needed to protect against core dumps from invalid lines.
      if ( (position+1) >= numFields )
      {
        // will core dump if this is true, so error out now
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Invalid Measure Line";
        return false;
      } 
      else if ( (parsed_line[(position+1)].string_ == "=") && ( (position+2) >= numFields) )
      {
        // will core dump if this is true
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Invalid Measure Line";
        return false;
      }
      else if( ((position+1) < numFields) || ((parsed_line[(position+1)].string_ == "=") && ((position+2) < numFields)) )
      {
        // value is in next or next + 1 position
        int valPosition = ++position;
        if( parsed_line[valPosition].string_ == "=" )
        {
          valPosition = ++position;
        }
        const std::string & value = parsed_line[valPosition].string_;
        if( (value[0]=='{') && (value[value.size()-1]=='}') )
        {
          // value is an expression, in the preferred Xyce syntax
          parameter.set(currentWord, Util::Expression(circuit_block.getExpressionGroup(),value));
        }
        else if ( (value == "PAR") && ((position + 3) < numFields ) && 
                  (parsed_line[position+1].string_ == "(") &&
                  (parsed_line[position+3].string_==")") )
        {
          // also handle HSPICE expression syntax of PAR({expString})
          std::string possExp = parsed_line[position+2].string_;
          if ( (possExp[0]=='{') && (possExp[possExp.size()-1]=='}') )
          {
            parameter.set(currentWord, Util::Expression(circuit_block.getExpressionGroup(),possExp));
            option_block.addParam(parameter);
            position += 3;
          }
        }
        else if ( (value== "(") && ((position + 2) < numFields) && 
                  (parsed_line[position+2].string_==")") )
        {
          // also handle HSPICE expression syntax of ({expString}) without the PAR
          std::string possExp = parsed_line[position+1].string_;
	  if ( (possExp[0]=='{') && (possExp[possExp.size()-1]=='}') )
	  {
            parameter.set(currentWord, Util::Expression(circuit_block.getExpressionGroup(),possExp));
            option_block.addParam(parameter);
            position += 2;
          } 
        }
        else if( Util::isValue(value) )
        {
          parameter.set(currentWord, Util::Value(value) );
        }
        else if( Util::isInt(value) )
        {
          parameter.set(currentWord, Util::Ival(value) );
        }
        else
        {
          Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Incomplete .MEASURE line.  "
				<< "TD, GOAL, WEIGHT, MINVAL, AT, TO, FROM, ON, OFF, IGNORE, YMIN, YMAX must be followed by a value";
          return false;
        }
        option_block.addParam(parameter);
      }
      else
      {
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Incomplete .MEASURE line.  "
	                        << "TD, GOAL, WEIGHT, MINVAL, AT, TO, FROM, ON, OFF, IGNORE, YMIN, YMAX must be followed by a value";
        return false;
      }
    }
    else
    {
      if( (currentWord[0]=='{') && (currentWord[currentWord.size()-1]=='}') )
      {
        // handle preferred Xyce syntax for expressions
        parameter.set(currentWord, std::string(currentWord));
        option_block.addParam(parameter);

	// skip over VAL= if this is a TRIG or TARG block 
        if (!handleValQualifier(netlist_filename, numFields, parsedType, parsed_line, position))
          return false;
      }
      else if ( currentWord== "PAR" ) 
      {
        // handle HSPICE expression syntax of PAR({expString})
        bool foundExp = false;
        if ( ((position + 3) < numFields ) && (parsed_line[position+1].string_ == "(") 
                && (parsed_line[position+3].string_==")") )
        {
          std::string possExp = parsed_line[position+2].string_;
          if ( (possExp[0]=='{') && (possExp[possExp.size()-1]=='}') )
          {
            parameter.set(possExp, std::string(possExp));
            option_block.addParam(parameter);
            position += 3;
            foundExp = true;
          } 
        }
        
        if (!foundExp)
	{
          Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Invalid PAR() expression syntax on .MEASURE line";
          return false;
        }

        // skip over VAL= if this is a TRIG or TARG block 
        if (!handleValQualifier(netlist_filename, numFields, parsedType, parsed_line, position))
          return false;
      }
      else if ( currentWord== "(" ) 
      {
        // handle HSPICE expression syntax of ({expString}) without the PAR
        bool foundExp = false;
        if ( ((position + 2) < numFields) && (parsed_line[position+2].string_==")") )
        {
          std::string possExp = parsed_line[position+1].string_;
	  if ( (possExp[0]=='{') && (possExp[possExp.size()-1]=='}') )
	  {
            parameter.set(possExp, std::string(possExp));
            option_block.addParam(parameter);
            position += 2;
            foundExp = true;
          } 
        }
        
        if (!foundExp)
	{
          Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Invalid () expression syntax on .MEASURE line";
          return false;
        }

        // skip over VAL= if this is a TRIG or TARG block
        if (!handleValQualifier(netlist_filename, numFields, parsedType, parsed_line, position))
          return false;
      }                      
      else
      {
      // The last form of legal syntax is out_var, out_var=val or out_var1=out_var2
      // where they are not within expression delimiters.  First, figure out how many 
      // positions between here and the next legal keyword or the end of the line.
      int endPosition=position;
      bool endNotFound=true;
      ExtendedString nextWord("");
      while( endNotFound )
      {
        endPosition++;
        if( endPosition == numFields)
        {
          endNotFound=false;
        }
        else
        {
          nextWord = parsed_line[endPosition].string_;
          nextWord.toUpper();
          if( (keywords.find( nextWord ) != keywords.end()) || (typeSet.find( nextWord ) != typeSet.end()) )
          {
            // found a keyword or type was found so we're at the end of space to
            // use in finding out_var, out_var, out_var=val or out_var1=out_var2
            endNotFound=false;
          }
        }
      }

      // ok, now package up the output vars.
      // It will be in the form ( [] indicate optional items )
      // V or I nodeName [nodeName] [ val | V or I nodeName [nodeName] ]
      // or { } expression delimited.  It might also use N, P or W.

      while( position < endPosition )
      {
        nextWord = parsed_line[position].string_;
        nextWord.toUpper();

        // check for ill-delimited expressions, without a matching '}' character.  
        // Any valid expression would have been processed above.
        if ( std::count(nextWord.begin(), nextWord.end(),'{') != std::count(nextWord.begin(), nextWord.end(),'}') )
	{
          Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Unmatched expression delimiter on .MEASURE line";
          return false;
	}

        // The second part of this if clause is to ensure we don't catch keywords that start
        // with I, V, N, P or W, and mistake them for Ixx( ) or Vxx().  Also need to not
        // mistake PAR for P(), and to allow for the HSPICE syntax of () around an expression.
        else if ( (position+1 < numFields && parsed_line[position+1].string_ == "(") &&
                  !((simpleKeywords.find(nextWord) != simpleKeywords.end()) || (nextWord == "PAR")
                    || (position+2 < numFields && parsed_line[position+2].string_[0]== '{')) )
	{
          // the error location and message returned by extractOperatorData() is currently not used here.
          int p_err=0;                           
          std::ostringstream msg;

	  if ( !IO::extractOperatorData(parsed_line, position, option_block, msg, p_err) )
          {
            Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Error in .MEASURE line.  "
	      << "Could not parse measured variable";
            return false;              
          }

          // extractOperatorData() advances to next token beyond the closing ).  Move back one, to
          // be compatible with the rest of the processing flow in this function.
          position -= 1;

          // now skip over VAL= if this is a TRIG or TARG block 
          if ( !handleValQualifier(netlist_filename, numFields, parsedType, parsed_line, position) )
            return false;
        }
        else if( nextWord.isValue() )
        {
          std::string objValue("OBJVAL");
          parameter.set( objValue, nextWord.Value());
          option_block.addParam( parameter );
        }
        else if( nextWord.isInt() )
        {
          std::string objValue( "OBJVAL");
          parameter.set( objValue, nextWord.Ival());
          option_block.addParam( parameter );
        }
        else if( (nextWord[0]=='{') && (nextWord[nextWord.size()-1]=='}') )
        {
          // value is an expression
          std::string objValue( "OBJVAL");
          parameter.set(objValue, Util::Expression(circuit_block.getExpressionGroup(),nextWord));
          option_block.addParam( parameter );
        }
        else if ( (nextWord == "PAR") && ((position + 3) < numFields ) && 
                  (parsed_line[position+1].string_ == "(") &&
                  (parsed_line[position+3].string_==")") )
        {
          // also handle HSPICE expression syntax of PAR({expString})
          std::string possExp = parsed_line[position+2].string_;
          if ( (possExp[0]=='{') && (possExp[possExp.size()-1]=='}') )
          {
            std::string objValue( "OBJVAL");
            parameter.set(objValue, Util::Expression(circuit_block.getExpressionGroup(),possExp));
            option_block.addParam(parameter);
            position += 3;
          }
        }
        else if ( (currentWord== "(") && ((position + 2) < numFields) && 
                  (parsed_line[position+2].string_==")") )
        {
          // also handle HSPICE expression syntax of ({expString}) without the PAR
          std::string possExp = parsed_line[position+1].string_;
	  if ( (possExp[0]=='{') && (possExp[possExp.size()-1]=='}') )
	  {
            std::string objValue( "OBJVAL");
            parameter.set(objValue, Util::Expression(circuit_block.getExpressionGroup(),possExp));
            option_block.addParam(parameter);
            position += 2;
          } 
        }
        else if( nextWord.possibleParam() )
        {
          // could be a param or measure name that will be evaluated
          // during the simulation.
          std::string objValue( "OBJVAL");
          parameter.set( objValue, std::string(nextWord));
          option_block.addParam( parameter );
        }

        position++;
      }
      // reset the position indicator to the end - 1 because
      // we're going to increment it at the end of the while loop
      position = endPosition - 1;
      }
    }
    position++;
  }

  circuit_block.addOptions(option_block);
  
  // success if we've reached here
  return true;
}

//-----------------------------------------------------------------------------
// Function      : populateMetadata
// Purpose       : 
// Special Notes : option blocks with the name MEASURE will come from .OPTIONS MEASURE
//                 lines.  option blocks from .MEASURE lines will have the name
//                 DOT_MEASURE_LINE.  The default values for MEASDGT and MEASFAIL
//                 match HSPICE, but are then overwritten by the values on the
//                 .OPTIONS MEASURE line.
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 09/28/2018
//-----------------------------------------------------------------------------
void populateMetadata(IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("MEASURE");

    parameters.insert(Util::ParamMap::value_type("MEASPRINT", Util::Param("MEASPRINT", "ALL")));
    parameters.insert(Util::ParamMap::value_type("MEASOUT", Util::Param("MEASOUT", 1)));
    parameters.insert(Util::ParamMap::value_type("MEASDGT", Util::Param("MEASDGT", 4)));
    parameters.insert(Util::ParamMap::value_type("MEASFAIL", Util::Param("MEASFAIL", 1)));
    parameters.insert(Util::ParamMap::value_type("DEFAULT_VAL", Util::Param("DEFAULT_VAL", -1)));
    parameters.insert(Util::ParamMap::value_type("USE_CONT_FILES", Util::Param("USE_CONT_FILES", 1)));
    parameters.insert(Util::ParamMap::value_type("USE_LTTM", Util::Param("USE_LTTM", 0)));
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : isComplexCurrentOp
// Purpose       : Determine if name is a "complex current operator" of the
//                 form IR, II, IM, IP or IDB
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/14/2021
//-----------------------------------------------------------------------------
bool isComplexCurrentOp(const std::string& name, int parenIdx)
{
  return ( ((parenIdx == 2) && ((name[1] == 'R') || (name[1] == 'I') || (name[1] == 'M') || (name[1] == 'P')))
           || ((parenIdx == 3) && (name.substr(1,2) == "DB")) );
}

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(Manager &manager, PkgOptionsMgr &options_manager)
{
  populateMetadata(options_manager);

  options_manager.addCommandParser(".MEAS", extractMEASUREData);
  options_manager.addCommandParser(".MEASURE", extractMEASUREData);

  // this handles .MEASURE lines
  options_manager.addCommandProcessor("DOT_MEASURE_LINE", new MeasureOptionsReg(manager));

  // this handles .OPTIONS MEASURE lines
  options_manager.addOptionsProcessor("MEASURE", 
      IO::createRegistrationOptions(manager, &Manager::registerMeasureOptions));

  return true;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

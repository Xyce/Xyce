//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose       : This file contains the functions to manage fourier objects
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_DEV_DeviceMgr.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_NetlistImportTool.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

// added for sensitivities:
#include <N_IO_Op.h>

#include <N_ERH_ErrorMgr.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_Op.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_SaveIOSState.h>

#include <Teuchos_ScalarTraits.hpp>

namespace Xyce {
namespace IO {

  // this is copied from the output manager:
namespace SensitivityOptions {

enum {
  DIRECT   = 0x01,
  ADJOINT  = 0x02,
  SCALED   = 0x04,
  UNSCALED = 0x08
};

}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::FourierMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
FourierMgr::FourierMgr(const CmdParse &cp)
  : commandLine_(cp),
    sensitivityOptions_(0),
    sensitivityRequested (false),
    numFreq_(10),
    gridSize_(200),
    calculated_(false)
{
  // Initialize first output variable pointer.  This vector is used to determine
  // which output variables are associated with which fundamental frequency.
  outputVarsPtr_.push_back( 0 );

  sensitivityOptions_ |= SensitivityOptions::ADJOINT;
  sensitivityOptions_ |= SensitivityOptions::UNSCALED;
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::FourierMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
FourierMgr::~FourierMgr()
{
  for (Util::Op::OpList::iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::notify
// Purpose       : Reset the Fourier analyses at the start of a STEP iteration,
//                 and output the Fourier results at the end of each STEP iteration.
//                 Output for the non-step case is currently handled by
//                 outputMacroResults().
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 6/01/2021
//-----------------------------------------------------------------------------
void FourierMgr::notify( const Analysis::StepEvent & step_event)
{
  switch (step_event.state_)
  {
    case Analysis::StepEvent::INITIALIZE:
      break;

    case Analysis::StepEvent::STEP_STARTED:
      reset();
      break;

    case Analysis::StepEvent::STEP_COMPLETED:
      outputResultsToFourFile(step_event.count_);
      break;

    case Analysis::StepEvent::FINISH:
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::reset
// Purpose       : Resets the object at the start of a .STEP loop
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 06/01/2021
//-----------------------------------------------------------------------------
void FourierMgr::reset()
{
  calculated_ = false;
  time_.clear();
  outputVarsValues_.clear();
  newTime_.clear();
  newValues_.clear();
  prdStart_.clear();
  lastPrdStart_.clear();

  mag_.clear();
  phase_.clear();
  nmag_.clear();
  nphase_.clear();
  thd_.clear();
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::fixupSensFourierParameters
//
// Purpose       : This function aguments the fixupFourerParameters function,
//                 but specifically deals with .FOUR for sensitivities.  The
//                 main complication is that (as of this writing) .FOUR for 
//                 sensitivities is specified as ".FOUR SENS" in the netlist, 
//                 and the actual sensitivity fields have to be inferred.
//
//                 So, the code below pulls the names of sensitivities out of
//                 the output manager (which is the only part of the code that 
//                 knows what they are) and then creates ops for them, following 
//                 the same patterns as used by the sensitivity outputters.
//
//                 Once it has done that, it aguments the existing, original data 
//                 structures such as  freqSolVarMap_, freqNamesMap_ 
//                 and freqNumOutputVarsMap_, so that the .FOUR machinery can be
//                 applied directly to the sensitivity outputs.
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 1/15/2019
//-----------------------------------------------------------------------------
void FourierMgr::fixupSensFourierParameters(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager) 
{
  if (sensitivityRequested)
  {


    // For each *unique* .FOUR SENS frequency, all the sensitivity outputs will be processed,
    for (std::map<double,int>::iterator it=sensFreqNumOutputVarsMap_.begin(); it!=sensFreqNumOutputVarsMap_.end(); ++it)
    {

      // find where each frequency begins and ends in the output Vars vector
      std::vector<int> tmpOutputVarsPtr(1,0);
      std::vector<double> tmpFreqVector;
      int ovPtrPos=0;
      for (std::map<double,int>::iterator it2=freqNumOutputVarsMap_.begin(); it2!=freqNumOutputVarsMap_.end(); ++it2, ++ovPtrPos)
      {
        tmpFreqVector.push_back(it2->first);
        tmpOutputVarsPtr.push_back( tmpOutputVarsPtr[ovPtrPos] + it2->second );
      }

      // temporary.  re-allocated each time thru the loop so each frequency gets unique sensitivity Ops
      Util::Op::OpList sensOutputVars_; 

      Util::Op::makeOps(comm, op_builder_manager, NetlistLocation(), 
          sensitivityVariableList_.begin(),
          sensitivityVariableList_.end(),
          std::back_inserter(sensOutputVars_));

      // copy the sensitivity information into the original .FOUR data structures
      double freq = it->first;
      it->second = sensOutputVars_.size();// don't know if I'll ever use this
      freqNumOutputVarsMap_[freq] += it->second;

      std::vector<double>::iterator foundFreq = std::find(tmpFreqVector.begin(), tmpFreqVector.end(), freq);

      if (foundFreq != tmpFreqVector.end()) 
      {
        int index = std::distance(tmpFreqVector.begin(), foundFreq); 
        int insertIndex = tmpOutputVarsPtr[index+1];
        outputVars_.insert( (outputVars_.begin() + insertIndex), sensOutputVars_.begin(), sensOutputVars_.end());
      }
      else
      {
        outputVars_.insert(outputVars_.end(), sensOutputVars_.begin(), sensOutputVars_.end());
      }

      for (Util::Op::OpList::iterator sensNamesIt = sensOutputVars_.begin(); sensNamesIt != sensOutputVars_.end(); ++sensNamesIt)
      {
        std::string name = (*sensNamesIt)->getName();
        freqNamesMap_.insert(std::make_pair(freq,name));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::fixupFourierParameters
// Purpose       : This function extracts info from the maps (freqSolVarMap_,
//                 freqNamesMap_ and freqNumOutputVarsMap_) created in 
//                 addFourierAnalysis(), in order to support multiple .FOUR
//                 lines that may have different fundamental frequencies. It
//                 creates the vectors (names_, freqVector_ and depSolVarIterVector_) 
//                 used by the existing code.  Simplifying this function likely 
//                 involves re-writing the rest of the .FOUR code.
//
// Special Notes : 
// Scope         : public
// Creator       : Not Sure (but likely Dave Baur)
// Creation Date : Not Sure
//-----------------------------------------------------------------------------
void FourierMgr::fixupFourierParameters(Parallel::Machine comm,
  const Util::Op::BuilderManager &op_builder_manager,
  const double endSimTime) 
{
  // this copy (from multi-map to vector) is done in order to support multiple .FOUR lines
  // The multi-map allows us to keep a sorted list of which solution variables are associated with
  // which fundamental frequency.
  for(std::multimap<double,Util::Param>::iterator it = freqSolVarMap_.begin(); it != freqSolVarMap_.end(); ++it ) 
  {
    // this vector is used by the makeOps function call below
    depSolVarIterVector_.push_back( it->second );   
  }

  if(!(depSolVarIterVector_.empty())) // no point in calling this if depSolVarIterVector is empty
  {
    Util::Op::makeOps(comm, op_builder_manager, NetlistLocation(), depSolVarIterVector_.begin(), depSolVarIterVector_.end(), std::back_inserter(outputVars_));
  }

  // this has to be called after the non-sensitivity ops are created (just above), 
  // but before the outputVarsPtr and names_, etc are set up.
  fixupSensFourierParameters(comm, op_builder_manager);

  // These copies (from set (or map) to vector) were done so that the uses of freqVector_ and names_
  // in the rest of the calculations and output functions did not have to be changed.  The 
  // freqNumOutputVarsMap_ contains the set of fundamental frequencies, and the number of output
  // variables at each frequency.  
  int ovPtrPos=0;
  for (std::map<double,int>::iterator it=freqNumOutputVarsMap_.begin(); it!=freqNumOutputVarsMap_.end(); ++it, ++ovPtrPos)
  {
    freqVector_.push_back(it->first);
    outputVarsPtr_.push_back( outputVarsPtr_[ovPtrPos] + it->second );
  }

  // check the AT values
  double startOfLastPeriod;
  for (int i=0; i<freqVector_.size(); ++i)
  {
    startOfLastPeriod = (freqVector_[i]*endSimTime - 1.0)/freqVector_[i];

    if ( !( (Teuchos::ScalarTraits<double>::magnitude(startOfLastPeriod) < Teuchos::ScalarTraits<double>::eps()) ||
            (startOfLastPeriod > 0) ) )
    {
      Xyce::Report::UserError0() << "The period (1/AT) requested on .FOUR line is greater than the length "
                                 << "of the time simulation for AT=" << freqVector_[i];
    }
  }

  // The freqNamesMap_ has a listing of "names" for the printResult_() function. 
  // It is needed in case the user intermixes .FOUR lines such as below. 
  //     .FOUR 1 V(1)
  //     .FOUR 2 V(2)
  //     .FOUR 1 V(3)
  for (std::multimap<double,std::string>::iterator it=freqNamesMap_.begin(); it!=freqNamesMap_.end(); ++it)
  {
    names_.push_back(it->second);
  } 

  // some error checking that the vectors have the right sizes and values
  if ( (names_.size() != outputVars_.size()) || ( outputVarsPtr_[outputVarsPtr_.size()-1] != outputVars_.size()) )
  {
    Report::UserFatal0() << "Unknown error while parsing .FOUR lines";
  }
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::getSensVars
// Purpose       : registers set of variables for .SENS.
//
// Special Notes : Copied from N_IO_OutputMgr::registerSens.  The main purpose
//                 is to set up the list of sensitivity variables.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/16/2019
//-----------------------------------------------------------------------------
bool FourierMgr::getSensVars (const Util::OptionBlock &option_block)
{
  bool bsuccess = true;

  std::vector<std::string> functions;
  std::vector<std::string> parameters;

  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    std::string tag = (*it).uTag();
    if ( std::string( (*it).uTag() ,0,7) == "OBJFUNC") // this is a vector
    {
      functions.push_back((*it).stringValue());
    }
    else if ( std::string( (*it).uTag() ,0,7) == "OBJVARS") // this is a vector
    {
      // do nothing. this is an option for AC
    }
    else if (std::string((*it).uTag(), 0, 5) == "PARAM")
    {
      parameters.push_back((*it).stringValue());
    }
    else if ( std::string( (*it).uTag() ,0,9) == "ACOBJFUNC") // this is a vector
    {
       // do nothing for now
    }
    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  int index = 0;
  int index2 = 0;
  for (std::vector<std::string>::const_iterator it1 = functions.begin(), end1 = functions.end(); it1 != end1; ++it1,++index2)
  {
    const std::string &function = (*it1);

    {
      Util::Marshal mout;
      mout << function << std::string("OBJECTIVEFUNCTION") << Util::Op::identifier<SensitivityObjFunctionOp>() << index2;
      sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
    }

    for (std::vector<std::string>::const_iterator it2 = parameters.begin(), end2 = parameters.end(); it2 != end2; ++it2, ++index)
    {
      const std::string &parameter = (*it2);

      if (sensitivityOptions_ & SensitivityOptions::DIRECT)
      {
        if (sensitivityOptions_ & SensitivityOptions::UNSCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpDirectOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
        if (sensitivityOptions_ & SensitivityOptions::SCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpDirectScaledOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
      }

      if (sensitivityOptions_ & SensitivityOptions::ADJOINT)
      {
        if (sensitivityOptions_ & SensitivityOptions::UNSCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpAdjointOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
        if (sensitivityOptions_ & SensitivityOptions::SCALED)
        {
          Util::Marshal mout;
          mout << function << parameter << Util::Op::identifier<SensitivitydOdpAdjointScaledOp>() << index;
          sensitivityVariableList_.push_back(Util::Param("SENS", mout.str()));
        }
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::registerSensOptions
// Purpose       : registers set of variables to set for .OPTIONS SENSITIVITY
//
// Special Notes : The only purpose of this is to set up the sensitivityOptions_ 
//                 variable, which is needed to know exactly which sensitivities 
//                 are being computed.
//
//                 Copied from N_IO_OutputMgr::registerSensOptions
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/16/2019
//-----------------------------------------------------------------------------
bool FourierMgr::registerSensOptions(const Util::OptionBlock &option_block)
{
  sensitivityOptions_ = 0;

  bool adjointGiven=false;
  bool outputUnscaledGiven=false;

  sensitivityOptions_ =0;

  Util::ParamList::const_iterator it = option_block.begin(); 
  Util::ParamList::const_iterator end = option_block.end();
  for ( ; it != end; ++it)
  {
    if ((*it).uTag() == "ADJOINT")
    {
      adjointGiven=true;
      if ((*it).getImmutableValue<bool>())
      {
        sensitivityOptions_ |= SensitivityOptions::ADJOINT;
      }
    }
    else if ((*it).uTag() == "DIRECT" && (*it).getImmutableValue<bool>())
    {
      sensitivityOptions_ |= SensitivityOptions::DIRECT;
    }
    else if ((*it).uTag() == "OUTPUTSCALED" && (*it).getImmutableValue<bool>())
    {
      sensitivityOptions_ |= SensitivityOptions::SCALED;
    }
    else if ((*it).uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledGiven=true;
      if ((*it).getImmutableValue<bool>())
      {
        sensitivityOptions_ |= SensitivityOptions::UNSCALED;
      }
    }
  }

  // default behavior is for the code to assume adjoint sensitivites are wanted
  // even if ADJOINT=1 wasn't specified.
  if (!adjointGiven)
  {
    sensitivityOptions_ |= SensitivityOptions::ADJOINT;
  }

  // default behavior is for the code to assume unscaled sensitivities are wanted
  // for output.
  if (!outputUnscaledGiven)
  {
    sensitivityOptions_ |= SensitivityOptions::UNSCALED;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::addFourierAnalysis
// Purpose       : Entry point when .four lines are pased in the netlist
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool FourierMgr::addFourierAnalysis(const Util::OptionBlock & fourierBlock)
{
  // based on what's in the option block passed in, we
  // create the needed fourier instance
  double freq;
  Util::ParamList variableList;   // Used to help register lead current requests with device manager.

  if (DEBUG_IO)
    Xyce::dout() << "In FourierMgr::addFourier" << std::endl
                 << ".FOUR line passed was: " << std::endl << fourierBlock << std::endl;

  const Util::Param *param = Util::findParameter(fourierBlock.begin(), fourierBlock.end(), "FREQ");
  if (!param)
  {
    // this shouldn't happen, but catch it if does
    Report::DevelFatal0() << "Missing FREQ in .FOUR line!";
  }

  for (Util::ParamList::const_iterator currentParamIt = fourierBlock.begin(); currentParamIt != const_cast<Util::OptionBlock &>(fourierBlock).end(); ++currentParamIt)
  {
    std::string name = "";
    std::string tag = currentParamIt->tag();

    if( tag == "FREQ" )
    {
      freq = currentParamIt->getImmutableValue<double>();
      if ( freq < 0 )
      {
        // a negative value for freq will cause a core dump later in the code.
        Report::DevelFatal0() << "Illegal FREQ value in .FOUR line!";
      }
      else
      {
        // a map is used here, in order to support multiple .FOUR lines with the same
        // fundamental frequency.  It stores a sorted version of the values of how many 
        // output variables are requested at each fundamental frequency
        freqNumOutputVarsMap_.insert(std::make_pair(freq,0));
      }
    }
    else if( (tag == "V") || (tag[0] == 'I') || (tag == "P") || (tag == "W") || (tag == "N") )
    {
      // tag[0] is used for I because branch currents for transistors can have two
      // characters.  An example is IS for the M-Device.
      int nodes = currentParamIt->getImmutableValue<int>();
      Util::Param aParam;
      aParam.set( tag, nodes );
      name += tag;
      name += "(";

      // add to list of variables.  This will be used later in netlist parsing 
      // to enable lead currents in the device manager.
      variableList.push_back(aParam);

      // here we just store the needed parts of V(a) or v(a,b) or I(device).
      // only the v(a,b) case will need an extra node in the outputVars_ array.
      // use list::insert() to keep an iterator pointing to this spot
    
      if ( freqNumOutputVarsMap_.size() == 0 )
      {
        // expected that a freq would have been found at this point, since the FREQ value 
        // comes before the ov variables on a .FOUR line.
        Report::DevelFatal0() << "Unable to parse .FOUR line";
      }
      else
      {
        // map is used here, in order to support multiple .FOUR lines, potentially with
        // different fundamental frequencies.
        freqSolVarMap_.insert(std::make_pair(freq,aParam));
        ++freqNumOutputVarsMap_[freq];
        for( int i=0; i<nodes; i++ )
        {
          currentParamIt++;
          aParam.set( currentParamIt->tag(), currentParamIt->getImmutableValue<double>() );
          freqSolVarMap_.insert(std::make_pair(freq,aParam));
          name += currentParamIt->tag(); 
          if (i != nodes-1 && nodes > 1) name += ",";     

          // add to list of variables.  This will be used later in netlist parsing 
          // to enable lead currents in the device manager.
          variableList.push_back(aParam);
        }
        name += ")";
      }
    }
    else if( Xyce::Util::hasExpressionTag(tag) )
    {
      if ( freqNumOutputVarsMap_.size() == 0 )
      {
        // expected that a freq would have been found at this point, since the FREQ value 
        // comes before the ov variables on a .FOUR line.
        Report::DevelFatal0() << "Unable to parse .FOUR line";
      }
      else
      {
        freqSolVarMap_.insert(std::make_pair(freq,*currentParamIt));
        ++freqNumOutputVarsMap_[freq];
        name = tag;

        // add to list of variables.  This will be used later in netlist parsing 
        // to enable lead currents in the device manager.
        Util::Param aParam;
        aParam.set( currentParamIt->tag(), currentParamIt->tag() );
        variableList.push_back(aParam);
      }
    }
    else if( (tag == "SENS") )
    {
      // don't know at this stage how many different sensitivities there are.    
      // So, sort it out later.
      sensitivityRequested = true;
      ++sensFreqNumOutputVarsMap_[freq];
    }

    // Save voltage or current variable name.
    if (name != "")
    {
      // a map is used here, in order to support multiple .FOUR lines, potentially with
      // different fundamental frequencies.
      freqNamesMap_.insert(std::make_pair(freq,name));
    }
  }
  
  // Used to help register lead current requests with device manager.
  getLeadCurrentDevices(variableList, devicesNeedingLeadCurrents_);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::updateFourierData
// Purpose       : Called during the simulation to update the fourier objects
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void FourierMgr::updateFourierData(Parallel::Machine comm, const double circuitTime, const Linear::Vector *solnVec, 
  const Linear::Vector *stateVec, const Linear::Vector * storeVec,
  const Linear::Vector *lead_current_vector, const Linear::Vector *junction_voltage_vector, 
  const Linear::Vector *lead_current_dqdt_vector,
  const std::vector<double> &         objectiveVec,
  const std::vector<double> &         dOdpVec,
  const std::vector<double> &         dOdpAdjVec,
  const std::vector<double> &         scaled_dOdpVec,
  const std::vector<double> &         scaled_dOdpAdjVec
  )
{
  // Save the time.
  if (outputVars_.size())
    time_.push_back(circuitTime);

  int vecIndex = 0;
  for (std::vector<Util::Op::Operator *>::const_iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
  {
    double retVal = getValue(comm, *(*it), Util::Op::OpData(vecIndex, solnVec, 0, stateVec, storeVec, 0, lead_current_vector, 0, junction_voltage_vector, 0, &objectiveVec, &dOdpVec, &scaled_dOdpVec, &dOdpAdjVec, &scaled_dOdpAdjVec)).real();

    outputVarsValues_.push_back(retVal);
    vecIndex++;
  }

}

//-----------------------------------------------------------------------------
// Function      : FFTMgr::outputResultsToFourFile
// Purpose       : Output all of the Fourier results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 06/01/2021
//-----------------------------------------------------------------------------
void FourierMgr::outputResultsToFourFile(int stepNumber)
{
  int numOutVars = outputVars_.size();

  if ( numOutVars && !time_.empty() && !calculated_ )
  {
    std::string filename = IO::makeOutputFileNameWithStepNum(commandLine_, ".four", stepNumber);
    std::ofstream outputFileStream;
    outputFileStream.open( filename.c_str() );
    outputResults(outputFileStream);
    outputFileStream.close();
  }
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::outputResults
// Purpose       : Output fourier results at end of simulation
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
void FourierMgr::outputResults(std::ostream& outputStream)
{
  // Only calculate something if a .four line was encountered and transient data was collected.
  int numOutVars = outputVars_.size();

  if ( numOutVars && !time_.empty() && !calculated_ )
  {
    // Calculated the fourier coefficients for the given nodes.
    getLastPeriod_();
    interpolateData_();
    calculateFT_();
    calculated_ = true;
  }

  // Output the information to the outputStream 
  printResult_( outputStream );

}


//-----------------------------------------------------------------------------
// Function      : FourierMgr::getLastPeriod_()
// Purpose       : finds the indices to access the last period of simulation
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void FourierMgr::getLastPeriod_()
{
  // We want to do the analysis on only the last period of the transient window. So here we find the indices
  // to access the endpoints of that interval.
  int numPoints = time_.size();
  int prdEnd = numPoints - 1;
  double endTime = time_[prdEnd];

  int nFreq = freqVector_.size();
  lastPrdStart_.resize( nFreq );
  prdStart_.resize( nFreq );

  for (int i=0; i<nFreq; ++i)
  {
    // Use this formulation to compute lastPrdStart to avoid issues with 32-bit arithmetic.
    lastPrdStart_[i] = (freqVector_[i]*endTime - 1.0)/freqVector_[i];

    if ( Teuchos::ScalarTraits<double>::magnitude( lastPrdStart_[i] ) < 
         Teuchos::ScalarTraits<double>::eps() )
    {
      lastPrdStart_[i] = 0.0;
      prdStart_[i] = 0;
    } 
    else if (lastPrdStart_[i] > 0)
    {
      // Initialize prdStart_ to be the index of the last element in time_.
      // Then scan until time_[i] <= endTime_ - period_.
      prdStart_[i] = prdEnd;
      while (time_[prdStart_[i]] > lastPrdStart_[i])
      {
        prdStart_[i] -= 1;
      }
    }
    else
    {
      std::string msg = "Error: The period is greater than the length of the time simulation. Exiting.";
      Report::UserFatal() << msg;
    }
  }
} 

//-----------------------------------------------------------------------------
// Function      : FourierMgr::interpolateData_()
// Purpose       : evaluates interpolating polynomial at equidistant time pts
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
bool FourierMgr::interpolateData_()
{
  double A, B, C;

  int numOutVars = outputVars_.size();
  int numPoints = time_.size();
  int nFreq = freqVector_.size();

  newTime_.resize(nFreq * gridSize_,0.0);
  newValues_.resize(numOutVars * gridSize_,0.0);

  for (int j=0; j<nFreq; ++j)
  {
    // Get the number of solution variables associated with this frequency
    int numCurrVars = outputVarsPtr_[j+1] - outputVarsPtr_[j];

    // Get the number of time points in the final period
    int nData = numPoints-prdStart_[j];
    double period = 1.0/freqVector_[j];
    std::vector<double> h(nData-1, 0.0);
    std::vector<double> b(nData-1, 0.0);
    std::vector<double> u(nData-1, 0.0);
    std::vector<double> v(nData-1, 0.0);
    std::vector<double> z(nData, 0.0);

    // Compute new, equally spaced time points.
    double step = period/gridSize_;

    newTime_[gridSize_*j] = lastPrdStart_[j];
    for (int i = 1; i < gridSize_; i++) 
    {
      newTime_[gridSize_*j + i] = newTime_[gridSize_*j + i-1] + step;
    }
  
    // Loop over all the data from the output variables associated with this frequency 
    for (int k=0; k<numCurrVars; k++)
    { 
      int offset = outputVarsPtr_[j] + k;

      // Cubic spline interpolation. We first need to find the z's.
      for (int i = 0; i < nData-1; i++)
      {
        h[i] = time_[i+1+prdStart_[j]]-time_[i+prdStart_[j]];
        b[i] = (6/h[i])*(outputVarsValues_[(i+1+prdStart_[j])*numOutVars + offset]
                         - outputVarsValues_[(i+prdStart_[j])*numOutVars + offset]);
      }

      u[1] = 2*(h[0]+h[1]);
      v[1] = b[1]-b[0];

      for (int i=2; i < nData-1; i++)
      {
        u[i] = 2*(h[i]+h[i-1])-((h[i-1])*(h[i-1]))/u[i-1];
        v[i] = b[i]-b[i-1]-(h[i-1]*v[i-1])/u[i-1];
      }

      z[nData-1] = 0;
      for (int i=nData-2; i > 0; i--)
      {
        z[i] = (v[i]-h[i]*z[i+1])/u[i];
      }
      z[0] = 0;

      // Calculate the new values at the new time points.
      for (int i = 0; i < gridSize_; i++)
      { 
        int idx = nData-1;
        while ((newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]]) < 0)
        {
          idx--;
        }
        A = (z[idx+1]-z[idx])/(6*h[idx]);
        B = z[idx]/2;
        C = -(h[idx]/6)*z[idx+1]-(h[idx]/3)*z[idx]+(outputVarsValues_[(idx+1+prdStart_[j])*numOutVars + offset]
                                                    - outputVarsValues_[(idx+prdStart_[j])*numOutVars + offset])/h[idx];
        newValues_[gridSize_*offset + i] = outputVarsValues_[(idx+prdStart_[j])*numOutVars + offset] + (newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]])*(C+(newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]])*(B+(newTime_[gridSize_*j + i]-time_[idx+prdStart_[j]])*A));
      }
    } 
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::calculateFT_()
// Purpose       : performs fourier analysis
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void FourierMgr::calculateFT_()
{ 
  int numOutVars = outputVars_.size();
  int nFreq = freqVector_.size();

  mag_.resize(numFreq_*numOutVars, 0.0);
  phase_.resize(numFreq_*numOutVars, 0.0);
  nmag_.resize(numFreq_*numOutVars, 0.0);
  nphase_.resize(numFreq_*numOutVars, 0.0);
  freq_.resize(numFreq_*nFreq, 0.0);
  thd_.resize(numOutVars, 0.0);

  // Compute frequencies 
  for (int j=0; j<nFreq; ++j)
  {
    for (int i=0; i < numFreq_; i++)
    {
      freq_[j*numFreq_ + i] = i * freqVector_[j];
    }
  }
   
  // Perform Fourier analysis for all the output variables. 
  for (int k=0; k < numOutVars; k++)
  { 
    for (int i=0; i < gridSize_; i++)
    { 
      for (int j=0; j < numFreq_; j++)
      {
        mag_[numFreq_*k + j] += (newValues_[gridSize_*k+i]*sin(j*2.0*M_PI*i/((double) gridSize_)));
        phase_[numFreq_*k + j] += (newValues_[gridSize_*k+i]*cos(j*2.0*M_PI*i/((double) gridSize_)));
      }
    }
  
    mag_[numFreq_*k] = phase_[numFreq_*k]/gridSize_;
    phase_[numFreq_*k] = 0;
    thd_[k] = 0; 

    double convRadDeg = 180.0/M_PI;

    for(int i = 1; i < numFreq_ ; i++)
    { 
      double tmp = mag_[numFreq_*k+i]*2.0 /gridSize_;
      phase_[numFreq_*k+i] *= 2.0/gridSize_;
      mag_[numFreq_*k+i] = sqrt(tmp*tmp+phase_[numFreq_*k+i]*phase_[numFreq_*k+i]);
      phase_[numFreq_*k+i] = atan2(phase_[numFreq_*k+i],tmp)*convRadDeg;
      nmag_[numFreq_*k+i] = mag_[numFreq_*k+i]/mag_[numFreq_*k+1];
      nphase_[numFreq_*k+i] = phase_[numFreq_*k+i]-phase_[numFreq_*k+1];
      if(i>1) thd_[k] += nmag_[numFreq_*k+i]*nmag_[numFreq_*k+i];
    }
    thd_[k] = 100*sqrt(thd_[k]);
  }
}

//-----------------------------------------------------------------------------
// Function      : FourierMgr::printResult_( std::ostream& os )
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
std::ostream& FourierMgr::printResult_( std::ostream& os )
{
  basic_ios_all_saver<std::ostream::char_type> save(os);

  // Compute measure.
  if (calculated_)
  { 
    for (unsigned int i=0; i<freqVector_.size(); ++i)
    {
      for (int j=outputVarsPtr_[i]; j<outputVarsPtr_[i+1]; j++)
      {
        int colWidth1=12, colWidth2 = 16, precision = 6;
        os << "Fourier analysis for " << names_[j] << ":" << std::endl;
        os << "  No. Harmonics: " << numFreq_ << ", THD: " 
           << std::scientific << std::setprecision(precision)
           << thd_[j] << ", Gridsize: " << gridSize_
           << ", Interpolation Type: Cubic Spline" << std::endl;
  
        os << std::setw(colWidth1) << "Harmonic" << std::setw(colWidth2) << "Frequency"
           << std::setw(colWidth2) << "Magnitude" << std::setw(colWidth2) << "Phase"
           << std::setw(colWidth2) << "Norm. Mag" << std::setw(colWidth2) << "Norm. Phase" << std::endl;
        for (int k = 0; k < numFreq_; k++)
        {
           os << std::setw(colWidth1) << k << std::setw(colWidth2) << freq_[numFreq_*i+k]
              << std::setw(colWidth2) << mag_[numFreq_*j+k]
              << std::setw(colWidth2) << phase_[numFreq_*j+k]
              << std::setw(colWidth2) << nmag_[numFreq_*j+k]
              << std::setw(colWidth2) << nphase_[numFreq_*j+k] << std::endl;
        }
        os << std::endl;
      }  
    }
  }
  return os;
}


namespace {

//-----------------------------------------------------------------------------
// Function      : extractFOURIERData
// Purpose       : Convert a .four line to an options block
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 06/03/2013
//-----------------------------------------------------------------------------
bool extractFOURIERData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("FOUR", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  Util::Param parameter;
  ExtendedString nextWord("");

  if(parsed_line.size() < 3)
  {
    Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
      << "Error: .FOUR line needs at least 3 arguments '.FOUR freq ov1 <ov2 ...>'";
    return false;
  }

  parameter.setTag("FREQ");
  parameter.setVal(parsed_line[1].string_ );
  option_block.addParam( parameter );

  int position = 2;
  int endPosition = parsed_line.size();
  while (position < endPosition)
  {

    nextWord = parsed_line[position].string_;
    nextWord.toUpper();
    // nextWord[0] is used for I because branch currents for transistors can have two
    // characters.  An example is IS for the M-Device.
    if( (nextWord[0] == 'I') || (nextWord == "V") || (nextWord == "P") || (nextWord == "W") || (nextWord == "N") )
    {
      // need to do a bit of look ahead here to see if this is a V(a) or V(a,b)
      int numNodes = 1;
      if( ((position+3) < endPosition) && (parsed_line[position+3].string_ == ",") )
      {
        numNodes = 2;
      }
      parameter.set(nextWord, numNodes);
      option_block.addParam( parameter );

      if( (position+2) < endPosition )
      {
        position+=2;
        nextWord = parsed_line[position].string_;
        nextWord.toUpper();
        parameter.set( nextWord, 0.0 );
        option_block.addParam( parameter );

        if( ((position+2) < endPosition) && (parsed_line[position+1].string_ == ",") )
        {
          // may have a voltage difference request
          position+=2;
          nextWord = parsed_line[position].string_;
          nextWord.toUpper();
          parameter.set( nextWord, 0.0 );
          option_block.addParam( parameter );
        }
      }
      else
      {
        Report::UserError0().at(netlist_filename, parsed_line[position].lineNumber_) << "Could not parse .FOUR variable";
      }
      position+=2;
    }
    else if ( (nextWord[0]=='{') && (nextWord[nextWord.size()-1]=='}') )
    {
      parameter.set(nextWord, std::string(nextWord));
      option_block.addParam(parameter);
      position++;
    }
    else if ( (nextWord =="SENS") )
    {
      parameter.set(nextWord, std::string("ALL"));
      option_block.addParam(parameter);
      position++;
    }
    else
    {
      Report::UserError().at(netlist_filename, parsed_line[position].lineNumber_) << "Could not parse .FOUR variable " << nextWord;
      position+=2;
    }
  }

  circuit_block.addOptions(option_block);

  return true;
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(FourierMgr &fourier_manager, PkgOptionsMgr &options_manager)
{
  options_manager.addCommandParser(".FOUR", extractFOURIERData);

  options_manager.addCommandProcessor("FOUR", IO::createRegistrationOptions(fourier_manager, &FourierMgr::addFourierAnalysis));
  options_manager.addCommandProcessor("SENS", IO::createRegistrationOptions(fourier_manager, &FourierMgr::getSensVars));

  options_manager.addOptionsProcessor("SENSITIVITY", IO::createRegistrationOptions(fourier_manager, &FourierMgr::registerSensOptions));

  return true;
}

} // namespace IO
} // namespace Xyce

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
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_MeasureError.h>
#include <N_IO_OutputPrn.h>
#include <N_IO_OutputCsd.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Interpolators.h>
#include <Teuchos_SerialDenseVector.hpp>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Error::Error()
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
Error::Error(const Manager &measureMgr, const Util::OptionBlock & measureBlock ):
  Base(measureMgr, measureBlock)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() and update AC() are likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();

  // Explicitly check that both column indexes are >= 0, since the default values for 
  // independentVarColumn_ and dependentVarColumn_ were set to -1 in the base object
  // constructor.  This handles incorrect values from the user.  For TRAN and AC
  // measures, the error conditions are:
  // 
  //   a) either the INDEPVARCOL or DEPVARCOL qualifiers were not specified on the 
  //      .MEASURE instance line.  
  //
  //   b) the INDEPVARCOL or DEPVARCOL qualifiers have the same values. 
  //
  //   c) The values in the INDEPVARCOL are not monotonically increasing.
  //  
  // For DC measures, only the DEPVARCOL qualifier must be specified on the .MEASURE
  // line, and be valid.  The INDEPVARCOL qualifer is silently ignored if it is specifed.  
  //
  // Finally, the case of either column index referring to a non-existent column in the 
  // data file is also handled below, after the dataFile is read in.
  if (mode_ == "DC")
  {
    if ( dependentVarColumn_ < 0 )
    {
      Report::UserError0() << "In measure " << name_ << ", missing or negative value for DEPVARCVOL";
      return;
    }
  }
  else
  { 
    // check TRAN and AC measures
    if ( (independentVarColumn_ < 0) || (dependentVarColumn_ < 0) )
    {
      Report::UserError0() << "In measure " << name_ << ", missing or negative value for INDEPVARCOL or DEPVARCVOL";
      return;
    }
    else if ( independentVarColumn_ == dependentVarColumn_ )
    {
      Report::UserError0() << "In measure " << name_ << ", identical values for INDEPVARCOL and DEPVARCVOL";
      return;
    }
  }

  // Some error checking on the comparison file name specified by FILE=  
  // Will check for a valid file extension (e.g., .PRN) below.
  if ( !(Util::checkIfValidFile(dataFileName_)) )
  {
    // Error out if the user-specified comparison file does not exist, cannot be opened,
    // or is a directory name rather than a file name.  See SON Bugs 730 
    // and 785 for more details.
    Report::UserError0() << "In measure " << name_ << ", could not find comparison file";
    return;
  }
  else if ( dataFileName_.length() < 4 )
  {
    // the ExtendedString fileExt() line below gives a "cryptic error message" without this check
    Report::UserError0() << "In measure " << name_ 
        << ", comparison filename must end in .PRN, .CSV or .CSD";
    return;
  }

  // Use the same functions for reading the .PRN, .CSD and .CSV files as re-measure.
  // The function logic is to first test for a valid extension. If the extension is invalid 
  // then the function errors out.  If the extension is valid, then make the correct object 
  // type and try to open the file.
  OutputFileBase * comparisonFile;
  ExtendedString fileExt(dataFileName_.substr(dataFileName_.length()-4));

  // The Xyce format for a .CSV file is basically the same as a .PRN file made with the
  // -delim COMMA command line option.  The only difference is that the .PRN file has something
  // like "End of Xyce(TM) Simulation" as its last line.  So, reading-in of .PRN and .CSV files 
  // currently use the same function.
  if ( (fileExt.toUpper() == ".PRN")  || (fileExt.toUpper() == ".CSV") )
  {
    comparisonFile = new OutputPrn();
  }
  else if (fileExt.toUpper() == ".CSD")
  {
    comparisonFile = new OutputCsd();
  }
  else
  {
    // unsupported file format.  Report error and exit.
    Report::UserError0() << "For measure " + name_ + ", ERROR measure only supports .PRN, .CSV or .CSD formats";
    return;
  }

  if (!comparisonFile->openFileForRead(dataFileName_))
  {
    // open failed.  Report error and exit.
    Report::UserError0() << "Could not open file in ERROR measure " + name_ + ". Failed filename was = " + dataFileName_ ;
    return;
  }

    // load data-names
  std::vector<std::string> fileVarNames;
  if (!(comparisonFile->getOutputVarNames(fileVarNames)))
  {
    // reading var names failed.  report error and exit remeasure routine.
    Report::UserError0() << "Problem reading variable names in file for measure " + name_ ;
    return;
  }

  // used for error measure generation below, during reading-in of comparison file
  double minIndepVarValues = 0.0;
  int idx=0;

  // This code section run though the lines in the comparison file, and makes the 
  // the vectors for the independent (e.g., time or frequency) and dependent values.
  int reading=1; 
  std::vector<double> varValuesVec;

  while (reading==1)
  {
    reading = comparisonFile->getOutputNextVarValuesSerial(&varValuesVec);

    if( reading == 1 )
    {
      int varValuesVecSize = varValuesVec.size();
      if (mode_ == "DC")
      {
        // DC measures only use the DEPVARCOL qualifier
        if( dependentVarColumn_ >=  varValuesVecSize )
        {
          Report::UserError0() << "In measure " << name_  << ", using data from file " << dataFileName_ 
                               << ". Requested column for dependent variable " 
                               << dependentVarColumn_ << " does not exist in the data file for entry " 
                               << idx;

          // exit after the first line with this error.
          return;
        }  
        else
        {
          dataValues_.push_back(varValuesVec[dependentVarColumn_]);
        }
      }
      else
      {
        // AC and TRAN measures use both the INDEPVARCOL or DEPVARCOL qualifiers.
        // Throw an error if either of the requested columns do not exist in the data file.
        // This will only report once such error, for each ERROR measure instance line, before exiting.
        int maxVarColumn = ( independentVarColumn_ > dependentVarColumn_ ) ? independentVarColumn_ : dependentVarColumn_;
        if( maxVarColumn >=  varValuesVecSize )
        {
          Report::UserError0() << "In measure " << name_  << ", using data from file " << dataFileName_ 
                               << ". Requested column for independent variable " 
                               << independentVarColumn_ << " or dependent variable " 
                               << dependentVarColumn_ << " does not exist in the data file for entry " 
                               << idx;

          // exit after the first line with this error.
          return;
        }  
        else
        {
          indepVarValues_.push_back(varValuesVec[independentVarColumn_]); 
          dataValues_.push_back(varValuesVec[dependentVarColumn_]);

          // check that time (or frequency) is monotonically increasing in the comparison file, 
          // for TRAN (or AC) mode measures, since interpolation is done for those cases.
          if ( indepVarValues_[idx] < minIndepVarValues )
          {
            Report::UserError0() << "In measure " << name_ << ", using data from file " << dataFileName_ 
                               << ". Independent variables are not sorted in monotonically increasing order at entry " 
                               << idx;
          }
          minIndepVarValues = indepVarValues_[idx];
          idx++;
        }
      } 

      varValuesVec.clear();
    }
  }

  // results from the simulation to compare to a column in dataValues.
  // simulationDataVals_ is the dependent variable.  simulationIndepVarVals_ is
  // then the independent variable, which is time for TRAN mode and frequency for 
  // AC mode.  simulationIndepVarVals_ is not used for DC mode.
  simulationDataVals_.reserve( dataValues_.size() );
  if (mode_ != "DC") simulationIndepVarVals_.reserve( dataValues_.size() );

  // clean up file object, which is local to this constructor
  comparisonFile->closeFileForRead();
  delete comparisonFile;
  comparisonFile = 0;
}

//-----------------------------------------------------------------------------
// Function      : Error::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void Error::prepareOutputVariables()
{
  // this measurement should have only one dependent variable.
  // Error out if it doesn't
  numOutVars_ = outputVars_.size();
  
  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for ERROR measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }

  outVarValues_.resize( numOutVars_, 0.0 );
}

//-----------------------------------------------------------------------------
// Function      : Error::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Error::reset() 
{
  resetBase();

  // need to reset both vectors used in updateTran() and updateAC()
  simulationDataVals_.clear();
  simulationIndepVarVals_.clear();
}


//-----------------------------------------------------------------------------
// Function      : Error::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void Error::updateTran(
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
  if( !calculationDone_ )
  {
    simulationIndepVarVals_.push_back( circuitTime );

    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i],
                                        solnVec, stateVec, storeVec, 0,
                                        lead_current_vector,
                                        junction_voltage_vector,
                                        lead_current_dqdt_vector, 0, 0, 0, 0);
      simulationDataVals_.push_back( outVarValues_[i] );
    }
    initialized_ = true;
  }
}


//-----------------------------------------------------------------------------
// Function      : Error::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
void Error::updateDC(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  // The dcParamsVec will be empty if the netlist has a .OP statement without a .DC statement.
  // In that case, a DC MEASURE will be reported as FAILED.
  if( !calculationDone_ && ( dcParamsVec.size() > 0 ) )
  {   
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i],
                                        solnVec, stateVec, storeVec, 0,
                                        lead_current_vector,
                                        junction_voltage_vector,
                                        lead_current_dqdt_vector, 0, 0, 0, 0);
      simulationDataVals_.push_back( outVarValues_[i] );
    }
    initialized_ = true;
    sweepVar_= getDCSweepVarName(dcParamsVec); // used in descriptive output to stdout
  }
 
}

//-----------------------------------------------------------------------------
// Function      : Error::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2014
//-----------------------------------------------------------------------------
void Error::updateAC(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{
  if( !calculationDone_ )
  {
    simulationIndepVarVals_.push_back( frequency );

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, 0, 0,
                                        imaginaryVec, 0, 0, 0, 0, 0, 0, RFparams );
      simulationDataVals_.push_back( outVarValues_[i] );
    }
    initialized_ = true;
  }
}

void Error::updateNoise(
  Parallel::Machine comm,
  double frequency,
  double fStart,
  double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  double totalOutputNoiseDens,
  double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec)
{
  if( !calculationDone_ )
  {
    simulationIndepVarVals_.push_back( frequency );

    // update our outVarValues_ vector
    for( int i=0; i< numOutVars_; i++ )
    {
      outVarValues_[i] = getOutputValue(comm, outputVars_[i], solnVec, 0, 0,
                                        imaginaryVec, 0, 0, 0,
                                        totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, 0);
      simulationDataVals_.push_back( outVarValues_[i] );
    }
    initialized_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : Error::printMeasureResult()
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes : ERROR measure needs its own printMeasureResult() function
//                 because it uses parallel comms.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 8/03/2016
//-----------------------------------------------------------------------------
std::ostream& Error::printMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if ( !initialized_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
    {
      // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 is given in the
      // netlist and this is a failed measure.
      os << name_ << " = FAILED" << std::endl;
    }
    else 
    {
      os << name_ << " = " << this->getMeasureResult() << std::endl;
    }

    return os;
}

//-----------------------------------------------------------------------------
// Function      : Error::printVerboseMeasureResult()
// Purpose       : used to print the "verbose" (more descriptive) measurement
//                 result to an output stream object, which is typically stdout
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 2/22/2015
//-----------------------------------------------------------------------------
std::ostream& Error::printVerboseMeasureResult(std::ostream& os)
{
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);

    if (initialized_)
    {
      os << name_ << " = " << this->getMeasureResult() << std::endl ;
    }
    else
    { 
      os << name_ << " = FAILED" << std::endl;
    }

    // also print out the COMP_FUNCTION used, and whether the default value was used
    if ( usedDefaultComparisonFunc_ && (comparisonFunctionName_ == "L2NORM") )
    {
      // explicitly state if the default value was used, because invalid values
      // default to the L2NORM.  See SON Bug 839.
      os << "COMP_FUNCTION defaulted to L2NORM" << std::endl; 
    }
    else
    {
      os << "COMP_FUNCTION was " << comparisonFunctionName_ << std::endl; 
    }     

    return os;
}

//-----------------------------------------------------------------------------
// Function      : Error::printMeasureWindow
// Purpose       : prints information related to time/frequency window used.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/8/2016
//-----------------------------------------------------------------------------
std::ostream& Error::printMeasureWindow(std::ostream& os, double endSimTime,
				        double startSweepVal, double endSweepVal) const
{
  // The measure window info is not printed for DC mode.
  if (mode_ != "DC")
  {
    basic_ios_all_saver<std::ostream::char_type> save(os);
    os << std::scientific << std::setprecision(precision_);
    
    // modeStr is "Time" for TRAN mode and "Freq" for AC mode.
    std::string modeStr = setModeStringForMeasureWindowText();

    // the "measurement window" for the ERROR measure is set equal to the 
    // begin and end times (or frequencies) of the comparison waveform.
    os << "Measure Start " << modeStr << "= " << indepVarValues_.front() 
       << "\tMeasure End " << modeStr << "= " << indepVarValues_.back() << std::endl; 
  }
  return os;
}

//-----------------------------------------------------------------------------
// Function      : Error::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
double Error::getMeasureResult()
{
  // Only need to execute these calculations on the first pass through this 
  // function.  This function will be called twice (once for the output to 
  // the .mtX file and a second time for the output to stdout). This can improve
  // performance if the comparision file has many data points in it.
  if (!gotMeasureResult_ && initialized_)
  {
    // make a vector to hold difference values
    Teuchos::SerialDenseVector<int,double> differenceVector(dataValues_.size());

    if (mode_ == "DC")
    {
      // No interpolation for DC mode measures
      for (unsigned int i=0 ; i<dataValues_.size(); i++)
      {
        differenceVector[i] = (simulationDataVals_[i]-dataValues_[i]);
      }
    }
    else
    {
      // Interpolate simulation data to input data values (in the comparison file)
      // for TRAN or AC mode measures.
      std::vector<double> interpDataVals( dataValues_.size(), 0.0 );
      if (simulationIndepVarVals_.size() > 0)
      {
        Util::akima<double> interp;
        interp.init( simulationIndepVarVals_, simulationDataVals_ );
        for (unsigned int i=0; i < dataValues_.size(); i++)
        {
          interp.eval( simulationIndepVarVals_, simulationDataVals_, indepVarValues_[i], interpDataVals[i] );
        }
      }
        
      // load the difference vector
      for (unsigned int i=0 ; i<dataValues_.size(); i++)
      {
        differenceVector[i] = (interpDataVals[i]-dataValues_[i]);
      }
    }

    // default to L2NORM for any string that is not L1NORM or INFNORM
    if( comparisonFunctionName_ == "L1NORM" )
    {
      calculationResult_ = differenceVector.normOne();
    }
    else if( comparisonFunctionName_ == "INFNORM" )
    {
      calculationResult_ = differenceVector.normInf();
    }
    else
    {
      calculationResult_ = differenceVector.normFrobenius();
    }
  }
  gotMeasureResult_=true;

  return calculationResult_;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

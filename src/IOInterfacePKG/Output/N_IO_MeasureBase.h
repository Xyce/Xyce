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
//
// Purpose        : Base class for Measure types
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 03/10/2009
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_MeasureBase_h
#define Xyce_N_IO_MeasureBase_h

#include <string>
#include <list>
#include <vector>

#include <N_ANP_SweepParam.h>
#include <N_IO_fwd.h>
#include <N_IO_MeasureManager.h>
#include <N_LAS_Vector.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Param.h>
#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-------------------------------------------------------------------------
// Class         : MeasureBase
// Purpose       : Base class for common analysis functions
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Base
{
  public:
    Base(const Manager &measureMgr, const Util::OptionBlock & measureBlock);

    virtual ~Base();

    virtual void prepareOutputVariables() = 0;
    virtual void reset() {initialized_ = false;}

    virtual void updateTran(
              Parallel::Machine comm,
              const double circuitTime,
              const Linear::Vector *solnVec,
              const Linear::Vector *stateVec,
              const Linear::Vector *storeVec,
              const Linear::Vector *lead_current_vector,
              const Linear::Vector *junction_voltage_vector,
              const Linear::Vector *lead_current_dqdt_vector) {}

    virtual void updateDC(
              Parallel::Machine comm,
              const std::vector<Analysis::SweepParam> & dcParamsVec,
              const Linear::Vector *solnVec,
              const Linear::Vector *stateVec,
              const Linear::Vector *storeVec,
              const Linear::Vector *lead_current_vector,
              const Linear::Vector *junction_voltage_vector,
              const Linear::Vector *lead_current_dqdt_vector) {}

    virtual void updateAC(
              Parallel::Machine comm,
              const double frequency,
              const double fStart,
              const double fStop,
              const Linear::Vector *solnVec,
              const Linear::Vector *imaginaryVec,
              const Util::Op::RFparamsData *RFparams) {}

    virtual void updateNoise(
              Parallel::Machine comm,
              const double frequency,
              const double fStart,
              const double fStop,
              const Linear::Vector *solnVec,
              const Linear::Vector *imaginaryVec,
              const double totalOutputNoiseDens,
              const double totalInputNoiseDens,
              const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec) {}

protected:
  // used by individual measure classes to update the output variables 
  // on which they depend 
    void updateOutputVars(
      Parallel::Machine comm,
      std::vector<double> & outputVarVec,
      const double circuitTime,
      const Linear::Vector *solnVec,
      const Linear::Vector *stateVec,
      const Linear::Vector *storeVec,
      const Linear::Vector *imaginaryVec,
      const Linear::Vector *lead_current_vector,
      const Linear::Vector *junction_voltage_vector,
      const Linear::Vector *lead_current_dqdt_vector,
      const double totalOutputNoiseDens,
      const double totalInputNoiseDens,
      const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec,
      const Util::Op::RFparamsData *RFparams);

    void resetBase();

public:
  bool finishedCalculation() {
      return calculationDone_;
    }

    void makeMeasureOps(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager);

    // these functions implement measurement window criteria
    // such as TD (delay time), TO/FROM windows and Rise/Fall/Cross counts.
    // the function withinFromToWindow() is used with both time (.tran) and 
    // frequency (.ac)
    bool withinTimeWindow( double time ); // used with TRAN
    bool withinFreqWindow( double freq ); // used with AC
    bool withinRiseFallCrossWindow( double measureVal, double crossVal  );
    bool newRiseFallCrossWindowforLast();
    void setRFCValueAndFlag( Util::ParamList::const_iterator currentParamIt, int &rfcVal, bool &rfcFlag ); 
    bool withinDCsweepFromToWindow( double sweepValue ); //used with DC
    bool withinMinMaxThresh( double value);

    // used to call the output manager's getPrgetImmutableValue<int>()
    double getOutputValue(
      Parallel::Machine comm,
      Util::Op::Operator *op,
      const Linear::Vector *solnVec,
      const Linear::Vector *stateVec,
      const Linear::Vector *storeVec,
      const Linear::Vector *imaginaryVec,
      const Linear::Vector *lead_current_vector,
      const Linear::Vector *junction_voltage_vector,
      const Linear::Vector *lead_current_dqdt_vector,
      const double totalOutputNoiseDens,
      const double totalInputNoiseDens,
      const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec,
      const Util::Op::RFparamsData *RFparams);

    // used to get the measurement result
    virtual double getMeasureResult() {
      return calculationResult_;
    }

    // used to check measure line for errors that will later cause core dumps in updateTran()
    virtual bool checkMeasureLine(); 

    // used to get the variable that controls where the measure output appears
    std::string getMeasurePrintOption() { return measurePrintOption_; }

    // used to print warnings about measurement time window, etc.
    virtual void printMeasureWarnings(const double endSimTime, const double startSweepVal,
                                      const double endSweepVal);

    std::string getDCSweepVarName(const std::vector<Analysis::SweepParam> & dcParamsVec);

    bool isInvalidTimeFreqWindow(double startSimVal, double endSimVal);
    bool isInvalidDCsweepWindow(double startSweepVal, double endSweepVal);

    // used to print message about measurement time window, etc.
    virtual std::ostream& printMeasureWindow(std::ostream& os, const double endSimTime,
                                             const double startSweepVal, const double endSweepVal);
    std::string setModeStringForMeasureWindowText();
    std::string setModeStringForMeasureResultText();

    // used to print the measurement result to an output stream object, which
    // is typically the mt0, ma0 or ms0 file
    virtual std::ostream& printMeasureResult(std::ostream& os)
    {
        basic_ios_all_saver<std::ostream::char_type> save(os);
        if ( !initialized_ && measureMgr_.isMeasFailGiven() && measureMgr_.getMeasFail() )
	{
          // output FAILED to .mt file if .OPTIONS MEASURE MEASFAIL=1 is given in the
          // netlist and this is a failed measure.
          os << name_ << " = FAILED" << std::endl;
        }
        else
	{
          os << name_ << " = " << std::scientific << std::setprecision(precision_) << this->getMeasureResult() << std::endl;
        }
        return os;
    }

    // used to print the "verbose" (more descriptive) measurement result to an output stream
    // object, which is typically stdout
    virtual std::ostream& printVerboseMeasureResult(std::ostream& os)
    {
        basic_ios_all_saver<std::ostream::char_type> save(os);
        if (initialized_)
        {
          os << name_ << " = " << std::scientific << std::setprecision(precision_) << this->getMeasureResult() << std::endl;
        }
        else
        { 
          os << name_ << " = FAILED" << std::endl;
        }
        return os;
    }

    // used to print information about RFC window
    virtual std::ostream& printRFCWindow(std::ostream& os);

    // allows measureManager to access an individual measure's name and mode (e.g., AC, DC or TRAN), 
    // mainly for the purposes of error checking and reporting.
    std::string getMeasureName() { return name_; }
    std::string getMeasureMode() { return mode_;}

    const Manager & measureMgr_;

    // this is the user defined name for this measurement
    std::string name_;
    
    // this is the mode under which the measurement is active (DC, AC or TRAN)
    std::string mode_;

    // This is the type of measurement to be done (e.g., MIN or MAX).  A text string is stored here.
    std::string type_;

    // this bool is set by the constructor of the derived class.  So those that
    // are not supported (i.e. not fully implemented) can warn the user and
    // then the measure manager can based on this flag not add them to the active list.
    bool typeSupported_;
    
    // this is a flag indicating if the measure is set up to start doing it's 
    // calculation.  If it's reset to "false" then the measure should start over 
    bool initialized_;

    // this is the list of output variables needed in the measure.
    // it could be one or more values.  Since each var in the Util::Param list
    // can take up multiple spots on the list, keep a count of the number and
    // a vector of iterators pointing to the start of each var in the list
    int numDepSolVars_;
    Util::ParamList depSolVarIterVector_;
    Util::Op::OpList outputVars_;
    double outputValueTarget_;
    bool  outputValueTargetGiven_;
    double lastOutputValue_;

    // many controls on how the measure calculation is done are set via keyword=val
    // we'll parse those out and hold them in the base class:
    double td_;
    bool tdGiven_;
    double goal_;
    double weight_;
    double minval_;
    double at_;
    bool atGiven_;
    double from_;
    bool fromGiven_;
    double to_;
    bool toGiven_;
    double ymin_;
    double ymax_;
    int rise_;
    bool riseGiven_;
    int fall_;
    bool fallGiven_;
    int cross_;
    bool crossGiven_; 
    int actualRise_;
    bool isRising_;
    
    // used when Rise/Fall/Cross uses an absolute value
    bool rfcLevelGiven_;
    double rfcLevel_;

    // flag indicates that LAST was specified for Rise,Fall or Cross measure 
    bool measureLastRFC_; 
    // max and min threshold are used to set upper and lower bounds on a value
    // prior to computation.
    double maxThresh_;
    bool maxThreshGiven_;
    double minThresh_;
    bool minThreshGiven_;

    int actualFall_;
    bool isFalling_;
    int actualCross_;

    // these are used by the Duty, On and Off measures to give values as to
    // when a variable is "on" or "off"
    double onValue_;
    bool onValueGiven_;
    double offValue_;
    bool offValueGiven_;

    // in the Rise/fall/delay measure the trigger and target can be nodes of the
    // circuit so they are of the type Util::Param  The TRIG and TARG signals can
    // also use different signals.  So, all of these variables had to be define for
    // trig and targ
    Util::Param trig_;
    Util::Param targ_;
    double trigOutputValueTarget_;
    bool trigOutputValueTargetGiven_;
    double targOutputValueTarget_;
    bool targOutputValueTargetGiven_;
    double trigFracMax_;  // fraction of the maxima for the trigger value
    bool trigFracMaxGiven_;
    double targFracMax_;  // fraction of the maxima for the target value
    bool targFracMaxGiven_;

    // separate Rise/Fall/Cross values can be given for TRIG and TARG,
    // so these are separate from the riseGiven_, etc. variables defined above
    bool trigRiseGiven_; 
    bool targRiseGiven_;
    bool trigFallGiven_;
    bool targFallGiven_;
    bool trigCrossGiven_;
    bool targCrossGiven_;
    int trigRise_;
    int targRise_;
    int trigFall_;
    int targFall_;
    int trigCross_;
    int targCross_;

    // flag used to indicate whether the processing for the first step in the 
    // measurement window has been done.  
    // It was added so that the Rise/Fall/Cross count is not incremented at t=0 if the waveform 
    // has a DC offset at t=0.  Also used to record start/end Sweep Values for DC measures.
    bool firstStepInMeasureWindow_;

    // flags used to indicate that a new rise, fall or cross window has occurred.
    // Used to help make the LAST qualifier work with RISE/FALL/CROSS
    bool newRiseWindow_, newFallWindow_, newCrossWindow_;

    // flag used to indicate whether the processing for the first step in the 
    // Rise/Fall/Cross window has been done.
    // It is used to help record the start time of the Rise/Fall/Cross window.
    bool firstStepInRfcWindow_;

    // variables used to record the start-time, or start-end times of the Rise/Fall/Cross window
    // if a valid one was found.
    bool rfcWindowFound_;
    double rfcWindowStartTime_;
    double rfcWindowEndTime_;

    // variables used for evaluating FROM-TO window for DC measures, and for printing
    // out the measure window info for AC, DC and NOISE measures.
    ExtendedString sweepVar_;
    bool firstSweepValueFound_;
    bool dcSweepAscending_;

    // used for error checking (find without a when qualifer) and for differentiating
    // between a find and a find-when measure
    bool findGiven_;
    bool whenGiven_;

    // this is used by max to measure the time when a requested fraction of the
    // maximum is reached.  as in 90% of max value.  Also applies to min as well.
    double fractionToExtrema_;
    bool fractionToExtremaGiven_;

    // this is used by fourier to determine how many harmonics to compute and sample grid size.
    int numFreq_;
    int gridSize_;

    // many measurements will finish before the end of a simulation.  This flag is used
    // to indicate if the measurement is done, and the measure no longer needs to be active
    bool calculationDone_;

    // this flag is used for compatibility with the LAST keyword for RISE/FALL/CROSS.  It indicates
    // that a measure value was found, if the initialized_ flag is not used by a given measure.  
    bool resultFound_;

    // This is where the results are stored.  calculationInstant_ will be a time value for
    // TRAN mode, a frequency value for AC mode and a Sweep Variable value for DC mode.
    double calculationResult_;
    double calculationInstant_;

    // this flag is used in getMeasureResult(), if the calculations done in that function
    // (e.g., for the ERROR measure) are complicated.  It is a performance optimization.
    bool gotMeasureResult_;

    // This is the default value of the calculation
    double calculationDefaultVal_;

    // controls the precision of the displayed results
    int precision_;

    // controls where the measure output appears.  Both .mt0 file and stdout, 
    // only stdout, or neither.
    std::string measurePrintOption_;

    // controls whether the measure output (to the .mt0 file) is the measure's
    // value or the measure's calculation time.  This is only supported for
    // MIN and MAX measures for now. 
    std::string measureOutputOption_;

    //  used in the Error measure where least squared error fit to data is measured.
    std::string dataFileName_;           // name of file that may hold data
    std::string comparisonFunctionName_; // what type of comparison will be done as in L1NORM, L2NORM or INFNORM
    bool usedDefaultComparisonFunc_;     // used in descriptive output, default is true
    int independentVarColumn_;            // if the data file has more then two columns, this is the indpendent var column (defaults to 0)
    int dependentVarColumn_;              // if there are more then two columns in the data file, this is the dependent var (defaults to 1);
};

double getDCSweepVal(const std::vector<Analysis::SweepParam> & dcParamsVec);

} // namespace Measure
} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_MeasureBase_h


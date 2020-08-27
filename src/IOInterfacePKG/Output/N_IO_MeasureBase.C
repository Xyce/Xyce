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
// Purpose       : Base class measure functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iterator>

#include <N_IO_MeasureManager.h>
#include <N_IO_MeasureBase.h>
#include <N_IO_Op.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : MeasureBase::MeasureBase
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Base::Base( const Manager &measureMgr, const Util::OptionBlock & measureBlock)
  : measureMgr_(measureMgr),
    name_(""),
    mode_(""),
    type_(""),
    typeSupported_(false),
    initialized_(false),
    numDepSolVars_(0),
    outputValueTarget_(0.0),
    outputValueTargetGiven_(false),
    lastOutputValue_(0.0),
    td_(0.0),
    tdGiven_(false),
    goal_(0.0),
    weight_(1.0),
    minval_(1.0e-12),
    at_(0.0),
    atGiven_(false),
    from_(0.0),
    fromGiven_(false),
    to_(0.0),
    toGiven_(false),
    ymin_(1.0e-15),
    ymax_(1.0e+15),
    rise_(0),
    riseGiven_(false),
    fall_(0),
    fallGiven_(false),
    cross_(0),
    crossGiven_(false),
    actualRise_(0),
    isRising_(false),
    rfcLevelGiven_(false),
    rfcLevel_(0.0),
    maxThresh_(0.0),
    maxThreshGiven_(false),
    minThresh_(0.0),
    minThreshGiven_(false),
    actualFall_(0),
    isFalling_(false),
    actualCross_(0),
    onValue_(0.0),
    onValueGiven_(false),
    offValue_(0.0),
    offValueGiven_(false),
    trigOutputValueTarget_(0.0),
    trigOutputValueTargetGiven_(false),
    targOutputValueTarget_(0.0),
    targOutputValueTargetGiven_(false),
    trigFracMax_(0.0),
    trigFracMaxGiven_(false),
    targFracMax_(0.0),
    targFracMaxGiven_(false),
    trigRiseGiven_(false),
    targRiseGiven_(false),
    trigFallGiven_(false),
    targFallGiven_(false),
    trigCrossGiven_(false),
    targCrossGiven_(false),
    trigRise_(0),
    targRise_(0),
    trigFall_(0),
    targFall_(0),
    trigCross_(0),
    targCross_(0),
    firstStepInMeasureWindow_(false),
    newRiseWindow_(0),
    newFallWindow_(0),
    newCrossWindow_(0),
    firstStepInRfcWindow_(false),
    rfcWindowFound_(false),
    rfcWindowStartTime_(0.0),
    rfcWindowEndTime_(0.0),
    sweepVar_(""),
    firstSweepValueFound_(false),
    dcSweepAscending_(true),
    findGiven_(false),
    whenGiven_(false),
    fractionToExtrema_(0.0),
    fractionToExtremaGiven_(false),
    numFreq_(10),
    gridSize_(200),
    calculationDone_(false),
    resultFound_(false),
    calculationResult_(-1.0),
    calculationInstant_(0.0),
    gotMeasureResult_(false),
    calculationDefaultVal_(-1.0),
    precision_(6),
    measurePrintOption_("ALL"),
    measureOutputOption_("VALUE"),
    comparisonFunctionName_("L2NORM"),
    usedDefaultComparisonFunc_(true),
    independentVarColumn_(-1),
    dependentVarColumn_(-1)
{
  // since many of the measure types share the use of keywords (like TD=<delay time>) we'll
  // parse those out here and store them.  We'll also pull out the out_var[= val | out_var2]
  // here since once we eliminate the keyword from the measureBlock, that's all that should left.

  // these are used to mark if we're in the trigger or target sections of a rise/fall/delay
  // measure.  This measure is different from the others in that it is essentiality
  // a compound measure which requires two objectives and thus two sets of qualifiers
  // on the objective like frac_max=, weight= etc.
  bool inTrigBlock = false;
  bool inTargBlock = false;

  // apply any "global settings" from .OPTIONS MEASURE lines, which will take
  // precedence over any qualifier values on the individual .MEASURE lines.
  if ( measureMgr_.isMeasDgtGiven() ) 
  { 
    // .OPTIONS MEASURE MEASDGT=<val> overrides PRECISION qualifier,
    // if it is given in the netlist.
    precision_ = measureMgr_.getMeasDgt();
  }

  if ( measureMgr_.isMeasFailGiven() ) 
  { 
    // .OPTIONS MEASURE MEASFAIL=<val> overrides the DEFAULT_VAL qualifier
    // and .OPTIONS MEASURE DEFAULT_VAL=<val>, if it is given in the netlist.  
    // Substitute 0 for the measure's default value.  The actual value output 
    // to the .mt0 file will then be 0 if MEASFAIL=FALSE and "FAILED" if 
    // MEASFAIL==TRUE. 
    calculationDefaultVal_ = 0;  
    calculationResult_ = calculationDefaultVal_;
  }
  else if ( measureMgr_.isMeasGlobalDefaultValGiven() ) 
  {
    // .OPTIONS MEASURE DEFAULT_VAL=<val> overrides the DEFAULT_VAL qualifier,
    // if it is given in the netlist
    calculationDefaultVal_ = measureMgr_.getMeasGlobalDefaultVal();
    calculationResult_ = calculationDefaultVal_;
  }

  for (Util::ParamList::const_iterator it = measureBlock.begin(), end = measureBlock.end(); it != end; ++it)
  {
    const std::string &tag = (*it).tag();
    //Xyce::dout() << " in measure base setup: tag = \"" << tag << "\"" << std::endl;
    if( tag == "NAME" )
    {
      name_ = (*it).stringValue();
    }
    else if( tag == "MODE" )
    {
      mode_ = (*it).stringValue();
    }
    else if( tag == "TYPE" )
    {
      type_ = (*it).stringValue();
      // if type is TRIG or TARG then the next Util::Param in the list the
      // node or expression that is the Trigger or Target.  We'll need to
      // catch and save this.  This oddity arises because all of the other
      // measures have one objective to read while the TRIG/TARG of the Rise
      // fall, delay measure has two and they are thus named.  Also, some of
      // the qualifiers like value=, weight= and frac_max= apply to the
      // last TRIG ar TARG statement so we need to remember if we just
      // passed one of those.
      if( type_ == "TRIG" )
      {
        inTrigBlock = true;
        inTargBlock = false;
      }
      else if( type_ == "TARG" )
      {
        inTrigBlock = false;
        inTargBlock = true;
      }
      else if( type_ == "FIND")
      {
        findGiven_ = true;
      }
      else if( type_ == "WHEN")
      {
        whenGiven_ = true;
      }
    }
    else if( tag == "TD" )
    {
      td_ = (*it).getImmutableValue<double>();
      tdGiven_ = true;
    }
    else if( ( tag == "GOAL" ) || ( tag == "VALUE" ) )
    {
      goal_ = (*it).getImmutableValue<double>();
    }
    else if( tag == "WEIGHT" )
    {
      weight_ = (*it).getImmutableValue<double>();
    }
    else if( tag == "MAX_THRESH" )
    {
      maxThresh_ = (*it).getImmutableValue<double>();
      maxThreshGiven_ = true;
    }
    else if( tag == "MIN_THRESH" )
    {
      minThresh_ = (*it).getImmutableValue<double>();
      minThreshGiven_ = true;
    }
    else if( tag == "MINVAL" )
    {
      minval_ = (*it).getImmutableValue<double>();
    }
    else if( tag == "AT" )
    {
      at_ = (*it).getImmutableValue<double>();
      atGiven_ = true;
      outputValueTargetGiven_ = true;
      if ( inTargBlock )
      {
        Report::UserError0() << "AT keyword not allowed in TARG block for measure " << name_ ;
      }
    }
    else if( tag == "FROM" )
    {
      from_ = (*it).getImmutableValue<double>();
      fromGiven_ = true;
    }
    else if( tag == "TO" )
    {
      to_ = (*it).getImmutableValue<double>();
      toGiven_ = true;
    }
    else if( (tag == "IGNORE") || (tag == "IGNOR") || (tag == "YMIN") )
    {
      ymin_ = (*it).getImmutableValue<double>();
      if (ymin_ < 0)
      {
        Report::UserError0() << "YMIN or IGNORE keyword for measure " << name_ << " must be non-negative";
      }

    }
    else if( tag == "YMAX" )
    {
      ymax_ = (*it).getImmutableValue<double>();
      if (ymax_ <= 0)
      {
        Report::UserError0() << "YMAX keyword for measure " << name_ << " must be positive";
      }
    }
    else if( tag == "RISE" )
    {
      if ( inTrigBlock )
      {
        setRFCValueAndFlag(it, trigRise_, trigRiseGiven_);
      }
      else if (inTargBlock)
      {
        setRFCValueAndFlag(it, targRise_, targRiseGiven_);
      }
      else
      {
        setRFCValueAndFlag(it, rise_, riseGiven_);
      }
    }
    else if( tag == "FALL" )
    {
      if ( inTrigBlock )
      {
        setRFCValueAndFlag(it, trigFall_, trigFallGiven_);
      }
      else if (inTargBlock)
      {
        setRFCValueAndFlag(it, targFall_, targFallGiven_);
      }
      else
      {
        setRFCValueAndFlag(it, fall_, fallGiven_);
      }
    }
    else if( tag == "CROSS" )
    {
      if ( inTrigBlock )
      {
        setRFCValueAndFlag(it, trigCross_, trigCrossGiven_);
      }
      else if (inTargBlock)
      {
        setRFCValueAndFlag(it, targCross_, targCrossGiven_);
      }
      else
      {
        setRFCValueAndFlag(it, cross_, crossGiven_);
      }
    }
    else if( tag == "RFC_LEVEL" )
    {
      rfcLevel_ = (*it).getImmutableValue<double>();
      rfcLevelGiven_ = true;
    }
    else if( tag == "FRAC_MAX" )
    {
      if( inTrigBlock )
      {
        trigFracMax_ = (*it).getImmutableValue<double>();
        trigFracMaxGiven_ = true;
      }
      else if( inTargBlock )
      {
        targFracMax_ = (*it).getImmutableValue<double>();
        targFracMaxGiven_ = true;
      }
    }
    else if( Xyce::Util::hasExpressionTag(tag) )
    {
      numDepSolVars_++;
      depSolVarIterVector_.push_back(*it);
    }
    else if( tag == "OBJVAL" )
    {
      if( (*it).getType() == Xyce::Util::INT )
      {
        outputValueTarget_ = (*it).getImmutableValue<int>();
        outputValueTargetGiven_ = true;
      }
      else if( (*it).getType() == Xyce::Util::DBLE )
      {
        outputValueTarget_ = (*it).getImmutableValue<double>();
        outputValueTargetGiven_ = true;
      }
      else if( (*it).getType() == Xyce::Util::STR )
      {
        // a bare string name that we will have to resovle
        Util::Param aParam;
        aParam.set( (*it).stringValue(), 0 );
        numDepSolVars_++;
        depSolVarIterVector_.push_back(aParam);
      }
      else if ( (*it).getType() == Xyce::Util::EXPR )
      {
        Util::Param aParam;
        aParam.set( '{' + (*it).stringValue() + '}', 0 );
        numDepSolVars_++;
        depSolVarIterVector_.push_back(aParam);
      }

      if( inTrigBlock )
      {
        trigOutputValueTarget_ = outputValueTarget_;
        trigOutputValueTargetGiven_ = true;
        outputValueTargetGiven_ = false;
      }
      else if( inTargBlock )
      {
        targOutputValueTarget_ = outputValueTarget_;
        targOutputValueTargetGiven_ = true;
        outputValueTargetGiven_ = false;
      }
    }
    else if( tag == "ON" )
    {
      onValue_ = (*it).getImmutableValue<double>();
      onValueGiven_ = true;
    }
    else if( tag == "OFF" )
    {
      offValue_ = (*it).getImmutableValue<double>();
      offValueGiven_ = true;
    }
    else if( tag == "NUMFREQ" )
    {
      numFreq_ = (*it).getImmutableValue<int>();
    }
    else if( tag == "GRIDSIZE" )
    {
      gridSize_ = (*it).getImmutableValue<int>();
    }
    else if( tag == "FILE" )
    {
      dataFileName_ = (*it).stringValue();
    }
    else if( tag == "COMP_FUNCTION" )
    {
      // default value for comparisonFunctionName_ is "L2NORM", and it is set
      // in the constructor. The other two recognized strings are "L1NORM" and 
      // "INFNORM". Any other string defaults to "L2NORM".
      ExtendedString tmpStr((*it).stringValue());
      if (tmpStr.toUpper() == "L1NORM")
      {
        comparisonFunctionName_ = "L1NORM";
        usedDefaultComparisonFunc_ = false;
      }
      else if (tmpStr.toUpper() == "INFNORM")
      {
        comparisonFunctionName_ = "INFNORM";
        usedDefaultComparisonFunc_ = false;
      }
      else if (tmpStr.toUpper() == "L2NORM")
      {
        // L2NORM is the default value, but the usedDefaultComparisonFunc_
        // flag is used in the descriptive output for the ERROR measure.  It
        // is set to false here because L2NORM was explicitly specifed on the
        // measure instance line.
        usedDefaultComparisonFunc_ = false;
      }
    }
    else if( tag == "INDEPVARCOL" )
    {
      independentVarColumn_ = (*it).getImmutableValue<int>();
    }
    else if( tag == "DEPVARCOL" )
    {
      dependentVarColumn_ = (*it).getImmutableValue<int>();
    }
    else if( tag == "DEFAULT_VAL" )
    {
      // use the value from the DEFAULT_VAL qualifier on the .MEASURE line,
      // unless .OPTION MEASURE MEASFAIL=<val> or .OPTION MEASURE DEFAULT_VAL=<val>
      // are in the netlist
      if ( !(measureMgr_.isMeasFailGiven() || measureMgr_.isMeasGlobalDefaultValGiven()) )
      {
        calculationDefaultVal_ = (*it).getImmutableValue<double>();
        calculationResult_ = calculationDefaultVal_;
      }
    }
    else if( tag == "PRECISION" ) 
    {
      // use the value from the PRECISION qualifier on the .MEASURE line,
      // unless .OPTION MEASURE MEASDGT=<val> is in the netlist
      if ( !measureMgr_.isMeasDgtGiven() )
      {
        precision_ = (*it).getImmutableValue<int>();
      }
    }
    else if( tag == "PRINT" )
    {
      // default value for measurePrintOption_ is "ALL", and it is set in the constructor.
      // The other two recognized strings are "NONE" and "STDOUT".
      // NONE suppresses the measure output in both .mt0 and standard output.
      // STDOUT suppress the measure output in just .mt0.
      // Any other string defaults to "ALL".
      ExtendedString tmpStr((*it).stringValue());
      if ( tmpStr.toUpper() == "STDOUT")
      {
        measurePrintOption_ = "STDOUT";
      }
      else if (tmpStr.toUpper() == "NONE")
      {
        measurePrintOption_ = "NONE";
      }
    }
    else if( tag == "OUTPUT" )
    {
      // default value for measureOutput_ is "VALUE", and it is set in the constructor.
      // The other recognized strings are "TIME", "FREQ" and "SV", which are basically 
      // synonyms that are "human-readable" for TRAN, AC and DC modes respectively.
      // Any other string defaults to "VALUE".
      ExtendedString tmpStr((*it).stringValue());
      if ( tmpStr.toUpper() == "TIME")
      {
        measureOutputOption_ = "TIME";
      }
      else if ( tmpStr.toUpper() == "FREQ")
      {
        measureOutputOption_ = "FREQ";
      }
      else if ( tmpStr.toUpper() == "SV")
      {
        measureOutputOption_ = "SV";
      }
    }
    else if( tag[0]=='V' || (tag[0]=='I' && tag != "INOISE") || tag[0]=='N' ||
             tag[0]=='P' || tag[0]=='W' || tag[0]=='D')
    {
      // this if clause must come last because we are only checking the 
      // first letter and don't with to get confused with kewords 
      // that happen to start with V, I, N, P, W or D.
      int nodes = (*it).getImmutableValue<int>();
      Util::Param aParam;
      aParam.set( tag, nodes );

      if( inTrigBlock )
      {
        trig_ = aParam;
      }
      else if( inTargBlock )
      {
        targ_ = aParam;
      }

      // at this point trig and targ hold the beginning of the argument and how
      // many nodes it depends on (typically 1 or 2)  That's all we need there
      // to later use outputVars_ and depSolVarIterVector_ to get their values.

      // here we just store the needed parts of V(a) or v(a,b) or I(device).
      // only the v(a,b) case will need an extra node in the outputVars_ array.

      numDepSolVars_++;

      depSolVarIterVector_.push_back(aParam);
      for( int i=0; i<nodes; i++ )
      {
        it++;
        aParam.set( (*it).tag(), (*it).getImmutableValue<double>() );
        depSolVarIterVector_.push_back( aParam );
      }
    }
    else if( tag[0]=='S' || tag[0]=='Y' || tag[0]=='Z')
    {
      // Turn the S, Y or Z operators into an expression-valued parameter.
      // We need to do this because for SR(1,1) the values 1 and 1 are
      // indices rather than (for example) nodes as in V(1,2).  So, the
      // dependent solution variable becomes the expression {SR(1,1)}
      // in this case.
      std::string expString('{'+tag+'(');
      int nodes = (*it).getImmutableValue<int>();

      for( int i=0; i< nodes; i++)
      {
        it++;
        expString += (*it).tag();
        if (i == 0) expString += ',';
      }
      expString += ")}";

      Util::Param aParam;
      aParam.set(expString,expString);

      numDepSolVars_++;
      depSolVarIterVector_.push_back(aParam);
    }
    else
    {
      Xyce::Report::UserWarning() << "Unknown tag in measure statement: " << tag << ", ignoring";
    }
  }

  // set flag related to whether LAST was specified for Measure for RISE, FALL or CROSS
  measureLastRFC_ = ((riseGiven_ && rise_ < 0) || (fallGiven_ && fall_ < 0) || 
	             (crossGiven_ && cross_ < 0)) ? true : false;
}

Base::~Base()
{
  for (Util::Op::OpList::iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
    delete *it;
}

void
Base::makeMeasureOps(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager) 
{
  makeOps(comm, op_builder_manager, NetlistLocation(), depSolVarIterVector_.begin(), depSolVarIterVector_.end(), std::back_inserter(outputVars_));

  prepareOutputVariables();
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinTimeWindow
// Purpose       : Checks if current time is within TD and FROM/TO windows.
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/8/2014
//-----------------------------------------------------------------------------
bool Base::withinTimeWindow( double time )
{
  bool retVal = true;
  if ( ( tdGiven_ && (time < td_)) || (fromGiven_ && (time < from_ )) || (toGiven_ && (time > to_)) )
  {
    retVal = false;
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinFreqWindow
// Purpose       : Checks if current frequency is within FROM/TO windows.
// Special Notes : The minVal_ tolerance is used as a fudge factor on the TO
//                 window because of numerical errors in the sweep values
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2019
//-----------------------------------------------------------------------------
bool Base::withinFreqWindow( double freq )
{
  bool retVal = true;
  if ( (fromGiven_ && (freq < from_ )) || (toGiven_ && (freq > (1+minval_)*to_)) )
  {
    retVal = false;
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinDCsweepFromToWindow
// Purpose       : Checks if the current value of the DC Sweep variable is
//                 within measurement window
// Special Notes : Used just for DC mode, since the first sweep variable (on a .DC
//                 line) can be either monotonically increasing or monotonically 
//                 decreasing.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 04/26/2017
//-----------------------------------------------------------------------------
bool Base::withinDCsweepFromToWindow(double sweepValue)
{
  // function used for DC mode
  bool retVal = true;

  if (fromGiven_ && toGiven_)
  {
    if( to_ >= from_ )
    {
      if ( (sweepValue < from_ ) || (sweepValue > to_) )
      {
        retVal = false;
      }
    }
    else
    {
      if( (sweepValue > from_ ) || (sweepValue < to_) )
      {
        retVal = false;
      }
    }
  }
  else if (toGiven_)
  {
    // handle both ascending and descending sweeps
    if ( (dcSweepAscending_ && (sweepValue > to_)) || (!dcSweepAscending_ && (sweepValue < to_)) )
      retVal = false;
  }
  else if (fromGiven_)
  {
    // handle both ascending and descending sweeps
    if ( (dcSweepAscending_ && (sweepValue < from_)) || (!dcSweepAscending_ && (sweepValue > from_)) )
      retVal = false;
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinRiseFallCrossWindow
// Purpose       : Checks if current value is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Base::withinRiseFallCrossWindow( double measureVal, double crossVal )
{
  // return true if neither rise, fall or cross is given.
  bool retVal = true;

  // used to enable the LAST keyword.  Reset to false, each time 
  // this function is called.
  newRiseWindow_=false;
  newFallWindow_=false;
  newCrossWindow_=false;

  // update rise/fall/cross counts
  if( riseGiven_ || fallGiven_ || crossGiven_ )
  {
    // default is to return false if a rise, fall or cross was given.
    retVal = false;

    // first check if we need to adjust rise/fall/cross counts.

    // use the target value in the WHEN clause (for the DERIV, FIND-WHEN and WHEN measures),
    // or the RFC_LEVEL if one was specified.  The former is for HSpice compatibility.
    if (whenGiven_ || rfcLevelGiven_ )
    {
      // for HSpice compatibility, a rising (or falling) waveform that equals the cross level
      // is considered to a rise (or fall).
      if ( ((measureVal-crossVal) >= 0.0) && ((lastOutputValue_-crossVal) < 0.0) )
      {
        actualRise_++;
        newRiseWindow_=true;
      }
      else if( ((measureVal-crossVal) <= 0.0) && ((lastOutputValue_-crossVal) > 0.0) )
      {   
        actualFall_++;
        newFallWindow_=true;
      }
    }
    else
    {
      // sense "absolute" rise and fall, otherwise.  This was the method supported
      // in Xyce 6.4, and earlier, for all measures.
      if( (measureVal > lastOutputValue_) && !isRising_ )
      {
        // we've started a rise
        isRising_= true;
        isFalling_ = false;
        actualRise_++;
        newRiseWindow_=true;
      }
      if( (measureVal < lastOutputValue_) && !isFalling_ )
      {
        // we've started a fall
        isRising_ = false;
        isFalling_ = true;
        actualFall_++;
        newFallWindow_=true;
      }
    }
    
    // CROSS qualifier always uses level-crossing approach. For HSpice compatibility, 
    // a rising (or falling) waveform that equals the cross level is considered to be
    // a cross.
    if( (((measureVal-crossVal) <= 0.0) && ((lastOutputValue_-crossVal) > 0.0)) 
     || (((measureVal-crossVal) >= 0.0) && ((lastOutputValue_-crossVal) < 0.0)) )
    {
      // we've crossed measureVal-crossVal == 0 
      actualCross_++;
      newCrossWindow_=true;
    }

    // now check if we're in the right window
    // this could be compressed to one statement, but this looks clearer.
    // Also note (as an example) that the RISE=LAST syntax sets riseGiven_ to true and
    // the rise_ parameter to -1.  So, we only want to return true in that case
    // if at least one rise has been found.  Similar logic applies to FALL=LAST and
    // CROSS=LAST.
    if( riseGiven_ && ( ((rise_ < 0) && (actualRise_ > 0)) || (rise_ == actualRise_) ) ) 
    {
      retVal=true;
    }
    else if( fallGiven_ && ( ((fall_ < 0) && (actualFall_ > 0)) || (fall_ == actualFall_) ) )
    {
      retVal=true;
    }
    else if( crossGiven_ && ( ((cross_< 0) && (actualCross_ > 0)) || (cross_ == actualCross_) ) )
    {
      retVal=true;
    }
    lastOutputValue_=measureVal;
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::newRiseFallCrossWindowForLast
// Purpose       : Determines if this is the start of a new
//                 Rise, Fall or Cross window.  It is used
//                 if the LAST keyword was specified.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 01/19/2016
//-----------------------------------------------------------------------------
bool Base::newRiseFallCrossWindowforLast()
{
  bool retVal=false;

  // the variables newRiseWindow_, newFallWindow_ and newCrossWindow_
  // were set in Base::withinRiseFallCrossWindow().
  if ( (rise_ < 0 && riseGiven_ && newRiseWindow_) || 
       (fall_ < 0 && fallGiven_ && newFallWindow_) || 
       (cross_ < 0 && crossGiven_ && newCrossWindow_) ) 
  {
    retVal=true;
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinMinMaxThresh
// Purpose       : Check if value is within MIN_THRESHOLD and MAX_THRESHOLD
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Base::withinMinMaxThresh( double value)
{
  bool returnValue = true;
  if( (minThreshGiven_ && (value < minThresh_)) )
    returnValue = false;
  if( (maxThreshGiven_ && (value > maxThresh_)) )
    returnValue = false;

  return returnValue;
}



//-----------------------------------------------------------------------------
// Function      : MeasureBase::updateOutputVars
// Purpose       : Call's the N_UTL_Op's getValue() function to update 
//                 the objects in Util::ParamList outputVars_;
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 11/01/2013
//-----------------------------------------------------------------------------
void Base::updateOutputVars(
  Parallel::Machine comm,
  std::vector<double> & outputVarVec,
  const double circuitTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector * storeVec,
  const Linear::Vector *imaginaryVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector,
  const double totalOutputNoiseDens,
  const double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec,
  const Util::Op::RFparamsData *RFparams)
{
  int vecIndex = 0;
  for (std::vector<Util::Op::Operator *>::const_iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
  {
    outputVarVec[vecIndex] = getValue(comm, *(*it), Util::Op::OpData(vecIndex, solnVec, imaginaryVec, stateVec, storeVec, 0, lead_current_vector, 0, junction_voltage_vector, 0, 0, 0, 0, 0, 0, totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, RFparams)).real();
    vecIndex++;
  }
}


//-----------------------------------------------------------------------------
// Function      : MeasureBase::resetBase
// Purpose       : When a measure is reset during .step, some actions need to
//                 be done in the base class to.  Put them here.
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 09/10/2014
//-----------------------------------------------------------------------------
void Base::resetBase()
{
  initialized_=false;
  calculationDone_=false;
  resultFound_=false;
  gotMeasureResult_=false;

  // reset any vars that were accumulators.
  actualRise_=0;
  actualFall_=0;
  actualCross_=0;
  isRising_=false;
  isFalling_=false;

  // reset values associated with sweep vector values
  firstSweepValueFound_=false;
  dcSweepAscending_=true;

  // reset default values
  calculationResult_ = calculationDefaultVal_;
  calculationInstant_=0;
  lastOutputValue_=0;
  firstStepInMeasureWindow_ = false;
  firstStepInRfcWindow_ = false;
  newRiseWindow_=false;
  newFallWindow_=false;
  newCrossWindow_=false;
  rfcWindowFound_ = false;
  rfcWindowStartTime_ = 0.0;
  rfcWindowEndTime_ = 0.0;
}


//-----------------------------------------------------------------------------
// Function      : MeasureBase::getOutputValue
// Purpose       : Call's the N_UTL_Op's getValue function to get sol. vars.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
double Base::getOutputValue(
  Parallel::Machine comm,
  Util::Op::Operator *op,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector * storeVec,
  const Linear::Vector *imaginaryVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector,
  const double totalOutputNoiseDens,
  const double totalInputNoiseDens,
  const std::vector<Xyce::Analysis::NoiseData*> *noiseDataVec,
  const Util::Op::RFparamsData *RFparams)
{
  double retVal = getValue(comm, *op, Util::Op::OpData(0, solnVec, imaginaryVec, stateVec, storeVec, 0, lead_current_vector, 0, junction_voltage_vector,0, 0, 0, 0, 0, 0, totalOutputNoiseDens, totalInputNoiseDens, noiseDataVec, RFparams)).real();
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::printMeasureWarnings
// Purpose       : prints error message related to invalid time windows, etc.
//                 This function currently only applies to TRAN, TRAN_CONT,
//                 AC, AC_CONT, NOISE and NOISE_CONT modes.  It does not do
//                 error checking for DC or DC_CONT modes.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2015
//-----------------------------------------------------------------------------
void Base::printMeasureWarnings(const double endSimTime, const double startSweepVal,
                                const double endSweepVal)
{
  if ( (calculationResult_ == calculationDefaultVal_) &&
       ( (mode_ == "TRAN") || (mode_ == "TRAN_CONT") || (mode_ == "AC") || (mode_ == "AC_CONT") ||
         (mode_ == "NOISE") || (mode_ == "NOISE_CONT") ) )
  {
    // print warning if time window or AT value was non-sensensical 
    if ( fromGiven_ && !tdGiven_ && toGiven_ && to_ < from_)
    {
      Xyce::Report::UserWarning() << name_ << " failed. TO value < FROM value";
    }
    else if ( tdGiven_ && toGiven_ && td_ > to_)
    {
      Xyce::Report::UserWarning() << name_ << " failed. TD value > TO value";
    } 
    else if ( toGiven_ && to_ <= 0.0 ) 
    {
      Xyce::Report::UserWarning() << name_ << " failed. TO value <= 0";
    }
    else if ( (riseGiven_ && actualRise_ < rise_) || (fallGiven_ && actualFall_ < fall_) ||
              (crossGiven_ && actualCross_ < cross_))
    {
        Xyce::Report::UserWarning() << name_ << " failed. Measured Rise,Fall,Cross=(" << 
          actualRise_ << "," << actualFall_ << "," << actualCross_ << ")";
    }    
    else if ( (mode_ == "TRAN") || (mode_ == "TRAN_CONT") )
    {
      if ( ( fromGiven_ && from_ >= endSimTime ) || ( tdGiven_ && td_ >= endSimTime ) )
      {
        Xyce::Report::UserWarning() << name_ << " failed. FROM or TD value > sim end time";
      }
      else if ( atGiven_ && (at_ < 0 || at_ > endSimTime) ) 
      {
        Xyce::Report::UserWarning() << name_ << " failed. AT value outside sim window";
      }
      else if ( atGiven_ && (at_ < from_ || at_ > to_) )
      {
        Xyce::Report::UserWarning() << name_ << " failed. AT value outside measurement window";
      }
    }
    else if ( (mode_ == "AC") || (mode_ == "AC_CONT") || (mode_ == "NOISE") || (mode_ == "NOISE_CONT") )
    {
      if ( ( fromGiven_ && from_ >= endSweepVal ) || ( tdGiven_ && td_ >= endSweepVal ) )
      {
        Xyce::Report::UserWarning() << name_ << " failed. FROM value > highest frequency value";
      }
      else if ( atGiven_ && (at_ < startSweepVal || at_ > endSweepVal) )
      {
        Xyce::Report::UserWarning() << name_ << " failed. AT value outside frequency sweep window";
      }
      else if ( atGiven_ && (at_ < from_ || at_ > to_) )
      {
        Xyce::Report::UserWarning() << name_ << " failed. AT value outside measurement window";
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::getDCSweepVarName
// Purpose       : Used to get the name of the DC sweep variable that is
//                 used for the measurement window.
// Special Notes : For TABLE-based sweeps (for .DC data=table), it is the row
//                 index for that table.  For all other DC sweep types, it
//                 is the name of the first variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 06/16/2020
//-----------------------------------------------------------------------------
std::string Base::getDCSweepVarName(const std::vector<Analysis::SweepParam> & dcParamsVec)
{
  ExtendedString sweepVarName("");

  if (dcParamsVec[0].type == "TABLE")
    sweepVarName = "Table Row";
  else
  {
    sweepVarName = dcParamsVec[0].name;
    sweepVarName.toUpper();
  }

  return sweepVarName;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::isInvalidTimeWindow
// Purpose       : returns true if the specified from_ and to_ values do NOT
//                 form a valid measurement window, based on the start/stop
//                 simulation times.
// Special Notes : Assumes that all transient simulations start at t=0
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/18/2020
//-----------------------------------------------------------------------------
bool Base::isInvalidTimeWindow(double endSimTime)
{
  return ( (fromGiven_&& toGiven_ && (from_ > to_)) ||
           (fromGiven_&& (from_ > endSimTime)) || (toGiven_&& (to_ < 0.0)) );
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::isInvalidFreqWindow
// Purpose       : returns true if the specified from_ and to_ values do NOT
//                 form a valid measurement window, based on the start/stop
//                 simulation frequencies.  So, this function can be used by
//                 AC and NOISE measure modes.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 08/06/2020
//-----------------------------------------------------------------------------
bool Base::isInvalidFreqWindow(double fStart, double fStop)
{
  return ( (fromGiven_&& toGiven_ && (from_ > to_)) ||
           (fromGiven_&& (from_ > fStop)) || (toGiven_&& (to_ < fStart)) );
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::isInvalidDCsweepWindow
// Purpose       : returns true if the specified from_ and to_ values do NOT
//                 form a valid measurement window, based on the sweep direction
//                 (ascending or descending) and the start/end DC sweep values.
// Special Notes : This can only happen if both FROM and TO are given, and that
//                 FROM-TO range does not overlap with the DC sweep range.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 08/06/2020
//-----------------------------------------------------------------------------
bool Base::isInvalidDCsweepWindow(double startSweepVal, double endSweepVal)
{
  bool retVal =false;

  if (fromGiven_ && toGiven_)
  {
    // this if-else block could be collapsed into one if statement, but this is clearer.
    if (dcSweepAscending_)
    {
      retVal = ( ((from_ > endSweepVal) && (to_ > endSweepVal)) ||
                 ((from_ < startSweepVal) && (to_ < startSweepVal)) );
    }
    else
    {
      // the descending case
      retVal = ( ((from_ < endSweepVal) && (to_ < endSweepVal)) ||
                 ((from_ > startSweepVal) && (to_ > startSweepVal)) );
    }
  }

  return retVal;

}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::printMeasureWindow
// Purpose       : prints information related to time, frequency or DC sweep
//                 window used.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2015
//-----------------------------------------------------------------------------
std::ostream& Base::printMeasureWindow(std::ostream& os, const double endSimTime,
				       const double startSweepVal, const double endSweepVal)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  double startOfWindow=0;
  double endOfWindow=0;

  if ( (mode_ == "TRAN") || (mode_ == "TRAN_CONT") )
  {
    if (tdGiven_ || fromGiven_) 
    {
      startOfWindow = (td_ > from_) ? td_ : from_;
    }

    endOfWindow = (toGiven_) ? to_ : endSimTime;
  }
  else if ( (mode_== "AC") || (mode_ == "AC_CONT") || (mode_== "NOISE") || (mode_ == "NOISE_CONT") )
  {
    // handle AC and NOISE cases
    if (initialized_)
    {
      startOfWindow = (fromGiven_) ? std::max(from_,startSweepVal) : startSweepVal;
      endOfWindow = (toGiven_) ? std::min(to_,endSweepVal) : endSweepVal;
    }
    else
    {
      startOfWindow = startSweepVal;
      endOfWindow = endSweepVal;
    }
  }
  else if ( (mode_ == "DC") || (mode_ == "DC_CONT") )
  {
    if (initialized_)
    {
      // adjust the FROM and TO values, printed to stdout, based on the direction
      // (ascending or descending) of the DC sweep.
      double fromVal, toVal;
      if ( (fromGiven_ && toGiven_ && dcSweepAscending_ && (from_ > to_)) ||
           (fromGiven_ && toGiven_ && !dcSweepAscending_ && (from_ < to_)) )
      {
        fromVal = to_;
        toVal = from_;
      }
      else
      {
        fromVal = from_;
        toVal = to_;
      }

      // account for both sweep directions, relative to FROM and TO values (if given).
      if (dcSweepAscending_)
      {
        startOfWindow = (fromGiven_) ? std::max(fromVal,startSweepVal) : startSweepVal;
        endOfWindow = (toGiven_) ? std::min(toVal,endSweepVal) : endSweepVal;
      }
      else
      {
        startOfWindow = (fromGiven_) ? std::min(fromVal,startSweepVal) : startSweepVal;
        endOfWindow = (toGiven_) ? std::max(toVal,endSweepVal) : endSweepVal;
      }
    }
    else
    {
      startOfWindow = startSweepVal;
      endOfWindow = endSweepVal;
    }
  }

  // modeStr is "Time" for TRAN and TRAN_CONT modes, "Freq" for AC and AC_CONT modes and
  // "<sweep variable> Value" for DC and DC_CONT modes.
  if ( (mode_ == "AC") || (mode_ == "NOISE") || (mode_ == "TRAN") || ((mode_ == "DC") && firstSweepValueFound_) ||
       (mode_ == "AC_CONT") || (mode_ == "NOISE_CONT") || (mode_ == "TRAN_CONT") ||
       ((mode_ == "DC_CONT") && firstSweepValueFound_))
  {
    std::string modeStr = setModeStringForMeasureWindowText();
    os << "Measure Start " << modeStr << "= " << startOfWindow 
       << "\tMeasure End " << modeStr << "= " << endOfWindow << std::endl; 
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::setModeStringForMeasureWindowText()
// Purpose       : set text string used in various printMeasureWindow() functions. 
// Special Notes : modeStr is "Time" for TRAN or TRAN_CONT mode, "Freq" for AC or
//                 AC_CONT mode and "<sweep variable> Value" for DC or DC_CONT mode.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::string Base::setModeStringForMeasureWindowText()
{
  std::string modeStr;
  if ( (mode_ == "TRAN") || (mode_ == "TRAN_CONT") )
  {
    modeStr = "Time";
  }
  else if ( (mode_ == "AC")  || (mode_ == "AC_CONT") || (mode_ == "NOISE")  || (mode_ == "NOISE_CONT"))
  {
    modeStr = "Freq";
  }
  else
  {
    // DC or DC_CONT case
    modeStr = sweepVar_ + " Value";
  }

  return modeStr;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::setModeStringForMeasureResultText()
// Purpose       : set text string used in various printMeasureResult() functions. 
// Special Notes : modeStr is "time" for TRAN or TRAN_CONT mode, "freq" for AC
//                 or AC_CONT mode and "<sweep variable> value" for DC or DC_CONT
//                 mode.
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::string Base::setModeStringForMeasureResultText()
{
  std::string modeStr;
  if ( (mode_ == "TRAN") || (mode_ == "TRAN_CONT") )
  {
    modeStr = "time";
  }
  else if ( (mode_ == "AC") || (mode_ == "AC_CONT") || (mode_ == "NOISE") || (mode_ == "NOISE_CONT"))
  {
    modeStr = "freq";
  }
  else
  {
    // DC case or DC_CONT case
    modeStr = sweepVar_ + " value";
  }

  return modeStr;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::printRFCWindow()
// Purpose       : print informaiton about the start time of the RISE, FALL or CROSS
//                 window, if a valid one was found.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystems Modeling
// Creation Date : 09/21/2015
//-----------------------------------------------------------------------------
std::ostream& Base::printRFCWindow(std::ostream& os)
{
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os << std::scientific << std::setprecision(precision_);

  // Printing of the RFC window information may not work for all
  // measures if LAST is specified.
  if ( ( (riseGiven_ && (actualRise_ > 0) ) || ( fallGiven_ && (actualFall_ > 0 ) )
	 || (crossGiven_ && (actualCross_ > 0) ) )  && rfcWindowFound_ )
  {
    if (riseGiven_)
    {
      rise_ < 0 ? ( os << "Last Rise" ) : ( os << "Rise " << rise_ ) ;
    }
    else if (fallGiven_)
    {
      fall_ < 0 ? ( os << "Last Fall" ) :( os << "Fall " << fall_ ) ;
    }
    else if (crossGiven_)
    {
      cross_ < 0 ? ( os << "Last Cross" ) : ( os << "Cross " << cross_) ;
    }
    os << ": Start Time= " << rfcWindowStartTime_ << "\tEnd Time= " << rfcWindowEndTime_ << std::endl; 
  }
  else
  { 
    // no op 
  }     

  return os;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::setRFCValueAndFlag
// Purpose       : sets the value and given_ flag for a rise, fall, cross count.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/11/2015
//-----------------------------------------------------------------------------
void Base::setRFCValueAndFlag( Util::ParamList::const_iterator currentParamIt, int &rfcVal, bool &rfcFlag )
{ 
  if( currentParamIt->getType() == Xyce::Util::STR )
  {
    ExtendedString lastStr(currentParamIt->stringValue());
    if (lastStr.toUpper() == "LAST")
    {
      // user requested LAST rise/fall/cross in simulation
      // so measure all of them and just report the last one.
      rfcVal = -1;
    }
    else
    {
      Report::UserError0() << "Invalid value for RISE, FALL or CROSS for measure " << name_ ;
    }
  }
  else
  {
    rfcVal = currentParamIt->getImmutableValue<int>();
  }
  rfcFlag = true;
}

//-----------------------------------------------------------------------------
// Function      : MeasureBase::checkMeasureLine
// Purpose       : check .MEASURE line for errors that will cause cause dumps
//               : later
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 09/01/2015
//-----------------------------------------------------------------------------
bool Base::checkMeasureLine()
{
  bool bsuccess = true;
  // incorrect number of dependent solution variables will cause core dumps in
  // updateTran() function
  if (numDepSolVars_ == 0)
  {
    bsuccess = false;
    Report::UserError0() << name_ << " has incomplete MEASURE line";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : getDCSweepVal
// Purpose       : Used to get the current value of the DC sweep variable that
//                 will be used for the measurement window.
// Special Notes : For TABLE-based sweeps (for .DC data=table), it is the current
//                 row index within that table where that index starts at 1.  For
//                 all other DC sweep types, it is the current value of the first
//                 variable in the DC sweep vector.
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 06/16/2020
//-----------------------------------------------------------------------------
double getDCSweepVal(const std::vector<Analysis::SweepParam> & dcParamsVec)
{
  double sweepVal;
  if (dcParamsVec[0].type == "TABLE")
    sweepVal = (dcParamsVec[0].count%dcParamsVec[0].maxStep) +1 ;
  else
    sweepVal = dcParamsVec[0].currentVal;

 return sweepVal;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

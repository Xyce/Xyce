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
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 9/4/04
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_ANP_SweepParam.h>
#include <N_ANP_SweepParamFreeFunctions.h>
#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LOA_CktLoader.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_MachDepParams.h>
#include <N_UTL_Math.h>
#include <N_UTL_Param.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_OptionBlock.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : SweepParam::updateCurrentVal
//
// Purpose       : Updates the values of the parameters used in a sweep.
//
// Special Notes : This is very similar to the "update" function in the
//                 class Device::SweepData.  (which no longer exists).
//
// Scope         : public
// Creator       : Eric R. Keiter,SNL, Computational Sciences
// Creation Date : 10/31/2003
//-----------------------------------------------------------------------------
bool SweepParam::updateCurrentVal (int stepNumberArg)
{
  // ERK. failsafes:
  if (maxStep==0) maxStep=1;
  if (interval==0) interval=1;

  outerStepNumber  = stepNumberArg/interval;
  int inum            = outerStepNumber/maxStep;
  int localStepNumber = outerStepNumber - inum*maxStep;

  // We must keep track of whether we're at the first step of our sweep,
  // because some device manager features need to know that.
  // It is important that we only set this when localStepNumber first becoms
  // zero, not every time localStepNumber *is* zero, because an outer loop
  // might remain at 0 for quite some time.

  if (localStepNumber == 0 && localStepNumber != lastLocalStepNumber_)
  {
    sweepResetFlag_=true;
  } 
  else
  {
    sweepResetFlag_=false;
  }


  // This flag determines if the setParam function will be called for this parameter.
  // Some parameters, like TEMP are expensive to set, so don't do it unless you have to.
  bool callSetParam=true; 
  if (type == "LIN")
  {
    if(localStepNumber != lastLocalStepNumber_)
    {
      currentVal = startVal + static_cast<double>(localStepNumber)*stepVal;
      callSetParam=true;
    }
    else
    {
      callSetParam=false;
    }
    ++count;
  }
  else if (type == "DEC" || type == "OCT")
  {
    if(localStepNumber != lastLocalStepNumber_)
    {
      currentVal = startVal*pow(stepMult, static_cast<double>(localStepNumber) );
      callSetParam=true;
    }
    else
    {
      callSetParam=false;
    }
    ++count;
  }
  else if (type == "LIST")
  {
    if(localStepNumber != lastLocalStepNumber_)
    {
      int size=  valList.size();
      int index = (localStepNumber < size)?localStepNumber:(size-1);
      currentVal = valList[index];
      callSetParam=true;
    }
    else
    {
      callSetParam=false;
    }
    ++count;
  }
  else if (type == "TABLE")
  {
    currentVal = valList[stepNumberArg];
    ++count;
  }
  else if (type=="NORMAL" || type=="UNIFORM" || type=="GAMMA")
  {
    // Doesn't do anything
    // The sampling values are set elsewhere
    // This class is used for sampling as well as sweep loops, although all the machinery
    // for determing sample points is in the UQSupport functions.
    ++count; // needed?
  }
  else
  {
    Report::DevelFatal0().in("SweepParam::updateCurrentVal") << "Unsupported type " << type << " specified";
  }


  if (DEBUG_TIME)
  {
#if 0
    Xyce::dout() << std::endl
                 << Xyce::subsection_divider << std::endl
                 << "updateCurrentVal" << std::endl
                 << "  name             = " << name << std::endl
                 << "  stepNumberArg    = " << stepNumberArg<< std::endl
                 << "  interval         = " << interval  << std::endl
                 << "  outerStepNumber  = " << outerStepNumber << std::endl
                 << "  localStepNumber  = " << localStepNumber << std::endl
                 << "  inum             = " << inum      << std::endl
                 << "  sweepResetFlag   = " << sweepResetFlag_ << std::endl
                 << "  currentVal       = " << currentVal << std::endl
                 << Xyce::subsection_divider << std::endl;
#else
    Xyce::dout() << "  " << name << "  "<< outerStepNumber << "  "<< localStepNumber << "  "<< lastLocalStepNumber_ << std::endl;
#endif
  }

  lastLocalStepNumber_=localStepNumber;
  return callSetParam;
}

//-----------------------------------------------------------------------------
// Function      : SweepParam::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
std::ostream &
operator<<(std::ostream & os, const SweepParam & sp)
{
  os << "\tname            = " << sp.name
     << "\tcurrentVal      = " << sp.currentVal
     << std::endl;
  return os;
}

//-----------------------------------------------------------------------------
// Function      : parseSweepParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Oct  1 09:19:28 2014
//-----------------------------------------------------------------------------
SweepParam parseSweepParams(Xyce::Util::ParamList::const_iterator first, Xyce::Util::ParamList::const_iterator last)
{
  if (DEBUG_ANALYSIS)
  {
    Xyce::dout() << std::endl << section_divider << std::endl
                 << "parseSweepParam" << std::endl;

    for (Xyce::Util::ParamList::const_iterator it = first, end = last; it != end; ++it)
    {
      Xyce::dout() << (*it).uTag() << "\t";
      if ((*it).uTag() == "PARAM" || 
          (*it).uTag() == "TYPE" || 
          (*it).uTag() == "DATASET" )
      {
        Xyce::dout() << (*it).stringValue();
      }
      else
      {
        Xyce::dout() << (*it).getImmutableValue<double>();
      }
      Xyce::dout() << std::endl;
    }
  }

  SweepParam sweep_param;
  
  Xyce::Util::ParamList::const_iterator it_param = last;
  for (Xyce::Util::ParamList::const_iterator it = first, end = last; it != end; ++it)
  {
    if ((*it).uTag() == "TYPE")
    {
      sweep_param.type = (*it).stringValue();
    }
    else if ((*it).uTag() == "PARAM")
    {
      it_param = it;
      sweep_param.name = (*it).stringValue();
    }
    else if ((*it).uTag() == "DATASET")
    {
      sweep_param.dataSetName = (*it).stringValue();
      Xyce::Util::toUpper( sweep_param.dataSetName );
    }
  }

  if (it_param != last)
  {
    if (sweep_param.type == "LIN") // default
    {
      sweep_param.startVal = (*++it_param).getImmutableValue<double>();
      sweep_param.stopVal  = (*++it_param).getImmutableValue<double>();
      sweep_param.stepVal  = (*++it_param).getImmutableValue<double>();
    }
    else if (sweep_param.type == "DEC" || sweep_param.type == "OCT")
    {
      sweep_param.startVal = (*++it_param).getImmutableValue<double>();
      sweep_param.stopVal  = (*++it_param).getImmutableValue<double>();

      ExtendedString numStepsStr((*++it_param).stringValue());
      sweep_param.numSteps = (*it_param).getImmutableValue<int>();

      // Error checking
      if ( !numStepsStr.isInt() )
        Report::UserError() << "Points parameter on .DC or .STEP line must be an integer";
      if ( sweep_param.numSteps < 1 )
        Report::UserError() << "Points parameter on .DC or .STEP line must be >= 1";
    }
    else if (sweep_param.type == "LIST")
    {
      for (Xyce::Util::ParamList::const_iterator it = ++it_param, end = last; it != end; ++it)
      {
        sweep_param.valList.push_back((*it).getImmutableValue<double>());
      }
    }
    else if (sweep_param.type == "TABLE")
    {
      for (Xyce::Util::ParamList::const_iterator it = ++it_param, end = last; it != end; ++it)
      {
        sweep_param.valList.push_back((*it).getImmutableValue<double>());
      }
    }
    else
    {
      Report::DevelFatal().in("parseSweepParam") << "Unsupported DC type";
    }
  }

  return sweep_param;
  
}

//-----------------------------------------------------------------------------
// Function      : setupSweepLoop
// Purpose       : Processes sweep parameters.
// Special Notes : Used for DC and STEP analysis classes.
// Scope         : public
// Creator       : Eric R. Keiter, SNL.
// Creation Date : 08/21/04
//-----------------------------------------------------------------------------
int setupSweepLoop(Parallel::Machine comm, Loader::Loader &loader, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end)
{
  // loop over the param containers, and check that all the params exist.
  // (the device package will complain (e.g, cause a netlist parsing error) if
  // it can't find the param)
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    SweepParam &sweep_param = (*it);

    loader.getParamAndReduce(comm, sweep_param.name);
  }

  // now set various values in the SweepParams, if all of the params exist
  // in all of them.
  int pinterval = setSweepLoopVals(begin, end);

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return pinterval;
}

//-----------------------------------------------------------------------------
// Function      : setSweepLoopVals
// Purpose       : Initialize various values in the SweepParam containers.  
// Special Notes : Used for DC and STEP analysis classes.
//                 This is broken out into its own function, so it can be 
//                 called by remeasure without the check by the device manager
//                 (loader.getParamAndReduce()) done in setupSweepLoop().  
//                 This code used to be part of setupSweepLoop().
// Scope         : public
// Creator       : Eric R. Keiter, SNL.
// Creation Date : 08/21/04
//-----------------------------------------------------------------------------
int setSweepLoopVals(std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end)
{
  double pinterval = 1.0;
  double pcount = 0.0, pstart, pstop, pstep;

  bool table=false;
  bool data=false;

  // loop over the param containers:
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    SweepParam &sweep_param = (*it);

    // set interval:
    sweep_param.interval = static_cast<int> (pinterval);

    // This stuff should probably be moved up into the SweepParam class.
    // obtain next pinterval:
    if (sweep_param.type=="LIN")
    {
      pstart = sweep_param.startVal;
      pstop  = sweep_param.stopVal;
      pstep  = sweep_param.stepVal;
      // ----------
      // pcount = floor(((pstop - pstart)/pstep) + 1.0);
      // The computation of "pcount" above is notoriously prone to roundoff
      // error, especially if the "floor" function is an in-line function
      // and is subject to high levels of optimization on x86 processors.
      // The symptom is that the last step of a DC sweep (or other sweep)
      // gets lost for very specific combinations of start/stop/step.
      // The next few lines are an attempt to mitigate this roundoff issue,
      // which was present in Xyce for years, and was inherited from SPICE3F5,
      // from which the above expression was taken verbatim.

      // Compute the number of steps of size pstep between pstart and pstop
      pcount = floor(((pstop - pstart)/pstep));
      // Here we're checking that adding one more step doesn't pass pstop
      // by more than machine precision.  If we're within machine precision
      // of pstop by taking one more step, that must mean we undercounted
      // due to roundoff in the division --- if we hadn't undercounted, we'd
      // exceed pstop by (nearly) a full pstep.
      if ( fabs(pstop-(pstart+(pcount+1.0)*pstep)) < 2.0*Xyce::Util::MachineDependentParams::MachinePrecision())
      {
        pcount += 1.0;
      }

      // Pcount is now the exact number of steps of size pstep between pstart
      // and pstop, with roundoff handled mostly cleanly.

      // finally, because our actual loop does a loop from zero to maxStep-1,
      // we have to pad maxStep (as is done in the original pcount expression
      // above) to get the full range.
      pcount += 1.0;

      // done this way, we should no longer miss final steps of DC sweeps.
      // Fixed 31 Jul 2012.  This was bug 695 in Bugzilla, and had plagued
      // us since Xyce was first ported to Linux with GCC.
      // ----------

      sweep_param.maxStep = static_cast<int>(pcount);

      if (sweep_param.maxStep < 1)
      {
        Report::UserWarning0() << "Linear DC or STEP parameters for sweep over "
                               << sweep_param.name
                               << " would result in no steps taken.  Check sign of step value to assure stop value can be reached from start value.  Start="
                               << pstart << " Stop=" << pstop
                               << " step=" << pstep
                               << ".  Only a single computation will be done at the start value.";
        pcount=1;
        sweep_param.maxStep = 1;
      }
    }
    else if(sweep_param.type=="DEC")
    {
      double numSteps = static_cast<double>(sweep_param.numSteps);
      // stepMult could also be calculated as pow(10,(1/numSteps))
      double stepMult = exp(log(10.0)/numSteps);
      sweep_param.stepMult = stepMult;

      pstart   = sweep_param.startVal;
      pstop    = sweep_param.stopVal;
      pcount   = floor(fabs(log10(pstart) - log10(pstop)) * numSteps + 1);
      // This check is to make sure that pstart is always smaller or equal to
      // pstop.  Only throw a warning if pstart exceeds pstop by more than
      // machine precision.  The case where they're equal should already
      // have been handled by the formula for pcount.
      if (pstart > pstop + 2.0*Xyce::Util::MachineDependentParams::MachinePrecision())
      {
        pcount=1;
        sweep_param.maxStep = 1;
        Report::UserWarning0() << "Decade DC or STEP parameters for sweep over "
                               << sweep_param.name
                               << " request stop value smaller than start value.  Xyce will compute only one step at the start value.  If decade sweep over values is desired, set start value smaller than stop value.  Start=" << pstart << " Stop="<< pstop << ".";
      }
      else
      {
        sweep_param.maxStep = static_cast<int>(pcount);
      }
    }
    else if(sweep_param.type=="OCT")
    {
      double numSteps = static_cast<double>(sweep_param.numSteps);
      // stepMult could also be calculated as pow(2,1/(numSteps))
      double stepMult = exp(log(2.0)/numSteps);

      // changed to remove dependence on "log2" function, which apparently
      // doesn't exist in the math libraries of FreeBSD or the mingw
      // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
      double ln2=log(2.0);

      sweep_param.stepMult = stepMult;
      pstart   = sweep_param.startVal;
      pstop    = sweep_param.stopVal;
      pcount   = floor(fabs(log(pstart) - log(pstop))/ln2 * numSteps + 1);
      // This check is to make sure that pstart is always smaller or equal to
      // pstop.  Only throw a warning if pstart exceeds pstop by more than
      // machine precision.  The case where they're equal should already
      // have been handled by the formula for pcount.
      if (pstart > pstop + 2.0*Xyce::Util::MachineDependentParams::MachinePrecision())
      {
        pcount=1;
        sweep_param.maxStep = 1;
        Report::UserWarning0() << "Octave DC or STEP parameters for sweep over "
                               << sweep_param.name
                               << " request stop value smaller than start value.  Xyce will compute only one step at the start value.  If decade sweep over values is desired, set start value smaller than stop value.  Start=" << pstart << " Stop="<< pstop << ".";
      }
      else
      {
        sweep_param.maxStep = static_cast<int>(pcount);
      }
    }
    else if(sweep_param.type=="LIST")
    {
      pcount = sweep_param.valList.size();
      sweep_param.maxStep = sweep_param.valList.size();
    }
    else if (sweep_param.type=="TABLE")
    {
      pcount = sweep_param.valList.size();
      sweep_param.maxStep = sweep_param.valList.size();
      table=true;
    }
    else if (sweep_param.type=="DATA")
    {
      pcount = 1;
      sweep_param.maxStep = 1;
      data=true;
    }
    else if (sweep_param.type=="NORMAL" || sweep_param.type=="UNIFORM" || sweep_param.type=="GAMMA")
    {
      // ERK. this is a no-op, to avoid the error message.  This class is used for
      // sampling as well as sweep loops.
    }
    else
    {
      Report::UserError0() << " Unsupported STEP type";
    }
    pinterval *= pcount;
  }

  if(table) // this is a bit of a quick hack for now
  {
    pinterval = pcount;

    // double check that the valLists for each param are same size, 
    // as for TABLE, all params change simultaneously
    for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
    {
      SweepParam &sweep_param = (*it);

      if(pcount != sweep_param.valList.size())
      {
        Report::UserError0() << "For the TABLE all parameters must have matching number of points.";
      }
    }
  }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return static_cast<int>(pinterval);
}

//-----------------------------------------------------------------------------
// Function      : updateSweepParams
// Purpose       : Update parameters either for DC or STEP sweeps
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/7/18
//-----------------------------------------------------------------------------
bool updateSweepParams(int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end)
{
  bool reset = false;

  // set parameter(s)
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    (*it).updateCurrentVal(step_count);
    reset = reset || (*it).getSweepResetFlag();
  }

  return reset;
}

//-----------------------------------------------------------------------------
// Function      : updateSweepParams
// Purpose       : Update parameters either for DC or STEP sweeps
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/26/04
//-----------------------------------------------------------------------------
bool updateSweepParams(Loader::Loader &loader, int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end, bool overrideOriginal)
{
  bool reset = false;

  // set parameter(s)
  for (std::vector<SweepParam>::iterator it = begin; it != end; ++it)
  {
    bool setParamFlag = (*it).updateCurrentVal(step_count);
    reset = reset || (*it).getSweepResetFlag();
    if (setParamFlag)
    {
      loader.setParam((*it).name, (*it).currentVal, overrideOriginal);
    }
    if (DEBUG_ANALYSIS)
    {
      Xyce::dout() << "Updating parameter " << (*it).name << " to " << (*it).currentVal << std::endl;
    }
  }

  // ERK. this call is necessary to ensure sweeps and homotopy play well together.
  loader.updateDependentParams();

  return reset;
}

//-----------------------------------------------------------------------------
// Function      : processDataStatements
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool processDataStatements(
    const Xyce::Util::OptionBlock & paramsBlock,
    std::map< std::string, std::vector<std::string> > & dataNamesMap,
    std::map< std::string, std::vector< std::vector<double> > > & dataTablesMap
    )
{
  // transpose the data statement numbers.
  //
  // Data statement-based .STEP is equivalent to the TABLE type.  However,
  // the data is expressed in a single line (with carriage returns, a 2D table)
  // with each column containing the values of the samples for a single param.

   std::string dataSetName;
   std::vector< std::string > params;
   std::vector<double> compressedRowData;

   Xyce::Util::ParamList::const_iterator iter = paramsBlock.begin();
   Xyce::Util::ParamList::const_iterator end  = paramsBlock.end();

   for(;iter!=end;++iter) 
   { 
     std::string tag = iter->tag();
     Xyce::Util::toUpper(tag);
     if (tag=="PARAM")
     {
       params.push_back(iter->stringValue());
     }
     else if (tag=="NAME")
     {
       dataSetName=iter->stringValue();
     }
     else if(tag=="VAL")
     {
       double tmp = iter->getImmutableValue<double>();
       compressedRowData.push_back(tmp);
     }
     else
     {
       Report::UserError0() << ".DATA line not formatted correctly.";
       return false;
     }
   }

   // check sizes:
   int paramSize = params.size();
   int CRdataSize = compressedRowData.size();

   // Must have both parameters and data rows.  The total number of data entries
   // must be an integer multiple of the number of parameters.
   if ( (paramSize == 0) || (CRdataSize == 0) || (CRdataSize%paramSize != 0) )
   {
     Report::UserError0() << ".DATA line " << dataSetName << " not formatted correctly.";
     return false;
   }

   int samples = CRdataSize / paramSize;
   int crIndex=0;
   std::vector< std::vector<double> > twoDimData(samples);
   for (int is=0;is<samples;++is)
   {
     twoDimData[is].resize(params.size(),0.0);
     for (int ip=0;ip<paramSize;++ip)
     {
       twoDimData[is][ip] = compressedRowData[crIndex++];
     }
   }

   dataNamesMap [ dataSetName ] = params;
   dataTablesMap[ dataSetName ] = twoDimData;

   return true;
}

//-----------------------------------------------------------------------------
// Function      : convertData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/5/18
//-----------------------------------------------------------------------------
bool convertData(
    SweepVector & stepSweepVector,
    const std::map< std::string, std::vector<std::string> > & dataNamesMap,
    const std::map< std::string, std::vector< std::vector<double> > > & dataTablesMap
    )
{
  if (dataTablesMap.empty() || dataNamesMap.empty())
  {
    Report::UserError0() << "Invalid sweep parameter name.  Netlist may lack any valid .DATA statements";
    return false;
  }

  std::vector<std::string> usedDataList;

  for (int is=0;is<stepSweepVector.size();++is)
  {
    if ( stepSweepVector[is].type == "DATA")
    {
      usedDataList.push_back(stepSweepVector[is].dataSetName);
    }
  }

  if (!(usedDataList.empty()))
  {
    stepSweepVector.clear();
  }

  // each column of the 2D data array is for a separate parameter
  // So, each column will create a "TABLE" style sweep parameter.
  for (int in=0;in<usedDataList.size();++in)
  {
    SweepParam  tmpSweepParam;

    const std::string dataSetName = usedDataList[in];
    const std::map< std::string, std::vector< std::vector<double> > >::const_iterator iterData = dataTablesMap.find( dataSetName );
    const std::map< std::string, std::vector<std::string> >::const_iterator iterParam = dataNamesMap.find (dataSetName);

    if ( iterData != dataTablesMap.end() && iterParam != dataNamesMap.end() )
    {
      const std::vector< std::vector<double> > & tmpData = iterData->second;
      const std::vector< std::string > & params = iterParam->second;

      int paramSize = params.size();
      int numSamples = tmpData.size();

      for (int ip=0;ip<paramSize;++ip)
      {
        tmpSweepParam.name=params[ip];
        tmpSweepParam.type="TABLE";
        tmpSweepParam.valList.clear();

        for (int is=0;is<numSamples;++is)
        {
          tmpSweepParam.valList.push_back( tmpData[is][ip] );
        }

        stepSweepVector.push_back(tmpSweepParam);
      }
    }
    else
    {
      Report::UserError0() << "Invalid table name " << dataSetName << " from .DATA line used as sweep variable";
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : isDataSpecified
// Purpose       : determine if a DATA=<val> was used on analysis line,
//                 like .AC, .DC, .NOISE or .STEP
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 6/6/19
//-----------------------------------------------------------------------------
bool isDataSpecified(const Xyce::Util::OptionBlock & paramsBlock)
{
  for (Xyce::Util::ParamList::const_iterator it = paramsBlock.begin(), end = paramsBlock.end(); it != end; ++it)
  {
    std::string tag = (*it).uTag();
    std::string val = (*it).stringValue();
    Xyce::Util::toUpper(tag);
    Xyce::Util::toUpper(val);
    if (tag == "TYPE" && val == "DATA")
    {
      return true;
    }
  }

  return false;
}

} // namespace Analysis

//-----------------------------------------------------------------------------
// Function      : SweepParam::packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/13/2020
//-----------------------------------------------------------------------------
template<>
int
Pack<Analysis::SweepParam>::packedByteCount(const Analysis::SweepParam &param)
{
  int byteCount = 0;

  //std::string name;
  byteCount += param.name.length() + sizeof(int); 
  //std::string opName;
  byteCount += param.opName.length() + sizeof(int); 
  //std::string baseName; 
  byteCount += param.baseName.length() + sizeof(int); 
  //std::string type;
  byteCount += param.type.length() + sizeof(int); 

  //double startVal;
  byteCount += sizeof(double);
  //double stopVal;
  byteCount += sizeof(double);
  //double stepVal;
  byteCount += sizeof(double);
  //double stepMult;
  byteCount += sizeof(double);

  // for normal distribution
  //double mean;
  byteCount += sizeof(double);
  //double stdDev;
  byteCount += sizeof(double);

  // for gamma distribution
  //double alpha;
  byteCount += sizeof(double);
  //double beta;
  byteCount += sizeof(double);

  //double upper_bound;
  byteCount += sizeof(double);
  //double lower_bound;
  byteCount += sizeof(double);

  //bool upper_boundGiven;
  byteCount += sizeof(int);
  //bool lower_boundGiven;
  byteCount += sizeof(int);

  //double currentVal;
  byteCount += sizeof(double);

  //int numSteps;
  byteCount += sizeof(int);

  //int count;
  byteCount += sizeof(int);
  //int maxStep;
  byteCount += sizeof(int);
  //int interval;
  byteCount += sizeof(int);

  //int outerStepNumber;
  byteCount += sizeof(int);

  //std::vector<double> valList;
  byteCount += sizeof(double) * param.valList.size();

  //std::string dataSetName;
  byteCount += param.dataSetName.length() + sizeof(int); 

  //int astOpIndex;
  byteCount += sizeof(int);
  //enum Util::astRandTypes astType;
  byteCount += sizeof(int);

  //bool sweepResetFlag_;
  byteCount += sizeof(int);
  //int lastLocalStepNumber_;
  byteCount += sizeof(int);

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : SweepParam::pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/13/2020
//-----------------------------------------------------------------------------
template<>
void
Pack<Analysis::SweepParam>::pack(const Analysis::SweepParam &param, char * buf, int bsize, int & pos, Parallel::Communicator* comm )
{
  int length;
  std::string tmp;
#if 0
  //pack tag
  length = param.tag_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( param.tag_.c_str(), length, buf, bsize, pos );

  //pack type
  int enum_type = param.data_->enumType();

  comm->pack( &enum_type, 1, buf, bsize, pos );

  //pack value
  switch( param.data_->enumType() )
  {
    case -1:
      break;

    case Util::STR:
      length = param.getValue<std::string>().length();
      comm->pack( &length, 1, buf, bsize, pos );
      comm->pack( param.getValue<std::string>().c_str(), length, buf, bsize, pos );
      break;

    case Util::DBLE:
      comm->pack( &(param.getValue<double>()), 1, buf, bsize, pos );
      break;

    case Util::INT:
      comm->pack( &(param.getValue<int>()), 1, buf, bsize, pos );
      break;

    case Util::BOOL:
    {
      int i;
      if (param.getValue<bool>())
        i = 1;
      else
        i = 0;
      comm->pack( &i, 1, buf, bsize, pos );
    }
      break;

    case Util::LNG:
      comm->pack( &(param.getValue<long>()), 1, buf, bsize, pos );
      break;

    case Util::EXPR:
      tmp = param.getValue<Util::Expression>().get_expression();
      length = tmp.length();
      comm->pack( &length, 1, buf, bsize, pos );
      comm->pack( tmp.c_str(), length, buf, bsize, pos );
      break;

    case Util::STR_VEC:
    {
      const std::vector<std::string> &string_vector = param.getValue<std::vector<std::string> >();
      length = (int) string_vector.size();
      comm->pack( &length, 1, buf, bsize, pos );
      for (int i=0; i < (int) string_vector.size(); i++)
      {
        length = string_vector[i].length();
        comm->pack( &length, 1, buf, bsize, pos );
        comm->pack( string_vector[i].c_str(), length, buf, bsize, pos );
      }
    }

      break;

    case Util::DBLE_VEC:
    {
      const std::vector<double> &double_vector = param.getValue<std::vector<double> >();
      length = (int) double_vector.size();
      comm->pack( &length, 1, buf, bsize, pos );
      comm->pack( &double_vector[0], length, buf, bsize, pos );
    }

      break;

    default:   Report::DevelFatal() << "SweepParam::pack: unknown type " << param.data_->enumType();
  }
#else
  //pack name
  length = param.name.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( param.name.c_str(), length, buf, bsize, pos );

  //pack opName
  length = param.opName.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( param.opName.c_str(), length, buf, bsize, pos );

  //pack baseName
  length = param.baseName.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( param.baseName.c_str(), length, buf, bsize, pos );

  //pack type
  length = param.type.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( param.type.c_str(), length, buf, bsize, pos );

  // pack startVal
  comm->pack( &(param.startVal), 1, buf, bsize, pos );
  // pack stopVal
  comm->pack( &(param.stopVal), 1, buf, bsize, pos );
  // pack stepVal
  comm->pack( &(param.stepVal), 1, buf, bsize, pos );
  // pack stepMult
  comm->pack( &(param.stepMult), 1, buf, bsize, pos );
  // pack mean
  comm->pack( &(param.mean), 1, buf, bsize, pos );
  // pack stdDev
  comm->pack( &(param.stdDev), 1, buf, bsize, pos );
  // pack alpha
  comm->pack( &(param.alpha), 1, buf, bsize, pos );
  // pack beta
  comm->pack( &(param.beta), 1, buf, bsize, pos );
  // pack upper_bound
  comm->pack( &(param.upper_bound), 1, buf, bsize, pos );
  // pack lower_bound
  comm->pack( &(param.lower_bound), 1, buf, bsize, pos );

  // pack upper_boundGiven
  {
    int i;
    if (param.upper_boundGiven){ i = 1; }
    else { i = 0; }
    comm->pack( &i, 1, buf, bsize, pos );
  }

  // pack lower_boundGiven
  {
    int i;
    if (param.lower_boundGiven) { i = 1; }
    else { i = 0; }
    comm->pack( &i, 1, buf, bsize, pos );
  }

  // pack currentVal
  comm->pack( &(param.currentVal), 1, buf, bsize, pos );

  // pack numSteps
  comm->pack( &(param.numSteps), 1, buf, bsize, pos );
  // pack count
  comm->pack( &(param.count), 1, buf, bsize, pos );
  // pack maxStep
  comm->pack( &(param.maxStep), 1, buf, bsize, pos );
  // pack interval
  comm->pack( &(param.interval), 1, buf, bsize, pos );
  // pack outerStepNumber
  comm->pack( &(param.outerStepNumber), 1, buf, bsize, pos );

  // pack valList
  {
    const std::vector<double> &double_vector = param.valList;
    length = (int) double_vector.size();
    comm->pack( &length, 1, buf, bsize, pos );
    comm->pack( &double_vector[0], length, buf, bsize, pos );
  }

  //pack dataSetName
  length = param.dataSetName.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( param.dataSetName.c_str(), length, buf, bsize, pos );

  // pack astOpIndex
  comm->pack( &(param.astOpIndex), 1, buf, bsize, pos );

  // pack astType
  {
    int i = static_cast<int> (param.astType);
    comm->pack( &i, 1, buf, bsize, pos );
  }

  // pack sweepResetFlag_
  {
    int i;
    if (param.sweepResetFlag_){ i = 1;}
    else{ i = 0;}
    comm->pack( &i, 1, buf, bsize, pos );
  }

  // pack lastLocalStepNumber_
  comm->pack( &(param.lastLocalStepNumber_), 1, buf, bsize, pos );
#endif
}

//-----------------------------------------------------------------------------
// Function      : SweepParam::unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/13/2020
//-----------------------------------------------------------------------------
template<>
void
Pack<Analysis::SweepParam>::unpack(Analysis::SweepParam &param, char * pB, int bsize, int & pos, Parallel::Communicator* comm )
{
  int length = 0;
  int vector_size = 0;

#if 0
  //unpack tag
  comm->unpack( pB, bsize, pos, &length, 1 );

  param.tag_ = std::string( (pB+pos), length );
  pos += length;

  //unpack type
  int enum_type = -1;
  comm->unpack( pB, bsize, pos, &enum_type, 1 );

  switch (enum_type)
  {
    case -1:
      break;

    case Util::STR:
      comm->unpack( pB, bsize, pos, &length, 1 );
      param.setVal(std::string( (pB+pos), length ));
      pos += length;
      break;

    case Util::DBLE: 
    {
      double d = 0.0;
      
      comm->unpack( pB, bsize, pos, &d, 1 );
      param.setVal(d);
    }
    break;

    case Util::INT: 
    {
      int i = 0;
      comm->unpack( pB, bsize, pos, &i, 1 );
      param.setVal(i);
    }
    break;

    case Util::BOOL:
    {
      int i = 0;
      comm->unpack( pB, bsize, pos, &i, 1 );
      if (i == 0)
        param.setVal(false);
      else
        param.setVal(true);
    }
    break;

    case Util::LNG:
    {
      long l = 0;
      comm->unpack( pB, bsize, pos, &l, 1 );
      param.setVal(l);
    }
    break;

    case Util::EXPR:
      comm->unpack( pB, bsize, pos, &length, 1 );
      // ERK.  This expression is allocated iwth the base group, which is easy 
      // to create, but doesn't do anything.  But, this doesn't matter.  .param 
      // and .global_param are not usually evaluated stand-alone, as they are 
      // usually attached to other expressions. The expressions they are 
      // attached to are the ones that need a meaningful group class.
      param.setVal(Util::Expression(Teuchos::rcp(new Xyce::Util::baseExpressionGroup()), std::string( (pB+pos), length )));
      pos += length;
      break;

    case Util::STR_VEC:
    {

      comm->unpack( pB, bsize, pos, &vector_size, 1 );
      param.setVal(std::vector<std::string>());
      std::vector<std::string> &x = param.getValue<std::vector<std::string> >();
      x.reserve(vector_size);

      for (int i=0; i< vector_size; i++)
      {
        comm->unpack( pB, bsize, pos, &length, 1 );
        x.push_back(std::string( (pB+pos), length ));
        pos += length;
      }
    }

    break;

    case Util::DBLE_VEC:
    {
      comm->unpack( pB, bsize, pos, &vector_size, 1 );
      param.setVal(std::vector<double>());
      std::vector<double> &x = param.getValue<std::vector<double> >();
      x.resize(vector_size, 0.0);
      comm->unpack( pB, bsize, pos, &(x[0]), vector_size );
    }

    break;

    default:
      Report::UserFatal() << "SweepParam::unpack: unknown type";

  }
#else

  //unpack name 
  comm->unpack( pB, bsize, pos, &length, 1 );
  param.name = std::string( (pB+pos), length );
  pos += length;

  //unpack opName 
  comm->unpack( pB, bsize, pos, &length, 1 );
  param.opName = std::string( (pB+pos), length );
  pos += length;

  //unpack baseName 
  comm->unpack( pB, bsize, pos, &length, 1 );
  param.baseName = std::string( (pB+pos), length );
  pos += length;

  //unpack type 
  comm->unpack( pB, bsize, pos, &length, 1 );
  param.type = std::string( (pB+pos), length );
  pos += length;

  // unpack startVal
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.startVal = d;
  }
  // unpack stopVal
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.stopVal = d;
  }
  // unpack stepVal
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.stepVal = d;
  }
  // unpack stepMult
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.stepMult = d;
  }
  // unpack mean
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.mean = d;
  }
  // unpack stdDev
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.stdDev = d;
  }
  // unpack alpha
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.alpha = d;
  }
  // unpack upper_bound
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.upper_bound = d;
  }
  // unpack lower_bound
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.lower_bound = d;
  }
  // unpack upper_boundGiven
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    if (i == 0) { param.upper_boundGiven = false; }
    else { param.upper_boundGiven = true; }
  }
  // unpack lower_boundGiven
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    if (i == 0) { param.lower_boundGiven = false; }
    else { param.lower_boundGiven = true; }
  }
  // unpack currentVal
  {
    double d = 0.0;
    comm->unpack( pB, bsize, pos, &d, 1 );
    param.currentVal = d;
  }
  // unpack numSteps
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.numSteps = i;
  }
  // unpack count
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.count = i;
  }
  // unpack maxStep
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.maxStep = i;
  }
  // unpack interval
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.interval = i;
  }
  // unpack outerStepNumber
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.outerStepNumber = i;
  }
  // unpack valList
  {
    comm->unpack( pB, bsize, pos, &vector_size, 1 );
    std::vector<double> &x = param.valList;
    x.resize(vector_size, 0.0);
    comm->unpack( pB, bsize, pos, &(x[0]), vector_size );
  }
  //unpack dataSetName 
  comm->unpack( pB, bsize, pos, &length, 1 );
  param.dataSetName = std::string( (pB+pos), length );
  pos += length;
  // unpack astOpIndex
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.astOpIndex = i;
  }
  // unpack astType
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.astType = i;
  }
  // unpack sweepResetFlag_
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    if (i == 0) { param.sweepResetFlag_ = false; }
    else { param.sweepResetFlag_ = true; }
  }
  // unpack lastLocalStepNumber_
  {
    int i = 0;
    comm->unpack( pB, bsize, pos, &i, 1 );
    param.lastLocalStepNumber_ = i;
  }
#endif
}

} // namespace Xyce


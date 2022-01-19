//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterTimeProbe.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>
#include <N_ANP_AnalysisManager.h>

namespace Xyce {
namespace IO {
namespace Outputter {


//-----------------------------------------------------------------------------
// Class         : TimeProbe
// Purpose       : Outputter class for transient runs, Probe output
//                 format (PSpice-compatibility output)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeProbe::TimeProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeProbe::TimeProbe(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    printCount_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0),
    printStepHeader_(false),
    opList_(),
    opInitialTime_(output_manager.getOpBuilderManager().createOp("ANALYSIS_INITIAL_TIME")),
    opFinalTime_(output_manager.getOpBuilderManager().createOp("ANALYSIS_FINAL_TIME"))
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csd";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::~TimeProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeProbe::~TimeProbe()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
  delete opInitialTime_;
  delete opFinalTime_;
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::timeHeader(Parallel::Machine comm)
{
  if (os_)
  {
    std::ostream &os = *os_;

    // count the number of output variables.
    printCount_ = 0;
    for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
    {
      if ((*it)->id() != Util::Op::identifier<Util::Op::UndefinedOp>())
        ++printCount_;
    }

    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << Util::Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetlistFilename() << "'" << std::endl;

    os.setf(std::ios::scientific);
    os.precision(1); // streamPrecision_);
    if (outputManager_.getStepSweepVector().empty())
    {
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      os << "SUBTITLE='Step param";
      for (std::vector<Analysis::SweepParam>::const_iterator it = outputManager_.getStepSweepVector().begin(); it != outputManager_.getStepSweepVector().end(); ++it)
      {
        os << " " << it->name << " = " << it->currentVal;
      }
    }

    os << "'" << std::endl;

    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
    os << "'" << std::endl;

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      os << "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
    }
    else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC)
    {
      // handle .PRINT AC_IC lines   
      os << "ANALYSIS='AC Sweep' SERIALNO='12345'" <<  std::endl;
    }
    else
    {
      os << "ANALYSIS='DC Sweep' " << "SERIALNO='12345'" <<  std::endl;
    }

    os << "ALLVALUES='NO' COMPLEXVALUES='NO' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    // handle .PRINT TRAN and .PRINT AC_IC lines
    if ( (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT) ||
         (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC) )
    {
      os << "SWEEPVAR='Time' SWEEPMODE='VAR_STEP'" <<
        std::endl;
    }
    else
    {
      os << "SWEEPVAR='";
      os << outputManager_.getPRINTDCname();
      os << "' SWEEPMODE=";
      if ((outputManager_.getDCSweepVector().size() > 0) && (outputManager_.getDCSweepVector()[0].type == "LIST"))
      {
        os << "'LIST'" << std::endl;
      }
      else
      {
        os << "'VAR_STEP'" << std::endl;
      }
    }

    // handle .PRINT TRAN and .PRINT AC_IC lines
    if ( (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT) ||
         (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC) )
    {
      double initial_time = getValue(comm, *(opInitialTime_), Util::Op::OpData()).real();
      double final_time = getValue(comm, *(opFinalTime_), Util::Op::OpData()).real();
      os << "XBEGIN='" << initial_time << "' XEND='" << final_time << "'" << std::endl;
    }
    else
    {
      os << "XBEGIN='" << outputManager_.getPRINTDCstart()
         << "' XEND='" << outputManager_.getPRINTDCstop() << "'" << std::endl;
    }

    os << "FORMAT='0 VOLTSorAMPS;EFLOAT : "
       << "NODEorBRANCH;NODE  '  " << std::endl;
    os << "DGTLDATA='NO'";

    int dcSize = outputManager_.getDCSweepVector().size();
    if (dcSize > 1)
    {
      os << "  ";
      for (int idc=1;idc<dcSize;++idc)
      {
        os << "SWEEP" << idc+1 << "PARM='";
        os << outputManager_.getDCSweepVector()[idc].name;
        os << "' ";
        os << "SWEEP" << idc+1 << "VALUE='";
        os.setf(std::ios::scientific);
        os.precision(printParameters_.streamPrecision_);
        os << outputManager_.getDCSweepVector()[idc].currentVal;
        os << "' ";
        os << std::endl;
      }
    }
    else
    {
      os << std::endl;
    }

    // The PSpice A/D waveform viewer apparently has an undocumented limit on the number
    // of characters that can be on a single line in the #N section of the header.  It
    // seemed to be about 430 characters in PSpice 16.6.  Because of that, we limit the 
    // number of characters on a single line in #N section of the Xyce .csd file to maxLineLength.  
    // Also, for readability, we won't put more than maxVarsPerLine variable names on 
    // a single line.
    os << "#N" << std::endl;

    int maxLineLength = 128;
    int maxVarsPerLine = 12;
    int i = 0;  // used to count # of variables in the outString below
    std::string outString, currentVarString;

    for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
    {
      // pad the currentVarString, so that V(1) is 'V(1)' followed by a space.
      // all subsequent operations then use that padded string.
      currentVarString = "'" + (*it)->getName() +"' ";

      if ( (outString.size() + currentVarString.size()) >= maxLineLength )
      {
        // handle cases where the currentVarString is too big
        if (currentVarString.size() >= maxLineLength)
	{
          if (outString.empty())
	  {
            // outString is empty, so just output the currentVarString.
            os << currentVarString << std::endl;
            i=0; // reset i and outString to 0 and "empty"
            outString="";
          }
          else
	  {
            // outString is not empty.  So output both strings, separated by a 
            // line break.
            os << outString << std::endl << currentVarString << std::endl;
            i=0; // reset i and outString to 0 and "empty"
            outString="";
          }
        }
        else
	{
          // output existing outString, and start a new outString containing
          // only the currentVarString.
          os << outString << std::endl;
          outString = currentVarString;
          i = 1;
        }  
      }
      else
      {
        // add curentVarString to outString, and increment i (# of padded variable
        // names in outString)
        outString += currentVarString;
        i++;
      }

      // for readability, also limit the number of variables on a given line to maxVarsPerLine.
      if ( i >= maxVarsPerLine )
      { 
        os << outString << std::endl;
        i = 0; // reset i and outString to 0 and "empty"
        outString = "";
      }
    }

    // Done, iterating through the variable names. So, output the last line, 
    // if it's not empty
    if (!outString.empty())
    {
      os << outString << std::endl;
    }

    os.flush();

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(std::ios::left, std::ios::adjustfield);

  }
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeProbe::doOutputTime(
  Parallel::Machine     comm,
  const Linear::Vector &  solnVecPtr,
  const Linear::Vector &  stateVecPtr,
  const Linear::Vector &  storeVecPtr,
  const Linear::Vector &  lead_current_vector,
  const Linear::Vector &  junction_voltage_vector)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    // the IF opens up the output file and prints out a header block if .STEP is not being used or
    // if this is step 0.
    outFilename_ = outputFilename(printParameters_.filename_, 
                                  printParameters_.defaultExtension_,
                                  printParameters_.suffix_+outputManager_.getFilenameSuffix(), 
                                  outputManager_.getNetlistFilename(),
                                  printParameters_.overrideRawFilename_,
                                  printParameters_.formatSupportsOverrideRaw_,
                                  printParameters_.dashoFilename_,
                                  printParameters_.fallback_);
    os_ = outputManager_.openFile(outFilename_);
    timeHeader(comm);
  }
  else if (Parallel::rank(comm) == 0 && printStepHeader_)
  {
    // the ELSE IF prints out a header block at the start of steps 1,2,3, ...
    // this flag is set in doStartStep().  This approach did not
    // require changing the signature of doStartStep() for all of
    // the other outputters.
    (*os_) << "#;" << std::endl;
    timeHeader(comm);
    printStepHeader_ = false;
  }

  // Each time entry begins with the line:
  // #C circuitTime printCount 
  // 
  // where print count is the number of values output at each time step.  For 
  // readability, the data for each  time step will then be broken up into rows 
  // with at most maxColsPerRow columns per row.  An example is as follows.  
  // Each data point is in the form result:colIdx
  //
  // #C 0.000e+00 5
  // 0.000e+00:1   0.000e+00:2   0.000e+00:3   0.000e+00:4
  // 0.000e+00:5   
  // #C 1.000e-06 5
  // 6.283e-04:1   5.236e-04:2   4.189e-04:3   3.142e-04:4
  // 2.094e-04:5  
  //
  if (os_) 
  {
    if ( (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)  ||
         (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC) )
    {
      *os_ << "#C " << outputManager_.getCircuitTime() << " " << printCount_ << std::endl;
    }
    else
    {
      *os_ << "#C " << outputManager_.getPRINTDCvalue() << " " << printCount_ << std::endl;
    }
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0, &lead_current_vector, 0, &junction_voltage_vector), result_list);

  int colIdx = 0;  // used as a column counter local to this function 
  int maxColsPerRow = 4; // max number of columns output per row for each time value
  for (int i = 0; i < result_list.size(); ++i)
  {
    colIdx++;
    if (os_)
      *os_ << result_list[i].real() << ":" << colIdx << (colIdx%maxColsPerRow == 0 ? "\n" : "   ");
  }

  // add a new line if the last outputted row had less than maxColsPerRow entries on it
  if (os_ && colIdx%maxColsPerRow != 0)
    *os_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      (*os_) << "#;" << std::endl;
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeProbe::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
  if ( current_step > 0)
  {
    // need to print out a header block at the start of steps 1,2,3, ...
    // this flag will be reset in doOutputTime()
    printStepHeader_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeProbe::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doSteppingComplete()
{
  if (os_)
  {
    (*os_) << "#;" << std::endl;
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

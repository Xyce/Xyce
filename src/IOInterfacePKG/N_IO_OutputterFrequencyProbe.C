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

#include <N_IO_OutputterFrequencyProbe.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {
//-----------------------------------------------------------------------------
// Class         : FrequencyProbe
// Purpose       : Outputter class for frequency-domain runs, Probe output
//                 format (PSpice-compatibility output)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::FrequencyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyProbe::FrequencyProbe(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    printCount_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0),
    printStepHeader_(false)

{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csd";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::~FrequencyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyProbe::~FrequencyProbe()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::frequencyHeader(Parallel::Machine comm,
                                     double fStart,
                                     double fStop)
{
  if (os_) {
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

    os << "ANALYSIS='AC Sweep' " <<
      "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='YES' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    os << "SWEEPVAR='";
    std::string varName = outputManager_.getPRINTDCname();
    if( varName == "" )
    {
      varName="FREQ";
    }
    os << varName;
    os << "' SWEEPMODE=";
    if ( (!outputManager_.getDCSweepVector().empty()) && (outputManager_.getDCSweepVector()[0].type == "LIST") )
    {
      os << "'LIST'" << std::endl;
    }
    else
    {
      os << "'VAR_STEP'" << std::endl;
    }


    os << "XBEGIN='" << fStart << 
          "' XEND='" << fStop << "'" << std::endl;

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
// Function      : FrequencyProbe::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyProbe::doOutputFrequency(
  Parallel::Machine     comm,
  double                frequency,
  double                fStart,
  double                fStop,
  const Linear::Vector &  real_solution_vector,
  const Linear::Vector &  imaginary_solution_vector,
  const Util::Op::RFparamsData & RFparams)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_, 
                                  printParameters_.defaultExtension_,
                                  printParameters_.suffix_+outputManager_.getFilenameSuffix(), 
                                  outputManager_.getNetlistFilename(),
                                  printParameters_.overrideRawFilename_,
                                  printParameters_.formatSupportsOverrideRaw_,
                                  printParameters_.dashoFilename_,
                                  printParameters_.fallback_);
    os_ = outputManager_.openFile(outFilename_);
    frequencyHeader(comm, fStart, fStop);
  }
  else if (Parallel::rank(comm) == 0 && printStepHeader_)
  {
    // This prints out the step footer (#;) between steps 0 and 1, steps 1 and 2, etc.  
    // It also prints out the header block at the start of steps 1,2,3, ...
    // The printStepHeader_ flag is set in doStartStep().  This approach did
    // not require changing the signature of doStartStep() for all of
    // the other outputters.
    (*os_) << "#;" << std::endl;
    frequencyHeader(comm, fStart, fStop);
    printStepHeader_ = false;
  }

  // Each frequency entry begins with the line:
  // #C frequency printCount 
  // 
  // where print count is the number of values output at each time step.  For 
  // readability, the data for each frequency step will then be broken up into rows 
  // with at most maxColsPerRow columns per row.  An example is as follows.  
  // Each data point is in the form result:colIdx
  //
  // #N 'V(1)' 'V(2)' 'V(3)' 'V(4)' 'V(5)' 
  // #C 1.000e+02 5
  // 1.000e+00/0.000e+00:1   7.896e-06/2.513e-03:2   5.922e-06/1.885e-03:3   3.948e-06/1.257e-03:4
  // 1.974e-06/6.283e-04:5   
  // #C 1.778e+02 5
  // 1.000e+00/0.000e+00:1   2.497e-05/4.469e-03:2   1.873e-05/3.352e-03:3   1.248e-05/2.235e-03:4
  // 6.242e-06/1.117e-03:5  
  //   
  if (os_)
    *os_ << "#C " << frequency << " " << printCount_ << std::endl;

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(0, &real_solution_vector, &imaginary_solution_vector, 0, 0, 0), result_list);

  // used as a column-counter local to this function 
  int colIdx = 0;
  // max number of columns output per row for each frequency value.  Note: If this value
  // is changed then the regression tests for SON Bug 646 may have to change. 
  int maxColsPerRow = 4; 
  for (int i = 0; i < result_list.size(); ++i)
  {
    colIdx++;
    if (os_)
      *os_ << result_list[i].real() << "/" << result_list[i].imag() << ":" << colIdx << (colIdx%maxColsPerRow == 0 ? "\n" : "   ");
  }

  // add a new line if the last outputted row had less than maxColsPerRow entries on it
  if (os_ && colIdx%maxColsPerRow != 0)
    *os_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doFinishOutput
// Purpose       : Close the stream if there is no .STEP loop.  This function 
//               : is also called after each step, if there is a .STEP loop, 
//               : but currently does nothing in that case.  The footers 
//               : between steps 0 and 1, 1 and 2, ... are added in the 
//               : function FrequencyProbe::doOutputFrequency().  This is done
//               : because this function is actually called twice, if .OP or
//               : .STEP is used in the netlist with .AC.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doFinishOutput()
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
// Function      : FrequencyProbe::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyProbe::doStartStep(
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
// Function      : FrequencyProbe::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyProbe::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doSteppingComplete
// Purpose       : Output footer (#;) and close the stream  when a .STEP loop 
//               : is used. This footer comes after the last step.  The
//               : footers between steps 0 and 1, 1 and 2, ... are added in
//               : the function FrequencyProbe::doOutputFrequency().
// Special Notes : 
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doSteppingComplete()
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

//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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

#include <N_IO_OutputterHomotopyProbe.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : HomotopyProbe
// Purpose       : Outputter class for homotopy output, Probe output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::HomotopyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyProbe::HomotopyProbe(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    currentStep_(0),
    numberOfSteps_(0),
    index_(0),
    printCount_(0),
    printStepHeader_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".HOMOTOPY.csd";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::~HomotopyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyProbe::~HomotopyProbe()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::homotopyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::homotopyHeader(
  Parallel::Machine comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const Linear::Vector &                  solution_vector)
{
  if (os_)
  {
    std::ostream &os = *os_;

    index_ = 0;

    // first homotopy parameter name will be used as the independent variable
    // (x-axis) in the file.  So, the expression for printCount_ has -1 in it.
    printCount_ = parameter_names.size() + opList_.size() - 1;

    // make a combined name list of the continuation parameter names and the output
    // variable names.  This allows the code for printing the variable list in the 
    // #N block to be used once, and to be basically be the same as for the .TRAN and .AC
    // Probe output code.  The first variable in the list is also needed for the
    // SWEEPVAR variable
    std::vector<std::string> nameList;
    std::string sweepVarName; // used to populate SWEEPVAR variable in header block
   
    // get the continutation parameter names:
    for (std::vector<std::string>::const_iterator iter_name = parameter_names.begin(); iter_name != parameter_names.end(); ++iter_name)
    {
      if ( iter_name == parameter_names.begin() ) 
      { 
        // first homotopy variable will be used as the independent variable.
        // It is used as the SWEEPVAR below.
        sweepVarName = *iter_name; 
      }
      else
      {
        // the rest of the homotopy variables (excluding the first one) will 
        // appear in the #N list.
        nameList.push_back(*iter_name);
      }
    }      

    // get the output variable names
    for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
    {
      nameList.push_back((*it)->getName());
    }

    // start making the header block
    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << Util::Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetlistFilename() << "'" << std::endl;

    os.setf(std::ios::scientific);
    os.precision(1); // streamPrecision_);

    // Make the SUBTITLE line 
    if (outputManager_.getStepSweepVector().empty())
    {
      // no .STEP line in the netlist
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      // The SUBTITLE data should include the values of the stepped parameters when .STEP is used.  
      // For the tecplot outputters, the step parameter information is not passed in as an argument.
      // However, the output manager has a copy, so get it from there.
      std::vector<Analysis::SweepParam> localStepSweepVector = outputManager_.getStepSweepVector();

      // Update elements of the local copy of stepSweepVector to have the correct 
      // values for currentStep_ and use those values for the SUBTITLE data.
      os << "SUBTITLE='Step param";
      for (std::vector<Analysis::SweepParam>::iterator it = localStepSweepVector.begin(); 
	   it != localStepSweepVector.end(); ++it)
      {
         os << " " << it->name << " = " << it->currentVal;
      }
    }

    os << "'" << std::endl;  // add trailing ' to the SUBTITLE line
    
    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
    os << "'" << std::endl;

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
      os << "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
    else
      os << "ANALYSIS='DC Sweep' " << "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='NO' " << 
       "NODES='" << printCount_ << "'" << std::endl;

    os << "SWEEPVAR='" << sweepVarName << "' SWEEPMODE='VAR_STEP'" << std::endl;
 
    // The next line assumes that we're doing a homotopy that goes from 0 to 1.
    // That assumption may be incorrect, and these two values should eventually
    // be filled in with the actual lower and upper limits on the first homotopy 
    // parameter, if those values are known, when this
    // header information is output.  This assumption does not appear to affect the
    // "viewability" of the resultant .csd file though.
    os << "XBEGIN='0.0'  XEND='1.0'"<<std::endl;

    os << "FORMAT='0 VOLTSorAMPS;EFLOAT : " << "NODEorBRANCH;NODE  '  " << std::endl;

    os << "DGTLDATA='NO'" << std::endl;

    // fill in the list of variable names in the #N section of the header block
    os << "#N" << std::endl;
    int maxLineLength = 128;
    int maxVarsPerLine = 12;
    int i = 0;  // used to count # of variables in the outString below
    std::string outString, currentVarString;

    for (std::vector<std::string>::const_iterator it = nameList.begin() ; it != nameList.end(); ++it)
    {
      // pad the currentVarString, so that V(1) is 'V(1)' followed by a space.
      // all subsequent operations then use that padded string.
      currentVarString = "'" + (*it) +"' ";

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
  ++numberOfSteps_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const Linear::Vector &                  solution_vector)
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
    homotopyHeader(comm,parameter_names, parameter_values, solution_vector);
  }
  else if (Parallel::rank(comm) == 0 && printStepHeader_)
  {
    // this prints out a header block at the start of steps 1,2,3, ...
    // this flag is set in doStartStep().  This approach did not
    // require changing the signature of doStartStep() for all of
    // the other outputters.
    homotopyHeader(comm, parameter_names, parameter_values, solution_vector);
    printStepHeader_ = false;
  }


  if (os_)
  {
    *os_ << "#C " << parameter_values[0] << " " << printCount_ << std::endl;
  }
    
  int colIdx = 0;  // used as a column counter local to this function 
  int maxColsPerRow = 4; // max number of columns output per row for each time value

  //-------------------------------------
  // HOMOTOPY PARAM VALUE OUTPUT GOES HERE.  Note: that the first homotopy
  // parameter is the SWEEPVAR.  So, it does not appear as a dependent variable
  // value
  //-------------------------------------
  for (int iparam=1;iparam < parameter_values.size(); ++iparam)
  {
    colIdx++;
    if (os_)
      *os_ << parameter_values[iparam] << ":" << colIdx << (colIdx%maxColsPerRow == 0 ? "\n" : "   ");
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(0, &solution_vector, 0, 0, 0, 0), result_list);
  
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
// Function      : HomotopyProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyProbe::doFinishOutput()
{
  if (os_)
  {
    (*os_) << "#;" << std::endl;

    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doStartStep(
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
// Function      : HomotopyProbe::doResetIndex
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doSteppingComplete()
{
  // close the homotopy file.
  //if (os_)
  //{
  //  (*os_) << "#;" << std::endl;
  //}
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

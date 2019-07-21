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

#include <N_IO_OutputterTimeRaw.h>

#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : TimeRaw
// Purpose       : Outputter class for transient output, rawfile output
//                 format
// Special Notes : Invoked by "FORMAT=raw" on .print line, not -r on command
//                 line.  -r is handled by the "OverrideRaw" classes.
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeRaw::TimeRaw
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeRaw::TimeRaw(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(""),
    numPoints_( 0),
    numPointsPos_( 0),
    os_( NULL),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::~TimeRaw
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeRaw::~TimeRaw()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeRaw::timeHeader()
{
  if (os_)
  {
    std::ostream &os = *os_;

    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_ = true;

      os << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      os << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepSweepVector().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      for (std::vector<Analysis::SweepParam>::const_iterator it = outputManager_.getStepSweepVector().begin();
          it != outputManager_.getStepSweepVector().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
      }
    }

    // while there is a doOutputHeaderAC(), AC can call this function
    // if it outputting a DC operating point.  Thus it's included
    // in this if statement.
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      plotName << "Transient Analysis";
    }
    else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC)
    {
      plotName << "DC operating point";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    os << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("real");
    os << "Flags: " << flags << std::endl;

    // prepare header for partial dump
    int numVars = opList_.size();

    os << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = os_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the
      // in the raw file header.  Optionally let one output the version
      // if whatever program is going to read the file expects it
      os << "Version: " << Util::Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    os << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;
    for (Util::Op::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      std::string varName = (*it)->getName();
      // Set the type.  Also change FREQ to FREQUENCY in the raw file header's
      // variable list.
      if (Util::hasExpressionTag((*it)->getName())) { tmpType = "expression"; }
      else if (varName == "INDEX")  { }
      else if (varName == "TIME")  { tmpType = "time"; }
      else if (varName == "FREQ")  { tmpType = "frequency"; varName = "FREQUENCY";}
      else if (varName == "sweep") { tmpType = "voltage"; }
      else if (varName[0] == 'I')  { tmpType = "current";    }
      else if (varName[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << varName
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeRaw::doOutputTime(
  Parallel::Machine     comm,
  const Linear::Vector &  solnVecPtr,
  const Linear::Vector &  stateVecPtr,
  const Linear::Vector &  storeVecPtr,
  const Linear::Vector &  lead_current_vector,
  const Linear::Vector &  junction_voltage_vector)
{
  if (Parallel::rank(comm) == 0 && os_ == 0)
  {
    outFilename_ = outputFilename(printParameters_.filename_, 
                                  printParameters_.defaultExtension_,
                                  printParameters_.suffix_+outputManager_.getFilenameSuffix(), 
                                  outputManager_.getNetlistFilename(),
                                  printParameters_.overrideRawFilename_,
                                  printParameters_.formatSupportsOverrideRaw_,
                                  printParameters_.dashoFilename_,
                                  printParameters_.fallback_);
    os_ = outputManager_.openBinaryFile(outFilename_);

    numPoints_ = 0;
  }

  if (numPoints_ == 0)
    timeHeader();

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(numPoints_, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0, &lead_current_vector, 0, &junction_voltage_vector), result_list);

  // select values to write from .PRINT line if FORMAT=RAW
  for (int i = 0; i < result_list.size(); ++i)
  {
    if (os_)
    {
      double x = result_list[i].real();
      os_->write((char *) &x, sizeof(x));
    }
  }

  // keep track of number of datapoints
  ++numPoints_;
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeRaw::doFinishOutput()
{
  if (os_)
  {
    if (numPoints_ != 0)
    {
      // need to move file pointer back to header and
      // write out the number of points.
      long currentFilePos = os_->tellp();

      // locate the position for number of points
      os_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*os_) << numPoints_;

      // move file pointer to the end again.
      os_->seekp( currentFilePos);
    }
  }

  // reset numPoints_ as it is used as a flag to print the header.
  numPoints_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

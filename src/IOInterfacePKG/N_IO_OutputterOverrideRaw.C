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

#include <N_IO_OutputterOverrideRaw.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : OverrideRaw
// Purpose       : Outputter class for rawfile output format
// Special Notes : Invoked by -r on command line.  FORMAT=RAW on print line
//                 is handled by the TimeRaw* and FrequencyRaw* classes
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : OverrideRaw::OverrideRaw
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OverrideRaw::OverrideRaw(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    numPoints_( 0),
    numPointsPos_( 0),
    os_( NULL),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::~OverrideRaw
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OverrideRaw::~OverrideRaw()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void OverrideRaw::timeHeader(Parallel::Machine comm)
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
      for (std::vector<Analysis::SweepParam>::const_iterator it = outputManager_.getStepSweepVector().begin(); it != outputManager_.getStepSweepVector().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
      }
    }

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      plotName << "Transient Analysis";
    }
    else if ( (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC) ||
              (printParameters_.analysisMode_ == Analysis::ANP_MODE_NOISE) )
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

  } // end proc0 check

  // prepare header for partial dump
  // get node types
  const std::vector< char > &typeRefs = outputManager_.getVarTypes();

  typedef std::vector<std::pair<std::string, char> > NodeNameTypeVector;
  NodeNameTypeVector local;

  // store local names and type
  if (solutionNodeNameMap_.empty())
    solutionNodeNameMap_.insert(outputManager_.getSolutionNodeMap().begin(), outputManager_.getSolutionNodeMap().end());
  for (OrderedNodeNameMap::const_iterator it = solutionNodeNameMap_.begin(), end = solutionNodeNameMap_.end(); it != end; ++it)
  {
    local.push_back(NodeNameTypeVector::value_type((*it).first, typeRefs[(*it).second]));
  }

  Util::Marshal mout;
  mout << local;

  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  if (os_)
  {
    std::ostream &os = *os_;

    NodeNameTypeVector global;
    for (int p = 0; p < Parallel::size(comm); ++p) {
      Util::Marshal min(dest[p]);

      std::vector<std::pair<std::string, char> > x;
      min >> x;
      global.insert(global.end(), x.begin(), x.end());
    }

    os << "No. Variables: " << global.size() + 1 << std::endl;

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

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpType;
    // write the variable information
    os << "Variables:" << std::endl;

    // add timestep header info(Spice3f5 style)
    os << "\t" << 0 << "\t";

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      os << "time\ttime\n";
    }
    else if ( (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC) ||
              (printParameters_.analysisMode_ == Analysis::ANP_MODE_NOISE) )
    {
      os << "frequency\tfrequency\n";
    }
    else
    {
      os << "sweep\tvoltage\n";
    }

    for (int i = 0; i < global.size(); ++i)
    {
      std::string tmpNodeName = global[i].first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" +(global[i]).first + ")";
      }

      std::string::size_type uPos = tmpNodeName.rfind( "_", tmpNodeName.size());
      if (uPos != std::string::npos)
      {
        tmpNodeName.replace( uPos, 1, "#");
      }

      os << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if (global[i].second == 'I' || global[i].second == 'i')
      {
        os << "current\n" ;
      }
      else
      {
        os << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  } // end proc0 check
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void OverrideRaw::frequencyHeader(Parallel::Machine comm)
{
  if (os_)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_=true;

      (*os_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*os_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepSweepVector().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      for (std::vector<Analysis::SweepParam>::const_iterator it = outputManager_.getStepSweepVector().begin(); it != outputManager_.getStepSweepVector().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
      }
    }
    if (!outputManager_.getDCSweepVector().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<Analysis::SweepParam>::const_iterator currItr = outputManager_.getDCSweepVector().begin();
      std::vector<Analysis::SweepParam>::const_iterator endItr = outputManager_.getDCSweepVector().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      plotName << "Transient Analysis";
    }
    else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC)
    {
      plotName << "AC Analysis";
    }
     else if (printParameters_.analysisMode_ == Analysis::ANP_MODE_NOISE)
    {
      plotName << "Noise Analysis";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*os_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    (*os_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  // get node types
  const std::vector< char > &typeRefs = outputManager_.getVarTypes();

  typedef std::vector<std::pair<std::string, char> > NodeNameTypeVector;
  NodeNameTypeVector local;

  // store local names and type
  if (solutionNodeNameMap_.empty())
    solutionNodeNameMap_.insert(outputManager_.getSolutionNodeMap().begin(), outputManager_.getSolutionNodeMap().end());
  for (OrderedNodeNameMap::const_iterator it = solutionNodeNameMap_.begin(), end = solutionNodeNameMap_.end(); it != end; ++it)
  {
    local.push_back(NodeNameTypeVector::value_type((*it).first, typeRefs[(*it).second]));
  }

  Util::Marshal mout;
  mout << local;

  std::vector<std::string> dest;
  Parallel::GatherV(comm, 0, mout.str(), dest);

  // format number of internal and external variables included here + time
  if (os_)
  {
    NodeNameTypeVector global;
    for (int p = 0; p < Parallel::size(comm); ++p) {
      Util::Marshal min(dest[p]);

      std::vector<std::pair<std::string, char> > x;
      min >> x;
      global.insert(global.end(), x.begin(), x.end());
    }

    (*os_) << "No. Variables: " << global.size() + 1 << std::endl;

    // format total number of data points & remember the streampos
    (*os_) << "No. Points: ";
    numPointsPos_ = os_->tellp();
    (*os_) << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the
      // in the raw file header.  Optionally let one output the version
      // if whatever program is going to read the file expects it
      (*os_) << "Version: " << Util::Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*os_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    // add timestep header info(Spice3f5 style)
    (*os_) << "\t" << 0 << "\t";

    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      (*os_) << "time\ttime\n";
    }
    else if ( (printParameters_.analysisMode_ == Analysis::ANP_MODE_AC) ||
              (printParameters_.analysisMode_ == Analysis::ANP_MODE_NOISE) )
    {
      (*os_) << "frequency\tfrequency\n";
    }
    else
    {
      (*os_) << "sweep\tvoltage\n";
    }

    std::string::size_type uPos;
    for ( int i = 0; i < global.size(); ++i)
    {
      std::string tmpNodeName = global[i].first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" + global[i].first + ")";
      }

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size());
      if (uPos != std::string::npos)
      {
        tmpNodeName.replace( uPos, 1, "#");
      }

      (*os_) << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if (global[i].second == 'I' || global[i].second == 'i')
      {
        (*os_) << "current\n" ;
      }
      else
      {
        (*os_) << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    (*os_) << "Binary:" << std::endl;
  } // end proc0 check
}


//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
OverrideRaw::doOutputTime(
  Parallel::Machine     comm,
  const Linear::Vector &  solnVecPtr,
  const Linear::Vector &  stateVecPtr,
  const Linear::Vector &  storeVecPtr,
  const Linear::Vector &  lead_current_vector,
  const Linear::Vector &  junction_voltage_vector)
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
    os_ = outputManager_.openBinaryFile(outFilename_);

    // set output value characteristics
    os_->setf( std::ios::scientific);
    os_->precision(8);
    os_->setf( std::ios::left, std::ios::adjustfield);

    numPoints_ = 0;
  }

  if (numPoints_ == 0)
    timeHeader(comm);

  // file IO on proc 0 only
  if (os_)
  {
    double tmp = 0.0;

// output the time/step values
    if (printParameters_.analysisMode_ == Analysis::ANP_MODE_TRANSIENT)
    {
      tmp = outputManager_.getCircuitTime();
      tmp *= printParameters_.outputTimeScaleFactor_;
    }
    else
    {
      tmp = outputManager_.getPRINTDCvalue();
    }
    os_->write((char *)&tmp , sizeof( double));
  }

  std::vector<double> local;
  for (OrderedNodeNameMap::const_iterator it = solutionNodeNameMap_.begin(), end = solutionNodeNameMap_.end(); it != end; ++it)
  {
    double result = solnVecPtr[(*it).second];
    result = filter(result, printParameters_.filter_);
    local.push_back(result);
  }

  std::vector<double> dest;
  Parallel::GatherV(comm, 0, local, dest);

  if (os_)
  {
    for (std::vector<double>::const_iterator it = dest.begin(), end = dest.end(); it != end; ++it)
    {
      double result = *it;

      // write binary data to rawfile
      os_->write((char *)&result , sizeof( double));
    }
  }

  // keep track of number of datapoints
  ++numPoints_;
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void OverrideRaw::doFinishOutput()
{
  if (os_)
  {
    // need to move file pointer back to header and
    // write out the number of points.
    long currentFelePost = os_->tellp();

    // locate the position for number of points
    os_->seekp( numPointsPos_);

    // overwrite blank space with value
    (*os_) << numPoints_;

    // move file pointer to the end again.
    os_->seekp( currentFelePost);
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
OverrideRaw::doOutputFrequency(
  Parallel::Machine     comm,
  double                frequency,
  double                fStart,
  double                fStop,
  const Linear::Vector &  real_solution_vector,
  const Linear::Vector &  imaginary_solution_vector,
  const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams)
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

    os_ = outputManager_.openBinaryFile(outFilename_);

    // set output value characteristics
    os_->setf( std::ios::scientific);
    os_->precision(8);
    os_->setf( std::ios::left, std::ios::adjustfield);

    numPoints_ = 0;
  }

  if (numPoints_ == 0)
    frequencyHeader(comm);

  double result = frequency;
  // file IO only on proc 0
  if (os_)
  {
    os_->write((char *) &result, sizeof(double));
    result = 0.0;
    os_->write((char *) &result, sizeof(double));
  }

  std::vector<complex> local;
  for (OrderedNodeNameMap::const_iterator it = solutionNodeNameMap_.begin(), end = solutionNodeNameMap_.end(); it != end; ++it)
  {
    complex result = complex(real_solution_vector[(*it).second], imaginary_solution_vector[(*it).second]);
    local.push_back(result);
  }

  std::vector<complex> dest;
  Parallel::GatherV(comm, 0, local, dest);

  if (os_)
  {
    for (std::vector<complex>::const_iterator it = dest.begin(), end = dest.end(); it != end; ++it)
    {
      complex result = *it;

      // write binary data to rawfile
      os_->write((char *)&result , sizeof(complex));
    }
  }

  // keep track of number of datapoints
  ++numPoints_;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

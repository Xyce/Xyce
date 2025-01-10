//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Outputter class for prn files for Embedded Sampling info.
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 7/26/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterEmbeddedSampling.h>
#include <N_IO_OutputterEmbeddedSamplingPrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::EmbeddedSamplingPrn
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
EmbeddedSamplingPrn::EmbeddedSamplingPrn(
    Parallel::Machine comm,
    OutputMgr &output_manager,
    const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
  {
    printParameters_.defaultExtension_ = ".ES.prn";
  }

  // Add columns for Index and TIME.  The columns for the PCE info will be added in
  // doOutputEmbeddedSampling().  Any variables on the .PRINT ES line will be ignored.
  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);

}

//-----------------------------------------------------------------------------
// Function      :EmbeddedSamplingPrn::~EmbeddedSamplingPrn
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
EmbeddedSamplingPrn::~EmbeddedSamplingPrn()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::EmbeddedSamplingHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingPrn::EmbeddedSamplingHeader()
{
  int column_index = 0;

  Table::ColumnList::const_iterator it = printParameters_.table_.columnList_.begin();
  Table::ColumnList::const_iterator end = printParameters_.table_.columnList_.end();
  for ( ; it != end; ++it, ++column_index)
  {
    if (it != printParameters_.table_.columnList_.begin())
    {
      *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
    }
    printHeader(*os_, (*it));
  }

  *os_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::doOutputEmbeddedSampling
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingPrn::doOutputEmbeddedSampling(
  Parallel::Machine comm,
  bool              regressionPCEenable,
  bool              projectionPCEenable,
  int               numSamples,
  const std::vector<std::string> & regressionPCEcoeffs,
  const std::vector<std::string> & projectionPCEcoeffs,
  const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec)
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

    // Generate names for the  additional header columns needed for PCE output
    std::vector<std::string> colNames;
    makeEmbeddedSamplingColumnNames(printParameters_, colNames, regressionPCEenable, projectionPCEenable,
                                    numSamples, regressionPCEcoeffs, projectionPCEcoeffs, outFuncDataVec);

    // add those additional columns to the printParameters_.table_
    fixupColumnsFromStrVec(comm, printParameters_, colNames);

    if (outputManager_.getPrintHeader())
    {
      // output the column names to the output file.
      EmbeddedSamplingHeader();
    }
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, 0, 0, 0, 0, 0), result_list);

  if (os_)
  {
    // Output the STEPNUM, Index and TIME values.
    for (int i = 0; i < result_list.size(); ++i)
    {
      if (os_)
        printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
    }

    // Now output the PCE values directly from the elements of the outFuncDataVec.
    // The next output column starts after the columns generated from the opList_
    outputEmbeddedSamplingData(printParameters_, os_, result_list,
			       regressionPCEenable, projectionPCEenable,
                               numSamples, outFuncDataVec);

    // send end-of-line character
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::doFinishOutput
// Purpose       : Output the footer, and close the stream if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingPrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        (*os_) << "End of Xyce(TM) Embedded Sampling Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingPrn::doStartStep( int current_step, int number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;

  // If using Format::GNUPLOT then add two blank lines before the output for
  // steps 1, 2, ... if there is a .STEP loop.  (Note: currentStep_ goes from
  // 0 to numberOfSteps_-1 if there is a .STEP loop.)
  //
  // If using Format::SPLOT then add a single blank line before the output for
  // steps 1, 2, ... if there is a .STEP loop.
  if (os_)
  {
    if ( (printParameters_.addGnuplotSpacing_) && (currentStep_ > 0 ) )
    {
      *os_ << std::endl << std::endl;
    }
    else if ( (printParameters_.addSplotSpacing_) && (currentStep_ > 0 ) )
    {
      *os_ << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingPrn::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingPrn::doSteppingComplete
// Purpose       : Output footer and close the stream  when a .STEP loop
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingPrn::doSteppingComplete()
{
  // close the embedded sampling file.
  if (os_)
  {
    // this end-of-simulation footer is used if there is a .STEP loop
    if ( outputManager_.getPrintFooter ())
    {
      (*os_) << "End of Xyce(TM) Embedded Sampling Simulation" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

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
// Purpose        : Outputter class for csv files for PCE info.
//
// Special Notes  :
//
// Creator        : Pete Sholander
//
// Creation Date  : 9/3/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterPCE.h>
#include <N_IO_OutputterPCECSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : PCECSV::PCECSV
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
PCECSV::PCECSV(
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
    printParameters_.defaultExtension_ = ".PCE.csv";
  }

  // Add a column for TIME.  The columns for the PCE info will be added in doOutputPCE().
  // Any variables on the .PRINT PCE line will be ignored.
  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : PCECSV::~PCECSV
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
PCECSV::~PCECSV()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : PCECSV::PCEHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCECSV::PCEHeader()
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
// Function      : PCECSV::doOutputPCE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCECSV::doOutputPCE(
  Parallel::Machine comm,
  int               numQuadPoints,
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
    makePCEColumnNames(printParameters_, colNames, numQuadPoints, outFuncDataVec);

    // add those additional columns to the printParameters_.table_
    fixupColumnsFromStrVec(comm, printParameters_, colNames);

    // output the column names to the output file.
    PCEHeader();
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, 0, 0, 0, 0, 0), result_list);

  if (os_)
  {
    // Output the TIME values.
    for (int i = 0; i < result_list.size(); ++i)
    {
      if (os_)
        printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
    }

    // Now output the PCE values directly from the the elements of the outFuncDataVec.
    // The next output column starts after the columns generated from the opList_
    outputPCEData(printParameters_, os_, result_list, numQuadPoints, outFuncDataVec);

    // send end-of-line character
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : PCECSV::doFinishOutput
// Purpose       : Close the stream if there is no .STEP loop.  This function
//                 is also called after each step, if there is a .STEP loop,
//                 but currently does nothing in that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCECSV::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PCECSV::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCECSV::doStartStep( int current_step, int number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : PCECSV::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCECSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : PCECSV::doSteppingComplete
// Purpose       : Close the stream when a .STEP loop is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCECSV::doSteppingComplete()
{
  // close the sensitivity file.
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

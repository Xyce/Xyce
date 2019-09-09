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
// Purpose        : Outputter class for Tecplot files for PCE info.
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 9/3/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterPCE.h>
#include <N_IO_OutputterPCETecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_IO_Tecplot.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : PCETecplot::PCETecplot
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
PCETecplot::PCETecplot(
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
    printParameters_.defaultExtension_ = ".PCE.dat";
  }

  // Add column for TIME variable.  The columns for the PCE info will be added in
  // doOutputPCE().  Any variables on the .PRINT PCE line will be ignored.
  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::~PCETecplot
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
PCETecplot::~PCETecplot()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::PCEHeader
// Purpose       :
// Special Notes : Adapted from tecplotFreqHeader function in N_IO_Tecplot.C
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCETecplot::PCEHeader()
{
  std::ostream &os = *os_;
  bool print_title =  (currentStep_ == 0);
  const std::string title = outputManager_.getNetlistFilename() + " - " + outputManager_.getTitle();
  const Table::ColumnList & columnList_ = printParameters_.table_.columnList_;

  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);
  os.precision(2);

  if (print_title)
  {
    std::string localTitle(title);
    std::size_t loc=0;
    while( (loc = localTitle.find_first_of( '"',loc ) ) != std::string::npos )
    {
      localTitle.insert(loc, 1, '\\');
      loc+=2;
    }
    os << " TITLE = \"" << localTitle << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    Table::ColumnList::const_iterator it2 = columnList_.begin();
    Table::ColumnList::const_iterator end2 = columnList_.end();
    for ( ; it2 != end2; ++it2)
    {
      os << "\" ";
      if (it2 != columnList_.begin())
      {
        *os_ << printParameters_.delimiter_;
      }
      printHeader(*os_, (*it2));
      os << "\" " << std::endl;
    }

    os << "DATASETAUXDATA ";
    os << getTecplotTimeDateStamp();
    os << std::endl;

    if (!outputManager_.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << outputManager_.getCircuitTemp() << " \"" << std::endl;
    }

    // this outputs the AuxData for Step 0
    outputAuxData(os);
  }
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::outputAuxData
// Purpose       : output some AUXDATA
// Special Notes : Adapted from tecplotFreqHeader function in N_IO_Tecplot.C
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCETecplot::outputAuxData(std::ostream &os)
{
  os << "ZONE F=POINT  ";

  if (outputManager_.getStepSweepVector().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<Analysis::SweepParam>::const_iterator it = outputManager_.getStepSweepVector().begin(); it != outputManager_.getStepSweepVector().end(); ++it)
    {
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  // put in the various sweep parameters as auxdata:
  if (!outputManager_.getStepSweepVector().empty())
  {
    for (std::vector<Analysis::SweepParam>::const_iterator iterParam = outputManager_.getStepSweepVector().begin();
    iterParam != outputManager_.getStepSweepVector().end();
    ++iterParam)
    {
      // convert any ":", "%" or "!" in the name to a "_", so as not to confuse tecplot.
      std::string tmpName(iterParam->name);
      replace(tmpName.begin(), tmpName.end(), '%', '_');
      replace(tmpName.begin(), tmpName.end(), ':', '_');
      replace(tmpName.begin(), tmpName.end(), '!', '_');
      os << "AUXDATA " << tmpName << " = " << "\" " << iterParam->currentVal << "\" ";
    }
    os << std::endl;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::doOutputPCE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCETecplot::doOutputPCE(
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

    if (outputManager_.getPrintHeader())
    {
      // Generate names for the  additional header columns needed for PCE output
      std::vector<std::string> colNames;
      makePCEColumnNames(printParameters_, colNames, numQuadPoints, outFuncDataVec);

      // add those additional columns to the printParameters_.table_
      fixupColumnsFromStrVec(comm, printParameters_, colNames);

      // output the column names to the output file.
      PCEHeader();
    }
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

    // Now output the PCE values directly from the elements of the outFuncDataVec.
    // The next output column starts after the columns generated from the opList_
    outputPCEData(printParameters_, os_, result_list, numQuadPoints, outFuncDataVec);

    // send end-of-line character
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::doFinishOutput
// Purpose       : Output the footer, and close the stream if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCETecplot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        (*os_) << "End of Xyce(TM) PCE Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes : This output the AUXDATA for steps 1,2,...)
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void PCETecplot::doStartStep( int current_step, int number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;

  // Output the AUXDATA if this is not Step 0.  The AUXDATA for Step 0
  // was output in the function PCEHeader()
  if (os_ && current_step != 0)
  {
    std::ostream &os = *os_;
    outputAuxData(os);
  }
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void PCETecplot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : PCETecplot::doSteppingComplete
// Purpose       : Output footer and close the stream  when a .STEP loop
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void PCETecplot::doSteppingComplete()
{
  // close the PCE file.
  if (os_)
  {
    // this end-of-simulation footer is used if there is a .STEP loop
    if ( outputManager_.getPrintFooter ())
    {
      (*os_) << "End of Xyce(TM) PCE Simulation" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

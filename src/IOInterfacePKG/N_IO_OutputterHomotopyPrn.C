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

#include <N_IO_OutputterHomotopyPrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {



//-----------------------------------------------------------------------------
// Class         : HomotopyPrn
// Purpose       : Outputter class for homotopy output, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::HomotopyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyPrn::HomotopyPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0),
    homotopyParamStartIndex_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".HOMOTOPY.prn";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);

  // adjust where the homotopy params start, in the output column list, based
  // on whether the STEPNUM and Index columns are output
  if (printParameters_.printStepNumColumn_)
    ++homotopyParamStartIndex_;
  if (printParameters_.printIndexColumn_)
    ++homotopyParamStartIndex_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::~HomotopyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyPrn::~HomotopyPrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::homotopyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::homotopyHeader(
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const Linear::Vector &                  solution_vector)
{
  if (columnList_.empty())
  {
    Table::Justification justification = printParameters_.delimiter_.empty() ?
      Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator it = parameter_names.begin();
        it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific,
            printParameters_.streamWidth_, printParameters_.streamPrecision_, justification));
    }
  }

  index_ = 0;

  if (currentStep_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator it = printParameters_.table_.columnList_.begin();
        it != printParameters_.table_.columnList_.end(); ++it, ++column_index)
    {
      if (it != printParameters_.table_.columnList_.begin())
      {
        *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }

      if (column_index == homotopyParamStartIndex_)
      {
        for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2)
        {
          if (it2 != columnList_.begin())
          {
            *os_ << printParameters_.delimiter_;
          }
          printHeader(*os_, (*it2));
        }
      }

      printHeader(*os_, (*it));
    }

    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyPrn::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const Linear::Vector &                  solution_vector)
{
  double tmpTime = outputManager_.getCircuitTime(); // outputManager_.getAnaIntPtr()->getTime();

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
        homotopyHeader(parameter_names, parameter_values, solution_vector);
      }
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, &solution_vector, 0, 0, 0, 0), result_list);

  if (Parallel::rank(comm) == 0)
  {
    for (int i = 0; i < result_list.size(); ++i) {
      if (i == homotopyParamStartIndex_)
        for (int j = 0; j < parameter_values.size(); ++j)
          printValue(*os_, columnList_[j], printParameters_.delimiter_, 1, parameter_values[j]);
       
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
    }
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doFinishOutput
// Purpose       : Output the footer, and close the stream if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in 
//               : that case.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyPrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
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
// Function      : HomotopyPrn::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyPrn::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doSteppingComplete
// Purpose       : Output footer and close the stream when a .STEP loop 
//               : is used.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doSteppingComplete()
{
  // close the homotopy file.
  if (os_)
  {
    // this end-of-simulation footer is used if there is a .STEP loop
    if ( outputManager_.getPrintFooter() )
    {
      (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}


} // namespace Outputter
} // namespace IO
} // namespace Xyce

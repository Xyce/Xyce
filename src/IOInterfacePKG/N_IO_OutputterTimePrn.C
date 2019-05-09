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

#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_PDS_MPI.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : TimePrn::TimePrn
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 06/07/2013
//-----------------------------------------------------------------------------
TimePrn::TimePrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".prn";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::~TimePrn
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 06/07/2013
//-----------------------------------------------------------------------------
TimePrn::~TimePrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doOutputTime
// Purpose       : Output the current data at a time point
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void
TimePrn::doOutputTime(
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
    os_ = outputManager_.openFile(outFilename_);

    if (outputManager_.getPrintHeader())
    {
      printHeader(*os_, printParameters_);
    }
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0, &lead_current_vector, 0, &junction_voltage_vector), result_list);

  for (int i = 0; i < result_list.size(); ++i) {
    result_list[i] = complex(filter(result_list[i].real(), printParameters_.filter_), 0.0);

    if (os_)
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
  }

  if (os_)
    *os_ << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doFinishOutput
// Purpose       : Output the footer, and close the stream if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in 
//               : that case.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        *os_ << "End of Xyce(TM) Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimePrn::doStartStep(
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
// Function      : TimePrn::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimePrn::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doSteppingComplete
// Purpose       : Output footer and close the stream  when a .STEP loop 
//               : is used.
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doSteppingComplete()
{
  if (os_)
  {
    // this end-of-simulation footer is used if there is a .STEP loop
    if (outputManager_.getPrintFooter())
    {
      *os_ << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : TimeCSV
// Purpose       : Outputter class for transient runs and "CSV"
//                 (comma separarated variable) output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeCSV::TimeCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeCSV::TimeCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csv";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::~TimeCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeCSV::~TimeCSV()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : TimeCSV::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::timeHeader(Parallel::Machine comm)
{
  if (Parallel::rank(comm) == 0)
  {
    printHeader(*os_, printParameters_);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputPRINT_
// Purpose       : .PRINT output
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void TimeCSV::doOutputTime(
  Parallel::Machine     comm,
  const Linear::Vector &  solnVecPtr,
  const Linear::Vector &  stateVecPtr,
  const Linear::Vector &  storeVecPtr,
  const Linear::Vector &  lead_current_vector,
  const Linear::Vector &  junction_voltage_vector)
{
  double time = outputManager_.getCircuitTime();

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
    timeHeader(comm);
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0, &lead_current_vector, 0, &junction_voltage_vector), result_list);

  for (int i = 0; i < result_list.size(); ++i)
  {
    result_list[i] = filter(result_list[i].real(), printParameters_.filter_);

    if (os_)
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
  }

  if (os_)
    *os_ << std::endl;

  ++index_;
}

void TimeCSV::doFinishOutput()
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
// Function      : TimeCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeCSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeCSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::doSteppingComplete()
{
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

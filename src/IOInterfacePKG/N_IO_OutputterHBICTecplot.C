//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Outputter class for tecplot files for HB_IC info.
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 11/19/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <fstream>

#include <N_IO_OutputterHBICTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : HBICTecplot::HBICTecplot
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
HBICTecplot::HBICTecplot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    tmpOutStream_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".hb_ic.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot::~HBICTecplot
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
HBICTecplot::~HBICTecplot()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
void HBICTecplot::doOutputTime(
  Parallel::Machine           comm,
  const Linear::Vector &        solnVecPtr,
  const Linear::Vector &        stateVecPtr,
  const Linear::Vector &        storeVecPtr,
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
    os_->setf(std::ios::scientific);
    os_->precision(printParameters_.streamPrecision_);
    os_->setf(std::ios::left, std::ios::adjustfield);
  }

  // If transient-assisted HB is used with .STEP then only the last (accepted) transient run
  // should be written to the output file.  The data points from the intermediate transient
  // runs are written to a tmp file first.  See SON Bug 928 for more details.
  if (Parallel::rank(comm) == 0 && !tmpOutStream_)
  {
    if ((outputManager_.getTaHBSpecified()) && (outputManager_.getStepSweepVector().size() > 0))
    {
      tmpOutStream_ = outputManager_.openFile(outFilename_ + ".tmp");
    }
  }

  if (tmpOutStream_ && index_ == 0)
    tecplotTimeHeader(*tmpOutStream_, currentStep_ == 0, outputManager_.getNetlistFilename() + " - " + outputManager_.getTitle(), opList_, outputManager_);
  else if (os_ && index_ == 0)
    tecplotTimeHeader(*os_, currentStep_ == 0, outputManager_.getNetlistFilename() + " - " + outputManager_.getTitle(), opList_, outputManager_);

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(0, &solnVecPtr, 0, &stateVecPtr, &storeVecPtr, 0, &lead_current_vector, 0, &junction_voltage_vector), result_list);

  for (int i = 0; i < result_list.size(); ++i)
  {
    result_list[i] = filter(result_list[i].real(), printParameters_.filter_);

    if (tmpOutStream_)
      (*tmpOutStream_) << std::setw(printParameters_.streamWidth_) << result_list[i].real() << " ";
    else if (os_)
      (*os_) << std::setw(printParameters_.streamWidth_) << result_list[i].real() << " ";
  }

  if (tmpOutStream_)
    (*tmpOutStream_) << std::endl;
  else if (os_)
    (*os_) << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
void HBICTecplot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter())
      {
        (*os_) << "End of Xyce(TM) Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot:::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
void HBICTecplot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
void HBICTecplot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot:: doSteppingComplete()
// Purpose       : Output footer and close the stream  when a .STEP loop
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
void HBICTecplot::doSteppingComplete()
{
  if (os_)
  {
    if (outputManager_.getPrintFooter())
    {
      (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot::reopenTmpFile
// Purpose       : close and re-open tmp file used when transient-assisted HB
//               : is used with .STEP
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/18/2019
//-----------------------------------------------------------------------------
void HBICTecplot::reopenTmpFile()
{
  if ((outputManager_.getTaHBSpecified()) && (outputManager_.getStepSweepVector().size() > 0)
      && tmpOutStream_ )
  {
    outputManager_.closeFile(tmpOutStream_);
    tmpOutStream_ = outputManager_.openFile(outFilename_ + ".tmp");
  }
}

//-----------------------------------------------------------------------------
// Function      : HBICTecplot:copyTmpFileToOutputFile
// Purpose       : copy tmp file to actual output file, for the case when
//               : transient-assisted HB is used with .STEP
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/18/2019
//-----------------------------------------------------------------------------
void HBICTecplot::copyTmpFileToOutputFile()
{
  if ( tmpOutStream_)
  {
    outputManager_.closeFile(tmpOutStream_);

    if ((outputManager_.getTaHBSpecified()) && (outputManager_.getStepSweepVector().size() > 0)
       && os_)
    {
      // re-open tmp file for input
      std::string tmpFilename = outFilename_ + ".tmp";
      std::ifstream infile;
      infile.open(tmpFilename.c_str(),std::ifstream::in);

      if (infile.good())
        *os_ << infile.rdbuf();

      // close and delete tmp file
      infile.close();
      std::remove(tmpFilename.c_str());
    }

    tmpOutStream_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

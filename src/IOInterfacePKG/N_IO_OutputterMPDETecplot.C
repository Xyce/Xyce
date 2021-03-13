//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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

#include <N_IO_OutputterMPDETecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_LAS_BlockVector.h>
#include <N_MPDE_Manager.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : MPDETecplot
// Purpose       : Outputter class for MPDE runs, Tecplot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDETecplot::MPDETecplot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
MPDETecplot::MPDETecplot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    n1_(0),
    n2_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".MPDE.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);

}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::~MPDETecplot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDETecplot::~MPDETecplot()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doOutputHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::mpdeHeader()
{
  (*os_) << " TITLE = \" Xyce MPDE data, " << outputManager_.getNetlistFilename() << "\", " << std::endl
	           << "\tVARIABLES = \"TIME1 \", \"TIME2\", " << std::endl;

  // output the user-specified solution vars:
  for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
  {
    (*os_) << "\" "<< (*it)->getName() << "\" " << std::endl;
  }

  // output some AUXDATA
  (*os_) << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl
                   << "ZONE I=" << n2_ + 1 << ", " << " J=" << n1_ << ", " << " F=POINT\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doOutputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points, 
  const Linear::BlockVector &     solution_block_vector)
{
  int blockCount = solution_block_vector.blockCount();
  n2_ = blockCount; // fast time points.
  ++n1_;            // slow time points.  increments by one each call.

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

  if (os_ && index_ == 0)
  {
    mpdeHeader();
  }

  // Loop over the fast time points of the Linear::BlockVecor:
  for (int iblock=0;iblock<n2_+1;++iblock)
  {
    const Linear::Vector &solnVecPtr = (iblock == n2_ ? solution_block_vector.block(0) : solution_block_vector.block(iblock));

    if (Parallel::rank(comm) == 0)
    {
      //-------------------------------------
      // Get the 2 time values first.
      //-------------------------------------
      double first  = 0.0;
      double second = 0.0;

      second = fast_time_points[iblock];
      first  = time;

      // time 1:
      (*os_) << std::setw(printParameters_.streamWidth_) << first << " " << std::setw(printParameters_.streamWidth_) << second;
    }

    std::vector<complex> result_list;
    getValues(comm, opList_, Util::Op::OpData(0, &solnVecPtr, 0, 0, 0, 0), result_list);

    for (int i = 0; i < result_list.size(); ++i)
    {
      if (os_)
      {
        (*os_) << result_list[i].real();
      }
    }

    if (os_)
    {
      (*os_) << "\n";
    }
  }

  if (os_)
  {
    (*os_) << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
      {
        (*os_) << "End of Xyce(TM) Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecplot::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecplot::doSteppingComplete()
{
  if (os_)
  {
    if (outputManager_.getPrintFooter ())
    {
      (*os_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

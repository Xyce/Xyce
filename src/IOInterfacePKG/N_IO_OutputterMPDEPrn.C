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

#include <N_IO_OutputterMPDEPrn.h>
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
// Class         : MPDEPrn
// Purpose       : Outputter class for MPDE runs, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDEPrn::MPDEPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDEPrn::MPDEPrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
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
    printParameters_.defaultExtension_ = ".MPDE.prn";

  printParameters_.table_.addColumn("TIME1", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_CENTER);
  printParameters_.table_.addColumn("TIME2", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_CENTER);

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::~MPDEPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDEPrn::~MPDEPrn()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::mpdeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::mpdeHeader()
{}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doOutputMPDE(
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

    if (outputManager_.getPrintHeader())
    {
      printHeader(*os_, printParameters_);
    }
  }

  // Loop over the fast time points of the Linear::BlockVecor:
  for (int iblock=0;iblock<n2_+1;++iblock)
  {
    const Linear::Vector &solution_vector = (iblock == n2_ ? solution_block_vector.block(0) : solution_block_vector.block(iblock));

    if (os_)
    {
      //-------------------------------------
      // Get the 2 time values first.
      //-------------------------------------
      double first  = time;
      double second = fast_time_points[iblock]; //genie 121913. This should be a bug since fast_time_points has size = n2_ so its index is at most n2_-1

      // time 1:
      printValue(*os_, printParameters_.table_.columnList_[0], printParameters_.delimiter_, 0, first);

      // time 2:
      printValue(*os_, printParameters_.table_.columnList_[1], printParameters_.delimiter_, 1, second);
    }

    std::vector<complex> result_list;
    getValues(comm, opList_, Util::Op::OpData(0, &solution_vector, 0, 0, 0, 0), result_list);

    for (int i = 0; i < result_list.size(); ++i)
    {
      if (os_)
      {
        printValue(*os_, printParameters_.table_.columnList_[i + 2], printParameters_.delimiter_, i + 2, result_list[i].real());
      }
    }

    if (os_)
    {
      (*os_) << std::endl;
    }
  } // fast time scale loop.

  if (os_)
  {
    if ( (printParameters_.addGnuplotSpacing_) || (printParameters_.addSplotSpacing_) )
    {
      *os_ << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        (*os_) << "End of Xyce(TM) Simulation" << std::endl;
      }
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
MPDEPrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
MPDEPrn::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doSteppingComplete()
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

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
// Creator        : Pete Sholander, Electrical Models and Simulation
//
// Creation Date  : 11/29/2018
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterMPDECSV.h>
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
// Class         : MPDECSV
// Purpose       : Outputter class for MPDE runs, CSV output
//                 format
// Special Notes :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDECSV::MPDECSV
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
MPDECSV::MPDECSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
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
    printParameters_.defaultExtension_ = ".MPDE.csv";

  printParameters_.table_.addColumn("TIME1", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_NONE);
  printParameters_.table_.addColumn("TIME2", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_NONE);

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::~MPDECSV
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
MPDECSV::~MPDECSV()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::mpdeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void MPDECSV::mpdeHeader()
{}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void MPDECSV::doOutputMPDE(
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
    printHeader(*os_, printParameters_);
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
}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void MPDECSV::doFinishOutput()
{
  outputManager_.closeFile(os_);
  os_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void
MPDECSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void
MPDECSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDECSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 10/29/2018
//-----------------------------------------------------------------------------
void MPDECSV::doSteppingComplete()
{
  // close the file.
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

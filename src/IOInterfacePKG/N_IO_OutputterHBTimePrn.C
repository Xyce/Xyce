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
// Purpose        : The outputter class for .PRINT HB_TD for FORMAT=STD
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL, Electrical Models and Simulation
//
// Creation Date  : 6/25/2018
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHBTimePrn.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::HBTimePrn
// Purpose       : Constructor for outputter class for .PRINT HB_TD output for
//                 STD output format
// Special Notes :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
HBTimePrn::HBTimePrn(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &time_print_parameters)
  : outputManager_(output_manager),
    timePrintParameters_(time_print_parameters),
    timeFilename_(),
    tos_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (timePrintParameters_.defaultExtension_.empty())
    timePrintParameters_.defaultExtension_ = ".HB.TD.prn";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), timePrintParameters_, timeOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::~HBTimePrn
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
HBTimePrn::~HBTimePrn()
{
  outputManager_.closeFile(tos_);

  deleteList(timeOpList_.begin(), timeOpList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::doOutputHB_TD
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBTimePrn::doOutputHB_TD(
  Parallel::Machine             comm,
  const std::vector<double> &   timePoints,
  const Linear::BlockVector &     timeDomainSolutionVec,
  const Linear::BlockVector &     timeDomainLeadCurrentVec,
  const Linear::BlockVector &     timeDomainJunctionVoltageVec)
{
  int blockCount = timeDomainSolutionVec.blockCount();

  if (Parallel::rank(comm) == 0 && !tos_)
  {
    timeFilename_ = outputFilename(timePrintParameters_.filename_, 
                                   timePrintParameters_.defaultExtension_,
                                   timePrintParameters_.suffix_+outputManager_.getFilenameSuffix(), 
                                   outputManager_.getNetlistFilename(),
                                   timePrintParameters_.overrideRawFilename_,
                                   timePrintParameters_.formatSupportsOverrideRaw_,
                                   timePrintParameters_.dashoFilename_,
                                   timePrintParameters_.fallback_); 

    tos_ = outputManager_.openFile(timeFilename_);

    if (outputManager_.getPrintHeader())
    {
      printHeader(*tos_, timePrintParameters_);
    }
  }

  // Loop over the time points of the Linear::BlockVecor:
  for (int iblock = 0; iblock < blockCount; ++iblock)
  {
    outputManager_.setCircuitTime(timePoints[iblock]);
    Linear::Vector * time_solution_vector = &timeDomainSolutionVec.block(iblock);
    Linear::Vector * time_lead_current_vector = &timeDomainLeadCurrentVec.block(iblock);

    Linear::Vector * time_junction_voltage_vector = &timeDomainJunctionVoltageVec.block(iblock);

    { // periodic time-domain steady-state output
      std::vector<complex> result_list;
      getValues(comm, timeOpList_, Util::Op::OpData(index_, time_solution_vector, 0, 0, 0, 0, time_lead_current_vector, 0, time_junction_voltage_vector), result_list);
      for (int i = 0; i < result_list.size(); ++i)
        if (tos_)
          printValue(*tos_, timePrintParameters_.table_.columnList_[i], timePrintParameters_.delimiter_, i, result_list[i].real());
    }

    if (tos_)
      *tos_ << std::endl;

    ++index_;
  }

}

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::doFinishOutput
// Purpose       : Output the footers, and close the streams if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in 
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBTimePrn::doFinishOutput()
{
  if (tos_)
  {
    if (numberOfSteps_ == 0)
    {
      // this end-of-simulation footer is used if there is no .STEP loop
      if (outputManager_.getPrintFooter())
      {
        *tos_ << "End of Xyce(TM) Simulation" << std::endl;
      }

      outputManager_.closeFile(tos_);
      tos_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void
HBTimePrn::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;

  // If using Format::GNUPLOT then add two blank lines before the output for
  // steps 1, 2, ... if there is a .STEP loop.  (Note: currentStep_ goes from 
  // 0 to numberOfSteps_-1 if there is a .STEP loop.)  Do this for both the
  // .TD.prn and .FD.prn files.
  //
  // If using Format::SPLOT then add a single blank line before the output for
  // steps 1, 2, ... if there is a .STEP loop.  
  if (tos_)
  {
    if ( (timePrintParameters_.addGnuplotSpacing_) && (currentStep_ > 0 ) )
    {
      *tos_ << std::endl << std::endl;
    }
    else if ( (timePrintParameters_.addSplotSpacing_) && (currentStep_ > 0 ) )
    {
      *tos_ << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void
HBTimePrn::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBTimePrn::doSteppingComplete
// Purpose       : Output footers and close the streams when a .STEP loop 
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBTimePrn::doSteppingComplete()
{ 
 if (outputManager_.getPrintFooter()) 
 {
    // this end-of-simulation footer is used if there is a .STEP loop
    if (tos_)
      (*tos_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(tos_);
  tos_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

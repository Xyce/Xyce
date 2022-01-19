//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : The outputter class for .PRINT HB_TD for FORMAT=TECPLOT
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

#include <N_IO_OutputterHBTimeTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::HBTimeTecPlot
// Purpose       : Constructor for outputter class for .PRINT HB_TD output for
//                 TECPLOT output
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
HBTimeTecPlot::HBTimeTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &time_print_parameters)
  : outputManager_(output_manager),
    timePrintParameters_(time_print_parameters),
    timeFilename_(),
    tos_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (timePrintParameters_.defaultExtension_.empty())
    timePrintParameters_.defaultExtension_ = ".HB.TD.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), timePrintParameters_, timeOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::~HBTimeTecPlot
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
HBTimeTecPlot::~HBTimeTecPlot()
{
  outputManager_.closeFile(tos_);

  deleteList(timeOpList_.begin(), timeOpList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::doOutputHB_TD
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBTimeTecPlot::doOutputHB_TD(
  Parallel::Machine             comm,
  const std::vector<double> &   timePoints,
  const Linear::BlockVector &   timeDomainSolutionVec,
  const Linear::BlockVector &   timeDomainLeadCurrentVec,
  const Linear::BlockVector &   timeDomainJunctionVoltageVec)
{
  int index = 0;

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

  }

  if (tos_ && index_ == 0)
    tecplotTimeHeader(*tos_, currentStep_ == 0, outputManager_.getNetlistFilename(), timeOpList_, outputManager_);

  // Loop over the time points of the Linear::BlockVecor:
  for (int iblock = 0; iblock < blockCount; ++iblock)
  {
    outputManager_.setCircuitTime(timePoints[iblock]);

    Linear::Vector *solution_vector = &(timeDomainSolutionVec.block(iblock));
    Linear::Vector * time_lead_current_vector = &timeDomainLeadCurrentVec.block(iblock);
    Linear::Vector * time_junction_voltage_vector = &timeDomainJunctionVoltageVec.block(iblock);

    { // periodic time-domain steady-state output
      std::vector<complex> result_list;
      getValues(comm, timeOpList_, Util::Op::OpData(index_, solution_vector, 0, 0, 0, 0, time_lead_current_vector, 0, time_junction_voltage_vector), result_list);

      for (int i = 0; i < result_list.size(); ++i)
        if (tos_)
          printValue(*tos_, timePrintParameters_.table_.columnList_[i], timePrintParameters_.delimiter_, i, result_list[i].real());
    }
    if (tos_)
      *tos_ << std::endl;

    ++index;
  }
  
}

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::doFinishOutput
// Purpose       : Output the footers, and close the streams if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in 
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBTimeTecPlot::doFinishOutput()
{
  if (tos_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
        *tos_ << "End of Xyce(TM) Simulation" << std::endl;

      outputManager_.closeFile(tos_);
      tos_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void
HBTimeTecPlot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : This function is executed at the start of each step.
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void
HBTimeTecPlot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBTimeTecPlot::doSteppingComplete
// Purpose       : Output footers and close the streams when a .STEP loop 
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBTimeTecPlot::doSteppingComplete()
{
  if (outputManager_.getPrintFooter ()) 
  {
    if (tos_)
      (*tos_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(tos_);
  tos_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

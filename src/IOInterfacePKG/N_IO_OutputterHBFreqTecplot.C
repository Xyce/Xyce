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
// Purpose        : The outputter class for .PRINT HB_FD for FORMAT=TECPLOT
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

#include <N_IO_OutputterHBFreqTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::HBFreqTecPlot
// Purpose       : Constructor for outputter class for .PRINT HB_FD output for
//                 TECPLOT output
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
HBFreqTecPlot::HBFreqTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &freq_print_parameters)
  : outputManager_(output_manager),
    freqPrintParameters_(freq_print_parameters),
    freqFilename_(),
    fos_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (freqPrintParameters_.defaultExtension_.empty())
    freqPrintParameters_.defaultExtension_ = ".HB.FD.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), freqPrintParameters_, freqOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::~HBFreqTecPlot
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
HBFreqTecPlot::~HBFreqTecPlot()
{
  outputManager_.closeFile(fos_);

  deleteList(freqOpList_.begin(), freqOpList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::doOutputHB_FD
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBFreqTecPlot::doOutputHB_FD(
  Parallel::Machine             comm,
  const std::vector<double> &   freqPoints,
  const Linear::BlockVector &   freqDomainSolutionVecReal,
  const Linear::BlockVector &   freqDomainSolutionVecImaginary,
  const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
  const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary)
{
  int index = 0;

   if (Parallel::rank(comm) == 0 && !fos_)
  {
    freqFilename_ = outputFilename(freqPrintParameters_.filename_,
                                   freqPrintParameters_.defaultExtension_,
                                   freqPrintParameters_.suffix_+outputManager_.getFilenameSuffix(), 
                                   outputManager_.getNetlistFilename(),
                                   freqPrintParameters_.overrideRawFilename_,
                                   freqPrintParameters_.formatSupportsOverrideRaw_,
                                   freqPrintParameters_.dashoFilename_,
                                   freqPrintParameters_.fallback_); 
    
    fos_ = outputManager_.openFile(freqFilename_);

  }

  if (fos_ && index_ == 0)
    tecplotFreqHeader(*fos_, currentStep_ == 0, outputManager_.getNetlistFilename(), freqOpList_, outputManager_);

  int fdblockCount = freqDomainSolutionVecReal.blockCount();
  index_ = 0;
  for (int iblock = 0; iblock < fdblockCount; ++iblock)
  {
    outputManager_.setCircuitFrequency(freqPoints[iblock]);

    Linear::Vector *real_solution_vector = &(freqDomainSolutionVecReal.block(iblock));
    Linear::Vector *imaginary_solution_vector = &(freqDomainSolutionVecImaginary.block(iblock));
    Linear::Vector * real_lead_current_vector = &freqDomainLeadCurrentVecReal.block(iblock);
    Linear::Vector * imaginary_lead_current_vector = &freqDomainLeadCurrentVecImaginary.block(iblock);
    Linear::Vector * real_junction_voltage_vector = &freqDomainJunctionVoltageVecReal.block(iblock);
    Linear::Vector * imaginary_junction_voltage_vector = &freqDomainJunctionVoltageVecImaginary.block(iblock);
    { // Fourier coefficient output
      std::vector<complex> result_list;

      getValues(comm, freqOpList_, Util::Op::OpData(index_, real_solution_vector, imaginary_solution_vector, 0, 0, 0, real_lead_current_vector, imaginary_lead_current_vector, real_junction_voltage_vector, imaginary_junction_voltage_vector), result_list);

      for (int i = 0; i < result_list.size(); ++i)
        if (fos_)
          printValue(*fos_, freqPrintParameters_.table_.columnList_[i], freqPrintParameters_.delimiter_, i, result_list[i].real());
    }

    if (fos_)
      *fos_ << std::endl;

    ++index_;
  }
}

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::doFinishOutput
// Purpose       : Output the footers, and close the streams if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in 
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBFreqTecPlot::doFinishOutput()
{
  if (fos_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
        *fos_ << "End of Xyce(TM) Simulation" << std::endl;

      outputManager_.closeFile(fos_);
      fos_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void
HBFreqTecPlot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void
HBFreqTecPlot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBFreqTecPlot::doSteppingComplete
// Purpose       : Output footers and close the streams when a .STEP loop 
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 6/25/2018
//-----------------------------------------------------------------------------
void HBFreqTecPlot::doSteppingComplete()
{
  if (outputManager_.getPrintFooter ()) 
  {
    if (fos_)
      (*fos_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(fos_);
  fos_ = 0;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

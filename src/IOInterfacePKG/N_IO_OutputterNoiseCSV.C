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
// Purpose        : Generate noise output in CSV format
//
// Special Notes  :
//
// Creator        : Pete Sholander, Electrical Models and Simulation
//
// Creation Date  : 8/24/2017
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterNoiseCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : NoiseCSV
// Purpose       : Outputter class for noise output, (CSV) output
//                 format
// Special Notes :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 8/24/2017
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : NoiseCSV::NoiseCSV
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 8/24/2017
//-----------------------------------------------------------------------------
NoiseCSV::NoiseCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".NOISE.csv";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : NoiseCSV::~NoiseCSV
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 8/24/2017
//-----------------------------------------------------------------------------
NoiseCSV::~NoiseCSV()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : NoiseCSV::noiseHeader
// Purpose       : Print out the header line in the .csv file
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 8/24/2017
//-----------------------------------------------------------------------------
void NoiseCSV::noiseHeader()
{
  if (os_ && currentStep_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator
        it = printParameters_.table_.columnList_.begin();
        it != printParameters_.table_.columnList_.end();
        ++it, ++column_index)
    {
      if (it != printParameters_.table_.columnList_.begin())
      {
        *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }
      printHeader(*os_, (*it));
    }

    for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2)
    {
      if (it2 != columnList_.begin())
      {
        *os_ << printParameters_.delimiter_;
      }
      printHeader(*os_, (*it2));
    }
    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : NoiseCSV::doOutputNoise
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date :
//-----------------------------------------------------------------------------
void NoiseCSV::doOutputNoise(
  Parallel::Machine   comm,
  double              frequency,
  const Linear::Vector &real_solution_vector, 
  const Linear::Vector &imaginary_solution_vector,
  double              totalOutputNoiseDens_, 
  double              totalInputNoiseDens_, 
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_)
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

    printHeader(*os_, printParameters_);
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, &real_solution_vector, &imaginary_solution_vector, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, totalOutputNoiseDens_, totalInputNoiseDens_, &noiseDataVec_  ), result_list);

  for (int i = 0; i < result_list.size(); ++i)
  {
    result_list[i] = complex(filter(result_list[i].real(), printParameters_.filter_), 0.0);

    if (os_)
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
  }

  if (os_)
    *os_ << std::endl;

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : NoiseCSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date :
//-----------------------------------------------------------------------------
void NoiseCSV::doFinishOutput()
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
// Function      : NoiseCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void NoiseCSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : NoiseCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void NoiseCSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : NoiseCSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void NoiseCSV::doSteppingComplete()
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

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
// Purpose        : Generate homotopy output in CSV format
//
// Special Notes  :
//
// Creator        : Pete Sholander, Electrical Models and Simulation
//
// Creation Date  : 7/24/2017
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHomotopyCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Probe.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {
namespace Outputter {



//-----------------------------------------------------------------------------
// Class         : HomotopyCSV
// Purpose       : Outputter class for homotopy output, (CSV) output
//                 format
// Special Notes :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::HomotopyCSV
// Purpose       : xonstructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
HomotopyCSV::HomotopyCSV(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".HOMOTOPY.csv";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::~HomotopyCSV
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
HomotopyCSV::~HomotopyCSV()
{
  outputManager_.closeFile(os_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::homotopyHeader
// Purpose       : Print out the header line in the .csv file
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
void HomotopyCSV::homotopyHeader(
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           param_values,
  const Linear::Vector &                  solution_vector)
{
  if (columnList_.empty())
  {
    Table::Justification justification = printParameters_.delimiter_.empty() ?
      Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator it = parameter_names.begin();
        it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific,
            printParameters_.streamWidth_, printParameters_.streamPrecision_, justification));
    }
  }

  index_ = 0;

  if (currentStep_ == 0)
  {
    int column_index = 0;
    // columns for the homotopy parameters
    for (Table::ColumnList::const_iterator it = columnList_.begin(); it != columnList_.end(); ++it)
    {
      printHeader(*os_, (*it));
      *os_ << printParameters_.delimiter_;
    }

    // columns for the variables on the .PRINT line, plus variables like TIME
    for (Table::ColumnList::const_iterator it2 = printParameters_.table_.columnList_.begin();
        it2 != printParameters_.table_.columnList_.end(); ++it2)
    {
      if (it2 != printParameters_.table_.columnList_.begin())
      {
        *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }
      printHeader(*os_, (*it2));
    }

    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::doOutputHomotopy
// Purpose       : Print out data rows in the .csv file
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
void
HomotopyCSV::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const Linear::Vector &                  solution_vector)
{
  double tmpTime = outputManager_.getCircuitTime(); // outputManager_.getAnaIntPtr()->getTime();

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

      homotopyHeader(parameter_names, parameter_values, solution_vector);
  }

  // Since CSV format does not have the Index column, use 0 for start of homotopy params
  int homotopyParamStartIndex=0;

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, &solution_vector, 0, 0, 0, 0), result_list);

  if (Parallel::rank(comm) == 0)
  {
    for (int i = 0; i < result_list.size(); ++i) {
      if (i == homotopyParamStartIndex)
        for (int j = 0; j < parameter_values.size(); ++j)
          printValue(*os_, columnList_[j], printParameters_.delimiter_, j, parameter_values[j]);
       
      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i+parameter_values.size(), result_list[i].real());
    }
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
void HomotopyCSV::doFinishOutput()
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
// Function      : HomotopyCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
void
HomotopyCSV::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
void
HomotopyCSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyCSV::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL, Electrical Models and Simulation
// Creation Date : 7/24/2017
//-----------------------------------------------------------------------------
void HomotopyCSV::doSteppingComplete()
{
  // close the homotopy file.
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

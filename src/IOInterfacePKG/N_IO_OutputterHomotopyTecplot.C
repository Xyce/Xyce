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
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterHomotopyTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Tecplot.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : HomotopyTecPlot
// Purpose       : Outputter class for homotopy output, TecPlot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::HomotopyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyTecPlot::HomotopyTecPlot(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".HOMOTOPY.dat";

  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::~HomotopyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyTecPlot::~HomotopyTecPlot()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::homotopyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::homotopyHeader(
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

  std::ostream &os = *os_;

  index_ = 0;
  if (currentStep_ == 0)
  {
    os << " TITLE = \" Xyce homotopy data, " << outputManager_.getNetlistFilename() << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the continuation parameters:
    std::vector<std::string>::const_iterator iter_name;
    for (iter_name = parameter_names.begin(); iter_name!= parameter_names.end(); ++iter_name)
    {
      os << "\" ";
      os << *iter_name;
      os << "\" " << std::endl;
    }

    // output the user-specified solution vars:
    for (Util::Op::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
    {
      os << "\" " << (*it)->getName() << "\" " << std::endl;
    }
  }

  // output some AUXDATA
  os << "DATASETAUXDATA ";
  os << getTecplotTimeDateStamp();
  os << std::endl;

  os << "ZONE F=POINT";

  if (outputManager_.getStepSweepVector().empty())
  {
    // no .STEP line in the netlist
    os << " T=\"Xyce data\" ";
  }
  else
  {
    // The ZONE data should include the values of the stepped parameters when .STEP is used.  
    // For the tecplot outputters, the step parameter information is not passed in as an argument.
    // However, the output manager has a copy, so get it from there.
    std::vector<Analysis::SweepParam> localStepSweepVector = outputManager_.getStepSweepVector();

    os << " T= \" "; 
    
    // Update elements of the local copy of stepSweepVector to have the correct 
    // values for currentStep_ and use those values for the SUBTITLE data.
    for (std::vector<Analysis::SweepParam>::iterator it = localStepSweepVector.begin(); 
	 it != localStepSweepVector.end(); ++it)
    {
      static const int tecplotHeaderPrecision = 2;
      os.setf(std::ios::scientific);
      os.precision(tecplotHeaderPrecision);
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doOutputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const Linear::Vector &                  solution_vector)
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

  if (os_ && index_ == 0)
    homotopyHeader(parameter_names, parameter_values, solution_vector);

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(0, &solution_vector, 0, 0, 0, 0), result_list);

  if (os_) {
    for (int i = 0; i < result_list.size(); ++i) {
      if (i == 0)
      {
        for (int i = 0; i < parameter_values.size(); ++i)
        {
          printValue(*os_, columnList_[i], printParameters_.delimiter_, 1, parameter_values[i]);
        }
      }

      printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
    }
  }

  if (os_)
  {
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter ())
      {
        (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doStartStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : void HomotopyTecPlot::doResetIndex
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doSteppingComplete
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doSteppingComplete()
{
  if (os_)
  {
    if (outputManager_.getPrintFooter () )
    {
      (*os_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
    }
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

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

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputterSensitivity.h>
#include <N_IO_OutputterSensitivityACPrn.h>
#include <N_IO_OutputterSensitivityACCSV.h>
#include <N_IO_OutputterSensitivityACTecplot.h>
#include <N_IO_OutputterSensitivityPrn.h>
#include <N_IO_OutputterSensitivityCSV.h>
#include <N_IO_OutputterSensitivityTecplot.h>
#include <N_IO_OutputterSensitivityDakota.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableSensitivityOutput
// Purpose       : Enable sensitivity output for analysis modes other than AC
// Special Notes :
// Scope         :
// Creator       : Dave Baur
// Creation Date : June 23, 2014
//-----------------------------------------------------------------------------
void enableSensitivityOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::SENS);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) 
    {
      PrintParameters sensitivity_print_parameters = (*it);

      if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      if (sensitivity_print_parameters.printIndexColumn_)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (sensitivity_print_parameters.printStepNumColumn_)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      output_manager.fixupPrintParameters(comm, sensitivity_print_parameters);

      Outputter::Interface *outputter;
      if (sensitivity_print_parameters.format_ == Format::STD) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
        outputter = new Outputter::SensitivityPrn(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::CSV) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.csv";
        outputter = new Outputter::SensitivityCSV(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::TECPLOT) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.dat";
        outputter = new Outputter::SensitivityTecPlot(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::DAKOTA) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".SENS.txt";
        outputter = new Outputter::SensitivityDakota(comm, output_manager, sensitivity_print_parameters);
      }
      else if ( (sensitivity_print_parameters.format_ == Format::RAW) ||
                (sensitivity_print_parameters.format_ == Format::RAW_ASCII) ||
                (sensitivity_print_parameters.format_ == Format::PROBE) ||
                (sensitivity_print_parameters.format_ == Format::TS1) ||
                (sensitivity_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "Sensitivity output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        sensitivity_print_parameters.format_ = Format::STD;
        outputter = new Outputter::SensitivityPrn(comm, output_manager, sensitivity_print_parameters); 
      }
      else
      {
        Report::UserWarning0()
          << "Sensitivity output cannot be written in requested format, using standard format";
        sensitivity_print_parameters.format_ = Format::STD;
        sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
        outputter = new Outputter::SensitivityPrn(comm, output_manager, sensitivity_print_parameters);
      }

      output_manager.addOutputter(PrintType::SENS, outputter);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : enableSensitivityACOutput
// Purpose       : Enable sensitivity output for AC analysis mode.  This will
//                 include both the direct and adjoint information.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 4/16/2019
//-----------------------------------------------------------------------------
void enableSensitivityACOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::SENS);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it)
    {
      PrintParameters sensitivity_print_parameters = (*it);
      sensitivity_print_parameters.expandComplexTypes_ = true;

      sensitivity_print_parameters.variableList_.push_front(Util::Param("FREQ", 0.0));
      if (sensitivity_print_parameters.printIndexColumn_)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (sensitivity_print_parameters.printStepNumColumn_)
        sensitivity_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      output_manager.fixupPrintParameters(comm, sensitivity_print_parameters);

      Outputter::Interface *outputter;
      if (sensitivity_print_parameters.format_ == Format::STD)
      {
        sensitivity_print_parameters.defaultExtension_ = ".FD.SENS.prn";
        outputter = new Outputter::SensitivityACPrn(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::CSV)
      {
        sensitivity_print_parameters.defaultExtension_ = ".FD.SENS.csv";
        outputter = new Outputter::SensitivityACCSV(comm, output_manager, sensitivity_print_parameters);
      }
      else if (sensitivity_print_parameters.format_ == Format::TECPLOT) 
      {
        sensitivity_print_parameters.defaultExtension_ = ".FD.SENS.dat";
        outputter = new Outputter::SensitivityACTecplot(comm, output_manager, sensitivity_print_parameters);
      }
      else if ( (sensitivity_print_parameters.format_ == Format::RAW) ||
                (sensitivity_print_parameters.format_ == Format::RAW_ASCII) ||
                (sensitivity_print_parameters.format_ == Format::PROBE) ||
                (sensitivity_print_parameters.format_ == Format::DAKOTA) ||
                (sensitivity_print_parameters.format_ == Format::TS1) ||
                (sensitivity_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "Sensitivity output cannot be written in PROBE, RAW, Dakota or Touchstone format, using standard format instead";
        sensitivity_print_parameters.format_ = Format::STD;
        outputter = new Outputter::SensitivityACPrn(comm, output_manager, sensitivity_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "AC Sensitivity output can only be written in standard format";
        sensitivity_print_parameters.format_ = Format::STD;
        sensitivity_print_parameters.defaultExtension_ = ".FD.SENS.prn";
        outputter = new Outputter::SensitivityACPrn(comm, output_manager, sensitivity_print_parameters);
      }

      output_manager.addOutputter(PrintType::SENS, outputter);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : enableAdjointSensitivityOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/12/2016
//-----------------------------------------------------------------------------
void enableAdjointSensitivityOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::TRANADJOINT);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) 
    {
      PrintParameters transientAdjoint_print_parameters = (*it);

      if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
      {
        transientAdjoint_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      }

      if (transientAdjoint_print_parameters.printIndexColumn_)
      {
        transientAdjoint_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      }

      output_manager.fixupPrintParameters(comm, transientAdjoint_print_parameters);

      Outputter::Interface *outputter;
      if (transientAdjoint_print_parameters.format_ == Format::STD) 
      {
        transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.prn";
        outputter = new Outputter::SensitivityPrn(comm, output_manager, transientAdjoint_print_parameters);
      }
      else if (transientAdjoint_print_parameters.format_ == Format::CSV) 
      {
        transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.csv";
        outputter = new Outputter::SensitivityCSV(comm, output_manager, transientAdjoint_print_parameters);
      }
      else if (transientAdjoint_print_parameters.format_ == Format::TECPLOT) 
      {
        transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.dat";
        outputter = new Outputter::SensitivityTecPlot(comm, output_manager, transientAdjoint_print_parameters);
      }
      else if (transientAdjoint_print_parameters.format_ == Format::DAKOTA) 
      {
        transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.txt";
        outputter = new Outputter::SensitivityDakota(comm, output_manager, transientAdjoint_print_parameters);
      }
      else if ( (transientAdjoint_print_parameters.format_ == Format::RAW) ||
                (transientAdjoint_print_parameters.format_ == Format::RAW_ASCII) ||
                (transientAdjoint_print_parameters.format_ == Format::PROBE) )
      {
        Report::UserWarning0()
          << "Transient adjoint output cannot be written in PROBE or RAW format, using standard format instead";
        transientAdjoint_print_parameters.format_ = Format::STD;
        outputter = new Outputter::SensitivityPrn(comm, output_manager, transientAdjoint_print_parameters); 
      }

      else
      {
        Report::UserWarning0()
          << "Sensitivity output cannot be written in requested format, using standard format";
        transientAdjoint_print_parameters.format_ = Format::STD;
        transientAdjoint_print_parameters.defaultExtension_ = ".TRADJ.prn";
        outputter = new Outputter::SensitivityPrn(comm, output_manager, transientAdjoint_print_parameters);
      }

      output_manager.addOutputter(PrintType::TRANADJOINT, outputter);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

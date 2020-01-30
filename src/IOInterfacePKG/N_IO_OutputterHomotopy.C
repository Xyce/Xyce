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

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputterHomotopy.h>
#include <N_IO_OutputterHomotopyCSV.h>
#include <N_IO_OutputterHomotopyPrn.h>
#include <N_IO_OutputterHomotopyTecplot.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableHomotopyOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Baur
// Creation Date : June 23, 2014
//-----------------------------------------------------------------------------
void enableHomotopyOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::HOMOTOPY);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it) {
      PrintParameters homotopy_print_parameters = (*it);

      if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
        homotopy_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      if (homotopy_print_parameters.printIndexColumn_)
        homotopy_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (homotopy_print_parameters.printStepNumColumn_)
        homotopy_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      output_manager.fixupPrintParameters(comm, homotopy_print_parameters);

      Outputter::Interface *outputter;
      if (homotopy_print_parameters.format_ == Format::STD) {
        outputter = new Outputter::HomotopyPrn(comm, output_manager, homotopy_print_parameters);
      }
      else if (homotopy_print_parameters.format_ == Format::TECPLOT) {
        outputter = new Outputter::HomotopyTecPlot(comm, output_manager, homotopy_print_parameters);
      }
      else if (homotopy_print_parameters.format_ == Format::CSV) {
        outputter = new Outputter::HomotopyCSV(comm, output_manager, homotopy_print_parameters);
      }
      else if ( (homotopy_print_parameters.format_ == Format::RAW) ||
              (homotopy_print_parameters.format_ == Format::RAW_ASCII) ||
              (homotopy_print_parameters.format_ == Format::PROBE) ||
              (homotopy_print_parameters.format_ == Format::TS1) ||
              (homotopy_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "Homotopy output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        homotopy_print_parameters.format_ = Format::STD;
        outputter = new Outputter::HomotopyPrn(comm, output_manager, homotopy_print_parameters);  
      }
      else
      {
        Report::UserWarning0()
          << "Homotopy output cannot be written in requested format, using standard format instead";
        homotopy_print_parameters.format_ = Format::STD;
        outputter = new Outputter::HomotopyPrn(comm, output_manager, homotopy_print_parameters);
      }

      output_manager.addOutputter(PrintType::HOMOTOPY, outputter);
    }
  }
}


} // namespace Outputter
} // namespace IO
} // namespace Xyce

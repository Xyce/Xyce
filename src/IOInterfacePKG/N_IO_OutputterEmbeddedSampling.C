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
// Purpose        : Outputter class for Embedded Sampling
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 7/26/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputterEmbeddedSampling.h>
#include <N_IO_OutputterEmbeddedSamplingPrn.h>
#include <N_IO_OutputterEmbeddedSamplingCSV.h>
#include <N_IO_OutputterEmbeddedSamplingTecplot.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableEmbeddedSamplingOutput
// Purpose       : Enable EmbeddedSampling output
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : July 26, 2019
//-----------------------------------------------------------------------------
void enableEmbeddedSamplingOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::ES);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) 
    {
      PrintParameters es_print_parameters = (*it);

      if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
        es_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      if (es_print_parameters.printIndexColumn_)
        es_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (es_print_parameters.printStepNumColumn_)
        es_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      output_manager.fixupPrintParameters(comm, es_print_parameters);

      Outputter::Interface *outputter;
      if (es_print_parameters.format_ == Format::STD)
      {
        es_print_parameters.defaultExtension_ = ".ES.prn";
        outputter = new Outputter::EmbeddedSamplingPrn(comm, output_manager, es_print_parameters);
      }
      else if (es_print_parameters.format_ == Format::CSV)
      {
        es_print_parameters.defaultExtension_ = ".ES.csv";
        outputter = new Outputter::EmbeddedSamplingCSV(comm, output_manager, es_print_parameters);
      }
      else if (es_print_parameters.format_ == Format::TECPLOT)
      {
        es_print_parameters.defaultExtension_ = ".ES.dat";
        outputter = new Outputter::EmbeddedSamplingTecplot(comm, output_manager, es_print_parameters);
      }
      else if ( (es_print_parameters.format_ == Format::RAW) ||
                (es_print_parameters.format_ == Format::RAW_ASCII) ||
                (es_print_parameters.format_ == Format::PROBE) ||
                (es_print_parameters.format_ == Format::TS1) ||
                (es_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "Embedded sampling output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        es_print_parameters.format_ = Format::STD;
        outputter = new Outputter::EmbeddedSamplingPrn(comm, output_manager, es_print_parameters);
      }
      else
      {
        Report::UserWarning0()
          << "Embedded Sampling output cannot be written in requested format, using standard format";
        es_print_parameters.format_ = Format::STD;
        es_print_parameters.defaultExtension_ = ".ES.prn";
        outputter = new Outputter::EmbeddedSamplingPrn(comm, output_manager, es_print_parameters);
      }

      output_manager.addOutputter(PrintType::ES, outputter);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

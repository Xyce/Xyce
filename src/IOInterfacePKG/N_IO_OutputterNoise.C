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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  :
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterNoise.h>
#include <N_IO_OutputterNoisePrn.h>
#include <N_IO_OutputterNoiseCSV.h>
#include <N_IO_OutputterNoiseTecplot.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableNoiseOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 01/26/2015
//-----------------------------------------------------------------------------
void enableNoiseOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::NOISE);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it) 
    {
      PrintParameters noise_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, noise_print_parameters);

      Outputter::Interface *outputter;
      if (noise_print_parameters.format_ == Format::STD) 
      {
        outputter = new Outputter::NoisePrn(comm, output_manager, noise_print_parameters);
      }
      else if (noise_print_parameters.format_ == Format::CSV) 
      {
        outputter = new Outputter::NoiseCSV(comm, output_manager, noise_print_parameters);
      }
      else if (noise_print_parameters.format_ == Format::TECPLOT) 
      {
        outputter = new Outputter::NoiseTecPlot(comm, output_manager, noise_print_parameters);
      }
      else if ( (noise_print_parameters.format_ == Format::RAW) ||
                (noise_print_parameters.format_ == Format::RAW_ASCII) ||
                (noise_print_parameters.format_ == Format::PROBE) ||
                (noise_print_parameters.format_ == Format::TS1) ||
                (noise_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "Noise output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        noise_print_parameters.format_ = Format::STD;
        outputter = new Outputter::NoisePrn(comm, output_manager, noise_print_parameters); 
      }
      else
      {
        Report::UserWarning0() << "Noise output cannot be written in " << noise_print_parameters.format_ << " format, using standard format";
        noise_print_parameters.format_ = Format::STD;
        outputter = new Outputter::NoisePrn(comm, output_manager, noise_print_parameters);
      }

      output_manager.addOutputter(PrintType::NOISE, outputter);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

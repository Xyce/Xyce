//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Creator        : Pete Sholander
//
// Creation Date  : 3/25/19
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterSParam.h>
#include <N_IO_OutputterSParamTS1.h>
#include <N_IO_OutputterSParamTS2.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableSParamOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2019
//-----------------------------------------------------------------------------
void enableSParamOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::SPARAM);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it) 
    {
      PrintParameters sparam_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, sparam_print_parameters);

      Outputter::Interface *outputter;
      if (sparam_print_parameters.format_ == Format::TS1) 
      {
        outputter = new Outputter::SParamTS1(comm, output_manager, sparam_print_parameters);
      }
      else if (sparam_print_parameters.format_ == Format::TS2) 
      {
        outputter = new Outputter::SParamTS2(comm, output_manager, sparam_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "SParam output can only be written in Touchstone format, using Touchstone2 format";
        sparam_print_parameters.format_ = Format::TS2;
        outputter = new Outputter::SParamTS2(comm, output_manager, sparam_print_parameters);
      }

      output_manager.addOutputter(PrintType::SPARAM, outputter);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

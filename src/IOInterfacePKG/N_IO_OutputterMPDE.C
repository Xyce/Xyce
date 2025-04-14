//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Outputter for MPDE
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
#include <N_IO_OutputterTransient.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterFrequencyPrn.h>
#include <N_IO_OutputterFrequencyCSV.h>
#include <N_IO_OutputterMPDEPrn.h>
#include <N_IO_OutputterMPDECSV.h>
#include <N_IO_OutputterMPDETecplot.h>
#include <N_IO_OutputterTimeTecplot.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableMPDEOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Baur
// Creation Date : June 23, 2014
//-----------------------------------------------------------------------------
void enableMPDEOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  // set up outputter for MPDE
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::MPDE);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters mpde_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, mpde_print_parameters);

      Outputter::Interface *outputter_mpde;
      if (mpde_print_parameters.format_ == Format::STD) 
      {
        outputter_mpde = new Outputter::MPDEPrn(comm, output_manager, mpde_print_parameters);
      }
      else if (mpde_print_parameters.format_ == Format::CSV)
      {
        outputter_mpde = new Outputter::MPDECSV(comm, output_manager, mpde_print_parameters);
      }
      else if (mpde_print_parameters.format_ == Format::TECPLOT) 
      {
        outputter_mpde = new Outputter::MPDETecplot(comm, output_manager, mpde_print_parameters);
      }
      else if ( (mpde_print_parameters.format_ == Format::RAW) ||
                (mpde_print_parameters.format_ == Format::RAW_ASCII) ||
                (mpde_print_parameters.format_ == Format::PROBE) ||
                (mpde_print_parameters.format_ == Format::TS1) ||
                (mpde_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "MPDE output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        mpde_print_parameters.format_ = Format::STD;
        outputter_mpde = new Outputter::MPDEPrn(comm, output_manager, mpde_print_parameters); 
      }
      else
      {
        Report::UserWarning0() << "MPDE output cannot be written in " << mpde_print_parameters.format_ << " format, using standard format";
        mpde_print_parameters.format_ = Format::STD;
        outputter_mpde = new Outputter::MPDEPrn(comm, output_manager, mpde_print_parameters);
      }

      output_manager.addOutputter(PrintType::MPDE, outputter_mpde);
    }
  }
  
  // setup outputter for MPDE_IC
  result = output_manager.findOutputParameter(OutputType::MPDE_IC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters mpde_ic_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, mpde_ic_print_parameters);

      Outputter::Interface *outputter_mpde_ic;
      if (mpde_ic_print_parameters.format_ == Format::STD) 
      {
        outputter_mpde_ic = new Outputter::TimePrn(comm, output_manager, mpde_ic_print_parameters);
      }
      else if (mpde_ic_print_parameters.format_ == Format::CSV) 
      {
        outputter_mpde_ic = new Outputter::TimeCSV(comm, output_manager, mpde_ic_print_parameters);
      }
      else if (mpde_ic_print_parameters.format_ == Format::TECPLOT) 
      {
        outputter_mpde_ic = new Outputter::TimeTecplot(comm, output_manager, mpde_ic_print_parameters);
      }
      else if ( (mpde_ic_print_parameters.format_ == Format::RAW) ||
                (mpde_ic_print_parameters.format_ == Format::RAW_ASCII) ||
                (mpde_ic_print_parameters.format_ == Format::PROBE) ||
                (mpde_ic_print_parameters.format_ == Format::TS1) ||
                (mpde_ic_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "MPDE_IC output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        mpde_ic_print_parameters.format_ = Format::STD;
        outputter_mpde_ic = new Outputter::TimePrn(comm, output_manager, mpde_ic_print_parameters); 
      }
      else
      {
        Report::UserWarning0() << "MPDE_IC output cannot be written in " << mpde_ic_print_parameters.format_ << " format, using standard format";
        mpde_ic_print_parameters.format_ = Format::STD;
        outputter_mpde_ic = new Outputter::TimePrn(comm, output_manager, mpde_ic_print_parameters);
      }

      output_manager.addOutputter(PrintType::MPDE_IC, outputter_mpde_ic);
    }
  }

  // set up outputter for MPDE_STARTUP
  result = output_manager.findOutputParameter(OutputType::MPDE_STARTUP);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters mpde_startup_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, mpde_startup_print_parameters);

      Outputter::Interface *outputter_startup;
      if (mpde_startup_print_parameters.format_ == Format::STD) 
      {
        outputter_startup = new Outputter::TimePrn(comm, output_manager, mpde_startup_print_parameters);
      }
      else if (mpde_startup_print_parameters.format_ == Format::CSV) 
      {
        outputter_startup = new Outputter::TimeCSV(comm, output_manager, mpde_startup_print_parameters);
      }
      else if (mpde_startup_print_parameters.format_ == Format::TECPLOT) 
      {
        outputter_startup = new Outputter::TimeTecplot(comm, output_manager, mpde_startup_print_parameters);
      }
      else if ( (mpde_startup_print_parameters.format_ == Format::RAW) ||
                (mpde_startup_print_parameters.format_ == Format::RAW_ASCII) ||
                (mpde_startup_print_parameters.format_ == Format::PROBE) ||
                (mpde_startup_print_parameters.format_ == Format::TS1) ||
                (mpde_startup_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "MPDE_STARTUP output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        mpde_startup_print_parameters.format_ = Format::STD;
        outputter_startup = new Outputter::TimePrn(comm, output_manager, mpde_startup_print_parameters); 
      }
      else
      {
        Report::UserWarning0() << "MPDE_STARTUP output cannot be written in " << mpde_startup_print_parameters.format_ << " format, using standard format";
        mpde_startup_print_parameters.format_ = Format::STD;
        outputter_startup = new Outputter::TimePrn(comm, output_manager, mpde_startup_print_parameters);
      }

      output_manager.addOutputter(PrintType::MPDE_STARTUP, outputter_startup);
    }
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

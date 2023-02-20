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
#include <N_IO_OutputterHBTimePrn.h>
#include <N_IO_OutputterHBFreqPrn.h>
#include <N_IO_OutputterHBICPrn.h>
#include <N_IO_OutputterHBTimeCSV.h>
#include <N_IO_OutputterHBFreqCSV.h>
#include <N_IO_OutputterHBICCSV.h>
#include <N_IO_OutputterHBTimeTecplot.h>
#include <N_IO_OutputterHBFreqTecplot.h>
#include <N_IO_OutputterHBICTecplot.h>
#include <N_IO_OutputterTimeTecplot.h>
#include <N_IO_OutputterExternal.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableHBOutput
// Purpose       : This function sets up both the frequency-domain and
//                 time-domain output for a .HB analysis.  It is used by
//                 .PRINT HB_FD, .PRINT HB_TD, .PRINT HB_IC and .PRINT 
//                 HB_STARTUP lines.
// Special Notes :
// Scope         :
// Creator       : Dave Baur
// Creation Date : June 23, 2014
//-----------------------------------------------------------------------------
void enableHBOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  // generate HB_FD outputters
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::HB_FD);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end();
         it != end; ++it) {
      PrintParameters freq_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, freq_print_parameters);

      Outputter::Interface *outputter_hb;
      if (freq_print_parameters.format_ == Format::STD) {
        outputter_hb = new Outputter::HBFreqPrn(comm, output_manager, freq_print_parameters);
      }
      else if (freq_print_parameters.format_ == Format::CSV) {
        outputter_hb = new Outputter::HBFreqCSV(comm, output_manager, freq_print_parameters);
      }
      else if (freq_print_parameters.format_ == Format::TECPLOT) {
        outputter_hb = new Outputter::HBFreqTecPlot(comm, output_manager, freq_print_parameters);
      }
      else if ( (freq_print_parameters.format_ == Format::RAW) ||
             (freq_print_parameters.format_ == Format::RAW_ASCII) ||
             (freq_print_parameters.format_ == Format::PROBE) ||
             (freq_print_parameters.format_ == Format::TS1) ||
             (freq_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0() << "HB_FD output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        freq_print_parameters.format_ = Format::STD;
        outputter_hb = new Outputter::HBFreqPrn(comm, output_manager, freq_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB_FD output cannot be written in requested format, using standard format instead";
        freq_print_parameters.format_ = Format::STD;
        outputter_hb = new Outputter::HBFreqPrn(comm, output_manager, freq_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB_FD, outputter_hb);
    }
  }

  // generate HB_TD outputters
  result = output_manager.findOutputParameter(OutputType::HB_TD);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end();
         it != end; ++it) {
      PrintParameters time_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, time_print_parameters);

      Outputter::Interface *outputter_hb;
      if (time_print_parameters.format_ == Format::STD) {
        outputter_hb = new Outputter::HBTimePrn(comm, output_manager, time_print_parameters);
      }
      else if (time_print_parameters.format_ == Format::CSV) {
        outputter_hb = new Outputter::HBTimeCSV(comm, output_manager, time_print_parameters);
      }
      else if (time_print_parameters.format_ == Format::TECPLOT) {
        outputter_hb = new Outputter::HBTimeTecPlot(comm, output_manager, time_print_parameters);
      }
      else if ( (time_print_parameters.format_ == Format::RAW) ||
             (time_print_parameters.format_ == Format::RAW_ASCII) ||
             (time_print_parameters.format_ == Format::PROBE) ||
             (time_print_parameters.format_ == Format::TS1) ||
             (time_print_parameters.format_ == Format::TS2) )
      {
        Report::UserWarning0() << "HB_TD output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        time_print_parameters.format_ = Format::STD;
        outputter_hb = new Outputter::HBTimePrn(comm, output_manager, time_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB_TD output cannot be written in requested format, using standard format instead";
        time_print_parameters.format_ = Format::STD;
        outputter_hb = new Outputter::HBTimePrn(comm, output_manager, time_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB_TD, outputter_hb);
    }
  }

  // generate HB_IC outputters.  Note that this output uses the output classes for .PRINT TRAN output,
  // such as Outputter::TimePrn
  result = output_manager.findOutputParameter(OutputType::HB_IC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters hb_ic_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, hb_ic_print_parameters);

      Outputter::Interface *outputter_init;
      if (hb_ic_print_parameters.format_ == Format::STD) {
        outputter_init = new Outputter::HBICPrn(comm, output_manager, hb_ic_print_parameters);
      }
      else if (hb_ic_print_parameters.format_ == Format::CSV) {
        outputter_init = new Outputter::HBICCSV(comm, output_manager, hb_ic_print_parameters);
      }
      else if (hb_ic_print_parameters.format_ == Format::TECPLOT) {
        outputter_init = new Outputter::HBICTecplot(comm, output_manager, hb_ic_print_parameters);
      }
      else if ( (hb_ic_print_parameters.format_ == Format::RAW) ||
                (hb_ic_print_parameters.format_ == Format::RAW_ASCII) ||
                (hb_ic_print_parameters.format_ == Format::PROBE) ||
                (hb_ic_print_parameters.format_ == Format::TS1) ||
                (hb_ic_print_parameters.format_ == Format::TS2))

      {
        Report::UserWarning0() << "HB_IC output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        hb_ic_print_parameters.format_ = Format::STD;
        outputter_init = new Outputter::TimePrn(comm, output_manager, hb_ic_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB_IC output cannot be written in requested format, using standard format";
        hb_ic_print_parameters.format_ = Format::STD;
        outputter_init = new Outputter::HBICPrn(comm, output_manager, hb_ic_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB_IC, outputter_init);
    }
  }

  // generate HC_STARTUP outputters.   Note that this output uses the output classes for .PRINT TRAN output,
  // such as Outputter::TimePrn
  result = output_manager.findOutputParameter(OutputType::HB_STARTUP);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it) {
      PrintParameters hb_startup_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, hb_startup_print_parameters);

      Outputter::Interface *outputter_startup;
      if (hb_startup_print_parameters.format_ == Format::STD) {
        outputter_startup = new Outputter::TimePrn(comm, output_manager, hb_startup_print_parameters);
      }
      else if (hb_startup_print_parameters.format_ == Format::CSV) {
        outputter_startup = new Outputter::TimeCSV(comm, output_manager, hb_startup_print_parameters);
      }
      else if (hb_startup_print_parameters.format_ == Format::TECPLOT) {
        outputter_startup = new Outputter::TimeTecplot(comm, output_manager, hb_startup_print_parameters);
      }
      else if ( (hb_startup_print_parameters.format_ == Format::RAW) ||
                (hb_startup_print_parameters.format_ == Format::RAW_ASCII) ||
                (hb_startup_print_parameters.format_ == Format::PROBE) ||
                (hb_startup_print_parameters.format_ == Format::TS1) ||
                (hb_startup_print_parameters.format_ == Format::TS2))

      {
        Report::UserWarning0() << "HB_STARTUP output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        hb_startup_print_parameters.format_ = Format::STD;
        outputter_startup = new Outputter::TimePrn(comm, output_manager, hb_startup_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "HB_STARTUP output cannot be written in requested format, using standard format";
        hb_startup_print_parameters.format_ = Format::STD;
        outputter_startup = new Outputter::TimePrn(comm, output_manager, hb_startup_print_parameters);
      }

      output_manager.addOutputter(PrintType::HB_STARTUP, outputter_startup);
    }
  }

  // We now need to set up any external outputters, which will only be "HB_TD"
  // or "HB_FD," and can only do one type of output each.
  std::pair<ExternalOutputWrapperMap::const_iterator, bool> resultHB_FD = output_manager.findExternalOutputWrapper(OutputType::HB_FD);
  if  (resultHB_FD.second)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it = (*resultHB_FD.first).second.begin(), end = (*resultHB_FD.first).second.end(); it != end; ++it)
    {
      ExternalOutputWrapper * theWrapperPtr = (*it);
      output_manager.fixupOutputVariables(comm, theWrapperPtr->getParamList());
      Outputter::Interface *outputter;
      outputter = new Outputter::OutputterExternal(comm,output_manager,
                                                   theWrapperPtr);
      output_manager.addOutputter(PrintType::HB_FD, outputter);
    }
  }
  std::pair<ExternalOutputWrapperMap::const_iterator, bool> resultHB_TD = output_manager.findExternalOutputWrapper(OutputType::HB_TD);
  if  (resultHB_TD.second)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it = (*resultHB_TD.first).second.begin(), end = (*resultHB_TD.first).second.end(); it != end; ++it)
    {
      ExternalOutputWrapper * theWrapperPtr = (*it);
      output_manager.fixupOutputVariables(comm, theWrapperPtr->getParamList());
      Outputter::Interface *outputter;
      outputter = new Outputter::OutputterExternal(comm,output_manager,
                                                   theWrapperPtr);
      output_manager.addOutputter(PrintType::HB_TD, outputter);
    }
  }

  // We may also have HB_IC and HB_STARTUP,
  // which are just transient outputters at heart.

  std::pair<ExternalOutputWrapperMap::const_iterator, bool> resultHB_IC = output_manager.findExternalOutputWrapper(OutputType::HB_IC);
  if  (resultHB_IC.second)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it = (*resultHB_IC.first).second.begin(), end = (*resultHB_IC.first).second.end(); it != end; ++it)
    {
      ExternalOutputWrapper * theWrapperPtr = (*it);
      output_manager.fixupOutputVariables(comm, theWrapperPtr->getParamList());
      Outputter::Interface *outputter;
      outputter = new Outputter::OutputterExternal(comm,output_manager,
                                                   theWrapperPtr);
      output_manager.addOutputter(PrintType::HB_IC, outputter);
    }
  }
  std::pair<ExternalOutputWrapperMap::const_iterator, bool> resultHB_STARTUP = output_manager.findExternalOutputWrapper(OutputType::HB_STARTUP);
  if  (resultHB_STARTUP.second)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it = (*resultHB_STARTUP.first).second.begin(), end = (*resultHB_STARTUP.first).second.end(); it != end; ++it)
    {
      ExternalOutputWrapper * theWrapperPtr = (*it);
      output_manager.fixupOutputVariables(comm, theWrapperPtr->getParamList());
      Outputter::Interface *outputter;
      outputter = new Outputter::OutputterExternal(comm,output_manager,
                                                   theWrapperPtr);
      output_manager.addOutputter(PrintType::HB_STARTUP, outputter);
    }
  }

}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

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
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputterAC.h>
#include <N_IO_OutputterExternal.h>
#include <N_IO_OutputterFrequencyCSV.h>
#include <N_IO_OutputterFrequencyPrn.h>
#include <N_IO_OutputterFrequencyProbe.h>
#include <N_IO_OutputterFrequencyRaw.h>
#include <N_IO_OutputterFrequencyRawASCII.h>
#include <N_IO_OutputterFrequencyTecplot.h>
#include <N_IO_OutputterTimeCSV.h>
#include <N_IO_OutputterTimePrn.h>
#include <N_IO_OutputterTimeProbe.h>
#include <N_IO_OutputterTimeRaw.h>
#include <N_IO_OutputterTimeRawASCII.h>
#include <N_IO_OutputterTimeTecplot.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enableACOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Baur
// Creation Date : June 23, 2014
//-----------------------------------------------------------------------------
void enableACOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::AC_IC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it)
    {
      PrintParameters ac_ic_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, ac_ic_print_parameters);

      Outputter::Interface *outputter;
      if (ac_ic_print_parameters.format_ == Format::STD)
      {
        ac_ic_print_parameters.defaultExtension_ = ".TD.prn";
        outputter = new Outputter::TimePrn(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::CSV)
      {
        ac_ic_print_parameters.defaultExtension_ = ".TD.csv";
        outputter = new Outputter::TimeCSV(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::PROBE)
      {
        ac_ic_print_parameters.defaultExtension_ = ".TD.csd";
        outputter = new Outputter::TimeProbe(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::TECPLOT)
      {
        ac_ic_print_parameters.defaultExtension_ = ".TD.dat";
        outputter = new Outputter::TimeTecplot(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::RAW)
      {
        ac_ic_print_parameters.defaultExtension_ = ".raw";
        outputter = new Outputter::TimeRaw(comm, output_manager, ac_ic_print_parameters);
      }
      else if (ac_ic_print_parameters.format_ == Format::RAW_ASCII)
      {
        ac_ic_print_parameters.defaultExtension_ = ".raw";
        outputter = new Outputter::TimeRawAscii(comm, output_manager, ac_ic_print_parameters);
      }
      else if ( (ac_ic_print_parameters.format_ == Format::TS1) ||
                (ac_ic_print_parameters.format_ == Format::TS2) )
      {
        Report::UserWarning0() << "AC output cannot be written in Touchstone format, using standard format";
        ac_ic_print_parameters.format_ = Format::STD;
        outputter = new Outputter::TimePrn(comm, output_manager, ac_ic_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "AC output cannot be written in " << ac_ic_print_parameters.format_ << " format, using standard format";
        ac_ic_print_parameters.format_ = Format::STD;
        outputter = new Outputter::TimePrn(comm, output_manager, ac_ic_print_parameters);
      }

      output_manager.addOutputter(PrintType::AC_IC, outputter);
    }
  }

  result = output_manager.findOutputParameter(OutputType::AC);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end2 = (*result.first).second.end(); it != end2; ++it)
    {
      PrintParameters ac_print_parameters = (*it);

      output_manager.fixupPrintParameters(comm, ac_print_parameters);

      Outputter::Interface *outputter;
      if (ac_print_parameters.format_ == Format::STD)
      {
        outputter = new Outputter::FrequencyPrn(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::CSV)
      {
        outputter = new Outputter::FrequencyCSV(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::PROBE)
      {
        outputter = new Outputter::FrequencyProbe(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::TECPLOT)
      {
        outputter = new Outputter::FrequencyTecplot(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::RAW)
      {
          outputter = new Outputter::FrequencyRaw(comm, output_manager, ac_print_parameters);
      }
      else if (ac_print_parameters.format_ == Format::RAW_ASCII)
      {
        outputter = new Outputter::FrequencyRawAscii(comm, output_manager, ac_print_parameters);
      }
       else if ( (ac_print_parameters.format_ == Format::TS1) ||
                 (ac_print_parameters.format_ == Format::TS2) )
      {
        Report::UserWarning0() << "AC output cannot be written in Touchstone format, using standard format";
        ac_print_parameters.format_ = Format::STD;
        outputter = new Outputter::FrequencyPrn(comm, output_manager, ac_print_parameters);
      }
      else
      {
        Report::UserWarning0() << "AC output cannot be written in " << ac_print_parameters.format_ << " format, using standard format";
        ac_print_parameters.format_ = Format::STD;
        outputter = new Outputter::FrequencyPrn(comm, output_manager, ac_print_parameters);
      }

      output_manager.addOutputter(PrintType::AC, outputter);
    }
  }


  // Now do pretty much the same thing, but using the External output stuff
  // handle both AC_IC and AC in exactly the same way, because the outputter
  // knows how to do both
  std::pair<ExternalOutputWrapperMap::const_iterator, bool> result2 = output_manager.findExternalOutputWrapper(OutputType::AC_IC);
  if  (result2.second)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it = (*result2.first).second.begin(), end = (*result2.first).second.end(); it != end; ++it)
    {
      ExternalOutputWrapper * theWrapperPtr = (*it);
      output_manager.fixupOutputVariables(comm, theWrapperPtr->getParamList());
      Outputter::Interface *outputter;
      outputter = new Outputter::OutputterExternal(comm,output_manager,
                                                   theWrapperPtr);
      output_manager.addOutputter(PrintType::AC_IC, outputter);
    }
  }

  result2 = output_manager.findExternalOutputWrapper(OutputType::AC);
  if  (result2.second)
  {
    for (std::vector<ExternalOutputWrapper *>::const_iterator it = (*result2.first).second.begin(), end = (*result2.first).second.end(); it != end; ++it)
    {
      ExternalOutputWrapper * theWrapperPtr = (*it);
      output_manager.fixupOutputVariables(comm, theWrapperPtr->getParamList());
      Outputter::Interface *outputter;
      outputter = new Outputter::OutputterExternal(comm,output_manager,
                                                   theWrapperPtr);
      output_manager.addOutputter(PrintType::AC, outputter);
    }
  }

}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

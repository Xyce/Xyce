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
// Purpose        : Outputter class for PCE (Polynomial Chaos Expansion)
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 9/3/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_OutputterPCE.h>
#include <N_IO_OutputterPCEPrn.h>
#include <N_IO_OutputterPCECSV.h>
#include <N_IO_OutputterPCETecplot.h>
#include <N_IO_OutputMgr.h>

#include <N_ANP_UQSupport.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : enablePCEOutput
// Purpose       : Enable PCE output
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : September 3, 2019
//-----------------------------------------------------------------------------
void enablePCEOutput(Parallel::Machine comm, OutputMgr &output_manager, Analysis::Mode analysis_mode)
{
  std::pair<OutputParameterMap::const_iterator, bool> result = output_manager.findOutputParameter(OutputType::PCE);
  if (result.second)
  {
    for (std::vector<PrintParameters>::const_iterator it = (*result.first).second.begin(), end = (*result.first).second.end(); it != end; ++it)
    {
      PrintParameters pce_print_parameters = (*it);

      if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
        pce_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
      if (pce_print_parameters.printIndexColumn_)
        pce_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
      if (pce_print_parameters.printStepNumColumn_)
        pce_print_parameters.variableList_.push_front(Util::Param("STEPNUM", 0.0));

      output_manager.fixupPrintParameters(comm, pce_print_parameters);

      Outputter::Interface *outputter;
      if (pce_print_parameters.format_ == Format::STD)
      {
        pce_print_parameters.defaultExtension_ = ".PCE.prn";
        outputter = new Outputter::PCEPrn(comm, output_manager, pce_print_parameters);
      }
      else if (pce_print_parameters.format_ == Format::CSV)
      {
        pce_print_parameters.defaultExtension_ = ".PCE.csv";
        outputter = new Outputter::PCECSV(comm, output_manager, pce_print_parameters);
      }
      else if (pce_print_parameters.format_ == Format::TECPLOT)
      {
        pce_print_parameters.defaultExtension_ = ".PCE.dat";
        outputter = new Outputter::PCETecplot(comm, output_manager, pce_print_parameters);
      }
      else if ( (pce_print_parameters.format_ == Format::RAW) ||
                (pce_print_parameters.format_ == Format::RAW_ASCII) ||
                (pce_print_parameters.format_ == Format::PROBE) ||
                (pce_print_parameters.format_ == Format::TS1) ||
                (pce_print_parameters.format_ == Format::TS2))
      {
        Report::UserWarning0()
          << "PCE output cannot be written in PROBE, RAW or Touchstone format, using standard format instead";
        pce_print_parameters.format_ = Format::STD;
        outputter = new Outputter::PCEPrn(comm, output_manager, pce_print_parameters);
      }
      else
      {
        Report::UserWarning0()
          << "PCE output cannot be written in requested format, using standard format";
        pce_print_parameters.format_ = Format::STD;
        pce_print_parameters.defaultExtension_ = ".PCE.prn";
        outputter = new Outputter::PCEPrn(comm, output_manager, pce_print_parameters);
      }

      output_manager.addOutputter(PrintType::PCE, outputter);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : makePCEColumnNames
// Purpose       : This helps makes the strings used for the column headers
//                 for the PCEPrn and PCECSV classes.  It helps make the
//                 strings for the VARIABLES list for the PCETecplot class.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : Septmeber 6, 2019
//-----------------------------------------------------------------------------
void makePCEColumnNames(
  const PrintParameters&       printParameters,
  std::vector<std::string>&    colNames,
  int                          numQuadPoints,
  const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec)
{
  for (int iout=0;iout<outFuncDataVec.size();++iout)
  {
    Xyce::Analysis::UQ::outputFunctionData & outFunc = *(outFuncDataVec[iout]);

    if (printParameters.outputPCEsampleStats_)
    {
      colNames.push_back(outFunc.outFuncString + "_mean");
      colNames.push_back(outFunc.outFuncString + "_meanPlus");
      colNames.push_back(outFunc.outFuncString + "_meanMinus");

      colNames.push_back(outFunc.outFuncString + "_stddev");
      colNames.push_back(outFunc.outFuncString + "_variance");
    }

#ifdef Xyce_STOKHOS_ENABLE
    colNames.push_back(outFunc.outFuncString + "_quad_pce_mean");
    colNames.push_back(outFunc.outFuncString + "_quad_pce_meanPlus");
    colNames.push_back(outFunc.outFuncString + "_quad_pce_meanMinus");

    colNames.push_back(outFunc.outFuncString + "_quad_pce_stddev");
    colNames.push_back(outFunc.outFuncString + "_quad_pce_variance");
#endif

    if (printParameters.outputAllPCEsamples_)
    {
      for(int i=0;i<numQuadPoints; ++i)
      {
        colNames.push_back(outFunc.outFuncString + "_"+std::to_string(i));
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : outputPCEData
// Purpose       : This helps output data for everything but the STEPNUM,
//                 Index and Time columns/variables for the PCEPrn, PCECSV
//                 PCETecplot classes.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : September 6, 2019
//-----------------------------------------------------------------------------
void outputPCEData(
  const PrintParameters&       printParameters,
  std::ostream *               os,
  const std::vector<complex>&  result_list,
  int                          numQuadPoints,
  const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec)
{
  int colIdx = result_list.size();
  for (int iout=0;iout<outFuncDataVec.size();++iout)
  {
    Xyce::Analysis::UQ::outputFunctionData & outFunc = *(outFuncDataVec[iout]);

    if (printParameters.outputPCEsampleStats_)
    {
      printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, outFunc.sm.mean);
      colIdx++;
      printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx,outFunc.sm.mean+outFunc.sm.stddev);
      colIdx++;
      printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, outFunc.sm.mean-outFunc.sm.stddev);
      colIdx++;
      printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, outFunc.sm.stddev);
      colIdx++;
      printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, outFunc.sm.variance);
      colIdx++;
    }

#ifdef Xyce_STOKHOS_ENABLE
    Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > & projectionPCE = outFunc.projectionPCE;

    double pce_mean = projectionPCE.mean();
    double pce_stddev = projectionPCE.standard_deviation();
    double pce_variance = pce_stddev*pce_stddev;

    if ( std::isinf(pce_mean) || std::isnan(pce_mean) )
    {
      pce_mean = 0.0;
    }

    if ( std::isinf(pce_stddev) || std::isnan(pce_stddev) )
    {
      pce_stddev = 0.0;
    }

    if ( std::isinf(pce_variance) || std::isnan(pce_variance) )
    {
      pce_variance = 0.0;
    }

    printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, pce_mean);
    colIdx++;
    printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, pce_mean+pce_stddev);
    colIdx++;
    printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, pce_mean-pce_stddev);
    colIdx++;
    printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, pce_stddev);
    colIdx++;
    printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, pce_variance);
    colIdx++;
#endif

    // output individual samples
    if (printParameters.outputAllPCEsamples_)
    {
      for(int i=0;i<numQuadPoints; ++i, colIdx++)
      {
        printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, outFunc.sampleOutputs[i]);
      }
    }
  }

  return;
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

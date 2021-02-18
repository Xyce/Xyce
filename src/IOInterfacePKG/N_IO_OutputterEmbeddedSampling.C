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

#include <N_ANP_UQSupport.h>

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

//-----------------------------------------------------------------------------
// Function      : makeEmbeddedSamplingColumnNames
// Purpose       : This helps makes the strings used for the column headers
//                 for the EmbeddedSamplingPrn and EmbeddedSamplingCSV classes.
//                 It helps make the strings for the VARIABLES list for the
//                 EmbeddedSamplingTecplot class.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : September 6, 2019
//-----------------------------------------------------------------------------
void makeEmbeddedSamplingColumnNames(
  const PrintParameters&       printParameters,
  std::vector<std::string>&    colNames,
  bool                         regressionPCEenable,
  bool                         projectionPCEenable,
  int                          numSamples,
  const std::vector<std::string> & regressionPCEcoeffs,
  const std::vector<std::string> & projectionPCEcoeffs,
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
    if (regressionPCEenable)
    {
      colNames.push_back(outFunc.outFuncString + "_regr_pce_mean");
      colNames.push_back(outFunc.outFuncString + "_regr_pce_meanPlus");
      colNames.push_back(outFunc.outFuncString + "_regr_pce_meanMinus");

      colNames.push_back(outFunc.outFuncString + "_regr_pce_stddev");
      colNames.push_back(outFunc.outFuncString + "_regr_pce_variance");

      if (printParameters.outputPCECoeffs_)
      {
        std::vector<std::string>::const_iterator it;
        for (it=regressionPCEcoeffs.begin();it!=regressionPCEcoeffs.end();++it)
        {
          colNames.push_back(outFunc.outFuncString + *it);
        }
      }
    }

    if (projectionPCEenable)
    {
      colNames.push_back(outFunc.outFuncString + "_quad_pce_mean");
      colNames.push_back(outFunc.outFuncString + "_quad_pce_meanPlus");
      colNames.push_back(outFunc.outFuncString + "_quad_pce_meanMinus");

      colNames.push_back(outFunc.outFuncString + "_quad_pce_stddev");
      colNames.push_back(outFunc.outFuncString + "_quad_pce_variance");

      if (printParameters.outputPCECoeffs_)
      {
        std::vector<std::string>::const_iterator it;
        for (it=projectionPCEcoeffs.begin();it!=projectionPCEcoeffs.end();++it)
        {
          colNames.push_back(outFunc.outFuncString + *it);
        }
      }
    }
#endif

    if (printParameters.outputAllPCEsamples_)
    {
      for(int i=0;i<numSamples; ++i)
      {
        colNames.push_back(outFunc.outFuncString + "_"+std::to_string(i));
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : outputEmbeddedSamplingData
// Purpose       : This helps output data for everything but the STEPNUM,
//                 Index and Time columns/variables for the EmbeddedSamplingPrn,
//                 EmbeddedSamplingCSV and EmbeddedSamplingTecplot classes.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : September 6, 2019
//-----------------------------------------------------------------------------
void outputEmbeddedSamplingData(
  const PrintParameters&       printParameters,
  std::ostream *               os,
  const std::vector<complex>&  result_list,
  bool                         regressionPCEenable,
  bool                         projectionPCEenable,
  int                          numSamples,
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
    if (regressionPCEenable)
    {
      Stokhos::OrthogPolyApprox<int,double> & regressionPCE = outFunc.regressionPCE;

      double pce_mean = regressionPCE.mean();
      double pce_stddev = regressionPCE.standard_deviation();
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

      if (printParameters.outputPCECoeffs_)
      {
        int NN=regressionPCE.size();
        for (int ii=0;ii<NN;ii++, colIdx++)
        {
          printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, regressionPCE[ii]);
        }
      }
    }

    if (projectionPCEenable)
    {
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

      if (printParameters.outputPCECoeffs_)
      {
        int NN=projectionPCE.size();
        for (int ii=0;ii<NN;ii++, colIdx++)
        {
          printValue(*os, printParameters.table_.columnList_[colIdx], printParameters.delimiter_, colIdx, projectionPCE.fastAccessCoeff(ii));
        }
      }
    }
#endif

    // output individual samples
    if (printParameters.outputAllPCEsamples_)
    {
      for(int i=0;i<numSamples; ++i, colIdx++)
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

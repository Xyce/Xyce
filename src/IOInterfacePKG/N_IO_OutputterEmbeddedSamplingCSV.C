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
// Purpose        : Outputter class for csv files for Embedded Sampling info.
//
// Special Notes  :
//
// Creator        : Pete Sholander
//
// Creation Date  : 4/23/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterEmbeddedSamplingCSV.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::EmbeddedSamplingCSV
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
EmbeddedSamplingCSV::EmbeddedSamplingCSV(
    Parallel::Machine comm,
    OutputMgr &output_manager,
    const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0)
{
  if (printParameters_.defaultExtension_.empty())
  {
    printParameters_.defaultExtension_ = ".ES.csv";
  }

  // Add columns for Index and TIME.  The columns for the PCE info will be added in
  // doOutputEmbeddedSampling().  Any variables on the .PRINT ES line will be ignored.
  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);

}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::~EmbeddedSamplingCSV
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
EmbeddedSamplingCSV::~EmbeddedSamplingCSV()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::EmbeddedSamplingHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingCSV::EmbeddedSamplingHeader()
{
  int column_index = 0;

  Table::ColumnList::const_iterator it = printParameters_.table_.columnList_.begin();
  Table::ColumnList::const_iterator end = printParameters_.table_.columnList_.end();
  for ( ; it != end; ++it, ++column_index)
  {
    if (it != printParameters_.table_.columnList_.begin())
    {
      *os_ << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
    }
    printHeader(*os_, (*it));
  }

  *os_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::doOutputEmbeddedSampling
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingCSV::doOutputEmbeddedSampling(
  Parallel::Machine comm,
  bool              regressionPCEenable,
  bool              projectionPCEenable,
  int               numSamples,
  const std::vector<std::string> & regressionPCEcoeffs,
  const std::vector<std::string> & projectionPCEcoeffs,
  const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec)
{
  if (Parallel::rank(comm) == 0 && !os_)
  {
    outFilename_ = outputFilename(printParameters_.filename_,
                                  printParameters_.defaultExtension_,
                                  printParameters_.suffix_+outputManager_.getFilenameSuffix(),
                                  outputManager_.getNetlistFilename(),
                                  printParameters_.overrideRawFilename_,
                                  printParameters_.formatSupportsOverrideRaw_,
                                  printParameters_.dashoFilename_,
                                  printParameters_.fallback_);

    os_ = outputManager_.openFile(outFilename_);

    if (outputManager_.getPrintHeader())
    {
      // Generate names for the  additional header columns needed for PCE output
      std::vector<std::string> colNames;

      for (int iout=0;iout<outFuncDataVec.size();++iout)
      {
	Xyce::Analysis::UQ::outputFunctionData & outFunc = *(outFuncDataVec[iout]);

        if (printParameters_.outputPCEsampleStats_)
        {
          colNames.push_back(outFunc.outFuncString + "_mean");
          colNames.push_back(outFunc.outFuncString + "_meanPlus");
          colNames.push_back(outFunc.outFuncString + "_meanMinus");

          colNames.push_back(outFunc.outFuncString + "_stddev");
          colNames.push_back(outFunc.outFuncString + "_variance");
        }

#if Xyce_STOKHOS_ENABLE
        if (regressionPCEenable)
        {
          colNames.push_back(outFunc.outFuncString + "_regr_pce_mean");
          colNames.push_back(outFunc.outFuncString + "_regr_pce_meanPlus");
          colNames.push_back(outFunc.outFuncString + "_regr_pce_meanMinus");

          colNames.push_back(outFunc.outFuncString + "_regr_pce_stddev");
          colNames.push_back(outFunc.outFuncString + "_regr_pce_variance");

          if (printParameters_.outputPCECoeffs_)
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

          if (printParameters_.outputPCECoeffs_)
          {
            std::vector<std::string>::const_iterator it;
            for (it=projectionPCEcoeffs.begin();it!=projectionPCEcoeffs.end();++it)
            {
              colNames.push_back(outFunc.outFuncString + *it);
            }
          }
        }
#endif
        if (printParameters_.outputAllPCEsamples_)
        {
          for(int i=0;i<numSamples; ++i)
          {
#if __cplusplus>=201103L
            colNames.push_back(outFunc.outFuncString + "_"+std::to_string(i));
#else
            std::stringstream ss;
            ss << i;
            colNames.push_back(outFunc.outFuncString + "_"+ss.str());
#endif
          }
        }
      }

      // add those additional columns to the printParameters_.table_
      fixupColumnsFromStrVec(comm, printParameters_, colNames);

      // output the column names to the output file.
      EmbeddedSamplingHeader();
    }
  }

  std::vector<complex> result_list;
  getValues(comm, opList_, Util::Op::OpData(index_, 0, 0, 0, 0, 0), result_list);

  if (os_)
  {
    // Output the Index and TIME values.
    for (int i = 0; i < result_list.size(); ++i)
    {
      if (os_)
        printValue(*os_, printParameters_.table_.columnList_[i], printParameters_.delimiter_, i, result_list[i].real());
    }

    // Now output the PCE values directly from the the elements of the outFuncDataVec.
    // The next output column starts after the columns generated from the opList_
    int colIdx = result_list.size();
    for (int i=0; i< outFuncDataVec.size(); ++i)
    {
      Xyce::Analysis::UQ::outputFunctionData & outFunc = *(outFuncDataVec[i]);

      if (printParameters_.outputPCEsampleStats_)
      {
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, outFunc.sm.mean);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx,outFunc.sm.mean+outFunc.sm.stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, outFunc.sm.mean-outFunc.sm.stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, outFunc.sm.stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, outFunc.sm.variance);
        colIdx++;
      }

 #if Xyce_STOKHOS_ENABLE
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

        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_mean);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_mean+pce_stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_mean-pce_stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_variance);
        colIdx++;

        if (printParameters_.outputPCECoeffs_)
        {
          int NN=regressionPCE.size();
          for (int ii=0;ii<NN;ii++, colIdx++)
          {
            printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, regressionPCE[ii]);
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

        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_mean);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_mean+pce_stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_mean-pce_stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_stddev);
        colIdx++;
        printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, pce_variance);
        colIdx++;

        if (printParameters_.outputPCECoeffs_)
        {
          int NN=projectionPCE.size();
          for (int ii=0;ii<NN;ii++, colIdx++)
          {
            printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, projectionPCE.fastAccessCoeff(ii));
          }
        }
      }
#endif

      // output individual samples
      if (printParameters_.outputAllPCEsamples_)
      {
        for(int i=0;i<numSamples; ++i, colIdx++)
        {
          printValue(*os_, printParameters_.table_.columnList_[colIdx], printParameters_.delimiter_, colIdx, outFunc.sampleOutputs[i]);
        }
      }
    }

    // send end-of-line character
    *os_ << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::doFinishOutput
// Purpose       : Close the stream if there is no .STEP loop.  This function
//                 is also called after each step, if there is a .STEP loop,
//                 but currently does nothing in that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingCSV::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingCSV::doStartStep( int current_step, int number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingCSV::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingCSV::doSteppingComplete
// Purpose       : Close the stream when a .STEP loop is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingCSV::doSteppingComplete()
{
  // close the sensitivity file.
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

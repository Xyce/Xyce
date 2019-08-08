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
// Purpose        : Outputter class for Embedded Sampling info.
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

#include <N_IO_OutputterEmbeddedSamplingTecplot.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_IO_Tecplot.h>
#include <N_UTL_DeleteList.h>
#include <N_UTL_SaveIOSState.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::EmbeddedSamplingTecplot
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
EmbeddedSamplingTecplot::EmbeddedSamplingTecplot(
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
    printParameters_.defaultExtension_ = ".ES.dat";
  }

  // Add columns for Index and TIME.  The columns for the PCE info will be added in
  // doOutputEmbeddedSampling().  Any variables on the .PRINT ES line will be ignored.
  fixupColumns(comm, outputManager_.getOpBuilderManager(), printParameters_, opList_);

}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::~EmbeddedSamplingTecplot
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
EmbeddedSamplingTecplot::~EmbeddedSamplingTecplot()
{
  outputManager_.closeFile(os_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::EmbeddedSamplingHeader
// Purpose       :
// Special Notes : Adapted from tecplotFreqHeader function in N_IO_Tecplot.C
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::EmbeddedSamplingHeader()
{
  std::ostream &os = *os_;
  bool print_title =  (currentStep_ == 0);
  const std::string title = outputManager_.getNetlistFilename() + " - " + outputManager_.getTitle();
  const Table::ColumnList & columnList_ = printParameters_.table_.columnList_;

  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);
  os.precision(2);

  if (print_title)
  {
    std::string localTitle(title);
    std::size_t loc=0;
    while( (loc = localTitle.find_first_of( '"',loc ) ) != std::string::npos )
    {
      localTitle.insert(loc, 1, '\\');
      loc+=2;
    }
    os << " TITLE = \"" << localTitle << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    Table::ColumnList::const_iterator it2 = columnList_.begin();
    Table::ColumnList::const_iterator end2 = columnList_.end();
    for ( ; it2 != end2; ++it2)
    {
      os << "\" ";
      if (it2 != columnList_.begin())
      {
        *os_ << printParameters_.delimiter_;
      }
      printHeader(*os_, (*it2));
      os << "\" " << std::endl;
    }

    os << "DATASETAUXDATA ";
    os << getTecplotTimeDateStamp();
    os << std::endl;

    if (!outputManager_.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << outputManager_.getCircuitTemp() << " \"" << std::endl;
    }

    // this outputs the AuxData for Step 0
    outputAuxData(os);
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::outputAuxData
// Purpose       : output some AUXDATA
// Special Notes : Adapted from tecplotFreqHeader function in N_IO_Tecplot.C
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::outputAuxData(std::ostream &os)
{
  os << "ZONE F=POINT  ";

  if (outputManager_.getStepSweepVector().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<Analysis::SweepParam>::const_iterator it = outputManager_.getStepSweepVector().begin(); it != outputManager_.getStepSweepVector().end(); ++it)
    {
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  // put in the various sweep parameters as auxdata:
  if (!outputManager_.getStepSweepVector().empty())
  {
    for (std::vector<Analysis::SweepParam>::const_iterator iterParam = outputManager_.getStepSweepVector().begin();
    iterParam != outputManager_.getStepSweepVector().end();
    ++iterParam)
    {
      // convert any ":", "%" or "!" in the name to a "_", so as not to confuse tecplot.
      std::string tmpName(iterParam->name);
      replace(tmpName.begin(), tmpName.end(), '%', '_');
      replace(tmpName.begin(), tmpName.end(), ':', '_');
      replace(tmpName.begin(), tmpName.end(), '!', '_');
      os << "AUXDATA " << tmpName << " = " << "\" " << iterParam->currentVal << "\" ";
    }
    os << std::endl;
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::doOutputEmbeddedSampling
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::doOutputEmbeddedSampling(
  Parallel::Machine comm,
  bool              regressionPCEenable,
  bool              projectionPCEenable,
  int               numSamples,
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
        }

        if (projectionPCEenable)
        {
          colNames.push_back(outFunc.outFuncString + "_quad_pce_mean");
          colNames.push_back(outFunc.outFuncString + "_quad_pce_meanPlus");
          colNames.push_back(outFunc.outFuncString + "_quad_pce_meanMinus");

          colNames.push_back(outFunc.outFuncString + "_quad_pce_stddev");
          colNames.push_back(outFunc.outFuncString + "_quad_pce_variance");
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

    // Now output the PCE values directly from the elements of the outFuncDataVec.
    // The next output column starts after the columns generated from the opList_
    int colIdx = result_list.size();
        for (int iout=0;iout<outFuncDataVec.size();++iout)
    {
      Xyce::Analysis::UQ::outputFunctionData & outFunc = *(outFuncDataVec[iout]);

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
// Function      : EmbeddedSamplingTecplot::doFinishOutput
// Purpose       : Output the footer, and close the stream if there is no
//               : .STEP loop.  This function is also called after each step,
//               : if there is a .STEP loop, but currently does nothing in
//               : that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      if (outputManager_.getPrintFooter())
      {
        // this end-of-simulation footer is used if there is no .STEP loop
        (*os_) << "End of Xyce(TM) Embedded Sampling Simulation" << std::endl;
      }

      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes : This output the AUXDATA for steps 1,2,...)
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/26/2019
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::doStartStep( int current_step, int number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;

  // Output the AUXDATA if this is not Step 0.  The AUXDATA for Step 0
  // was output in the function EmbeddedSamplingHeader()
  if (os_ && current_step != 0)
  {
    std::ostream &os = *os_;
    outputAuxData(os);
  }
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : EmbeddedSamplingTecplot::doSteppingComplete
// Purpose       : Output footer and close the stream  when a .STEP loop
//               : is used.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void EmbeddedSamplingTecplot::doSteppingComplete()
{
  // close the embedded sampling file.
  if (os_)
  {
    // this end-of-simulation footer is used if there is a .STEP loop
    if ( outputManager_.getPrintFooter ())
    {
      (*os_) << "End of Xyce(TM) Embedded Sampling Simulation" << std::endl;
    }

    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

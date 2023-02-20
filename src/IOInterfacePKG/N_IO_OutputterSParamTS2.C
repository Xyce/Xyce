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
// Purpose        : Outputter for .LIN analysis.  This handles the output of
//                  S-, Y- and Z-parameters in Touchstone2 format.
// Special Notes  :
//
// Creator        : Pete Sholander
//
// Creation Date  : 3/25/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_OutputterSParamTS2.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : SParamT2S::SParamTS2
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
SParamTS2::SParamTS2(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    os_(0),
    index_(0),
    currentStep_(0),
    numberOfSteps_(0),
    numPorts_(1)
{
  // note the number of ports (e.g. 2) will be inserted into the default
  // extension in doOutputSParams()
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".sp";

  // Note: the usual call to fixupColumns() is deferred to the function
  // doOutputSParams() because the column names are inferred from the
  // size of Sparams matrix rather than being specified on the .LIN line.
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::~SParamTS2
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
SParamTS2::~SParamTS2()
{
  outputManager_.closeFile(os_);
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::sparamHeader
// Purpose       :
// Special Notes : The Sparams matrix contains the data for the parameter type
//                 (S, Y or Z) requested on the associcated .LIN line.
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS2::sparamHeader(
  Parallel::Machine comm,
  double            numFreq,
  std::vector<double> & Z0sVec,
  const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams)
{
  if (os_ && currentStep_ == 0)
  {
    std::ostream &os = *os_;
    std::string dataFormatStr;
    if (printParameters_.dataFormat_ == DataFormat::MA)
      dataFormatStr = "MA";
    else if (printParameters_.dataFormat_ == DataFormat::DB)
      dataFormatStr = "DB";
    else
      dataFormatStr = "RI";

    // Determine how many port impedances we need to print out.  Only
    // print out one Z0 value on the option line if they are all the
    // same value.
    bool printAllImpedances = false;
    for (int i=0; (i < numPorts_ && printAllImpedances == false); ++i)
    {
      for (int j = i; j < numPorts_; ++j)
      {
        if (Z0sVec[j] != Z0sVec[i])
	{
	  printAllImpedances = true;
          break;
        }
      }
    }

    os << "[Version] 2.0" << std::endl;
    os << "# Hz " << printParameters_.RFparamType_ << " " << dataFormatStr << " R" << " "  << Z0sVec[0];
    if (printAllImpedances && numPorts_ > 1)
      for (int i=1; i < numPorts_; ++i)
        os << " "  << Z0sVec[i];
    os << std::endl;
    os << "[Number of Ports] " << numPorts_ << std::endl;

    // This line is only used for two-port data and specifies the ordering of
    // S12 and S21 data.  The other option (for TS2) is 21_12, which reverses
    // their order in the output.
    if (numPorts_ == 2) { os << "[Two-Port Data Order] 12_21" << std::endl;}

    os << "[Number of Frequencies] " << numFreq << std::endl;
    os << "[Reference]";
    for (int i=0; i < numPorts_; ++i)
      os << " " << Z0sVec[i];
    os << std::endl;
    os << "[Network Data]" << std::endl;

    // add comment line with column header names
    os << "!";
    printHeader(os, printParameters_);
  }
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::doOutputSParams
// Purpose       :
// Special Notes : The RFparams map contains pointers to the matrices that contain
//                 the S-, Y- and Z-parameters.  The Sparams matrix is then selected
//                 from that map based on parameter type (S, Y or Z) requested on
//                 the associated .LIN line.  So, that matrix may end up containing
//                 S-, Y- or Z-parameter data.
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS2::doOutputSParams(
  Parallel::Machine   comm,
  double              frequency,
  double              numFreq,
  std::vector<double> & Z0sVec,
  const Util::Op::RFparamsData & RFparams)
{
  // Get a reference to the correct matrix (e.g., Sparams, Yparams or Zparams) in
  // the RFparams map based on the RFTYPE parameter on the .LIN line.
  Util::Op::RFparamsData::const_iterator it;
  it = RFparams.find(printParameters_.RFparamType_);
  const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams = *it->second;

  if (Parallel::rank(comm) == 0 && !os_)
  {
    numPorts_ = Z0sVec.size();
    // add number of ports to default extension
    std::ostringstream de;
    de << ".s" << numPorts_ << "p";
    printParameters_.defaultExtension_ = de.str();
    outFilename_ = outputFilename(printParameters_.filename_,
                                  printParameters_.defaultExtension_,
                                  printParameters_.suffix_+outputManager_.getFilenameSuffix(),
                                  outputManager_.getNetlistFilename(),
                                  printParameters_.overrideRawFilename_,
                                  printParameters_.formatSupportsOverrideRaw_,
                                  printParameters_.dashoFilename_,
                                  printParameters_.fallback_);
    os_ = outputManager_.openFile(outFilename_);

    // make the column names, based on the size of the Sparams matrix, and then set
    // the formatting of the output columns.
    std::vector<std::string> colNames;
    colNames.push_back("Freq");

    // The RF parameter type requested on the .LIN line.  "S" is the default.
    std::string paramType(printParameters_.RFparamType_);

    std::ostringstream ss;
    for (int i = 0; i < Sparams.numRows(); ++i)
    {
      for (int j = 0; j < Sparams.numRows(); ++j)
      {
        if (printParameters_.dataFormat_ == DataFormat::MA)
        {
          // Port numbering in Touchstone output starts at 1.
          // The Sparam matrix's indexing starts at 0.
          ss.str("");
          ss << "mag" << paramType << i+1 << j+1;
          colNames.push_back(ss.str());

          ss.str("");
          ss << "ang" << paramType << i+1 << j+1;
          colNames.push_back(ss.str());
        }
        else if (printParameters_.dataFormat_ == DataFormat::DB)
	{
          // Port numbering in Touchstone output starts at 1.
          // The Sparam matrix's indexing starts at 0.
          ss.str("");
          ss << "magdb" << paramType << i+1 << j+1;
          colNames.push_back(ss.str());

          ss.str("");
          ss << "ang" << paramType << i+1 << j+1;
          colNames.push_back(ss.str());
        }
        else
        {
          // default to DataFormat::RI
          ss.str("");
          ss << "Re" <<paramType << i+1 << j+1;
          colNames.push_back(ss.str());

          ss.str("");
          ss << "Im" << paramType << i+1 << j+1;
          colNames.push_back(ss.str());
        }
      }
    }
    // This call allows the PRECISION and WIDTH qualifiers to be used on the .LIN
    // line, just like on a .PRINT line.
    fixupColumnsFromStrVec(comm, printParameters_, colNames);

    // make the header lines
    sparamHeader(comm, numFreq, Z0sVec, Sparams);
  }

  if (os_)
  {
    // output the data lines
    printValue(*os_, printParameters_.table_.columnList_[0], printParameters_.delimiter_, 0, frequency);

    int count=0;
    double val1, val2;
    for (int i = 0; i < Sparams.numRows(); ++i)
    {
      for (int j = 0; j < Sparams.numRows(); ++j)
      {
        count++;
        // calculate the pair of values that will be output, for each complex value, based on the requested format
        if (printParameters_.dataFormat_ == DataFormat::MA)
	{
          val1 = std::abs(Sparams(i, j));
          val2 = 180*(arg(Sparams(i, j)))/M_PI;
        }
        else if (printParameters_.dataFormat_ == DataFormat::DB)
	{
          val1 = 20.0*std::log10(std::abs(Sparams(i, j)));
          val2 = 180*(arg(Sparams(i, j)))/M_PI;
        }
        else
        {
          // default to RI
          val1 = Sparams(i, j).real();
          val2 = Sparams(i, j).imag();
        }
        printValue(*os_, printParameters_.table_.columnList_[count], printParameters_.delimiter_, count, val1);
        printValue(*os_, printParameters_.table_.columnList_[count+1], printParameters_.delimiter_, count+1, val2);
      }
    }

    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::doFinishOutput
// Purpose       : Close the stream if there is no .STEP loop.  This function
//               : is also called after each step, if there is a .STEP loop,
//               : but currently does nothing in that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS2::doFinishOutput()
{
  if (os_)
  {
    if (numberOfSteps_ == 0)
    {
      *os_ << "[End]" << std::endl;
      outputManager_.closeFile(os_);
      os_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS2::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS2::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : SParamTS2::doSteppingComplete
// Purpose       : Close the stream  when a .STEP loop is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void SParamTS2::doSteppingComplete()
{
  // close the file.
  if (os_)
  {
    *os_ << "[End]" << std::endl;
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

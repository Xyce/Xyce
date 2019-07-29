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
// Purpose        : Outputter for .LIN analysis.  This handles the output of
//                  S-, Y- and Z-parameters in Touchstone1 format.
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

#include <N_IO_Outputter.h>
#include <N_IO_OutputterSParamTS1.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_Op.h>
#include <N_UTL_DeleteList.h>

namespace Xyce {
namespace IO {
namespace Outputter {

//-----------------------------------------------------------------------------
// Function      : SParamTS1::SParamTS1
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
SParamTS1::SParamTS1(Parallel::Machine comm, OutputMgr &output_manager, const PrintParameters &print_parameters)
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
// Creator       : Pete Sholander
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
SParamTS1::~SParamTS1()
{
  outputManager_.closeFile(os_);
}

//-----------------------------------------------------------------------------
// Function      : SParamTS1::sparamHeader
// Purpose       :
// Special Notes : The Sparams matrix contains the data for the parameter type
//                 (S, Y or Z) requested on the associcated .LIN line.
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS1::sparamHeader(
  Parallel::Machine comm,
  double            numFreq,
  std::vector<double> & Z0sVec,
  const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams)
{
  if (os_ && currentStep_ == 0)
  {
    // variable for wrapping the columns in the output file
    int wrapCount = (numPorts_ <= 3) ? numPorts_ : 4;

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

    os << "# Hz " << printParameters_.RFparamType_ << " " << dataFormatStr << " R"  << " "  << Z0sVec[0];
    if (printAllImpedances && numPorts_ > 1)
      for (int i=1; i < numPorts_; ++i)
        os << " "  << Z0sVec[i];
    os << std::endl;

    // add optional comment lines that indicates what the data fields are.  Each comment
    // line starts with a ! character.
    os << "!";
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

      // Account for line wrapping if there are more than two ports.  The formatting matches
      // what is in the Touchstone v2 standard's examples, including a "blank header column" at
      // the start of header rows 2,3,4,...
      if (numPorts_ >= 3)
      {
	if (( (column_index !=0) && ( column_index%(2*wrapCount) == 0 ) && ( column_index != (numColumns_-1) ) ))
	  os << std::endl << "!" << std::setw(printParameters_.table_.columnList_[0].width_) << " " ;
      }
    }
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : SParamTS1::doOutputSParams
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
void SParamTS1::doOutputSParams(
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
    // extra column is the Freq column, and each entry in the Sparams matrix
    // is a complex number
    numColumns_ = 2*Sparams.numRows()*Sparams.numRows() + 1;
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

    if (numPorts_ == 2)
    {
      // Touchstone1 format has a special ordering for 2 ports.  It it column-row order
      if (printParameters_.dataFormat_ == DataFormat::MA)
      {
        colNames.push_back("mag" + paramType + "11");
        colNames.push_back("ang" + paramType + "11");
        colNames.push_back("mag" + paramType + "21");
        colNames.push_back("ang" + paramType + "21");
        colNames.push_back("mag" + paramType + "12");
        colNames.push_back("ang" + paramType + "12");
        colNames.push_back("mag" + paramType + "22");
        colNames.push_back("ang" + paramType + "22");
      }
      else if (printParameters_.dataFormat_ == DataFormat::DB)
      {
        colNames.push_back("magdb" + paramType + "11");
        colNames.push_back("ang" + paramType + "11");
        colNames.push_back("magdb" + paramType + "21");
        colNames.push_back("ang" + paramType + "21");
        colNames.push_back("magdb" + paramType + "12");
        colNames.push_back("ang" + paramType + "12");
        colNames.push_back("magdb" + paramType + "22");
        colNames.push_back("ang" + paramType + "22");
      }
      else
      {
        // default to DataFormat::RI
        colNames.push_back("Re" + paramType + "11");
        colNames.push_back("Im" + paramType + "11");
        colNames.push_back("Re" + paramType + "21");
        colNames.push_back("Im" + paramType + "21");
        colNames.push_back("Re" + paramType + "12");
        colNames.push_back("Im" + paramType + "12");
        colNames.push_back("Re" + paramType + "22");
        colNames.push_back("Im" + paramType + "22");
      }
    }
    else
    {
      // handle every other value of numPorts_, which use row-column order
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
            ss << "Re" << paramType << i+1 << j+1;
            colNames.push_back(ss.str());

            ss.str("");
            ss << "Im" << paramType << i+1 << j+1;
            colNames.push_back(ss.str());
          }
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

    // output the frequency column
    printValue(*os_, printParameters_.table_.columnList_[0], printParameters_.delimiter_, 0, frequency);

    int count=0;
    double val1, val2;
    if (numPorts_ == 2)
    {
      // Touchstone1 format has a special ordering for 2 ports.  So, this prints out the columns
      // in the Sparams matrix.
      for (int j = 0; j < Sparams.numRows(); ++j)
      {
        for (int i = 0; i < Sparams.numRows(); ++i)
        {
          count++;
          // calculate the pair of values that will be output, for each complex value, based on the requested format
          if (printParameters_.dataFormat_ == DataFormat::MA)
	  {
            val1 = abs(Sparams(i, j));
            val2 = 180*(arg(Sparams(i, j)))/M_PI;
          }
          else if (printParameters_.dataFormat_ == DataFormat::DB)
	  {
            val1 = 20.0*std::log10(abs(Sparams(i, j)));
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
    }
    else
    {
      // handle every other value of numPorts_.  This just prints out each row in the Sparams matrix.
      int wrapCount = (numPorts_ <= 3) ? numPorts_ : 4;
      for (int i = 0; i < Sparams.numRows(); ++i)
      {
        for (int j = 0; j < Sparams.numRows(); ++j)
        {
          count++;
          // calculate the pair of values that will be output, for each complex value, based on the requested format
          if (printParameters_.dataFormat_ == DataFormat::MA)
	  {
            val1 = abs(Sparams(i, j));
            val2 = 180*(arg(Sparams(i, j)))/M_PI;
          }
          else if (printParameters_.dataFormat_ == DataFormat::DB)
	  {
            val1 = 20.0*std::log10(abs(Sparams(i, j)));
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
          // Only allowed to have four data values on a line.  Also add a "spacer" equal to the width of
          // the Freq column to the start of the 2nd,3rd,4th, ... lines of data.
          if ( (count%wrapCount == 0) && ((2*count) != numColumns_-1) ) *os_ << std::endl <<
                   " " << std::setw(printParameters_.table_.columnList_[0].width_);
        }
      }
    }

    *os_ << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : SParamTS1::doFinishOutput
// Purpose       : Close the stream if there is no .STEP loop.  This function
//               : is also called after each step, if there is a .STEP loop,
//               : but currently does nothing in that case.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS1::doFinishOutput()
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
// Function      : SParamTS1::doStartStep
// Purpose       : This function is executed at the start of each step.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS1::doStartStep(
  int                           current_step,
  int                           number_of_steps)
{
  index_ = 0;
  currentStep_ = current_step;
  numberOfSteps_ = number_of_steps;
}

//-----------------------------------------------------------------------------
// Function      : SParamTS1::doResetIndex
// Purpose       : Reset the value for the Index column to zero
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS1::doResetIndex()
{
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : SParamTS1::doSteppingComplete
// Purpose       : Close the stream  when a .STEP loop is used.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void SParamTS1::doSteppingComplete()
{
  // close the file.
  if (os_)
  {
    outputManager_.closeFile(os_);
    os_ = 0;
  }
}

} // namespace Outputter
} // namespace IO
} // namespace Xyce

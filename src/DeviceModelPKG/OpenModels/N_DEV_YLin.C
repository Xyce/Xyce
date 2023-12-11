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

//----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 08/19/2019
//
//
//
//
//----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_DEV_YLin.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Message.h>
#include <N_IO_OutputFileBase.h>
#include <N_IO_OutputPrn.h>
#include <N_IO_ParsingHelpers.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_RFparams.h>
#include <N_UTL_Math.h>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>

//#include <N_UTL_Interpolators.h>

namespace Xyce {
namespace Device {
namespace YLin {

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::initializeJacobianStamp
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Ting Mei, SNL
// Creation Date : 08/19/2019
//-----------------------------------------------------------------------------

void Instance::initializeJacobianStamp()
{
  if (jacStamp.empty())
  {
    jacStamp.resize(numExtVars);

    for ( int i=0; i < numExtVars; i++)
    {    
      jacStamp[i].resize(numExtVars);

      for ( int j=0; j< numExtVars; j++)
        jacStamp[i][j] = j;

    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Traits::loadInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 08/19/2019
//-----------------------------------------------------------------------------
void Traits::loadInstanceParameters(ParametricData<YLin::Instance> &p)
{
  p.addPar("M", 1.0, &YLin::Instance::multiplicityFactor)
    .setUnit(U_NONE)
    .setDescription("Multiplicity Factor");

  p.addPar("DTEMP", 0.0, &YLin::Instance::dtemp)
    .setGivenMember(&YLin::Instance::dtempGiven)
    .setUnit(U_DEGC)
    .setDescription("For compatibility only. Parameter is NOT used");

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Traits::loadModelParameters
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Pete Sholander, SNL
// Creation Date : 8/19/2019
//-----------------------------------------------------------------------------
void Traits::loadModelParameters(ParametricData<YLin::Model> &p)
{
  // Create parameter definitions for parameter member variables
  p.addPar("TSTONEFILE", "", &YLin::Model::TSFileName_)
    .setGivenMember(&YLin::Model::TSFileNameGiven_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("Touchstone File Name");

  p.addPar("ISC_FD", false, &YLin::Model::IscFD_)
    .setUnit(U_LOGIC)
    .setCategory(CAT_NONE)
    .setDescription("Touchstone file contains frequency-domain short-circuit current data");

  p.addPar("ISC_TD_FILE", "", &YLin::Model::ISC_TD_FileName_)
    .setGivenMember(&YLin::Model::ISC_TD_FileNameGiven_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("ISC Time Domain File Name");

  p.addPar("ISC_TD_FILE_FORMAT", "STD", &YLin::Model::ISC_TD_FileFormat_)
    .setGivenMember(&YLin::Model::ISC_TD_FileFormatGiven_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("Format of ISC Time Domain File");

  p.addPar("INTERPOLATION", 1, &YLin::Model::interpolation_ )
//    .setGivenMember(&YLin::Model::interpGiven_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("Interpolation method");


  p.addPar("HIGHPASS",  1, &YLin::Model::extrapolationHigh_ )
//    .setGivenMember(&YLin::Model::interpGiven_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("method to extrapolate higher frequency points" );

  p.addPar("LOWPASS",  1, &YLin::Model::extrapolationLow_ )
//    .setGivenMember(&YLin::Model::interpGiven_)
    .setUnit(U_NONE)
    .setCategory(CAT_NONE)
    .setDescription("method to extrapolate lower frequency points");
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::Instance
// Purpose       : Instance constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/19/2019
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    R(0.0),
    multiplicityFactor(0.0),
    dtemp(0.0),
    dtempGiven(false),
    G(0.0),
    i0(0.0),
    li_Pos(-1),
    li_Neg(-1),
    li_branch_data(0),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ,
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0)
  // ,
#endif

{
  // Initialize DeviceInstance values
  numIntVars   = 0;    // Initialize number if internal nodes in DeviceInstance
  numExtVars   = instance_block.numExtVars;    // Initialize number if external nodes in DeviceInstance
  numStateVars = 0;    // Initialize number if state variables in DeviceInstance
  setNumStoreVars(0);  // Initialize number if store variables in DeviceInstance

  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

//  initializeJacobianStamp();

  // Set params to constant default values from parameter definition
  setDefaultParams();

  // Set params according to instance line and constant defaults from metadata
  setParams(instance_block.params);

  // Calculate any parameters specified as expressions
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors
  processParams();

  if (numExtVars != 2 * model_.numPorts_)
  {
    UserError(*this) << "The number of nodes, "<< numExtVars << ", must be twice the number of ports, " << model_.numPorts_;
  }

  initializeJacobianStamp();

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/19/2019
//-----------------------------------------------------------------------------
bool Instance::processParams()
{

  // THIS IS A HACK.  Set resistance to first element of Z0Vec_ populated
  // from model's Touchstone file.
  if (model_.Z0Vec_.size() == 0)
  {
    R = 1000.0;
    UserWarning(*this) << "Z0Vec_ is empty, setting to the default, " << R << " ohms";
  }
  else
  {
    R=model_.Z0Vec_[0];
  }

  // M must be non-negative
  if (multiplicityFactor <= 0)
  {
    UserError(*this) << "Multiplicity Factor (M) must be non-negative" << std::endl;
  }

  // Apply the multiplicityFactor (M) from the instance line
  double factor = 1.0/multiplicityFactor;

  if (R*factor != 0.0)
    G = 1.0/(R * factor);
  else
    G = 0.0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs(
  const std::vector<int> & intLIDVecRef,
  const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy the local ID lists.
  intLIDVec = intLIDVecRef;                           // Set the internal local IDs in DeviceInstance
  extLIDVec = extLIDVecRef;                           // Set the external local IDs in DeviceInstance

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    dout() << getName() << " LIDs"
      //<< Util::push << std::endl
      << std::endl
           << "li_Pos_ = " << li_Pos << std::endl
           << "li_Neg_ = " << li_Neg << std::endl
           //<< Util::pop << std::endl;
           << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}


//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::registerBranchDataLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec)
{
  // Let DeviceInstance do its work.
  DeviceInstance::registerJacLIDs(jacLIDVec);

  // Store offsets of the components of the Jacobian of this instance
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];

  LIDVec_ = jacLIDVec;

}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  i0 = (solVec[li_Pos]-solVec[li_Neg])*G;

  fVec[li_Pos] += i0;
  fVec[li_Neg] += -i0;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_data] = i0;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Instance::loadDAEdFdx
// Purpose       : Loads the F-vector contributions for a single
//                 YLin  instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;
  return true;
}


// Begin model functions

//-----------------------------------------------------------------------------
// Function      : Model::readTouchStoneFile
// Purpose       : This function reads in the first few lines of a Touchstone 1
//                 or Touchstone 2 formatted file.  After finding the first
//                 non-comment line that is not blank, processing then branches
//                 based on whether the file format appears to be Touchstone 1
//                 or Touchstone 2.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
bool Model::readTouchStoneFile()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Processing Touchstone input file " << TSFileName_
                 << " for model " << getName() << std::endl;
  }

  // where we are in the file, and whether parsing had an error
  int lineNum = 0;
  bool psuccess = false;

  // try to open the data file
  std::ifstream inputFile;

  // Error out if the user-specified file does not exist, cannot
  // be opened, or is a directory name rather than a file name.
  // See SON Bug 785 and SRN Bug 2100 for more details.
  if ( !(Util::checkIfValidFile(TSFileName_)) )
  {
    Report::UserError() << "Touchstone input file \"" << TSFileName_ << "\" for model "
                        << getName() << " could not be found.";
    return false;
  }

  // open the Touchstone formatted input file
  inputFile.open( TSFileName_.c_str(),std::ifstream::in);
  if( !inputFile.good() )
  {
    Report::UserError() << "Touchstone file \"" << TSFileName_ << "\" for model "
                        << getName() << " could not be opened.";
    return false;
  }

  // Start reading lines in the input file
  bool firstLineFound = false;  // found first line that is not blank or a comment
  bool voLineFound = false;     // potentially valid [Version] or Option line found
  ExtendedString aLine("");
  IO::TokenVector parsedLine;
  readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);

  while( !(inputFile.eof() || (aLine.substr(0,5) == "[END]") || (voLineFound)) )
  {
    // start processing with first non-blank or non-comment line
    if ( (aLine.size() > 0) && (aLine[0] != TSCommentChar_ ) )
    {
      // [Version] line must be the first non-blank and non-comment line in
      // a Touchstone 2 formatted file.  The Option line must be the first
      // such line in a Touchstone 1 file.
      if (firstLineFound)
      {
        Report::UserError() << "No valid [Version] or Option line in file " << TSFileName_
               << " for model " << getName() << " at line " << lineNum;
        return false;
      }

      if (aLine.substr(0,9) == "[VERSION]")
      {
	psuccess = readTouchStone2File(inputFile, aLine, lineNum);
        voLineFound = true;  // stops further processing of input file in this loop
      }
      else if (aLine[0] == '#')
      {
        psuccess = readTouchStone1File(inputFile, aLine, lineNum);
        voLineFound = true;
      }

      firstLineFound = true; // used for error handling
    }
    else
    {
      // skip over leading blank and comment lines
      readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
    }
  }

  inputFile.close();

  return psuccess;
}

//-----------------------------------------------------------------------------
// Function      : Model::readTouchStone2File
// Purpose       : This is a parser to read in a Touchstone 2 formatted file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
bool Model::readTouchStone2File(std::ifstream& inputFile, const ExtendedString& firstLine, int lineNum)
{
  // where we are in the file, and whether parsing had an error
  int numVersionLinesFound = 0;
  bool optionLineFound = false;
  int numPortsLinesFound = 0;
  int numTwoPortDataOrderLinesFound = 0;
  int numFreqLinesFound = 0;
  int numReferenceLinesFound = 0;
  int numMatrixFormatLinesFound = 0;
  int numNetworkDataLinesFound = 0;
  int numNetworkDataElementsFound = 0;

  bool skipReadNextLine = false;
  bool psuccess = true;

  // used to parse each line of the file
  ExtendedString aLine("");
  IO::TokenVector parsedLine;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Found [Version] line at lineNum " << lineNum << std::endl;
  }

  ++numVersionLinesFound;
  splitTouchStoneFileLine(firstLine,parsedLine);

  if ( parsedLine.size() < 2 )
  {
    Report::UserError() << "Invalid [Version] line in file " << TSFileName_
      << " for model " << getName() << " at line " << lineNum;
    return false;
  }
  else
  {
    if (parsedLine[1].string_ != "2.0")
    {
      Report::UserError() << "Unsupported [Version] " << parsedLine[1].string_ << " in file "
        << TSFileName_ << " for model " << getName() << " at line " << lineNum;
        return false;
    }
    else
    {
      TSVersion_ = parsedLine[1].string_;
    }
  }

  // skip over any comment and blank lines
  if (!inputFile.eof())
  {
    readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
  }
  while( (!inputFile.eof()) && isTouchStoneBlankOrCommentLine(aLine) )
  {
    readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
  }

  // now parse the Option line, which starts with #
  if (inputFile.eof() || (aLine[0] != '#') )
  {
    Report::UserError() << "Option line not found immediately after [Version] line in file "
      << TSFileName_ << " for model " << getName() << " at line " << lineNum;
      return false;
  }
  else
  {
    optionLineFound = processTouchStoneOptionLine(aLine, lineNum);
    if (!optionLineFound)
      return false;
  }

  while ( !(inputFile.eof() || (aLine.substr(0,5) == "[END]")) )
  {
    // subsequent Option lines are silently ignored, along with blank or comment lines
    if ( !isTouchStoneBlankCommentOrOptionLine(aLine) )
    {
      // There can only be one [Version] in a a Touchstone 2 formatted file.
      if (aLine.substr(0,9) == "[VERSION]")
      {
        Report::UserError() << "Invalid [Version] line in file " << TSFileName_
               << " for model " << getName() << " at line " << lineNum;
        return false;
      }
      else if (aLine.substr(0,17) == "[NUMBER OF PORTS]")
      {
        // This line is required and may only appear once.  It must have an
        // integer value > 0.
        ++numPortsLinesFound;
        splitTouchStoneFileLine(aLine,parsedLine);
        if ( (parsedLine.size() != 4) || (numPortsLinesFound > 1) )
        {
          Report::UserError() << "Invalid [Number of Ports] line in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          psuccess = false;
        }
        else
	{
          ExtendedString numPortStr(parsedLine[3].string_);
          numPorts_ = numPortStr.Value();
          if ( (numPorts_ < 1) || !numPortStr.isInt())
          {
            Report::UserError() << "Number of ports in file " << TSFileName_ << " for model "
              << getName() << " must be an integer > 0";
            psuccess = false;
          }

          // populate the Z0Vec_ vector for a default Option line.  This may be
          // overwritten later by info on the [Reference] line.
          if (defaultOptionLine_)
            for (int i=0; i < numPorts_; i++) {Z0Vec_.push_back(50.0);}

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Set numPorts_ = " << numPorts_ << " at lineNum " << lineNum << std::endl;
          }
        }
      }
      else if (aLine.substr(0,21) == "[TWO-PORT DATA ORDER]")
      {
        // This line is required if numPorts=2.  It is forbidden otherwise.
        // If required then it must appear after the [Number of Ports] line
        // and before the [Network Data] line.  It only has the allowed
        // string values "12_21" or "21_12".
        ++numTwoPortDataOrderLinesFound;
        splitTouchStoneFileLine(aLine,parsedLine);
        if ( (parsedLine.size() != 4 ) || (numTwoPortDataOrderLinesFound > 1) || (numPorts_ != 2) ||
             (numNetworkDataLinesFound > 0) )
        {
          Report::UserError() << "Invalid [Two-Port Data Order] line in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          psuccess = false;
        }
        else
	{
          twoPortDataOrder_ = parsedLine[3].string_;
          if ( !( (twoPortDataOrder_ == "12_21" ) || (twoPortDataOrder_ == "21_12" ) ))
          {
            Report::UserError() << "Invalid string " << twoPortDataOrder_ << " in file " << TSFileName_
	      << "for [Two-Port Data Order] for model " << getName();
            psuccess = false;
          }

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Set twoPortDataOrder_ = " << twoPortDataOrder_ << " at lineNum " << lineNum << std::endl;
          }
        }
      }
      else if (aLine.substr(0,23) ==  "[NUMBER OF FREQUENCIES]")
      {
        // This line is required and may only appear once.  It must have an
        // integer value > 0.  It must appear after the [Number of Ports] line
        // and before the [Network Data] line.
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "Found [Number of Frequencies] line at lineNum " << lineNum << std::endl;
        }

        ++numFreqLinesFound;
        splitTouchStoneFileLine(aLine,parsedLine);

        if ( (parsedLine.size() != 4 ) || (numFreqLinesFound > 1) || (numPorts_ < 1) ||
	     ( numNetworkDataLinesFound > 0) )
        {
          Report::UserError() << "Invalid [Number of Frequencies] line in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          psuccess = false;
        }
        else
	{
          ExtendedString numFreqStr(parsedLine[3].string_);
          numFreq_ = numFreqStr.Value();
          if ( (numFreq_ < 1) || !numFreqStr.isInt())
          {
            Report::UserError() << "Number of frequencies in file " << TSFileName_ << " for model "
              << getName() << " must be an integer > 0";
            psuccess = false;
          }

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Set numFreq_ = " << numFreq_ << " at lineNum " << lineNum << std::endl;
          }
        }
      }
      else if (aLine.substr(0,11) ==  "[REFERENCE]")
      {
        // this line type needs to look ahead at the next line
        skipReadNextLine = true;

        // This line is optional, but may only appear once if it included.
        // It must appear after the [Number of Ports] line and before the [Network Data] line.
        ++numReferenceLinesFound;
        splitTouchStoneFileLine(aLine,parsedLine);

        if ( (parsedLine.size() < 2 ) || (numReferenceLinesFound > 1) || (numPorts_ < 1) ||
             (numNetworkDataLinesFound > 0) )
        {
          Report::UserError() << "Invalid [Reference] line in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          return false;
        }
        else
	{
          Z0Vec_.clear();
          // Handle line 1 of [Reference] line block
          for (int i=1; i<parsedLine.size(); ++i)
	  {
            ExtendedString z0Str(parsedLine[i].string_);
            Z0Vec_.push_back(z0Str.Value());
          }

          // read next line to see if the [Reference] block spans multiple lines
          if (!inputFile.eof())
          {
	     readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
          }

          // this a continuation of the [Reference] block
          if (aLine[0] !='[')
	  {
            splitTouchStoneFileLine(aLine,parsedLine);
            for (int i=0; i<parsedLine.size(); ++i)
	    {
              ExtendedString z0Str(parsedLine[i].string_);
              Z0Vec_.push_back(z0Str.Value());
            }
            readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
          }

          if ( Z0Vec_.size() != numPorts_ )
          {
            Report::UserError() << "[Reference] line in file " << TSFileName_
              << " for model " << getName() << " lacks enough entries";
            psuccess = false;
          }

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
	    Xyce::dout() << "[Reference] values are:";
            for (std::vector<double>::const_iterator it=Z0Vec_.begin(); it!=Z0Vec_.end(); ++it)
	      Xyce::dout() << " " << *it;
	    Xyce::dout() << std::endl;
          }
        }
      }
      else if (aLine.substr(0,15) ==  "[MATRIX FORMAT]")
      {
        // This line is optional.  If included, it must appear after the [Number of Ports] line
        // and before the [Network Data] line.  It has the allowed string values of
        // "FULL", "LOWER" or "UPPER".  This parsing only supports "FULL" now.
        ++numMatrixFormatLinesFound;
        splitTouchStoneFileLine(aLine,parsedLine);

        if ( (parsedLine.size() != 3) || (numPorts_ < 1) || (numNetworkDataLinesFound > 0) )
        {
          Report::UserError() << "Invalid [Matrix Format] line(s) in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          return false;
        }
        else
	{
          ExtendedString tokenStr(parsedLine[2].string_);
          if (IscFD_ && (tokenStr != "FULL"))
	  {
            Report::UserError() << "Only [Matrix Format] = FULL is supported when ISC=TRUE is used for YLIN model";
	    return false;
          }
          else if (tokenStr == "FULL")
	  {
	    matrixFormat_ = MatrixFormat::FULL;
	  }
          else if (tokenStr == "UPPER")
	  {
	    matrixFormat_ = MatrixFormat::UPPER;
	  }
          else if (tokenStr == "LOWER")
          {
	    matrixFormat_ = MatrixFormat::LOWER;
	  }
          else
          {
            // standard Touchstone 2 file, without frequency-domain ISC data
            Report::UserError() << "File " << TSFileName_ << " for model " << getName()
               << " had invalid value on [Matrix Format] line";
            return false;
          }

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Found matrix format = " << matrixFormat_ << " at lineNum " << lineNum << std::endl;
          }
        }

      }
      else if (aLine.substr(0,23) ==  "[NETWORK DATA]")
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "Found [Network Data] line at lineNum " << lineNum << std::endl;
        }

        // this line type needs to look ahead at the next line
        skipReadNextLine = true;

        // This line is required and may only appear once.  It must have an
        // integer value > 0.  It must appear after the [Number of Ports] line, and
        // immediately before the actual network data.
        ++numNetworkDataLinesFound;
        int numDataLinesFound = 0;

        if ( (numPorts_ < 1) || (numFreq_ < 1 ) )
        {
          Report::UserError() << "Unable to parse [Network Data Section] in file " << TSFileName_
             << " for model " << getName() << ".  Number of ports or frequencies < 1";
          return false;
        }
        else if ( numNetworkDataLinesFound > 1 )
	{
           Report::UserError() << "Invalid [Network Data] line in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          return false;
        }

        // skip over any comment lines
        if (!inputFile.eof())
	{
	  readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
        }
        while( (!inputFile.eof()) && isTouchStoneBlankCommentOrOptionLine(aLine) )
	{
          readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
        }

        // Set expected number of data elements on first line, assuming "FULL", "UPPER" or
        // "LOWER" formats
        if (matrixFormat_ == MatrixFormat::FULL)
          expectedNumElementsPerNetworkDataLine_ = 2*(numPorts_*numPorts_) + 1;
        else
          expectedNumElementsPerNetworkDataLine_ = numPorts_*(numPorts_+1) + 1;

        // only format "FULL" is supported if frequency-domain ISC data is given
        if (IscFD_)
          expectedNumElementsPerNetworkDataLine_ += 2*numPorts_;

        parsedLine.clear();

	// now populate the matrix
	while ( (!inputFile.eof()) && ( aLine[0] !='[') )
	{
          // account for Touchstone 2 files where each row of network data may
          // be "wrapped" across several rows in the input file
          IO::TokenVector tempParsedLine;
          splitTouchStoneFileLine(aLine,tempParsedLine);
          numNetworkDataElementsFound += tempParsedLine.size();
          if ( tempParsedLine.size() <= expectedNumElementsPerNetworkDataLine_ )
	  {
            parsedLine.insert(parsedLine.end(),tempParsedLine.begin(),tempParsedLine.end());
          }

          // We have a complete line of network data if this conditional is true
          if (parsedLine.size() == expectedNumElementsPerNetworkDataLine_)
	  {
            ++numDataLinesFound;
            if ( !processTouchStoneNetworkDataLine(parsedLine) )
              return false;

            // clear parsedLine here, to support parsing network data that
            // is "wrapped" across multiple rows in the input file
            parsedLine.clear();
          }

          // read in next line
          readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
        }

        if ( numDataLinesFound != numFreq_)
        {
          Report::UserError() << "Number of lines in [Network Data] does not match [Number of Frequencies] "
	    << "in file " << TSFileName_ << " for model " << getName();
          psuccess = false;
        }
      }
    }

    // Some line types (like [Reference]) used look-ahead to already read the next line
    if (!skipReadNextLine)
      readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
    else
      skipReadNextLine=false;
  }

  if ( (Z0Vec_.size() == 0) && optionLineFound )
  {
    // Handle case where optional [Reference] line is missing and the Option
    // line did not have a R value.  Use the default R value of 50.
    for (int i=1; i<=numPorts_; ++i)
      Z0Vec_.push_back(50.0);
  }
  else if ( (Z0Vec_.size() == 1) && (numPorts_ > 1) )
  {
    // Duplicate the single R value from the Options line into the rest
    // of the vector.
    for (int i=2; i<=numPorts_; ++i)
      Z0Vec_.push_back(Z0Vec_[0]);
  }

  // Some final error checking that the Touchstone 2 file had all of the required
  // lines.
  if ( (numVersionLinesFound == 0) || !optionLineFound )
  {
    Report::UserError() << "File " << TSFileName_ << " for model " << getName()
      << " lacked valid [Version] and/or Option lines";
    psuccess = false;
  }

  if ( (numPortsLinesFound == 0) || (numPorts_ == 0) )
  {
    Report::UserError() << "No valid [Number of Ports] line found in file " << TSFileName_
      << " for model " << getName();
    psuccess = false;
  }

  if ( (numFreqLinesFound == 0) || (numFreq_ == 0) )
  {
    Report::UserError() << "No valid [Number of Frequencies] line found in file " << TSFileName_
      << " for model " << getName();
    psuccess = false;
  }

  if ( (numNetworkDataLinesFound == 0) || (inputNetworkDataVec_.size() == 0) )
  {
    Report::UserError() << "No valid [Network Data] block of lines found in file " << TSFileName_
      << " for model " << getName();
    psuccess = false;
  }

  if ( (numNetworkDataElementsFound != numFreq_*expectedNumElementsPerNetworkDataLine_) &&
       (numPorts_ != 0) && (numFreq_ != 0) )
  {
    Report::UserError() << "Incorrect number of entries in [Network Data] block found in file "
                        << TSFileName_ << " for model " << getName()
                        << ". Found " << numNetworkDataElementsFound
                        << ". Expected " << numFreq_*expectedNumElementsPerNetworkDataLine_;
    psuccess = false;
  }

  return psuccess;
}

//-----------------------------------------------------------------------------
// Function      : Model::readTouchStone1File
// Purpose       : This is a parser to read in a Touchstone 1 formatted file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 05/25/2020
//-----------------------------------------------------------------------------
bool Model::readTouchStone1File(std::ifstream& inputFile, const ExtendedString& firstLine, int lineNum)
{
  // assume this is a Touchstone 1 file, since the first non-comment line started with #
  TSVersion_="1.0";
  twoPortDataOrder_ = "21_12";

  // process Option line
  if ( !processTouchStoneOptionLine(firstLine, lineNum) )
    return false;
  if (defaultOptionLine_)
    Z0Vec_.push_back(50.0); // remaining port impedances are set below

  // read next line in file, and then begin parsing
  ExtendedString aLine("");
  IO::TokenVector tempParsedLine;
  IO::TokenVector parsedLine;

  // skip over comment and blank lines to find first network data line
  readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
  while (!inputFile.eof() && isTouchStoneBlankCommentOrOptionLine(aLine))
    readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);

  if (aLine[0] == '[')
  {
    // Touchstone 2 meta-data lines not allowed in Touchstone 1 format
    Report::UserError() << "Invalid Touchstone 1 format, or possible missing [Version] line in Touchstone 2 "
      << "format, in file " << TSFileName_ << " for model " << getName() << " at line " << lineNum;
    return false;
  }
  splitTouchStoneFileLine(aLine,parsedLine);

  // skip over comment and blank lines to find second network data line
  readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
  while (!inputFile.eof() && isTouchStoneBlankCommentOrOptionLine(aLine))
    readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);

  if (aLine[0] == '[')
  {
    // Touchstone 2 meta-data lines not allowed in Touchstone 1 format
    Report::UserError() << "Invalid Touchstone 1 format, or possible missing [Version] line in Touchstone 2 "
      << "format, in file " << TSFileName_ << " for model " << getName() << " at line " << lineNum;
    return false;
  }

  // loop over remaining lines
  while (!inputFile.eof())
  {
    if ( !isTouchStoneBlankCommentOrOptionLine(aLine) )
    {
      if (aLine[0] == '[')
      {
        // Touchstone 2 meta-data lines not allowed in Touchstone 1 format
        Report::UserError() << "Invalid Touchstone 1 format, or possible missing [Version] line in Touchstone 2 "
	   << "format, in file " << TSFileName_ << " for model " << getName() << " at line " << lineNum;
        return false;
      }

      splitTouchStoneFileLine(aLine,tempParsedLine);

      if (tempParsedLine.size()%2 == 0)
      {
        // lines with a even number of elements are "column wrapped".  So, add them to the
        // parsedLine.
        parsedLine.insert(parsedLine.end(),tempParsedLine.begin(),tempParsedLine.end());
      }
      else
      {
        // lines with an odd number of elements on them are the start of the next complete line
        // of network data. So, process the previous "unwrapped" line here.
        numFreq_++;

        // Can determine the number of ports, and populate the Z0 vector, once the first
	// complete data line has been found.
        if (numFreq_ == 1)
        {
          if ( !setVarsFromTouchStone1File(parsedLine) )
	    return false;
        }
        else if (parsedLine.size() != expectedNumElementsPerNetworkDataLine_)
	{
          Report::UserError() << "Invalid network data line in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          return false;
        }

        if ( !processTouchStoneNetworkDataLine(parsedLine) )
          return false;

        // start accumulating a new line of network data
        parsedLine.clear();
        parsedLine = tempParsedLine;
      }
    }

    readAndUpperCaseTouchStoneFileLine(inputFile,aLine,lineNum);
  }

  if (numFreq_ == 0)
  {
    Report::UserError() << "Unable to parse network data section in file " << TSFileName_
      << " for model " << getName() << ".  Number of frequencies < 1";
    return false;
  }

  // handle last complete line of network data, if there is one
  if (parsedLine.size() > 0)
  {
    numFreq_++;

    // Handle the pathological case of only one network data line
    if (numFreq_ == 1)
    {
      if ( !setVarsFromTouchStone1File(parsedLine) )
	return false;
    }
    else if (parsedLine.size() != expectedNumElementsPerNetworkDataLine_)
    {
      Report::UserError() << "Invalid network data line in file " << TSFileName_
        << " for model " << getName() << " at line " << lineNum;
     return false;
    }

    if ( !processTouchStoneNetworkDataLine(parsedLine) )
      return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::setVarsFromTouchStone1File
// Purpose       : Set model variables based on first complete line of network
//                 data in a Touchstone 1 file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
bool Model::setVarsFromTouchStone1File(const IO::TokenVector & parsedLine)
{
  // Assume that the first complete line of network data has the correct
  // number of elements on it.
  expectedNumElementsPerNetworkDataLine_ = parsedLine.size();

  if (!IscFD_)
    numPorts_ = sqrt(0.5*(parsedLine.size()-1));
  else
  {
    // solution of quadratic equation with a=b=2 and c=(1-parsedLine.size())
    numPorts_ = 0.25*(sqrt(4.0 + 8.0*(parsedLine.size()-1)) - 2.0);
  }

  if ( (numPorts_ < 1) || (std::floor(numPorts_) != numPorts_) )
  {
    // derived value for numPorts_ must be an integer >=1
    Report::UserError() << "Error determining number of ports from file " << TSFileName_
	      << "for model " << getName();
    return false;
  }

  // Populate Z0 vector also
  for (int i=1; i<numPorts_; ++i)
    Z0Vec_.push_back(Z0Vec_[0]);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::splitTouchStoneFileLine
// Purpose       : Handle in-line comments, and then split a line from a
//                 Touchstone 1 or Touchstone 2 formatted file into a TokenVector
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
void Model::splitTouchStoneFileLine(const ExtendedString& aLine, IO::TokenVector & parsedLine)
{
  // trim off any in-line comments
  std::string trimmedLine(aLine);
  size_t pos = trimmedLine.find(TSCommentChar_);
  if (pos != std::string::npos) trimmedLine.erase(trimmedLine.begin()+pos,trimmedLine.end());

  // now split the line into "tokens", and also trim out extra white space from those
  // tokens.
  IO::splitLineNoWS(trimmedLine,parsedLine);

  return;
}

//-----------------------------------------------------------------------------
// Function      : Model::readAndUpperCaseTouchStoneFileLine
// Purpose       : Read in a line from a Touchstone formatted file, upper-case
//                 it, remove any leading whitespace, and then also increment
//                 the lineNum counter.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
void Model::readAndUpperCaseTouchStoneFileLine(std::istream & inputFile, ExtendedString& aLine, int& lineNum)
{
  IO::readLine(inputFile,aLine);
  aLine.toUpper();

  // remove any leading whitespace
  const std::string nonid(" \t\n\r\0");
  size_t start=aLine.find_first_not_of(nonid);
  if (start == std::string::npos)
    aLine="";
  else
    aLine = aLine.substr(start);

  lineNum++;

  return;
}

//-----------------------------------------------------------------------------
// Function      : Model::isTouchStoneBlankOrCommentLine
// Purpose       :
// Special Notes : "Blank lines" have zero length, since any leading whitespace
//                 has been removed from aLine already.  So, the comment
//                 character will always be the first character on the line.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/08/2020
//-----------------------------------------------------------------------------
bool Model::isTouchStoneBlankOrCommentLine(const ExtendedString& aLine)
{
  return ( (aLine.size() == 0) || (aLine[0] == TSCommentChar_) );
}

//-----------------------------------------------------------------------------
// Function      : Model::isTouchStoneBlankCommentOptionLine
// Purpose       :
// Special Notes : Option lines begin with a single '#' character, since any
//                 leading whitespace has been removed from aLine already.
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 06/08/2020
//-----------------------------------------------------------------------------
bool Model::isTouchStoneBlankCommentOrOptionLine(const ExtendedString& aLine)
{
  return ( isTouchStoneBlankOrCommentLine(aLine) || (aLine[0] == '#') );
}

//-----------------------------------------------------------------------------
// Function      : Model::processTouchStoneOptionLine
// Purpose       : Process the option line read in from either a Touchstone 1
//                 or Touchstone 2 formatted file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 05/25/2020
//-----------------------------------------------------------------------------
bool Model::processTouchStoneOptionLine(const ExtendedString& aLine, int lineNum)
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Found Option line at lineNum " << lineNum << std::endl;
  }

  // Record if this is a default Option line with only a # character on it,
  // so that Z0Vec_ vector can be populated with the default of 50, for each
  // port, once the number of ports is known.  The other defaults are set
  // correctly in the Model::Model constructor
  if (aLine.size() == 1)
    defaultOptionLine_ = true;

  // The fields on the Option line can be in any order. Yuck!  Also an Option
  // line that only has a # is allowed.
  IO::TokenVector parsedLine;
  splitTouchStoneFileLine(aLine,parsedLine);

  for (int i=1; i<parsedLine.size(); ++i)
  {
    ExtendedString tokenStr(parsedLine[i].string_);
    if (tokenStr == "S")
    {
      paramType_ = ParamType::S;
    }
    else if (tokenStr == "Y")
    {
      paramType_ = ParamType::Y;
    }
    else if (tokenStr == "Z")
    {
      paramType_ = ParamType::Z;
    }
    else if (tokenStr == "RI")
    {
      dataFormat_ = DataFormat::RI;
    }
    else if (tokenStr == "MA")
    {
      dataFormat_ = DataFormat::MA;
    }
    else if (tokenStr == "DB")
    {
      dataFormat_ = DataFormat::DB;
    }
    else if (tokenStr == "HZ")
    {
      // There are four allowed frequency multipliers.  So, explicitly
      // check for each one and hard-code the conversion.
      freqMultiplier_ = 1.0;
    }
    else if (tokenStr == "KHZ")
    {
      freqMultiplier_ = 1.0e+3;
    }
    else if (tokenStr == "MHZ")
    {
      freqMultiplier_ = 1.0e+6;
    }
    else if (tokenStr == "GHZ")
    {
      freqMultiplier_ = 1.0e+9;
    }
    else if ( tokenStr == "R" )
    {
      if (i < parsedLine.size()-1)
      {
        ExtendedString RStr(parsedLine[i+1].string_);
        if (RStr.isValue())
        {
          Z0Vec_.push_back(RStr.Value());
        }
        else
        {
          Report::UserError() << "Invalid value for R on option line in file " << TSFileName_
            << " for model " << getName() << " at line " << lineNum;
          return false;
        }
      }
      else
      {
        Report::UserError() << "Missing value for R on option line in file " << TSFileName_
          << " for model " << getName() << " at line " << lineNum;
        return false;
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
     Xyce::dout() << "After processing Option line:" << std::endl
       << "  Parameter type is " << paramType_ << std::endl
       << "  Data format is " << dataFormat_ << std::endl
       << "  Frequency Mutliplier is " << freqMultiplier_ << std::endl
       << "  R value is " << Z0Vec_[0] << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::processTouchStoneNetworkDataLine
// Purpose       : Process a network data line read in from either a Touchstone 1
//                 or Touchstone 2 formatted file
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 05/25/2020
//-----------------------------------------------------------------------------
bool Model::processTouchStoneNetworkDataLine(const IO::TokenVector& parsedLine)
{
  std::vector<std::complex<double> > inputIscData;
  Teuchos::SerialDenseMatrix<int, std::complex<double> > inputNetworkData;
  inputNetworkData.shape(numPorts_, numPorts_);

  ExtendedString freqStr(parsedLine[0].string_);
  freqVec_.push_back(freqMultiplier_*freqStr.Value());

  // Offset into parsedLine, for the data values.  The frequency value is at offset=0
  int offset=1;

  // These values work, for FULL format, for the inner (column index) loop below
  int startIdx=0;
  int endIdx=numPorts_-1;

  // For 2-port networks, this assumes [Two-Port Data Order] = "12_21".  If it is
  // actually "21_12" then that is handled below.
  for (int i=0; i<=numPorts_-1; ++i)
  {
    // adjust the starting or ending column index, in the inner loop, for UPPER or
    // LOWER format respectively
    if (matrixFormat_ == MatrixFormat::UPPER)
       startIdx = i;
    else if (matrixFormat_ == MatrixFormat::LOWER)
       endIdx = i;

    for (int j=startIdx; j<=endIdx; ++j)
    {
      // data will be converted to RI format for internal use in
      // YLin model
      ExtendedString Str1(parsedLine[offset].string_);
      ExtendedString Str2(parsedLine[offset+1].string_);
      if (dataFormat_ == DataFormat::RI)
      {
        // Use inputNetworkData(i,j) format for accessing, in
        // order to use row-column indexing.  Also, indexing
        // starts at (0,0).
        inputNetworkData(i,j).real(Str1.Value());
        inputNetworkData(i,j).imag(Str2.Value());
      }
      else if (dataFormat_ == DataFormat::MA)
      {
        double mag = Str1.Value();
        double angle = M_PI*Str2.Value()/180.0;
        inputNetworkData(i,j).real(mag*cos(angle));
        inputNetworkData(i,j).imag(mag*sin(angle));
      }
      else if (dataFormat_ == DataFormat::DB)
      {
        double mag = pow(10.0,0.05*Str1.Value());
        double angle = M_PI*Str2.Value()/180.0;
        inputNetworkData(i,j).real(mag*cos(angle));
        inputNetworkData(i,j).imag(mag*sin(angle));
      }
      else
      {
        Report::UserError() << "Incorrect data format " << dataFormat_ << " for network data in file "
          << TSFileName_ << " for model " << getName();
        return false;
      }

      // move offset to start of next pair of real/imaginary values
      offset += 2;

      // internally, the data is stored in Full format.  So, populate the values
      // on the other side of the diagonal
      if ( (i!=j) && ((matrixFormat_ == MatrixFormat::UPPER) || (matrixFormat_ == MatrixFormat::LOWER)) )
        inputNetworkData(j,i) = inputNetworkData(i,j);
    }
  }

  if ( (twoPortDataOrder_ == "21_12") && (numPorts_ == 2) )
  {
    // indexing into inputNetworkData starts at (0,0)
    complex tempVal = inputNetworkData(0,1);
    inputNetworkData(0,1) = inputNetworkData(1,0);
    inputNetworkData(1,0) = tempVal;
  }

  // YLin model will use Y-parameter format internally
  if (paramType_== ParamType::S)
  {
    Teuchos::SerialDenseMatrix<int, std::complex<double> > YParams;
    YParams.shape(numPorts_, numPorts_);
    Util::stoy(inputNetworkData,YParams,Z0Vec_);
    inputNetworkDataVec_.push_back(YParams);
  }
  else if (paramType_== ParamType::Z)
  {
    Teuchos::SerialDenseMatrix<int, std::complex<double> > YParams;
    YParams.shape(numPorts_, numPorts_);
    Util::ztoy(inputNetworkData,YParams);
    inputNetworkDataVec_.push_back(YParams);
  }
  else
  {
    // input was in Y-parameter format
    inputNetworkDataVec_.push_back(inputNetworkData);
  }

  // Handle optional Isc data. There are two additional columns for
  // each port added to each row, if that data is included.  This code only
  // only works for Full matrix format.
  if (IscFD_)
  {
    inputIscData.clear();
    for (int i=2*numPorts_*numPorts_+1, pIdx=1; i < expectedNumElementsPerNetworkDataLine_; i=i+2, pIdx++)
    {
      ExtendedString Str1(parsedLine[i].string_);
      ExtendedString Str2(parsedLine[i+1].string_);
      if (dataFormat_ == DataFormat::RI)
      {
        inputIscData.push_back(std::complex<double>(Str1.Value(), Str2.Value()));
      }
      else if (dataFormat_ == DataFormat::MA)
      {
        double mag = Str1.Value();
        double angle = M_PI*Str2.Value()/180.0;
        inputIscData.push_back(std::complex<double>(mag*cos(angle), mag*sin(angle)));
      }
      else if (dataFormat_ == DataFormat::DB)
      {
        double mag = pow(10.0,0.05*Str1.Value());
        double angle = M_PI*Str2.Value()/180.0;
        inputIscData.push_back(std::complex<double>(mag*cos(angle), mag*sin(angle)));
      }

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "For Freq " << freqVec_[freqVec_.size()-1] << " Isc" << pIdx
                     << "=(" << inputIscData[pIdx-1].real() << ", "
                     << inputIscData[pIdx-1].imag() << ")" << std::endl;
      }

    }
    inputIscFDVec_.push_back(inputIscData);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::readISC_TD_File
// Purpose       : Read in file with time-domain short-circuit current data
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 04/10/2020
//-----------------------------------------------------------------------------
void Model::readISC_TD_File()
{
  // Some error checking on the file name specified by ISC_TD_FILE
  if ( !(Util::checkIfValidFile(ISC_TD_FileName_)) )
  {
    Report::UserError() << "ISC_TD_FILE \"" << ISC_TD_FileName_ << "\" for model "
                        << getName() << " could not be found.";
    return;
  }

  // check the ISC_TD_FILE_FORMAT model parameter, if given
  if ( ISC_TD_FileFormatGiven_ &&
       !((ISC_TD_FileFormat_ == "STD") || (ISC_TD_FileFormat_ == "NOINDEX") || (ISC_TD_FileFormat_ == "CSV")) )
  {
    Report::UserError() << "ISC_TD_FILE_FORMAT for model " << getName()
                        << " must be STD, NOINDEX or CSV.";
    return;
  }

  // Use the same functions for reading the .PRN and .CSV files as re-measure does.
  IO::OutputFileBase* inputFile;
  inputFile = new IO::OutputPrn();

  ExtendedString fileExt(ISC_TD_FileName_.substr(ISC_TD_FileName_.length()-4));
  int timeColumn;  // column with the time variable

  // Set the position of the TIME column, and also auto-sense .prn or .csv files
  // is the ISC_TD_FILE_FORMAT model parameter is not given.
  if (ISC_TD_FileFormatGiven_)
  {
    if (ISC_TD_FileFormat_ == "STD")
      timeColumn = 1;
    else if ( (ISC_TD_FileFormat_ == "CSV") || (ISC_TD_FileFormat_ == "NOINDEX") )
      timeColumn = 0;
  }
  else if (!ISC_TD_FileFormatGiven_)
  {
    if (fileExt.toUpper() == ".PRN")
    {
      ISC_TD_FileFormat_ = "STD";
      timeColumn = 1;
    }
    else if (fileExt.toUpper() == ".CSV")
    {
      ISC_TD_FileFormat_ = "CSV";
      timeColumn = 0;
    }
    else
    {
      Report::UserError() << "Problem determining ISC_TD_FILE_FORMAT for file \"" << ISC_TD_FileName_
                        << "\" for model " << getName();
      return;
    }
  }

  if (!inputFile->openFileForRead(ISC_TD_FileName_))
  {
    // open failed.  Report error and exit.
    Report::UserError() << "Could not open ISC_TD_FILE \"" << ISC_TD_FileName_ << "\" for model " << getName();
    return;
  }

  // load data-names, as a way of reading and discarding the first row in the input file
  std::vector<std::string> fileVarNames;
  if (!(inputFile->getOutputVarNames(fileVarNames)))
  {
    // reading var names failed
    Report::UserError() << "Problem reading variable names in ISC_TD_FILE \"" << ISC_TD_FileName_ << "\" for model " << getName();
    return;
  }

  // This code section run though the lines in the input file, and makes the
  // the vectors for the time values and per-port short-circuit current values.
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "For ISC_TD_FILE \"" << ISC_TD_FileName_ << "\" for model " << getName() << std::endl;
  }
  int reading=1;
  std::vector<double> varValuesVec;
  inputIscTDVec_.resize(numPorts_);

  while (reading==1)
  {
    reading = inputFile->getOutputNextVarValuesSerial(&varValuesVec);

    if( reading == 1 )
    {
      // number of entries on each line should be (number of ports)+1 (plus possibly the INDEX column)
      if ( ((ISC_TD_FileFormat_ == "STD") && (varValuesVec.size() != numPorts_+2)) ||
           (((ISC_TD_FileFormat_ == "CSV") || (ISC_TD_FileFormat_ == "NOINDEX")) && (varValuesVec.size() != numPorts_+1)) )
      {
        Report::UserError() << "Incorrect number of entries found in ISC_TD_FILE \"" << ISC_TD_FileName_ << "\" for model " << getName();
        return;
      }

      iscTDTimeVec_.push_back(varValuesVec[timeColumn]);
      // inputIscTDVec_ will have one vector for each port.  Each of those vectors contains the short-circuit
      // current values for each time point for that port.
      for (int i=0; i<numPorts_; i++)
        inputIscTDVec_[i].push_back(varValuesVec[timeColumn+1+i]);

      varValuesVec.clear();
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    for (int i=0; i<iscTDTimeVec_.size(); i++)
    {
      Xyce::dout() << "  For Time " <<  iscTDTimeVec_[i] << " Isc =";
      for (int j=0; j<numPorts_; j++)
        Xyce::dout() << " " << inputIscTDVec_[j][i];

      Xyce::dout() << std::endl;
    }
  }

  // clean up file object, which is local to this function
  inputFile->closeFileForRead();
  delete inputFile;
  inputFile = 0;

  return;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
///
/// Process model parameters
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   6/03/02
bool Model::processParams()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
///
/// Process the instance parameters of instance owned by this model
///
/// This method simply loops over all instances associated with this
/// model and calls their processParams method.
///
/// @return true
///
/// @author Dave Shirely, PSSI
/// @date   03/23/06

bool Model::processInstanceParams()
{
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    (*it)->processParams();
  }

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Model::N_DEV_YLinModel
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
///
/// Construct a YLin model from a "model block" that was created
/// by the netlist parser.
///
/// @param configuration
/// @param model_block
/// @param factory_block
///
/// @author Pete Sholander, SNL
/// @date   8/19/2019
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
    TSFileName_(""),
    TSFileNameGiven_(false),
    TSCommentChar_('!'),
    TSVersion_(""),
    freqUnit_("GHZ"),
    freqMultiplier_(1.0e9),
    paramType_(ParamType::S),
    matrixFormat_(MatrixFormat::FULL),
    dataFormat_(DataFormat::MA),
    numPorts_(0),
    twoPortDataOrder_(""),
    numFreq_(0),
    Z0Vec_(0),
    freqVec_(0),
    defaultOptionLine_(false),
    expectedNumElementsPerNetworkDataLine_(0),
    IscFD_(false), 
    interpolation_(1),
    extrapolationHigh_(1),
    extrapolationLow_(1),
    IscTD_(false),
    yInterpolator(0)
{
  // Set params to constant default values.
  setDefaultParams();

  // Set params according to .model line and constant defaults from metadata.
  setModParams(model_block.params);

  // Calculate any parameters specified as expressions.
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors.
  processParams();

  // read Touchstone 2 formatted input file
  bool TouchstoneFileRead=false;
  if (TSFileNameGiven_)
    TouchstoneFileRead = readTouchStoneFile();
  else
    UserError(*this) << "No Touchstone input file given for model " << getName();

  // read time-domain ISC data, if given
  if ( ISC_TD_FileNameGiven_ && (numPorts_ > 0) )
  {
    readISC_TD_File();
    IscTD_ = true;
  }

  // it the file was successfully converted it is in Y-parameters and RI format
  if (TouchstoneFileRead)
  {
    paramType_= ParamType::Y;
    dataFormat_= DataFormat::RI;
  }



  switch (interpolation_ )
  {
    case 1:
      yInterpolator = new Util::linear<double>();
      break;

    case 2:
      yInterpolator = new Util::akima<double>();
      break;

    default:
      UserFatal(*this) << "Unsupported interpolation method. ";
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Model::~N_DEV_YLinModel
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Destroy this model.
///
/// Also destroys all instances that use this model.
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   3/16/00
Model::~Model()
{
  // Destory all owned instances
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    delete (*it);
  }

  yInterpolator->clear();
  delete yInterpolator;

}

//-----------------------------------------------------------------------------
// Function      : N_DEV_YLinModel::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
///
/// Print instances associated with this model.
///
/// Used only for debugging
///
/// @param os output stream
///
/// @return reference to output stream
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   4/03/00
///
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  os << std::endl;
  os << "Number of YLin Instances: " << instanceContainer.size() << std::endl;
  os << "    name     model name  Parameters" << std::endl;

  int i = 0;
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    os << "  " << i << ": " << (*it)->getName() << "\t";
    os << getName();
    os << "\t\tR(Tnom) = " << (*it)->R;
    os << "\tG(T) = " << (*it)->G;
    os << std::endl;
    ++i;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
///
/// @param[in] op Operator to apply to all instances.
///
///
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}

             

 
void Model::interpLin( double freq,  Teuchos::SerialDenseMatrix<int, std::complex<double> > & result, std::vector<std::complex<double> >  & Iscvals  )
{
//  std::vector<std::complex<double> > 
//  result.resize ( freqData[0].yth_.size());

//  int fsize = freqData.size();

  // Perform linear interpolation or extrapolation

         
  if (Iscvals.empty())
    Iscvals.resize(numPorts_ );


  double fmin = freqVec_[0];
  double fmax = freqVec_[numFreq_ - 1];

  double f0;
  double f1;

  int n0;
  int n1;

  if (freq <= fmin)
  {
// Use fmin constant data
//    sparams = (ymat_[0]);

// if we want to use the same extra/interpolation
//    n0 = 0;
//    n1 = 1;

    result = inputNetworkDataVec_[0];

    if (IscFD_)
      Iscvals = inputIscFDVec_[0];

    return;
  }
  else
  {
    if (freq >= fmax)
    {
//      n0 = fsize - 2;
//      n1 = fsize - 1;

      result = inputNetworkDataVec_[numFreq_ - 1];


      if (IscFD_)
        Iscvals = inputIscFDVec_[numFreq_ - 1];

      return;

    }
    else
    {
      n0 = 0;

      for (int n = 0; n < numFreq_; n++)
        if ( freqVec_[n] >=freq)
	{
	  n0 = n - 1; break;
	}

      n1 = n0 + 1;
   
    }
  }


  f0 = freqVec_[n0];
  f1 = freqVec_[n1];

  std::complex<double> y0, y1;

  for (int i=0; i< numPorts_ ; i++)
  {

    for (int j=0; j< numPorts_ ; j++)
    {
      y0 = (inputNetworkDataVec_[n0])(i, j);
      y1 = (inputNetworkDataVec_[n1])(i, j);
      result(i, j)= y0 + (y1 - y0)*(freq - f0) / (f1 - f0);

    }

    if (IscFD_)
    {
      y0 = inputIscFDVec_[n0][i];
      y1 = inputIscFDVec_[n1][i];
      Iscvals[i] = y0 + (y1 - y0)*(freq - f0) / (f1 - f0);
    }
  }

}


//-----------------------------------------------------------------------------
// Function      : Model::interpData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 8/2023
//----------------------------------------------------------------------------
void Model::interpData( double freq,  Teuchos::SerialDenseMatrix<int, std::complex<double> > & result, std::vector<std::complex<double> >  & Iscvals  )
{

// Perform interpolation or extrapolation

         
  if (Iscvals.empty())
    Iscvals.resize(numPorts_ );


  double fmin = freqVec_[0];
  double fmax = freqVec_[numFreq_ - 1];

  result.shape(numPorts_, numPorts_);


  std::vector<double> yreal,  yimag;

  yreal.resize(numFreq_ );

  yimag.resize(numFreq_ );

  if (freq < fmin)
  {

    extrapolateData( freq, result, Iscvals, extrapolationLow_ );
  }
  else if (freq > fmax)
  {

    extrapolateData( freq, result, Iscvals , extrapolationHigh_);
  }
  else
  {

    double y0, y1;

    for (int i=0; i< numPorts_ ; i++)
    {

      for (int j=0; j< numPorts_ ; j++)
      {

        for (int n=0; n< numFreq_ ; n++)
        {
          yreal[n] = ((inputNetworkDataVec_[n])(i, j)).real();
          yimag[n] = ((inputNetworkDataVec_[n])(i, j)).imag();
        }

        yInterpolator->clear();
        yInterpolator->init( freqVec_,   yreal );
        yInterpolator->eval( freqVec_,  yreal,  freq, y0 );

        yInterpolator->clear();
        yInterpolator->init( freqVec_,   yimag );
        yInterpolator->eval( freqVec_,  yimag,  freq, y1 );

        result(i,j) = {  y0,  y1  };

      }

      if (IscFD_)
      {

        for (int n=0; n< numFreq_ ; n++)
        {
          yreal[n] = (inputIscFDVec_[n][i]).real();

          yimag[n] = (inputIscFDVec_[n][i]).imag();
        }
        yInterpolator->clear();
        yInterpolator->init( freqVec_,   yreal );
        yInterpolator->eval( freqVec_,  yreal,  freq, y0 );

        yInterpolator->clear();
        yInterpolator->init( freqVec_,   yimag );
        yInterpolator->eval( freqVec_,  yimag,  freq, y1 );

        Iscvals[i] = {  y0,  y1  };
      }
    }

  }

}


//-----------------------------------------------------------------------------
// Function      : Model:: extrapolateData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 8/2023
//----------------------------------------------------------------------------
void Model::extrapolateData( double freq,  Teuchos::SerialDenseMatrix<int, std::complex<double> > & result, std::vector<std::complex<double> >  & Iscvals, int extrapolation )
{

// Perform extrapolation


  if (Iscvals.empty())
    Iscvals.resize(numPorts_ );

  double fmin = freqVec_[0];
  double fmax = freqVec_[numFreq_ - 1];

  std::complex<double> y0, y1;

  int n0, n1, n ;
   
  double f0, f1;

  if (freq < fmin)
  {
    n = 0;

    n0 = 0;    

    n1 = 1;
  }
  else
  {
    n =  numFreq_ - 1;

    n0 =  numFreq_ - 2;

    n1 =  numFreq_ - 1;
  }


  f0 = freqVec_[n0];
  f1 = freqVec_[n1];

  result.shape(numPorts_, numPorts_);


  if ( extrapolation == 0)
  {
     // Use cut off

    result.putScalar ( 0.0 );

    if (IscFD_)
      Iscvals.assign( numPorts_ , 0.0 ); 

  }
  else if ( extrapolation == 1)
  {
    result = inputNetworkDataVec_[n];

    if (IscFD_)
      Iscvals = inputIscFDVec_[n];

  }
  else if ( extrapolation == 2)
  {

    for (int i=0; i< numPorts_ ; i++)
    {

      for (int j=0; j< numPorts_ ; j++)
      {
        y0 = (inputNetworkDataVec_[n0])(i, j);
        y1 = (inputNetworkDataVec_[n1])(i, j);
        result(i, j)= y0 + (y1 - y0)*(freq - f0) / (f1 - f0);

      }

      if (IscFD_)
      {
        y0 = inputIscFDVec_[n0][i];
        y1 = inputIscFDVec_[n1][i];
        Iscvals[i] = y0 + (y1 - y0)*(freq - f0) / (f1 - f0);
      }

    }

  }
  else
  {
    UserFatal(*this) << "Unsupported extrapolation method. ";
    return;
  }

}



//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Load DAE vectors of all YLin instances, regardless of model
///
/// @param solVec solution vector
/// @param fVec f vector
/// @param qVec q vector
/// @param leadF store lead current f vector
/// @param leadQ store lead current q vector
///
/// @return true on success
///
/// @note Because the YLin device re-implements the base-class
/// Master::loadDAEVectors, the Instance::loadDAEFVector method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all YLin instances.
///
/// @see Xyce::Device::YLin::Instance::loadDAEFVector
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec,
                             double * leadF, double * leadQ, double * junctionV, int loadType)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Load DAE matrices for all YLin instances, regardless of model
///
/// @param dFdx matrix of derivatives of F vector with respect to solution
/// @param dQdx matrix of derivatives of Q vector with respect to solution
///
/// @return true on success
///
/// @note Because the YLin device re-implements the base-class
/// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all YLin instances.
///
/// @see Xyce::Device::YLin::Instance::loadDAEdFdx
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEMatrices(Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
{
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Master::loadFreqDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 10/1/2019
//----------------------------------------------------------------------------
bool Master::loadFreqDAEVectors(double frequency, std::complex<double>* solVec,
                                std::vector<Util::FreqVecEntry>& fVec,
                                std::vector<Util::FreqVecEntry>& bVec)
{

  InstanceVector::const_iterator it, end;

  it = getInstanceBegin();
  end = getInstanceEnd();

  Util::FreqVecEntry tmpEntry;
     
  for ( ; it != end; ++it )
  { 
    Instance & inst = *(*it);

    std::vector<std::complex<double> > port_vals;

    for (size_t i = 0; i < inst.extLIDVec.size(); i += 2)
    {
      port_vals.push_back( solVec[inst.extLIDVec[i]] - solVec[inst.extLIDVec[i + 1]] );
    }

    Teuchos::SerialDenseVector<int,std::complex<double> > Fvec( inst.model_.numPorts_ );
    Teuchos::SerialDenseVector<int,std::complex<double> > Xvec( Teuchos::View, &port_vals[0], inst.model_.numPorts_ );

    inst.model_.interpData( frequency,  inst.yvals, inst.Iscvals );

//    inst.model_.interpLin( frequency,  inst.yvals, inst.Iscvals );
                                    

    Fvec.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, Teuchos::ScalarTraits<std::complex<double> >::one(),
                   inst.yvals, Xvec, Teuchos::ScalarTraits<std::complex<double> >::zero() );   


    for (size_t i = 0; i < inst.extLIDVec.size(); i += 2)
    {

      tmpEntry.val = Fvec[i/2];
      tmpEntry.lid = inst.extLIDVec[i];   
      fVec.push_back(tmpEntry);

      tmpEntry.val = inst.Iscvals[i/2];
      tmpEntry.lid = inst.extLIDVec[i];   
      bVec.push_back(tmpEntry);


      // Add RHS vector element for the negative circuit node KCL equ.
      tmpEntry.val = -Fvec[i/2];
      tmpEntry.lid = inst.extLIDVec[i+1];
      fVec.push_back(tmpEntry);


      tmpEntry.val = -inst.Iscvals[i/2];
      tmpEntry.lid = inst.extLIDVec[i+1];   
      bVec.push_back(tmpEntry);

    }


  }
  return true;
}


bool Master::loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                 std::vector<Util::FreqMatEntry>& dFdx)
{
  InstanceVector::const_iterator it, end;


  it =  getInstanceBegin();
  end = getInstanceEnd();

  Util::FreqMatEntry tmpEntry;

  for ( ; it != end; ++it )
  {
    Instance & inst = *(*it);

    for (int i=0; i< inst.model_.numPorts_ ; i++)
    {

      for (int j=0; j< inst.model_.numPorts_ ; j++)
      {
        std::complex<double> tmpVal = std::complex<double>( (inst.yvals)(i, j) );

        tmpEntry.val = tmpVal;
        tmpEntry.row_lid = inst.extLIDVec[2*i + 0];
        tmpEntry.col_lid = inst.LIDVec_[2*i + 0][ 2*j + 0];
        dFdx.push_back(tmpEntry);

        tmpEntry.row_lid = inst.extLIDVec[2*i + 1];
        tmpEntry.col_lid = inst.LIDVec_[2*i + 1][ 2*j + 1];
        dFdx.push_back(tmpEntry);

        tmpEntry.val = -tmpVal;
        tmpEntry.row_lid = inst.extLIDVec[2*i + 0];
        tmpEntry.col_lid = inst.LIDVec_[2*i + 0][2*j + 1];
        dFdx.push_back(tmpEntry);

        tmpEntry.row_lid = inst.extLIDVec[2*i + 1];
        tmpEntry.col_lid = inst.LIDVec_[2*i + 1][2*j + 0];
        dFdx.push_back(tmpEntry);
      }
    }


    // Add RHS vector element for the positive circuit node KCL equ.
/*    tmpEntry.val = std::complex<double>(ri.G, 0.0);
    tmpEntry.row_lid = ri.li_Pos;
    tmpEntry.col_lid = ri.APosEquPosNodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = ri.li_Neg;
    tmpEntry.col_lid = ri.ANegEquNegNodeOffset;
    dFdx.push_back(tmpEntry);

    // Add RHS vector element for the negative circuit node KCL equ.
    tmpEntry.val = std::complex<double>(-ri.G, 0.0);
    tmpEntry.row_lid = ri.li_Pos;
    tmpEntry.col_lid = ri.APosEquNegNodeOffset;
    dFdx.push_back(tmpEntry);

    tmpEntry.row_lid = ri.li_Neg;
    tmpEntry.col_lid = ri.ANegEquPosNodeOffset;
    dFdx.push_back(tmpEntry);      */
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Traits::factory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date :
//-----------------------------------------------------------------------------
///
/// Create a new instance of the YLin device.
///
/// @param configuration
/// @param factory_block
///
Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::registerDevice
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 8/19/2019
//-----------------------------------------------------------------------------
///
/// Define how to use the device in a netlist.
///
/// This method is called from the Xyce::Device::registerOpenDevices
///
/// The device is declared here to be an "YLIN" device, which is required
/// to have a model card of type "LIN".  This device will correspond to model
/// level 1 of LIN models.
void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  static bool initialized = false;

  if (!initialized && (deviceMap.empty() || (deviceMap.find("LIN")!=deviceMap.end())))
  {
    initialized = true;

    Config<Traits>::addConfiguration()
      .registerDevice("lin", 1)
      .registerModelType("lin", 1);
  }
}


} // namespace YLin
} // namespace Device
} // namespace Xyce

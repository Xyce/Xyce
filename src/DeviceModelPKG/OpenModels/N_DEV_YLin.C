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
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Message.h>
#include <N_IO_ParsingHelpers.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_RFparams.h>

namespace Xyce {
namespace Device {
namespace YLin {

///
/// Common Jacobian Stamp for all YLin devices.
/// Because all resistors have identical Jacobian stamps, this data is
/// declared static and is shared by all resistor instances.
///
std::vector<std::vector<int> > Instance::jacStamp;

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
// Purpose       : This is the initial version of a parser to read in a
//                 Touchstone 2 formatted file.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
bool Model::readTouchStoneFile()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Processing Touchstone2 input file " << TSFileName_
                 << " for model " << getName() << std::endl;
  }

  // where we are in the file, and whether parsing had an error
  int lineNum = 0;
  bool firstLineFound = false;
  int numVersionLinesFound = 0;
  int numOptionLinesFound = 0;
  int numPortsLinesFound = 0;
  int numTwoPortDataOrderLinesFound = 0;
  int numFreqLinesFound = 0;
  int numReferenceLinesFound = 0;
  int numMatrixFormatLinesFound = 0;
  int numNetworkDataLinesFound = 0;
  bool skipReadNextLine = false;
  bool psuccess = true;

  // try to open the data file
  std::ifstream inputFile;

  // Error out if the user-specified file does not exist, cannot
  // be opened, or is a directory name rather than a file name.
  // See SON Bug 785 and SRN Bug 2100 for more details.
  if ( !(Util::checkIfValidFile(TSFileName_)) )
  {
    Report::UserError() << "Touchstone2 input file \"" << TSFileName_ << "\" for model "
                        << getName() << " could not be found.";
    return false;
  }

  // open the Touchstone 2 formatted input file
  inputFile.open( TSFileName_.c_str(),std::ifstream::in);
  if( !inputFile.good() )
  {
    Report::UserError() << "Touchstone file \"" << TSFileName_ << "\" for model "
                        << getName() << " could not be opened.";
    return false;
  }

  // parse the file
  std::string aLine;
  IO::TokenVector parsedLine;
  readTouchStoneFileLine(inputFile,aLine,lineNum);

  while( (!inputFile.eof()) || (aLine.substr(0,5) == "[End]") )
  {
    if ( aLine[0] != TSCommentChar_ )
    {
      // [Version] line must be the first non-blank and non-comment line in
      // a Touchstone 2 formatted file.  This parser assumes a version 2.0 file.
      // So, the option line must be the next non-comment line after the [Version]
      // line.
      if (aLine.substr(0,9) == "[Version]")
      {
        if (!firstLineFound)
	{
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
	    Xyce::dout() << "Found [Version] line at lineNum " << lineNum << std::endl;
          }

          firstLineFound = true;
          ++numVersionLinesFound;
          splitTouchStoneFileLine(aLine,parsedLine);

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
        }
        else if (firstLineFound || numVersionLinesFound > 1)
	{
          Report::UserError() << "Invalid [Version] line in file " << TSFileName_
               << " for model " << getName() << " at line " << lineNum;
          return false;
        }

        // skip over any comment lines
        if (!inputFile.eof())
        {
	  readTouchStoneFileLine(inputFile,aLine,lineNum);
        }
        while( (!inputFile.eof()) && ( aLine[0] == TSCommentChar_) )
	{
          readTouchStoneFileLine(inputFile,aLine,lineNum);
        }

        // now parse the Option line, which starts with #
        if (inputFile.eof() || (aLine[0] != '#') )
	{
          Report::UserError() << "Option line not found immediately after [Version] line in file "
            << TSFileName_ << " for model " << getName() << " at line " << lineNum;
          return false;
        }
        else if ( numOptionLinesFound < 1 )
	{
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Found Option line at lineNum " << lineNum << std::endl;
          }

          // Only read the first Option line.  Any others are silently ignored.
          ++numOptionLinesFound;

          // The fields on the Option line can be in any order. Yuck!  Also an Option
          // line that only has a # is allowed.
          splitTouchStoneFileLine(aLine,parsedLine);

          for (int i=1; i<parsedLine.size(); ++i)
	  {
            ExtendedString tokenStr(parsedLine[i].string_);
            if ( (tokenStr == "S") || (tokenStr == "Y") || (tokenStr == "Z") )
	    {
              paramType_ = tokenStr[0];
            }
            else if ( (tokenStr == "RI") || (tokenStr == "MA") || ( tokenStr == "DB") )
	    {
              dataFormat_ = tokenStr;
            }
            else if ( (tokenStr.toUpper() == "HZ") || (tokenStr.toUpper() == "KHZ") ||
                      (tokenStr.toUpper() == "MHZ") || (tokenStr.toUpper() == "GHZ") )
	    {
              ExtendedString freqMultStr("1.0"+parsedLine[i].string_);
              freqMultiplier_ = freqMultStr.Value();
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
        }
      }
      else if (aLine.substr(0,17) == "[Number of Ports]")
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

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Set numPorts_ = " << numPorts_ << " at lineNum " << lineNum << std::endl;
          }
        }
      }
      else if (aLine.substr(0,21) == "[Two-Port Data Order]")
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
      else if (aLine.substr(0,23) ==  "[Number of Frequencies]")
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
      else if (aLine.substr(0,11) ==  "[Reference]")
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
          // Handle line 1 of [Reference} line block
          for (int i=1; i<parsedLine.size(); ++i)
	  {
            ExtendedString z0Str(parsedLine[i].string_);
            Z0Vec_.push_back(z0Str.Value());
          }

          // read next line to see if the [Reference] block spans multiple lines
          if (!inputFile.eof())
          {
	     readTouchStoneFileLine(inputFile,aLine,lineNum);
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
            readTouchStoneFileLine(inputFile,aLine,lineNum);
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
      else if (aLine.substr(0,15) ==  "[Matrix Format]")
      {
        // This line is optional.  If included, it must appear after the [Number of Ports] line
        // and before the [Network Data] line.  It has the allowed string values of
        // "Full", "Lower" or "Upper".  This parsing only supports "Full" now.
        ++numMatrixFormatLinesFound;
        splitTouchStoneFileLine(aLine,parsedLine);

        if ( (parsedLine.size() != 3) || (numPorts_ < 1) || (numNetworkDataLinesFound > 0) )
        {
          Report::UserError() << "Invalid [Matrix Format] line(s) in file " << TSFileName_
             << " for model " << getName() << " at line " << lineNum;
          psuccess = false;
        }
        else
	{
	  std::string matrixFormat(parsedLine[2].string_);
          if ( matrixFormat != "Full" )
          {
            Report::UserError() << "Only [Matrix Format] = Full is supported";
            psuccess = false;
          }

          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
            Xyce::dout() << "Found matrix format = " << matrixFormat << " at lineNum " << lineNum << std::endl;
          }
        }

      }
      else if (aLine.substr(0,23) ==  "[Network Data]")
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
	  readTouchStoneFileLine(inputFile,aLine,lineNum);
        }
        while( (!inputFile.eof()) && ( aLine[0] == TSCommentChar_) )
	{
          readTouchStoneFileLine(inputFile,aLine,lineNum);
        }

        Teuchos::SerialDenseMatrix<int, std::complex<double> > inputNetworkData;
        inputNetworkData.shape(numPorts_, numPorts_);

	// now populate the matrix assuming the input is in "Full" format
	while ( (!inputFile.eof()) && ( aLine[0] !='[') )
	{
          ++numDataLinesFound;
          splitTouchStoneFileLine(aLine,parsedLine);
          if ( parsedLine.size() != (2*(numPorts_*numPorts_)+1 ) )
	  {
            Report::UserError() << "Incorrect number of entries for network data on lineNum "
               << lineNum << " in file " << TSFileName_ << " for model " << getName();
            return false;
          }
          else
	  {
            ExtendedString freqStr(parsedLine[0].string_);
            freqVec_.push_back(freqMultiplier_*freqStr.Value());

            // for 2-port networks, this assumes [Two-Port Data Order] = "12_21"
	    for (int i=0; i<numPorts_; ++i)
	    {
              for (int j=0; j<numPorts_; ++j)
	      {
                // data will be converted to RI format for internal use in
                // YLin model
                ExtendedString Str1(parsedLine[2*(i*numPorts_+j)+1].string_);
	        ExtendedString Str2(parsedLine[2*(i*numPorts_+j)+2].string_);
                if (dataFormat_ == "RI")
		{
                  inputNetworkData[i][j].real(Str1.Value());
                  inputNetworkData[i][j].imag(Str2.Value());
                }
                else if (dataFormat_ == "MA")
		{
                  double mag = Str1.Value();
                  double angle = M_PI*Str2.Value()/180.0;
                  inputNetworkData[i][j].real(mag*cos(angle));
                  inputNetworkData[i][j].imag(mag*sin(angle));
                }
                else if (dataFormat_ == "DB")
		{
                  double mag = pow(10.0,0.05*Str1.Value());
                  double angle = M_PI*Str2.Value()/180.0;
                  inputNetworkData[i][j].real(mag*cos(angle));
                  inputNetworkData[i][j].imag(mag*sin(angle));
                }
                else
		{
                  Report::UserError() << "Incorrect data format " << dataFormat_ << " for network data in file "
                    << TSFileName_ << " for model " << getName();
                  return false;
                }
              }
            }

            if ( (twoPortDataOrder_ == "21_12") && (numPorts_ == 2) )
	    {
	      complex tempVal = inputNetworkData[1][2];
              inputNetworkData[1][2] = inputNetworkData[2][1];
              inputNetworkData[2][1] = tempVal;
	    }

            // YLin model will use Y-parameter format internally
            if (paramType_=='S')
	    {
              Teuchos::SerialDenseMatrix<int, std::complex<double> > YParams;
	      Util::stoy(inputNetworkData,YParams,Z0Vec_);
              inputNetworkDataVec_.push_back(YParams);
            }
            else if (paramType_=='Z')
	    {
              Teuchos::SerialDenseMatrix<int, std::complex<double> > YParams;
	      Util::ztoy(inputNetworkData,YParams);
              inputNetworkDataVec_.push_back(YParams);
            }
            else
	    {
              // input was in Y-parameter format
              inputNetworkDataVec_.push_back(inputNetworkData);
            }
          }

          // read in next line
          readTouchStoneFileLine(inputFile,aLine,lineNum);
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
      readTouchStoneFileLine(inputFile,aLine,lineNum);
    else
      skipReadNextLine=false;
  }

  if ( (Z0Vec_.size() == 0) && (numOptionLinesFound > 0) )
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
  if ( (numVersionLinesFound == 0) || (numOptionLinesFound == 0) )
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

  inputFile.close();

  return psuccess;
}


//-----------------------------------------------------------------------------
// Function      : Model::splitTouchStoneFileLine
// Purpose       : Handle in-line comments, and then split a line from a
//                 Touchstone 2 formatted file into a TokenVector
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
void Model::splitTouchStoneFileLine(const std::string& aLine, IO::TokenVector & parsedLine)
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
// Function      : Model::readTouchStoneFileLine
// Purpose       : Read in a line from a Touchstone 2 formatted file, and
//                 also increment the lineNum counter
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL, Electrical Models & Simulation
// Creation Date : 08/20/2019
//-----------------------------------------------------------------------------
void Model::readTouchStoneFileLine(std::istream & inputFile, std::string& aLine, int& lineNum)
{
  IO::readLine(inputFile,aLine);
  lineNum++;

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
    paramType_('S'),
    dataFormat_("MA"),
    numPorts_(0),
    twoPortDataOrder_(""),
    numFreq_(0),
    Z0Vec_(0),
    freqVec_(0)
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


  // it the file was successfully converted it is in Y-parameters and RI format
  if (TouchstoneFileRead)
  {
    paramType_='Y';
    dataFormat_="RI";
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
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
//  else if (loadType == LINEAR_FREQ)
//  {
//    it = linearInstances_.begin();
//    end = linearInstances_.begin();
//  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ri = *(*it);

    // Load RHS vector element for the positive circuit node KCL equ.
    ri.i0 = (solVec[ri.li_Pos]-solVec[ri.li_Neg])*ri.G;

    fVec[ri.li_Pos] += ri.i0;
    fVec[ri.li_Neg] += -ri.i0;
    if( ri.loadLeadCurrent )
    {
      leadF[ri.li_branch_data] = ri.i0;
      junctionV[ri.li_branch_data] = solVec[ri.li_Pos] - solVec[ri.li_Neg];
    }
  }
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
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
//  else if (loadType == LINEAR_FREQ)
//  {
//    // Don't do anything, this is a frequency domain load
//    it = linearInstances_.begin();
//    end = linearInstances_.begin();
//  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & ri = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(ri.f_PosEquPosNodePtr) += ri.G;
    *(ri.f_PosEquNegNodePtr) -= ri.G;
    *(ri.f_NegEquPosNodePtr) -= ri.G;
    *(ri.f_NegEquNegNodePtr) += ri.G;
#else
    dFdx[ri.li_Pos][ri.APosEquPosNodeOffset] += ri.G;
    dFdx[ri.li_Pos][ri.APosEquNegNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquPosNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquNegNodeOffset] += ri.G;
#endif
  }

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

/*  InstanceVector::const_iterator it, end;

  it = linearInstances_.begin();
  end = linearInstances_.end();

  Util::FreqVecEntry tmpEntry;
     
  for ( ; it != end; ++it )
  { 
    Instance & ri = *(*it);

//    int num_ports = ri.ports_.size();

    std::vector<std::complex<double> > port_vals;

    for (size_t i = 0; i < ri.extLIDVec.size(); i += 2)
    {
      port_vals.push_back( solVec[ri.extLIDVec[i]] - solVec[ri.extLIDVec[i + 1]] );
    }

    Teuchos::SerialDenseVector<int,std::complex<double> > Fvec( ri.numPorts_ );
    Teuchos::SerialDenseVector<int,std::complex<double> > Xvec( Teuchos::View, &port_vals[0], ri.numPorts_ );



//    interpLin( frequency, freqVec_,  inputNetworkDataVec_);
                                    
    Fvec.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, Teuchos::ScalarTraits<std::complex<double> >::one(),
                     ri.Ymatrix_, Xvec, Teuchos::ScalarTraits<std::complex<double> >::zero() );


    for (size_t i = 0; i < ri.extLIDVec.size(); i += 2)
    {
      tmpEntry.val = Fvec[i/2];
      tmpEntry.lid = ri.extLIDVec[i];   
      fVec.push_back(tmpEntry);

      // Add RHS vector element for the negative circuit node KCL equ.
      tmpEntry.val = -Fvec[i/2];
      tmpEntry.lid = ri.extLIDVec[i+1];
      fVec.push_back(tmpEntry);
    }


  }   */
  return true;
}


bool Master::loadFreqDAEMatrices(double frequency, std::complex<double>* solVec,
                                 std::vector<Util::FreqMatEntry>& dFdx)
{
  InstanceVector::const_iterator it, end;
/*
  it = linearInstances_.begin();
  end = linearInstances_.end();

  Util::FreqMatEntry tmpEntry;

  for ( ; it != end; ++it )
  {
    Instance & ri = *(*it);

    // Add RHS vector element for the positive circuit node KCL equ.
    tmpEntry.val = std::complex<double>(ri.G, 0.0);
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
    dFdx.push_back(tmpEntry);
  }
*/
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::YLin::Master::storeInstance
// Purpose       :
// Special Notes : The load lead current logic can be removed when lead currents
//                 are moved out of device loading.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 8/1/17
//-----------------------------------------------------------------------------
void Master::storeInstance( const FactoryBlock& factory_block, Instance* instance )
{
  Xyce::Device::DeviceMaster<Traits>::storeInstance( factory_block, instance );

  bool loadLeadCurrent = false;

  const std::set<std::string>& leadCurrentSet = factory_block.deviceManager_.getDevicesNeedingLeadCurrentLoads();
  std::string outputName = (instance->getName()).getEncodedName();
  if ( factory_block.deviceOptions_.calculateAllLeadCurrents ||
      leadCurrentSet.find(outputName) != leadCurrentSet.end() )
  {
    loadLeadCurrent = true;
  }

  if ( instance->isLinearDevice() && !loadLeadCurrent )
  {
    linearInstances_.push_back( instance );
  }
  else
  {
    nonlinearInstances_.push_back( instance );
  }
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

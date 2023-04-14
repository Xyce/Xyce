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

//-----------------------------------------------------------------------------
//
// Purpose        : Define the N_IO_ParameterBlock class instantiations of
//                  which are associated with netlist .MODEL lines.
//
// Special Notes  : ERK.  It seems that the name "N_IO_ModelBlock" would have been
//                  more appropriate and less confusing, or possibly
//                  N_IO_ModelParametersBlock.  Calling it "parameters block"
//                  makes it sound much more general than it really is.
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/10/2001
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <algorithm>
#include <cstring>

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_IO_CircuitMetadata.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_Report.h>
#include <N_IO_fwd.h>
#include <N_DEV_Param.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_PDS_Comm.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::ParameterBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
ParameterBlock::ParameterBlock(
    std::string const& fileName,
    TokenVector const& parsedInputLine)
: defaultApplied_(false)
{
  setLevel("1");
  modelData.setNetlistLocation(NetlistLocation(fileName, parsedInputLine[0].lineNumber_));

  extractModelData( parsedInputLine );
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::print
// Purpose       : Output the details of a device block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
void ParameterBlock::print()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << "Parameter Block Information" << std::endl;
  Xyce::dout() << "---------------------------" << std::endl;
  Xyce::dout() << "  name : " << getName() << std::endl;
  Xyce::dout() << "  type : " << getType() << std::endl;
  Xyce::dout() << "  level: " << getLevel() << std::endl;

  Xyce::dout() << "  parameters: " << std::endl;
  int numParameters = getNumberOfParameters();
  for (int i = 0; i < numParameters; ++i)
  {
    Xyce::dout() << "  " << getParameter(i).uTag() << " : ";
    Xyce::dout() << getParameter(i).stringValue();
    if ( getParameter(i).isTimeDependent() )
    {
      Xyce::dout() << "  time dependent";
    }
    Xyce::dout() << std::endl;
  }
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::extractModelData
// Purpose       : Extract model data from parsedInputLine using model metadata.
//
// Special Notes : ERK - this function doesn't appear to actually use the
//                 metadata object.
//
//                 Unfortunately, this model data gets extracted before
//                 the models defaults get determined.  That makes it much
//                 harder to deal with unusual cases like VECTOR and
//                 VECTOR-COMPOSITE.  Defaults are set up later.
//
//                 Instances are handled in the opposite order.
//
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
bool ParameterBlock::extractModelData(TokenVector const& parsedInputLine)
{
  const int numFields = parsedInputLine.size();

  // Check that the minimum required number of fields are on the line.
  if ( numFields < 3 )
  {
    Report::UserError().at(modelData.getNetlistLocation())
      << ".model line has too few fields";
    return false;
  }

  // Extract the model name and type from parsedInputLine.
  ExtendedString field ( parsedInputLine[1].string_ );
  field.toUpper();
  setName( field );
  field = parsedInputLine[2].string_;
  field.toUpper();
  setType( field );

  // Set the start and end position of the model parameters in parsedInputLine
  // accounting for optional enclosing parentheses.
  int paramStartPos;
  int offsetEq;
  bool levelSet = false;
  std::vector<Device::Param> inputParameters;
  if ( numFields > 3 )
  {
    int paramEndPos;
    bool beginParenFound = false;
    bool endParenFound = false;

    if ( parsedInputLine[3].string_ == "(" )
    {
      paramStartPos = 4;
      beginParenFound = true;
    }
    else
    {
      paramStartPos = 3;
    }
    
    // can't just look at the end of the list for a matching paren for the first one as
    // last item could be matched up with one from last "tag = value" as in tag = (xxxx)
    if((parsedInputLine[numFields - 1].string_ == ")") && (parsedInputLine[numFields - 3].string_ != "("))
    {
      endParenFound=true;
    }

    //if ( (beginParenFound && parsedInputLine[numFields - 1].string_ != ")") ||
    //     (!beginParenFound && parsedInputLine[numFields - 1].string_ == ")") )
    if( beginParenFound != endParenFound )
    {
      Report::UserError().at(modelData.getNetlistLocation())
        << "Unmatched parentheses for .model " << getName();
      return false;
    }
    else if ( parsedInputLine[numFields - 1].string_ == ")" )
    {
      paramEndPos = numFields - 2;
    }
    else
    {
      paramEndPos = numFields - 1;
    }

    // flag equal sign usage in parameter list
    if ( ( paramStartPos + 1 >= numFields ) ||
         ( parsedInputLine[paramStartPos + 1].string_ == "=" ) )
    {
      offsetEq = 0;
    }
    else
    {
      offsetEq = -1;
    }

    // Extract the parameters from parsedInputLine.
    Device::Param parameter("", "");

    // new way (ERK), which handles traditional name=value format params,
    // VECTOR params, and VECTOR-COMPOSITE params.
    for (int i = paramStartPos; i <= paramEndPos-2-offsetEq; i+=3+offsetEq )
    {
      bool paramDone(false);
      if ( parsedInputLine[i+1].string_ == "=" && parsedInputLine[i+2].string_ == "{" )
      {
        // This might be a VECTOR or a VECTOR-COMPOSITE
        //   1. Find the matching final brace.
        //   2. check for "=" signs.  If present between curly braces,
        //   then this is a vector-composite.
        int unmatchedLeftBraces = 1;
        int tmpIndex=i+3;
        int numEqualSigns=0;
        int numCommas=0;
        int finalLeftBraceIndex=i+3;

        for (; tmpIndex <= paramEndPos;++tmpIndex)
        {
          if ( parsedInputLine[tmpIndex].string_  == "}")
          {
            unmatchedLeftBraces--;
            if (unmatchedLeftBraces==0)
            {
              finalLeftBraceIndex=tmpIndex;
              break;
            }
          }

          if ( parsedInputLine[tmpIndex].string_  == "{" )
          {
            ++unmatchedLeftBraces;
          }

          if ( parsedInputLine[tmpIndex].string_ == "=" )
          {
            ++numEqualSigns;
          }

          if ( parsedInputLine[tmpIndex].string_ == "," )
          {
            ++numCommas;
          }
        }

        if (unmatchedLeftBraces > 0)
        {
          Report::UserError().at(modelData.getNetlistLocation())
            << "Unmatched curly braces in .model " << getName();
          return false;
        }

        bool composite(numEqualSigns > 0);

        if (composite)
        {
          if (DEBUG_IO)
          {
            Xyce::dout() 
              << " ParameterBlock::extractModelData: YES! I found a VECTOR-COMPOSITE!!!!" << std::endl;
          }

          int linePosition=i+3;
          int savedLinePosition=linePosition;
          int numBlocks = 0;
          int numComponents = 0;
          ExtendedString paramBase ( parsedInputLine[i].string_ );
          paramBase.toUpper();

          // count the number of blocks, to get strawman structure.
          // Note: The private function validLinePosition_() is used in several places below,
          // to stop segfaults caused by trying to read past the end of the parsedInputLine vector.
          std::vector<Device::Param> tmpParamVec;
          while (validLinePosition_(linePosition,numFields,paramBase) && (parsedInputLine[linePosition].string_ != "}") )
          {
            ExtendedString component ( parsedInputLine[linePosition].string_ );
            component.toUpper();
            ++numComponents;

            parameter.set(component, 0.0);
            tmpParamVec.push_back(parameter);
            if (DEBUG_IO)
              Xyce::dout() << "parameter = " << parameter << std::endl;

            linePosition += 2;
            
            int blockIndex=0;
            bool startVC = true;
            while ( validLinePosition_(linePosition,numFields,paramBase) &&
                   (startVC || parsedInputLine[linePosition].string_ == ",") )
            {
              if (!startVC) ++linePosition;
              startVC = false;

              ++linePosition;
              ++blockIndex;
            } // end of while loop.

            numBlocks = blockIndex;
          } // end of while loop.
          if (DEBUG_IO)
            Xyce::dout() << "numBlocks = " << numBlocks << std::endl;

          inputCompositeParamVecMap[paramBase].resize(numBlocks);
          int iblock=0;
          for (iblock=0;iblock<numBlocks;++iblock)
          {
            inputCompositeParamVecMap[paramBase][iblock] = tmpParamVec;
            if (DEBUG_IO)
            {
              Xyce::dout() << "paramBase = " << paramBase;
              Xyce::dout() << "  iblock = " << iblock << std::endl;
              for (size_t ieric=0;ieric<tmpParamVec.size();++ieric)
              {
                Xyce::dout() << "par["<<ieric<<"]="<<tmpParamVec[ieric]<<std::endl;
              }
            }
          }

          linePosition = savedLinePosition;

          // fill in the allocated data structure
          size_t componentIndex=0;
          while (validLinePosition_(linePosition,numFields,paramBase) && (parsedInputLine[linePosition].string_ != "}") )
          {
            ExtendedString component ( parsedInputLine[linePosition].string_ );
            component.toUpper();

            // In the instance version of vector-composite processing, there is
            // an error check at this point to see if the specified component
            // exists or not.  We can't do that here, because this function
            // is called prior to metadata being set up.

            linePosition += 2;

            // Mark the component as given in components. Later we will
            // add all of the components and their defaults that were
            // not given.
            //paramIter->setGiven(true);

            size_t blockIndex=0;
            bool startVC = true;
            while ( validLinePosition_(linePosition,numFields,paramBase) &&
                    (startVC || parsedInputLine[linePosition].string_ == ",") )
            {
              if (validLinePosition_(linePosition+1,numFields,paramBase) && !startVC) ++linePosition;
              startVC = false;

              std::string value = parsedInputLine[linePosition].string_;
              if (DEBUG_IO)
                Xyce::dout() << " ParameterBlock::extractModelData: value = " << value << std::endl;

              //
              // This line commented  out by TVR on 4 Feb 2009.  It is not
              // clear that it is necessary in any case to upper-case the
              // values of parameters as a matter of course, and it definitely
              // breaks things when the parameter is a file name and the
              // platform has a file system that honors case (e.g. UNIX as
              // opposed to Losedows or Mac).
              // It is probably better to let the code that *uses* a parameter
              // decide whether to force it to upper.
              //
              // value.toUpper();

              // Check if the insertion into inputCompositeParamVecMap will 
              // cause a segfault
              if ( (blockIndex > inputCompositeParamVecMap[paramBase].size() -1) ||
                   (componentIndex > inputCompositeParamVecMap[paramBase][blockIndex].size() -1) )
              {
                Report::UserError().at(modelData.getNetlistLocation())
                  << "Fatal error when parsing vector-composite parameter " << paramBase 
                  << " in .model statement " << getName();
                return false;
              }
                   
              inputCompositeParamVecMap[paramBase][blockIndex][componentIndex].setVal(value);
              inputCompositeParamVecMap[paramBase][blockIndex][componentIndex].setGiven( true );

              ++linePosition;
              ++blockIndex;
            } // end of while loop.

            ++componentIndex;
            //numComponents = blockIndex;
          } // end of while loop.

          paramDone=true;
          i=finalLeftBraceIndex-2;

          if (DEBUG_IO)
          {
            Xyce::dout() << "User-input .MODEL statement vector-composite for param = "<< paramBase << std::endl;
            std::vector<std::vector<Device::Param> > & tmpParamVecVec = inputCompositeParamVecMap[paramBase];
            for (size_t ivec=0;ivec<tmpParamVecVec.size();++ivec)
            {

              Xyce::dout() << "Composite for column (block) "<< ivec<< ":" <<std::endl;
              std::vector<Device::Param> & tmpVec = tmpParamVecVec[ivec];
              int vecSize = tmpVec.size();

              for (int ip=0; ip<vecSize; ++ip)
              {
                Xyce::dout() << "param["<<ip<<"] = " << tmpVec[ip];
              }
              Xyce::dout() << "-----------" << std::endl;
            }
          }

          continue;
        }
        else
        {
          if (DEBUG_IO)
          {
            Xyce::dout() << " ParameterBlock::extractModelData: YES! I found a VECTOR!!!!" << std::endl;
            for (int k = i; k<=finalLeftBraceIndex; k++)
              Xyce::dout() << "\"" << parsedInputLine[k].string_ << "\" ";
            Xyce::dout() << std::endl;
          }
          
          // Extract parameter name and value from parsedInputLine and add to
          // parameter list, check for LEVEL parameter treat as special case.
          field = parsedInputLine[i].string_;
          field.toUpper();
          std::vector<std::string> values;
          for( int k=i+3; k<finalLeftBraceIndex; k+=2 )
          {
            //Xyce::dout() << "adding \"" << parsedInputLine[k].string_ << "\"" << std::endl;
            values.push_back( parsedInputLine[k].string_ );
          }
          parameter.set( field, values );
          parameter.setGiven( true );
          inputParameters.push_back( parameter );
          paramDone=true;
          i=finalLeftBraceIndex-2;
          continue;
        }
      }

      // If paramDone=false, then it was not a VECTOR or a VECTOR-COMPOSITE.
      if ( ( parsedInputLine[i+1].string_ == "=" && !paramDone ) ||
           ( offsetEq == -1 && !paramDone ) )
      {
        // Extract parameter name and value from parsedInputLine and add to
        // parameter list, check for LEVEL parameter treat as special case.
        field = parsedInputLine[i].string_;
        field.toUpper();
        
        // carefully check for the right structure of typeQualifier ( string_value )
        //if( numFields > (i+5+offsetEq+1) )
        ExtendedString typeQualifierValue=parsedInputLine[i+2+offsetEq].string_;
        ExtendedString openingParen("");
        ExtendedString closingParen("");
        // use caution looking forward on list for opening and closing parens 
        if( (i+5+offsetEq)<numFields)
        {
          openingParen=parsedInputLine[i+3+offsetEq].string_;
          closingParen=parsedInputLine[i+5+offsetEq].string_;
        }
        typeQualifierValue.toUpper();
        if ((openingParen=="(") && (closingParen==")") &&
           ((typeQualifierValue == "TABLEFILE") || (typeQualifierValue == "STRING")))
        {
          std::string parameterValue=parsedInputLine[i+2+offsetEq].string_ + parsedInputLine[i+2+offsetEq+2].string_;
          parameter.set( field, parameterValue );
          // have to advance the counter i due to the extra parens in this syntax
          i+=3;
        }
        else if ((openingParen=="(") && (closingParen==")") )
        {
          Report::UserError().at(modelData.getNetlistLocation())
            << "Unknown string type qualifier.  Should be either TABLEFILE or STRING in model " << getName();
        }
        else
        {
          parameter.set( field, parsedInputLine[i+2+offsetEq].string_ );
        }
        if ( parameter.uTag() != "LEVEL" )
        {
          parameter.setGiven( true );
          inputParameters.push_back( parameter );
        }
        else
        {
          setLevel( parameter.stringValue() );
          levelSet = true;
        }
        paramDone=true;
      }

      if (!paramDone)
      {
        Report::UserError().at(modelData.getNetlistLocation())
          << "Parameter " << parsedInputLine[i].string_  << " not formatted correctly in .model " << getName();
        return false;
      }
    }
  }

  if (!levelSet)
  {
    setLevel("1");
  }

  if (DEBUG_IO)
  {
    Xyce::dout() << " ParameterBlock::extractModelData.  inputParameters: " << std::endl;

    int iIPSize = inputParameters.size();
    for (int ieric=0;ieric<iIPSize;++ieric)
    {
      Xyce::dout() << inputParameters[ieric];
    }
  }

  // For now, accept any parameter.  We will have to wait until the model metadata is
  // generated before checking that the parameters are valid.  This will not happen
  // until pass two.

  // Insert the vector of input parameters into the official parameter vector
  // for this model.
  addParameters(inputParameters);

  // ERK: 12/14/09 note.  at this stage the input parameters were conditionally added
  // to the expressionValuedParms object.  However, in order to include parameters
  // from vector-composites, this process is now done at the end of the
  // addDefaultModelParameters function (at which point vector-composites have been
  // added in to modelData.params object).

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::validLinePosition_
// Purpose       : Generate fatal error if parameter block parsing is
//               : about to walk off the end of the parsedLine vector.
//               : Otherwise return true.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 10/1/2017
//-----------------------------------------------------------------------------
bool ParameterBlock::validLinePosition_(int pos, int lineLength, ExtendedString paramName )
{
  if (pos > lineLength-1)
  {
    Report::UserFatal() << "Fatal error parsing vector-composite parameter " << paramName 
      << " in .model statement " << getName();
    return false;
  }
  
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addDefaultModelParameters
// Purpose       : Add the default parameters for a model from metadata.
// Special Notes : While not required, it is generally expected that
//                 the parameters given in the input netlist will have
//                 already been extracted via extractModelData().
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/17/2001
// Last Modified : 05/28/2002
//-----------------------------------------------------------------------------
void ParameterBlock::addDefaultModelParameters(CircuitMetadata & metadata )
{
  std::map<std::string,bool> paramMetadataExistMap;

  if (defaultApplied_)
    return;

  // get model metadata
  const std::vector<Device::Param> &modelParameterPtr = metadata.getModelParameters(getType(), getLevel());

  // Process the metadata.
  if (DEBUG_IO)
  {
    Xyce::dout() << "BEFORE adding final params:"<<std::endl;
    for (const auto &modelDataParam : modelData.params)
    {
      Xyce::dout() << modelDataParam;
    }
    Xyce::dout() << "End of BEFORE parameter List"<<std::endl;
  }

  for (const auto &param : modelParameterPtr)
  {
    std::string paramSVal(param.stringValue());

    paramMetadataExistMap[param.tag()] = true;

    // ERK:  If this is a vector-composite param, we have to do some special stuff.
    // The metadata is stored in a different structure, and the user-specified
    // stuff is in a different structure as well.  Ultimately, however, the
    // vector-composite components has to be converted to the param structure.
    if (paramSVal == "VECTOR-COMPOSITE")
    {
      modelData.params.push_back(param);
      addDefaultCompositeModelParameters_(metadata, param, paramMetadataExistMap);
    }

    if (paramSVal == "VECTOR")
    {
      Xyce::dout() << "ParameterBlock::addDefaultModelParameters() Parameter type is VECTOR" << std::endl;
    }
  }

  if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
  {
    Xyce::dout() << " ParameterBlock::addDefaultModelParameters: metadata modelParameterPtr[] (list of defaults):"<<std::endl;
    for (const auto &param : modelParameterPtr)
    {
      Xyce::dout() << param;
    }
  }

  // loop over user-specified params.  At this point the vector composites (if they
  // exist), including defaults, have already been added the params vector.
  int pos = 0;
  unordered_map<std::string,std::vector<Device::Param>::const_iterator> paramDefaultProcessed;
  for (std::vector<Device::Param>::iterator paramIter = modelData.params.begin(), paramIterEnd = modelData.params.end(); paramIter != paramIterEnd; )
  {
    const Device::Param & param = (*paramIter);

    if (param.tag() == "INDEPENDENT;PARAM")
    {
      std::vector<Device::Param> defaults;
      for (std::vector<Device::Param>::const_iterator pparamIter = modelParameterPtr.begin(), pparamIterEnd = modelParameterPtr.end() ; pparamIter != pparamIterEnd ; ++pparamIter)
      {
        if (paramDefaultProcessed.find((*pparamIter).tag()) == paramDefaultProcessed.end())
        {
          defaults.push_back(*pparamIter);
        }
      }
      paramDefaultProcessed.clear();

      modelData.params.insert(paramIter, defaults.begin(), defaults.end());
      pos += defaults.size()+1;
      paramIter = modelData.params.begin() + pos;
      paramIterEnd = modelData.params.end();
    }
    else if (paramMetadataExistMap.find(param.tag()) == paramMetadataExistMap.end())
    {
      Report::UserWarning().at(modelData.getNetlistLocation())
        << "No model parameter " << param.tag() << " found for model " << getName() << " of type " << getType() << ", parameter ignored.";

      paramIter = modelData.params.erase(paramIter);
      paramIterEnd = modelData.params.end();
    }
    else
    {
      paramDefaultProcessed[(*paramIter).tag()] = paramIter;
      ++paramIter;
      ++pos;
    }
  }

  // Any parameters, that have not been pushed back on the params vector, but
  // which do exist in metadata, should be pushed back at this point.  Vector-composite
  // defaults are already here, but other defaults may not be.
  // for (paramIter=modelParameterPtr.begin() ; paramIter != modelParameterPtr.end() ; ++paramIter)
  for (std::vector<Device::Param>::const_iterator paramIter = modelParameterPtr.begin(), paramIterEnd = modelParameterPtr.end(); paramIter != paramIterEnd; ++paramIter)
  {
    if (paramDefaultProcessed.find((*paramIter).tag()) == paramDefaultProcessed.end())
    {
      modelData.params.push_back(*paramIter);
    }
  }

  // Add to the vector of expression parameter pointers.
  // ERK. 12/14/09 note:  this procedure used to be done
  // at the end of the extractModelData function.  I moved
  // it here to include vector-composite params.
  int numParameters = modelData.params.size();
  for ( int i = 0; i < numParameters; ++i )
  {
    Device::Param & parameter = modelData.params[i];
    if (parameter.hasExpressionValue() || parameter.isQuoted() || parameter.isTableFileTypeQuoted() || parameter.isStringTypeQuoted() )
    {
      expressionValuedParams_.push_back(parameter);
    }
  }

  if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
  {
    Xyce::dout() << "AFTER adding final default params:"<<std::endl;
    for (const auto &modelDataParam : modelData.params)
    {
      Xyce::dout() << modelDataParam;
    }
    Xyce::dout() << "End of AFTER parameter List"<<std::endl;
  }

  defaultApplied_ = true;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addDefaultCompositeModelParameters_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 05/15/2008
//-----------------------------------------------------------------------------
void
ParameterBlock::addDefaultCompositeModelParameters_(
  CircuitMetadata &             metadata,
  const Device::Param &         baseParam ,
  std::map<std::string,bool> &  paramMetadataExistMap)
{
  size_t icomp(0), ic(0), iNumCols(0);
  std::string typ(getType());
  int lev(getLevel());

  // Get the defaultComponents from metadata.
  std::vector<Device::Param> defaultComponents;
  std::string baseTag(baseParam.uTag());
  //metadata.getModelCompositeComponents(string(""),
  metadata.getModelCompositeComponents(getType(), baseTag, lev, defaultComponents);

  size_t iDefaultCompSize=defaultComponents.size();

  if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
  {
    Xyce::dout() << "FOUND A COMPOSITE! "
                 << "ParameterBlock::addDefaultCompositeModelParameters_: " << baseParam;

    for (icomp=0;icomp<iDefaultCompSize;++icomp)
    {
      Xyce::dout() << "component["<<icomp<<"] : ";
      Xyce::dout() << defaultComponents[icomp];
    }
  }

  if (inputCompositeParamVecMap.find(baseTag) == inputCompositeParamVecMap.end())
  {
    // do nothing.  If user did not specify, then don't include it in params.
  }
  else
  {

    // push defaults into the inputCompositeParamVecMap, as needed.
    for (icomp=0;icomp<iDefaultCompSize;++icomp)
    {
      Device::Param & paramDefault = defaultComponents[icomp];
      std::string utag1(defaultComponents[icomp].uTag());

      bool found=false;
      size_t iInputCompSize =  inputCompositeParamVecMap[baseTag][0].size();
      for (size_t itmp=0;itmp<iInputCompSize;++itmp)
      {
        std::string utag2(inputCompositeParamVecMap[baseTag][0][itmp].uTag());
        if (utag1 == utag2)
        {
          found=true;
          break;
        }
      }

      // if not found, then user did not specify.
      // Vector-composite is always specified as a matrix.  So, if the parameter
      // wasn't specified for one column, it wasn't specified for any of them.
      if (!found)
      {
        iNumCols = inputCompositeParamVecMap[baseTag].size();
        for (ic=0;ic<iNumCols;++ic)
        {
          inputCompositeParamVecMap[baseTag][ic].push_back(paramDefault);
        }
      }
    }

    size_t iInputCompSize =  inputCompositeParamVecMap[baseTag][0].size();
    if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
    {
      Xyce::dout() << "ParameterBlock::addDefaultCompositeModelParameters_: input composite for " << baseTag << ": " <<std::endl;
      iNumCols = inputCompositeParamVecMap[baseTag].size();
      Xyce::dout() << "Column: " ;
      for (ic=0;ic<iNumCols;++ic)
      {
        Xyce::dout() << "\t\t" << ic;
      }
      Xyce::dout() << std::endl;
      for (icomp=0;icomp<iInputCompSize;++icomp)
      {
        Xyce::dout().width(14);
        Xyce::dout() << inputCompositeParamVecMap[baseTag][0][icomp].tag();
        for (ic=0;ic<iNumCols;++ic)
        {
          if (icomp >= inputCompositeParamVecMap[baseTag][ic].size())
          {
            std::string msg("ParameterBlock::addDefaultCompositeModelParameters_: inputCompositeParamVecMap overrun!");
            Report::UserFatal() << msg;
          }

          Device::Param & tmpPar = inputCompositeParamVecMap[baseTag][ic][icomp];
          Xyce::dout() << "\t";
          Xyce::dout().width(20);
          Xyce::dout() << tmpPar.stringValue();
          Xyce::dout() << "(given=";
          if (tmpPar.given())
            Xyce::dout() << "TRUE) ";
          else
            Xyce::dout() << "FALSE)";
        }
        Xyce::dout() << std::endl;
      }
      Xyce::dout() << std::endl;
    }

    // convert from composite components to conventional parameters.
    std::ostringstream paramName;
    ExtendedString paramBase ( baseTag );
    paramBase.toUpper();
    std::vector<Device::Param> convertedCompositeParams;

    iInputCompSize =  inputCompositeParamVecMap[baseTag][0].size();
    iNumCols = inputCompositeParamVecMap[baseTag].size();
    for (icomp=0;icomp<iInputCompSize;++icomp)
    {
      for (ic=0;ic<iNumCols;++ic)
      {
        Device::Param newParam = inputCompositeParamVecMap[baseTag][ic][icomp];
        paramName.str("");
        paramName << paramBase << ic << "." << newParam.uTag();
        newParam.setTag(paramName.str());

        convertedCompositeParams.push_back(newParam);
        paramMetadataExistMap[paramName.str()] = true;
      }
    }

    modelData.params.insert
      (modelData.params.end(),
       convertedCompositeParams.begin(),
       convertedCompositeParams.end());

    if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
    {
      Xyce::dout() << "After adding composite params:"<<std::endl;
      for (const auto &modelDataParam : modelData.params)
      {
        Xyce::dout() << modelDataParam;
      }
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setParameterValues
// Purpose       : Look for expression valued parameters in the parameter
//                 list, evaluate expressions found and reset the parameter
//                 value accordingly.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/19/2001
//-----------------------------------------------------------------------------
void ParameterBlock::setParameterValues(CircuitContext* contextPtr)
{
  Device::Param* parameterPtr;
  int numParameters = expressionValuedParams_.size();
  int i;

  for ( i = 0; i < numParameters; ++i )
  {
    parameterPtr = findParameter(expressionValuedParams_[i]);
    if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
    {
      Xyce::dout() << " setParameterValues, processing expressionValuedParams[" << i << "]" << std::endl;
      if (parameterPtr == NULL)
      {
        Xyce::dout() << "ParameterBlock::setParameterValues.  expressionValuedParams_["<<i<<"].uTag = " <<
          expressionValuedParams_[i].uTag() << " is NOT found in model data!" << std::endl;
        Report::DevelFatal() << expressionValuedParams_[i].uTag() << " is NOT found in model data!" << std::endl;
      }
    }
    
    if (parameterPtr != NULL)
    {
#if 0
      if (!contextPtr->resolveParameter(*parameterPtr))
#else
      resolveStatus paramResolveStatus;
      contextPtr->resolveParameter(*parameterPtr,paramResolveStatus);
      if (!(paramResolveStatus.success))
#endif
      {
        Util::Expression & expr = parameterPtr->getValue<Util::Expression>();

        Report::UserFatal message;
        message.at(modelData.getNetlistLocation());
        message << "Parameter " << parameterPtr->uTag() << " for model " << getName();

        if ( !( expr.getLeadCurrents().empty()) ) // ERK how does this make sense?
        {
          message << " contains illegal use of lead current: ";
        }
        else
        {
          message << " contains unrecognized symbols: ";
        }
        message << expr.get_expression();
      }
    }
    if (DEBUG_IO && isActive(Diag::IO_DEVICE_PARAMETERS))
    {
      Xyce::dout() << " after resolution, "  << std::endl;
      Xyce::dout() << "   Tag is " << parameterPtr->uTag() << std::endl;
      if (parameterPtr->hasExpressionValue())
      {
        Xyce::dout() << " Expression did not fully resolve, still has value "
                     << parameterPtr->getValue<Util::Expression>().get_expression() << std::endl;
      }
      else
      {
        Xyce::dout() << "   value is " << parameterPtr->stringValue() << std::endl;
      }
    }
  }

  std::vector<Device::Param> & modelDataParams = modelData.params;
  int modelDataSize = modelDataParams.size();
  for (int iparam=0;iparam<modelDataSize;++iparam)
  {
    Device::Param & param = modelDataParams[iparam];

    if (param.getType() == Xyce::Util::STR && !param.isNumeric())
    {
      ExtendedString paramNameOrig(param.stringValue());

      // since we're going to upcase paramNameOrig to try resolving it as an
      // expression, save the original just in case we *can't* resolve it
      // and need to restore it.  This is necessary because the string value
      //  might later be used as a file name, for example, and upcasing the
      // value is Just Wrong.
      ExtendedString paramNameSave(paramNameOrig);

      paramNameOrig.toUpper();
      if (paramNameOrig.possibleParam())
      {
        param.setVal(std::string("{" + paramNameOrig + "}"));

        // try to resolve this string as a simple parameter.  If it can't
        // resolve, restore it to its original value (it's real original
        // value, not its upcased value)
#if 0 
        if (!contextPtr->resolveParameter(param))
#else
        resolveStatus paramResolveStatus;
        contextPtr->resolveParameter(param,paramResolveStatus);
        if (!(paramResolveStatus.success))
#endif
          param.setVal(std::string(paramNameSave));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::findParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
Device::Param* ParameterBlock::findParameter( Device::Param const& parameter )
{
  std::vector<Device::Param>::iterator paramIter = std::find_if( modelData.params.begin(), modelData.params.end(), Util::EqualParam(parameter.tag()));
  if ( paramIter != modelData.params.end() )
  {
    return &(*paramIter);
  }
  else
  {
    return NULL;
  }
}

NetlistLocation ParameterBlock::netlistLocation() const {
  return NetlistLocation(modelData.getNetlistLocation());
}

} // namespace IO

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
int Pack<IO::ParameterBlock>::packedByteCount(
  const IO::ParameterBlock &    parameter_block)
{
  int byteCount = 0;
  int size, j;

  // count modelData
  byteCount += Xyce::packedByteCount(parameter_block.modelData);

  // count defaultApplied_
  byteCount += sizeof( int );

  // count expressionValuedParams_
  size = parameter_block.expressionValuedParams_.size();
  byteCount += sizeof( int );
  for( j = 0; j < size; ++j )
  {
    byteCount += Pack<Device::Param>::packedByteCount(parameter_block.expressionValuedParams_[j]);
  }

  // count the input composite size
  byteCount += sizeof( int ); // size

  // count the internals of the composite, if not empty.
  if ( !(parameter_block.inputCompositeParamVecMap.empty()) )
  {
    std::map< std::string, std::vector<std::vector<Device::Param> > >::const_iterator  iter;
    std::map< std::string, std::vector<std::vector<Device::Param> > >::const_iterator  begin = parameter_block.inputCompositeParamVecMap.begin();
    std::map< std::string, std::vector<std::vector<Device::Param> > >::const_iterator  end = parameter_block.inputCompositeParamVecMap.end();
    for (iter=begin;iter!=end;++iter)
    {
      // count paramBase
      std::string paramBase = iter->first;
      byteCount += sizeof( int ); // size
      byteCount += paramBase.length();

      // count size of paramVecVec.
      const std::vector<std::vector<Device::Param> > & tmpParamVecVec = iter->second;
      int sizeTmpParamVec = tmpParamVecVec.size();
      byteCount += sizeof( int ); // size

      for (int ivec=0;ivec< sizeTmpParamVec; ++ivec)
      {
        const std::vector<Device::Param> & tmpVec = tmpParamVecVec[ivec];
        int vecSize = tmpVec.size();
        byteCount += sizeof( int ); // size

        for (int ip=0; ip<vecSize; ++ip)
        {
          byteCount += Pack<Device::Param>::packedByteCount(tmpVec[ip]);
        }
      }
    }
  }

  return byteCount;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::pack
// Purpose       : Packs parameter block into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
void
Pack<IO::ParameterBlock>::pack(const IO::ParameterBlock &parameter_block, char * buf, int bsize, int & pos, Parallel::Communicator* comm )
{
  int size, j;
  int def;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+Pack<IO::ParameterBlock>::packedByteCount( parameter_block );
#endif

  // pack modelData;
  Xyce::pack(parameter_block.modelData, buf, bsize, pos, comm );

  if (parameter_block.defaultApplied_)
    def = 1;
  else
    def = 0;
  comm->pack( &def, 1, buf,  bsize, pos );

  // pack expressionValuedParams_
  size = parameter_block.expressionValuedParams_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for( j = 0; j < size; ++j )
  {
    Pack<Device::Param>::pack(parameter_block.expressionValuedParams_[j], buf, bsize, pos, comm );
  }

  // pack inputCompositeParamVecMap
  int sizeComposite=parameter_block.inputCompositeParamVecMap.size();
  comm->pack( &sizeComposite, 1, buf, bsize, pos );

  if ( !(parameter_block.inputCompositeParamVecMap.empty()) )
  {
    std::map< std::string, std::vector<std::vector<Device::Param> > >::const_iterator  iter;
    std::map< std::string, std::vector<std::vector<Device::Param> > >::const_iterator  begin = parameter_block.inputCompositeParamVecMap.begin();
    std::map< std::string, std::vector<std::vector<Device::Param> > >::const_iterator  end = parameter_block.inputCompositeParamVecMap.end();
    for (iter=begin;iter!=end;++iter)
    {
      // pack paramBase name.
      std::string paramBase = iter->first;
      int baseLength = paramBase.length();
      comm->pack( &baseLength, 1, buf, bsize, pos );
      comm->pack( paramBase.c_str(), baseLength, buf, bsize, pos );

      // pack size of tmpParamVecVec
      const std::vector<std::vector<Device::Param> > & tmpParamVecVec = iter->second;
      int sizeTmpParamVec = tmpParamVecVec.size();
      comm->pack( &sizeTmpParamVec, 1, buf, bsize, pos );

      // pack contents of tmpParamVecVec
      for (int ivec=0;ivec< sizeTmpParamVec; ++ivec)
      {
        const std::vector<Device::Param> & tmpVec = tmpParamVecVec[ivec];
        int vecSize = tmpVec.size();
        comm->pack( &vecSize, 1, buf, bsize, pos );

        for (int ip=0; ip<vecSize; ++ip)
        {
          Pack<Device::Param>::pack(tmpVec[ip], buf, bsize, pos, comm);
        }
      }
    }
  }

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    Report::DevelFatal() << "ParameterBlock::pack - predicted pos " << predictedPos
                         << ") does not match actual pos (" << pos << ")";
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : ParameterBlock::unpack
// Purpose       : Unpacks parameter block from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
template<>
void Pack<IO::ParameterBlock>::unpack(
  IO::ParameterBlock &          parameter_block,
  char *                        pB,
  int                           bsize,
  int &                         pos,
  Parallel::Communicator *      comm)
{
  int size = 0;
  int length = 0;
  int def;

  // unpack optionData
  Xyce::unpack(parameter_block.modelData, pB, bsize, pos, comm );

  comm->unpack( pB, bsize, pos, &def, 1 );
  if (def != 0)
    parameter_block.defaultApplied_=true;
  else
    parameter_block.defaultApplied_=false;

  // unpack expressionValuedParams_
  comm->unpack( pB, bsize, pos, &size, 1 );
  for (int j = 0; j < size; ++j )
  {
    parameter_block.expressionValuedParams_.push_back( Device::Param() );
    Xyce::unpack(parameter_block.expressionValuedParams_[j], pB, bsize, pos, comm );
  }

  // unpack inputCompositeParamVecMap
  int sizeComposite=0;
  comm->unpack( pB, bsize, pos, &size, 1 );
  sizeComposite=size;

  if ( sizeComposite>0 )
  {
    for (int ic=0;ic<sizeComposite;++ic)
    {
      // unpack paramBase name.
      comm->unpack( pB, bsize, pos, &length, 1 );
      std::string paramBase = std::string( ( pB + pos ), length );
      pos += length;

      // unpack size of paramVecVec
      int sizeTmpParamVec(0);
      comm->unpack( pB, bsize, pos, &sizeTmpParamVec, 1 );
      std::vector<std::vector<Device::Param> > tmpParamVecVec;
      tmpParamVecVec.resize(sizeTmpParamVec);

      for (int ivec=0;ivec< sizeTmpParamVec; ++ivec)
      {
        // unpack size of paramVec
        comm->unpack( pB, bsize, pos, &length, 1 );
        int vecSize = length;
        tmpParamVecVec[ivec].resize(vecSize);

        for (int ip=0; ip<vecSize; ++ip)
        {
          Xyce::unpack(tmpParamVecVec[ivec][ip], pB, bsize, pos, comm );
        }
      }
      parameter_block.inputCompositeParamVecMap[paramBase] = tmpParamVecVec;
    }
  }
}

} // namespace Xyce

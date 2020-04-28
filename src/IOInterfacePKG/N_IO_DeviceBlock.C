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

//-----------------------------------------------------------------------------
//
// Purpose        : Define the N_IO_DeviceBlock class instantiations of which
//                  are associated with netlist device lines.
//
// Special Notes  : ERK.  It seems that the name "N_IO_InstanceBlock" would have been
//                  more appropriate and less confusing, or possibly
//                  N_IO_InstanceParametersBlock.  Calling it "device block"
//                  makes it sound much more general than it really is.
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/09/2001
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <sstream>
#include <iostream>

#include <N_DEV_SourceData.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CircuitMetadata.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_Report.h>
#include <N_TOP_NodeDevBlock.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace IO {

//--------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 09/09/2001
//--------------------------------------------------------------------------
DeviceBlock::DeviceBlock(
  CircuitContext & cc,
  CircuitMetadata & md )
  : circuitContext_(cc),
    metadata_(md)
{
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::DeviceBlock
// Purpose       : constructor
// Special Notes : This constructor keeps a copy of the parsedLine_, which is only
//               : needed for the mutual inductors at this time.  No other device
//               : should use this.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
DeviceBlock::DeviceBlock(
  CircuitContext & cc,
  CircuitMetadata & md,
  std::string const& fileName,
  TokenVector const& parsedInputLine)
  : circuitContext_(cc),
    metadata_(md),
    parsedLine_(parsedInputLine),
    subcircuitInstance_(false),
    extracted_(false)
{
  deviceData_.getDevBlock().setNetlistLocation(NetlistLocation(fileName, parsedInputLine[0].lineNumber_));
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::DeviceBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
DeviceBlock::DeviceBlock(DeviceBlock const& rhsDB)
  : circuitContext_(rhsDB.circuitContext_),
    metadata_(rhsDB.metadata_),
    parsedLine_(rhsDB.parsedLine_),
    netlistType_(rhsDB.netlistType_),
    deviceData_(rhsDB.deviceData_),
    subcircuitInstance_(rhsDB.subcircuitInstance_),
    extracted_(rhsDB.extracted_)
{}


//-----------------------------------------------------------------------------
// Function      : DeviceBlock::print
// Purpose       : Output the details of a device block to standard out.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
void DeviceBlock::print()
{
  Xyce::dout() << std::endl << Xyce::section_divider << std::endl
               << "Device Information" << std::endl
               << "------------------" << std::endl
               << "device line:" << std::endl;
  int numFields = parsedLine_.size();
  int i;
  for ( i = 0; i < numFields; ++i )
  {
    Xyce::dout() << "  " << parsedLine_[i].string_;
  }
  Xyce::dout() << std::endl;
  Xyce::dout() << "  name: " << getInstanceName() << std::endl;

  Xyce::dout() << "  nodes: ";

  std::vector<std::string>::const_iterator paramIter;
  std::vector<std::string>::const_iterator paramEnd=getNodeValues().end();
  i=1;
  for (paramIter=getNodeValues().begin(); paramIter!=paramEnd; ++paramIter)
  {
    Xyce::dout() << "Node " << i << ": ";
    Xyce::dout() << *paramIter << " ";
    ++i;
  }

  Xyce::dout() << std::endl;

  if ( getModelName() != "" )
  {
    Xyce::dout() << "  model name: " << getModelName() << std::endl;
  }
  Xyce::dout() << std::endl;

  if ( getNumberOfInstanceParameters() > 0 )
  {
    Xyce::dout() << "  Instance Parameters:" << std::endl;
    size_t numParams = getNumberOfInstanceParameters();
    for ( size_t k = 0; k < numParams; ++k )
    {
      Xyce::dout() << "    " << getInstanceParameter(k).uTag();
      Xyce::dout() << "    " << getInstanceParameter(k).stringValue();

      if ( getInstanceParameter(k).given() )
      {
         Xyce::dout() << "    given";
      }
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;
  }

  Xyce::dout() << std::endl << Xyce::section_divider << std::endl;
  Xyce::dout() << std::endl;
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::clear
// Purpose        : Reset all attributes to their defaults.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/06/2003
//----------------------------------------------------------------------------
void DeviceBlock::clear()
{
  parsedLine_.clear();
  netlistType_ = "";
  subcircuitInstance_ = false;

  setModelName("");
  deviceData_.clear();
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractData
// Purpose       : Extract the device data from parsed line. Use device
//                 metadata to determine device type, number of nodes, etc.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractData( std::string const& fileName,
                               TokenVector const& parsedInputLine,
                               bool resolveParams,
                               bool modelBinning)
{
  // Set the device name and netlist type.
  ExtendedString deviceName(parsedInputLine[0].string_);
  deviceName.toUpper();
  setName(deviceName);
  setNetlistType(deviceName[0]);
  deviceData_.getDevBlock().setNetlistLocation(NetlistLocation(fileName, parsedInputLine[0].lineNumber_));

  if (!isValidDeviceType(getNetlistDeviceType()))
  {
    Report::UserError().at(fileName, parsedInputLine[0].lineNumber_)
      << "Invalid device type for device " << deviceName;
    return false;
  }

  bool result = false;

  // Invoke appropriate method to handle the rest of the line. Independent
  // sources and mutual inductances need to be handled differently than
  // other devices.
  if ( getNetlistDeviceType() == "V" || getNetlistDeviceType() == "I" || getNetlistDeviceType() =="P" )
  {
              
    if ( getNetlistDeviceType() =="P" )
      setNetlistType ( 'V' );

    result = Device::extractSourceData(parsedInputLine, *this, metadata_.getPrimaryParameter(getNetlistDeviceType(), -1));
  }
  else if ( getNetlistDeviceType() == "E" || getNetlistDeviceType() == "G" ||
            getNetlistDeviceType() == "F" || getNetlistDeviceType() == "H" )
  {
    result = extractBehavioralDeviceData( parsedInputLine );
  }
  else if ( getNetlistDeviceType() == "K" )
  {
    result = extractMutualInductanceData( parsedInputLine );
  }
  else if ( getNetlistDeviceType() == "S" || getNetlistDeviceType() == "W" )
  {
    result = extractSwitchDeviceData( parsedInputLine );
  }
  else if ( getNetlistDeviceType() == "X" )
  {
    result = extractSubcircuitInstanceData( parsedInputLine );
  }
  else if ( getNetlistDeviceType() == "Y" )
  {
    result = extractYDeviceData( parsedInputLine );
  }
  else if ( getNetlistDeviceType() == "U" )
  {
    result = extractUDeviceData( parsedInputLine );
  }
  else
  {
    result = extractBasicDeviceData( parsedInputLine, resolveParams, modelBinning );
  }

  // Check the status and extracting device data
  // and return if something went wrong.
  if ( !result )
    return false;

  //allows testing later to avoid "reextracting" data
  extracted_ = true;

  if (resolveParams)
  {
    // Now that the data has been extracted, given the circuit context,
    // the parameter values can be set.
    setParameterValues();
  }

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractSubcircuitInstanceData
// Purpose       : Extract the subcircuit instance data given on a netlist 'X'
//                 line.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/28/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractSubcircuitInstanceData( const TokenVector & parsedInputLine )
{
  // Set the flag indicating this device is a subcircuit instance.
  subcircuitInstance_ = true;

  // Set the device name and netlist type.
  ExtendedString fieldES( parsedInputLine[0].string_ );
  fieldES.toUpper();
  setName( fieldES );
  setNetlistType( fieldES[0] );

  // Set the number of fields on the line.
  size_t numFields = parsedInputLine.size();

  // Look for "PARAMS:" on the line, set the parameter start position and
  // subcircuit name position accordingly.
  int parameterPosition = 0;
  int subcircuitNamePosition = 0;
  ExtendedString field("");
  for ( size_t i = 1; i < numFields; ++i )
  {
    field = parsedInputLine[i].string_;
    field.toUpper();
    if ( field == "PARAMS:" )
    {
      parameterPosition = i+1;
      subcircuitNamePosition = i-1;
      break;
    }
    else if ( i < numFields-1 && parsedInputLine[i+1].string_ == "=" )
    {
      parameterPosition = i;
      subcircuitNamePosition = i-1;
      break;
    }
  }

  // If the "PARAMS:" keyword was not found, then set the
  // subcircuitNamePos to the last field on the line.
  if ( parameterPosition == 0 )
  {
    subcircuitNamePosition = numFields - 1;
  }

  // Use the DeviceBlock getModelName() field to hold the subcircuit name field.
  fieldES = parsedInputLine[ subcircuitNamePosition ].string_;
  fieldES.toUpper();
  setModelName( fieldES );

  // Set list of subcircuit external nodes, the last node
  // on the line precedes the subcircuit name.
  for ( int i = 1; i < subcircuitNamePosition; ++i )
  {
    addNodeValue( ExtendedString(parsedInputLine[i].string_).toUpper() );
  }

  // Collect the parameters on the line if there are any.
  if ( parameterPosition > 0 )
  {
    size_t i = parameterPosition;
    Device::Param parameter("","");
    while ( i < numFields )
    {
      fieldES =  parsedInputLine[i].string_;
      fieldES.toUpper();
      if (!fieldES.possibleParam())
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[i].lineNumber_)
          << "Parameter name " << fieldES << " contains illegal character(s)";
        return false;
      }
      parameter.setTag( fieldES );

      if ( parsedInputLine[i+1].string_ != "=" )
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[i].lineNumber_)
          << "Equal sign required between parameter and value in subcircuit parameter list for subcircuit instance " << getInstanceName();
        return false; // prevents core dump from next increment of i
      }

      i+=2; // Advance past "=" sign

      fieldES =  parsedInputLine[i].string_;
      fieldES.toUpper();

      if (fieldES.possibleParam())
      {
        fieldES = "{" + parsedInputLine[i].string_ + "}";
      }
      parameter.setVal(  std::string(fieldES) );

      addInstanceParameter( parameter );

      // Advance to next field.
      ++i;
    }
  }

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractBasicDeviceData
// Purpose       : Extract the device data from parsedLine for devices other
//                 than independent sources, mutual inductances and "Y" devices.
//                 These devices are exected to have input lines of the
//                 following format:
//
//    dname node node [node]* [Value] [model-name] [param=Value]* [Area]
//    +     [ON|OFF] [IC = value [,value]*] [TEMP=Value]
//
//                  where fields enclosed in [] indicates are optional, *
//                  indicates a field that may be reapeated 0 or more times,
//                  fields given in all lower case a string is expected, fields
//                  starting with a capital letter indicate a numerical value
//                  is expected, and fields in all upper case are expected
//                  literally.
//
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/24/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractBasicDeviceData( const TokenVector & parsedInputLine , bool failIfUnresolved, bool modelBinning)
{
  bool result;
  int i, n_start, n_end, n_req, n_opt, n_fill;

  // The device name has been extracted by extractData. The
  // remaining data in parsedLine is extracted here.

  int numFields = parsedInputLine.size();

  // Extract the model name from parsedLine if one exists. If
  // a model name was found, find its type.
  int modelLevel, modelNamePosition;
  std::string modelType;
  ModelFoundState modelFound = extractModelName( parsedInputLine, modelType, modelLevel, modelNamePosition, modelBinning );
  if (modelFound == MODEL_NOT_FOUND)
    return false;
  else if (modelFound == MODEL_NOT_SPECIFIED)
    modelLevel = -1;

  // Some devices require a model, check that a model was found if it is
  // required for this device.
  if (metadata_.isModelRequired(getNetlistDeviceType(), modelLevel) && modelFound == MODEL_NOT_SPECIFIED)
  {
    Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
      << "Model is required for device " << getInstanceName() << " and no valid model card found.";
    return false;
  }

  std::string netlistModelType("");
  if ( modelFound == MODEL_FOUND && modelNamePosition < numFields)
  {
    // Query the metadata for the validity of the device model.
    if (!metadata_.isModelTypeValid(getNetlistDeviceType(), modelType, modelLevel))
    {
      Report::UserError().at(getNetlistFilename(), parsedInputLine[modelNamePosition].lineNumber_)
        << "Model type \"" << modelType << "\" not valid for device " << getInstanceName() << " of type " << getNetlistDeviceType();
      return false;
    }

    netlistModelType = modelType;
  }

  if ( netlistModelType == "" )
  {
    netlistModelType = "default";
  }

  // Extract the device nodes from parsedInputLine.
  result = extractNodes(parsedInputLine, modelLevel, modelNamePosition);
  if (!result)
    return result;

  // Process optional nodes
  n_req = metadata_.getNumberOfNodes(getNetlistDeviceType(), modelLevel);
  int nodesFound = n_req;
  n_opt = metadata_.getNumberOfOptionalNodes(getNetlistDeviceType(), modelLevel);
  std::string primaryDeviceParameter(
    metadata_.getPrimaryParameter(getNetlistDeviceType(), modelLevel));

  if (n_opt > 0)
  {
    n_start = n_req+1;
    n_end = n_start;
    n_fill = metadata_.getNumberOfFillNodes(getNetlistDeviceType(), modelLevel);
    if (modelFound == MODEL_FOUND && modelNamePosition > n_req)
    {
      n_end = modelNamePosition;
    }
    else if (modelFound == MODEL_NOT_SPECIFIED && primaryDeviceParameter == "")
    {
      for (i=n_start ; i<numFields ; ++i)
      {
        if (parsedInputLine[i].string_ == "=")
        {
          n_end = i - 1;
          break;
        }
      }
    }
    if (n_end > n_start)
    {
      for (i=n_start ; i<n_end ; ++i)
      {
        addNodeValue( ExtendedString(parsedInputLine[i].string_).toUpper());
        ++nodesFound;
      }
    }
    for (i=n_end ; i<=n_req+n_fill ; ++i) {
      addNodeValue( "0" );
    }
  }
  deviceData_.getDevBlock().numExtVars = nodesFound;

  // Check for value field and set position of start of instance
  // parameters on the line.
  Device::Param primaryParameter( "", "" );
  int parameterStartPosition;
  if ( primaryDeviceParameter != "" )
  {
    // test for un-tagged parameter after model name and use as offset
    int offsetUntaggedParamAfterModel = 0;
    if ( modelFound == MODEL_FOUND && modelNamePosition + 1 < numFields &&
     !metadata_.isDeviceParameter( getNetlistDeviceType(), modelLevel,
     parsedInputLine[ modelNamePosition + 1 ].string_ ) )
    {
      offsetUntaggedParamAfterModel = 1;
    }

    if ( ((getNumberOfNodes() + 1 < numFields) && modelFound == MODEL_NOT_SPECIFIED) ||
         (modelFound == MODEL_FOUND && (getNumberOfNodes()+1 < modelNamePosition)) ||
         (modelFound == MODEL_FOUND && (getNumberOfNodes()+1 == modelNamePosition) &&
         bool(offsetUntaggedParamAfterModel)) )
    {
      if (metadata_.isDeviceParameter( getNetlistDeviceType(), modelLevel,
            parsedInputLine[getNumberOfNodes()+1+offsetUntaggedParamAfterModel].string_) &&
          !(numFields < getNumberOfNodes()+3+offsetUntaggedParamAfterModel ||
           parsedInputLine[getNumberOfNodes()+2+offsetUntaggedParamAfterModel].string_ != "="))
      {
        // If the field following the node list is in the instance
        // parameter list then it is not the value field.
        parameterStartPosition = getNumberOfNodes() + 1+offsetUntaggedParamAfterModel;
      }
      else
      {
        // The field following the node list is the value field.
        primaryParameter.set(
                       primaryDeviceParameter,
                       parsedInputLine[getNumberOfNodes()+1+offsetUntaggedParamAfterModel].string_ );

        // The value field must either be a valid number or an expression.
        if (primaryParameter.getType() == Xyce::Util::STR && !primaryParameter.isNumeric())
        {
          ExtendedString p_orig(primaryParameter.stringValue());
          p_orig.toUpper();
          if (p_orig.possibleParam())
          {
            primaryParameter.setVal(std::string("{" + p_orig + "}"));
            if (!circuitContext_.resolveParameter(primaryParameter))
              primaryParameter.setVal(std::string(p_orig));
          }
        }
        if ( failIfUnresolved && !primaryParameter.isNumeric() &&
             !primaryParameter.hasExpressionValue() )
        {
          Report::UserError().at(getNetlistFilename(), parsedInputLine[getNumberOfNodes()+1+offsetUntaggedParamAfterModel].lineNumber_)
            << "Illegal value found for device " <<  getInstanceName();
        }

        if ( modelFound == MODEL_FOUND)
        {
          parameterStartPosition = getNumberOfNodes() + 3;
        }
        else
        {
          parameterStartPosition = getNumberOfNodes() + 2;
        }
      }
    }
    else if ( modelFound == MODEL_FOUND && (getNumberOfNodes()+1 == modelNamePosition) )
    {
      // No device type parameter on line.
      parameterStartPosition = modelNamePosition + 1;
    }
    else if ( getNumberOfNodes()+1 >= numFields )
    {
      Report::UserWarning().at(getNetlistFilename(), parsedInputLine[getNumberOfNodes()].lineNumber_)
        << "Expected value field for device " << getInstanceName() << ", continuing with value of 0";

      // There are no instance parameters given for this device,
      // set this variable for consistency.
      parameterStartPosition = getNumberOfNodes() + 1;
    }
  }
  else
  {
    // The device does not have a value field.
    if ( modelFound == MODEL_FOUND)
    {
      parameterStartPosition = modelNamePosition + 1;
    }
    else
    {
      parameterStartPosition = getNumberOfNodes() + 1;
    }
  }

  // Add the device instance parameters and their default values
  // to instanceParameters and check parsedLine for instance
  // parameters.
  int parameterEndPosition;
  extractInstanceParameters( parsedInputLine,
                             parameterStartPosition,
                             parameterStartPosition,
                             parameterEndPosition, "break",
                             modelLevel);

  // If the primary device parameter was found, reset its value.
  if ( primaryParameter.tag() != "" )
  {
    Device::Param *parameterPtr = findInstanceParameter( primaryParameter );
    if (parameterPtr)
    {
      // parameterPtr->setVal( static_cast<Util::Param &>(primaryParameter));
      setParamValue(*parameterPtr, primaryParameter);
      parameterPtr->setGiven( true );
    }
    else
    {
      if ( getNetlistDeviceType() == "L")
      {
        addInstanceParameter( Device::Param (primaryParameter.tag(), parsedInputLine[modelNamePosition-1].string_));
      }
    }
  }

  // Check for AREA, OFF, IC and TEMP parameters.

  // Check for Area parameter.
  int linePosition = parameterEndPosition + 1;
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedInputLine[linePosition].string_ );
    field.toUpper();
    if ( field != "ON" && field != "OFF" && field != "IC" && field != "TEMP"  )
    {
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("AREA", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( parsedInputLine[linePosition].string_ );
        parameterPtr->setGiven( true );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
         // Either the device does not have an area parameter, or there is
         // some unrecognized parameter on the line. Either way, flag an error
         // on the line and stop.
         issueUnrecognizedParameterError(field);
      }
    }
  }

  // Check for ON parameter.
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedInputLine[linePosition].string_ );
    field.toUpper();
    if ( field == "ON" )
    {
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("ON", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( 1.0 );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
        // There is no ON parameter for this device.
        issueUnrecognizedParameterError("ON");
      }
    }
  }

  // Check for OFF parameter.
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedInputLine[linePosition].string_ );
    field.toUpper();
    if ( field == "OFF" )
    {
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("OFF", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( 1.0 );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
        // There is no OFF parameter for this device.
        issueUnrecognizedParameterError("OFF");
      }
    }
  }


  // Check for temperature parameter.
  if ( linePosition < numFields )
  {
    ExtendedString field ( parsedInputLine[linePosition].string_ );
    field.toUpper();
    if ( field == "TEMP" )
    {
      linePosition += 2; // advance past "=" sign
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("TEMP", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( parsedInputLine[linePosition].string_ );
        parameterPtr->setGiven( true );
        ++linePosition; // advance to next field in parsedLine
      }
      else
      {
         // This device does not have a TEMP parameter.
         issueUnrecognizedParameterError("TEMP");
      }
    }
  }

  // Issue fatal error if there are more fields left on the line.
  if ( linePosition < numFields )
  {
    Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition].lineNumber_)
      << "Unrecognized fields for device " << getInstanceName();
    return false;
  }

  return true; // Only get here on success.
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::extractBehavioralDeviceData
// Purpose        : Extract the device data for a "E", "F", "G", and "H" devices. If
//                  the usage of these devices is as a PSpice style behavioral
//                  device (i.e. a dependent source device), convert the device
//                  to a Xyce B-source device.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 12/09/2002
//----------------------------------------------------------------------------
bool DeviceBlock::extractBehavioralDeviceData( const TokenVector & parsedInputLine )
{
  bool result;

  int numFields = parsedInputLine.size();

  // Check the device line for a "VALUE", "TABLE", or "POLY" field. If
  // none of these fields are found and the device is an E or G device,
  // it can be parsed by the usual device rules (extractBasicDeviceData).
  // If there is a "VALUE", "TABLE", or "POLY" field it should be the
  // fourth field on the line.
  ExtendedString typeField("");
  if ( numFields >= 4 )
  {
    typeField = parsedInputLine[3].string_;
    typeField.toUpper();
  }

  // For HSPICE compatibility, VOL is a synonym for VALUE for the E source.
  // CUR is a synonym for VALUE for the G-source.
  if ( ((getNetlistDeviceType() == "E") && (typeField != "VALUE" && typeField != "VOL" && typeField != "TABLE" && typeField != "POLY")) ||
       ((getNetlistDeviceType() == "G") && (typeField != "VALUE" && typeField != "CUR" && typeField != "TABLE" && typeField != "POLY")) )
  {
    bool resolveParams=true;
    bool modelBinning=false;
    result = extractBasicDeviceData( parsedInputLine, resolveParams, modelBinning);
    return result;
  }

  // For all other cases, convert the device to a B-source.

  // Get the device type.
  std::string deviceType(getNetlistDeviceType());

  // Change the type of the device to a B-source.
  setNetlistType ( 'B' );
  deviceData_.getDevBlock().bsourceFlag = true;

  // Extract the device nodes from parsedInputLine.
  result = extractNodes(parsedInputLine, -1, 0);
  if (!result)
    return result;

  std::string expression("");
  if ((deviceType == "E" || deviceType == "G") && typeField != "POLY")
  {
    // Get the expression from the line. The expression may need to be
    // transformed by adding braces where Xyce expects them.
    if ( (typeField == "VALUE") || (deviceType == "E" && typeField == "VOL") ||
         (deviceType == "G" && typeField == "CUR") )
    {
      // For HSPICE compatibility, VOL is a synonym for VALUE for the E source.
      // CUR is a synonym for VALUE for the G-source.
      int exprStart = 4;
      if ( parsedInputLine[exprStart].string_ == "=" ) ++exprStart;

      for ( int i = exprStart; i < numFields; ++i )
      {
//        expression += parsedInputLine[i].string_;

        ExtendedString field = parsedInputLine[i].string_;
        field.toUpper();


        if (field != "SMOOTHBSRC" &&  field != "RCCONST" )
          expression += parsedInputLine[i].string_;
        else
        {

          int parameterStartPosition = i;
          int parameterEndPosition;
          extractInstanceParameters( parsedInputLine,
                             parameterStartPosition,
                             parameterStartPosition,
                             parameterEndPosition, "break",
                             -1 );
          break;
        }


      }
    }
    else if ( typeField == "TABLE" )
    {
      expression = "TABLE";

      // Find the next "=" sign.
      int equalPos;
      for ( equalPos = 5; equalPos < numFields; ++equalPos )
      {
        if ( parsedInputLine[equalPos].string_ == "=" ) break;
      }

      if ( equalPos == numFields )
      {
        // Problem, did not find expected "=" sign.
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Required \"=\" sign missing after the expression in the TABLE for device " << getInstanceName();
        return false;
      }

      // Get the table expression.
      std::string tableExpression("");
      for ( int i = 4; i < equalPos; ++i )
      {
        tableExpression += parsedInputLine[i].string_;
      }

      // Add braces to tableExpression if needed.
      if ( tableExpression[0] != '{' &&
        tableExpression[tableExpression.size()-1] != '}' )
      {
        tableExpression = "{" + tableExpression + "}";
      }

      expression += " " + tableExpression + " = ";

      // Add the table data to the expression.
      for ( int i = equalPos+1; i < numFields; ++i )
      {
        expression += " " + parsedInputLine[i].string_;
      }
    }
  }
  else if (typeField == "POLY")
  {
    // Convert an E, F, G, or H device with a POLY function to the equivalent
    // B-source device.

    // Get the dimension of the polynomial.
    if (parsedInputLine.size() < 6)
    {
      Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
        << "Not enough fields on input line for device " << getInstanceName();
      return false;
    }

    Device::Param dimensionParam("dim", parsedInputLine[5].string_);
    int dimension = dimensionParam.getImmutableValue<int>();

    // Build up the B-source POLY expression, some fields need to
    // be manipulated depending on the device type.
    expression = "POLY(";
    expression += dimensionParam.stringValue() + ")";

    int linePosition = 7;
    if (deviceType == "E" || deviceType == "G")
    {
      if (parsedInputLine.size() < 7+2*dimension)
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Not enough fields on input line for device " << getInstanceName();
        return false;
      }

      // Get the POLY controlling nodes, there should be as many pairs
      // of these as the dimension of the polynomial.
      for (int i = 0; i < dimension; ++i, linePosition += 2)
      {
        std::string cnode1(parsedInputLine[linePosition].string_);
        std::string cnode2(parsedInputLine[linePosition+1].string_);

        if (cnode2 == "0")
        {
          expression += " V(" + cnode1 + ")";
        }
        else
        {
          expression += " V(" + cnode1 + "," + cnode2 + ")";
        }
      }
    }
    else
    {
      // This is an F or an H device.
      if (parsedInputLine.size() < 7+dimension)
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Not enough fields on input line for device " << getInstanceName();
        return false;
      }

      // Get the controlling sources, there should be as many of these as
      // as the dimenstion of the polynomial.
      for (int i = 0; i < dimension; ++i, ++linePosition)
      {
        std::string csource(parsedInputLine[linePosition].string_);

        expression += " I(" + csource + ")";
      }
    }

    // Add all remaining fields to the expression.
    for ( ;linePosition < parsedInputLine.size(); ++linePosition)
    {
      expression += " " + parsedInputLine[linePosition].string_;
    }
  }
  else if (deviceType == "F" || deviceType == "H")
  {
    // Linear F or H conversion to B-source.
    if (parsedInputLine.size() < 5)
    {
      Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
        << "Not enough fields on input line for device " << getInstanceName();
      return false;
    }

    // Get the source and value.
    std::string source(parsedInputLine[3].string_);
    std::string value(parsedInputLine[4].string_);

    // If the user has specified the gain as an expression, it might have
    // braces around it.  Replace them with parens
    if (value[0] == '{' && value[value.size()-1] == '}')
    {
      value[0] = '(';
      value[value.size()-1]=')';
    }
    expression = value + " * I(" + source + ")";
  }
  else
  {
    // Have something unrecognizeable on the line.
  }

  // Added braces to expression if needed.
  if ( expression[0] != '{' &&
      expression[expression.size()-1] != '}' )
  {
    expression = "{" + expression + "}";
  }

  // Set the device parameter
  Device::Param parameter( "", "" );
  if ( deviceType == "E" || deviceType == "H" )
  {

    
      Device::Param* parameterPtr =
        findInstanceParameter( Device::Param("V", "") );
      if ( parameterPtr != NULL )
      {
        parameterPtr->setVal( expression );
        parameterPtr->setGiven( true );
      }
      else
      {

      parameter.setTag( "V" );
      parameter.setVal( expression );
      parameter.setGiven( true );
      addInstanceParameter( parameter );

      parameter.setTag( "I" ); // This B-source parameter is required but
                             // won't be used.
      parameter.setVal( "" );
      parameter.setGiven( false );
      addInstanceParameter( parameter );
    }
    
  }
  else if ( deviceType == "F" || deviceType == "G" )
  {
    parameter.setTag( "I" );
    parameter.setVal( expression );
    parameter.setGiven( true );
    addInstanceParameter( parameter );

    parameter.setTag( "V" ); // This B-source parameter is required but
                             // won't be used.
    parameter.setVal( "" );
    parameter.setGiven( false );
    addInstanceParameter( parameter );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractYDeviceData
// Purpose       : Extract the device data from parsedLine for "Y" devices.
//                 These devices are exected to have input lines of the following
//                 format:
//
//                 Ytype name node [node]* [Value] [model-name] [param=Value]*
//
//                 where fields enclosed in [] indicates are optional, *
//                 indicates a field that may be reapeated 0 or more times,
//                 fields given in all lower case a string is expected, fields
//                 starting with a capital letter indicate a numerical value
//                 is expected, and fields in all upper case are expected
//                 literally.
//
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/24/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractYDeviceData( const TokenVector & parsedInputLine )
{
  bool result;

  // Copy the parsed line since it will be modified.
  parsedLine_ = parsedInputLine;

  int numFields = parsedLine_.size();

  if ( numFields < 2 )
  {
    Report::UserError().at(getNetlistFilename(), parsedLine_[0].lineNumber_)
      << "Not enough fields on input line for device " << getInstanceName();
    return false;
  }

  // Reset the device name and type which were set by extractData based on
  // the standard rules for netlist devices.
  ExtendedString deviceName(parsedLine_[1].string_);
  deviceName.toUpper();
  ExtendedString deviceType ( parsedLine_[0].string_.substr(1,parsedLine_[0].string_.size()) );
  deviceType.toUpper();
  setNetlistType( deviceType );
  setName( "Y" + deviceType + "!" + deviceName );

  // drop the first element of the parsed line vector
  // and treat this like a basic device
  parsedLine_.erase(parsedLine_.begin());

  // special handling for new mutual inductor
  if( deviceType == "MIL" || deviceType == "MIN" )
  {
    result = extractMIDeviceData( parsedLine_ );
  }
  else
  {
    bool resolveParams=true;
    bool modelBinning=false;
    result = extractBasicDeviceData( parsedLine_, resolveParams, modelBinning );
  }

  return result;

}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractUDeviceData
// Purpose       : Extract the device data from parsedLine for "U" devices.
//                 These devices are exected to have input lines of the following
//                 format:
//
//                 Uname type(N) node [node]* [Value] [model-name] [param=Value]*
//
//                 where fields enclosed in [] indicates are optional, *
//                 indicates a field that may be reapeated 0 or more times,
//                 fields given in all lower case a string is expected, fields
//                 starting with a capital letter indicate a numerical value
//                 is expected, and fields in all upper case are expected
//                 literally.  Some digital devices (e.g., NAND, AND, OR and
//                 NOR) may have a variable number of inputs.  In that case,
//                 the number of inputs is denoted by (for example) AND(N),
//                 where N is an integer.  Other gates (e.g, INV) have a fixed
//                 number of inputs, and their type is (for example) INV.
//
// Special Notes : The Xyce U device is intended to be similar to the PSpice U 
//                 digital devices, rather than the SPICE3F5 U device which is
//                 the "uniform RC transmission line".  
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 02/14/2014
//-----------------------------------------------------------------------------
bool DeviceBlock::extractUDeviceData( const TokenVector & parsedInputLine )
{
  bool noError = true;

  // Copy the parsed line since it will be modified.
  parsedLine_ = parsedInputLine;

  int numFields = parsedLine_.size();

  if ( numFields < 2 )
  {
    Report::UserError().at(getNetlistFilename(), parsedLine_[0].lineNumber_)
      << "Not enough fields on input line for U device " << getInstanceName();
    noError = false;
  }

  // Parse AND(N) syntax where AND is the gate type, and N is the number
  // of inputs.  Error handling for N not being an integer, N missing or
  // trailing ')' missing.
  std::string numInputs = "";
  if ( numFields >=5)
  {
    if ( parsedLine_[2].string_ == "(" )
    {
      if ( parsedLine_[4].string_ == ")" )
      {
        // found possible (N) syntax
        numInputs = parsedLine_[3].string_;
        if (!Util::isInt(numInputs))
	{
          Report::UserError().at(getNetlistFilename(), parsedLine_[0].lineNumber_)
            << "Found (N) syntax in U device " << getInstanceName() << " but N not an integer";
          noError = false;
        }
        else
	{
          // found valid (N) syntax.  Remove (N) from parsed line so that
          // extractBasicDeviceData() function works correctly
          parsedLine_.erase(parsedLine_.begin()+2,parsedLine_.begin()+5);
	}
      }
      else
      {
        Report::UserError().at(getNetlistFilename(), parsedLine_[0].lineNumber_)
          << "Apparent error parsing number of inputs, (N) syntax, for U device " << getInstanceName();
        noError = false;
      }
    } 
  }

  // Reset the device name and type which were set by extractData based on
  // the standard rules for netlist devices.  Note that the order of the
  // deviceType and deviceName are reversed from the Y-device syntax
  ExtendedString deviceType(parsedLine_[1].string_);
  deviceType.toUpper();
  ExtendedString deviceName( parsedLine_[0].string_.substr(1,parsedLine_[0].string_.size()) );
  deviceName.toUpper();
  // U devices are not allowed to have a '!' character in them because it hoses
  // subsequent processing in the digital-device constructor.
  if ( (deviceType.find_first_of('!') != std::string::npos) || 
       (deviceName.find_first_of('!') != std::string::npos) )
  {
    Report::UserError().at(getNetlistFilename(), parsedLine_[0].lineNumber_)
      << "U device name: " << deviceName << " of type: " << deviceType <<  " should not have ! character.";
    noError = false;  
  }
  setNetlistType( deviceType );
  if (numInputs == "")
  {
    // name string format for U devices with a fixed number of inputs
    setName( "U" + deviceType + "!" + deviceName );
  }
  else
  {
    // name string format for U devices with a variable number of inputs
    setName( "U" + deviceType + "!" + deviceName + "!" + numInputs);
  }

  // drop the second element of the parsed line vector, which is the
  // gate type, and treat this like a basic device.  Previous statement 
  // already had error message: "not enough fields on input line for device"
  if ( numFields >= 2 )
  {
    parsedLine_.erase(parsedLine_.begin()+1);
  }

  // remove the U from the device name, before extracting the device data.
  // this mimics the Y device parsing
  parsedLine_[0].string_ = parsedLine_[0].string_.substr(1,parsedLine_[0].string_.size()); 

  bool resolveParams=true;
  bool modelBinning=false;
  return (noError == true) ? extractBasicDeviceData( parsedLine_, resolveParams, modelBinning ) : false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractMutualInductanceData
// Purpose       : Extract the device data from parsedLine for mutual
//                 inductances.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/25/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractMutualInductanceData( const TokenVector & parsedInputLine )
{
  int numFields = parsedInputLine.size();
  bool kequals = false;

  // Extract the inductor names.
  int numInductors = 0;
  Device::Param parameter( "", "" );
  while ( parsedInputLine[numInductors+1].string_[0] == 'L' ||
          parsedInputLine[numInductors+1].string_[0] == 'l' ||
          parsedInputLine[numInductors+1].string_[0] == 'X' ||
          parsedInputLine[numInductors+1].string_[0] == 'x'
          )
  {
    //Have to add some additional checking if the inductor is part of a
    //subcircuit instance.

    if (parsedInputLine[numInductors+1].string_[0] == 'X' ||
          parsedInputLine[numInductors+1].string_[0] == 'x'
        )
    {
      //For now, we're just not going to allow this.  Fix later with a new
      //parser, hopefullY:
      Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
        << "Subcircuit calls ('X' devices) are not allowed in mutual inductor definitions.";
      return false;
    }
    else
    {
      ++numInductors;

      parameter.setTag( parsedInputLine[numInductors].string_ );
      parameter.setVal( numInductors );
      addInstanceParameter( parameter );
    }
  }

  parameter.setTag( "COUPLING" );

  //Have to check whether coupling is present in the form "k = {coupling}"

  if ((parsedInputLine[numInductors + 1].string_ == "K") ||
      ( (parsedInputLine[numInductors + 1].string_ == "k")
        && (parsedInputLine[numInductors + 2].string_ == "=")
        && (numInductors + 3 < numFields) ) )
  {
    kequals = true;
    parameter.setVal( parsedInputLine[numInductors + 3].string_ );
  }
  else
  {
    parameter.setVal( parsedInputLine[numInductors + 1].string_ );
  }

  addInstanceParameter( parameter );

  int modelLevel, modelNamePosition;
  std::string modelType;
  ModelFoundState modelFound = extractModelName( parsedInputLine, modelType, modelLevel, modelNamePosition );

  // check format for errors; bogus Lnames are handled later
  if ( modelFound == MODEL_FOUND)
  {
    if( kequals )
    {
      if ( ( modelNamePosition != ( numInductors + 4 ) ) ||
           ( modelNamePosition != ( numFields - 1 ) ) )
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Malformed line for device " << getInstanceName();
      }
    }
    else
    {
      if ( ( modelNamePosition != ( numInductors + 2 ) ) ||
           ( modelNamePosition != ( numFields - 1 ) ) )
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Malformed line for device " << getInstanceName();
      }
    }
  }
  else 
  {
    if (kequals)
    {
      if ( ( numInductors + 3 ) != ( numFields - 1 ) )
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Specified model not found for device " << getInstanceName();
      }
    }
    else
    {
      if ( ( numInductors + 1 ) != ( numFields - 1 ) )
      {
        Report::UserError().at(getNetlistFilename(), parsedInputLine[0].lineNumber_)
          << "Specified model not found for device " << getInstanceName();
      }
    }
  }

  return true;  // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractSwitchDeviceData
// Purpose       : Convert SPICE switch formats to general expression controlled
//                 switch
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/16/04
//-----------------------------------------------------------------------------
bool DeviceBlock::extractSwitchDeviceData( const TokenVector & parsedInputLine )
{
  bool result;
  int modelLevel, modelNamePosition, controlPosition;
  std::string modelType, expression;
  int numFields = parsedInputLine.size();
  int i, parameterStartPosition, parameterEndPosition;
  bool is_w=false;

  if ( getNetlistDeviceType() == "W" )
  {
    is_w = true;
    setNetlistType ( 'S' );
  }

  ModelFoundState modelFound = extractModelName( parsedInputLine, modelType, modelLevel, modelNamePosition );
  if (modelFound != MODEL_FOUND)
  {
    Report::UserError().at(getNetlistFilename(), getLineNumber())
      << "No model found for switch device " << getInstanceName();
    return false;
  }
  controlPosition = 0;
  for (i=modelNamePosition+1 ; i<numFields ; ++i)
  {
    if (ExtendedString(parsedInputLine[i].string_).toUpper() == "CONTROL")
    {
      controlPosition = i;
      break;
    }
  }

  if ( is_w )
  {
    if ( modelNamePosition != 4)
    {
      Report::UserError().at(getNetlistFilename(), getLineNumber())
          << "Wrong number of nodes for switch device " << getInstanceName();
      return false;
    }
    else if (controlPosition != 0)
    {
      Report::UserError().at(getNetlistFilename(), getLineNumber())
          << "CONTROL param found for W switch device " << getInstanceName();
      return false;
    }
      
    setNetlistType ( 'S' );
    expression = "{I(" + parsedInputLine[modelNamePosition-1].string_ + ")}";
  }
  else
  {
    if ((modelNamePosition != 3 && controlPosition > 0) ||
        (modelNamePosition != 5 && controlPosition == 0))
    {
      Report::UserError().at(getNetlistFilename(), getLineNumber())
          << "Wrong number of nodes for switch device " << getInstanceName();
      return false;
    }
    if (modelNamePosition == 5)
      expression = "{V(" + parsedInputLine[3].string_ + ")-V(" + parsedInputLine[4].string_ + ")}";
  }

  addNodeValue( ExtendedString(parsedInputLine[1].string_).toUpper());
  addNodeValue( ExtendedString(parsedInputLine[2].string_).toUpper());
  parameterStartPosition = modelNamePosition+1;
  parameterEndPosition = numFields;
  extractInstanceParameters( parsedInputLine,
                             parameterStartPosition,
                             parameterStartPosition,
                             parameterEndPosition,
                             "break");

  if (controlPosition == 0)
  {
    Device::Param parameter;
    int numParameters = getNumberOfInstanceParameters();
    for ( i = 0; i < numParameters; ++i )
    {
      if ( getInstanceParameter(i).uTag() == "CONTROL")
      {
        parameter = getInstanceParameter(i);
        parameter.setVal( expression );
        parameter.setGiven( true );
        setInstanceParameter( i, parameter );
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/27/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::extractNodes(const TokenVector & parsedInputLine,
                               int modelLevel, int modelNamePosition)
{
  size_t numFields = parsedInputLine.size();
  int numNodes;

  numNodes = metadata_.getNumberOfNodes(getNetlistDeviceType(), modelLevel);
  if (numNodes == -1)
    return false;

  const size_t nodeStartPos = 1;
  const size_t nodeEndPos = nodeStartPos + numNodes - 1;

  if ( modelNamePosition > 0 && modelNamePosition <= nodeEndPos)
  {
    if (metadata_.isModelTypeValid(getNetlistDeviceType(), parsedInputLine[modelNamePosition].string_, modelLevel))
      Report::UserError().at(getNetlistFilename(), getLineNumber())
        << "Insufficient nodes specified or node name '" << parsedInputLine[modelNamePosition].string_ << "' matches one of this device's model name";
    else
      Report::UserError().at(getNetlistFilename(), getLineNumber())
      << "Insufficient number of nodes specified";

    return false;
  }

  if ( numFields < nodeEndPos + 1)
  {
    Report::UserError().at(getNetlistFilename(), getLineNumber())
      << "Not enough fields on input line for device " << getInstanceName();
    return false;
  }

  std::vector<std::string> nodeValues;
  for ( int i = nodeStartPos; i <= nodeEndPos; ++i )
  {
    nodeValues.push_back( ExtendedString(parsedInputLine[i].string_).toUpper() );
  }
  setNodeValues(nodeValues);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function       : getLandW
// Purpose        : 
// Special Notes  : Added to support model binning
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//-----------------------------------------------------------------------------
bool getLandW (
  const TokenVector & parsedInputLine,
  CircuitContext & circuitContext_, 
  double & L,
  double & W
    )
{
  size_t model_search_begin_index = 1;
  size_t model_search_end_index = parsedInputLine.size();

  bool LWfound=false, foundL=false, foundW=false;
  for (size_t LWsearchFieldno = model_search_begin_index; LWsearchFieldno < model_search_end_index; ++LWsearchFieldno )
  {
    std::string tmp = parsedInputLine[LWsearchFieldno].string_;
    Util::toUpper(tmp);
    if (tmp == "L" && LWsearchFieldno+2<model_search_end_index)
    {
      Device::Param par(std::string("L"), parsedInputLine[LWsearchFieldno+2].string_);
      foundL = circuitContext_.fullyResolveParam(par,L);
    }
    if (tmp == "W" && LWsearchFieldno+2<model_search_end_index)
    {
      Device::Param par(std::string("W"), parsedInputLine[LWsearchFieldno+2].string_);
      foundW = circuitContext_.fullyResolveParam(par,W);
    }
  }

  if ( foundL && foundW )
  { 
    LWfound=true;
  }

  return LWfound;
}

//-----------------------------------------------------------------------------
// Function       : getLandNFIN
// Purpose        : 
// Special Notes  : Added to support model binning
// Scope          :
// Creator        : Eric R. Keiter
// Creation Date  : 10/2/2018
//-----------------------------------------------------------------------------
bool getLandNFIN (
  const TokenVector & parsedInputLine,
  CircuitContext & circuitContext_, 
  double & L,
  double & NFIN
    )
{
  size_t model_search_begin_index = 1;
  size_t model_search_end_index = parsedInputLine.size();

  bool LNFINfound=false, foundL=false, foundNFIN=false;
  for (size_t LNFINsearchFieldno = model_search_begin_index; LNFINsearchFieldno < model_search_end_index; ++LNFINsearchFieldno )
  {
    std::string tmp = parsedInputLine[LNFINsearchFieldno].string_;
    Util::toUpper(tmp);
    if (tmp == "L" && LNFINsearchFieldno+2<model_search_end_index)
    {
      Device::Param par(std::string("L"), parsedInputLine[LNFINsearchFieldno+2].string_);
      foundL = circuitContext_.fullyResolveParam(par,L);
    }
    if (tmp == "NFIN" && LNFINsearchFieldno+2<model_search_end_index)
    {
      Device::Param par(std::string("NFIN"), parsedInputLine[LNFINsearchFieldno+2].string_);
      foundNFIN = circuitContext_.fullyResolveParam(par,NFIN);
    }
  }

  if ( foundL && foundNFIN )
  { 
    LNFINfound=true;
  }

  return LNFINfound;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractModelName
// Purpose       : Check parsedLine for existance of model name. If a model
//                 name exists, set its position in parsedLine and the device
//                 model name and return true, otherwise return false.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/27/2001
//-----------------------------------------------------------------------------
ModelFoundState
DeviceBlock::extractModelName(
  const TokenVector & parsedInputLine,
  std::string & modelType,
  int &         modelLevel,
  int &         modelNamePosition,
  bool modelBinning )
{
  modelType = "";
  modelNamePosition = 0;
  modelLevel = 0;

  const std::string &device_type = getNetlistDeviceType();

  // Start looking for model name right after device name.
  size_t model_search_begin_index = 1;

  size_t model_search_end_index = parsedInputLine.size();

  ModelFoundState modelFound = MODEL_NOT_SPECIFIED;


  // Model Binning stuff:
  bool LWfound=false;
  bool LNFINfound=false;
  double L = 0.0, W = 0.0, NFIN=0.0;
  if (modelBinning && device_type == "M")
  {
    LWfound = getLandW (parsedInputLine, circuitContext_, L, W);
    LNFINfound = getLandNFIN (parsedInputLine, circuitContext_, L, NFIN);
  }

  for (size_t fieldno = model_search_begin_index; fieldno < model_search_end_index; ++fieldno )
  {
    const std::string &model_name = parsedInputLine[fieldno].string_;
    if (device_type == "K") 
    {
      if (model_name == "=" && !equal_nocase(parsedInputLine[fieldno - 1].string_, "K"))
        break;
    }
    else 
    {
      if (model_name == "=") // equal sign implies parameter is being parsed, we've gone to far
        break;
    }

    ParameterBlock* modelPtr;

    bool foundTheModel = circuitContext_.findModel(model_name, modelPtr);

    // Model Binning stuff:
    std::string final_model_name = model_name;
    bool foundTheBinnedModel = false;
    if ( modelBinning && !foundTheModel && (LWfound || LNFINfound))
    {
      if (device_type == "M")
      {
        std::string binNumber;
        foundTheBinnedModel = circuitContext_.findBinnedModel(model_name, modelPtr, LWfound, LNFINfound, L, W, NFIN, binNumber);
        final_model_name = model_name + "." + binNumber;
      }
    }

    if (foundTheModel || foundTheBinnedModel)
    {
      if ( device_type != "K")
      {
        modelType = modelPtr->getType();
        modelLevel = modelPtr->getLevel();

        const DeviceMetadata &model_metadata = metadata_.getDeviceMetadata(device_type, modelLevel);

        // Only consider models that are of valid type for this device type!
        if (model_metadata.isModelLevelValid() &&
            model_metadata.isModelTypeValid(modelType))
        {          
          // Now make absolutely sure that it is even legal for a model
          // of this type/level to appear at this position!
          
          const DeviceMetadata &model_level_metadata = metadata_.getDeviceMetadata(device_type,modelLevel);

          if (fieldno >= model_level_metadata.numNodes+1)
          {

            modelPtr->addDefaultModelParameters(metadata_);
            modelFound = MODEL_FOUND;
            modelNamePosition = fieldno;
            setModelName(final_model_name);
            break;
          }
          else
          {
            // Ignore it --- it is not a candidate for this device's model
            // due to appearing too early on the line
          }
        }
        else
        {
          // Ignore it --- it is not a candidate for this device's model
          // due to not having the right model type or level.  But if is
          // the right type and no valid level, error out.
          if (model_metadata.isModelTypeValid(modelType) &&
              !model_metadata.isModelLevelValid() )
          {
            Report::UserError().at(getNetlistFilename(), getLineNumber())
              << "Model type \"" << modelType << "\" does not have level " << modelLevel << " defined";
            return MODEL_NOT_FOUND;
          }
        }
      }
      else
      {
        // Apparently we just accept any model name here for K devices
        modelFound = MODEL_FOUND;
        modelType = modelPtr->getType();
        modelLevel = modelPtr->getLevel();
        modelNamePosition = fieldno;
        setModelName(model_name);
      }
    }
  }

  return modelFound;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractInstanceParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 09/27/2001
//-----------------------------------------------------------------------------
void DeviceBlock::extractInstanceParameters( const TokenVector & parsedInputLine,
                                             int searchStartPosition,
                                             int & parameterStartPosition,
                                             int & parameterEndPosition,
                                             std::string const& action,
                                             int modelLevel )
{
  addDefaultInstanceParameters(modelLevel);

  // Check for instance parameters on the device line, reset the parameter
  // value in instanceParameters to the value given.
  int linePosition = searchStartPosition;
  bool foundParameters = false;
  int numFields = parsedInputLine.size();
  parameterStartPosition = numFields; // Default in case there are no
                                      // parameters on the line.
  parameterEndPosition = numFields; // Default in case there are no
                                      // parameters on the line.
  Device::Param paramToFind;
  while (linePosition < numFields)
  {
    // Look for the parameter in the instance parameter list.
    paramToFind.setTag(parsedInputLine[linePosition].string_);
    Device::Param* parameterPtr = findInstanceParameter(paramToFind);

    if ( parameterPtr != NULL )
    {

      ExtendedString es(parsedInputLine[linePosition].string_);
      es.toUpper();

      if ( !foundParameters)
      {
        if ( es != "PORT" && es != "Z0")
        {
          foundParameters = true;
          parameterStartPosition = linePosition;
        }
      }

      if ( es == "ON" || es == "OFF")
      {
        // If the ON or OFF field is untagged set its value to 1 (true)
        // and continue to next field.
        if ( !(linePosition + 1 < numFields &&
             parsedInputLine[linePosition + 1].string_ == "=") )
        {
          parameterPtr->setVal( 1.0 );
          parameterPtr->setGiven( true );
          ++linePosition;
          continue;
        }
      }

      ++linePosition;    // advance to next field in parsedLine

      if ( linePosition >= numFields )
      {
        // Hit end-of-line unexpectedly.  Report error, and set the parameterEndPosition
        // variable correctly, so that subsequent parsing is correct.
        Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition-1].lineNumber_)
          << "Unexpectedly reached end of line while looking for parameters for device " << getInstanceName();
        parameterEndPosition = linePosition - 1;
        return;
      }

      if ( parsedInputLine[linePosition].string_ == "=" )
      {
        ++linePosition; // if field is "=", advance to next field
      }

      if (parameterPtr->stringValue() != "VECTOR" &&
          parameterPtr->stringValue() != "VECTOR-COMPOSITE")
      {
        if ( linePosition >= numFields )
        {
          // Hit end-of-line unexpectedly.  Report error, and set the parameterEndPosition
          // variable correctly, so that subsequent parsing is correct.
          Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition-1].lineNumber_)
            << "Unexpectedly reached end of line while looking for parameters for device " << getInstanceName();
          parameterEndPosition = linePosition - 1;
          return;
        }

        if (DEBUG_IO)
          Xyce::dout() << " Setting parameter " << parameterPtr->uTag()
                       << "to value " << parsedInputLine[linePosition].string_ << std::endl;

        // Set the parameter value.
        if (parameterPtr->getType() == Xyce::Util::DBLE)
        {
          const std::string & tmpStr = (parsedInputLine[linePosition].string_);
          if (Util::possibleParam(tmpStr))
          {
            parameterPtr->setVal( "{" + parsedInputLine[linePosition].string_ + "}" );
          }
          else if (Util::isValue(tmpStr))
          {
            parameterPtr->setVal( Util::Value(tmpStr) );
          }
          else
          {
            parameterPtr->setVal( parsedInputLine[linePosition].string_ );
          }
        }
        else if (parameterPtr->getType() == Xyce::Util::INT)
        {
          const std::string & tmpStr = (parsedInputLine[linePosition].string_);

          if (Util::possibleParam(tmpStr))
          {
            parameterPtr->setVal( "{" + parsedInputLine[linePosition].string_ + "}");
          }
          else if (Util::isInt(tmpStr))
          {
            parameterPtr->setVal( Util::Ival(tmpStr) );
          }
          else
          {
            parameterPtr->setVal( parsedInputLine[linePosition].string_ );
          }
        }
        else if(parameterPtr->isTableFileTypeQuoted())
        {
          circuitContext_.resolveTableFileType(*parameterPtr);
        }
        else
        {
          parameterPtr->setVal( parsedInputLine[linePosition].string_ );
        }
        parameterPtr->setGiven( true );

        ++linePosition; // Advance to next field.
      }
      else if (parameterPtr->stringValue() == "VECTOR")
      {
        std::ostringstream paramName;
        Device::Param parameter;
        std::string paramBaseName(parameterPtr->uTag());
        int j = 1;
        paramName << paramBaseName << j;
        // check linePosition before trying to set the parameter
        if ( linePosition >= numFields )
        {
          // Hit end-of-line unexpectedly.  Report error, and set the parameterEndPosition
          // variable correctly, so that subsequent parsing is correct.
          Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition-1].lineNumber_)
            << "Unexpectedly reached end of line while looking for parameters for device " << getInstanceName();
          parameterEndPosition = linePosition - 1;
          return;
        }
        parameter.set(paramName.str(), parsedInputLine[linePosition].string_);
        ExtendedString p_orig(parameter.stringValue());
        p_orig.toUpper();
        if (p_orig.possibleParam())
        {
          parameter.setVal(std::string("{" + p_orig + "}"));
          if (!circuitContext_.resolveParameter(parameter))
            parameter.setVal(std::string(p_orig));
        }
        parameter.setGiven( true );
        addInstanceParameter(parameter);
        ++linePosition;

        while ( (linePosition < numFields) && (parsedInputLine[linePosition].string_ == ",") )
        {
          paramName.str("");
          ++j;
          ++linePosition;
          paramName << paramBaseName << j;
          // check linePosition before trying to set the parameter
          if ( linePosition >= numFields )
          {
            // Hit end-of-line unexpectedly.  Report error, and set the parameterEndPosition
            // variable correctly, so that subsequent parsing is correct.
            Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition-1].lineNumber_)
              << "Unexpectedly reached end of line while looking for parameters for device " << getInstanceName();
            parameterEndPosition = linePosition - 1;
            return;
          }
          parameter.set(paramName.str(), parsedInputLine[linePosition].string_);
          ExtendedString p_orig(parameter.stringValue());
          p_orig.toUpper();
          if (p_orig.possibleParam())
          {
            parameter.setVal(std::string("{" + p_orig + "}"));
            if (!circuitContext_.resolveParameter(parameter))
              parameter.setVal(std::string(p_orig));
          }
          parameter.setGiven( true );
          addInstanceParameter(parameter);
          ++linePosition;
        }
      }
      else if (parameterPtr->stringValue() == "VECTOR-COMPOSITE")
      {
        std::ostringstream paramName;
        // Get the components from metadata.
        std::vector<Device::Param> components;
        metadata_.getInstanceCompositeComponents(
            getNetlistDeviceType(),
            parameterPtr->uTag(), modelLevel,
            components);

        if (DEBUG_IO) {
          Xyce::dout() << " Processing composite for parameter " << parameterPtr->uTag() << std::endl;

          for (int ieric=0; ieric < components.size(); ++ieric)
          {
            Xyce::dout() << "tag[" << ieric << "] = " << components[ieric].uTag() << "  val = " << components[ieric].stringValue() << std::endl;
          }
        }

        // Do some error checking.
        // find the right curly brace, if it exists.
        // This is not a perfect test.  It just finds the first one.
        bool foundRightB = false;
        int rightBracePosition = linePosition;
        int ipos1;

        if (DEBUG_IO)
          Xyce::dout() << "    Doing error checking for vector-composite..." << std::endl;

        for (ipos1=linePosition;ipos1<numFields;++ipos1)
        {
          if (DEBUG_IO)
            Xyce::dout() << "    position..." << ipos1 << " string is "
                         << parsedInputLine[ipos1].string_ ;

          if (parsedInputLine[ipos1].string_ == "}")
          {
            foundRightB = true;
            rightBracePosition = ipos1;
            if (DEBUG_IO)
              Xyce::dout() << "    found right brace!" ;
          }
          if (DEBUG_IO)
            Xyce::dout() << std::endl;
        }

        if (DEBUG_IO)
          Xyce::dout() << "String at line position is " << parsedInputLine[linePosition].string_ << std::endl ;

        if (parsedInputLine[linePosition].string_ != "{" || !foundRightB)
        {
          // Expect components of VECTOR-COMPOSITE to be enclosed in braces.
          Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition].lineNumber_)
            << "Composite parameter " << parameterPtr->uTag() << " for device "  << getInstanceName()
            << " must be enclosed in braces {}";
          return;
        }
        ++linePosition;

        Device::Param parameter( "", "" );
        std::string paramBase(parameterPtr->uTag());
        int numComponents = 0;
        while (parsedInputLine[linePosition].string_ != "}")
        {
          ExtendedString component ( parsedInputLine[linePosition].string_ );
          component.toUpper();
          std::vector<Device::Param>::iterator paramIter =
            std::find_if(components.begin(), components.end(), Util::EqualParam(component));
          if (paramIter == components.end())
          {
            Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition].lineNumber_)
              << "Found unexpected component \"" << component << "\" for composite parameter " << paramBase
              << " for device "  << getInstanceName()
              << ". Instance line may be improperly formatted";
            return;
          }

          if (DEBUG_IO)
            Xyce::dout() << " Found component " << paramIter->uTag() << std::endl;

          linePosition += 2;

          // Mark the component as given in components. Later we will
          // add all of the components and their defaults that were
          // not given.
          paramIter->setGiven(true);

          int j = 0;
          bool startVC = true;
          while (startVC || parsedInputLine[linePosition].string_ == ",")
          {
            if (!startVC) ++linePosition;
            startVC = false;
            paramName << paramBase << j << "." << component;

            ExtendedString value ( parsedInputLine[linePosition].string_ );
            // Commented out by TVR on 9 Nov 2010
            // See similar comment in ParameterBlock::extractModelData
            // It is unreasonable to up-case parameter values here, because
            // the possibility exists that the string parameter value is
            // intended as a file name, and this will break on any system with
            // a case-dependent file system, such as Linux or BSD or any other
            // *nix other than Mac.
            // The responsibility for up-casing a parameter value should be
            // left to the user of the parameter if necessary, never done here.
            //
            //value.toUpper();

            if (value == "DEFAULT")
            {
              // parameter.set(paramName.str(), static_cast<Util::Param &>(*paramIter));
              setParam(parameter, paramName.str(), *paramIter);
            }
            else if ( (value == ",") || (value == "}") )
	    {
              // Catch the user error of inserting extra commas into the composite
              // parameter definition
              Report::UserError().at(getNetlistFilename(), parsedInputLine[linePosition].lineNumber_)
                << "Invalid syntax (likely spurious comma) for composite parameter " << paramBase
                << " for device "  << getInstanceName();
              return;
            }
            else
            {
              parameter.set(paramName.str(), std::string(value));
            }

            parameter.setGiven(true);
            addInstanceParameter(parameter);
            paramName.str("");

            ++linePosition;
            ++j;
          } // end of while loop.

          numComponents = j;
        } // end of while loop.

        // Now find the components not given in the netlist and
        // add them with their defaults.
        parameter.setGiven(false);
        std::vector<Device::Param>::iterator paramIter = components.begin();
        std::vector<Device::Param>::iterator paramEnd = components.end();
        for (; paramIter != paramEnd; ++paramIter)
        {
          if (!paramIter->given())
          {
            for (int j = 0; j < numComponents; ++j)
            {
              paramName << paramBase << j << "." << paramIter->uTag();
              // parameter.set(paramName.str(), static_cast<Util::Param &>(*paramIter));
              setParam(parameter, paramName.str(), *paramIter);
              addInstanceParameter(parameter);
              paramName.str("");
            }
          }
        } // end of paramIter loop
        // Advance past the right bracket, so that it will not be an
        // unrecognized instance parameter.
        ++linePosition;

      } // end of VECTOR-COMPOSITE if statement.
    }
    else //  parameterPtr != NULL
    {
      if ( action == "break" )
      {
        // Break out of loop and return if a field is encountered
        // that is not an instance parameter. Only care about
        // parameterEndPosition in this case.
        parameterEndPosition = linePosition - 1;
        break;
      }
      else if ( action == "continue" )
      {
        // Continue to the next field if a field is encountered that
        // is not an instance parameter.
        ++linePosition;
      }
    }
  }

  // Set parameterEndPosition if we get here with action = "break".
  if ( action == "break" )
  {
    parameterEndPosition = linePosition - 1;
  }
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::addDefaultInstanceParameters
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 08/11/2003
//----------------------------------------------------------------------------
void DeviceBlock::addDefaultInstanceParameters(int modelLevel)
{
  addInstanceParameters(metadata_.getInstanceParameters(getNetlistDeviceType(), modelLevel));
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::issueUnrecognizedParameterError
// Purpose       : Use the Xyce error manager to issue an error message to
//                 the user and exit the program.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 10/01/2001
//-----------------------------------------------------------------------------
void DeviceBlock::issueUnrecognizedParameterError(
   std::string const& parameterName)
{
  Report::UserError().at(getNetlistFilename(), getLineNumber())
    << "Unrecognized parameter " << parameterName << " for device " << getInstanceName();
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::setParameterValues
// Purpose       : Look for expression valued parameters in the parameter
//                 list, evaluate expressions found and reset the parameter
//                 value accordingly.
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 12/31/2001
//-----------------------------------------------------------------------------
bool DeviceBlock::setParameterValues()
{
  Device::Param parameter( "", "" );
  int numParameters = getNumberOfInstanceParameters();
  int i;
  std::vector<std::string> strings;
  std::vector<std::string> funcs;

  for ( i = 0; i < numParameters; ++i )
  {
    parameter = getInstanceParameter(i);
    if ( parameter.hasExpressionValue() || parameter.isQuoted() || parameter.isTableFileTypeQuoted()  || parameter.isStringTypeQuoted()  )
    {
      if (!circuitContext_.resolveParameter(parameter))
      {
        std::ostringstream msg;
        msg << "Parameter " << parameter.uTag() << " for device "
            << getInstanceName().getEncodedName() << " contains unrecognized symbol";
        if (parameter.getType() == Xyce::Util::EXPR)
        {
          Util::Expression &e = parameter.getValue<Util::Expression>();

          strings.clear();
          funcs.clear();
          e.get_names(XEXP_STRING, strings);
          e.get_names(XEXP_FUNCTION, funcs);
          if (strings.size() + funcs.size() == 1)
            msg << ":";
          else if (strings.size() + funcs.size() > 1)
            msg << "s:";
          for (std::vector<std::string>::iterator s = strings.begin(), s_end =strings.end(); s != s_end; ++s)
            msg << " " << *s;
          for (std::vector<std::string>::iterator s = funcs.begin(), s_end =funcs.end(); s != s_end; ++s)
            msg << " " << *s << "()";
        }
        else
        {
          msg << "(s): " + parameter.stringValue();
        }
        if (strings.size() + funcs.size() > 0)
        {
          Report::UserError().at(getNetlistFilename(), getLineNumber())
            << msg.str();
        }
      }
      setInstanceParameter( i, parameter );
    }
    else
    {
      if ( parameter.getType() == Xyce::Util::STR && !(parameter.isNumeric()) )
      {
        if (Util::possibleParam(parameter.stringValue()))
        {
          ExtendedString p_orig(parameter.stringValue()); p_orig.toUpper();

          parameter.setVal(std::string("{" + p_orig + "}"));
          if (!circuitContext_.resolveParameter(parameter))
            parameter.setVal(std::string(p_orig));
        }
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::findInstanceParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2001
//-----------------------------------------------------------------------------
Device::Param* DeviceBlock::findInstanceParameter( Util::Param const& parameter )
{
  std::vector<Device::Param>::iterator paramIter;

  paramIter = find_if(deviceData_.getDevBlock().params.begin(), deviceData_.getDevBlock().params.end(), Util::EqualParam(parameter.tag()));
  if ( paramIter != deviceData_.getDevBlock().params.end() )
  {
    return &(*paramIter);
  }
  else
  {
    return NULL;
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceBlock::findInstanceParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2001
//-----------------------------------------------------------------------------
Device::Param* DeviceBlock::findInstanceParameter(
      std::string const& parameterName )
{
  Device::Param parameter( parameterName, "" );

  return findInstanceParameter( parameter );
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::getInstanceParameters
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/10/2003
//----------------------------------------------------------------------------
void DeviceBlock::getInstanceParameters(
    std::vector<Device::Param>& parameters)
{
  parameters.reserve(deviceData_.getDevBlock().params.size());

  parameters.insert(parameters.end(),
                    deviceData_.getDevBlock().params.begin(),
                    deviceData_.getDevBlock().params.end());
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::extractMIDeviceData
// Purpose       : Extract the device data from parsedLine for YMIL/YMIN device
//
//                  Ytype!name inductor_list coupling_list model params
//
//                  type :- MIL | MIN
//                  for linear and nonlinear coupling respectively
//
//                  name :- K1_K2_...KN
//                  linear couplings with shared inductors have Knames
//                  concatenated w/underscores otherwise the Kname is used
//
//                  inductor_list :- Lname T1 T2 I
//                  a list of one or more inductors comprising the name,
//                  terminals, and inductance
//
//                  coupling_list :- L1 L2 ... LN C
//                  a list of inductor names followed by coupling value
//                  linear couplings with shared inductors will have several
//                  of these
//
//                  model :- the model name for nonlinear coupling
//
//                  params
//
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceBlock::extractMIDeviceData( const TokenVector & parsedInputLine )
{
  // The device name has been extracted by extractData. The
  // remaining data in parsedLine is extracted here.
  int numFields = parsedInputLine.size();

  // Extract the model name from parsedLine if one exists. If
  // a model name was found, find its type.
  int modelLevel, modelNamePosition;
  std::string modelType;
  ModelFoundState modelFound = extractModelName( parsedInputLine, modelType, modelLevel, modelNamePosition );
  if (modelFound == MODEL_NOT_FOUND)
    return false;
  else if (modelFound == MODEL_NOT_SPECIFIED)
    modelLevel = -1;

  // Some devices require a model, check that a model was found if it is
  // required for this device.
  if (metadata_.isModelRequired(getNetlistDeviceType(), modelLevel) && modelFound == MODEL_NOT_SPECIFIED)
  {
    Report::UserError().at(getNetlistFilename(), getLineNumber())
      << "Did not find required model for device " << getInstanceName();
    return false;
  }

  if (modelFound == MODEL_FOUND)
  {
    // Query the metadata for the validity of the device model.  Note that the
    // YMIL/YMIN model metadata is extracted from the K model metadata
    if( !metadata_.isModelTypeValid(getNetlistDeviceType(), modelType, modelLevel ) )
    {
      Report::UserError().at(getNetlistFilename(), getLineNumber())
        << "Model type \"" << modelType << "\" not valid for device " << getInstanceName() << " of type " << getNetlistDeviceType();
      return false;
    }
  }

  std::set< std::string > inductors;
  Device::Param param;

  // parameters names differ between YMIL and YMIN
  std::string tmpDiff;
  if (modelFound == MODEL_FOUND)
  {
    param.setTag( "NONLINEARCOUPLING" );
    param.setGiven( true );
    param.setVal( 1 );
    addInstanceParameter( param );
    tmpDiff = "COUPLEDMutIndNonLin";
  }
  else
  {
    tmpDiff = "COUPLEDMutIndLin";
  }

  // Extract the device nodes from parsedInputLine.
  std::vector<std::string> nodeValues;
  int nodesFound = 0, i = 1;

  // stop when first inductor name is seen again == start of coupling list
  std::string couplingListFlag(parsedInputLine[i].string_);

  do
  {
    // add inductor to mutual inductance
    param.setTag( tmpDiff );
    param.setVal( parsedInputLine[i].string_ );
    param.setGiven(true);
    addInstanceParameter( param );
    inductors.insert( parsedInputLine[i].string_ );

    // store the two node values and add to parameter list
    nodeValues.push_back( ExtendedString( parsedInputLine[i + 1].string_ ).toUpper() );
    param.setTag( "NODE1" );
    param.setVal( parsedInputLine[i + 1].string_ );
    param.setGiven(true);
    addInstanceParameter( param );

    nodeValues.push_back( ExtendedString( parsedInputLine[i + 2].string_ ).toUpper() );
    param.setTag( "NODE2" );
    param.setVal( parsedInputLine[i + 2].string_ );
    param.setGiven(true);
    addInstanceParameter( param );

    nodesFound += 2;

    // add inductance value for this inductor
    param.setTag( "COUPLEDINDUCTANCE" );
    param.setVal( ( parsedInputLine[i + 3].string_ ).c_str() );
    // The value field must either be a valid number or an expression.
    if (param.getType() == Xyce::Util::STR && !param.isNumeric())
    {
      ExtendedString p_orig(param.stringValue());
      p_orig.toUpper();
      if (p_orig.possibleParam())
      {
        param.setVal(std::string("{" + p_orig + "}"));
        if (!circuitContext_.resolveParameter(param))
          param.setVal(std::string(p_orig));
      }
    }
    if( !param.isNumeric() && !param.hasExpressionValue() )
    {
      Report::UserError()
        << "Illegal value found for device " <<  getInstanceName();
    }
    else
    {
      // At this point, the inductance really should be resolved
      if (param.hasExpressionValue())
      {
        if (!circuitContext_.resolveParameter(param))
        {
          Report::UserError()
            << "Could not resolve parameter " 
            << param.uTag() << " value of " 
            << param.stringValue()
            << " for device " <<  getInstanceName();
        }
      }
      else
      {
        double value=param.getImmutableValue<double>();
        param.setVal(value);
      }
    }
    param.setGiven(true);
    addInstanceParameter( param );
    
    param.setTag( "IC" );
    param.setVal( ( parsedInputLine[i + 4].string_ ).c_str() );
    param.setGiven(true);
    addInstanceParameter( param );
    
    // move to next inductor in the list
    i += 5;
  } while ( ( i < numFields ) && ( parsedInputLine[i].string_ != couplingListFlag ) );

  // shouldn't need to do this so I'll comment it out.
  // remove duplicate node names
  // nodeValues.unique();

  // set device node names
  setNodeValues( nodeValues );
  deviceData_.getDevBlock().numExtVars = nodesFound;

  bool moreCouplings = ( parsedInputLine[i].string_ == couplingListFlag );
  int curr = i;
  ++i;

  // retrieve coupling list
  while( moreCouplings && ( i < numFields ) )
  {
    if( inductors.find( parsedInputLine[i].string_ ) == inductors.end() )
    {
      // add coupling list
      param.setTag( "COUPLING" );
      const std::string &  coup = ((parsedInputLine[i].string_));
      if (Util::isValue(coup))
      {
        param.setVal(Util::Value(coup));
      }
      else
      {
        std::string exp("{"+ parsedInputLine[i].string_ + "}");
        param.setVal(exp);
      }

      param.setGiven(true);
      addInstanceParameter( param );

      for( ; curr < i; ++curr )
      {
        // add inductor to mutual coupling list
        param.setTag( "COUPLEDINDUCTOR" );
        param.setVal( parsedInputLine[curr].string_ );
        param.setGiven(true);
        addInstanceParameter( param );
      }

      // move indices to next list coupling value
      ++curr;
      ++i;

      // check for beginning of another coupling list
      moreCouplings = ( ( i + 1 < numFields ) &&
       ( inductors.find( parsedInputLine[curr].string_ ) != inductors.end() ) );
    }

    // move to next item on the list
    ++i;
  }


  // Add the device instance parameters and their default values
  // to instanceParameters and check parsedLine for instance
  // parameters.
  int pStart = (modelFound == MODEL_FOUND) ? i : curr;
  int pEnd = numFields - 1 ;

  // extract remaining instance parameters on the device line
  extractInstanceParameters( parsedInputLine, pStart, pStart, pEnd, "break", modelLevel );

  // Issue error if there are more fields on the line.
  if( pEnd + 1 < numFields )
  {
    Report::UserError().at(getNetlistFilename(), getLineNumber())
      << "Unrecognized fields for device " << getInstanceName();
    return false;
  }

  return true;  // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::isValidDeviceType
// Purpose       : Given a device type (a letter) return false
//                 if that letter represents  an  illegal, unimplemented
//                 device type
// Special Notes : Having added this function, there is now one extra place
//                 that has to be updated when adding a new device letter.
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 06/20/2008
//-----------------------------------------------------------------------------
bool DeviceBlock::isValidDeviceType(const std::string & deviceType)
{
  bool retcode;
  if (deviceType == "A" || deviceType == "N" )
//      deviceType == "N" ||
//      deviceType == "P" ) 
  {
    retcode=false;
  }
  else
  {
    retcode=true;
  }

  return (retcode);
}

//----------------------------------------------------------------------------
// Function       : DeviceBlock::setNodeValues
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 08/12/2003
//----------------------------------------------------------------------------
void DeviceBlock::setNodeValues(std::vector<std::string> const& nodeValues)
{
  deviceData_.set_NodeList(nodeValues);
}

//-----------------------------------------------------------------------------
// Function      : DeviceBlock::setInstanceParamter
// Purpose       :
// Special Notes : It assumed that getNumberOfInstanceParameters was called
//                 to ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/04/2002
//-----------------------------------------------------------------------------
void DeviceBlock::setInstanceParameter(int i, Device::Param & parameter)
{
  // deviceData_.getDevBlock().params[i].setVal(static_cast<Util::Param &>(parameter));
  setParamValue(deviceData_.getDevBlock().params[i], parameter);
  deviceData_.getDevBlock().params[i].setGiven( parameter.given() );
}

//-----------------------------------------------------------------------------
// Function      : extractSubcircuitModelName
// Purpose       : Extract ONLY the model name from a parsed input line.
// Special Notes : This non-member method is for extracting only the model name
//               : from the parsed input line without creating a DeviceModel object.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 06/16/2015
//-----------------------------------------------------------------------------
bool
extractSubcircuitModelName( const TokenVector& parsedInputLine, std::string& modelName )
{
  // Set the number of fields on the line.
  size_t numFields = parsedInputLine.size();
  
  // Look for "PARAMS:" on the line, set the parameter start position and
  // subcircuit name position accordingly.
  int parameterPosition = 0;
  int subcircuitNamePosition = 0;
  ExtendedString field("");
  for ( size_t i = 1; i < numFields; ++i )
  {
    field = parsedInputLine[i].string_;
    field.toUpper();
    if ( field == "PARAMS:" )
    {
      parameterPosition = i+1;
      subcircuitNamePosition = i-1;
      break;
    }
    else if ( i < numFields-1 && parsedInputLine[i+1].string_ == "=" )
    {
      parameterPosition = i;
      subcircuitNamePosition = i-1;
      break;
    }
  }
  
  // If the "PARAMS:" keyword was not found, then set the
  // subcircuitNamePos to the last field on the line.
  if ( parameterPosition == 0 )
  {
    subcircuitNamePosition = numFields - 1;
  }
  
  ExtendedString fieldES = parsedInputLine[ subcircuitNamePosition ].string_;
  fieldES.toUpper();
  modelName = fieldES;
 
  return true; 
}

} // namespace IO
} // namespace Xyce

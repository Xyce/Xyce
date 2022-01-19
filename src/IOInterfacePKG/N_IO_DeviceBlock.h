//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : Declare the N_IO_DeviceBlock class instantiations of which
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

#ifndef N_IO_DEVICEBLOCK_H
#define N_IO_DEVICEBLOCK_H

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_DEV_Param.h>

#include <N_TOP_NodeDevBlock.h>

namespace Xyce {
namespace IO {

enum ModelFoundState {MODEL_FOUND, MODEL_NOT_SPECIFIED, MODEL_NOT_FOUND};

class DeviceBlock
{
private:
  DeviceBlock();

public:
  DeviceBlock( CircuitContext & cc, CircuitMetadata & md );

  DeviceBlock(
     CircuitContext &                                                cc,
     CircuitMetadata &                                               md,
     const std::string &                                             fileName,
     const TokenVector &                                             parsedInputLine);

  DeviceBlock(const DeviceBlock &rhsDB );

private:
  DeviceBlock &operator=(const DeviceBlock &rhsDB );

public:
  ~DeviceBlock()
  {}

  // Setters and Getters
  const std::string & getNetlistFilename() const
  {
    return deviceData_.getDevBlock().getNetlistLocation().getFilename();
  }

  int getLineNumber() const
  {
    return deviceData_.getDevBlock().getNetlistLocation().getLineNumber();
  }

  void setName(const std::string &name )
  {
    deviceData_.getDevBlock().setInstanceName(Device::InstanceName(name));
  }

  void setNetlistType( char type )
  {
    netlistType_ = type;
  }

  void setNetlistType(const std::string &type)
  {
    netlistType_ = type;
  }

  void addNodeValue( std::string const& nodeValue )
  {
    deviceData_.addNode(nodeValue);
  }

  void setModelName( std::string const& modelName )
  {
    deviceData_.getDevBlock().setModelName(modelName);
    deviceData_.getDevBlock().modelFlag = (modelName != "");
  }

  void addInstanceParameter(const Device::Param & parameter )
  {
    deviceData_.getDevBlock().params.push_back( parameter );
  }

  void addInstanceParameters(const std::vector<Device::Param> &parameters )
  {
    deviceData_.getDevBlock().params.insert(deviceData_.getDevBlock().params.end(), parameters.begin(), parameters.end());
  }

  const Device::InstanceName &getInstanceName() const 
  {
    return deviceData_.getDevBlock().getInstanceName();
  }

  const Device::ModelName &getModelName() const 
  {
    return deviceData_.getDevBlock().getModelName();
  }

  const std::string getNetlistDeviceType() const 
  {
    return netlistType_;
  }

  int getNumberOfNodes() const 
  {
    return deviceData_.get_NodeList().size();
  }

  const std::vector<std::string> & getNodeValues() const 
  {
    return deviceData_.get_NodeList();
  }

  std::vector<std::string> & getNodeValues() 
  {
    return deviceData_.get_NodeList();
  }

  int getNumberOfInstanceParameters() const 
  {
    return deviceData_.getDevBlock().params.size();
  }

  const Device::Param &getInstanceParameter(int i) const 
  {
    return deviceData_.getDevBlock().params[i];
  }

  bool isSubcircuitInstance() const 
  {
    return subcircuitInstance_;
  }

  bool isExtracted() const 
  {
    return extracted_;
  }

  const Topo::NodeDevBlock &getDeviceData() const 
  {
    return deviceData_;
  }

  void setInstanceParameter(int i, Device::Param &parameter );
  void setNodeValues( std::vector<std::string> const& nodeValues );

  Device::Param* findInstanceParameter( Util::Param const& parameter );
  Device::Param* findInstanceParameter( std::string const& parameter );
  void getInstanceParameters(std::vector<Device::Param>& parameters);

  // Print the details of a device to standard out.
  void print();

  // Clear the device, reset all attributes to their defaults.
  void clear();

  bool extractData( const std::string & fileName, 
                    const TokenVector & parsedInputLine,
                    bool resolveParams,
                    bool modelBinning,
                    double scale);

  bool extractData(bool resolveParams, bool modelBinning,double scale)
  {
    return extractData( getNetlistFilename(), parsedLine_, resolveParams, modelBinning, scale ); 
  } 

  // Extract the subcircuit instance data given on a netlist 'X' line.
  bool extractSubcircuitInstanceData( const TokenVector & parsedInputLine );

  bool extractNodes(const TokenVector & parsedInputLine,
                    int modelLevel, int modelNamePosition);

  void extractInstanceParameters( const TokenVector & parsedInputLine,
                                  int searchStartPosition,
                                  int & parameterStartPosition,
                                  int & parameterEndPosition,
                                  std::string const& action,
                                  int modelLevel = -1 );

  void issueUnrecognizedParameterError(std::string const& parameterName);

private:
  bool extractMutualInductanceData( const TokenVector & parsedInputLine );

  //bool extractBasicDeviceData( bool failIfUnresolved=true)
  bool extractBasicDeviceData( bool failIfUnresolved, bool modelBinning, double scale)
  {
    //return extractBasicDeviceData( parsedLine_ ,failIfUnresolved);
    return extractBasicDeviceData( parsedLine_ ,failIfUnresolved, modelBinning, scale);
  }

  //bool extractBasicDeviceData( const TokenVector & parsedInputLine , bool failIfUnresolved=true);
  bool extractBasicDeviceData( const TokenVector & parsedInputLine , bool failIfUnresolved, bool modelBinning, double scale);

  bool extractBehavioralDeviceData( const TokenVector & parsedInputLine );

  bool extractYDeviceData( const TokenVector & parsedInputLine );

  bool extractUDeviceData( const TokenVector & parsedInputLine );

  bool extractMIDeviceData( const TokenVector & parsedInputLine );

  bool extractSwitchDeviceData( const TokenVector & parsedInputLine );

  int checkIfModelValid (const std::string & modelType, int modelLevel, int fieldNo);

  ModelFoundState extractModelName( const TokenVector & parsedInputLine,
                                    std::string& modelType, 
                                    int & modelLevel, 
                                    int & modelNamePosition,
                                    bool modelBinning,
                                    double scale);

  void addDefaultInstanceParameters(int modelLevel);


  // Look for expression valued parameters in the parameter
  // list, evaluate expressions found and reset the parameter
  // value accordingly.
  bool setParameterValues();

  bool isValidDeviceType(const std::string & deviceType);

private:
  CircuitContext &              circuitContext_;
  CircuitMetadata &             metadata_;

  TokenVector                   parsedLine_;

  std::string                   netlistType_;
  Topo::NodeDevBlock            deviceData_;

  bool                          subcircuitInstance_;
  bool                          extracted_;
};

// Extract only the model name from the subcircuit line.
// NOTE:  This is for the first pass of the parser to avoid creating and destroying DeviceModel objects.
bool extractSubcircuitModelName( const TokenVector& parsedInputLine, std::string& modelName );

} // namespace IO
} // namespace Xyce

#endif

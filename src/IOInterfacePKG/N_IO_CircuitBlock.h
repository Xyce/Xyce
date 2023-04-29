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
// Purpose        : Declare the circuit level containers for holding netlist
//                  circuit data and the associated circuit level methods.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/06/2001
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_CircuitBlock_h
#define Xyce_N_IO_CircuitBlock_h

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

// trilinos includes
#include <Teuchos_RCP.hpp>

#include <N_IO_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_DeviceBlock.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : CircuitBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/06/2001
//-----------------------------------------------------------------------------

class CircuitBlock
{
public:
  CircuitBlock(
     const std::string &                                                fileName,
     const CmdParse &                                                   command_line,
     HangingResistor &                                                  hanging_resistor, 
     CircuitMetadata &                                                  md,
     unordered_set<std::string> &                                       mn,
     std::map<std::string,FileSSFPair> &                                ssfm,
     std::map<std::string,IncludeFileInfo> &                            iflm,
     CircuitContext &                                                   cc,
     Topo::Topology &                                                   topology,
     Device::DeviceMgr &                                                device_manager,
     unordered_set<std::string> &                                       dNames,
     unordered_set<std::string> &                                       nNames,
     AliasNodeMap &                                                     alias_node_map,
     const std::vector< std::pair< std::string, std::string> > &        externalNetlistParams,
     Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group);

  CircuitBlock(
     const std::string &                                                fileName,
     const CmdParse &                                                   command_line,
     HangingResistor &                                                  hanging_resistor, 
     CircuitMetadata &                                                  md,
     unordered_set<std::string> &                                       mn,
     std::map<std::string,FileSSFPair> &                                ssfm,
     std::map<std::string,IncludeFileInfo> &                            iflm,
     CircuitContext &                                                   cc,
     CircuitBlock *                                                     mainCircPtr,
     CircuitBlock *                                                     parentCircPtr,
     Topo::Topology &                                                   topology,
     Device::DeviceMgr &                                                device_manager,
     unordered_set<std::string> &                                       dNames,
     unordered_set<std::string> &                                       nNames,
     AliasNodeMap &                                                     alias_node_map,
     const std::vector< std::pair< std::string, std::string> > &        externalNetlistParams,
     Teuchos::RCP<Xyce::Util::baseExpressionGroup> &                    group,
     std::vector<bool> &                                                pFilter,
     bool                                                               removeRedundant,
     bool                                                               modelBinning,
     double                                                             scale);

  ~CircuitBlock();

private:
  // Copy Constructor.
  CircuitBlock(CircuitBlock const& rhsCB);
  CircuitBlock &operator=(const CircuitBlock& rhsCB);

public:

  // Get the name of the circuit, if this is the root
  const std::string &getTitle() const
  {
    return title_;
  }

  // Get the name of the subcircuit, if this is a subcircuit  
  const std::string& getName() const
  {
    return name_;
  }

  // Set the name of the subcircuit, if this is a subcircuit
  void setName(std::string const& name)
  {
    name_ = name;
  }

  // Get the name of the netlist file
  const std::string& getNetlistFilename() const
  {
    return netlistFilename_;
  }

  // Get the path to the top-level netlist file, which may be either
  // absolute or relative to the execution subdirectory.
  const std::string& getTopLevelPath() const
  {
    return topLevelPath_;
  }

  AliasNodeMap &getAliasNodeMap()
  {
    return aliasNodeMap_;
  }

  const unordered_set< std::string > &getAliasNodeMapHelper() const 
  {
    return aliasNodeMapHelper_;
  }

  const std::vector<bool>& getPreprocessFilter() const
  {
    return preprocessFilter_;
  }

  const HangingResistor &getHangingResistor() const
  {
    return hangingResistor_;
  }

  HangingResistor &getHangingResistor()
  {
    return hangingResistor_;
  }

  CircuitContext* getCircuitContextPtr()
  {
    return &circuitContext_;
  }

  bool isSubcircuit() const
  {
    return parentCircuitPtr_ != NULL;
  }

  // This function parses the netlist file and fills in the
  // details of the circuit. This is phase 1 of netlist parsing.
  // The devices cannot be completely handled in this phase.
  bool parseNetlistFilePass1(PkgOptionsMgr &options_manager);
  bool parseNetlistFilePass1(PkgOptionsMgr &options_manager, const std::string &libSelect, 
                             std::vector< std::string >& libInside);

  // Perform special pass for mutual inductances
  bool parseMutualInductances();

  // Set data_->ssfPtr_ .
  void setSSFPtr( SpiceSeparatedFieldTool* ssfPtr )
  {
    ssfPtr_ = ssfPtr;
  }

  SpiceSeparatedFieldTool* getSSFPtr()
  {
    return ssfPtr_;
  }

  std::map<std::string,FileSSFPair>& getSSFMap()
  {
    return ssfMap_;
  }

  void setFilePosition(std::streampos const& position);
  void setLinePosition(int const& position);
  void setStartPosition();
  void setEndPosition();

  const std::streampos getStartPosition() const { 
    return fileStartPosition_; 
  }
  const std::streampos getEndPosition() const { 
    return fileEndPosition_; 
  }
  int getLineStartPosition() const { 
    return lineStartPosition_; 
  }
  int getLineEndPosition() const { 
    return lineEndPosition_; 
  }

  bool getSimpleSingleDevice() const {
    return simpleSingleDevice_;
  }
  const std::streampos getDevicePosition() const { 
    return devicePosition_; 
  }
  int getDeviceLine() const { 
    return deviceLine_;
  }

  // At the top level reset the file and line positions and skip the comment line.
  void resetSSFPtr();

  // Add a mutual inductor to the circuit.
  void addMutualInductor( DeviceBlock& device, CircuitContext* context );

  // Add a model to the circuit.
  void addModel(const ParameterBlock * model, std::string const& modelPrefix);

  // Add a device to the circuit.
  void addTableData(DeviceBlock & device );

  // Add a set of options corresponding to a .OPTIONS netlist line
  // to the circuit.
  void addOptions(const Util::OptionBlock &options);

  void addParams(const Util::OptionBlock &options);
  void addGlobalParams(const Util::OptionBlock &options);

  void registerGlobalParams(Util::UParamList & globalParams);

  // Search the subcircuitInstanceTable of the current circuit block for the
  // subcircuit of the given name. If it is not found, recursively
  // search each parent subcircuit. Return a pointer to the circuit
  // block if it is found, otherwise return NULL.
  CircuitBlock* findSubcircuit( std::string const& subcircuitName );

  // Change netlist file name
  void setFileName ( const std::string & fileNameIn );

  // resolve expressions in optionBlocks like .print
  bool resolveExpressionsInOptionBlocks();

  // update the aliasNodeMapHelper_ map to include "subcircuit
  // interface node names" that were embedded within expressions
  void updateAliasNodeMapHelper();

  // write out a copy of the netlist if it is requested.
  void writeOutNetlist();

  // Print the contents of CircuitBlock.
  void print();

  // Get information from handleAnalysis method.
  void setAnalysisName(const std::string analysisName)
  {
    analysisName_ = analysisName;
  }
  std::string& getAnalysisName() { return analysisName_; }
  bool getMORFlag() { return morFlag_; }

  const std::list<Util::OptionBlock>& getOptionsTable() const
  {
    return optionsTable_;
  }

  std::list<Util::OptionBlock>& getOptionsTable()
  {
    return optionsTable_;
  }

  CircuitMetadata& getMetadata()
  {
    return metadata_;
  }

  void setModelBinningFlag (bool modelBinning)
  {
    model_binning_flag_ = modelBinning;
  }

  bool getModelBinningFlag ()
  {
    return model_binning_flag_;
  }

  void setLengthScale (double scale)
  {
    lengthScale_ = scale;
  }

  double getLengthScale ()
  {
    return lengthScale_;
  }

  void setLevelSet( const std::set<int>& levelSet )
  {
    levelSet_ = levelSet;
  }

  const std::set<int>& getLevelSet() const
  {
    return levelSet_;
  }

  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & getExpressionGroup() { return expressionGroup_; }

public:

  // Lookup table for initial conditions
  std::map< std::string, TokenVector > initCondIndex;

private:
  std::string                           netlistFilename_;
  std::string                           topLevelPath_;     // path of the top-level netlist. May be absolute, or relative to execution dir
  std::string                           title_;                   // For top level circuit, given by first line of netlist
  std::string                           name_;                    // For subcircuits
  std::string                           analysisName_;

  unordered_set<std::string> &          nodeNames_;
  unordered_set<std::string> &          modelNames_;
  ModelMap                              modMap_;
  unordered_set<std::string> &          deviceNames_;
  std::list<Util::OptionBlock>          optionsTable_;
  std::set<int>                         levelSet_;                // Collect level numbers for later

  unordered_map<std::string, CircuitBlock *> circuitBlockTable_;

  // keep track of where include files are being called from.
  std::map<std::string, IncludeFileInfo> & includeFileLocation_; 
  
  // keep track of K lines that need extracted
  std::multimap< CircuitContext *, DeviceBlock > rawMIs_;

  const CmdParse &                      commandLine_;
  HangingResistor &                     hangingResistor_;
  
  // Circuit Context object to hold context information.
  CircuitContext & circuitContext_;

  CircuitMetadata & metadata_;

  std::vector< std::pair< std::string, std::string> > externalNetlistParams_;

  // Checking the netlist syntax.
  bool netlistSave_;

  // Check if MOR is requested.
  bool morFlag_;

  int devProcessedNumber_;

  // This is a map of node aliases.  Interface nodes to a subcircuit are removed as
  // the subcircuit is expanded into the netlist. We'll store the names of the interface
  // nodes in case the user accesses them elsewhere (as in a print statement)
  // The keys are the alias names and the values are the real circuit node names
  AliasNodeMap &                aliasNodeMap_;
  unordered_set< std::string >  aliasNodeMapHelper_;

  std::ifstream* netlistIn_;

  SpiceSeparatedFieldTool* ssfPtr_;

  std::streampos fileStartPosition_;
  std::streampos fileEndPosition_;
  int lineStartPosition_;
  int lineEndPosition_;

  // Data for single-device subcircuit parsing
  std::streampos devicePosition_;
  int            deviceLine_;
  int            numDevices_;
  bool           simpleSingleDevice_;

  void setDevicePosition()
  {
    devicePosition_ = ssfPtr_->getFilePosition();
    deviceLine_ = ssfPtr_->getLineNumber();
  }

  // Top level circuit pointer
  CircuitBlock* mainCircuitPtr_;

  // For subcircuits, points to the circuitBlock instance that contains
  // this subcircuit. NULL for top level.
  CircuitBlock* parentCircuitPtr_; 

  std::map<std::string,FileSSFPair> & ssfMap_;

  DeviceBlock device_;
  ParameterBlock tmpModel_;

  // This vector keeps track of any filters that are enabled to check for removal of 
  // "redundant" devices (where all device nodes are the same) or replace any occurrence 
  // of "GND", "GND!", or "GROUND" with "0" (so that all four terms are synonyms)
  std::vector<bool> preprocessFilter_;
  bool remove_any_redundant_;

  // model binning enable/disable and other parser settings
  bool model_binning_flag_;
  double lengthScale_;

  Topo::Topology &              topology_;
  Device::DeviceMgr  &          deviceManager_;

  Teuchos::RCP<Xyce::Util::baseExpressionGroup> expressionGroup_; ///< required for setting up expressions

  //This function preprocesses the netlist file to provide the user the
  //option of removing "redundant" devices (devices where all the nodes are
  //the same.  The info gathered here affects the phase 1 and phase 2 parse.
  bool parsePreprocess();

  //This function will reproduce a copy of the netlist file under the name
  //netlistfilename_copy.cir (to be used, in the end, to produce netlist files
  //which contain large resistors connecting "dangling" nodes to ground.
  void produceUnflattenedNetlist();

  // Handle a netlist line, determine the line type and take the
  // appropriate action.
  bool handleLinePass1( bool & result,
                        PkgOptionsMgr &options_manager, 
                        std::map<std::string,int> & fun, ModelMap & modMap,
                        const std::string & libSelect, 
                        std::vector< std::string > & libInside );

  // Handle a line for Mutual Inductance Pass
  bool getLinePassMI();

  // Post process analysis commands.
  bool handleAnalysis();

  // Parse the given include file adding the contents
  // to the current CircuitBlock.
  bool parseIncludeFile(
    PkgOptionsMgr &options_manager, 
    const std::string & includeFile, const std::string & libSelect, 
    std::map<std::string,int> & fun, ModelMap & modMap);

  // helper function for parseIncludeFile()
  void restorePrevssfInfo(
    SpiceSeparatedFieldTool* oldssfPtr,
    const std::string& old_netlistFilename,
    int oldFilePos,
    int oldLineNumber);

  // Retrieve separate IC= data from line or external file and temporarily
  // store in CircuitBlock
  void handleInitCond(TokenVector const& parsedLine );

  // Recurse through circuitBlockTable_ and get any .IC or .NODESET option blocks for
  // the usedSubcircuits.
  void getICNodesetList( std::vector<std::string>& usedSubcircuits, 
                         std::vector<Util::OptionBlock>& icNodesetList );
};

//-----------------------------------------------------------------------------
// Function      : packCircuitOptions
// Purpose       : send option blocks to all procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int packCircuitOptions(const std::list<Util::OptionBlock>& options, char* char_buffer,
                       int char_buffer_size, Xyce::Parallel::Communicator* pds_comm_ptr);

//-----------------------------------------------------------------------------
// Function      : unpackCircuitOptions
// Purpose       : unpack option blocks from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool unpackCircuitOptions(std::list<Util::OptionBlock>& options, char* char_buffer,
                          int bsize, Xyce::Parallel::Communicator* pds_comm_ptr);

//-----------------------------------------------------------------------------
// Function      : packAliasNodeMap
// Purpose       : send alias node map to all procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int packAliasNodeMap(const AliasNodeMap& alias_node_map, char* char_buffer,
                     int char_buffer_size, Xyce::Parallel::Communicator* pds_comm_ptr);

//-----------------------------------------------------------------------------
// Function      : unpackAliasNodeMap
// Purpose       : unpack alias node map from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool unpackAliasNodeMap(AliasNodeMap& alias_node_map, char* char_buffer,
                        int bsize, Xyce::Parallel::Communicator* pds_comm_ptr);


} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_CircuitBlock_h

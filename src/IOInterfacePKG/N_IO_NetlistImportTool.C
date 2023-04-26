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
// Filename      : NetlistImportTool.C
//
// Purpose       : Implement the interface to to read and parse a netlist for
//                 an electrical circuit.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 07/28/00
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <N_UTL_fwd.h>

#include <N_DEV_DeviceMgr.h>
#include <N_DEV_RegisterDevices.h>
#include <N_DEV_RegisterOpenDevices.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_DistToolFactory.h>
#include <N_IO_MORAnalysisTool.h>
#include <N_IO_NetlistImportTool.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_ParsingMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_WildcardSupport.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Marshal.h>
#include <N_UTL_Expression.h>

#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

#include <expressionGroup.h>

namespace Xyce {

namespace Util {

Marshal &operator<<(Marshal &mout, const IO::UndefinedName &undefined_name)
{
  return mout << undefined_name.getName() << undefined_name.getNetlistLocation();
}

Marshal &operator>>(Marshal &min, IO::UndefinedName &undefined_name)
{
  std::string name;
  NetlistLocation netlist_location;
  min >> name >> netlist_location;
  undefined_name.setName(name);
  undefined_name.setNetlistLocation(netlist_location);
  return min;
}

}


namespace IO {

namespace {

//-----------------------------------------------------------------------------
// Class         : STEPOptionsReg
// Purpose       : functor for registering STEP options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct STEPOptionsReg : public PkgOptionsReg
{
  STEPOptionsReg( NetlistImportTool &netlist_import_tool )
    : netlistImportTool_(netlist_import_tool)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return netlistImportTool_.registerSTEPOptions(options.begin(), options.end());
  }

  NetlistImportTool &netlistImportTool_;
};

//-----------------------------------------------------------------------------
// Class         : SAMPLINGOptionsReg
// Purpose       : functor for registering SAMPLING options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : Eric Keiter, SNL
// Creation Date : 11/1/2017
//-----------------------------------------------------------------------------
struct SAMPLINGOptionsReg : public PkgOptionsReg
{
  SAMPLINGOptionsReg( NetlistImportTool &netlist_import_tool )
    : netlistImportTool_(netlist_import_tool)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return netlistImportTool_.registerSAMPLINGOptions(options.begin(), options.end());
  }

  NetlistImportTool &netlistImportTool_;
};

//-----------------------------------------------------------------------------
// Class         : EMBEDDEDSAMPLINGOptionsReg
// Purpose       : functor for registering EMBEDDEDSAMPLING options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : Eric Keiter, SNL
// Creation Date : 6/3/2019
//-----------------------------------------------------------------------------
struct EMBEDDEDSAMPLINGOptionsReg : public PkgOptionsReg
{
  EMBEDDEDSAMPLINGOptionsReg( NetlistImportTool &netlist_import_tool )
    : netlistImportTool_(netlist_import_tool)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return netlistImportTool_.registerEMBEDDEDSAMPLINGOptions(options.begin(), options.end());
  }

  NetlistImportTool &netlistImportTool_;
};

//-----------------------------------------------------------------------------
// Class         : PCEOptionsReg
// Purpose       : functor for registering PCE options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : Eric Keiter, SNL
// Creation Date : 6/3/2019
//-----------------------------------------------------------------------------
struct PCEOptionsReg : public PkgOptionsReg
{
  PCEOptionsReg( NetlistImportTool &netlist_import_tool )
    : netlistImportTool_(netlist_import_tool)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return netlistImportTool_.registerPCEOptions(options.begin(), options.end());
  }

  NetlistImportTool &netlistImportTool_;
};

//-----------------------------------------------------------------------------
// Class         : DCOptionsReg
// Purpose       : functor for registering DC options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct DCOptionsReg : public PkgOptionsReg
{
  DCOptionsReg( NetlistImportTool &netlist_import_tool )
    : netlistImportTool_(netlist_import_tool)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return netlistImportTool_.registerDCOptions(options.begin(), options.end());
  }

  NetlistImportTool &netlistImportTool_;
};

//-----------------------------------------------------------------------------
// Class         : DISTOptionsReg
// Purpose       : functor for registering Dist options
// Special Notes : Used by package manager addOptionsProcessor method
// Creator       : hkthorn
// Creation Date : 01/11/2018
//-----------------------------------------------------------------------------
struct DISTOptionsReg : public PkgOptionsReg
{
  DISTOptionsReg( NetlistImportTool &netlist_import_tool )
    : netlistImportTool_(netlist_import_tool)
  {}

  bool operator()( const Util::OptionBlock & options )
  {
    return netlistImportTool_.setDISTOptions( options );
  }

  NetlistImportTool &netlistImportTool_;
};

//-----------------------------------------------------------------------------
// Function      : errorUndefinedParameters
// Purpose       : used to print out error message if there are undefined
//                 parameters on the .PRINT line
// Special Notes :
// Creator       : dgbaur
// Creation Date : 03/23/2015
//-----------------------------------------------------------------------------
void errorUndefinedParameters(
  const UndefinedNameSet &      undefined_parameters)
{
  std::map<NetlistLocation, std::string> msg_map;
  int count = 0;
  for (UndefinedNameSet::const_iterator it = undefined_parameters.begin(), end = undefined_parameters.end(); it != end; ++it) {
    std::string &msg = msg_map[(*it).getNetlistLocation()];
    if (!msg.empty())
      msg += ", ";
    msg += (*it).getName();
    ++count;
  }

  for (std::map<NetlistLocation, std::string>::const_iterator it = msg_map.begin(), end = msg_map.end(); it != end; ++it) {
    Report::UserError0().at((*it).first)
      << "There " << (count == 1 ? "was " : "were ") << count << " undefined symbol"
      << (count == 1 ? "" : "s") << " in .PRINT command: " << (*it).second;
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(NetlistImportTool &netlist_import_tool, PkgOptionsMgr &options_manager)
{
  NetlistImportTool::populateMetadata( options_manager );

  options_manager.addCommandProcessor("DC", new DCOptionsReg(netlist_import_tool));
  options_manager.addCommandProcessor("STEP", new STEPOptionsReg(netlist_import_tool));
  options_manager.addCommandProcessor("DIST", new DISTOptionsReg(netlist_import_tool));
  options_manager.addCommandProcessor("SAMPLING", new SAMPLINGOptionsReg(netlist_import_tool));
  options_manager.addCommandProcessor("EMBEDDEDSAMPLING", new EMBEDDEDSAMPLINGOptionsReg(netlist_import_tool));
  options_manager.addCommandProcessor("PCE", new PCEOptionsReg(netlist_import_tool));

  return true;
}


//-------------------------------------------------------------------------
// Name          : NetlistImportTool::NetlistImportTool
// Purpose       : Default constructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::NetlistImportTool(
  Util::Op::BuilderManager &    op_builder_manager,
  const ParsingMgr &            parsing_manager,
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group)
  : parsingMgr_(parsing_manager),
    mainCircuitBlock_(NULL),
    distributionTool_(NULL),
    currentContextPtr_(NULL),
    metadata_(),
    circuitContext_(op_builder_manager, parsing_manager, contextList_, currentContextPtr_),
    useMOR_(false),
    expressionGroup_(group)
{}

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::~NetlistImportTool
// Purpose       : Destructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
NetlistImportTool::~NetlistImportTool()
{
  // Destroy the SpiceSeparatedFieldTool map, close the input files.
  for (std::map<std::string, FileSSFPair>::iterator it = ssfMap_.begin(), end = ssfMap_.end(); it != end; ++it)
  {
    delete (*it).second.first;
    delete (*it).second.second;
  }

  delete mainCircuitBlock_;
  delete distributionTool_;
}

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::populateMetadata
// Purpose       : 
// Special Notes :
// Creator       : 
// Creation Date : 
//-------------------------------------------------------------------------
void NetlistImportTool::populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("DIST");

  parameters.insert(Util::ParamMap::value_type("STRATEGY", Util::Param("STRATEGY", 0)));
}

//-------------------------------------------------------------------------
// Function      : NetlistImportTool::constructCircuitFromNetlist
//
// Purpose       : Construct a circuit from a netlist.  
//
// Special Notes : This function contains two passes of the netlist.
//
//                 The first is entirely on proc0, and gathers device 
//                 counts and global information.  The global info is 
//                 broadcast to all procs.
//
//                 The second pass distributes devices in parallel 
//                 and expands subcircuits.  How they are distributed 
//                 depends on which distribution manager is selected.
//
// Creator       : Lon Waters
// Creation Date : 07/28/2000
//-------------------------------------------------------------------------
int NetlistImportTool::constructCircuitFromNetlist(
  const CmdParse &                                              command_line,
  HangingResistor &                                             hanging_resistor,
  const std::string &                                           netlist_filename,
  const std::vector< std::pair< std::string, std::string> > &   externalNetlistParams,
  Topo::Topology &                                              topology,
  Parallel::Communicator &                                      pds_comm,
  PkgOptionsMgr &                                               options_manager,
  OutputMgr &                                                   output_manager,
  Device::DeviceMgr &                                           device_manager,
  Measure::Manager &                                            measure_manager,
  FourierMgr &                                                  fourier_manager,
  FFTMgr &                                                      fft_manager,
  unordered_set<std::string> &device_names )
{
  Parallel::Machine comm = pds_comm.comm();

  // Create a circuitBlock instance to hold netlist circuit data.
  mainCircuitBlock_ = new IO::CircuitBlock(
    netlist_filename,
    command_line,
    hanging_resistor,
    metadata_,
    modelNames_ ,
    ssfMap_,
    iflMap_,
    circuitContext_,
    topology,
    device_manager,
    device_names,
    nodeNames_,
    aliasNodeMap_,
    externalNetlistParams,
    expressionGroup_
    );

  // Register the mutual inductors here, they are handled differently than every other 
  // device right now.
  Device::registerMutualInductors();

  // Parse the netlist file.
  if (DEBUG_IO)
    Xyce::dout() << "Starting netlist parsing." << std::endl;

  {
    Stats::StatTop _parseStat("Parse Context");
    Stats::TimeBlock _parseTimer(_parseStat);

    if (Parallel::rank(comm) == 0)
    {
      // Perform initial pass through the circuit file to generate hierarchical context object.
      // Read in the circuit context and circuit options on the root processor.
      mainCircuitBlock_->parseNetlistFilePass1(options_manager); 
#if  0
      // ERK.  Comments to possibly remove later, related to issue 24.
      //
      // My work on issue 24 has to do with resolving parameters.  These are resolved in each circuit 
      // context via the CircuitContext::resolve function call.
      //
      // For the top level context, CircuitContext::resolve is called on proc 0 during the first pass.
      // The top level contains lots of global information that will be needed on all processors, so a
      // global broadcast is performed after this.
      //
      // One of the things included in that global broadcast is the list of global parameters, which until 
      // recently were only allowed at the top level of the netlist.  This list of global parameters will 
      // later be passed over to the device manager.  non-global parameters are not saved at all.
#endif
    }

    // Just in case an error was reported in parsing the netlist.
    N_ERH_ErrorMgr::safeBarrier(comm);
  }
 
  // If we are only performing a device count, we can exit here.
  if (command_line.argExists("-count"))
  {
    return 1;
  }

  // If MOR is requested, analyze the circuit before performing full simulation
  if (mainCircuitBlock_->getMORFlag() && mainCircuitBlock_->getAnalysisName() != "MOR")
  {
    IO::MORAnalysisTool morAnalyzer(circuitContext_, mainCircuitBlock_->getOptionsTable(),
                                    command_line, ssfMap_, iflMap_);

     useMOR_ = morAnalyzer.reduceLinearSubcircuits();

     // If model-order reduction was used, there is no need to run the original circuit.
     if (useMOR_) {
       return 1;
     }
     else {
       // Remove .MOR and .OPTIONS MOR_OPTS from option table before continuing simulation
       // MOR was determined to not be viable.
       morAnalyzer.removeMOROptions();
     }
  }

  // register distribution options before creating distribution tool
  registerDistOptions(options_manager, mainCircuitBlock_->getOptionsTable());

  distributionTool_ = DistToolFactory::create(&pds_comm, distOptions_, 
                                              *mainCircuitBlock_, ssfMap_, 
                                              iflMap_, externalNetlistParams, parsingMgr_);

  // Distribute the circuit context and circuit options to all other processors.
  distributionTool_->broadcastGlobalData();
#if 0
  // ERK.  Comments to possibly remove later, related to issue 24.
  // From here (the distributionTool_->broadcastGlobalData() call)
  // all the global parameters (among other things) from the top level context are broadcast to other processors.
  //
  // Note, the total device *count* is determined during the first pass, I think.   If this is the case, that 
  // means that a skeleton of the subcircuit heirarchy must be traversed.  The count is needed for the distribution.
#endif

  // Now we know what devices and models are used in the netlist, so create the device configuration.
  // NOTE:  This configuration is created without the mutual inductors and inductor, hence the 'false' argument.
  Device::registerDevices(circuitContext_.getDeviceCountMap(), mainCircuitBlock_->getLevelSet(), false);

  // register the global parameters
  registerGlobalParams(device_manager, circuitContext_.getGlobals().begin(), circuitContext_.getGlobals().end());

  // register circuit options
  // NOTE:  This method does not process distribution options, since they were processed above.
  registerCircuitOptions(options_manager, mainCircuitBlock_->getOptionsTable());

  mainCircuitBlock_->setModelBinningFlag( parsingMgr_.getModelBinningFlag() );
  mainCircuitBlock_->setLengthScale( parsingMgr_.getLengthScale() );


  // Perform second pass
  {
    Stats::StatTop _distributeStat("Distribute Devices");
    Stats::TimeBlock _distributeTimer(_distributeStat);

    // Distribute and instantiate the circuit devices.
    distributionTool_->distributeDevices();   

#if 0
    // ERK.  Comments to possibly remove later, related to issue 24.
    //
    // all the subcircuit "CircuitContext::resolve" calls will happen here, under the distributeDevices call. 
    // Each subckt instance is handled like a device and expanded as those lines are encountered.  
    // They are handled by the function DistToolDefault::expandSubcircuitInstance, which subsequently calls CircuitContext::resolve.
    //
    // This means that for the top-level context, "CircuitContext::resolve" is called in the first pass on proc0, and can 
    // be broadcast to all procs, via the "broadcastGlobalData" function call. 
    //
    // For all subordinate contexts, "CircuitContext::resolve" is called in the second pass, but also on proc 0, like before.
    // So, the .param and .global_param containers are complete at this point, and on proc zero.  
    // But .global_params aren't broadcast in parallel at this stage.  They need to be.
    //
    // Other stuff (not .global_params) from this second pass gets sent around in parallel, but the focus for that parallel communication is on sending device metadata.
    //
    // So, any .params that needed to be converted to .global_params down in subordinate subcircuits during ::resolve
    // will need to be broadcast in parallel, much like above.  
    //
    // All processors ultimately need a single master list of global parameters in the device package.  All non-global params get thrown away once parsing is done.
    //
    // So, this is one practical reason that global_params have previously only been allowed in the top level netlist.  
    // The machinery for handling subcircuit global_params doesn't exist.
#endif
  }

#if 1
  // register additional global parameters, which have been resolved in subcircuits
  Util::UParamList & newGlobals = distributionTool_->getAdditionalGlobalParams();
  registerGlobalParams(device_manager, newGlobals.begin(), newGlobals.end());
#endif

  // register additional circuit options after distribution
  registerCircuitOptions(options_manager, distributionTool_->getAdditionalOptions());

  // Check for name collisions between devices
  checkNodeDevConflicts(device_names, pds_comm);

  if (DEBUG_IO)
    Xyce::dout() << "Completed netlist parsing. ";

  // Write out preprocessed netlist, if requested.
  if (Parallel::rank(comm) == 0)
  {
    mainCircuitBlock_->writeOutNetlist();
  }

  output_manager.setTitle(mainCircuitBlock_->getTitle());

  // printLineDiagnostics scans all output parameters requested by print lines
  // and flags any errors or inconsistencies.
  printLineDiagnostics(comm, output_manager.getOutputParameterMap(),
                       device_manager, measure_manager, nodeNames_,
                       stepParams_,
                       samplingParams_,
                       embeddedSamplingParams_,
                       pceParams_,
                       dcParams_, device_names, aliasNodeMap_,
                       deferredUndefinedParameters_);

  return 1;
}

namespace {

bool findNode(const std::string &name, const unordered_set<std::string> &node_set, const AliasNodeMap &alias_map)
{
  bool ret = true;
  if (!node_set.count( name )) {
    AliasNodeMap::const_iterator alias_node_it = alias_map.find(name);
    if (alias_node_it != alias_map.end())
      ret = (node_set.count((*alias_node_it).second) > 0);
    else
      ret = false;
  }

  return ret;
}

} // end of unnamed

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDCOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/21/06
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerDCOptions(
  Util::ParamList::const_iterator       it,
  Util::ParamList::const_iterator       end)
{
  for (; it != end; ++it)
  {
    if (it->tag() == "PARAM")
    {
      dcParams_.push_back((*it).stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSTEPOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/21/06
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerSTEPOptions(
  Util::ParamList::const_iterator       it,
  Util::ParamList::const_iterator       end)
{
  for (; it != end; ++it)
  {
    if ((*it).tag() == "PARAM")
    {
      stepParams_.push_back((*it).stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setDISTOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 01/11/2018
//-----------------------------------------------------------------------------
bool NetlistImportTool::setDISTOptions( const Util::OptionBlock& distOptions )
{
  distOptions_ = distOptions; 

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSAMPLINGOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/1/2017
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerSAMPLINGOptions(
  Util::ParamList::const_iterator       it,
  Util::ParamList::const_iterator       end)
{
  for (; it != end; ++it)
  {
    if ((*it).tag() == "PARAM")
    {
      samplingParams_.push_back((*it).stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerEMBEDDEDSAMPLINGOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/3/2019
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerEMBEDDEDSAMPLINGOptions(
  Util::ParamList::const_iterator       it,
  Util::ParamList::const_iterator       end)
{
  for (; it != end; ++it)
  {
    if ((*it).tag() == "PARAM")
    {
      embeddedSamplingParams_.push_back((*it).stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerPCEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 6/3/2019
//-----------------------------------------------------------------------------
bool NetlistImportTool::registerPCEOptions(
  Util::ParamList::const_iterator       it,
  Util::ParamList::const_iterator       end)
{
  for (; it != end; ++it)
  {
    if ((*it).tag() == "PARAM")
    {
      pceParams_.push_back((*it).stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : printLineDiagnostics
// Purpose       :
// Special Notes : erkeite:  8/10/2007
//
// For commentary and some explanation of this function, see the comments
// for the deferredParameterDiagnostics function.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/18/06
//-----------------------------------------------------------------------------
void printLineDiagnostics(
  Parallel::Machine                             comm,
  const OutputParameterMap &                    output_parameter_map,
  const Device::DeviceMgr &                     device_manager,
  const Measure::Manager &                      measure_manager,
  const unordered_set<std::string> &            node_names,
  const std::vector<std::string> &              step_params,
  const std::vector<std::string> &              sampling_params,
  const std::vector<std::string> &              embeddedSampling_params,
  const std::vector<std::string> &              pce_params,
  const std::vector<std::string> &              dc_params,
  const unordered_set<std::string> &    device_map,
  const IO::AliasNodeMap &                      alias_node_map,
  UndefinedNameSet &                            deferred_undefined_parameters)
{
  UndefinedNameSet undefined_parameters;

  for (OutputParameterMap::const_iterator it1 = output_parameter_map.begin(), end1 = output_parameter_map.end(); it1 != end1; ++it1)
  {
    const OutputParameterMap::mapped_type &parameter_vector = (*it1).second;

    for (OutputParameterMap::mapped_type::const_iterator it2 = parameter_vector.begin(), end2 = parameter_vector.end(); it2 != end2; ++it2)
    {
      const PrintParameters &print_parameters = (*it2);

      // Populate stringStat, nodeStat and instanceState with variable names from print line
      std::map<std::string, bool> stringStat;
      std::map<std::string, bool> nodeStat;
      std::map<std::string, bool> instanceStat;

      for (Util::ParamList::const_iterator it3 = print_parameters.variableList_.begin(), end3 = print_parameters.variableList_.end(); it3 != end3; ++it3)
      {
        const Util::Param &parameter = (*it3);

        std::vector<std::string> nodes;
        std::vector<std::string> instances;
        std::vector<std::string> leads;
        std::vector<std::string> strings;
        std::vector<std::string> special;
        if (Util::hasExpressionTag(parameter) )
        {
          // parameter starts with "{" but may not have been been parsed into an expression.
          // check if there is an underlying expression object with the parameter
          if (parameter.getType() == Util::EXPR)
          {
            parameter.getValue<Util::Expression>().getVoltageNodes(nodes);
            parameter.getValue<Util::Expression>().getDeviceCurrents(instances);
            parameter.getValue<Util::Expression>().getLeadCurrents(leads);
            parameter.getValue<Util::Expression>().getUnresolvedParams(strings);
            parameter.getValue<Util::Expression>().getSpecials(special);
          }
          else if ( ((parameter.getType()) == Util::DBLE) || 
                    ((parameter.getType()) == Util::INT) || 
                    ((parameter.getType()) == Util::CMPLX) )
          {
          }
          instances.insert(instances.end(), leads.begin(), leads.end());  // and lead currents to instances
          strings.insert(strings.end(), special.begin(), special.end());  // add specials to strings
        }
        else
        {
          const std::string &varType = parameter.tag();

          // Note: for I and V, the check for (varType.size() == 2 || varType.size() == 3) && varType[0] == ...)
          // is needed to support operators like IR, II, IM, IP and IDB.  
          if ((varType == "I" || ((varType.size() == 2 || varType.size() == 3) && varType[0] == 'I')) && parameter.getImmutableValue<int>() > 0)
          {
            // any devices found in this I(xxx) structure need to be communicated to the device manager
            // so that the lead currents can be calculated
            if (parameter.getImmutableValue<int>() != 1)
            {
              Report::UserError() << "Only one device argument allowed in I() in .PRINT command";
            }
            else
            {
              ++it3;
              if (varType.size() == 2)
              {
                instances.push_back((*it3).tag() + "{" + varType[1] + "}");
              }
              else
              {
                instances.push_back((*it3).tag());
              }
            }
          }
          else if ((varType == "V" || ((varType.size() == 2 || varType.size() == 3) && varType[0] == 'V')) && parameter.getImmutableValue<int>() > 0)
          {
            int numIndices = parameter.getImmutableValue<int>();

            if (numIndices < 1 || numIndices > 2)
            {
              Report::UserError0() << "Only one or two node arguments allowed in V() in .PRINT command";
            }
            else
            {
              for (; numIndices > 0 ; --numIndices)
              {
                ++it3;
                nodes.push_back((*it3).tag());
              }
            }
          }
          else if (varType == "N" && parameter.getImmutableValue<int>() > 0)
          {
            // Don't attempt to check this type of variable.
            //
            // N-variables can only be fully resolved once the devices are allocated,
            // and variables which are internal to the devices have been assigned names.
            //(N-variables are the same as V-variables, except that they can include
            // internal vars).
            if (parameter.getImmutableValue<int>() != 1)
              Report::UserError() << "Only one device argument allowed in N() in .PRINT command";
            else
              ++it3;
          }
          else if (varType.size() == 3 && varType[0] == 'D' && parameter.getImmutableValue<int>() > 0)
          {
            int numIndices = parameter.getImmutableValue<int>();
            if (numIndices < 1 || numIndices > 2)
            {
              Report::UserError() << "Only one or two arguments allowed in DNI() or DNO() in .PRINT command";
            }
            else 
            {
              // DNO or DNI operators come in two forms, DNO(deviceName) or DNO(deviceName,noiseType)
              // We only want to push the deviceName onto instances, and not the noiseType
              ++it3;
              instances.push_back((*it3).tag());
              if (numIndices == 2) ++it3;
            }

          }
          else if ( ((varType[0] == 'S') || (varType[0] == 'Y') || (varType[0] == 'Z')) &&
                    varType.size() <= 3 && parameter.getImmutableValue<int>() > 0)
          {
            // Don't trigger this clause if varType == "SENS". Only for S(), SR(), SI(),
            // SP(), SM() and SDB(), and the corresponding Y and Z operators.
            int numIndices = parameter.getImmutableValue<int>();
            if (numIndices != 2)
            {
              Report::UserError() << "S(), Y() and Z() must have two arguments in .PRINT command";
            }
            else
            {
	      ++it3; ++it3;
            }
          }
          else if (((varType == "P") || (varType == "W")) && parameter.getImmutableValue<int>() > 0)
          {
            // any devices found in this P(xxx) or W(xxx) structure need to be communicated to the device manager
            // so that the lead currents and power can be calculated
            if (parameter.getImmutableValue<int>() != 1)
            {
              Report::UserError() << "Only one device argument allowed in P() or W() in .PRINT command";
            }
            else
            {
              ++it3;
              instances.push_back((*it3).tag());
            }
          }
          else
          {
            strings.push_back(parameter.tag());
          }
        }

        for (std::vector<std::string>::const_iterator it4 = strings.begin(), end4 = strings.end(); it4 != end4; ++it4)
        {
          const std::string &name = *it4;

          bool done = false;
          if (name == "TEMP" || name == "TEMPER" || name == "TIME" || name == "FREQ" || name == "HERTZ" || name == "INDEX" ||
              name == "STEPNUM" || name == "OBJFUNC" || name == "OBJVARS" || name == "SENS" || name == "NOISE" || name == "sweep" ||
              name == "ONOISE" || name == "INOISE")
          {
            done = true;
          }

          if (!done)
          {
            // check if this is a dc sweep value.
            for (std::vector<std::string>::const_iterator it = dc_params.begin(), end = dc_params.end(); it != end ; ++it)
            {
              if (*it == name)
              {
                done = true;
                break;
              }
            }
          }

          if (!done)
          {
            // check if this is a step sweep value.
            for (std::vector<std::string>::const_iterator it = step_params.begin(), end = step_params.end(); it != end ; ++it)
            {
              if (*it == name)
              {
                done = true;
                break;
              }
            }
          }

          if (!done)
          {
            // check if this is a sampling value.
            for (std::vector<std::string>::const_iterator it = sampling_params.begin(), end = sampling_params.end(); it != end ; ++it)
            {
              if (*it == name)
              {
                done = true;
                break;
              }
            }

            // check if this is an embedded sampling value.
            for (std::vector<std::string>::const_iterator it = embeddedSampling_params.begin(), end = embeddedSampling_params.end(); it != end ; ++it)
            {
              if (*it == name)
              {
                done = true;
                break;
              }
            }

            // check if this is a PCE value.
            for (std::vector<std::string>::const_iterator it = pce_params.begin(), end = pce_params.end(); it != end ; ++it)
            {
              if (*it == name)
              {
                done = true;
                break;
              }
            }

          }

          double mValue;
          bool mFound = measure_manager.getMeasureValue(name, mValue);
          if (mFound)
            done = true;

          if (!done)
          {
            done = device_manager.parameterExists(comm, name);

            // Note:  when this is called, the devices haven't been allocated
            // yet.  So, this particular check isn't reliable.  If the device and/or
            // parameter is not found by the device package via getParam, then
            // simply save it for a later diagnostic.
            if (!done)
            {
              deferred_undefined_parameters.insert(UndefinedName(name, print_parameters.netlistLocation_));
            }
            done = true;
          }
          stringStat[name] = done;
        }

        for (std::vector<std::string>::const_iterator iter_s = nodes.begin(); iter_s != nodes.end(); ++iter_s)
        {
          const std::string &name= *iter_s;

          // Allow * to pass.  This also assumes that something like V(1*) is a wildcard request
          // for all nodal voltages that start with '1', rather than the voltage at node 1*.
          nodeStat[name] = (name == "*") || findWildCardMatch(name,node_names) ||
	                   (findNode(name, node_names, alias_node_map));
        }

        for (std::vector<std::string>::iterator iter_s= instances.begin(); iter_s != instances.end(); ++iter_s)
        {
          const std::string &name = *iter_s;

          if (name.substr(name.size() - 1, 1) == "}")
          {
            // Allow * to pass.  This also assumes that something like I(R1*) is a wildcard
            // request for all R devices nodes that start with "R1".
            instanceStat[name] = (name[0] == '*') || findWildCardMatch(name.substr(0, name.size() - 3),device_map) ||
                                 (device_map.find(name.substr(0, name.size() - 3)) != device_map.end());
          }
          else
          {
            // Allow * to pass.  This also assumes that something like I(R1*) is a wildcard
            // request for all R devices nodes that start with "R1".
            instanceStat[name] = (name == "*") || findWildCardMatch(name,device_map) ||
                                 (device_map.find(name) != device_map.end());
          }
        }
      }

      // Sum uses of names in nodeStat and instanceStat
      std::vector<int> stat;
      for (std::map<std::string, bool>::const_iterator it = nodeStat.begin(), end = nodeStat.end(); it != end; ++it)
      {
        stat.push_back((*it).second ? 1 : 0);
      }

      for (std::map<std::string, bool>::const_iterator it = instanceStat.begin(), end = instanceStat.end(); it != end; ++it)
      {
        stat.push_back((*it).second ? 1 : 0);
      }

      Parallel::AllReduce(comm, MPI_MAX, stat);

// Generate message
      std::ostringstream oss;

      int stat_index = 0;
      for (std::map<std::string, bool>::const_iterator it = nodeStat.begin(), end = nodeStat.end(); it != end; ++it)
      {
        if (stat[stat_index++] == 0)
        {
          std::ostringstream oss;
          oss << "node " << (*it).first;
          undefined_parameters.insert(UndefinedName(oss.str(), print_parameters.netlistLocation_));
        }
      }

      for (std::map<std::string, bool>::const_iterator it = instanceStat.begin(), end = instanceStat.end(); it != end; ++it)
      {
        if (stat[stat_index++] == 0)
        {
          // Insert into undefined_parameters, except for devices (L,K,N,U and Y) that may not
          // support lead currents or power.  For those device types, error checking occurs during 
          // operator creation.  For the other devices, that support lead currents and power, we will
          // only get here if the node name is invalid.
          if ( !deferErrorCheckUntilOpCreation((*it).first) )
	  {
            std::ostringstream oss;
            oss << "device " << (*it).first;
            undefined_parameters.insert(UndefinedName(oss.str(), print_parameters.netlistLocation_));
          }
        }
      }

      for (std::map<std::string, bool>::const_iterator it= stringStat.begin(); it != stringStat.end(); ++it)
      {
        if (!(*it).second)
        {
          std::ostringstream oss;
          oss << (*it).first;
          undefined_parameters.insert(UndefinedName(oss.str(), print_parameters.netlistLocation_));
        }
      }
    }

  }

  if (!undefined_parameters.empty())
    errorUndefinedParameters(undefined_parameters);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr:deferErrorCheckUntilOpCreation
// Purpose       : Some devices do not support lead currents and/or power
//                 calculations.  In addition, inductors (L) can be part of a 
//                 mutual inductor (K).  For those devices, error checking (for
//                 either non-existent devices or for "non-support" of the I, 
//                 P or W operators) is deferred until Operator creation.  This 
//                 is its own function because it is used in multiple places.
//                 It's also useful to explain what it does.
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 02/14/17
//-----------------------------------------------------------------------------
bool deferErrorCheckUntilOpCreation(std::string devStr)
{
  bool bsuccess = ((startswith_nocase(devStr.c_str(), "L")) || (devStr.find( ":L") != std::string::npos ) || 
                (devStr.find( ":l") != std::string::npos ) || 
                (startswith_nocase(devStr.c_str(), "K")) || (devStr.find( ":K") != std::string::npos ) || 
                (devStr.find( ":k") != std::string::npos ) || 
                (startswith_nocase(devStr.c_str(), "O")) || (devStr.find( ":O") != std::string::npos ) || 
                (devStr.find( ":o") != std::string::npos ) ||
                (startswith_nocase(devStr.c_str(), "U")) || (devStr.find( ":U") != std::string::npos ) || 
                (devStr.find( ":u") != std::string::npos ) ||
                (startswith_nocase(devStr.c_str(), "Y")) || (devStr.find( ":Y") != std::string::npos ) || 
		(devStr.find( ":y") != std::string::npos ));

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::deferredParameterDiagnostics
// Purpose       : Yet another function that checks the .print line.
//
// Special Notes : erkeite.  8/10/2007
//
//                 This function has to be called after the devices in
//                 the device package have all been allocated.  The function
//                 printLineDiagnostics is called earlier than that, so
//                 it doesn't handle device-specific parameter outputs.
//
//                 We have 3 completely different functions that do different
//                 aspects of checking the print line:
//
//                (1) printLineDiagnostics:  happens really early, as soon as
//                 topology is constructed, but before devivces are allocated.
//
//                (2) delayedPrintLineDiagnostics:  happens just a little later, to
//                 check .print line fields that are specific to device entities(ie
//                 aren't just I and V solution vars).
//
//                (3) check_output:  Happens at the beginning of the simulation,
//                 after everything is completely set up, but before any solves
//                 happen.  The advantage of this one is that it actually uses
//                 all the same function calls as the regular output call that
//                 will happen later.  So, it is probably the most thorough test
//                 of all.  However, for really big simulations, it happens relatively
//                 late in the setup.  For really large simulations, setup can take
//                 a long time, and it is often useful to get diagnostics earlier
//                 than that.
//
//                 Note from TVR 10/22/2013:  "check_output", referred to
//                 above, has been broken apart and renamed "prepareOutput".
//                 The notes from ERK above are therefore somewhat out of date.
//                 The refactor, though, is not complete, and the replacement
//                 code is not well documented.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
void
deferredParameterDiagnostics(
  Parallel::Machine                             comm,
  const UndefinedNameSet &                      deferred_undefined_parameters,
  const Device::DeviceMgr &                     device_manager)
{
  std::ostringstream oss;

  Util::Marshal mout;
  mout << deferred_undefined_parameters;

  std::vector<std::string> dest;
  Parallel::AllGatherV(comm, mout.str(), dest);

  UndefinedNameSet common_parameter_set;
  for (int p = 0; p < Parallel::size(comm); ++p)
  {
    Util::Marshal min(dest[p]);
    std::vector<UndefinedName> x;
    min >> x;
    std::copy(x.begin(), x.end(), std::inserter(common_parameter_set, common_parameter_set.begin()));
  }

  UndefinedNameSet undefined_parameters;
  for (std::set<UndefinedName>::const_iterator it = common_parameter_set.begin(), end = common_parameter_set.end(); it != end; ++it)
  {
    // double result = 0.0;
    // if (!getParamAndReduce(comm, device_manager, (*it).getName(), result))
    if (!device_manager.parameterExists(comm, (*it).getName()))
      undefined_parameters.insert(*it);
  }

  if (!undefined_parameters.empty())
    errorUndefinedParameters(undefined_parameters);
}

//-----------------------------------------------------------------------------
// Function      : checkNodeDevConflicts
// Purpose       : Check for name collisions between devices
// Special Notes : Non-member helper function
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 07/28/14
//-----------------------------------------------------------------------------
void
checkNodeDevConflicts(
  const unordered_set<std::string> &     device_names,
  Parallel::Communicator &               pds_comm)
{
  Parallel::Machine comm = pds_comm.comm();

  int proc_size = Parallel::size(comm);
  int proc_rank = Parallel::rank(comm);

  if (proc_size > 1)
  {
    int byteCount = 0;

    // Count the number of bytes needed to send this buffer.
    if (proc_rank > 0)
    {
      for (unordered_set<std::string>::const_iterator it = device_names.begin(), end = device_names.end(); it != end; ++it)
      {
        byteCount += sizeof(int) + (*it).size();
      }
    }

    N_ERH_ErrorMgr::safeBarrier(comm);     // All procs call (5)

    int maxByteCount = 0;
    Parallel::AllReduce(comm, MPI_MAX, &byteCount, &maxByteCount, 1);

    std::vector<char> sendBuffer, recvBuffer;
    unordered_set<std::string> tmp_device_names;

    if (proc_rank == 0)
    {
      recvBuffer.resize(maxByteCount);
      tmp_device_names = device_names;
    }

    for (int proc = 1; proc < proc_size; ++proc)
    {
      int bsize = 0;
      if (proc_rank == 0)
      {
        // Get the size of the buffer and then get the buffer.
        pds_comm.recv( &bsize, 1, proc );
        pds_comm.recv( &recvBuffer[0], bsize, proc ); 

        int pos = 0;
        while (pos < bsize)
        { 
          int length = 0;
          
          pds_comm.unpack(&recvBuffer[0], bsize, pos, &length, 1);
        
          // find duplicates across procs
          std::pair<unordered_set<std::string>::iterator, bool> result = tmp_device_names.insert(std::string((&recvBuffer[0])+pos, length));
          if (!result.second)
          {
            Report::UserError() << "Duplicate device " << std::string((&recvBuffer[0])+pos, length);
          }
          pos += length;
        }
      }
      else if (proc_rank == proc) 
      {
        // Pack the buffer and send it
        sendBuffer.resize(byteCount);

        int pos = 0;
        for (unordered_set<std::string>::const_iterator it = device_names.begin(), end = device_names.end(); it != end; ++it)
        {
          int length = (*it).size();
          pds_comm.pack(&length, 1, &sendBuffer[0], byteCount, pos);
          pds_comm.pack((*it).c_str(), length, &sendBuffer[0], byteCount, pos);
        }

        pds_comm.send( &byteCount, 1, 0 );
        pds_comm.send( &sendBuffer[0], byteCount, 0 );
      }
    }

    N_ERH_ErrorMgr::safeBarrier(comm);          // All procs call (5+N)
  }
}

//----------------------------------------------------------------------------
// Function       : getLeadCurrentDevices
// Purpose        : 
// Special Notes  :
// Creator        : 
// Creation Date  : 
//----------------------------------------------------------------------------
void getLeadCurrentDevices(const Util::ParamList &variable_list, std::set<std::string> &devicesNeedingLeadCurrents)
{
  for (Util::ParamList::const_iterator iterParam = variable_list.begin() ; iterParam != variable_list.end(); ++iterParam)
  {
    const std::string &varType = iterParam->tag();

    if (Util::hasExpressionTag(*iterParam))
    {

      // this is a do-nothing group
      Teuchos::RCP<Xyce::Util::baseExpressionGroup> exprGroup = 
        Teuchos::rcp(new Xyce::Util::baseExpressionGroup());

      Util::Expression exp(exprGroup, iterParam->tag());

      const std::vector<std::string> & leads = exp.getLeadCurrents();

      // any lead currents found in this expression need to be communicated to the device manager.
      // Multi terminal devices have an extra designator on the name as in name{lead_name}
      // need to remove any {} in the name.
      for (std::vector<std::string>::const_iterator currLeadItr = leads.begin(); currLeadItr != leads.end(); ++currLeadItr)
      {
        size_t leadDesignator = currLeadItr->find_first_of("{");
        devicesNeedingLeadCurrents.insert( currLeadItr->substr(0, leadDesignator));
      }
    }
    else
    {
      if ( ((varType == "I" || (varType.size() == 2 && varType[0] == 'I')) || (varType == "P") || (varType == "W")) &&  (iterParam->getImmutableValue<int>() > 0))
      {
        // any devices found in this I(xxx) structure need to be communicated to the device manager
        // so that the lead currents can be calculated
        if (iterParam->getImmutableValue<int>() != 1)
        {
          Report::UserError0() << "Only one device argument allowed in I(), W() or P() in .print";
        }
        else
        {
          ++iterParam;
          if (varType.size() == 2)
          {
            devicesNeedingLeadCurrents.insert(iterParam->tag());
          }
          else
          {
            devicesNeedingLeadCurrents.insert(iterParam->tag());
          }
        }
      }
    }
  }

  // the list of devices that need lead currents.
  if (DEBUG_IO && !devicesNeedingLeadCurrents.empty())
  {
    std::set<std::string>::iterator currentDeviceNameItr = devicesNeedingLeadCurrents.begin();
    std::set<std::string>::iterator endDeviceNameItr = devicesNeedingLeadCurrents.end();
    Xyce::dout() << "Devices for which lead currents were requested: ";
    while ( currentDeviceNameItr != endDeviceNameItr)
    {
      Xyce::dout() << *currentDeviceNameItr << "  ";
      currentDeviceNameItr++;
    }
  }
}

//----------------------------------------------------------------------------
// Function       : getWildCardLeadCurrentDevices
// Purpose        : Determine which devices need to have lead currents enable
//                  because they match a wildcard specification in an I(), P()
//                  or W() operator.
// Special Notes  :
// Creator        : Pete Sholander, SNL
// Creation Date  : 9/25/2019
//----------------------------------------------------------------------------
void getWildCardLeadCurrentDevices(
  const Util::ParamList &variable_list,
  const unordered_set<std::string> & device_names,
  std::set<std::string> &devicesNeedingLeadCurrents)
{
  for (Util::ParamList::const_iterator iterParam = variable_list.begin() ; iterParam != variable_list.end(); ++iterParam)
  {
    const std::string &varType = iterParam->tag();

    if ( !Util::hasExpressionTag(*iterParam) && ((varType == "I" || (varType.size() == 2 && varType[0] == 'I')) ||
         (varType == "P") || (varType == "W")) &&  (iterParam->getImmutableValue<int>() > 0))
    {
      // any devices that match this I(xx), P(xx) or W(xx) wildcard request need to be
      // communicated to the device manager so that their lead currents can be calculated
      if (iterParam->getImmutableValue<int>() != 1)
      {
        Report::UserError0() << "Only one device argument allowed in I(), W() or P() in .print";
      }
      else
      {
        ++iterParam;
	std::vector<std::string> matches;
	findAllWildCardMatches(iterParam->tag(), device_names, matches);
        devicesNeedingLeadCurrents.insert(matches.begin(),matches.end());
      }
    }
  }

  // the list of devices that need lead currents.
  if (DEBUG_IO && !devicesNeedingLeadCurrents.empty())
  {
    std::set<std::string>::iterator currentDeviceNameItr = devicesNeedingLeadCurrents.begin();
    std::set<std::string>::iterator endDeviceNameItr = devicesNeedingLeadCurrents.end();
    Xyce::dout() << "Devices for which lead currents were requested: ";
    while ( currentDeviceNameItr != endDeviceNameItr)
    {
      Xyce::dout() << *currentDeviceNameItr << "  ";
      currentDeviceNameItr++;
    }
  }
}

  struct X
  {
    const char *first;
    int second;
  };

  static X x[] = {
    {"GLOBAL", -1},
    {"PARALLEL", 1},
    {"SENS", 2},
    {"PRINT", 3}
  };

struct Sorter
{
  using result_type = bool;
  using first_argument_type = Util::OptionBlock;
  using second_argument_type = Util::OptionBlock;

  bool operator()(const Util::OptionBlock &s0, const Util::OptionBlock &s1) const 
  {
    int s0_order = 0;
    int s1_order = 0;

    for (int i = 0; i < sizeof(x)/sizeof(x[0]); ++i) {
      if (compare_nocase(s0.getName().c_str(), x[i].first) == 0)
        s0_order = x[i].second;
      if (compare_nocase(s1.getName().c_str(), x[i].first) == 0)
        s1_order = x[i].second;
    }
    return s0_order < s1_order;
  }
  bool operator()(const Util::OptionBlock *s0, const Util::OptionBlock *s1) const 
  {
    int s0_order = 0;
    int s1_order = 0;

    for (int i = 0; i < sizeof(x)/sizeof(x[0]); ++i) {
      if (compare_nocase(s0->getName().c_str(), x[i].first) == 0)
        s0_order = x[i].second;
      if (compare_nocase(s1->getName().c_str(), x[i].first) == 0)
        s1_order = x[i].second;
    }
    return s0_order < s1_order;
  }
};

//-----------------------------------------------------------------------------
// Function      : registerCircuitOptions
// Purpose       : register options
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
registerCircuitOptions(
  PkgOptionsMgr &                       options_manager,
  std::list<Util::OptionBlock> &        option_block_list)
{
  option_block_list.sort( Sorter() );
  for (std::list<Util::OptionBlock>::iterator it = option_block_list.begin(), end = option_block_list.end(); it != end; ++it)
  {
    if (it->getName() != "DIST")
      options_manager.submitOptions(*it);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : registerDistOptions
// Purpose       : register distribution options
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
registerDistOptions(
  PkgOptionsMgr &                       options_manager,
  std::list<Util::OptionBlock> &        option_block_list)
{
  option_block_list.sort( Sorter() );
  for (std::list<Util::OptionBlock>::iterator it = option_block_list.begin(), end = option_block_list.end(); it != end; ++it)
  {
    if (it->getName() == "DIST")
     options_manager.submitOptions(*it);
  }

  return true;
}

} // namespace IO
} // namespace Xyce

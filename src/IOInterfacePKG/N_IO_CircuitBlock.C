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
// Purpose        : Define the circuit level containers for holding netlist
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

#include <Xyce_config.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include <ctime>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CircuitContext.h>
#include <N_IO_CmdParse.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_ParameterBlock.h>
#include <N_IO_ParsingHelpers.h>
#include <N_TOP_Topology.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Stats.h>
#include <N_UTL_Expression.h>
#include <N_PDS_PackTraits.h>

#include <N_IO_fwd.h>

#include <expressionGroup.h>

namespace Xyce {
namespace IO {

//--------------------------------------------------------------------------
// Function      : CircuitBlock::CircuitBlock
// Purpose       : Constructor
// Special Notes : This constructor is the one used for subcircuiting.
//                 Changed on 10/10/2007 by KRS to accept boolean remove
//                 variables.
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlock::CircuitBlock(
  std::string const &                                           fileName,
  const CmdParse  &                                             command_line,
  HangingResistor &                                             hanging_resistor, 
  CircuitMetadata  &                                            md,
  unordered_set<std::string> &                                  mn,
  std::map<std::string,FileSSFPair>  &                          ssfm,
  std::map<std::string,IncludeFileInfo>  &                      iflm,
  CircuitContext  &                                             cc,
  CircuitBlock*                                                 mainCircPtr,
  CircuitBlock*                                                 parentCircPtr,
  Topo::Topology &                                              topology,
  Device::DeviceMgr &                                           device_manager,
  unordered_set<std::string>  &                                 dNames,
  unordered_set<std::string>  &                                 nNames,
  AliasNodeMap &                                                alias_node_map,
  const std::vector< std::pair< std::string, std::string> >  &  externalNetlistParams,
  std::vector<bool> &                                           pFilter,
  bool                                                          removeRedundant,
  bool                                                          modelBinning,
  double                                                        scale)
: netlistFilename_(fileName),
  title_(""),
  name_(""),
  analysisName_(""),
  nodeNames_(nNames),
  modelNames_(mn),
  deviceNames_(dNames),
  optionsTable_(),
  includeFileLocation_(iflm),
  commandLine_(command_line),
  hangingResistor_(hanging_resistor),
  circuitContext_(cc),
  metadata_(md),
  externalNetlistParams_(externalNetlistParams),
  netlistSave_(true),
  morFlag_(false),
  devProcessedNumber_(0),
  aliasNodeMap_(alias_node_map),
  netlistIn_(0),
  ssfPtr_(0),
  fileStartPosition_(0),
  fileEndPosition_(0),
  lineStartPosition_(1),
  lineEndPosition_(1),
  mainCircuitPtr_(mainCircPtr),
  parentCircuitPtr_(parentCircPtr),
  ssfMap_(ssfm),
  device_(cc,md),
  preprocessFilter_(pFilter),
  remove_any_redundant_(removeRedundant),
  model_binning_flag_(modelBinning),
  lengthScale_(scale),
  topology_(topology),
  deviceManager_(device_manager)
{
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::CircuitBlock
// Purpose       : Constructor
// Special Notes : Use of this constructor coincides with the top-level or
//                 main circuit, the pointer to the main circuit is
//                 set here.  This constructor is only called from
//                 IO_NetlistImportTool.  It is never called from inside of
//                 CircuitBlock.
//
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlock::CircuitBlock(
  const std::string &                                           netlistFilename_In,
  const CmdParse &                                              command_line,
  HangingResistor &                                             hanging_resistor, 
  CircuitMetadata &                                             md,
  unordered_set<std::string> &                                  mn,
  std::map<std::string,FileSSFPair> &                           ssfm,
  std::map<std::string,IncludeFileInfo> &                       iflm,
  CircuitContext &                                              cc,
  Topo::Topology &                                              topology,
  Device::DeviceMgr &                                           device_manager,
  unordered_set<std::string> &                                  dNames,
  unordered_set<std::string> &                                  nNames,
  AliasNodeMap &                                                alias_node_map,
  const std::vector< std::pair< std::string, std::string> > &   externalNetlistParams)
: netlistFilename_(netlistFilename_In),
  title_(""),
  name_(""),
  analysisName_(""),
  nodeNames_(nNames),
  modelNames_(mn),
  deviceNames_(dNames),
  includeFileLocation_(iflm),
  commandLine_(command_line),
  hangingResistor_(hanging_resistor),
  circuitContext_(cc),
  metadata_(md),
  externalNetlistParams_(externalNetlistParams),
  netlistSave_(true),
  morFlag_(false),
  devProcessedNumber_(0),
  aliasNodeMap_(alias_node_map),
  netlistIn_(0),
  ssfPtr_(0),
  fileStartPosition_(0),
  fileEndPosition_(0),
  lineStartPosition_(1),
  lineEndPosition_(1),
  mainCircuitPtr_(this),
  parentCircuitPtr_(NULL),
  ssfMap_(ssfm),
  device_(cc,md),
  preprocessFilter_(PreprocessType::NUM_PREPROCESS, false),
  remove_any_redundant_(false),
  model_binning_flag_(false),
  topology_(topology),
  deviceManager_(device_manager)
{
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::~CircuitBlock
// Purpose       : Destructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlock::~CircuitBlock()
{
  // The original instance of the class is being destroyed.
  unordered_map< std::string, CircuitBlock * >::iterator itcbt = circuitBlockTable_.begin();
  for ( ; itcbt != circuitBlockTable_.end(); ++itcbt )
  {
    delete itcbt->second;
  }
  circuitBlockTable_.clear();
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::parseNetlistFilePass1
// Purpose       : Top level entry for parsing netlist This is where the
//                 library context is initialized.
//
// Special Notes :
//
// Creator       : Dave Shirley, PSSI
//
// Creation Date : 10/08/2009
//--------------------------------------------------------------------------
bool CircuitBlock::parseNetlistFilePass1(PkgOptionsMgr &options_manager)
{
  std::string libSelect;
  std::vector< std::string > libInside;
  return parseNetlistFilePass1(options_manager, libSelect, libInside);
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::parseNetlistFilePass1
// Purpose       : Parse the netlist file. Netlists parsing is a two-phase
//                 operation since the models are needed when determining
//                 how to handle the device lines. In the first phase,
//                 the netlist is read, the type of each line is determined,
//                 and an object is instantiated corresponding to the type
//                 of the line. If the line is not a device line it can
//                 be fully treated in this phase. Device lines are completed
//                 in the second phase.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
bool CircuitBlock::parseNetlistFilePass1(
  PkgOptionsMgr &              options_manager,
  const std::string &          libSelect,
  std::vector< std::string >&  libInside )
{
  bool result = true;

  // If this is the parent circuit, open the netlist file register an
  // instance of SpiceSeparatedFieldTool for the netlist input stream and
  // get the title line. Note: this instance of SpiceSeparatedFieldTool
  // will be shared by subcircuits (if any) that recursively call this
  // method.

  if ( parentCircuitPtr_ == NULL )
  {
    netlistIn_ = new std::ifstream;
    // Open the netlist file.  Using binary to avoid issues with compiler/plat
    // *fstream differences in implementation
    netlistIn_->open( netlistFilename_.c_str(), std::ios::in | std::ios::binary );

    if ( !netlistIn_->is_open() )
    {
      Report::UserError0() << "Could not open netlist file " << netlistFilename_;
      return false;
    }

    ssfPtr_ = new SpiceSeparatedFieldTool(*netlistIn_, netlistFilename_, externalNetlistParams_);
    ssfMap_[netlistFilename_] = FileSSFPair(netlistIn_, ssfPtr_);

    Xyce::IO::readLine( *netlistIn_, title_ );

    // Increment the line number in the SSF object for the netlist to
    // account for the title line.
    ssfPtr_ ->changeCursorLineNumber(1);

    NetlistLocation topLevelLocation( netlistFilename_, lineStartPosition_ );
    circuitContext_.setLocation( topLevelLocation );

    //Added 10/8/07, KRS.  Adding a preprocessing phase to determine which, if
    //any, devices will be selected to remove "redundancy" (i.e., if "resistor"
    //is selected as an option in a .PREPROCESS statement, Xyce will ignore
    //all resistors for which both nodes are the same, i.e.
    // "R1 3 3 1" will not be added to the netlist.  This first phase just
    //extracts the parameters (diodes, capacitors, resistors, etc.) for which
    //we want Xyce to eliminate redundancy.
    //
    // 12/10/07 KRS:  Also using this preprocess phase to detect if we want to
    // add resistors of specified value between ground and "dangling" nodes
    parsePreprocess();

    // Reset the location/line number (parsepreprocess brings us to end of file).
    resetSSFPtr();

    if (commandLine_.argExists("-syntax") || commandLine_.argExists("-count"))
    {
      netlistSave_ = false;
    }
  }

  if (DEBUG_IO) {
    if ( parentCircuitPtr_ == NULL )
    {
      Xyce::dout() << "Pass 1 parsing of netlist file: " << netlistFilename_ << std::endl;
    }
  }

  std::map<std::string,int> fun;
  for (;;) {
    bool line_parsed = true;

    if (handleLinePass1( line_parsed, options_manager, fun, modMap_, libSelect, libInside ) )
      result = result && line_parsed;
    else
      break;
  }

  // Check if the parser thinks it is still inside a library, even though it is done.
  // That means that a user has tried to substitute .LIB for .INC, where the rest of the
  // netlist is ignored in the pursuit of finding .ENDL.  See SON bugs 387 and 980.
  if (( parentCircuitPtr_ == NULL) && libInside.size())
  {
    Report::UserError() << "Could not find .ENDL statement for '.LIB " << libInside.front() 
                        << "'.  Maybe '.LIB " << libInside.front() << " <library_name>' or '.INC " 
                        << libInside.front() << "' was intended.";
    return false;
  }

  if (!result)
    return result;

  // if K lines found, collect coupled inductance data for parsing later
  if( ( parentCircuitPtr_ == NULL ) && ( !(rawMIs_.empty()) ) )
  {
    std::multimap< CircuitContext *, DeviceBlock >::iterator mm =
     rawMIs_.begin();

    if (DEBUG_IO)
      Xyce::dout() << "Total K lines found:  " << rawMIs_.size() << std::endl;

    for( ; mm != rawMIs_.end(); ++mm )
    {
      circuitContext_.setContext( ( *mm ).first );

      // Parses the K line
      // Note: the optional "false" parameter to extractData tells it not
      // to bother trying to do full resolution of expression-valued
      // coupling coefficients.  This is unnecessary at this stage (we'll
      // do that later, when the K line has been turned into a Y), and
      // is actually harmful for coil gun applications, where K may be
      // unresolvable at this stage in processing.
      bool resolveParams=false;
      bool modelBinning=false;
      double scale=1.0;
      ( ( *mm ).second ).extractData(resolveParams,modelBinning,scale);

      // Add mutual inductance to circuit context
      circuitContext_.addMutualInductance( ( *mm ).second );

      if (DEBUG_IO)
        Xyce::dout() << "In Pass 1:  adding: "
                     << ( ( *mm ).second ).getInstanceName() << " with model "
                     << ( ( *mm ).second ).getModelName() << " to "
                     << circuitContext_.getCurrentContextName() << std::endl;

      circuitContext_.restorePreviousContext();
    }

    // temporary blocks no longer needed; free memory
    rawMIs_.clear();
  }

  modMap_.clear();

  // Clean up device counts and unused subcircuits and resolve expressions.

  if ( parentCircuitPtr_ == NULL )
  {
    // Finish resolving Mutual Inductances if necessary
    // This means finding all the assoc. inductors and including their inductance in MI.

    // Note that the parsing of mutual inductances includes converting coupled sets of
    // inductors into devices like YMIL (linear coupling) and YMIN (non-linear
    // coupling) that contain the coupling plus the inductors.  This creates an extra
    // pass between pass1 and pass2 which occurs on proc zero in parallel.

    if( circuitContext_.totalMutualInductanceCount() )
    {
      parseMutualInductances();
    }

    // Get the total device count for the circuit.
    circuitContext_.getTotalDeviceCount();
    
    // Collect the device types from the circuit context if we are only performing syntax checking.
    if ( !netlistSave_ )
    {
      deviceManager_.addDevicesToCount( circuitContext_.getDeviceCountMap() );
    }
 
    // Get the name of all the subcircuit instances that were used in computing the device count. 
    std::vector<std::string>& usedSubcircuits = circuitContext_.getUsedContextList();
    std::sort( usedSubcircuits.begin(), usedSubcircuits.end() );
    usedSubcircuits.erase(std::unique(usedSubcircuits.begin(), usedSubcircuits.end()), usedSubcircuits.end());

    // Prune unused subcircuits from the CircuitContext list (so they don't get communicated to other processors)
    circuitContext_.pruneContexts( usedSubcircuits );

    // Collect any IC or NODESET option blocks that exist in the used subcircuits.
    if (usedSubcircuits.size() > 0)
    {
      // Recurse through the circuitBlockTable_ to retrieve any .ICs or .NODESETs for
      // the usedSubcircuits.
      std::vector<Util::OptionBlock> icNodesetOB;
      getICNodesetList( usedSubcircuits, icNodesetOB );

      optionsTable_.insert( optionsTable_.end(), icNodesetOB.begin(), icNodesetOB.end() );
    }

    // Resolve current context parameters.
    std::vector<Device::Param> params;
    result = circuitContext_.resolve(params);
    if (!result)
    {
      return result;
    }

    // resolve any functions and parameters in expression on the print line
    // at this point (or for that matter any functions/parameters in
    // the optionsTable data (so ".OP" ".OPTIONS" ".OUTPUT"
    // ".PRINT"".TRAN"  ".STEP"  ".RESULT" ".IC" ".DCVOLT"
    //  ".NODESET" ".SAVE" ".LOAD" ".MPDE" ".HB"  ".AC" ".MEASURE" ".MEAS")
    // This could happen later, as in the actual classes that handle the above
    // functions, however at this stage we have all the contextual information
    // to resolve this without duplicating code elsewhere.
    resolveExpressionsInOptionBlocks();

    // update the aliasNodeMapHelper_ map to include "subcircuit interface node names"
    // that were embedded within expressions
    updateAliasNodeMapHelper();

    // Sanity check the analysis and print line(s)
    result = handleAnalysis();
    if (!result)
      return result;

    if (DEBUG_IO) {
      Xyce::dout() << "Done with pass 1 netlist file parsing" << std::endl;
    }
  }

  return true; // Only get here on success.
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::writeOutNetlist
// Purpose        : Write out netlist, if requested.
// Special Notes  :
// Scope          :
// Creator        : Heidi Thornquist
// Creation Date  : 07/25/2014
//----------------------------------------------------------------------------
void CircuitBlock::writeOutNetlist()
{
  // Here's where we call netlist copy stuff
  if (hangingResistor_.getNetlistCopy())
  {
    produceUnflattenedNetlist();

    if (DEBUG_IO)
      Xyce::dout() << "Done writing preprocessed netlist" << std::endl;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::parseMutualInductances
// Purpose        : Special Pass for Mutual Inductances
// Special Notes  :
// Scope          :
// Creator        : Rob Hoekstra
// Creation Date  : 08/27/04
//----------------------------------------------------------------------------
bool CircuitBlock::parseMutualInductances()
{
  if (DEBUG_IO) 
  {
    if (parentCircuitPtr_ == NULL)
    {
      Xyce::dout() << "Pass MI parsing of netlist file: " << netlistFilename_ << std::endl;
    }
  }

  // Set the start location of the circuit or subcircuit in its associated file.
  // If this is the main circuit, skip over the title line and continue.
  resetSSFPtr();

  if( circuitContext_.haveMutualInductances() )
  {
    while( getLinePassMI() ) {}

    // retrieve tables and MI references from current circuit context
    std::vector<CircuitContext::MutualInductance> & MIs =
     circuitContext_.getMutualInductances();
    std::vector< std::set< std::string > > & iTable = circuitContext_.getSharedInductorTable();
    std::vector< std::vector< int > > & mTable = circuitContext_.getAllIndexedMIs();
    std::set< std::string > & cTable = circuitContext_.getAllCoupledInductors();
    std::map<std::string,std::string>::iterator nIter;
    std::map<std::string,std::string>::iterator nIter_end;
    int numMIs = MIs.size();
    int doneKey = 1;
    bool done = false;
    int imin = 0;


    mTable.push_back( std::vector< int >() );
    iTable.push_back( std::set< std::string >() );

    while (!done )
    {
      std::set<std::string> indUsed;

      mTable.push_back( std::vector< int >() );
      iTable.push_back( std::set< std::string >() );

      //Add the information for imin'th mutual inductor to all of the tables:
      MIs[imin].sharedKey = doneKey;
      mTable[doneKey].push_back(imin);
      nIter = MIs[imin].inductors.begin();
      nIter_end = MIs[imin].inductors.end();
      for( ; nIter != nIter_end; ++nIter )
      {
        indUsed.insert((*nIter).first);
        cTable.insert((*nIter).first);
        iTable[doneKey].insert((*nIter).first);
      }

      for( int i = imin + 1; i < numMIs; ++i )
      {
        if( MIs[i].model == "" && MIs[i].sharedKey == 0)
        {
          bool addit = false;

          //Search to see if the i'th mutual inductor contains inductors that
          //are already contained in the imin'th mutual inductor.  If so, add
          //these inductances to the imin'th mutual inductor
          nIter = MIs[i].inductors.begin();
          nIter_end = MIs[i].inductors.end();
          for( ; nIter != nIter_end; ++nIter )
          {
            if (indUsed.find((*nIter).first) != indUsed.end())
              addit = true;
          }

          if (addit)
          {
            MIs[i].sharedKey = doneKey;
            mTable[doneKey].push_back(i);
            nIter = MIs[i].inductors.begin();
            nIter_end = MIs[i].inductors.end();
            for( ; nIter != nIter_end; ++nIter )
            {
              indUsed.insert((*nIter).first);
              cTable.insert((*nIter).first);
              iTable[doneKey].insert((*nIter).first);
            }
            i = imin;
          }
        }
          //Whenever we add a new mutual inductor to the imin'th mutual
          //inductor, we have to start the loop all over again to make sure we
          //haven't missed anything.  That's why we set i=imin after we go
          //through the "addit" process, so that we start the loop again at
          //i = imin + 1.
      }

      done=true;

      //If there are no more mutual inductors with sharedKey tag = 0, then
      //we're done.  But, otherwise, we have to repeat this process.  imin is
      //now set to the first mutual inductor index which has a sharedKey tag of
      //0

      for (int i = 0; i < numMIs; i++)
      {
        if (MIs[i].sharedKey == 0)
        {
          imin = i;
          doneKey++;
          done=false;
          break;
        }
      }
    }

    //Now that all the MIs have been broken up into Y-devices, we have to
    //make a correction to the total number of devices.  Before we get to this
    //point, Xyce computes the total number of devices as simply being the
    //total number of device lines it encounters in parsing the netlist file,
    //but when we lump the coupled inductors and multiple K-devices into big
    //Y devices, the number of total devices decreases.  Essentially, what we
    //do here is the following:  we take the total device count, subtract the
    //total number of coupled inductors, subtract the total number of K-lines
    //found in the netlist file, and add on the number of Y-devices that we
    //formed in the above loop.

    int totalCoupledIs = 0;
    int ilength = iTable.size();

    for (int i=0; i < ilength; i++)
    {
      totalCoupledIs += iTable[i].size();
    }

    circuitContext_.augmentTotalDeviceCount(numMIs,totalCoupledIs,doneKey);

    // package MIs for distribution later
    circuitContext_.bundleMIs();
  }

  unordered_map<std::string, CircuitBlock*>::iterator itcbt = circuitBlockTable_.begin();
  for( ; itcbt != circuitBlockTable_.end(); ++itcbt )
  {
    CircuitBlock * subcircuitPtr = itcbt->second;

    // Locate the subcircuit in the netlist file. It can either be in
    // the file currently being read, or in a separate include file.
    if(subcircuitPtr->netlistFilename_ != netlistFilename_)
    { // The subcircuit is in an include file.
      // Get SSF from Pass 1's ssf map
      if( ssfMap_.count( subcircuitPtr->netlistFilename_ ) )
        subcircuitPtr->setSSFPtr( ssfMap_[subcircuitPtr->netlistFilename_].second );
      else
      {
        Report::UserError() << "Can't find include file " << subcircuitPtr->netlistFilename_;
        return false;
      }
    }
    else
      subcircuitPtr->setSSFPtr( ssfPtr_ );

    // Set the position of the subcircuit in its file.
    subcircuitPtr->setFilePosition(subcircuitPtr->getStartPosition());
    subcircuitPtr->setLinePosition( subcircuitPtr->getLineStartPosition() );

    // switch context
    circuitContext_.setContext( subcircuitPtr->getName() );

    // Parse subckt's MIs
    subcircuitPtr->parseMutualInductances();

    // restore context
    circuitContext_.restorePreviousContext();
  }

  if (DEBUG_IO) {
    print();

    if (parentCircuitPtr_ == NULL)
    {
      Xyce::dout() << "Done with pass MI netlist file parsing" << std::endl;
    }
  }

  return true; // Only get here on success.
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::print
// Purpose       : Print the circuit.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/08/2001
//--------------------------------------------------------------------------
void CircuitBlock::print()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << std::endl;
  Xyce::dout() << std::endl << Xyce::section_divider << std::endl;
  Xyce::dout() << "CircuitBlock::print" << std::endl;
  if ( parentCircuitPtr_ == NULL )
  {
    Xyce::dout() << "Circuit Title: " << title_ << std::endl;
  }
  else 
  {
    Xyce::dout() << "Subcircuit Name:  " << name_ << std::endl;
  }
  Xyce::dout() << std::endl;

  if ( !optionsTable_.empty() )
  {
    Xyce::dout() << "Options: " << std::endl;
    std::list<Util::OptionBlock>::iterator optionIter = optionsTable_.begin();
    std::list<Util::OptionBlock>::iterator optionIterEnd = optionsTable_.end();
    for ( ; optionIter != optionIterEnd; ++optionIter )
    {
      Xyce::dout() << std::endl
                   << "Option Information" << std::endl
                   << "------------------" << std::endl
                   << std::endl
                   << "  name: " << optionIter->getName() << std::endl;

      Xyce::dout() << "  parameters: " << std::endl;
      Util::ParamList::const_iterator paramIter = optionIter->begin();
      Util::ParamList::const_iterator paramIterEnd = optionIter->end();
      for ( ; paramIter != paramIterEnd; ++paramIter )
      {
        Xyce::dout() << "  " << paramIter->tag() << "  ";
        Xyce::dout() << paramIter->stringValue() << std::endl;
      }
    }
    Xyce::dout() << std::endl << std::endl;
  }

  if ( !circuitBlockTable_.empty() )
  {
    Xyce::dout() << "Subcircuits: " << std::endl;
    unordered_map< std::string, CircuitBlock * >::iterator itcbt = circuitBlockTable_.begin();
    for ( ; itcbt != circuitBlockTable_.end(); ++itcbt )
    {
      itcbt->second->print();
    }
    Xyce::dout() << "End Subcircuits" << std::endl;

    Xyce::dout() << std::endl;
  }
  Xyce::dout() << std::endl << Xyce::section_divider << std::endl;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setStartPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 05/19/2003
//----------------------------------------------------------------------------
void CircuitBlock::setStartPosition()
{
  fileStartPosition_ = ssfPtr_->getFilePosition();
  lineStartPosition_ = ssfPtr_->getLineNumber();

  if (DEBUG_IO) {
    Xyce::dout() << "CircuitBlock::setStartPosition being called for file "
                 << ssfPtr_->getFileName() << std::endl
                 << "  start position  = " <<  fileStartPosition_ << std::endl
                 << "  start line  = " << lineStartPosition_ << std::endl;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setEndPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/25/2003
//----------------------------------------------------------------------------
void CircuitBlock::setEndPosition()
{
  fileEndPosition_ = ssfPtr_->getFilePosition();
  lineEndPosition_ = ssfPtr_->getLineNumber();

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlock::setEndPosition:" << std::endl
                 << "  Setting file end position: " << fileEndPosition_ << std::endl
                 << "  setting line end position: " << lineEndPosition_ << std::endl;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setFilePosition
// Purpose        : Set the location in the input file to the given position.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/21/2003
//----------------------------------------------------------------------------
void CircuitBlock::setFilePosition(std::streampos const& position)
{

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlock::setFilePosition: Setting file position to "
                 << position <<std::endl;

  ssfPtr_->setLocation(position);
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setLinePosition
// Purpose        : Set the location of the currentLine counter in the ssfPtr
// Special Notes  :
// Scope          : public
// Creator        : Eric Rankin
// Creation Date  : 10/13/2004
//----------------------------------------------------------------------------
void CircuitBlock::setLinePosition( int const& position )
{
  ssfPtr_->setLineNumber( position );
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::extractSubcircuitData
// Purpose       : Extract subcircuit data from parsedLine. The bulk of
//                 the subcircuit data is stored in the circuit context.
//                 A circuit block is created to represent the subcircuit
//                 but mainly exists to help with the recursive descent
//                 through the circuit during pass 2 to instantiate devices.
//                 All that is needed for this is the subcircuit name.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/21/2001
//--------------------------------------------------------------------------
bool CircuitBlock::extractSubcircuitData(std::string fileName, 
                                         TokenVector const& parsedLine)
{
  int numFields = parsedLine.size();

  if ( numFields < 3 )
  {
    Report::DevelFatal0().in("CircuitBlock::extractSubcircuitData").at(fileName, parsedLine[0].lineNumber_)
      << "This should have been detected earlier in parse";
    return false;
  }

  // Extract the subcircuit name.
  ExtendedString field ( parsedLine[1].string_ );
  field.toUpper();
  name_ = field;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addTableData
// Purpose       : Add a device to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addTableData( DeviceBlock & device )
{
  if (DEBUG_IO) {
    ++devProcessedNumber_;
    if (devProcessedNumber_%1000 == 0)
      Xyce::dout() << devProcessedNumber_ << "Devices processed" << std::endl;
  }

  const std::string& dName = device.getDeviceData().getDevBlock().getInstanceName().getEncodedName();

  std::pair<unordered_set<std::string>::iterator, bool> result = deviceNames_.insert(dName);
  if (!result.second) 
  {
    Report::UserError().at(device.getDeviceData().getDevBlock().getNetlistLocation()) << "Duplicate device "<< dName;
  }

  if (netlistSave_)
  {
    if (DEBUG_IO) {
      int lastColon = dName.find_last_of( ':' );
      Xyce::dout() << "Inserting device: " << dName
                   << " locally named " <<  dName.substr( lastColon + 1 )
                   << " from file: " << device.getDeviceData().getDevBlock().getNetlistLocation().getFilename()
                   << " line: " << device.getDeviceData().getDevBlock().getNetlistLocation().getLineNumber() << std::endl;
    }

    topology_.addDevice(deviceManager_, device.getDeviceData());
  }

  // save node names for syntax diagnostics
  const std::vector<std::string> &nodeList = device.getDeviceData().get_NodeList();
  for (std::vector<std::string>::const_iterator it = nodeList.begin(), end = nodeList.end(); it != end; ++it)
  {
    if (DEBUG_IO)
      Xyce::dout() << "  Node: " << *it << std::endl;

    nodeNames_.insert(*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addMutualInductor
// Purpose       : Add a mutual inductor device to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addMutualInductor( DeviceBlock& device, CircuitContext* context )
{
  rawMIs_.insert(std::pair< CircuitContext *, DeviceBlock >
                ( context, device ) );
}


//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addModel
// Purpose       : Add a model to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addModel( const ParameterBlock * modelPtr, std::string const& modelPrefix)
{
  std::string modelName(modelPtr->getName());
  if (modelPrefix != "")
  {
    modelName = modelPrefix + ":" + modelName;
  }

  std::pair<unordered_set<std::string>::iterator,bool> ret = modelNames_.insert(modelName);

  if ( ret.second )
  {
    // Copy the parameter block for this model.
    tmpModel_ = (*modelPtr);
    tmpModel_.setName(modelName);

    // Set the model parameter values.
    tmpModel_.setParameterValues(&circuitContext_);

    deviceManager_.addDeviceModel(tmpModel_.modelData);
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addParams
// Purpose       : Add param to circuitContext params (.PARAM)
// Special Notes :
// Creator       : David G. Baur
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addParams(const Util::OptionBlock & options)
{
  circuitContext_.addParams(options.begin(), options.end());
}

//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addGlobalParams
// Purpose       : Add param to circuitContext global params (.GLOBAL_PARAM)
// Special Notes :
// Creator       : David G. Baur
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addGlobalParams(const Util::OptionBlock & options)
{
  circuitContext_.addGlobalParams(options.begin(), options.end());
}

//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addOptions
// Purpose       : Add a set of options corresponding to a .OPTIONS netlist
//                 line to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addOptions(const Util::OptionBlock & options)
{
  const std::string &name = options.getName();
  
  // This handles .PRINT, .MEASURE, .FOUR, .IC and .NODESET lines.
  // Note that option blocks that have name="MEASURE" are actually 
  // from .OPTIONS MEASURE lines and not .MEASURE lines.
  if ( (name == "PRINT") || (name == "DOT_MEASURE_LINE") || (name == "FOUR") || 
       (name == "IC") || (name == "NODESET") )
  {
    // Update the aliasNodeMapHelper_ map with any parameters that start with X 
    // (subcircuit interface node names) or that are an expression.  The expressions 
    // will be parsed later, for embedded subcircuit interface node names, 
    // in updateAliasNodeMapHelper().  See SRN Bug 1962 for more details.
    for (Util::ParamList::const_iterator it = options.begin(), end = options.end(); it != end; ++it)
    {
      std::string aliasStr = std::string((*it).uTag());
      if ( (aliasStr[0]=='{') || (aliasStr[0]=='X') )
          aliasNodeMapHelper_.insert(aliasStr);
    }
  }
  else if (name == "RESULT") 
  {
    // Update the aliasNodeMapHelper_ map with the expression given on the .RESULT
    // line. It will be parsed later, for embedded subcircuit interface node names, 
    // in updateAliasNodeMapHelper().
    for (Util::ParamList::const_iterator it = options.begin(), end = options.end(); it != end; ++it)
    {
      if  ( (std::string((*it).uTag()) == "EXPRESSION") && ((*it).hasExpressionValue()) )
      {
        aliasNodeMapHelper_.insert((*it).stringValue());
      }
    }
  }
  else if (name == "SENS")
  {
    // Update the aliasNodeMapHelper_ map with the expressions given on the .SENS
    // line. They will be parsed later, for embedded subcircuit interface node names, 
    // in updateAliasNodeMapHelper().
    for (Util::ParamList::const_iterator it = options.begin(), end = options.end(); it != end; ++it)
    {
      if  ( (std::string((*it).uTag(),0,7) == "OBJFUNC") && ((*it).hasExpressionValue()) )
      {
        aliasNodeMapHelper_.insert((*it).stringValue());
      }
      else if  ( (std::string((*it).uTag(),0,7) == "OBJVARS") && ((*it).hasExpressionValue()) )
      {
        aliasNodeMapHelper_.insert((*it).stringValue());
      }
    }
  }
  else if (name == "OUTPUT-LINE")
  {
    // The parameters on a .OUTPUT line need to be added to the option block
    // with the name "OUTPUT" which was created by a prior .OPTIONS OUTPUT
    // line in the netlist. Find the option block, report an error if not
    // found.
    std::list<Util::OptionBlock>::iterator it = optionsTable_.begin();
    std::list<Util::OptionBlock>::iterator end = optionsTable_.end();
    for (; it != end; ++it)
      if (it->getName() == "OUTPUT")
        break;

    Util::ParamList::const_iterator paramIter = options.begin();

    if (it == optionsTable_.end())
    {
      // The line number of the .OUTPUT line was stored as the 3rd parameter.
      paramIter++;
      paramIter++;
      int lineNum = (*paramIter).getImmutableValue<int>();

      // Could not find required option block, report error.
      Report::UserError0().at(netlistFilename_, lineNum) << "A .OPTIONS OUTPUT line is required before any .OUTPUT line in the netlist";
    }
    else
    {
      // If we get here, all is well, add the parameters.
      it->addParam(*paramIter);
      paramIter++;
      it->addParam(*paramIter);
    }
    return;
  }

  optionsTable_.push_back(options);
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::handleLinePass1
// Purpose       : Determine the type of netlist line in parsedLine and
//                 handle appropriately.
//
// Special Notes : Returns false on end of file
//
// Creator       : Lon Waters
//
// Creation Date : 09/08/2001
//--------------------------------------------------------------------------
bool CircuitBlock::handleLinePass1(
  bool &                        result,
  PkgOptionsMgr &               options_manager,
  std::map<std::string, int> &  fun,
  ModelMap &                    modMap,
  const std::string &           libSelect,
  std::vector< std::string > &  libInside)
{
  result = true;

  char lineType;
  TokenVector line;
  ExtendedString ES1 ( " " ); 

  // For the first pass, we can peek at the first character of the next usable line
  // and see if we actually need to tokenize the line.
  int eof = ssfPtr_->peekAtNextLine( lineType );
  if (eof)
    return !eof;

  // If we are removing redundant devices, or this isn't a device line, get the line.
  if (remove_any_redundant_ || lineType < 'A' || lineType > 'Z' 
      || lineType == 'X' || lineType == 'K' || lineType == 'U' || lineType == 'Y')
  {
    // Get the next line of input.
    eof = !ssfPtr_->getLine(line,preprocessFilter_[PreprocessType::REPLACE_GROUND]); // Breaks the line into fields.
  
    // If blank line, nothing to do here.
    if ( line.empty() )
    {
      return !eof;
    }

    // Get the longer, uppercase, first word
    ES1 = line[0].string_;
    ES1.toUpper();
  } 
  else
  {
    // The lineType will be enough to process this line, so skip to the end.
    ssfPtr_->skipToEndOfLine();
 
    // Remove any continuation lines.
    char contChar;
    eof = ssfPtr_->peekAtNextLine( contChar );
    while (!eof && contChar == '+')
    {
      ssfPtr_->skipToEndOfLine();
      eof = ssfPtr_->peekAtNextLine( contChar );
    }
  }

  // Handle selecting lines from a library file.
  // If the lines are from a library that is not selected, then ignore them.
  std::string libInsideHere = "";
  if (libInside.size())
  {
    libInsideHere = libInside.back();
  }

  if (libSelect != libInsideHere && libInsideHere != "" && ES1 != ".ENDL")
  {
    return !eof;
  }

  // This is a device, only get the line if necessary.
  bool removecomponent = false;
  if (lineType >= 'A' && lineType <= 'Z')
  {
    if (lineType == 'C' || lineType == 'D' || lineType == 'I' ||
        lineType == 'L' || lineType == 'R' || lineType == 'V')
    {
      if (line.size() > 2 && remove_any_redundant_) //make sure that there are two nodes to check!
      {
        ExtendedString node1 ( line[1].string_ );
        ExtendedString node2 ( line[2].string_ );
        node1.toUpper();
        node2.toUpper();

        removecomponent = IO::removeTwoTerminalDevice(preprocessFilter_, lineType, node1, node2);
      }
    }
    else if (lineType == 'M' || lineType == 'Q')
    {
      if (line.size() > 3 && remove_any_redundant_) //make sure that there are three nodes to check!
      {
        ExtendedString node1 ( line[1].string_ );
        ExtendedString node2 ( line[2].string_ );
        ExtendedString node3 ( line[3].string_ );
        node1.toUpper();
        node2.toUpper();
        node3.toUpper();

        removecomponent = IO::removeThreeTerminalDevice(preprocessFilter_, lineType, node1, node2, node3);
      }
    }

    if (!removecomponent)
    {
      std::string ES2(1, lineType);
      if (lineType != 'X' && lineType != 'K')
      {
        // Get the device type as a single character or character string.
        if (lineType == 'Y')
        {
          ES2 = ES1.substr(1);
        }
        
        // Increment device count.
        circuitContext_.incrementDeviceCount( ES2 );

        // Get the configuration for this device type.
        const Device::Configuration *configuration = Device::Configuration::findConfiguration( ES2, 1 );
        if (configuration)
        {
          if (configuration->getLinearDevice())
          {
            circuitContext_.incrementLinearDeviceCount();
          }         
        }
      } 
      if (lineType == 'X')
      {
        std::string modelName;
        Xyce::IO::extractSubcircuitModelName( line, modelName ); 
        circuitContext_.addInstance(modelName,ES1,netlistFilename_,line[0].lineNumber_);
      }
      if (lineType == 'K')
      {
        circuitContext_.incrementDeviceCount( ES2 );
        DeviceBlock device(circuitContext_, metadata_, netlistFilename_, line);

        // save this information for later resolution
        mainCircuitPtr_->addMutualInductor( device, circuitContext_.getCurrentContextPtr() );
      }

      if (lineType == 'U')
      {
        if (line.size() < 2)
        {
          Report::UserError().at(netlistFilename_, line[0].lineNumber_)
            << "U device line specified with only one token: " <<  ES1 << ". Need at least two.";
          result = false;
        }
        else {
          if (line.size() < 3)
          {
            Report::UserWarning().at(netlistFilename_, line[0].lineNumber_)
              << "U device line (" << ES1 << ") specified with only two tokens. Likely need more.";
          }
        }
      }

      if (lineType == 'Y')
      {
        if (line.size() < 2)
        {
          Report::UserError().at(netlistFilename_, line[0].lineNumber_)
            << "Y device line specified with only one token: " <<  ES1 << ". Need at least two.";
          result = false;
        }
        else {
          if (line.size() < 3)
          {
            Report::UserWarning().at(netlistFilename_, line[0].lineNumber_)
              << "Y device line (" << ES1 << ") specified with only two tokens. Likely need more.";
          }
        }
      }
    }

    else if (DEBUG_IO)
    {
      Xyce::dout() << "Netlist Parse 1:  ";
      Xyce::dout() << "removing component " << ES1 << ".  All nodes on the device";
      Xyce::dout() << " are the same."  << std::endl;
    }
  }

  else if (lineType == '.')
  {
    if (extractData(options_manager, *this, netlistFilename_, line))
      ;
    else if (ES1 == ".GLOBAL")
    {
      if (parentCircuitPtr_ != NULL )
      {
        Report::UserError().at(netlistFilename_, line[0].lineNumber_)
          << "Attempt to assign global node inside of subcircuit";
      }
      if (line.size() != 2)
      {
        Report::UserError().at(netlistFilename_, line[0].lineNumber_)
          << "Syntax error in .global, should be .global <node>";
      }
      ExtendedString ES2(line[1].string_);
      ES2.toUpper();
      circuitContext_.addGlobalNode ( ES2 );
    }

    else if (ES1 == ".FUNC")
    {
      // Create a FunctionBlock for the .FUNC line, extract the
      // data from line, and add the function to the circuit.
      FunctionBlock function( netlistFilename_, line );
      circuitContext_.addFunction( function );
      ExtendedString F ( line[1].string_ );
      F.toUpper();
      if (fun[F]++ != 0)
      {
        Report::UserError().at(netlistFilename_, line[0].lineNumber_)
          << "Duplicate function definition " <<  F;
        result = false;
      }
    }

    else if (ES1 == ".PARAM")
    {
      return Xyce::IO::extractParamData( *this, netlistFilename_, line );
    }

    else if (ES1 == ".GLOBAL_PARAM")
    {
      return Xyce::IO::extractGlobalParamData( *this, netlistFilename_, line );
    }

    else if (ES1 == ".END")
    {
      std::string name = circuitContext_.getCurrentContextPtr()->getName();
      while (circuitContext_.endSubcircuitContext()) {
        Report::UserError0().at(netlistFilename_, line[0].lineNumber_) << "Subcircuit " << name << " missing .ENDS";
        setEndPosition();
        name = circuitContext_.getCurrentContextPtr()->getName();
      }
      return false;
    }

    else if (ES1 == ".ENDS")
    {
      // End the current subcircuit context.
      if (circuitContext_.endSubcircuitContext()) {
        setEndPosition();
        return false;
      }
      else
      {
        Report::UserWarning0().at(netlistFilename_, line[0].lineNumber_) << "Subcircuit .ENDS without .SUBCKT, ignoring";
      }
    }

    else if (ES1 == ".INCLUDE" || ES1 == ".INCL" || ES1 == ".INC" || ES1 == ".LIB")
    {
      // HSPICE documents .INC, .INCL and .INCLUDE as being a valid .INC line
      std::string includeFile, libSelect_new = libSelect, libInside_new;
      Xyce::IO::handleIncludeLine( netlistFilename_, line, 
                                   ES1, includeFile, libSelect_new, libInside_new );

      // Check for recursive .INC/.INCLUDE of same file, which will create an infinite loop.
      if (includeFile == netlistFilename_ && libSelect_new == "" && libInside_new == "")
      {
        Report::UserError().at(netlistFilename_, line[0].lineNumber_) << "Recursive inclusion of same file results in infinite loop";
        return false;
      }
      if (libInside_new != "")
      {
        libInside.push_back( libInside_new );
      }
 
      if (includeFile != "")
      {
        // Get status of context before the include file is parsed to see if we will need to look at the file later.
        std::pair<int,bool> beforeCount = circuitContext_.getCurrentContextPtr()->getDeviceCount();
        int beforeSubckts = (circuitContext_.getCurrentContextPtr())->getSubcktList().size();
        int beforeMod = modMap.size();
        int beforeSubcktDef = circuitBlockTable_.size();

        result = parseIncludeFile(options_manager, includeFile, libSelect_new, fun, modMap);
  
        // Load include file information into struct for future use.
        IncludeFileInfo info;
        info.numDevices = (circuitContext_.getCurrentContextPtr())->getDeviceCount().first - beforeCount.first;
        info.numSubckts = (circuitContext_.getCurrentContextPtr())->getSubcktList().size() - beforeSubckts;
        info.numModels = modMap.size() - beforeMod;
        info.numSUBCKTdefs = circuitBlockTable_.size() - beforeSubcktDef;
        info.inSUBCKT = isSubcircuit();
        info.location = NetlistLocation( netlistFilename_, line[0].lineNumber_ );
        info.parentSUBCKT = circuitContext_.getCurrentContextName();
        includeFileLocation_[ includeFile ] = info;

        if (DEBUG_IO) {
          std::cout << "Include file " << includeFile << " has been parsed : " << std::endl;
          std::cout << " Number of devices = " << info.numDevices << ", number of X lines " << info.numSubckts
                    << ", number of models = " << info.numModels << ", number of SUBCKTs " << info.numSUBCKTdefs << std::endl;
          if ( info.inSUBCKT )
            std::cout << "THIS FILE IS INCLUDED IN .SUBCKT DEFINITION " << info.parentSUBCKT << " !" << std::endl;
        }
      }
    }
    else if (ES1 == ".ENDL")
    {
      Xyce::IO::handleEndlLine ( netlistFilename_, line, libInsideHere );
      if (libInside.size())
        libInside.pop_back();
    }
    else if (ES1 == ".INITCOND" )
    {
      handleInitCond( line );
    }
    else if (ES1 == ".MODEL")
    {
      ParameterBlock* modelPtr =
        new ParameterBlock (netlistFilename_, line);

      ExtendedString M ( line[1].string_ );
      M.toUpper();
      std::map<std::string,ParameterBlock*,LessNoCase>::iterator mp = modMap.find(M);
      if (mp == modMap.end())
      {
        // Save for potential use later, like if there are data for other temperatures
        modMap[M] = modelPtr;
        // Add the model to the circuit context. Note, the circuit context
        // will handle deleting the model.
        circuitContext_.addModel(modelPtr);
        // Store the level number for device configuration
        (*mainCircuitPtr_).levelSet_.insert( modelPtr->getLevel() );
      }
      else
      {
        // A duplicate model name has been detected.  Hopefully this is a data
        // point at another temperature, or other independent parameter.
        ParameterBlock *pb = modMap[M];
        if (pb->getType() == modelPtr->getType() && pb->getLevel() == modelPtr->getLevel())
        {
          std::vector<Device::Param> addMP;
          addMP.push_back(Device::Param("INDEPENDENT;PARAM","TNOM"));
          pb->addParameters(addMP);
          pb->addParameters(modelPtr->getParams());
          delete modelPtr;
        }
        else
        {
          Report::UserError().at(netlistFilename_, line[0].lineNumber_)
            << "Duplicate model definition " << M;
          result = false;
        }
      }
    }
    else if (ES1 == ".SAVE")
    {
      Report::UserWarning0().at(netlistFilename_, line[0].lineNumber_)
        << ".SAVE line not handled properly, statement skipped";
    }
    else if (ES1 == ".SUBCKT")
    {
      // Create a new CircuitBlock to hold the subcircuit definition.
      // Set the parentCircuitPtr of the new CircuitBlock.
      CircuitBlock* subcircuitBlockPtr =
        new CircuitBlock(
          netlistFilename_,
          commandLine_,
          hangingResistor_,
          metadata_,
          modelNames_,
          ssfMap_,
          includeFileLocation_,
          circuitContext_,
          mainCircuitPtr_,
          this,
          topology_,
          deviceManager_,
          deviceNames_,
          nodeNames_,
          aliasNodeMap_,
          externalNetlistParams_,
          preprocessFilter_,
          remove_any_redundant_,
          model_binning_flag_,
          lengthScale_);

      // Start a new subcircuit context in the circuit context.
      result = circuitContext_.beginSubcircuitContext(netlistFilename_, line);

      // Subcircuits must use the ssfPtr_ that was set up by the
      // main circuit for parsing the input file.
      subcircuitBlockPtr->setSSFPtr(ssfPtr_);
      subcircuitBlockPtr->setStartPosition();

      // Extract the subcircuit data from line.
      if (result) {
        result = subcircuitBlockPtr->extractSubcircuitData(netlistFilename_, line)
                 && result;

        ExtendedString S ( line[1].string_ );
        S.toUpper();
        if ( circuitBlockTable_.find( S ) != circuitBlockTable_.end() )
        {
          Report::UserError().at(netlistFilename_, line[0].lineNumber_)
            << "Duplicate subcircuit definition detected: " <<  S;
          result = false;
        }

        circuitBlockTable_[S] = subcircuitBlockPtr;
      }

      result = subcircuitBlockPtr->parseNetlistFilePass1(options_manager, libSelect, libInside)
               && result;
    }

    else if (ES1 == ".PREPROCESS")
    {
      result=true; //ignore this line; it's job is done in the preprocess
      //phase
    }
    else if (ES1 == ".ENDDATA")
    {
      result=true; //ignore this line;  it isn't needed for Xyce to parse .DATA, 
      // but is part of HSpice syntax, so will often appear in netlists that use it.
    }
    else
    {
      // If we get here then we have an unrecognized "." line, flag it
      // with a warning, ignore it and continue.
      Report::UserWarning0().at(netlistFilename_, line[0].lineNumber_)
        <<  "Unrecognized dot line will be ignored";
    }
  }
  else if (lineType == '*' || lineType == ' ' || lineType == '\t')
  {
    result = true;
  }
  else
  {
    Report::UserError().at(netlistFilename_, line[0].lineNumber_)
      << "Unrecognized line";
    result = false;
  }
  return !eof;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getLinePassMI
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 08/01/2003
//----------------------------------------------------------------------------
bool CircuitBlock::getLinePassMI()
{
  int eof = 0;

  TokenVector line;

  while (!eof)
  {
    // Breaks the line into fields.
    eof = !ssfPtr_->getLine(line,preprocessFilter_[PreprocessType::REPLACE_GROUND]);
    //FLAG  Should I be using the optional boolean argument here?  I think so.
    if (DEBUG_IO) 
    {
      Xyce::dout() << "pass MI read netlist line: ";
      for (unsigned int i = 0; i < line.size(); ++i)
      {
        Xyce::dout() << line[i].string_ << " ";
      }
      Xyce::dout() << std::endl;
    }

    // better not try to do anything if getLine returned an empty line!
    if ( !(line.empty()) )
    {
      // Determine what to do with the parsed line.
      char lineType;
      ExtendedString ES1 ( line[0].string_ );
      ES1.toUpper();
      lineType = ES1[0];

      if (lineType == 'L')
      {
        // This is an inductor, check for assoc. MI
        std::vector<CircuitContext::MutualInductance> & MIs =
          circuitContext_.getMutualInductances();
        int numMIs = MIs.size();
        for( int i = 0; i < numMIs; ++i )
        {
          if( MIs[i].inductors.count( ES1 ) )
          {
            device_.clear();
            bool resolveParams=false;
            bool modelBinning=false;
            double scale=1.0;
            device_.extractData(netlistFilename_, line, resolveParams, modelBinning,scale);

            int numParams = device_.getNumberOfInstanceParameters();
            for( int j = 0; j < numParams; ++j )
            {
              Device::Param param = device_.getInstanceParameter(j);
              if( param.uTag() == "L" )
              {
                MIs[i].inductors[ES1] = param.stringValue();

                // store terminal names associated with this inductor
                std::vector<std::string>::const_iterator paramIter=device_.getNodeValues().begin();
                ( MIs[i].terminals[ES1] ).push_back( *paramIter );
                ++paramIter;
                ( MIs[i].terminals[ES1] ).push_back( *paramIter );
              }
              if( param.uTag() == "IC" )
              {
                ( MIs[i].otherParams[ES1] ).push_back( param );
              }
            }
          }
        }
        return true;
      }
      else if (lineType == '.')
      {
        // Jump to the end of the subcircuit; parseMutualInductances will
        // handle subckt directly
        if (ES1 == ".SUBCKT")
        {
          // Find the subcircuit corresponding to this instance.
          CircuitBlock* subcircuitPtr = 
            findSubcircuit(ExtendedString(line[1].string_).toUpper());

          // Set the end location of the subcircuit in its associated file.
          subcircuitPtr->setFilePosition(subcircuitPtr->getEndPosition());
          subcircuitPtr->setLinePosition( subcircuitPtr->getLineEndPosition() );

        }
        else if (ES1 == ".INCLUDE" || ES1 == ".INCL" || ES1 == ".INC")
        {
          // HSPICE documents .INC, .INCL and .INCLUDE as being a valid .INC line
          std::string includeFile;
          std::string libInside, libSelect;
          Xyce::IO::handleIncludeLine( netlistFilename_, line,
                                       ES1, includeFile, libSelect, libInside );
          if (includeFile != "")
          {
            // Save current ssfPtr_ and netlistFilename_.
            SpiceSeparatedFieldTool* oldssfPtr = ssfPtr_;

            // Save the old file name (parent file)
            std::string old_netlistFilename(netlistFilename_);

            // set the current netlist file to the name of this include file.
            netlistFilename_ = includeFile;

            // If this include file was found once, it will be found again. 
            ssfPtr_ = ssfMap_[includeFile].second;
            ssfPtr_->setLocation(0);
            ssfPtr_->setLineNumber(1);
 
            // Process all the lines in this include file.
            while( getLinePassMI() ) {}

            // Restore old ssfPtr_ and netlistFilename_.
            ssfPtr_ = oldssfPtr;

            // restore the parent file name
            netlistFilename_ = old_netlistFilename;
          }
        }
        else if (ES1 == ".ENDS" || ES1 == ".END")
          return parentCircuitPtr_ == 0;
      }
      else
      {
        return true;
      }
    }
  }

  // Only get here if end of file.
  return false;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::setFileName
// Purpose       : Change netlist file name 
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/2006
//--------------------------------------------------------------------------
void CircuitBlock::setFileName ( const std::string & fileNameIn )
{
  netlistFilename_ = fileNameIn;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::findSubcircuit
// Purpose       : Search the circuitBlockTable_ of the current circuit block
//                 for the subcircuit of the given name. If it is not found,
//                 recursively search each parent subcircuit. Return a
//                 pointer to the circuit block if it is found, otherwise
//                 return NULL.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/28/2001
//--------------------------------------------------------------------------
CircuitBlock* CircuitBlock::findSubcircuit( std::string const& subcircuitName)
{
  // Search this circuit blocks subcircuit list.
  if ( circuitBlockTable_.find( subcircuitName ) != circuitBlockTable_.end() )
  {
      return circuitBlockTable_.find( subcircuitName )->second;
  }
  else
  {
    // The subcircuit was not found in the current circuit's subcircuit list,
    // recursively search the parent circuit's subcircuit list.
    CircuitBlock* circuitBlockPtr = NULL;
    if ( parentCircuitPtr_ != NULL )
    {
      return circuitBlockPtr = parentCircuitPtr_->findSubcircuit( subcircuitName );
    }
    else
    {
      return NULL;
    }
  }
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::resolveExpressionsInOptionBlocks
// Purpose       : Scans the OptionBlock data held in the class var
//                 std::list<Util::OptionBlock> optionsTable for expressions
//                 and tries to resolve them so that a call on an
//                 expressions eval() method will return the right result.
//
// Special Notes : ERK: Note, this function was 100% empty until 9/23/2015, 
//                 when I used it to resolve bug 592 on joseki (Make .IC 
//                 and .NODESET be able to use expressions).
//
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 02/02/2012
//--------------------------------------------------------------------------
bool CircuitBlock::resolveExpressionsInOptionBlocks()
{
  std::list<Util::OptionBlock>::iterator  iter = optionsTable_.begin();
  std::list<Util::OptionBlock>::iterator  end = optionsTable_.end();

  for ( ; iter != end; ++iter)
  {
    // ERK.  for now, only enable this for IC and NODESET.  
    // Fix to use a find command.
    if (iter->getName()=="IC" || iter->getName()=="NODESET" ||
        iter->getName()=="DOT_MEASURE_LINE")
    {
      Util::ParamList::iterator iterPar = iter->begin();
      Util::ParamList::iterator endPar = iter->end();
      for (; iterPar != endPar; ++iterPar)
      {
        // If this is an IC or NODESET from a subcircuit, we can't resolve
        // expressions now.
        if (iterPar->tag() == "SUBCKT")
        {
          break;
        }
        
        circuitContext_.resolveParameter((*iterPar));
      }
    } 
  }

  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::updateAliasNodeMapHelper
// Purpose       : This functions updates the aliasNodeMapHelper_ map to
//                 include "subcircuit interface node names" that were embedded
//                 within expressions.  It also removes those expressions 
//                 from the aliasNodeMapHelper_ map.
// Special Notes : 
// Creator       : Pete Sholander, Electrical Models and Simulation
// Creation Date : 05/10/2018
//--------------------------------------------------------------------------
void CircuitBlock::updateAliasNodeMapHelper()
{
  // Needed for call to circuitContext_.resolveStrings() but this value is
  // then otherwise unused within this function.
  std::vector<std::string> exceptionStrings; 

  // Find the entries in aliasNodeMapHelper_ that are expressions (e.g.,
  // start with '}' and store them in a vector (expStrings)
  std::vector<std::string> expStrings; 
  unordered_set< std::string >::iterator itA=aliasNodeMapHelper_.begin();
  unordered_set< std::string >::iterator endA=aliasNodeMapHelper_.end();
  for ( ; itA != endA; ++itA)
  {
    if ( (*itA)[0] == '{' ) expStrings.push_back((*itA)); 
  }

  // Iterate through the contents of expStrings and add any subcircuit interface
  // node names found to aliasNodeMapHelper_.  Also remove the expression 
  // from aliasNodeMapHelper_ 
  std::vector<std::string>::iterator it=expStrings.begin();
  std::vector<std::string>::iterator end=expStrings.end();
  for ( ; it != end; ++it)
  {
    // Turn the string in expStrings into an expression, and parse it
    std::string expressionString;
    expressionString = (*it).substr(1, (*it).size()-2);
    Util::Expression expression(expressionGroup_,expressionString);

    if (expression.parsed())
    {
      // Resolve the strings in the expression. Unresolved strings
      // may be parameters defined in a .PARAM statement or global
      // parameters defined in .GLOBAL_PARAM statement or may
      // be due to function arguments if the expression is the
      // body of a function defined in a .FUNC statement.
      bool stringsResolved = circuitContext_.resolveStrings(expression, exceptionStrings);

      // Resolve functions in the expression.
      bool functionsResolved = circuitContext_.resolveFunctions(expression);

      // resolve variables in the function body
    
      const std::vector<std::string> & strings = expression.getUnresolvedParams();
      if ( !(strings.empty()) )
      //if ( expression.get_num(XEXP_STRING) > 0 )
      {
        circuitContext_.resolveStrings(expression, exceptionStrings);
      }

      if (stringsResolved && functionsResolved)
      {
        // Check the expression for nodes, and add any that start with X (which are 
        // subcircuit interface node names) to aliasNodeMapHelper_
        std::vector<std::string> nodes;
        expression.getVoltageNodes(nodes);
        for (std::vector<std::string>::iterator node_it = nodes.begin(), end = nodes.end(); node_it != end; ++node_it)
        {
          if ((*node_it)[0] == 'X')
             aliasNodeMapHelper_.insert((*node_it));
        }
      }
    }

    // Remove the expression from aliasNodeMapHelper_, since it's not useful in
    // that form.
    aliasNodeMapHelper_.erase((*it));
  }
        
  return;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::parseIncludeFile
// Purpose       : Parse each include file in includeFiles_ adding the
//                 contents to the current CircuitBlock.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 01/10/2001
//--------------------------------------------------------------------------
bool CircuitBlock::parseIncludeFile(
  PkgOptionsMgr &               options_manager,
  const std::string &           includeFile,
  const std::string &           libSelect,
  std::map<std::string, int> &  fun,
  ModelMap &                    modMap)
{
  // Save current ssfPtr_ and netlistFilename_.
  SpiceSeparatedFieldTool* oldssfPtr = ssfPtr_;

  // save the old file name (parent file)
  std::string old_netlistFilename(netlistFilename_);

  // set the current netlist file to the name of this include file.
  netlistFilename_ = includeFile;

  // get the location in the file just in case we are re-entering the file (only with .lib)
  int oldLineNumber = ssfPtr_->getLineNumber();
  std::streampos oldFilePos = ssfPtr_->getFilePosition();

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlock::parseIncludeFile: Parsing include file: " << includeFile << std::endl;

  if( !ssfMap_.count(includeFile) )
  {
    // Create a new SpiceSeparatedFieldTool for this include file.
    std::ifstream * includeIn = new std::ifstream;

    // Error out if the user-specified include file does not exist, cannot be opened,
    // or is a directory name rather than a file name.  See SON Bugs 730 
    // and 785 for more details.
    if ( !(Util::checkIfValidFile(includeFile)) )
    {
      Report::UserError0() << "Could not find include file " << includeFile;
      return false;
    }
     
    // Using binary to avoid issues with compiler/plat
    // *fstream differences in implementation
    includeIn->open( includeFile.c_str(), std::ios::in | std::ios::binary );
    if ( !includeIn->is_open() )
    {
      Report::UserError0() << "Could not open include file " << includeFile;
      return false;
    }

    ssfPtr_ = new SpiceSeparatedFieldTool(*includeIn, includeFile, externalNetlistParams_);

    ssfMap_[includeFile] = FileSSFPair( includeIn, ssfPtr_ );
  }
  else
  {
    // we already have an ssF for this file in the map, just pull it out,
    // rewind, and proceed.
    ssfPtr_ = ssfMap_[includeFile].second;

    if (DEBUG_IO) {
      Xyce::dout() << "  CircuitBlock::parseIncludeFile: found eisting ssFT " << std::endl
                   << " \t its file name is " << ssfPtr_->getFileName() << std::endl
                   << "\t its current location is " << ssfPtr_->getFilePosition() << std::endl
                   << "\t its current line number is " << ssfPtr_->getLineNumber() << std::endl
                   << "\t Rewinding to location 0 and line 1" << std::endl;
    }

    ssfPtr_->setLocation(0);
    ssfPtr_->setLineNumber(1);
  }

  // Handle the include file lines.
  bool result = false;
  std::vector< std::string > libInside;
  for (;;) {
    bool line_parsed = true;

    if (handleLinePass1(line_parsed, options_manager, fun, modMap, libSelect, libInside) )
      result = result && line_parsed;
    else
      break;
  }

  // Restore old ssfPtr_ and netlistFilename_.
  ssfPtr_ = oldssfPtr;

  // restore the parent file name
  netlistFilename_ = old_netlistFilename;

  // get the location in the file just in case we are re-entering the file (only with .lib)
  ssfPtr_->setLocation(oldFilePos);
  ssfPtr_->setLineNumber(oldLineNumber);

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlock::parseIncludeFile: finished with include file: " << includeFile << std::endl;

  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::handleInitCond
// Purpose       : Retrieve separate IC= data from line or external file and
//               : temporarily store in CircuitBlock
// Special Notes : Validation of .initcond lines is not done here.  Semantic
//               : errors are handled during device instantiation.
// Creator       :
// Creation Date :
//--------------------------------------------------------------------------
void CircuitBlock::handleInitCond(TokenVector const& parsedLine )
{
  // check for multiple .initcond lines
  if( !(mainCircuitPtr_->initCondIndex.empty()) )
  {
    Report::UserError0() << ".INITCOND line may appear only once.";
  }

  // check minimum line length
  if( parsedLine.size() < 3 )
  {
    Report::UserError0().at(netlistFilename_, parsedLine[0].lineNumber_ )
      << ".INITCOND line is missing information";
  }

  ExtendedString tmpType ( parsedLine[1].string_ );
  tmpType.toUpper();

  // read IC values from separate file:  .INITCOND FILE fileName|"filename"
  if( tmpType == "FILE" )
  {
    // Strip off the enclosing double quotes if they are present.
    std::string initCondFile(parsedLine[2].string_);
    if ( ( initCondFile[0] == '"' ) &&
     ( initCondFile[initCondFile.length() - 1] == '"' ) )
    {
      initCondFile = initCondFile.substr( 1, initCondFile.length() - 2 );
    }

    // open the file for reading
    std::ifstream initCondIn;
    initCondIn.open( initCondFile.c_str(), std::ios::in | std::ios::binary );
    if( !initCondIn.is_open() )
    {
      Report::UserError0() << "Could not open the .INITCOND file " << initCondFile;
      return;
    }

    // use parser to extract data from the file
    SpiceSeparatedFieldTool ssfICPtr( initCondIn, initCondFile, externalNetlistParams_ );
    TokenVector line;

    while( !initCondIn.eof() )
    {
      // tokenize lines; file ptr is advanced in getLine()
      ssfICPtr.getLine( line, preprocessFilter_[PreprocessType::REPLACE_GROUND] );

      // check for enough data on line
      if( line.size() < 4 )
      {
        Report::UserError0() << ".INITCOND file '" << initCondFile << "' is not formatted properly.";
      }
      else {

        // store tokenized line in the index
        tmpType = line[0].string_;
        tmpType.toUpper();
        mainCircuitPtr_->initCondIndex[tmpType] =
          TokenVector(
            line.begin() + 1, line.end() );
      }
    }
  }

  // read IC values from line:  .INITCOND ( fqDevName IC = val (, val)* )+
  else
  {
    TokenVector::const_iterator
      pIter, pNextIter, pEndIter;

    pIter = parsedLine.begin() + 1;
    pEndIter = parsedLine.end();

    // check for enough data on line
    if( distance( pIter, pEndIter ) < 2 )
    {
      Report::UserError0() << ".INITCOND line is not formatted properly.";
    }
    else {

      // find the next dev name
      while( pIter != pEndIter )
      {
        // point to beginning of IC=val1...valN list
        pNextIter = pIter + 3;

        // build var1..varN list
        while( pNextIter != pEndIter && (*pNextIter).string_ != "=")
        {
          ++pNextIter;
        }

        // keep checking if more data is on the line
        if( pNextIter != pEndIter )
        {
          // move back to end of list
          pNextIter -= 2;
        }

        // copy into map
        ExtendedString tmpType ( (*pIter).string_ );
        tmpType.toUpper();
        mainCircuitPtr_->initCondIndex[tmpType] =
          TokenVector(
            pIter + 1, pNextIter );

        // move to end of line
        pIter = pNextIter;
      }
    }
  }

  if (DEBUG_IO) {
    std::map< std::string, TokenVector >::const_iterator a, b;
    a=mainCircuitPtr_->initCondIndex.begin();
    b=mainCircuitPtr_->initCondIndex.end();
    Xyce::dout() << ".INITCOND line yields " << mainCircuitPtr_->initCondIndex.size() <<
      " parsed devices:  ";
    for(; a!=b;++a)
    {
      Xyce::dout() << (*a).first << " ";
      for(unsigned int i=0; i< mainCircuitPtr_->initCondIndex[(*a).first].size();++i)
        Xyce::dout() << ( (mainCircuitPtr_->initCondIndex[(*a).first])[i] ).string_;
      Xyce::dout() << " ";
    }
    Xyce::dout() << std::endl;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getICNodesetList
// Purpose        : Collect IC and NODESET statements from used subcircuits
// Special Notes  : This requires traversing the hierarchy to get the resolved node names
// Scope          :
// Creator        : Heidi K. Thornquist
// Creation Date  : 08/03/2016
//----------------------------------------------------------------------------
void CircuitBlock::getICNodesetList( std::vector<std::string>& usedSubcircuits, 
                                     std::vector<Util::OptionBlock>& icNodesetList )
{
   // Add any IC or NODESET statements from this subcircuit
   if (std::binary_search( usedSubcircuits.begin(), usedSubcircuits.end(), name_ ))
   {
     std::list<Util::OptionBlock>::iterator optionIter = optionsTable_.begin();
     for ( ; optionIter != optionsTable_.end(); ++optionIter )
     {
       // Only get the IC or NODESET option blocks.
       const std::string& optionName = optionIter->getName();
       if ( (optionName == "IC") || (optionName == "NODESET") )
       {
         icNodesetList.push_back( *optionIter );
       }
     }
   }

   // Check subcircuits of this subcircuit, recursion.
   unordered_map< std::string, CircuitBlock * >::iterator itcbt = circuitBlockTable_.begin();
   for ( ; itcbt != circuitBlockTable_.end(); ++itcbt )
   {
     itcbt->second->getICNodesetList( usedSubcircuits, icNodesetList );
   }
}



namespace {

struct GetNameEqual
{
  bool operator()(const Util::OptionBlock &op, const char *name)
  {
    return op.getName() == name;
  }
};

} // namespace <unnamed>

//----------------------------------------------------------------------------
// Function       : CircuitBlock::handleAnalysis
// Purpose        : Post process analysis statements
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/06/2002
//----------------------------------------------------------------------------
bool CircuitBlock::handleAnalysis()
{
  static const char *analysisOptions_[] = {"DC", "TRAN", "TR", "MPDE", "HB", "AC", "OP", "NOISE", "ROL"}; // TT
  static const char *printOptions_[] = {"PRINT"};

  std::list<Util::OptionBlock>::const_iterator op_analysis_it
    = std::find_first_of(optionsTable_.begin(), optionsTable_.end(),
                         &analysisOptions_[0], &analysisOptions_[sizeof(analysisOptions_)/sizeof(analysisOptions_[0])], GetNameEqual());

  // Check if Model Order Reduction (MOR) is requested by the user.
  std::list<Util::OptionBlock>::const_iterator op_analysis_it2 = optionsTable_.begin();
  for ( ; op_analysis_it2 != optionsTable_.end(); op_analysis_it2++)
  {
    if (op_analysis_it2->getName() == "MOR")
    {
      morFlag_ = true;
      break;
    }
  }

  if (op_analysis_it == optionsTable_.end())
  {
    if (morFlag_)
    {
      setAnalysisName( "MOR" );
    }
    else
    {
      // Problem, no analysis specified.
      Report::UserError0() << "No analysis specified.";
      return false;
    }
  }
  else
  {
    setAnalysisName( op_analysis_it->getName() );
  }

  // Add the DC sweep parameters to the PRINT options. First, find the
  // DC PRINT options.
  std::list<Util::OptionBlock>::const_iterator op_param_it
    = std::find_first_of(optionsTable_.begin(), optionsTable_.end(),
                         &printOptions_[0], &printOptions_[sizeof(printOptions_)/sizeof(printOptions_[0])], GetNameEqual());
  if (op_param_it == optionsTable_.end())
  {
    Report::UserWarning0() << "No print specified";
    return true;
  }

  Util::ParamList::const_iterator paramIter;
  paramIter = std::find_if(op_param_it->begin(), op_param_it->end(), Util::EqualParam("TYPE"));

  // Check for consistency between analysis type and print type.
  std::string usVal = paramIter->usVal();
  if (analysisName_ == "TR")
  {
    analysisName_ = "TRAN"; // TR is a synonym for TRAN
  }
  if (usVal == "TR")
    usVal = "TRAN";

  if (!( (analysisName_ == "TRAN" && usVal == "TRAN") ||
         (analysisName_ == "OP" && usVal == "TRAN") ||
         (analysisName_ == "OP" && usVal == "AC") ||
         (analysisName_ == "OP" && usVal == "AC_IC") ||
         (analysisName_ == "OP" && usVal == "DC") ||
         (analysisName_ == "OP" && usVal == "ES") ||
         (analysisName_ == "OP" && usVal == "HB") ||
         (analysisName_ == "OP" && usVal == "HB_TD") ||
         (analysisName_ == "OP" && usVal == "HB_FD") ||
         (analysisName_ == "OP" && usVal == "HB_IC") ||
         (analysisName_ == "OP" && usVal == "HB_STARTUP") ||
         (analysisName_ == "OP" && usVal == "NOISE") ||
         (analysisName_ == "OP" && usVal == "PCE") ||
         (analysisName_ == "TRAN" && usVal == "HOMOTOPY") ||
         (analysisName_ == "OP" && usVal == "HOMOTOPY") ||
         (analysisName_ == "DC" && usVal == "HOMOTOPY") ||
         (analysisName_ == "DC" && usVal == "DC") ||
         (analysisName_ == "DC" && usVal == "ES") ||
         (analysisName_ == "DC" && usVal == "PCE") ||
         (analysisName_ == "OP" && usVal == "SENS") ||
         (analysisName_ == "DC" && usVal == "SENS") ||
         (analysisName_ == "TRAN" && usVal == "ES") ||
         (analysisName_ == "TRAN" && usVal == "PCE") ||
         (analysisName_ == "TRAN" && usVal == "SENS") ||
         (analysisName_ == "MPDE" && usVal == "TRAN") ||
         (analysisName_ == "MPDE" && usVal == "MPDE") ||
         (analysisName_ == "MPDE" && usVal == "MPDE_IC") ||
         (analysisName_ == "MPDE" && usVal == "MPDE_STARTUP") ||
         (analysisName_ == "HB" && usVal == "HB") ||
         (analysisName_ == "HB" && usVal == "HB_TD") ||
         (analysisName_ == "HB" && usVal == "HB_FD") ||
         (analysisName_ == "HB" && usVal == "HB_IC") ||
         (analysisName_ == "HB" && usVal == "HB_STARTUP") ||
         (analysisName_ == "AC" && usVal == "AC") ||
         (analysisName_ == "AC" && usVal == "AC_IC") ||
         (analysisName_ == "AC" && usVal == "SENS") ||
         (analysisName_ == "AC" && usVal == "SPARAM") ||
         (analysisName_ == "NOISE" && usVal == "NOISE") ||
         (analysisName_ == "MOR" && usVal == "MOR")  ||
         (analysisName_ == "ROL" && usVal == "DC"))) // TT
  {
    // Problem, inconsistent analysis type and print type.
    Report::UserError0() << "Analysis type " << analysisName_ << " and print type " << usVal << " are inconsistent.";
    return false;
  }
  
  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::parsePreprocess
// Purpose       : This function introduces a preprocessing phase whereby
//                 it is determined whether Xyce should remove "redundant"
//                 devices (devices for which all of the nodes are the
//                 same).  This has to be caught before the initial parse
//                 to determine whether to add a device to the circuit or
//                 not.
//
// Special Notes : Updated 12/6/07 for additional detection:  we detect flags
//                 to create netlist files which contain resistors to ground
//                 for nodes which are connected only to one device terminal,
//                 and/or resistors to ground for nodes that have no DC path
//                 to ground.
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/05/2007
//--------------------------------------------------------------------------
bool CircuitBlock::parsePreprocess()
{
  if (DEBUG_IO && isActive(Diag::IO_PARSE))
    Xyce::dout() << "Preprocess of netlist file: " << netlistFilename_ << std::endl;

  //Get the first character of input.
  TokenVector line;
  char lineType;
  int eof = ssfPtr_->peekAtNextLine( lineType );
  int removecounter = 0;
  int replacecounter = 0;
  int onetermcounter = 0;
  int nodcpathcounter = 0;

  while (!eof)
  {
    if (lineType == '.')
    {
      eof = !ssfPtr_->getLine(line); //Breaks the line into fields.
      ExtendedString ES1 ( line[0].string_ );
      ES1.toUpper();

      if ( ES1 != ".PREPROCESS" )
      {
        //do nothing
      }
      else if (line.size() < 3)
      {
        Report::UserError().at(netlistFilename_, line[0].lineNumber_)
          << "Too few parameters specified in .PREPROCESS statement.";
      }
      else
      {
        ExtendedString preprocarg ( line[1].string_ );
        preprocarg.toUpper();

        if (preprocarg == "REMOVEUNUSED")
        {
          if (removecounter != 0)
          {
            Report::UserError().at(netlistFilename_, line[0].lineNumber_)
              << "Multiple .PREPROCESS REMOVEUNUSED statements.";
          }
          else
          {
            removecounter++;
            ExtendedString removeparam ( ES1 );
            bool anyparamsremoved = false;

            for (unsigned int i = 2; i < line.size(); ++i)
            {
              removeparam=line[i].string_;
              removeparam.toUpper();
              if (removeparam == "C")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_C] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "D")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_D] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "I")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_I] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "L")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_L] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "M")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_M] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "Q")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_Q] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "R")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_R] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "V")
              {
                preprocessFilter_[PreprocessType::REDUNDANT_V] = true;
                anyparamsremoved = true;
              }
              else if (removeparam == ",")
              {
                //skip commas
              }
              else
              {
                Report::UserError().at(netlistFilename_, line[0].lineNumber_)
                  << "Unknown argument type " << removeparam << " in .PREPROCESS REMOVEUNUSED statement.";
              }
            }
            if (anyparamsremoved)
            {
              remove_any_redundant_ = true;
            }
            else
            {
              //didn't find any parameters on the line
              Report::UserError().at(netlistFilename_, line[0].lineNumber_)
                << "No remove parameters specified in .PREPROCESS REMOVEUNUSED statement.";
            }
          }
        }
        else if (preprocarg == "REPLACEGROUND")
        {
          if (replacecounter != 0)
          {
            Report::UserError().at(netlistFilename_, line[0].lineNumber_)
              << "Multiple .PREPROCESS REPLACEGROUND statements.";
          }
          else
          {
            replacecounter++;
            if (line.size() > 3)
            {
              Report::UserWarning().at(netlistFilename_, line[0].lineNumber_)
                << "Additional parameters in .PREPROCESS REPLACEGROUND statement.  Ignoring.";
            }

            ExtendedString replaceparam(line[2].string_);
            replaceparam.toUpper();

            if (replaceparam == "TRUE")
            {
              preprocessFilter_[PreprocessType::REPLACE_GROUND]=true;
            }
            else if (replaceparam != "FALSE")
            {
              Report::UserError().at(netlistFilename_, line[0].lineNumber_)
                << "Unknown argument " << replaceparam << " in .PREPROCESS REPLACEGROUND statement.";
            }
          }
        }
        else if (preprocarg == "ADDRESISTORS")
        {
          hangingResistor_.setNetlistCopy(true);

          if (line.size() > 4)
          {
            Report::UserWarning().at(netlistFilename_, line[0].lineNumber_)
              << "Additional parameters in .PREPROCESS ADDRESISTORS statement.  Ignoring.";
          }
          else if (line.size() < 4)
          {
            Report::UserError().at(netlistFilename_, line[0].lineNumber_)
              << "Missing resistance value in .PREPROCESS ADDRESISTORS statement.";
          }

          ExtendedString netlistparam(line[2].string_);
          netlistparam.toUpper();
          ExtendedString resistanceparam(line[3].string_);
          resistanceparam.toUpper();

          if (netlistparam == "ONETERMINAL")
          {
            if (onetermcounter != 0)
            {
              Report::UserError().at(netlistFilename_, line[0].lineNumber_)
                << "Multiple .PREPROCESS ADDRESISTORS ONETERMINAL statements.";
            }
            else
            {
              onetermcounter++;
              hangingResistor_.setOneTerm(true);
              hangingResistor_.setOneTermRes(resistanceparam);
            }
          }
          else if (netlistparam == "NODCPATH")
          {

            if (nodcpathcounter != 0)
            {
              Report::UserError().at(netlistFilename_, line[0].lineNumber_)
                << "Multiple .PREPROCESS ADDRESISTORS NODCPATH statements.";
            }
            else
            {
              nodcpathcounter++;
              hangingResistor_.setNoDCPath(true);
              hangingResistor_.setNoDCPathRes(resistanceparam);
            }
          }
          else
          {
            Report::UserError().at(netlistFilename_, line[0].lineNumber_)
              <<  "Unknown argument " << netlistparam << " in .PREPROCESS ADDRESISTORS statement.";
          }
        }
        else
        {
          Report::UserError().at(netlistFilename_, line[0].lineNumber_)
            << "Unknown keyword " <<  preprocarg << " specified in .PREPROCESS statement.";
        }
      }
    }
    else
    {
      // Skip to the end.
      ssfPtr_->skipToEndOfLine();
    }
    eof = ssfPtr_->peekAtNextLine( lineType );
  }

  if (DEBUG_IO) {
    Xyce::dout() << std::endl << "Unused components to be removed.  (1 means remove "
                 << "redundancies,"  << std::endl <<" 0 means do not remove redundancies): "
                 << std::endl << std::endl

                 << "Remove Unused C:  " << preprocessFilter_[PreprocessType::REDUNDANT_C] << std::endl
                 << "Remove Unused D:  " << preprocessFilter_[PreprocessType::REDUNDANT_D] << std::endl
                 << "Remove Unused I:  " << preprocessFilter_[PreprocessType::REDUNDANT_I] << std::endl
                 << "Remove Unused L:  " << preprocessFilter_[PreprocessType::REDUNDANT_L] << std::endl
                 << "Remove Unused M:  " << preprocessFilter_[PreprocessType::REDUNDANT_M] << std::endl
                 << "Remove Unused Q:  " << preprocessFilter_[PreprocessType::REDUNDANT_Q] << std::endl
                 << "Remove Unused R:  " << preprocessFilter_[PreprocessType::REDUNDANT_R] << std::endl
                 << "Remove Unused V:  " << preprocessFilter_[PreprocessType::REDUNDANT_V] << std::endl << std::endl

                 << "Replace Ground Flag set to:  " << preprocessFilter_[PreprocessType::REPLACE_GROUND] << std::endl << std::endl
                 << "Netlist copy Flag set to:  "
                 << hangingResistor_.getNetlistCopy() << std::endl << std::endl
                 << "One terminal Flag set to:  "
                 << hangingResistor_.getOneTerm() << std::endl << std::endl
                 << "No DC Path Flag set to:  "
                 << hangingResistor_.getNoDCPath() << std::endl << std::endl
                 << "One terminal resistance:  "
                 << hangingResistor_.getOneTermRes() << std::endl << std::endl
                 << "No DC path resistance:  "
                 << hangingResistor_.getNoDCPathRes() << std::endl << std::endl

                 << "Done with preprocess netlist file parsing." << std::endl << std::endl;
  }

  return true;
}


//--------------------------------------------------------------------------
// Function      : CircuitBlock::produceUnflattenedNetlist
// Purpose       : Generates a copy of the current netlist in the file
//                 netlistFilename_copy.cir.  This is used to create netlist
//                 files that contain resistors that connect ground to
//                 "dangling" nodes (nodes which either don't have a dc path
//                 to ground, or which are only connected to one device
//                 terminal).
//
// Special Notes :
//
// Creator       : Keith Santarelli, Electrical and Microsystems Modeling
//
// Creation Date : 12/5/07
//--------------------------------------------------------------------------
void CircuitBlock::produceUnflattenedNetlist()
{

  if (DEBUG_IO)
    Xyce::dout() << "Producing copy of netlist file: " << netlistFilename_ << std::endl;

  // Reset the SSF pointer to the beginning of the file.
  ssfPtr_->setLocation(fileStartPosition_);
  ssfPtr_->setLineNumber( lineStartPosition_ );
  netlistIn_->clear();
  netlistIn_->seekg(0, std::ios::beg);
  ssfPtr_->changeCursorLineNumber( 1 );

  // Create the output file stream
  std::string netlistCopy(netlistFilename_);
  netlistCopy += "_xyce.cir";
  std::ofstream copyFile(netlistCopy.c_str());

  // Some error checking in case we can't open the file.
  if(copyFile.fail())
  {
    Report::UserError0() << ".PREPROCESS NETLISTCOPY cannot open output file " << netlistCopy;
    return;
  }


  //Create date/time stamp
  const time_t now = time( NULL );
  char timeDate[ 40 ];
  strftime( timeDate, 40, "TIME='%I:%M:%S %p' DATE='%b %d, %Y' ",
    localtime( &now ) );

  //Create title line (not the same as title line of original netlist file!)
  std::string header("XYCE-generated Netlist file copy:  ");
  header += timeDate;
  copyFile << header << std::endl;

  //Add the original title:
  copyFile << "*Original Netlist Title:  " << std::endl << std::endl;
  copyFile << "*";

  TokenVector separatedLine;
  int eof = !ssfPtr_->getLineWithComments(separatedLine);

  //Add the original title text:
  for (unsigned int i = 0; i < separatedLine.size(); ++i)
  {
    copyFile << separatedLine[i].string_;
  }

  copyFile << std::endl;

  eof = !ssfPtr_->getLineWithComments(separatedLine);

  bool endflag=false;
  bool addresistbool=false;
  ExtendedString firstarg("");
  ExtendedString addresistarg("");

  if ( !(separatedLine.empty()) )
  {
    firstarg=separatedLine[0].string_;
    firstarg.toUpper();
  }

  while (!eof && !endflag)
  {
    if (firstarg != ".PREPROCESS")
    {
      for (unsigned int i = 0; i < separatedLine.size(); ++i)
      {
        copyFile << separatedLine[i].string_;
      }
    }
    else
    {
      for (unsigned int i = 0; i < separatedLine.size(); ++i)
      {
        addresistarg = separatedLine[i].string_;
        addresistarg.toUpper();
        if (addresistarg == "ADDRESISTORS")
        {
          addresistbool = true;
        }
      }

      if (!addresistbool)
      {
        for (unsigned int i = 0; i < separatedLine.size(); ++i)
        {
          copyFile << separatedLine[i].string_;
        }
      }
      else
      {
        copyFile << "*";
        for (unsigned int i = 0; i < separatedLine.size(); ++i)
        {
          copyFile << separatedLine[i].string_;
        }
        copyFile << "*Xyce:  \".PREPROCESS ADDRESISTORS\" statement";
        copyFile << " automatically commented out in netlist copy.";
        copyFile << std::endl;
      }
    }

    //Get the next line of input
    eof = !ssfPtr_->getLineWithComments(separatedLine);

    if ( !(separatedLine.empty()) )
    {
      firstarg = separatedLine[0].string_;
    }
    else
    {
      firstarg="";
    }
    firstarg.toUpper();

    //We don't reproduce anything after a .END statement!
    if ( !(separatedLine.empty()) && firstarg == ".END")
    {
      endflag = true;
    }
  }

 copyFile.close();
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::resetTopLevelSSFPtr
// Purpose       : Reset the file and line positions and at the top level,
//               : skip the comment line
// Special Notes :
//
// Creator       : Heidi Thornquist, SNL
//
// Creation Date : 07/30/2015
//--------------------------------------------------------------------------
void CircuitBlock::resetSSFPtr()
{
  ssfPtr_->setLocation(fileStartPosition_);
  ssfPtr_->setLineNumber( lineStartPosition_ );

  // If this is the main circuit, skip over the title line and continue.
  if ( parentCircuitPtr_ == NULL )
  {
    std::string title("");
    netlistIn_->clear();
    netlistIn_->seekg(0, std::ios::beg);
    Xyce::IO::readLine( *netlistIn_, title );
    ssfPtr_->changeCursorLineNumber( 1 );
  }
}


} // namespace IO
} // namespace Xyce

//-----------------------------------------------------------------------------
// Function      : packCircuitOptions
// Purpose       : send option blocks to all procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Xyce::IO::packCircuitOptions(const std::list<Util::OptionBlock>& options, char* char_buffer,
                                 int char_buffer_size, Xyce::Parallel::Communicator* pds_comm_ptr)
{
  int bsize = 0;

  if (Parallel::is_parallel_run(pds_comm_ptr->comm()))
  {
    int pos = 0;

    // pack options
    int count = options.size();
    pds_comm_ptr->pack( &count, 1, char_buffer, char_buffer_size, pos );
    for (std::list<Util::OptionBlock>::const_iterator it = options.begin(), end = options.end(); it != end; ++it)
    {
      Pack<Util::OptionBlock>::pack(*it, char_buffer, char_buffer_size, pos, pds_comm_ptr);
    }

    bsize = pos;
  }

  return bsize;
}

//-----------------------------------------------------------------------------
// Function      : unpackCircuitOptions
// Purpose       : unpack option blocks from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Xyce::IO::unpackCircuitOptions(std::list<Util::OptionBlock>& options, char* char_buffer,
                                    int bsize, Xyce::Parallel::Communicator* pds_comm_ptr)
{
  if (Parallel::is_parallel_run(pds_comm_ptr->comm()))
  {
    int pos = 0;
    int size = 0;

    // unpack options
    pds_comm_ptr->unpack( char_buffer, bsize, pos, &size, 1 );
    for( int i = 0; i < size; ++i )
    {
      Util::OptionBlock anOptionBlock;
      Pack<Util::OptionBlock>::unpack(anOptionBlock, char_buffer, bsize, pos, pds_comm_ptr);
      options.push_back( anOptionBlock );
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : packAliasNodeMap
// Purpose       : send alias node map to all procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Xyce::IO::packAliasNodeMap(const AliasNodeMap& alias_node_map, char* char_buffer,
                               int char_buffer_size, Xyce::Parallel::Communicator* pds_comm_ptr)
{
  int bsize = 0;

  if (Parallel::is_parallel_run(pds_comm_ptr->comm()))
  {
    int pos = 0;

    // pack options
    int count = alias_node_map.size();
    pds_comm_ptr->pack( &count, 1, char_buffer, char_buffer_size, pos );
    for (AliasNodeMap::const_iterator it = alias_node_map.begin(), end = alias_node_map.end(); it != end; ++it)
    {
      Xyce::Parallel::PackTraits<std::string>::pack((*it).first, char_buffer, char_buffer_size, pos, *pds_comm_ptr);
      Xyce::Parallel::PackTraits<std::string>::pack((*it).second, char_buffer, char_buffer_size, pos, *pds_comm_ptr);
    }

    bsize = pos;
  }

  return bsize;
}

//-----------------------------------------------------------------------------
// Function      : unpackAliasNodeMap
// Purpose       : unpack alias node map from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Xyce::IO::unpackAliasNodeMap(AliasNodeMap& alias_node_map, char* char_buffer,
                                  int bsize, Xyce::Parallel::Communicator* pds_comm_ptr)
{
  if (Parallel::is_parallel_run(pds_comm_ptr->comm()))
  {
    int pos = 0;
    int size = 0;

    // unpack options
    pds_comm_ptr->unpack( char_buffer, bsize, pos, &size, 1 );
    for( int i = 0; i < size; ++i )
    {
      std::string key, value;
      Xyce::Parallel::PackTraits<std::string>::unpack(key, char_buffer, bsize, pos, *pds_comm_ptr);
      Xyce::Parallel::PackTraits<std::string>::unpack(value, char_buffer, bsize, pos, *pds_comm_ptr);
      alias_node_map[key] = value;
    }
  }

  return true;
}


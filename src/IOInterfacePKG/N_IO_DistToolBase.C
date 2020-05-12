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
// Purpose        : Defines the DistToolBase class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
// Special Notes  :
//
// Creator        : Eric Rankin, SNL
//
// Creation Date  : 03/12/2003
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_IO_fwd.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_ERH_Message.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_DistToolBase.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_ParsingHelpers.h>
#include <N_PDS_Comm.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Stats.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : DistToolBase::DistToolBase
// Purpose       : ctor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistToolBase::DistToolBase(
  Parallel::Communicator                 * pdsCommPtr,
  CircuitBlock                           & circuit_block,
  std::map<std::string,FileSSFPair>      & ssfMap
  )
  : pdsCommPtr_(pdsCommPtr),
    numProcs_(pdsCommPtr->numProc()),
    circuitContextReady_(false),
    circuitOptionsReady_(false),
    checkSubcktICNODESET_(false),
    resolveSubcktICNODESET_(false),
    charBufferSize_(250000),
    charBufferPos_(0),
    charBuffer_(0),
    circuitBlock_(circuit_block),
    circuitContext_(circuit_block.getCircuitContextPtr()),
    ssfMap_(ssfMap),
    options_(circuit_block.getOptionsTable()),
    netlistFilename_(),
    device_(*circuit_block.getCircuitContextPtr(), circuit_block.getMetadata()),
    mainCircuitPtr_(&circuit_block), 
    parentCircuitPtr_(NULL),
    currentCircuitPtr_(&circuit_block),
    preprocessFilter_(PreprocessType::NUM_PREPROCESS, false),
    ssfPtr_(0),
    remove_any_redundant_(false)
{
}

//-----------------------------------------------------------------------------
// Function      : DistToolBase::~DistToolBase
// Purpose       : dtor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistToolBase::~DistToolBase()
{
  delete[] charBuffer_;
}

//-----------------------------------------------------------------------------
// Function      : DistToolBase::broadcastGlobalData
// Purpose       : send bufsize, options, and context to procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
DistToolBase::broadcastGlobalData()
{
  // Get the preprocessor filter before broadcasting.
  preprocessFilter_ = circuitBlock_.getPreprocessFilter();

  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    // tx buffer size
    pdsCommPtr_->bcast( &charBufferSize_, 1, 0 );

    // check for halt request
    if (charBufferSize_ < 0)
    {
      return false;
    }

    // FIRST communicate filename vector to all processors
    const Filename::FilenameVector & files = Filename::getFilenameVector();
   
    // reserve memory; including offset
    charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

    // make sure data is ready to tx
    if (numProcs_ > 1)
    {
      int bsize = 0;
      int pos = 0;

      // pack filenames, broadcast buffer size needed for filenames
      bsize = Xyce::packedByteCount(files);
      pdsCommPtr_->bcast( &bsize, 1, 0 );
      char * filenameBuf = new char[ bsize ];

      // pack buffer
      if (pdsCommPtr_->procID() == 0)
        Pack<Filename::FilenameVector>::pack( files, filenameBuf, bsize, pos, pdsCommPtr_ );

      pdsCommPtr_->bcast( filenameBuf, bsize, 0 ); 

      if (pdsCommPtr_->procID() != 0)
      {
        Filename::FilenameVector & nonconst_files = const_cast<Filename::FilenameVector &>(files);
        Pack<Filename::FilenameVector>::unpack( nonconst_files, filenameBuf, bsize, pos, pdsCommPtr_ );
      }

      // delete buffer
      delete [] filenameBuf;

      // tx global data
      if (circuitOptionsReady_)
      {
        bsize = Xyce::IO::packCircuitOptions(options_, charBuffer_, charBufferSize_, pdsCommPtr_);
      }

      // broadcast options
      pdsCommPtr_->bcast( &bsize, 1, 0 );
      pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

      // receive global data
      if (!circuitOptionsReady_)
      {
        Xyce::IO::unpackCircuitOptions(options_, charBuffer_, bsize, pdsCommPtr_);
      }

      // tx global data
      if (circuitContextReady_)
      {
        bsize = Xyce::IO::packCircuitContext(circuitContext_, charBuffer_, charBufferSize_, pdsCommPtr_);
      }
      // broadcast context
      pdsCommPtr_->bcast( &bsize, 1, 0 );
      pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

      if (!circuitContextReady_)
      {
        Xyce::IO::unpackCircuitContext(circuitContext_, charBuffer_, bsize, pdsCommPtr_);
      }
    }

    HangingResistor &hanging_resistor = circuitBlock_.getHangingResistor();
 
    // transform booleans into integers and transmit them
    bool netlistcopy = hanging_resistor.getNetlistCopy();
    bool oneterm = hanging_resistor.getOneTerm();
    bool nodcpath = hanging_resistor.getNoDCPath();

    int nlcopyint = 0;
    int otint = 0;
    int nodcint = 0;

    if (netlistcopy)
      nlcopyint = 1;
    if (oneterm)
      otint = 1;
    if (nodcpath)
      nodcint = 1;

    // tx global data
    pdsCommPtr_->bcast(&nlcopyint,1,0);
    pdsCommPtr_->bcast(&otint,1,0);
    pdsCommPtr_->bcast(&nodcint,1,0);

    std::string onetermres(hanging_resistor.getOneTermRes());
    std::string nodcpathres(hanging_resistor.getNoDCPathRes());
    int otreslength = onetermres.size();
    int nodcreslength = nodcpathres.size();

    // tx global data
    pdsCommPtr_->bcast(&otreslength,1,0);
    onetermres.resize(otreslength);
    pdsCommPtr_->bcast(&onetermres[0],otreslength,0);
    pdsCommPtr_->bcast(&nodcreslength,1,0);
    nodcpathres.resize(nodcreslength);
    pdsCommPtr_->bcast(&nodcpathres[0],nodcreslength,0);

    if (nlcopyint == 1)
      hanging_resistor.setNetlistCopy(true);
    if (otint == 1)
      hanging_resistor.setOneTerm(true);
    if (nodcint == 1)
      hanging_resistor.setNoDCPath(true);

    hanging_resistor.setOneTermRes(onetermres);
    hanging_resistor.setNoDCPathRes(nodcpathres);

    N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());

    // tx preprocess filter
    for (int i=0; i<PreprocessType::NUM_PREPROCESS; i++)
    {
      int preprocint = 0;
      if ( preprocessFilter_[i] )
        preprocint = 1;

      pdsCommPtr_->bcast(&preprocint,1,0);

      if ( preprocint == 1 )
        preprocessFilter_[i] = true;
    }
  }

  if ( preprocessFilter_[PreprocessType::REDUNDANT_C]
       || preprocessFilter_[PreprocessType::REDUNDANT_D]
       || preprocessFilter_[PreprocessType::REDUNDANT_I]
       || preprocessFilter_[PreprocessType::REDUNDANT_L]
       || preprocessFilter_[PreprocessType::REDUNDANT_M]
       || preprocessFilter_[PreprocessType::REDUNDANT_Q]
       || preprocessFilter_[PreprocessType::REDUNDANT_R]
       || preprocessFilter_[PreprocessType::REDUNDANT_V] )
  {
    remove_any_redundant_ = true;
  }

  // Broadcast analysis type
  int len = mainCircuitPtr_->getAnalysisName().length();
  char* analysisType = 0;
  if (len > 0)
  {
    analysisType = new char[len+1];
    std::strcpy( analysisType, mainCircuitPtr_->getAnalysisName().c_str() );
  }
  pdsCommPtr_->bcast(&len,1,0);
  if (analysisType == 0)
  {
    analysisType = new char[len+1];
  }
  pdsCommPtr_->bcast(analysisType,len+1,0);
  std::string newAnalysisType( analysisType );
  mainCircuitPtr_->setAnalysisName( newAnalysisType );

  delete [] analysisType;

  // Broadcast level set
  if (numProcs_ > 1)
  {
    // Broadcast size of level set
    const std::set<int>& ls = mainCircuitPtr_->getLevelSet();
    int size = ls.size();
    pdsCommPtr_->bcast(&size,1,0);

    if (size)
    {
      // allocate buffer
      std::vector<int> levelSet(size);
   
      // pack buffer
      if (pdsCommPtr_->procID() == 0)
      {
        std::set<int>::const_iterator it_ls = ls.begin();
        std::set<int>::const_iterator it_ls_end = ls.end();
        for (int i=0; it_ls!=it_ls_end; ++i, ++it_ls)
          levelSet[i] = *it_ls;
      }

      // send
      pdsCommPtr_->bcast(&levelSet[0],size,0);

      // unpack buffer
      if (pdsCommPtr_->procID() != 0)
      {
        std::set<int> newls(levelSet.begin(), levelSet.end());
        mainCircuitPtr_->setLevelSet( newls );
      }
    }
  }    

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolBase::check_IC_NODESET_OptionBlock()
// Purpose       : check whether option block has .IC or .NODESET for subcircuits
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolBase::check_IC_NODESET_OptionBlock()
{
  if (!checkSubcktICNODESET_)
  {
    // Determine if there are any .IC or .NODESET statements in subcircuits to resolve
    std::list<Util::OptionBlock>::const_iterator it = options_.begin();
    for ( ; it != options_.end(); ++it )
    {
      if (it->getName()=="IC" || it->getName()=="NODESET")
      {
        if (it->begin()->tag() == "SUBCKT")
        {
          resolveSubcktICNODESET_ = true;
        }
      }
      if (resolveSubcktICNODESET_)
        break;
    }

    // We've now checked
    checkSubcktICNODESET_ = true;
  }

  return resolveSubcktICNODESET_;
}

//-----------------------------------------------------------------------------
// Function      : DistToolBase::setCircuitContext
// Purpose       : check buffer size wrt the circuit context size
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void
DistToolBase::setCircuitContext()
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if (pdsCommPtr_->procID() == 0)
    {
      // resize buffer if necessary
      int tmpSize = Xyce::packedByteCount(*circuitContext_) + sizeof( int );

      if ( tmpSize > charBufferSize_ )
      {
        charBufferSize_ = tmpSize;
      }

      // set flag to tx global data
      circuitContextReady_ = true;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DistToolBase::setCircuitOptions
// Purpose       : stage options for tx
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void
DistToolBase::setCircuitOptions()
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if (pdsCommPtr_->procID() == 0)
    {
      // resize buffer if necessary
      int tmpSize = sizeof(int);
      for (std::list<Util::OptionBlock>::const_iterator it = options_.begin(), end = options_.end(); it != end; ++it)
        tmpSize += Xyce::packedByteCount(*it);

      if ( tmpSize > charBufferSize_ )
        charBufferSize_ = tmpSize;

      // set flag to tx global data
      circuitOptionsReady_ = true;
    }
  }
}

//----------------------------------------------------------------------------
// Function       : DistToolBase::getLine
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 08/01/2003
//----------------------------------------------------------------------------
bool DistToolBase::getLine( TokenVector &line,
  const std::string &libSelect, std::vector<std::string> &libInside)
{
  int eof = 0;

  while (!eof)
  {
    {
      //Stats::StatTop _getLineStat("Tokenize Line");
      //Stats::TimeBlock _getLineTimer(_getLineStat);

      eof = !ssfPtr_->getLine(line,preprocessFilter_[PreprocessType::REPLACE_GROUND]); // Breaks the line into fields.
    }

    if (DEBUG_IO) 
    {
      Xyce::dout() << "Pass 2 read netlist line:  "<< std::endl;
      for (unsigned int i = 0; i < line.size(); ++i)
      {
        Xyce::dout() << line[i].string_ << " ";
      }
      Xyce::dout() << std::endl;

      if (eof)
      {
        Xyce::dout() << "We are at the end of file, size of line is " << line.size()
                     << std::endl;
      }
    }

    // better not try to do anything if getLine returned an empty line!
    if ( line.empty() )
      continue;

    // Determine what to do with the parsed line.
    ExtendedString ES1(line[0].string_);
    ES1.toUpper();

    std::string libInsideHere = "";
    if (libInside.size())
    {
      libInsideHere = libInside.back();
    }

    if ( !(libSelect != libInsideHere && libInsideHere != "" && ES1 != ".ENDL"))
    {
      char lineType;
      lineType = ES1[0];
      bool removecomponent=false;

      if (lineType >= 'A' && lineType <= 'Z')
      {

        if ( lineType == 'C' || lineType == 'D' || lineType == 'I' ||
             lineType == 'L' || lineType == 'R' || lineType == 'V')
        {
          if (remove_any_redundant_ && line.size() > 2) //make sure that there are two nodes to check!
          {
            ExtendedString node1 ( line[1].string_ );
            ExtendedString node2 ( line[2].string_ );
            node1.toUpper();
            node2.toUpper();

            removecomponent = IO::removeTwoTerminalDevice(preprocessFilter_,
                                                          lineType, node1, node2);
          }
        }
        else if (lineType == 'M' || lineType == 'Q')
        {
          if (remove_any_redundant_ && line.size() > 3) //make sure that there are three nodes to check!
          {
            ExtendedString node1 ( line[1].string_ );
            ExtendedString node2 ( line[2].string_ );
            ExtendedString node3 ( line[3].string_ );
            node1.toUpper();
            node2.toUpper();
            node3.toUpper();

            removecomponent = IO::removeThreeTerminalDevice(preprocessFilter_,
                                                            lineType, node1,
                                                            node2, node3);
          }
        }
        if (!removecomponent)
        {
          if (lineType != 'X' && lineType != 'K' && lineType != 'L')
          {
            // check for .INITCOND values
            std::string tmpName(circuitContext_->getPrefix());
            if (tmpName == "")
              tmpName = ES1;
            else
              tmpName = tmpName + ":" + ES1;

            if( mainCircuitPtr_->initCondIndex.find( tmpName ) !=
                mainCircuitPtr_->initCondIndex.end() )
            {
              // append the IC=val1...valN list to device line
              line.insert( line.end(),
                           ( mainCircuitPtr_->initCondIndex[tmpName] ).begin(),
                           ( mainCircuitPtr_->initCondIndex[tmpName] ).end() );
            }

            return true;
          }
          else if (lineType == 'L')
          {
            std::set< std::string > & cTable =
              circuitContext_->getAllCoupledInductors();

            if( cTable.find( ES1 ) == cTable.end() )
            {
              // treat as a normal device line if inductor is not coupled
              return true;
            }
            continue;
          }
          else if (lineType == 'X')
          {
            handleDeviceLine(line, libSelect, libInside);
          }
        }
        else
        {
          // This device is being removed, clear the line.
          line.clear();

          if (DEBUG_IO)
          {
            Xyce::dout() << "Netlist Parse 2:  ";
            Xyce::dout() << "removing component " << ES1 << ".  All nodes on the device";
            Xyce::dout() << " are the same."  << std::endl;
          }

          return true;
        }
      }
      else if (lineType == '.')
      {
        if (ES1 == ".SUBCKT")
        {
          // Jump to the end of the subcircuit.
          // Find the subcircuit corresponding to this instance.
          CircuitBlock* subcircuitPtr =
            currentCircuitPtr_->findSubcircuit(ExtendedString(line[1].string_).toUpper());

          ssfPtr_->setLocation(subcircuitPtr->getEndPosition());
          ssfPtr_->setLineNumber(subcircuitPtr->getLineEndPosition());

        }
        else if ( ES1 == ".INCLUDE" || ES1 == ".INCL" || ES1 == ".INC" || ES1 == ".LIB")
        {
          // HSPICE documents .INC, .INCL and .INCLUDE as being a valid .INC line
          std::string includeFile;
          std::string libSelect_new = libSelect, libInside_new;
          Xyce::IO::handleIncludeLine( netlistFilename_, line,
                                       ES1, includeFile, libSelect_new, libInside_new);
          if (includeFile != "")
          {
            int result = this->parseIncludeFile(includeFile, libSelect_new);
            if (!result)
              return result;
          }
          else if (libInside_new != "")
          {
            libInside.push_back( libInside_new );
          }
        }
        else if (ES1 == ".ENDL")
        {
          Xyce::IO::handleEndlLine ( netlistFilename_, line, libInside.back() );
          libInside.pop_back();
        }
        else if (ES1 == ".ENDS")
        {
          return parentCircuitPtr_ == 0;
        }
        else if (ES1 == ".END")
        {
          return false;
        }
      }
      else
      {
        // Do nothing for other types of lines
      }
    }
  }
  // Only get here if end of file.

  return false;
}

//----------------------------------------------------------------------------
// Function       : DistToolBase::handleDeviceLine
// Purpose        : Processor zero / serial processing of devices
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/21/2003
//----------------------------------------------------------------------------
bool DistToolBase::handleDeviceLine( TokenVector const& deviceLine,
    const std::string &libSelect, std::vector<std::string> &libInside)
{
  bool result;

  device_.clear();

  // Fully parse the device line and extract the relevant data.
  {
    //Stats::StatTop _extractData("Extract Device Data");
    //Stats::TimeBlock _extractDataTimer(_extractData);

    bool resolveParams=true;
    bool modelBinning=circuitBlock_.getModelBinningFlag();
    result = device_.extractData(netlistFilename_, deviceLine, resolveParams,modelBinning);
    if (!result)
      return result;
  }

  unordered_map<std::string, std::string> * nodeMapPtr = circuitContext_->getNodeMapPtr();
  result = instantiateDevice(device_, circuitContext_->getPrefix(), 
                             *nodeMapPtr, libSelect, libInside);

  return result;
}

//----------------------------------------------------------------------------
// Function       : DistToolBase::instantiateDevice
// Purpose        : Extract data from tokenized device and place all data about
//                  the device in Topology.
// Special Notes  : This method is misnamed as it does not actually instantiate
//                  the device.  That is done after partitioning (if parallel)
//                  and is called from Topology.
//                  This actually fills the DeviceBlock which is stored in
//                  Topology.
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 07/28/2003
//----------------------------------------------------------------------------
bool DistToolBase::instantiateDevice(
  DeviceBlock &                                 device,
  const std::string &                           prefix,
  const unordered_map<std::string,std::string> &nodeMap,
  const std::string &                           libSelect,
  std::vector<std::string> &                    libInside)
{
  bool result = true;

  unordered_map<std::string,std::string>::const_iterator nodeMapIter;
  unordered_map<std::string,std::string>::const_iterator nodeMapEnd = nodeMap.end();

  // Global params were discovered by proc 0 in pass 1 and stored in root level of
  // circuitContext_->  Then circuitContext was transmitted to all procs.

  // Now, before device parsing, move global params from circuitContext_ to deviceMgr

  // This is done here because the global parameters are needed for the parsing of
  // device lines, but are not known until after pass 1.  This point is immediately
  // before device line parsing on all processors.

  // Map device nodes.
  if (!nodeMap.empty())
  {
    //Stats::StatTop _instantiateStat("Determine Device Nodes");
    //Stats::TimeBlock _instantiateTimer(_instantiateStat);

    std::vector<std::string>::iterator nodeIter = device.getNodeValues().begin();
    std::vector<std::string>::iterator nodeIterEnd = device.getNodeValues().end();
    for (; nodeIter != nodeIterEnd; ++nodeIter)
    {
      nodeMapIter = nodeMap.find( *nodeIter );
      if (nodeMapIter != nodeMap.end())
      {
        // The device node is an external node, map it using nodeMap
        *nodeIter = nodeMapIter->second;
      }
      else if ( *nodeIter != "0" && nodeIter->substr(0,2) != "$G" &&
               !circuitContext_->globalNode( *nodeIter ) )
      {
        // The node is internal, prepend subcircuitPrefix. Note: the
        // ground node and global nodes (0 and $G*) are unchanged.
        *nodeIter = std::string( prefix + ":" + *nodeIter );
      }
    }
  }

  // If the device is not a subcircuit instance add its data to the
  // tables, otherwise, expand the instance.
  if (!device.isSubcircuitInstance())
  {
    //Stats::StatTop _instantiateStat("Process Device Data");
    //Stats::TimeBlock _instantiateTimer(_instantiateStat);

    // Once last check to verify that this is device should be instantiated, 
    // now that the node values have been resolved.
    // NOTE:  It is easier to NOT instantiate the device, than remove it from topology later!
    if (remove_any_redundant_)
    {
      if ( IO::removeDevice(preprocessFilter_, device) )
      {
        return result;
      }
    }

    // Map device name.
    if (prefix != "")
      device.setName(prefix + ":" + device.getInstanceName().getEncodedName());

    // If the device has a model, find it and instantiate it (if it has not
    // already been instantiated). Prepend the model name in the device and
    // the model with the model prefix.  Add the model to the circuit.
    if (device.getModelName() != "")
    {
      //Stats::StatTop _modelStat("Map Model Name");
      //Stats::TimeBlock _modelTimer(_modelStat);

      std::string modelPrefix;
      ParameterBlock* modelPtr;

      bool modelFound = circuitContext_->findModel(device.getModelName(), modelPtr, modelPrefix);
      if (modelFound)
      {
        // Add the model prefix to the device's model name.
        if (modelPrefix != "")
        {
          device.setModelName(modelPrefix + ":" + device.getModelName());
        }
      }
      else
      {
        Report::UserError0() << "Unable to find model " << device.getModelName()
                             << " for device " << device.getInstanceName();
        return false;
      }

      // Add the model to the circuit.
      {
        //Stats::StatTop _modelStat2("Add Model Circuit Block");
        //Stats::TimeBlock _modelTimer2(_modelStat2);
        mainCircuitPtr_->addModel(modelPtr, modelPrefix);
      }
    }

    // If the device is a Bsource, the nodes and sources appearing in
    // the expression must be appropriately prefixed.

    // THIS IS AN AWFUL HACK!  It means that *ANY* device that has
    // "SOLN_DEP" parameters must be hacked in here or they can't be used
    // in a subcircuit!

    if ((device.getNetlistDeviceType() == "B" || device.getNetlistDeviceType() == "S" || device.getNetlistDeviceType() == "C")&& (!nodeMap.empty()))
    {
      Device::Param* givenParameter = 0;
      if (device.getNetlistDeviceType() == "B")
      {
        Device::Param* I_parameter = device.findInstanceParameter( "I" );
        Device::Param* V_parameter = device.findInstanceParameter( "V" );
        if ( I_parameter->given() )
          givenParameter = I_parameter;
        else if ( V_parameter->given() )
          givenParameter = V_parameter;
      }
      else
      {
        Device::Param* C_parameter = 0;
        Device::Param* Q_parameter = 0;
        if (device.getNetlistDeviceType() == "S")
        {
            C_parameter = device.findInstanceParameter( "CONTROL" );
            if ( C_parameter->given() )
              givenParameter = C_parameter;
        }
        else if (device.getNetlistDeviceType() == "C")
        {
          C_parameter = device.findInstanceParameter( "C" );
          Q_parameter = device.findInstanceParameter( "Q" );
          // Note!  The capacitor could be a semiconductor capacitor, and
          // therefore might not even HAVE a C or Q parameter.
          // Must guard against  that.
          if ( C_parameter->given() )
            givenParameter = C_parameter;
          else if (Q_parameter->given() )
            givenParameter = Q_parameter;
          else
            givenParameter = 0;
        }
      }

      if (givenParameter && givenParameter->getType() == Xyce::Util::EXPR)
      {
        Util::Expression &expression = givenParameter->getValue<Util::Expression>();

        std::vector<std::string> names;
        std::vector<std::string> nodes;
        std::vector<std::string> instances;
        std::vector<std::string> leads;

        expression.getVoltageNodes(nodes);
        expression.getDeviceCurrents(instances);
        expression.getLeadCurrents(leads);

//ERK
#if 0
        if ( (expression.get_num( XEXP_NODE ) > 0) ||
            (expression.get_num( XEXP_INSTANCE ) > 0) ||
            (expression.get_num( XEXP_LEAD ) > 0) )
#else
        if ( (!(nodes.empty())) || (!(instances.empty())) || (!(leads.empty())) )
#endif
        {
          // If the expression has nodes or voltage source instances, get
          // the nodes and map them appropriately or add the subcircuit
          // prefix to them.
          if ( (!(nodes.empty()))     ) { names.insert( names.end(), nodes.begin(), nodes.end() ); }
          if ( (!(instances.empty())) ) { names.insert( names.end(), instances.begin(), instances.end() ); }
          if ( (!(leads.empty()))     ) { names.insert( names.end(), leads.begin(), leads.end() ); }

          std::vector<std::string> actualName;
          std::string newName, tmpName;
          unsigned int i;

          nodeMapEnd = nodeMap.end();

          for (i = 0; i < names.size(); ++i )
          {
            nodeMapIter = nodeMap.find(names[i]);
            if (nodeMapIter != nodeMapEnd)
            {
              newName = nodeMapIter->second;
            }
            else if (!(i < nodes.size() && ( names[i].substr(0,2) == "$G" ||
                     circuitContext_->globalNode(names[i]))))
            {
              newName = prefix + ":" + names[i];
            }
            else
            {
              // which leaves the one case, where we're a global node $G...,
              // in which case we just keep the original name.
              newName = names[i];
            }

            // Replace the old node name with the new node name
            // in the expression.
            if (names[i] != newName)
            {
              if (DEBUG_IO)
                Xyce::dout() << "Replacing " << names[i] << " with " << newName << std::endl;

              actualName.push_back(newName);
              tmpName = ";" + newName;
              expression.replace_name( names[i], tmpName );
              if (DEBUG_EXPRESSION) {
                Xyce::dout() << "DistToolBase::instantiateDevice:  After this replacement, get_expression returns "
                             << expression.get_expression() << std::endl
                             << " Parse Tree: " << std::endl;
                expression.dumpParseTree();
              }
            }
          }
          for (i=0 ; i<actualName.size() ; ++i)
          {
            if (DEBUG_IO)
              Xyce::dout() << "Replacing ;"<<actualName[i]<<" with " << actualName[i] << std::endl;

            newName = actualName[i];
            tmpName = ";" + newName;
            expression.replace_name( tmpName, newName);

            if (DEBUG_EXPRESSION) {
              Xyce::dout() << "DistToolBase::instantiateDevice:  After this replacement get_expression returns "
                           << expression.get_expression() << std::endl;
              expression.dumpParseTree();
            }
          }

          // Reset Bsource's expression.
          if (DEBUG_IO)
            Xyce::dout() << "DistToolBase::instantiateDevice:  After all expression handling, get_expression returns "
                         << expression.get_expression() << std::endl;
        }
      }
    }

    if (device.getNetlistDeviceType() == "L" && circuitContext_->haveMutualInductances() )
    {
      handleMutualInductance(device);
    }

    if (device.getNetlistDeviceType() == "IBIS")
    {
      handleIBISdevice(device);
    } 

    // Add device info to the main circuit device and node tables.
    {
      //Stats::StatTop _addTableStat("Add Table Data");
      //Stats::TimeBlock _addTableTimer(_addTableStat);

      circuitBlock_.addTableData(device);
    }
  }
  else
  {
    result = this->expandSubcircuitInstance(device, libSelect, libInside);
  }

  if (DEBUG_IO)
  {
    device.print();
  }

  return result;
}


//----------------------------------------------------------------------------
// Function       : DistToolBase::handleMutualInductance
// Purpose        : Post-process the mutual inductors in the current circuit,
//                  each inductor must get the list of inductors it is coupled
//                  to, the coupling coeefficient, and model information.
// Special Notes  :
// Scope          :
// Creator        : Rob Hoekstra
// Creation Date  : 08/28/04
//----------------------------------------------------------------------------
bool DistToolBase::handleMutualInductance( DeviceBlock & device )
{
  std::string subcircuitPrefix = circuitContext_->getPrefix();

  //Check for mutual inductance assoc. with this inductor
  std::vector<CircuitContext::MutualInductance> & MIs = circuitContext_->getMutualInductances();

  std::string name = device.getInstanceName().getEncodedName();
  std::string::size_type pos = name.find_last_of(":");
  if (pos != std::string::npos)
    name = name.substr( pos + 1, name.length() - (pos + 1));

  for( unsigned int i = 0; i < MIs.size(); ++i)
  {
    if( MIs[i].inductors.count( name ) )
    {
      // add mutual inductance info to this inductor
      if( MIs[i].inductors.begin()->first == name )
      {
        Device::Param param;

        param.setTag( "FIRSTINDUCTOR" );
        param.setGiven(true);
        param.setVal( 1 );
        device.addInstanceParameter( param );
      }

      if( MIs[i].model != "" )
      {
        std::string modelName(MIs[i].model);

        Device::Param param;

        param.setTag( "NONLINEARCOUPLING" );
        param.setVal( 1 );
        param.setGiven(true);
        device.addInstanceParameter( param );

        std::string modelPrefix;
        ParameterBlock* MI_modelPtr;
        bool modelFound = circuitContext_->
              findModel( modelName, MI_modelPtr, modelPrefix );

        if(modelFound)
        {
          // Check whether the current context is dependent on subcircuit
          // parameters (via a "params:" list on the .subckt line. This
          // will affect the model prefixing.
          if (circuitContext_->hasSubcircuitParams() &&
                MI_modelPtr->hasExpressionValuedParams())
          {
            // If there is a subcircuit parameter dependency, we will
            // assume the worst and add the subcircuit instance prefix
            // to the model prefix. This will insure that model parameters
            // that might depend on subcircuit parameters are assigned to
            // a unique model for each subcircuit instance.
            if( subcircuitPrefix != "" )
              modelPrefix = subcircuitPrefix + ":" + modelPrefix;
          }

          // Add the model prefix to the device's model name.
          if( modelPrefix != "" ) modelName = modelPrefix + ":" + modelName;

          device.setModelName( modelName );
        }
        else
        {
          Report::UserError0() << "Unable to find mutual inductor model " << modelName;
          return false;
        }

        // Add the model to the circuit.
        mainCircuitPtr_->addModel(MI_modelPtr, modelPrefix);
      }
      {
        Device::Param param;

        param.setTag( "COUPLING" );
        param.setVal( MIs[i].coupling );
        param.setGiven(true);
        device.addInstanceParameter( param );
      }

      // add the set of inductors to which it is coupled
      std::map<std::string,std::string>::iterator iterI = MIs[i].inductors.begin();
      std::map<std::string,std::string>::iterator  endI = MIs[i].inductors.end();
      for( ; iterI != endI; ++iterI )
      {
        std::string ci = iterI->first;

        if( name != ci )
        {
          Device::Param param;

          param.setTag( "COUPLEDINDUCTOR" );
          if( subcircuitPrefix != "" )
          {
            ci = subcircuitPrefix + ":" + ci;
          }

          param.setVal( ci );
          param.setGiven (true);
          device.addInstanceParameter( param );

          param.setTag("COUPLEDINDUCTANCE");

          if (DEBUG_IO)
            Xyce::dout() << "coupledinductance value from " << name << " to " << ci
                         << " is " << iterI->second << std::endl;

          param.setVal( iterI->second);
          param.setGiven (true);
          device.addInstanceParameter( param );
        }
      }
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : DistToolBase::handleIBISdevice
// Purpose        : The IBIS devices need the fully qualified names of their
//                  interface nodes in order to auto-generate their various
//                  I-V tables.  This function puts that list into the NODELIST
//                  instance parameter.
// Special Notes  : This approach is a bit of a "hack", but the same hack is
//                  used for the mutual inductor device.
// Scope          :
// Creator        : Pete Sholander, Electrical Models & Simulation
// Creation Date  : 08/28/18
//----------------------------------------------------------------------------
bool DistToolBase::handleIBISdevice( DeviceBlock & device )
{
  if (!device.getNodeValues().empty())
  {
    Device::Param param;
    std::string subcircuitPrefix = circuitContext_->getPrefix();
    std::vector<std::string> fullNodeNames;
    
    std::vector<std::string>::iterator nodeIter = device.getNodeValues().begin();
    std::vector<std::string>::iterator nodeIterEnd = device.getNodeValues().end();

    for (; nodeIter != nodeIterEnd; ++nodeIter)
    {
      if ( !subcircuitPrefix.empty() )
      {
        fullNodeNames.push_back(subcircuitPrefix + ":"+ *nodeIter);
      }
      else
      {
        fullNodeNames.push_back(*nodeIter);
      }
    }

    param.setTag("NODELIST");
    param.setVal(fullNodeNames);
    param.setGiven (true);
    device.addInstanceParameter( param );
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : DistToolBase::find_IC_NODESET_OptionBlock
// Purpose        : Process the option block to find any IC and NODESET
//                : statements that are attached to the current subcircuit
// Special Notes  : Check if any of the nodes are interface nodes.
// Scope          :
// Creator        : Heidi Thornquist
// Creation Date  : 08/04/2016
//----------------------------------------------------------------------------
void DistToolBase::find_IC_NODESET_OptionBlock(const std::string& modelName,
                                               const std::string& subcircuitPrefix,
                                               const std::vector<std::string>& subcircuitNodes,
                                               const std::vector<std::string>& subcircuitInstanceNodes)

{
  if (check_IC_NODESET_OptionBlock())
  {
    // Add any .IC or .NODESET statements for this subcircuit to the top level
    // option block and resolve any parameters in the current context.
    std::list<Util::OptionBlock>::const_iterator it = options_.begin();
    for ( ; it != options_.end(); ++it )
    {
      if (it->getName()=="IC" || it->getName()=="NODESET")
      {
        if (it->begin()->tag() == "SUBCKT")
        {
          if (it->begin()->stringValue() == modelName)
          {
            // Copy this option block and fill out the template (resolve nodes, params)
            Util::OptionBlock newOB( *it );
            newOB.removeParam( "SUBCKT" );
  
            Util::ParamList::iterator iterPar = newOB.begin();
            Util::ParamList::iterator endPar = newOB.end();
            for (; iterPar != endPar; ++iterPar)
            {
              // Fix up the node name.
              if (iterPar->tag() == "V")
              {
                // The next param is the node name.
                iterPar++;
                std::string newTag;
                std::vector<std::string>::const_iterator subIt = std::find( subcircuitNodes.begin(), subcircuitNodes.end(), iterPar->tag() );
                if ( subIt != subcircuitNodes.end() )
                {
                  int dist = std::distance( subcircuitNodes.begin(), subIt );
                  newTag = subcircuitInstanceNodes[ dist ];
                }
                else
                {
                  newTag = subcircuitPrefix + ":" + iterPar->tag();
                }
                iterPar->setTag( newTag );
              }
  
              // Resolve any local params.
              circuitContext_->resolveParameter((*iterPar));
            }
  
            // Add this new option block to the current list.
            addOptions_.push_back( newOB );
          }
        }  
      }
    }
  }
}

} // namespace IO
} // namespace Xyce

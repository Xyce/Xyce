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
// Purpose        : Defines the DistToolDefault class.  Distribution tool
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
#include <N_IO_DistToolDefault.h>
#include <N_IO_ParsingHelpers.h>
#include <N_PDS_Comm.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Stats.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : DistToolDefault::DistToolDefault
// Purpose       : ctor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistToolDefault::DistToolDefault(
  Parallel::Communicator *                 pdsCommPtr,
  CircuitBlock &                           circuit_block,
  std::map<std::string,FileSSFPair>      & ssfMap,
  std::map<std::string, IncludeFileInfo> & iflMap,
  const ParsingMgr                       & parsing_manager
  )
  : DistToolBase(pdsCommPtr, circuit_block, ssfMap,parsing_manager),
    currProc_(0),
    iflMap_(iflMap),
    procDeviceCount_(0),
    devices_(0),
    deviceLinesSent_(0),
    bufs_(),
    bufSize_(),
    subcircuitNames_(),
    subcircuitNodes_(),
    subcircuitPrefixes_(),
    subcircuitParams_()
{
  // Set the circuit context and options table.
  setCircuitContext();
  setCircuitOptions();

  // determine how many devices to send to each far proc
  devices_ = (circuit_block.getCircuitContextPtr())->getTotalDeviceCount();
  procDeviceCount_ = devices_ / numProcs_;

  // procID == 0 may not be true if we're running in a hierarchical parallel context
  // so check if numProc() == 1 too.
  if (pdsCommPtr_->procID() == 0)
  {
    currProc_ = numProcs_ == 1 ? 0 : 1;
  }
}

//-----------------------------------------------------------------------------
// Function      : DistToolDefault::circuitDeviceLine
// Purpose       : Send a circuit device line to current proc
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolDefault::circuitDeviceLine(TokenVector & deviceLine )
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if( currProc_ != 0 )
    {
      int size = 0;

      // count line type
      size += sizeof(char);

      // count line size
      size += sizeof(int);

      // count device line
      int deviceLineSize = deviceLine.size();
      for (int i = 0; i < deviceLineSize; ++i)
      {
        size += Xyce::packedByteCount(deviceLine[i]);
      }

      // flush buffer as needed
      send(size);

      // pack line type; "d"evice line
      char lineType = 'd';
      pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );

      // pack device line size
      pdsCommPtr_->pack( &deviceLineSize, 1, charBuffer_, charBufferSize_, charBufferPos_ );

      // pack device line
      for (int i = 0; i < deviceLineSize; ++i)
      {
        Xyce::pack(deviceLine[i],  charBuffer_, charBufferSize_, charBufferPos_, pdsCommPtr_ );
      }

      // increment device line counter
      deviceLinesSent_++;

      // manage proc change if n/p reached
      if (deviceLinesSent_ >= procDeviceCount_)
      {
        int minus1 = -1;

        // flush buffer
        send();

        // allow this proc to exit receive loop and begin processing
        pdsCommPtr_->send( &minus1, 1, currProc_ );

        // reset device line counter
        deviceLinesSent_ = 0;

        // move to next processor
        currProc_++;

        // switch to proc 0 if nearly complete
        if (currProc_ == numProcs_)
        {
          currProc_ = 0;
        }

        // otherwise, prepare currProc for next set of device lines
        else
        {
          int length = netlistFilename_.size();
          lineType = 'f';

          // flush buffer as needed
          send(sizeof(char) + sizeof(int) + length);

          // pack the filename
          pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
          pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
          pdsCommPtr_->pack( netlistFilename_.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

          // pack the subcircuit names
          int subcircuitNamesSize = subcircuitNames_.size();
          for (int i = 0; i < subcircuitNamesSize; ++i)
          {
            if( currProc_ != 0 )
            {
              packSubcircuitData(subcircuitNames_[i], subcircuitNodes_[i],
                                 subcircuitPrefixes_[i], subcircuitParams_[i]);
            }
          }
        }
      }

      return true;
    }
    else
      return false;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : DistToolDefault::endDeviceLines()
// Purpose       : Make sure other processors have been told that device processing
//                 is ended.  This can be a problem with small device count circuits
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/06/06
//-----------------------------------------------------------------------------
void DistToolDefault::endDeviceLines()
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if (currProc_ > 0)
    {
      int minus1 = -1;

      // flush buffer
      send();

      // end parsing on current node
      pdsCommPtr_->send( &minus1, 1, currProc_ );
      ++currProc_;

      // end parsing on remaining nodes
      for ( ; currProc_ < numProcs_ ; ++currProc_)
      {
        pdsCommPtr_->send( &minus1, 1, currProc_ );
      }

      // complete any remaining parser action on node 0
      currProc_ = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DistToolDefault::circuitStart
// Purpose       : change current subcircuit context after a new subcircuit is
//               : started
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolDefault::circuitStart( std::string const & subcircuitName,
                                     std::vector<std::string> const & nodes,
                                     std::string const & prefix,
                                     std::vector<Device::Param> const & params )
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if( currProc_ != 0 )
    {
      // save new context data on stacks
      subcircuitNames_.push_back( subcircuitName );
      subcircuitPrefixes_.push_back( prefix );
      subcircuitNodes_.push_back( nodes );
      subcircuitParams_.push_back( params );

      // pack data into buffer
      packSubcircuitData( subcircuitName, nodes, prefix, params );
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistToolDefault::circuitEnd
// Purpose       : change current subcircuit context to previous context
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolDefault::circuitEnd()
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if( currProc_ != 0 )
    {
      // remove latest context data on stacks
      subcircuitNames_.pop_back();
      subcircuitPrefixes_.pop_back();
      subcircuitNodes_.pop_back();
      subcircuitParams_.pop_back();

      // flush buffer as needed
      send(sizeof(char));

      // pack line type; subcircuit "e"nd
      char lineType = 'e';
      pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::setFileName
// Purpose       : Change name of netlist file
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DistToolDefault::setFileName(std::string const & fileNameIn)
{
  netlistFilename_ = fileNameIn;
  circuitBlock_.setFileName( fileNameIn );

  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    char lineType = 'f';
    int length = netlistFilename_.size();

    // flush buffer as needed
    send(sizeof(char) + sizeof(int) + length);

    pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    pdsCommPtr_->pack( netlistFilename_.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );
  }
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::distributeDevices()
// Purpose       : Distribute the devices according to the strategy.
// Special Notes : The strategy is first-come-first-served, flattening the circuit
//               : as you go through the netlist files.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 08/28/2014
//-----------------------------------------------------------------------------
void DistToolDefault::distributeDevices()
{
  Xyce::Parallel::Machine comm = pdsCommPtr_->comm();

  if (DEBUG_IO) {
    Xyce::dout() << "Pass 2 parsing of netlist file and distributing devices: " <<
      circuitBlock_.getNetlistFilename() << std::endl;
  }

  //Stats::StatTop _distDataStat("Distribute Data");
  {
  //Stats::TimeBlock _distDataTimer(_distDataStat);
  
  if (Parallel::rank(comm) == 0)
  {
    setFileName(circuitBlock_.getNetlistFilename());

    ssfPtr_ = ssfMap_[netlistFilename_].second;
    ssfPtr_->setLocation( circuitBlock_.getStartPosition() );
    ssfPtr_->setLineNumber( circuitBlock_.getLineStartPosition() );

    // If this is the main circuit, skip over the title line and continue.
    std::ifstream* netlistIn = ssfMap_[netlistFilename_].first;
    std::string title("");
    netlistIn->clear();
    netlistIn->seekg(0, std::ios::beg);
    Xyce::IO::readLine( *netlistIn, title ); 
    ssfPtr_->changeCursorLineNumber( 1 );

    // Instantiate all devices in the current context.
    std::string libSelect;
    std::vector<std::string> libInside;
    TokenVector line;

    while (getLine(line, libSelect, libInside))
    {
      if (!line.empty() && compare_nocase(line[0].string_.c_str(), ".ends") != 0)
      {
        // parse locally if distool does not distribute
        if( !circuitDeviceLine( line ) )
        {
          handleDeviceLine( line, libSelect, libInside );
        }
      }
    }

    // send MIs if present in top circuit
    if( circuitContext_->haveMutualInductances() )
    {
      // Distribute each YMI?|name device line, end of circuit
      int n = circuitContext_->getNumMILines();
      for( int i = 0; i < n; ++i )
      {
        // parse locally if not distributed to another processor
        if( !circuitDeviceLine( circuitContext_->getMILine( i ) ) )
        {
          handleDeviceLine( circuitContext_->getMILine( i ), libSelect, libInside );
        }
      }
    }

    endDeviceLines();

  }
  else
  {
    // Wait to receive the circuit data from processor 0.
    receiveCircuitData();
  }

  int bsize = 0;

  // Send any additional option blocks that have been found during device distribution.
  // NOTE:  These are .IC or .NODESET lines within subcircuits, etc.
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if (pdsCommPtr_->procID() == 0)
    {
      // resize buffer if necessary
      int tmpSize = sizeof(int);
      for (std::list<Util::OptionBlock>::const_iterator it = addOptions_.begin(), end = addOptions_.end(); it != end; ++it)
        tmpSize += Xyce::packedByteCount(*it);

      if ( tmpSize > charBufferSize_ )
      {
        charBufferSize_ = tmpSize;
        delete[] charBuffer_;
        charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];
      }

      bsize = Xyce::IO::packCircuitOptions(addOptions_, charBuffer_, charBufferSize_, pdsCommPtr_);
    }
   
    // broadcast options
    pdsCommPtr_->bcast( &bsize, 1, 0 );
    if ( bsize > charBufferSize_ )
    {
      charBufferSize_ = bsize;
      delete[] charBuffer_;
      charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];
    }
    pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

    // receive data
    if (pdsCommPtr_->procID() != 0)
    {
      Xyce::IO::unpackCircuitOptions(addOptions_, charBuffer_, bsize, pdsCommPtr_);
    }
  }

  // Send node alias map.
  if (Parallel::rank(comm) == 0)
  {  
    bsize = Xyce::IO::packAliasNodeMap(circuitBlock_.getAliasNodeMap(), charBuffer_, charBufferSize_, pdsCommPtr_);
  }

  // broadcast alias node map
  pdsCommPtr_->bcast( &bsize, 1, 0 );
  pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

  if (Parallel::rank(comm) != 0)
  {
    Xyce::IO::unpackAliasNodeMap(circuitBlock_.getAliasNodeMap(), charBuffer_, charBufferSize_, pdsCommPtr_);
  }

  if (DEBUG_IO)
    Xyce::dout() << "Done with pass 2 netlist file parsing and device distribution" << std::endl;

  // Just in case an error was reported in parsing the netlist.
  N_ERH_ErrorMgr::safeBarrier(comm);

  }
  //std::cout << "Processor " << pdsCommPtr_->procID() << " took " << _distDataStat.getTop().getMetric<Xyce::Stats::WallTime>().getAccumulatedLap() << " to distribute devices." << std::endl;
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::packSubcircuitData
// Purpose       : pack subcircuit data into buffer; helper for circuitStart()
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolDefault::packSubcircuitData( std::string const & subcircuitName,
                                           std::vector<std::string> const & nodes,
                                           std::string const & prefix,
                                           std::vector<Device::Param> const & params )
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    int size = 0;

    // count line type, name size, and name chars
    size += sizeof(char) + sizeof(int) + subcircuitName.length();

    // count # of nodes,
    size += sizeof(int);

    // count node sizes and node chars
    std::vector<std::string>::const_iterator it_sbL = nodes.begin();
    std::vector<std::string>::const_iterator it_seL = nodes.end();
    for( ; it_sbL != it_seL; ++it_sbL  )
    {
      size += sizeof(int) + it_sbL->length();
    }

    // count prefix size and prefix chars
    size += sizeof(int) + prefix.length();

    // count # of params
    size += sizeof(int);

    // count params
    std::vector<Device::Param>::const_iterator it_pbL = params.begin();
    std::vector<Device::Param>::const_iterator it_peL = params.end();
    for( ; it_pbL != it_peL; ++it_pbL  )
    {
      size += Pack<Device::Param>::packedByteCount(*it_pbL);
    }

    // flush buffer as needed
    send(size);


    // pack line type; subcircuit "s"tart
    char lineType = 's';
    pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack name size and name chars
    int length = subcircuitName.length();
    pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    pdsCommPtr_->pack( subcircuitName.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack # of nodes
    size = nodes.size();
    pdsCommPtr_->pack( &size, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack nodes name sizes and chars
    it_sbL = nodes.begin();
    for( ; it_sbL != it_seL; ++it_sbL  )
    {
      length = it_sbL->length();
      pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
      pdsCommPtr_->pack( it_sbL->c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );
    }

    // pack prefix size and prefix chars
    length = prefix.length();
    pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    pdsCommPtr_->pack( prefix.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack # of params
    size = params.size();
    pdsCommPtr_->pack( &size, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack params
    it_pbL = params.begin();
    for( ; it_pbL != it_peL; ++it_pbL  )
    {
      Pack<Device::Param>::pack(*it_pbL, charBuffer_, charBufferSize_, charBufferPos_, pdsCommPtr_ );
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::receiveCircuitData
// Purpose       : receive all data from proc 0 and store it in device buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolDefault::receiveCircuitData()
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    //Stats::StatTop _getLineStat("Receive Data");
    
    {
    //Stats::TimeBlock _getLineTimer(_getLineStat);
    int bsize;
    char *currBuffer;

    // receive the device lines first
    while (true)
    {
      pdsCommPtr_->recv( &bsize, 1, 0 );

      // check for halt request
      if (bsize < 0)
      {
        break;
      }

      currBuffer = new char[bsize];
      bufs_.push_back(currBuffer);
      bufSize_.push_back(bsize);
      //std::cout << "Processor " << pdsCommPtr_->procID() << " receiving buffer number " << bufs_.size() << std::endl;

      pdsCommPtr_->recv(currBuffer, bsize, 0);
    }
    }
    //std::cout << "Processor " << pdsCommPtr_->procID() << " took " << _getLineStat.getTop().getMetric<Xyce::Stats::WallTime>().getAccumulatedLap() << " seconds receiving " << bufs_.size() << " buffers of size " << bufSize_[bufs_.size()-2] << std::endl;

    bool ok = true;

    //Stats::StatTop _processDeviceStat("Process Devices");
    {
      //Stats::TimeBlock _processDeviceTimer(_processDeviceStat);
      ok = processDeviceBuffer();
    }
    //std::cout << "Processor " << pdsCommPtr_->procID() << " took " << _processDeviceStat.getTop().getMetric<Xyce::Stats::WallTime>().getAccumulatedLap() << " to process device buffers." << std::endl;

    return ok;
  }
  else
    return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::send
// Purpose       : send buffer from proc 0 and adjust buffer limits
// Special Notes : size == -1, and no args send(), will force a transmission
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void
DistToolDefault::send(
  int                   size)
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    // transmit buffer if next set of object will fill it completely, or forced
    if ((charBufferPos_ + size >= charBufferSize_) || (size == -1))
    {
      if (DEBUG_DISTRIBUTION)
        dout() << "node " << pdsCommPtr_->procID() << " packed "
               << charBufferPos_ << " bytes into " << charBufferSize_ << " byte buffer"
               << " [data size:  " << size << "]" << std::endl
               << "node " << currProc_ << " was sent "
               << deviceLinesSent_ << " of " << procDeviceCount_ << " devices "
               << (float)deviceLinesSent_ / procDeviceCount_ * 100.0
               << ( size == -1 ? "% complete  [flushed]" : "% complete"  ) << std::endl;

      // tx buffer size
      pdsCommPtr_->send( &charBufferPos_, 1, currProc_ );

      // tx buffer contents
      pdsCommPtr_->send( charBuffer_, charBufferPos_, currProc_ );

      // reset counters
      charBufferPos_ = 0;

      // adjust buffer size if objects larger than existing buffer
      if (size > charBufferSize_)
      {

        // free memory
        // delete [] charBuffer_;

        // set new buffer size
        charBufferSize_ = size;

        // reserve memory; including offset
        charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

        if (DEBUG_DISTRIBUTION)
          dout() << "node " << pdsCommPtr_->procID()
                 << " resized buffer to " << charBufferSize_
                 << " with offset " << sizeof(char) + sizeof(int) << std::endl;
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::processDeviceBuffer
// Purpose       : process all data in the buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
DistToolDefault::processDeviceBuffer()
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    int size, bsize, length, pos, i;

    char *currBuffer;
    unsigned int ib;

    CircuitContext* circuit_context = circuitBlock_.getCircuitContextPtr();

    // process received device lines, contexts info, and exit calls
    for (ib=0 ; ib<bufs_.size() ; ++ib)
    {
      currBuffer = bufs_[ib];
      bsize = bufSize_[ib];

      // get ready to process buffer
      pos = 0;

      // break apart buffered lines and parse
      while( pos < bsize )
      {
        char lineType = '\0';
        // get the linetype marker that indicates incoming line
        pdsCommPtr_->unpack( currBuffer, bsize, pos, &lineType, 1 );

        // process line according to type
        switch( lineType )
        {
          case 'd': // process a device line
          {
            // unpack the device line
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &size, 1 );
            TokenVector deviceLine( size );

            for( i = 0; i < size; ++i )
            {
              Xyce::unpack(deviceLine[i], currBuffer, bsize, pos, pdsCommPtr_ );
            }

            // hand to circuit block for processing
            std::string libSelect("");
            std::vector<std::string> libInside; 
            handleDeviceLine( deviceLine, libSelect, libInside );

            break;
          }

          case 's': // subcircuit start found
          {
            // get the subcircuit name
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
            std::string subcircuitName(std::string( ( currBuffer + pos ), length ));
            pos += length;

            // get the nodes for this subcircuit call
            std::vector<std::string> nodes;
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &size, 1 );
            for( i = 0; i < size; ++i )
            {
              pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
              nodes.push_back( std::string( ( currBuffer + pos ), length ) );
              pos += length;
            }

            // get the prefix
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
            std::string prefix(std::string( ( currBuffer + pos ), length ));
            pos += length;

            // get params
            std::vector<Device::Param> params;
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &size, 1 );
            for( i = 0; i < size; ++i )
            {
              Device::Param param;
              Pack<Device::Param>::unpack(param, currBuffer, bsize, pos, pdsCommPtr_ );
              params.push_back( param );
            }

            // send to cktblk
            circuit_context->setContext( subcircuitName, prefix, nodes );
            circuit_context->resolve( params );

            break;
          }

          case 'e': // subcircuit end found
          {
            // adjust the circuit context pointer
            circuit_context->restorePreviousContext();

            break;
          }

          case 'f': // change netlist file name
          {
            pdsCommPtr_->unpack( currBuffer, bsize, pos, &length, 1 );
            netlistFilename_ = std::string( ( currBuffer + pos ), length );
            pos += length;
            setFileName(netlistFilename_);

            break;
          }

          default:  // something went wrong
          {
            Report::DevelFatal().in("DistToolDefault::processDeviceBuffer")
              << "Node " << pdsCommPtr_->procID() << " received invalid message type \"" << lineType;
          }
        }
      }
      delete [] bufs_[ib];
    }

    bufs_.clear();
    bufSize_.clear();
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolDefault::broadcastGlobalData
// Purpose       : send bufsize, options, and context to procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
DistToolDefault::broadcastGlobalData()
{
  return DistToolBase::broadcastGlobalData();
}

//--------------------------------------------------------------------------
// Function      : DistToolDefault::parseIncludeFile
// Purpose       : Jump to include file for 2nd pass
// Special Notes :
// Creator       : Rob Hoekstra
// Creation Date : 08/02/2004
//--------------------------------------------------------------------------
bool DistToolDefault::parseIncludeFile(std::string const& includeFile,
             const std::string &libSelect)
{
  // Save current ssfPtr_ and netlistFilename_.
  SpiceSeparatedFieldTool* oldssfPtr = ssfPtr_;
  std::string old_netlistFilename(netlistFilename_);
  int oldLineNumber = ssfPtr_->getLineNumber();
  std::streampos oldFilePos = ssfPtr_->getFilePosition();

  setFileName(includeFile);

  if (DEBUG_IO)
    Xyce::dout() << "Parsing include file Pass 2: " << includeFile << std::endl;

  // Find SSF for file
  if( !ssfMap_.count( includeFile ) )
  {
    Report::UserError() << "Could not find include file SSF " << includeFile;
    restorePrevssfInfo(oldssfPtr, old_netlistFilename, oldFilePos, oldLineNumber);
    return false;
  }
  ssfPtr_ = ssfMap_[includeFile].second;

  //Save current position and move to beginning
  ssfPtr_->setLocation(0);
  ssfPtr_->setLineNumber( 1 );

  TokenVector line;
  std::vector<std::string> libInside;
  while (getLine(line, libSelect, libInside)) 
  {
    if (!line.empty() && compare_nocase(line[0].string_.c_str(), ".ends") != 0)
    {
      // parse locally if distool does not distribute
      if( !circuitDeviceLine( line ) )
      {
        handleDeviceLine( line, libSelect, libInside);
      }
    }
  }

  // Restore old ssfPtr_ and netlistFilename_.
  restorePrevssfInfo(oldssfPtr, old_netlistFilename, oldFilePos, oldLineNumber);

  if (DEBUG_IO)
    Xyce::dout() << "Done with include file Pass 2: " << includeFile << std::endl;

  return true; // Only get here on success.
}

//--------------------------------------------------------------------------
// Function      : DistToolDefault::restorePrevssfInfo
// Purpose       : This is a helper function for parseIncludeFile(). It
//                 restores the information about the previous file.  It
//                 should be called before each return statement in that
//                 function.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/11/2020
//--------------------------------------------------------------------------
void DistToolDefault::restorePrevssfInfo(
    SpiceSeparatedFieldTool* oldssfPtr,
    const std::string& old_netlistFilename,
    int oldFilePos,
    int oldLineNumber)
{
  // Restore old ssfPtr_ and netlistFilename_.
  ssfPtr_ = oldssfPtr;
  setFileName(old_netlistFilename);

  // get the location in the file just in case we are re-entering the file (only with .lib)
  ssfPtr_->setLocation(oldFilePos);
  ssfPtr_->setLineNumber(oldLineNumber);

  return;
}

//--------------------------------------------------------------------------
// Function      : DistToolDefault::expandSubcircuitInstance
// Purpose       : Expand a subcircuit instance by adding the devices and
//                 device models that compose the subcircuit to the main
//                 (top level) circuit.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/28/2001
//--------------------------------------------------------------------------
bool DistToolDefault::expandSubcircuitInstance(
  IO::DeviceBlock &     subcircuitInstance,
  const std::string &   libSelect,
  std::vector<std::string> & libInside)
{

  // Set subcircuitPrefix.
  std::string subcircuitPrefix;
  if ( circuitContext_->getPrefix() != "" )
  {
    subcircuitPrefix = circuitContext_->getPrefix() +
      parsingMgr_.getSeparator() + subcircuitInstance.getInstanceName().getEncodedName();
  }
  else
  {
    subcircuitPrefix = subcircuitInstance.getInstanceName().getEncodedName();
  }

  // Find the subcircuit corresponding to this instance.
  CircuitBlock* subcircuitPtr = currentCircuitPtr_->findSubcircuit(subcircuitInstance.getModelName());
  if ( subcircuitPtr == NULL )
  {
    Report::UserError0().at(netlistFilename_, subcircuitInstance.getLineNumber())
      << "Subcircuit " << subcircuitInstance.getModelName()
      << " has not been defined for instance " << subcircuitInstance.getInstanceName();
    return false;
  }

  // For hierarchical circuits, set the pointers accordingly. 
  parentCircuitPtr_ = currentCircuitPtr_;
  currentCircuitPtr_ = subcircuitPtr;
  CircuitBlock* oldParentCircuitPtr = parentCircuitPtr_;

  // get the list of X nodes
  const std::vector< std::string >& subcircuitInstanceNodes =
    subcircuitInstance.getNodeValues();

  // Set the context for this subcircuit instance.
  std::string subcircuitName(subcircuitInstance.getModelName());
  bool cresult = circuitContext_->setContext(subcircuitName, subcircuitPrefix, subcircuitInstanceNodes);
  if (!cresult)
  {
    Report::UserError0().at(netlistFilename_, subcircuitInstance.getLineNumber())
      << "Error invoking subcircuit " << subcircuitInstance.getModelName()
      << " instance " << subcircuitInstance.getInstanceName();
    return false;
  }

  // get the list of .SUBCKT nodes
  std::vector<std::string> subcircuitNodes = circuitContext_->getNodeList();

  // Make sure the subcircuit instance and subcircuit definition agree on the
  // number of nodes.
  if ( subcircuitInstanceNodes.size() != subcircuitNodes.size() )
  {
    Report::UserError0() << "Number of nodes for subcircuit instance " << subcircuitInstance.getInstanceName()
                         << " does not agree with number of nodes in subcircuit "
                         << circuitContext_->getCurrentContextName();
    return false;
  }

  // these iterators loop over the nodes in this context (i.e. the real node
  // connected to the interface nodes)
  std::vector<std::string>::const_iterator circuitInstanceNodeIt = subcircuitInstanceNodes.begin();
  std::vector<std::string>::const_iterator endCircuitInstanceNodeIt = subcircuitInstanceNodes.end();

  // these iterators loop over the subcircuits interface nodes
  std::vector<std::string>::iterator subcircuitNodeInterfaceIt = subcircuitNodes.begin();
  std::vector<std::string>::iterator endSubcircuitNodeInterfaceIt = subcircuitNodes.end();

  // add the interface nodes of this subcircuit to the aliasNodeMap so that
  // we can look up the aliases later if needed
  while( ( circuitInstanceNodeIt != endCircuitInstanceNodeIt ) &&
         ( subcircuitNodeInterfaceIt != endSubcircuitNodeInterfaceIt ) )
  {
    std::string key( subcircuitPrefix + parsingMgr_.getSeparator() + *subcircuitNodeInterfaceIt );
  
    if ( ( mainCircuitPtr_->getAliasNodeMapHelper() ).find( key ) !=
         ( mainCircuitPtr_->getAliasNodeMapHelper() ).end() )
    {
      (mainCircuitPtr_->getAliasNodeMap())[key] = *circuitInstanceNodeIt;

      if (DEBUG_IO)
        Xyce::dout() << "Found node alias:  " << key << " ==> " << *circuitInstanceNodeIt << std::endl;
    }

    circuitInstanceNodeIt++;
    subcircuitNodeInterfaceIt++;
  }

  // Save the current file pointer and netlist filename for later. 
  SpiceSeparatedFieldTool * oldssf = ssfPtr_;
  std::string oldNetlistFilename = netlistFilename_;

  // Locate the subcircuit in the netlist file. It can either be in
  // the file currently being read, or in a separate include file.
  SpiceSeparatedFieldTool * newssf = ssfPtr_;

  if (subcircuitPtr->getNetlistFilename() != netlistFilename_)
  { // The subcircuit is in an include file.
    // Get SSF from Pass 1's ssf map
    if( ssfMap_.count( subcircuitPtr->getNetlistFilename()) )
      newssf = ssfMap_[subcircuitPtr->getNetlistFilename()].second;
    else
    {
      Report::UserError().at(netlistFilename_, subcircuitInstance.getLineNumber())
        << "Can't find include file " << subcircuitPtr->getNetlistFilename();
    }

    netlistFilename_ = subcircuitPtr->getNetlistFilename();
  }

  // Store old positions in the parent netlist file.
  std::streampos oldLoc = oldssf->getFilePosition();
  int oldLine = oldssf->getLineNumber();

  // Set file pointer to location of subcircuit.
  // Set the position of the subcircuit in its file.
  ssfPtr_ = newssf;
  ssfPtr_->setLocation( subcircuitPtr->getStartPosition() );
  ssfPtr_->setLineNumber( subcircuitPtr->getLineStartPosition() );

  // Resolve parameters and functions in the current context.
  bool result;
  std::vector<Device::Param> subcircuitInstanceParams;
  subcircuitInstance.getInstanceParameters(subcircuitInstanceParams);

  result = circuitContext_->resolve(subcircuitInstanceParams);
  if (!result)
    return result;

  // Add any .IC or .NODESET statements for this subcircuit to the top level 
  // option block and resolve any parameters in the current context.
  // Pass in local nodes and instance nodes so that the IC/NODESET node is identified properly.
  find_IC_NODESET_OptionBlock( subcircuitInstance.getModelName(), subcircuitPrefix,
                               subcircuitNodes, subcircuitInstanceNodes );

  // Tell the distribution tool that the context has changed.
  circuitStart( subcircuitName, subcircuitInstanceNodes, subcircuitPrefix,
                subcircuitInstanceParams );

  // Instantiate the devices in this subcircuit instance.
  TokenVector line;

  // If this is a simple single-device subcircuit, jump to the device line directly
  if ( subcircuitPtr->getSimpleSingleDevice() )
  {
    //std::cout << "Subcircuit " <<  subcircuitInstance.getModelName() << " is a simple single device subcircuit!" << std::endl;
    ssfPtr_->setLocation( subcircuitPtr->getDevicePosition() );
    ssfPtr_->setLineNumber( subcircuitPtr->getDeviceLine() );
  }

  while (getLine(line, libSelect, libInside))
  {
    if (!line.empty() && compare_nocase(line[0].string_.c_str(), ".ends") != 0)
    {
      // parse locally if distool does not distribute
      if( !circuitDeviceLine( line ) )
      {
        handleDeviceLine( line, libSelect, libInside );
      }
    }
    // There was only one device to parse and distribute, so move on.
    if ( subcircuitPtr->getSimpleSingleDevice() )
      break;
  }

  // send MIs if present
  if( circuitContext_->haveMutualInductances() )
  {
    // Distribute each YMI?!name device line, end of circuit
    int n = circuitContext_->getNumMILines();
    for( int i = 0; i < n; ++i )
    {
      // parse locally if distool does not distribute; normally the
      // DistToolDefault::instantiateDevices() performs this step
      if( !circuitDeviceLine( circuitContext_->getMILine( i ) ) )
      {
        handleDeviceLine( circuitContext_->getMILine( i ), libSelect, libInside );
      }
    }
  }

  // Return to previous context and line numbers.
  circuitContext_->restorePreviousContext();

  // Return file pointer to parent file pointers and location.
  ssfPtr_ = oldssf;
  ssfPtr_->setLocation( oldLoc );
  ssfPtr_->setLineNumber( oldLine );
  netlistFilename_ = oldNetlistFilename;

  // For hierarchical circuits, reset the pointers accordingly. 
  currentCircuitPtr_ = parentCircuitPtr_;
  parentCircuitPtr_ = oldParentCircuitPtr;

  // Tell the distribution tool that the current context has ended and
  // to switch back to the previous context.
  circuitEnd();

  // Only get here on success.
  return true;
}


} // namespace IO
} // namespace Xyce

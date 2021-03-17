//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Defines the DistToolFlatRoundRobin class.  Distribution tool
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
#include <N_IO_DistToolFlatRoundRobin.h>
#include <N_IO_ParsingHelpers.h>
#include <N_PDS_Comm.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Stats.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::DistToolFlatRoundRobin
// Purpose       : ctor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistToolFlatRoundRobin::DistToolFlatRoundRobin(
  Parallel::Communicator *                 pdsCommPtr,
  CircuitBlock &                           circuit_block,
  std::map<std::string,FileSSFPair>      & ssfMap,
  std::map<std::string, IncludeFileInfo> & iflMap,
  const std::vector< std::pair< std::string, std::string> > & externalNetlistParams)
  : DistToolBase(pdsCommPtr, circuit_block, ssfMap),
    iflMap_(iflMap),
    procDeviceCount_(0),
    blockSize_(0),
    numBlocks_(20),
    devices_(0),
    myDeviceLines_(0),
    bufs_(),
    externalNetlistParams_(externalNetlistParams)
{
  myProc_ = pdsCommPtr_->procID();

  // Set the circuit context and options table.
  setCircuitContext();
  setCircuitOptions();

  devices_ = (circuit_block.getCircuitContextPtr())->getTotalDeviceCount();
  pdsCommPtr_->bcast( &devices_, 1, 0 );
  procDeviceCount_ = devices_ / numProcs_;

  blockSize_ = procDeviceCount_ / numBlocks_;

  // have 500 devices as the minimum size block.
  if (blockSize_ < 500)
  {
    numBlocks_ = procDeviceCount_ / 500;
    if (numBlocks_ == 0 || numBlocks_ == 1)
    {
      blockSize_ = procDeviceCount_;
      numBlocks_ = 1;
    }
    else
    {
      blockSize_ = procDeviceCount_ / numBlocks_ + 1;
    }
  }

  // Get each processors neighbor.
  // Any remaining devices should be processed by the first processor.
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    if (myProc_ == 0)
    {
      int rem = devices_ - numProcs_*procDeviceCount_;
      procDeviceCount_ += rem;
      if (numBlocks_ == 1)
        blockSize_ += rem;
    }
 
    fromProc_ = myProc_-1;
    if (fromProc_ < 0)
      fromProc_ = numProcs_-1; 
    toProc_ = (myProc_+1)%numProcs_;
  }
  else
  {
    numBlocks_ = 1;
    blockSize_ = devices_;
  }

  // Resize the buffer here, just in case the blocksize is changed above. 
  bufs_.resize( blockSize_ );

  //std::cout << "Processor " << myProc_ << ", procDeviceCount_ = " << procDeviceCount_ << ", blockSize_ = " << blockSize_ << std::endl;
}


//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::endDeviceLines()
// Purpose       : Make sure other processors have been told that device processing
//                 is ended.  This can be a problem with small device count circuits
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/06/06
//-----------------------------------------------------------------------------
void DistToolFlatRoundRobin::endDeviceLines(int stopProc)
{
  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    //std::cout << "Processor " << myProc_ << " is sending the finish signal to processor " << toProc_ << std::endl;
    // Tell neighbor processor to stop parsing.
    int minus2 = -2;
    pdsCommPtr_->send( &minus2, 1, toProc_ );

    // Send the neighbor processor which processor requested the parsing stop.
    pdsCommPtr_->send( &stopProc, 1, toProc_ );
  }
}


//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::setFileName
// Purpose       : Change name of netlist file
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DistToolFlatRoundRobin::setFileName(std::string const & fileNameIn)
{
  netlistFilename_ = fileNameIn;
  circuitBlock_.setFileName( fileNameIn );
}

//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::distributeDevices()
// Purpose       : Distribute the devices according to the strategy.
// Special Notes : The strategy performs a round-robin assignment of devices to
//               : processors to interleave communication and computation.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 08/28/2014
//-----------------------------------------------------------------------------
void DistToolFlatRoundRobin::distributeDevices()
{
  Xyce::Parallel::Machine comm = pdsCommPtr_->comm();

  if (DEBUG_IO) {
    Xyce::dout() << "Pass 2 parsing of netlist file and distributing devices: " <<
      circuitBlock_.getNetlistFilename() << std::endl;
  }

  bool ok = true;

  // Let processor zero take the first block of devices.
  if (myProc_ == 0)
  {
    setFileName(circuitBlock_.getNetlistFilename());

    ssfPtr_ = ssfMap_[netlistFilename_].second;
    ssfPtr_->setLocation( circuitBlock_.getStartPosition() );
    ssfPtr_->setLineNumber( circuitBlock_.getLineStartPosition() );

    // Skip over the title line and continue.
    std::ifstream* netlistIn = ssfMap_[netlistFilename_].first;
    std::string title("");
    netlistIn->clear();
    netlistIn->seekg(0, std::ios::beg);
    Xyce::IO::readLine( *netlistIn, title ); 
    ssfPtr_->changeCursorLineNumber( 1 );

    int startIdx = 0;

    // Process MIs if present in top circuit
    if( circuitContext_->haveMutualInductances() )
    {
      // Distribute each YMI?|name device line, end of circuit
      int n = circuitContext_->getNumMILines();
      if (n > blockSize_)
      {
        bufs_.resize( n );
      }
      for( int i = 0; i < n; ++i )
      {
        bufs_[i] = circuitContext_->getMILine( i );
/*        std::cout << "Parsed MI line on processor " << myProc_ << " ";
        TokenVector::iterator currField = bufs_[i].begin();
        TokenVector::iterator endField = bufs_[i].end();
        while( currField != endField )
        {
          std::cout << currField->string_ << " ";
          currField++;
        } 
        std::cout << std::endl;
*/
      }
      myDeviceLines_ += n;
    
      if ( n >= blockSize_ )
      {  
        processDeviceBuffer();
        bufs_.resize( blockSize_ );
      }
      else
      {
        startIdx = n;
      }
    }

    // Get the first block of devices from the netlist.
    std::string libSelect;
    std::vector<std::string> libInside;
    ok = bufferDeviceLines(netlistFilename_, libSelect, libInside,
                           ssfPtr_->getFilePosition(), ssfPtr_->getLineNumber(),
                           myDeviceLines_, startIdx);

    if (bufs_.size() > 0)
    {
      processDeviceBuffer();
    }
  }

  // Wait to receive the signal to buffer data from the ifstream.
  // NOTE: Every processor but zero will enter this function immediately.
  //       If processor 0 encounters an error and communicates that error
  //       in the bufferDeviceLines() it will not enter this function.
  if (ok)
  {
    waitYourTurn();
  }

  if (DEBUG_IO)
    Xyce::dout() << "Done with pass 2 netlist file parsing and device distribution" << std::endl;

  // Just in case an error was reported in parsing the netlist.
  N_ERH_ErrorMgr::safeBarrier(comm);
}

//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::bufferCircuitData
// Purpose       : place next block of data in device buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolFlatRoundRobin::bufferDeviceLines(std::string fileName,
                                               std::string libSelect, 
                                               std::vector<std::string>& libInside,
                                               long filePos, long lineNum,
                                               int totalDevCount, int startIdx)
{
  bool ok=true;

  //std::cout << "Processor " << myProc_ << " is now buffering device lines!" << std::endl;
  //std::cout << "totalDevCount = " << totalDevCount << " of " << devices_ << " and " << myDeviceLines_ << " of " << procDeviceCount_ << std::endl;
  if (totalDevCount < devices_ && myDeviceLines_ < procDeviceCount_)
  {

    // Open this file, if it is not open already.
    if( !ssfMap_.count( fileName ) )
    {
      std::ifstream * in = new std::ifstream;

      // Using binary to avoid issues with compiler/plat *fstream differences in implementation
      in->open( fileName.c_str(), std::ios::in | std::ios::binary );
      if ( !in->is_open() )
      {
        Report::UserError() << "Could not open file " << fileName << " on processor "
                            << pdsCommPtr_->procID() << std::endl;
        ok=false;
      }
      ssfMap_[netlistFilename_] = FileSSFPair( in, new SpiceSeparatedFieldTool(*in, netlistFilename_, externalNetlistParams_) );
      //std::cout << "Opened file " << netlistFilename_ << " on processor " << pdsCommPtr_->procID() << std::endl;
    }
    ssfPtr_ = ssfMap_[fileName].second;
    ssfPtr_->setLocation( filePos );
    ssfPtr_->setLineNumber( lineNum );
 
    int numLines = startIdx;
    while (ok && numLines < blockSize_)
    {
      ok = getLine(bufs_[numLines], libSelect, libInside);
      if (ok && !bufs_[numLines].empty() && compare_nocase((bufs_[numLines])[0].string_.c_str(), ".ends") != 0)
      {
/*
        TokenVector::iterator currField = bufs_[numLines].begin();
        TokenVector::iterator endField = bufs_[numLines].end();
        std::cout << "Parsed line on processor " << myProc_ << " ";
        while( currField != endField )
        {
          std::cout << currField->string_ << " ";
          currField++;
        } 
        std::cout << std::endl;
*/        
        numLines++;
      }
    }
    // Resize the buffer before processing.
    bufs_.resize( numLines );

    // Keep a total of the number of devices processed on this processor.
    myDeviceLines_ += (numLines-startIdx);
    totalDevCount += (numLines-startIdx);
    //std::cout << "Processor " << myProc_ << " has buffered " << numLines << " device lines!" << std::endl;
  }
  else
  {
    // This processor is not doing any work.
    bufs_.clear();
  }
 
  // If there was a file error or all the devices have been read, end the parsing.
  if (!ok || totalDevCount >= devices_)
  {
    //std::cout << "Processor " << myProc_ << ", ok = " << ok << ", totalDevCount = " << totalDevCount << ", devices_ = " << devices_ << std::endl;
    // Tell neighbor processor to stop parsing.
    endDeviceLines( myProc_ );

    // Make sure to process the current buffers on this processor.
    if (bufs_.size() > 0)
    {
      processDeviceBuffer();
    }
    
    // Now clear buffer.
    bufs_.clear();

    // Make sure the exit is clean.
    ok=false;
  }
  else
  {
    if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
    {
      // Everything looks good to continue.
      // Send the next processor the information.
      int minus1 = -1;
      pdsCommPtr_->send( &minus1, 1, toProc_ );

      //std::cout << "Processor " << myProc_ << " is sending a signal to processor " << toProc_ << " to buffer device lines." << std::endl;
      // Send the current filename, file position, line number, and total devices processed.
      // This enables the next processor to start processing from where this processor left off.
      int length = fileName.size();
      long pos = ssfPtr_->getFilePosition();
      long num = ssfPtr_->getLineNumber();
      //std::cout << "Processor " << myProc_ << " sending netlist name " << netlistFilename_ << " pos = " << pos << " num = " << num << " totalDevCount = " << totalDevCount << std::endl;
      pdsCommPtr_->send( &length, 1, toProc_ );
      pdsCommPtr_->send( netlistFilename_.c_str(), length, toProc_ );
      pdsCommPtr_->send( &pos, 1, toProc_ );
      pdsCommPtr_->send( &num, 1, toProc_ );
      pdsCommPtr_->send( &totalDevCount, 1, toProc_ );
    }
  }

  return ok;
}


//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::waitYourTurn
// Purpose       : receive all data from proc 0 and store it in device buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistToolFlatRoundRobin::waitYourTurn()
{
  bool ok = true;

  if (Parallel::is_parallel_run(pdsCommPtr_->comm()))
  {
    int signal = 0;
    //Stats::StatTop _processDeviceStat("Process Devices");

    // The outer loop controls the number of blocks of device buffers that are sent by proc 0.
    while (true)
    {
      // The inner loop controls the number of buffers necessary to send the current block of devices from proc 0.
      while (true)
      {
        // Wait until we get the a signal from our neighbor processor.
        pdsCommPtr_->recv( &signal, 1, fromProc_ );

        // If the signal is -1, then get the next block of device lines from the file.
        if (signal == -1)
        {
           //std::cout << "Processor " << myProc_ << " just received the signal from processor " << fromProc_ << " to buffer device lines." << std::endl;
           // Get the current filename, file position, line number, and total devices processed.
           // This enables the this processor to start processing from where the last processor left off.
           int length = 0, totalDevCount = 0;
           long pos = 0, num = 0;
           pdsCommPtr_->recv( &length, 1, fromProc_ );
           char * file = new char[ length ];
           pdsCommPtr_->recv( file, length, fromProc_ );
           netlistFilename_ = std::string( file, length );
           pdsCommPtr_->recv( &pos, 1, fromProc_ );
           pdsCommPtr_->recv( &num, 1, fromProc_ );
           pdsCommPtr_->recv( &totalDevCount, 1, fromProc_ );

           //std::cout << "Processor " << myProc_ << " recieving netlist name " << netlistFilename_ << " pos = " << pos << " num = " << num << " totalDevCount = " << totalDevCount << std::endl;
           // Now we can continue parsing.
           std::string libSelect;
           std::vector<std::string> libInside;
           ok = bufferDeviceLines( netlistFilename_, libSelect, libInside,
                                   pos, num, totalDevCount );
        
           // Return to the inner loop to process this buffer and wait for the next block of devices.
           break; 
        }

        // If the signal is -2, then netlist processing is complete.  This processor should exit the loop.  
        if (signal == -2)
        {
          int stopProc = 0;
          pdsCommPtr_->recv( &stopProc, 1, fromProc_ ); 
          //std::cout << "Processor " << myProc_ << " just received the stop signal from processor " << fromProc_ << " to end processing device lines." << std::endl;
 
          // Only forward the message if this processor didn't send it.
          if (stopProc != myProc_ && stopProc != toProc_)
          {
            endDeviceLines( stopProc );
          }

          // Processing has ended, resize the buffer to 0.
          bufs_.clear();
          break;
        }
      }  
   
      // Netlist parsing has concluded, time to exit. 
      if (!ok || signal == -2)
      {
        break;
      }

      {
        //Stats::TimeBlock _processDeviceTimer(_processDeviceStat);
        processDeviceBuffer();
      }
    
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::processDeviceBuffer
// Purpose       : process all data in the buffer
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
DistToolFlatRoundRobin::processDeviceBuffer()
{
  // process received device lines, contexts info, and exit calls
  for (unsigned int ib=0 ; ib<bufs_.size() ; ++ib)
  {
    // hand to circuit block for processing
    std::string libSelect("");
    std::vector<std::string> libInside; 
    handleDeviceLine( bufs_[ib], libSelect, libInside );
    bufs_[ib].clear();
  }
  //std::cout << "Processor " << myProc_ << " just processed " << bufs_.size() << " device lines." << std::endl;

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::broadcastGlobalData
// Purpose       : send bufsize, options, and context to procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool
DistToolFlatRoundRobin::broadcastGlobalData()
{
  return DistToolBase::broadcastGlobalData();
}

//--------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::parseIncludeFile
// Purpose       : Jump to include file for 2nd pass
// Special Notes :
// Creator       : Rob Hoekstra
// Creation Date : 08/02/2004
//--------------------------------------------------------------------------
bool DistToolFlatRoundRobin::parseIncludeFile(std::string const& includeFile,
             const std::string &libSelect)
{
/*
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
*/
  return true; // Only get here on success.
}

//--------------------------------------------------------------------------
// Function      : DistToolFlatRoundRobin::restorePrevssfInfo
// Purpose       : This is a helper function for parseIncludeFile(). It
//                 restores the information about the previous file.   It
//                 should be called before each return statement in that
//                 function.
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 11/11/2020
//--------------------------------------------------------------------------
void DistToolFlatRoundRobin::restorePrevssfInfo(
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

} // namespace IO
} // namespace Xyce

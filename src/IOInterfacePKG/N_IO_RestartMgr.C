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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 7/19/01
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <sstream>
#include <fstream>
#include <algorithm>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_RestartMgr.h>
#include <N_IO_RestartNode.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TOP_Topology.h>
#include <N_TIA_OneStep.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {

#ifdef Xyce_RESTART_NOPACK
static const bool RESTART_NOPACK = true;
#else
static const bool RESTART_NOPACK = false;
#endif

//-----------------------------------------------------------------------------
// Function      : RestartMgr::RestartMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
RestartMgr::RestartMgr()
  : restartFlag_(false),
    restartFileName_(""),
    restartJobName_(""),
    initialRestartInterval_(0.0),
    pack_(!RESTART_NOPACK),
    printOptions_(false)
{}

//-----------------------------------------------------------------------------
// Function      : RestartMgr::registerRestartOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool RestartMgr::registerRestartOptions(const Util::OptionBlock & option_block, int proc_size, int proc_rank)
{
  bool success = true;
  double start_time = 0.0;

  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    ExtendedString tag(it->tag());
    tag.toUpper();

    if (tag == "TIME")
    {
      double t = it->getImmutableValue<double>();
      ++it;
      if( it == option_block.end())
      {
        Report::UserError() << ".OPTIONS Restart: Badly formed INITIAL_INTERVAL\n";
        success = false;
      }
      double iv = it->getImmutableValue<double>();
      restartIntervalPairs_.push_back(Interval(t,iv));
    }

    else if (tag == "PACK")
    {
      pack_ = it->getImmutableValue<int>();

      if (!pack_ && proc_size != 1)
      {
        Report::UserError() << ".OPTIONS Restart: Restart Data Must Be Packed for Parallel Runs, Please remove 'pack=0' from .OPTIONS RESTART line";
        return false;
      }
    }

    else if (tag == "INITIAL_INTERVAL")
    {
      initialRestartInterval_ = it->getImmutableValue<double>();
    }

    else if (tag == "FILE")
    {
      restartFileName_ = it->stringValue();
      restartFlag_ = true;
    }

    else if (tag == "JOB")
    {
      restartJobName_ = it->stringValue();
    }

    else if (tag == "START_TIME")
    {
      start_time = it->getImmutableValue<double>();
      restartFlag_ = true;
    }

    else if (tag == "PRINT_TIMEINT_OPTIONS")
    {
      printOptions_ = it->getImmutableValue<int>();
    }
  }

  if (restartFlag_ && restartFileName_.empty())
  {
    std::ostringstream oss;
    oss << restartJobName_ << start_time;
    restartFileName_ = oss.str();
  }

  if (pack_ && RESTART_NOPACK)
  {
    pack_ = false;
    Report::UserWarning() << "Restart packing disabled for this implementation of Xyce";
  }

  if (!pack_  && proc_size != 1)
  {
    Report::UserError() << ".OPTIONS Restart: Parallel Restarts must be packed!";
    success = false;
  }

  if (DEBUG_RESTART && proc_rank == 0)
  {
    Xyce::dout() << Xyce::subsection_divider << std::endl
                 << "RESTART OPTIONS SETUP: " << std::endl
                 << "restartFileName: " << restartFileName_ << std::endl
                 << "restartJobName: " << restartJobName_ << std::endl
                 << "isRestart: " << restartFlag_ << std::endl
                 << "initial Interval: " << initialRestartInterval_ << std::endl
                 << "pack data: " << pack_ << std::endl
                 << "print time integration options: " << printOptions_ << std::endl;

    for (unsigned int i = 0; i < restartIntervalPairs_.size(); ++i)
      Xyce::dout() << restartIntervalPairs_[i].first << " " << restartIntervalPairs_[i].second << std::endl;

    Xyce::dout() << Xyce::subsection_divider << std::endl;
  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : RestartMgr::registerTimeintOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool RestartMgr::registerTimeintOptions(const Util::OptionBlock & option_block)
{
  // Save the time integrator option block
  savedTimeintOB_ = option_block;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : dumpRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool dumpRestartData(
  Parallel::Communicator &      comm,
  Topo::Topology &              topology,
  Analysis::AnalysisManager &   analysis_manager,
  Device::DeviceMgr &           device_manager,
  const std::string &           job_name,
  bool                          pack,
  double                        time)
{
  bool success = true;

  int procID = comm.procID();
  int numProcs = comm.numProc();

  //Get current time integration method
  int method = analysis_manager.getIntegrationMethod();

  //Get Restart Node Data and setup buffers
  //---------------------------------------
  std::vector<RestartNode *> nodeVec;

  // Get restart nodes and compute global node count
  topology.getRestartNodes(analysis_manager, nodeVec);

  int nodeCount = nodeVec.size();
  double nC = static_cast<double>(nodeCount);

  double gNC = 0;
  comm.sumAll( &nC, &gNC, 1);

  int globalNodeCount = static_cast<int>(gNC);

  // Compute maxSize to allocate data buffer.
  std::vector<int> dataSizeVec;

  int analysis_dataSize = 0;
  int device_dataSize = 0;
  int node_dataSize = sizeof(int);

  if (procID == 0)
  {
    analysis_dataSize = analysis_manager.getRestartDataSize(pack);
    device_dataSize = device_manager.restartDataSize(pack);
 
    dataSizeVec.push_back( analysis_dataSize );
    dataSizeVec.push_back( device_dataSize );
  }

  for( int i = 0; i < nodeCount; ++i)
    node_dataSize += Xyce::packedByteCount(*nodeVec[i]);

  int dataSize = node_dataSize;

  if (DEBUG_RESTART)
    Xyce::dout() << "fProcID: " << procID << std::endl
                 << "dataSize: " << node_dataSize << std::endl;
      
  if (procID == 0)
  {
    dataSizeVec.push_back( node_dataSize );
    std::vector<int>::iterator maxIter = std::max_element( dataSizeVec.begin(), dataSizeVec.end() );
    dataSize = *maxIter;
  }

  int maxSize = 0;
  int bsize;
  char * buf;
#ifdef Xyce_PARALLEL_MPI
  comm.maxAll( &dataSize, &maxSize, 1);
  if (procID == 0)
  {
    buf = new char[maxSize];
    bsize = maxSize;
  }
  else
#endif
  {
    buf = new char[dataSize];
    bsize = dataSize;
    maxSize = dataSize;
  }

  std::string versionString = Xyce::Util::Version::getShortVersionString();
  std::ofstream outStreamSt;
  std::string outName;

  int proc = 0;
  if (procID == 0)
  {
    std::ostringstream ost;
    ost << job_name << time;
    outName.assign(ost.str());
    outStreamSt.open( outName.c_str());

    if (!outStreamSt.is_open())
    {
      Report::UserWarning0() << "Cannot Open CheckPoint File: " << outName;

      return false;
    }

    if (DEBUG_RESTART)
      Xyce::dout() << Xyce::subsection_divider << std::endl
                   << "DUMPING RESTART: " << outName << " version: " << versionString << std::endl
                   << "numProcs: " << numProcs << " maxSize: " << maxSize << std::endl
                   << "proc: " << proc << " dataSize: " << dataSize << std::endl
                   << "nodeCount: " << globalNodeCount << " pack: " << pack << std::endl;

    int packed = pack?1:0;
    outStreamSt << numProcs << " " << maxSize << " " << packed << " " 
                << globalNodeCount << " " << method << " " << versionString << " ";
  }

  // 1) Add Time Int stuff from proc 0
  //------------------------------
  if (procID == 0)
  {
    int pos = 0;
    success &= analysis_manager.dumpRestartData( buf, bsize, pos, &comm, pack);

    outStreamSt << analysis_dataSize << " ";
    outStreamSt.write( buf, analysis_dataSize);
  }

  comm.barrier();

  // 2) Pack Restart Nodes
  if (pack)
  {
    int pos = 0;
    comm.pack( &nodeCount, 1, buf, bsize, pos);
    for( int i = 0; i < nodeCount; ++i)
    {
      Xyce::pack(*nodeVec[i], buf, bsize, pos, &comm);
    }

    if (procID == 0)
    {
      outStreamSt << procID << " " << node_dataSize << " ";
      outStreamSt.write( buf, node_dataSize);
    }
  }
  else
  {
    outStreamSt << " " << nodeCount << " ";
    for( int i = 0; i < nodeCount; ++i)
    {
      nodeVec[i]->dump(outStreamSt);
    }
  }

#ifdef Xyce_PARALLEL_MPI
  //PARALLEL, get and store nodes from other procs
  //----------------------------------------------
  comm.barrier();

  int size;
  for( proc = 1; proc < numProcs; ++proc)
  {
    if (procID == 0)
    {
      comm.recv( &size, 1, proc);
      comm.recv( buf, size, proc);

      outStreamSt << proc << " " << size << " ";
      outStreamSt.write( buf, size);
    }
    else if (procID == proc)
    {
      comm.send( &dataSize, 1, 0);
      comm.send( buf, dataSize, 0);
    }

    comm.barrier();
  }
#endif

  if (DEBUG_RESTART && procID == 0)
    Xyce::dout() << Xyce::subsection_divider << std::endl;

  // 3) Add Device Mgr stuff from proc 0
  //------------------------------
  if (procID == 0)
  {
    int pos = 0;
    success &= device_manager.dumpRestartData( buf, bsize, pos, &comm, pack);

    outStreamSt << device_dataSize << " ";
    outStreamSt.write( buf, device_dataSize );
  }

  comm.barrier();

  if (DEBUG_RESTART && procID == 0)
    Xyce::dout() << Xyce::subsection_divider << std::endl;

  int nodeSize = nodeVec.size();
  for (int i = 0; i < nodeSize; ++i)
    delete nodeVec[i];

  delete[] buf;
  bsize = 0;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : restoreRestartData
// Purpose       :
// Special Notes : version used with distributed parser
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool
RestartMgr::restoreRestartData(
  Parallel::Communicator &      comm,
  Topo::Topology &              topology,
  Analysis::AnalysisManager &   analysis_manager,
  Device::DeviceMgr &           device_manager,
  const std::string &           path)
{
  bool success = true;

  int procID = comm.procID();
  int numProcs = comm.numProc();

  char * buf = 0;
  int dataSize;
  int pos = 0;
  int bsize = 0;
  int fProcID;

  std::ifstream * inStream = 0;

  int nodeCount;
  int totNumNodes;
  int oldNumProcs, maxSize;
  std::string oldVersion;

  bool pack = RESTART_NOPACK ? false : numProcs > 1;

  // Get requested time integration method and set it in analysis manager
  // before restoring analysis manager information.
  // NOTE:  This might not be specified in the option block, so use default if it hasn't.
  int integrationMethod = TimeIntg::OneStep::type;
  for (Util::ParamList::const_iterator it = savedTimeintOB_.begin(), end = savedTimeintOB_.end(); it != end; ++it)
  {
    const Util::Param &param = (*it);

    if (param.uTag() == "METHOD")
    {
      if (param.isInteger())
        integrationMethod = param.getImmutableValue<int>();
      else
      {
        ExtendedString stringVal ( param.stringValue() );
        stringVal.toUpper();

        if (stringVal == "TRAP" || stringVal == "TRAPEZOIDAL")
          integrationMethod = 7;
        else if (stringVal == "GEAR")
          integrationMethod = 8;
        else
        {
          IO::ParamError(savedTimeintOB_, param) << "RESTART requested: Unsupported time integration method: " << stringVal;
          return false;
        }
      }
    }
  }
  analysis_manager.setIntegrationMethod(integrationMethod); 

  int oldIntMethod;

  //1. Proc 0 reads in analysis data and distribures to every processor
  if (procID == 0)
  {
    // Error out if the user-specified checkpoint file does not exist, cannot be opened,
    // or is a directory name rather than a file name.  See SON Bugs 730 
    // and 785 for more details.  
    if ( !(Util::checkIfValidFile(path)) )
    {
      Report::UserFatal0() << "Cannot Find CheckPoint File: " << path << " for Restart";
      return false;
    }

    inStream = new std::ifstream( path.c_str());

    if (!inStream->is_open())
    {
      Report::UserFatal0() << "Cannot Open CheckPoint File: " << path << " for Restart";
      return false;
    }

    int packed;
    (*inStream) >> oldNumProcs >> maxSize >> packed >> totNumNodes >> oldIntMethod >> oldVersion;
    pack = packed;

    if (DEBUG_RESTART)
      Xyce::dout() << Xyce::subsection_divider << std::endl
                   << "RESTORING RESTART DATA" << std::endl
                   << "oldNumProcs: " << oldNumProcs << std::endl
                   << "maxSize: " << maxSize << std::endl
                   << "packed: " << pack << std::endl
                   << "totNumNodes: " << totNumNodes << std::endl
                   << "oldIntMethod: " << oldIntMethod << std::endl
                   << "oldVersion: " << oldVersion << std::endl;

    // read in Time data
    (*inStream) >> dataSize;
    char dummy = 'x';
    inStream->read( &dummy, 1);

    if (DEBUG_RESTART)
      Xyce::dout() << Xyce::subsection_divider << std::endl
                   << "RESTORING TIME INT STUFF\n" << std::endl
                   << "dataSize: " << dataSize << std::endl
                   << Xyce::subsection_divider << std::endl;

    maxSize = dataSize;
    buf = new char[maxSize];
    bsize = maxSize;

    inStream->read( buf, dataSize);
  }

#ifdef Xyce_PARALLEL_MPI
  comm.barrier();
  comm.bcast( &oldIntMethod, 1, 0);
  comm.bcast( &dataSize, 1, 0);

  if (procID != 0)
  {
    buf = new char[dataSize];
    bsize = dataSize;
  }

  comm.bcast( buf, dataSize, 0);
#endif

  pos = 0;
  analysis_manager.restoreRestartData(buf, bsize, pos, &comm, pack);

  // If the time integration method has changed, treat restart as though it happened
  // at a breakpoint since history information is tied to integration method.
  //std::cout << "oldIntMethod: " << oldIntMethod << std::endl;

  if (oldIntMethod != analysis_manager.getIntegrationMethod())
  {
    analysis_manager.setBeginningIntegrationFlag( true );
  }

  //Table of restart nodes used by all procs
  std::vector<RestartNode*> nodeTable;

  //2. Proc 0 reads in data pushing NumDevices/NumProcs to every processor
  if (procID == 0)
  {
    char * buf1 = 0;
    int    bsize1 = 0;
    int    pos1   = 0;

    unsigned int DevsPerProc = totNumNodes/numProcs;
    if (totNumNodes % numProcs)
      DevsPerProc++;

    RestartNode * nodeP;
    int CurrProc = 1;
    for( int oldProc = 0; oldProc < oldNumProcs; ++oldProc)
    {
      if (pack)
      {
        (*inStream) >> fProcID >> dataSize;
        char dummy = 'x';
        inStream->read( &dummy, 1);

        if (DEBUG_RESTART)
          Xyce::dout() << "fProcID: " << fProcID << std::endl
                       << "dataSize: " << dataSize << std::endl;

        bsize = dataSize;
        buf = new char[bsize];

        inStream->read( buf, dataSize);

        pos = 0;
        comm.unpack( buf, bsize, pos, &nodeCount, 1);
      }
      else
        (*inStream) >> nodeCount;

      if (DEBUG_RESTART)
        Xyce::dout() << "nodeCount: " << nodeCount << std::endl;

      for( int i = 0; i < nodeCount; ++i)
      {
        nodeP = new RestartNode();

        if (pack)
          Xyce::unpack(*nodeP, buf, bsize, pos, &comm);
        else
          nodeP->restore(*inStream);

        nodeTable.push_back( nodeP);

        if ((nodeTable.size()==DevsPerProc) ||
             ( (i==nodeCount-1) && (oldProc==oldNumProcs-1)))
        {
          //Pack up data and send to CurrProc
          int Size = 0;
          int nTSize = nodeTable.size();
          Size += sizeof(int);
          for( int j = 0; j < nTSize; ++j)
            Size += Xyce::packedByteCount(*nodeTable[j]);

          if (Size > bsize1)
          {
            if (bsize1) delete [] buf1;
            bsize1 = Size;
            buf1 = new char[bsize1];
          }

          pos1 = 0;
          comm.pack( &nTSize, 1, buf1, bsize1, pos1);
          for( int j = 0; j < nTSize; ++j)
            Xyce::pack(*nodeTable[j], buf1, bsize1, pos1, &comm);

          if (CurrProc < numProcs)
          {
            if (DEBUG_RESTART)
              Xyce::dout() << "Sending restart nodes to processor " << CurrProc 
                           << ", number nodes sent " << i << " of " << nodeCount << std::endl;

            comm.send( &Size, 1, CurrProc);
            comm.send( buf1, Size, CurrProc);
            ++CurrProc;
          }

          for_each( nodeTable.begin(), nodeTable.end(), DeletePtr<RestartNode>());
          nodeTable.clear();
        }
      }
    }
   
    // If there are processors that weren't assigned nodes, send buffer size 0. 
    for ( int i=CurrProc; i<numProcs; ++i )
    { 
      int Size=0;
      comm.send( &Size, 1, i );
    }

    if (bsize) delete [] buf;
    buf = buf1;
    bsize = bsize1;
  }

#ifdef Xyce_PARALLEL_MPI
  //Everybody else receive their set of devices
  if (procID != 0)
  {
    comm.recv( &bsize, 1, 0);
    if (DEBUG_RESTART)
      Xyce::dout() << "Proc " << procID << " preparing to receive restart nodes, bsize: " << bsize << std::endl;
    
    if (bsize)
    {
      buf = new char[bsize];
      comm.recv( buf, bsize, 0);
    }
    else
    {
      // Place a zero in the buffer because this processor was not allocated nodes
      bsize = sizeof(int);
      buf = new char[bsize];
      int pos = 0;
      nodeCount = 0;
      comm.pack( &nodeCount, 1, buf, bsize, pos);
    } 
  }

  if (DEBUG_RESTART)
    Xyce::dout() << "Finished ReadIn of Restart Node Data on Proc: " << procID << std::endl;

  comm.barrier();

  //Set Buffers to largest necessary size
  int bsize2 = 0;
  comm.maxAll( &bsize, &bsize2, 1);
  char * buf2 = new char[bsize2];
  memcpy( buf2, buf, bsize);

  delete [] buf;
  bsize = bsize2;
  buf = new char[bsize];

  std::swap( buf, buf2);
#endif

  //Check devices in your buffer, then push around ring
  for( int i = 0; i < numProcs; ++i)
  {
    pos = 0;
    comm.unpack( buf, bsize, pos, &nodeCount, 1);
    nodeTable.resize(nodeCount);

    if (DEBUG_RESTART)
      Xyce::dout() << "Restart Ring Stage: " << i << " " << nodeCount << std::endl;

    for( int j = 0; j < nodeCount; ++j)
    {
      nodeTable[j] = new RestartNode();
      Xyce::unpack(*nodeTable[j], buf, bsize, pos, &comm);
    }

    topology.restoreRestartNodes(analysis_manager, nodeTable);

    for_each( nodeTable.begin(), nodeTable.end(), DeletePtr<RestartNode>());
    nodeTable.clear();

#ifdef Xyce_PARALLEL_MPI
    int sendProc = procID-1;
    if (sendProc == -1) sendProc = numProcs-1;
    int recvProc = procID+1;
    if (recvProc == numProcs) recvProc = 0;

    int MaxSendSize = 1e6;

    int NumSends = 1;
    int SendSize = bsize;
    if (bsize > MaxSendSize)
    {
      NumSends = bsize/MaxSendSize;
      if (bsize != NumSends*MaxSendSize) ++NumSends;
      SendSize = MaxSendSize;
    }

    if (DEBUG_RESTART)
      Xyce::dout() << "NumSends: " << NumSends << std::endl;

    for( int j = 0; j < NumSends; ++j)
    {
      int Offset = j*SendSize;

      if (j == (NumSends-1)) SendSize = bsize - Offset;

      comm.iRecv( buf2+Offset, SendSize, sendProc);

      comm.rSend( buf+Offset, SendSize, recvProc);

      comm.waitAll();
    }

    std::swap( buf, buf2);
#endif
  }

#ifdef Xyce_PARALLEL_MPI
  bsize2 = 0;
  delete [] buf2;
  buf2 = 0;
#endif

// read in Dev data
  if (procID == 0)
  {
    (*inStream) >> dataSize;
    char dummy = 'x';
    inStream->read( &dummy, 1);

    if (DEBUG_RESTART)
      Xyce::dout() << Xyce::subsection_divider << std::endl
                   << "RESTORING DEV MGR STUFF" << std::endl
                   << "dataSize: " << dataSize << std::endl
                   << Xyce::subsection_divider << std::endl;

    if (dataSize > bsize)
    {
      delete [] buf;
      bsize = dataSize;
      buf = new char[bsize];
    }

    inStream->read( buf, dataSize);

    delete inStream;
  }

#ifdef Xyce_PARALLEL_MPI
  comm.barrier();

  comm.bcast( &dataSize, 1, 0);
  if (procID != 0)
    if (dataSize > bsize)
    {
      delete [] buf;
      buf = new char[dataSize];
      bsize = dataSize;
    }

  comm.bcast( buf, dataSize, 0);
#endif

  pos = 0;
  device_manager.restoreRestartData( buf, bsize, pos, &comm, pack);
  delete [] buf;

  return success;
}

namespace {

struct OptionsReg : public PkgOptionsReg
{
  OptionsReg(RestartMgr &restart_manager, int proc_size, int proc_rank)
    : restartManager_(restart_manager),
      size_(proc_size),
      rank_(proc_rank)
  {}

  bool operator()(const Util::OptionBlock & options)
  {
    return restartManager_.registerRestartOptions(options, size_, rank_);
  }

  RestartMgr &                  restartManager_;
  const int                     size_;
  const int                     rank_;
};

struct TimeintOptionsReg : public PkgOptionsReg
{
  TimeintOptionsReg(RestartMgr &restart_manager)
    : restartManager_(restart_manager)
  {}

  bool operator()(const Util::OptionBlock & options)
  {
    return restartManager_.registerTimeintOptions(options);
  }

  RestartMgr &                  restartManager_;
};

void
populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("RESTART");

  parameters.insert(Util::ParamMap::value_type("PACK", Util::Param("PACK", 1)));
  parameters.insert(Util::ParamMap::value_type("PRINT_TIMEINT_OPTIONS", Util::Param("PRINT_TIMEINT_OPTIONS", 0)));
  parameters.insert(Util::ParamMap::value_type("JOB", Util::Param("JOB", "")));
  parameters.insert(Util::ParamMap::value_type("START_TIME", Util::Param("START_TIME", 0.0)));
  parameters.insert(Util::ParamMap::value_type("FILE", Util::Param("FILE", "")));
  parameters.insert(Util::ParamMap::value_type("INITIAL_INTERVAL", Util::Param("INITIAL_INTERVAL", 0.0)));
  parameters.insert(Util::ParamMap::value_type("TIME", Util::Param("TIME", 0.0)));
  parameters.insert(Util::ParamMap::value_type("INTERVAL", Util::Param("INTERVAL", 0.0)));
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
bool
registerPkgOptionsMgr(
  RestartMgr &          restart_manager,
  PkgOptionsMgr &       options_manager,
  int                   proc_size,
  int                   proc_rank)
{
  populateMetadata(options_manager);

  options_manager.addOptionsProcessor("RESTART", new OptionsReg(restart_manager, proc_size, proc_rank));
  options_manager.addOptionsProcessor("TIMEINT", new TimeintOptionsReg(restart_manager));

  return true;
}

} // namespace IO
} // namespace Xyce

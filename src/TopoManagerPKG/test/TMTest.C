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


#include <string>
#include <iostream>
#include <fstream>
#include <list>

#ifdef Xyce_PARALLEL_MPI
  #include <mpi.h>
#endif

#include <N_TOP_CktNode.h>
#include <N_TOP_CktNodeCreator.h>
#include <N_TOP_CktGraphSupport.h>
#include <N_TOP_CktGraphCreator.h>
#include <N_TOP_CktGraph.h>
#include <N_TOP_Topology.h>
#include <N_PDS_Comm.h>
#include <N_IO_DistribMgr.h>
#include <N_DEV_DeviceMgr.h>

int main (int argc, char* argv[]) {

#ifdef Xyce_PARALLEL_MPI
    MPI_Init( &argc, &argv );
#endif

//    N_PDS_Manager PDSMgr;

    N_PDS_Comm PDSComm;

#ifdef Xyce_PARALLEL_MPI
    N_ERH_ErrorMgr::registerComm( &PDSComm );
#endif

    N_DEV_DeviceMgr* DevMgrPtr = N_DEV_DeviceMgr::factory();

    N_TOP_Topology Topo( DevMgrPtr );

//    ParseSpiceNetlist psn ( argv[1], Topo, *DevMgrPtr, true );

//    int isuccess = psn.PSN_Read();

    N_IO_DistribMgr* DistMgrPtr = N_IO_DistribMgr::factory();
    DistMgrPtr->registerTopology( &Topo );
    DistMgrPtr->registerDeviceMgr( DevMgrPtr );

#ifdef Xyce_PARALLEL_MPI
    DistMgrPtr->registerParallelServices( &PDSComm );
#endif

    DistMgrPtr->generateParser( string( argv[1] ) );
    DistMgrPtr->run();

    cout << Topo << endl;

    Topo.OutputBFSGraphLists();

    DevMgrPtr->printOutLists();

#ifdef Xyce_PARALLEL_MPI
    MPI_Finalize();
#endif

    return 0;
}

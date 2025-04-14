//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Builder class for parallel/serial linear objects
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/22/01
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_fwd.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Graph.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_QueryUtil.h>
#include <N_LAS_Vector.h>
#include <N_LAS_SystemHelpers.h>
#include <N_PDS_Comm.h>
#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_Manager.h>
#include <N_PDS_ParHelpers.h>
#include <N_PDS_MPI.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Functors.h>

using std::max;
using Xyce::VERBOSE_LINEAR;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Function      : Builder::createMultiVector
// Purpose       : returns Soln/RHS sized multiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
MultiVector * Builder::createMultiVector( const int numVectors ) const
{
  return Xyce::Linear::createMultiVector( *(pdsMgr_->getParallelMap(Parallel::SOLUTION)),
                                          *(pdsMgr_->getParallelMap(Parallel::SOLUTION_OVERLAP_GND)),
                                          numVectors );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createStateMultiVector
// Purpose       : returns State sized multiVector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
MultiVector * Builder::createStateMultiVector( const int numVectors ) const
{
  return Xyce::Linear::createMultiVector( *(pdsMgr_->getParallelMap(Parallel::STATE)),
                                          numVectors );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createStoreMultiVector
// Purpose       : returns Store sized multiVector
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
MultiVector * Builder::createStoreMultiVector( const int numVectors ) const
{
  return Xyce::Linear::createMultiVector( *(pdsMgr_->getParallelMap(Parallel::STORE)),
                                          numVectors );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createVector
// Purpose       : returns Soln/RHS sized vector
// Special Notes : Takes an initial value argument.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 03/02/01
//-----------------------------------------------------------------------------
Vector * Builder::createVector() const
{
    return Xyce::Linear::createVector( *(pdsMgr_->getParallelMap(Parallel::SOLUTION)),
                                       *(pdsMgr_->getParallelMap(Parallel::SOLUTION_OVERLAP_GND)) );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createStateVector
// Purpose       : returns State sized vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
Vector * Builder::createStateVector() const
{
    return Xyce::Linear::createVector( *(pdsMgr_->getParallelMap(Parallel::STATE)) );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createStoreVector
// Purpose       : returns Store sized vector
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
Vector * Builder::createStoreVector() const
{
    return Xyce::Linear::createVector( *(pdsMgr_->getParallelMap(Parallel::STORE)) );
}


//-----------------------------------------------------------------------------
// Function      : Builder::createStoreVector
// Purpose       : returns Store sized vector
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
Vector * Builder::createLeadCurrentVector() const
{
    return Xyce::Linear::createVector( *(pdsMgr_->getParallelMap(Parallel::LEADCURRENT)) );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createMatrix
// Purpose       : returns Matrix initialized based on QueryUtil and ParMap
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/10/00
//-----------------------------------------------------------------------------
Matrix * Builder::createMatrix() const
{
  return Xyce::Linear::createMatrix( pdsMgr_->getMatrixGraph( Parallel::JACOBIAN_OVERLAP ),
                                     pdsMgr_->getMatrixGraph( Parallel::JACOBIAN ) );
}

//-----------------------------------------------------------------------------
// Function      : Builder::createSolnColoring
// Purpose       : Color Map representing variable types in solution vector
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/08/04
//-----------------------------------------------------------------------------
const std::vector<int> & Builder::createSolnColoring() const
{
  if (solnColoring_.empty())
  {
    const std::vector<char> & charColors = lasQueryUtil_->rowList_VarType();

    int size = charColors.size();
    solnColoring_.resize( size );

    for( int i = 0; i < size; ++i )
    {
      switch( charColors[i] )
      {
        case 'V': solnColoring_[i] = 0;
                  break;
        case 'I': solnColoring_[i] = 1;
                  break;
        default : solnColoring_[i] = 2;
                  break;
      }
    }
  }

  return solnColoring_;
}


//-----------------------------------------------------------------------------
// Function      : Builder::createInitialConditionColoring
// Purpose       :
// Special Notes : The .IC and .NODESET capabilities will use the variables
//                 which are colored 0.  This will be all voltage nodes not
//                 connected to independent sources.
// Scope         : Public
// Creator       : Eric R. Keiter,  SNL
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
const std::vector<int> & Builder::createInitialConditionColoring() const
{
  if (icColoring_.empty())
  {
    const std::vector<char> & charColors = lasQueryUtil_->rowList_VarType();
    const std::vector<int> & vsrcGIDColors = lasQueryUtil_->vsrcGIDVec();

    int size = charColors.size();
    int vsrcSize = vsrcGIDColors.size();
    icColoring_.resize( size );

    for( int i = 0; i < size; ++i )
    {
      switch( charColors[i] )
      {
        case 'V': icColoring_[i] = 0;
                  break;
        case 'I': icColoring_[i] = 1;
                  break;
        default : icColoring_[i] = 2;
                  break;
      }
    }

    Parallel::ParMap * solnMap = pdsMgr_->getParallelMap(Parallel::SOLUTION);
    for( int i=0; i < vsrcSize; ++i )
    {
      int vsrcID = vsrcGIDColors[i];
      // Convert the ID from local to global if it is valid and the build is parallel.
      if (vsrcID >= 0)
      {
        if (!pdsMgr_->getPDSComm()->isSerial()) 
          vsrcID = solnMap->globalToLocalIndex( vsrcGIDColors[i] );

        if (vsrcID < size && vsrcID >= 0)
          icColoring_[vsrcID] = 1;
      }
    }
  }

  return icColoring_;
}

//-----------------------------------------------------------------------------
// Function      : Builder::generateParMaps
// Purpose       : Creates parallel maps for SOLN, STATE, STORE, and LEAD_CURRENT
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/24/01
//-----------------------------------------------------------------------------
bool Builder::generateParMaps()
{
  int numLocalRows = lasQueryUtil_->numLocalRows();
  int numLocalStateVars = lasQueryUtil_->numLocalStateVars();
  int numLocalStoreVars = lasQueryUtil_->numLocalStoreVars();
  int numLocalLeadCurrentVars = lasQueryUtil_->numLocalLeadCurrentVars();

  if( pdsMgr_->getPDSComm()->isSerial() )
  {
    // Copy the soln GID list, since ground will be appended to it.
    std::vector<int> arrayGIDs = lasQueryUtil_->rowList_GID();
    const std::vector<int>& arrayStateGIDs = lasQueryUtil_->rowList_StateGID();
    const std::vector<int>& arrayStoreGIDs = lasQueryUtil_->rowList_StoreGID();
    const std::vector<int>& arrayLeadCurrentGIDs = lasQueryUtil_->rowList_LeadCurrentGID();
  
    Parallel::ParMap *solnMap = Parallel::createPDSParMap(numLocalRows, numLocalRows, arrayGIDs, 0, *(pdsMgr_->getPDSComm()));
    Parallel::ParMap *stateMap = Parallel::createPDSParMap(numLocalStateVars, numLocalStateVars, arrayStateGIDs, 0, *(pdsMgr_->getPDSComm()));
    Parallel::ParMap *storeMap = Parallel::createPDSParMap(numLocalStoreVars, numLocalStoreVars, arrayStoreGIDs, 0, *(pdsMgr_->getPDSComm()));
    Parallel::ParMap *leadCurrentMap = Parallel::createPDSParMap(numLocalLeadCurrentVars, numLocalLeadCurrentVars, 
                                                             arrayLeadCurrentGIDs, 0, *(pdsMgr_->getPDSComm()));

    // Create solution map with ground by appending ground to the end of the arrayGIDs vector.
    arrayGIDs.push_back(-1);
    numLocalRows++;
    Parallel::ParMap *overlapGndSolnMap = Parallel::createPDSParMap(numLocalRows, numLocalRows, arrayGIDs, -1, *(pdsMgr_->getPDSComm()));

    // Register serial maps. 
    pdsMgr_->addParallelMap( Parallel::SOLUTION, solnMap );
    pdsMgr_->addParallelMap( Parallel::STATE, stateMap );
    pdsMgr_->addParallelMap( Parallel::STORE, storeMap );
    pdsMgr_->addParallelMap( Parallel::LEADCURRENT, leadCurrentMap );
    pdsMgr_->addParallelMap( Parallel::SOLUTION_OVERLAP_GND, overlapGndSolnMap );

    // Link maps that are the same in serial and parallel because of the lack of overlap.
    pdsMgr_->linkParallelMap( Parallel::SOLUTION_OVERLAP, Parallel::SOLUTION );
  }
  else
  {
    // Create solution map, solution overlap map, and solution overlap map with ground nodes.
    std::vector<int> arrayGIDs = lasQueryUtil_->rowList_GID();

    int numGlobalRows = lasQueryUtil_->numGlobalRows();
    Parallel::ParMap *solnMap = Parallel::createPDSParMap(numGlobalRows, numLocalRows, arrayGIDs, 0, *(pdsMgr_->getPDSComm()));

    int procCnt = pdsMgr_->getPDSComm()->numProc();
    int numExternRows = lasQueryUtil_->numExternRows();
    int totalRows = numLocalRows + numExternRows;

    int size = lasQueryUtil_->rowList_ExternGID().size();
    for(int i = 0; i < size; ++i)
      arrayGIDs.push_back(lasQueryUtil_->rowList_ExternGID()[i].first);

    int numGlobalExternRows = lasQueryUtil_->numGlobalExternRows();
    int totalGlobalRows = numGlobalRows + numGlobalExternRows;
    Parallel::ParMap *overlapSolnMap = Parallel::createPDSParMap(totalGlobalRows, totalRows, arrayGIDs, 0, *(pdsMgr_->getPDSComm()));

    arrayGIDs.push_back(-1);
    totalGlobalRows += procCnt;
    totalRows++;
    Parallel::ParMap *overlapGndSolnMap = Parallel::createPDSParMap(totalGlobalRows, totalRows, arrayGIDs, -1, *(pdsMgr_->getPDSComm()));

    // Register the parallel maps.
    pdsMgr_->addParallelMap( Parallel::SOLUTION, solnMap );
    pdsMgr_->addParallelMap( Parallel::SOLUTION_OVERLAP, overlapSolnMap );
    pdsMgr_->addParallelMap( Parallel::SOLUTION_OVERLAP_GND, overlapGndSolnMap );

    // Output size of overlap.
    if (VERBOSE_LINEAR)
    {
      int procID = pdsMgr_->getPDSComm()->procID();

      std::cout << "Solution variable distribution, Processor: " << procID << " of " << procCnt << ": numLocalRows= " << numLocalRows << " of numGlobalRows= " << numGlobalRows << std::endl; 
      std::cout << "Solution variable overlap, Processor " << procID << " of " << procCnt << ": numExternRows= " << numExternRows << " of numGlobalExternRows= " << numGlobalExternRows << std::endl; 
    }

    // Create state map.
    int numGlobalStateVars = lasQueryUtil_->numGlobalStateVars();
    std::vector<int> arrayStateGIDs = lasQueryUtil_->rowList_StateGID();
    Parallel::ParMap *stateMap = Parallel::createPDSParMap(numGlobalStateVars, numLocalStateVars, 
                                                       arrayStateGIDs, 0, *(pdsMgr_->getPDSComm()));
    pdsMgr_->addParallelMap( Parallel::STATE, stateMap );

    // Create store map.
    int numGlobalStoreVars = lasQueryUtil_->numGlobalStoreVars();
    std::vector<int> arrayStoreGIDs = lasQueryUtil_->rowList_StoreGID();
    Parallel::ParMap *storeMap = Parallel::createPDSParMap(numGlobalStoreVars, numLocalStoreVars, 
                                                       arrayStoreGIDs, 0, *(pdsMgr_->getPDSComm()));
    pdsMgr_->addParallelMap( Parallel::STORE, storeMap );

    // Create lead current map and lead current overlap map.
    int numGlobalLeadCurrentVars = lasQueryUtil_->numGlobalLeadCurrentVars();
    std::vector<int> arrayLeadCurrentGIDs = lasQueryUtil_->rowList_LeadCurrentGID();
    Parallel::ParMap *leadCurrentMap = Parallel::createPDSParMap(numGlobalLeadCurrentVars, numLocalLeadCurrentVars,
                                                             arrayLeadCurrentGIDs, 0, *(pdsMgr_->getPDSComm()));
    pdsMgr_->addParallelMap( Parallel::LEADCURRENT, leadCurrentMap );
 
    // Add global accessors for parallel data migration.
    // NOTE:  This has to be after the map registration.
    Parallel::GlobalAccessor * solnGA = pdsMgr_->addGlobalAccessor( Parallel::SOLUTION );
    solnGA->registerExternGIDVector( lasQueryUtil_->rowList_ExternGID() );
    solnGA->generateMigrationPlan();
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Builder::generateGraphs
// Purpose       : Generation of Matrix Graphs, stored with ParMgr
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/02
//-----------------------------------------------------------------------------
bool Builder::generateGraphs()
{
  const std::vector< std::vector<int> > & rcData = lasQueryUtil_->rowList_ColList();
  const std::vector<int> & arrayNZs = lasQueryUtil_->rowList_NumNZs();

  int numLocalRows_Overlap = rcData.size();

  Parallel::ParMap * solnOvGMap = pdsMgr_->getParallelMap( Parallel::SOLUTION_OVERLAP_GND );
  Parallel::ParMap * solnOvMap = pdsMgr_->getParallelMap( Parallel::SOLUTION_OVERLAP );
  Parallel::ParMap * solnMap = pdsMgr_->getParallelMap( Parallel::SOLUTION );

 Graph * overlapGraph = Xyce::Linear::createGraph( *solnOvMap, *solnOvGMap, arrayNZs, rcData );

  pdsMgr_->addMatrixGraph( Parallel::JACOBIAN_OVERLAP, overlapGraph );

  if( pdsMgr_->getPDSComm()->isSerial() )
  {
    // The Jacobian graph and Jacobian overlap graph are the same in serial.
    pdsMgr_->linkMatrixGraph( Parallel::JACOBIAN, Parallel::JACOBIAN_OVERLAP );
  }
  else
  {
    Graph * solnGraph = overlapGraph->exportGraph( *solnMap );

    if (VERBOSE_LINEAR)
      Xyce::lout() << "Local Graph Transformed!\n"  << std::endl;

    pdsMgr_->addMatrixGraph( Parallel::JACOBIAN, solnGraph );
  }

  // Clean up GID arrays in lasQueryUtil
  lasQueryUtil_->cleanRowLists();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Builder::getSolutionMap
// Purpose       :
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<const Parallel::ParMap> Builder::getSolutionMap() const
{
  return(Teuchos::rcp(pdsMgr_->getParallelMap(Parallel::SOLUTION),false));
}

//-----------------------------------------------------------------------------
// Function      : Builder::getSolutionMap
// Purpose       :
// Special Notes : This is specifically for blockAnalysis types (like MPDE & HB)
// so we can get a valid map from the builder.
// Scope         : Public
// Creator       : Todd Coffey, 1414
// Creation Date : 09/05/08
//-----------------------------------------------------------------------------
RCP<Parallel::ParMap> Builder::getSolutionMap() 
{
  return(Teuchos::rcp(pdsMgr_->getParallelMap(Parallel::SOLUTION),false));
}

//-----------------------------------------------------------------------------
// Function      : Builder::vnodeGIDVec()
// Purpose       :
// Special Notes : This is overridden for blockAnalysis types (like MPDE & HB)
// Scope         : Public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 09/11/2019
//-----------------------------------------------------------------------------
const std::vector<int> & Builder::vnodeGIDVec() const
{
  return(lasQueryUtil_->vnodeGIDVec());
}

} // namespace Linear
} // namespace Xyce

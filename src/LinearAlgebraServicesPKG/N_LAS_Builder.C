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

  Graph * overlapGraph = Xyce::Linear::createGraph( *solnOvMap, arrayNZs );

  for( int i = 0; i < numLocalRows_Overlap; ++i )
  {
    std::vector<int>& rcData_i = const_cast<std::vector<int> &>(rcData[i]);
    if( solnOvGMap->localToGlobalIndex(i) != -1 && arrayNZs[i] )
    {
      if( rcData_i[0] == -1 )
        overlapGraph->insertGlobalIndices( solnOvMap->localToGlobalIndex(i), arrayNZs[i]-1, &rcData_i[1] );
      else
        overlapGraph->insertGlobalIndices( solnOvMap->localToGlobalIndex(i), arrayNZs[i], &rcData_i[0] );
    }
  }
  overlapGraph->fillComplete();
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
// Function      : Builder::setupSeparatedLSObjects()
// Purpose       :
// Special Notes : This is for the separation of the solution variables related
//                 to linear and nonlinear devices.
// Scope         : Public
// Creator       : Heidi Thornquist and Ting Mei, Electrical Simulation and Modeling
// Creation Date : 04/01/15
//-----------------------------------------------------------------------------
bool Builder::setupSeparatedLSObjects()
{
  // First create separate linear and nonlinear maps.

  // Get all of the local GIDs
  const std::vector<int>& arrayGIDs = lasQueryUtil_->rowList_GID();

  // Get only the local nonlin GIDs
  // NOTE:  These are not necessarily local in terms of the current solution map.
  const std::vector<int>& nonlinGIDs = lasQueryUtil_->nonlinGIDVec();

  // Collect the global nonlin GIDs
  int tnumNonlinGIDs = nonlinGIDs.size();
  int numNonlinGIDs = nonlinGIDs.size();
  pdsMgr_->getPDSComm()->sumAll( &numNonlinGIDs, &tnumNonlinGIDs, 1 );
  std::vector<int> globalNLGIDs( tnumNonlinGIDs );
  std::vector<int> finalNLGIDs( tnumNonlinGIDs );

  std::vector<int>::iterator nonlin_it;

#ifdef Xyce_PARALLEL_MPI
  int procCnt = pdsMgr_->getPDSComm()->numProc();
  int procID = pdsMgr_->getPDSComm()->procID();
  std::vector<int> numNLGIDs( procCnt, 0.0 );
  Parallel::AllGather(pdsMgr_->getPDSComm()->comm(), numNonlinGIDs, numNLGIDs);
  for (int i=1; i<procCnt; i++)
    numNLGIDs[i] += numNLGIDs[i-1];
  numNLGIDs.insert(numNLGIDs.begin(), 0);

  std::vector<int> tmp_allNLGIDs( tnumNonlinGIDs ); 
  for (int i=numNLGIDs[procID]; i<numNLGIDs[procID+1]; i++)
    tmp_allNLGIDs[ i ] = nonlinGIDs[ i - numNLGIDs[procID] ];

  pdsMgr_->getPDSComm()->sumAll( &tmp_allNLGIDs[0], &globalNLGIDs[0], tnumNonlinGIDs );
  std::sort( globalNLGIDs.begin(), globalNLGIDs.end() );
  std::vector<int>::iterator last = std::unique(globalNLGIDs.begin(), globalNLGIDs.end());
  globalNLGIDs.erase(last, globalNLGIDs.end());
  tnumNonlinGIDs = globalNLGIDs.size();

  // Get list of global nonlinear GIDs that are on this processor.
  // NOTE:  nonlinGIDs may have off-processor GIDs, so this is necessary in parallel.
  nonlin_it = std::set_intersection (globalNLGIDs.begin(), globalNLGIDs.end(),
                                     arrayGIDs.begin(), arrayGIDs.end(),
                                     finalNLGIDs.begin());
  finalNLGIDs.resize(nonlin_it-finalNLGIDs.begin());      
#else
  globalNLGIDs = nonlinGIDs;
  finalNLGIDs = nonlinGIDs;
#endif
  
  // Make linear GID vector.
  // NOTE: This is the initial cut at the linear GIDs, it will be pruned after
  // the graph is generated in the event that any linearGID is associated with
  // a nonlinear node 
  // (ex. a vsrc connected to a nonlinear device, branch equation uses nonlinear node)
  std::vector<int> linGIDs;

  nonlin_it = globalNLGIDs.begin();
  std::vector<int>::const_iterator lcl_it = arrayGIDs.begin();
  for ( ; lcl_it != arrayGIDs.end(); lcl_it++ )
  {
    std::vector<int>::iterator found_it = std::find( nonlin_it, globalNLGIDs.end(), *lcl_it );
    if (found_it == globalNLGIDs.end())
      linGIDs.push_back( *lcl_it );
    else 
      nonlin_it = found_it;
  }
  std::sort( linGIDs.begin(), linGIDs.end() );

  // First create the graphs, then correct the maps, remove linear GIDs if necessary.
  // The separation of the graph into linear and nonlinear Jacobian entries will
  // generate 4 graphs:  linear, nonlinear, linear->nonlinear, nonlinear->linear.
  //
  // Visually this is equivalent to separating a matrix into a Schur-complement form.
  // [ A  B ]
  // [ B' D ]
  // Where A is the linear graph, D is the nonlinear graph, and B/C is the connection
  // between linear and nonlinear portions of the graph.  The entries of the Jacobian
  // for A and B are static, where C and D can change over the course of a transient.
  //
  Parallel::ParMap * solnMap = pdsMgr_->getParallelMap(Parallel::SOLUTION);
  const Graph* graph = pdsMgr_->getMatrixGraph(Parallel::JACOBIAN);
  int numLocalRows = graph->numLocalEntities();

  std::vector<int> linArrayNZs(numLocalRows), nonlinArrayNZs(numLocalRows);
  std::vector<int> linNLArrayNZs(numLocalRows), nlLinArrayNZs(numLocalRows);
  std::vector< std::vector<int> > linArrayCols(numLocalRows), nonlinArrayCols(numLocalRows);
  std::vector< std::vector<int> > linNLArrayCols(numLocalRows), nlLinArrayCols(numLocalRows);
  int maxLinArrayNZs = 0, maxNonlinArrayNZs = 0, maxLinNLArrayNZs = 0, maxNlLinArrayNZs = 0;
 
  int numIndices, maxIndices = graph->maxNumIndices();
  std::vector<int> Indices(maxIndices); 

  // Check if there are any linear GID entries in linArrayCols.  If not, this linear variable
  // must have an equation that uses purely nonlinear nodes (ex. branch equations)
  int numFGIDs = finalNLGIDs.size();
  std::vector<int>::iterator lin_it = linGIDs.begin();
  for ( ; lin_it != linGIDs.end(); lin_it++ )
  {
    graph->extractGlobalRowCopy( *lin_it, maxIndices, numIndices, &Indices[0] );
    std::sort( Indices.begin(), Indices.begin()+numIndices );

    // There are no linear GIDs associated with this linear row, move it to a nonlinear GID.
    bool isEmpty = std::includes( globalNLGIDs.begin(), globalNLGIDs.end(),
                                  Indices.begin(), Indices.begin()+numIndices );
    if (isEmpty)
    { 
      int removedGID = *(lin_it);
      lin_it = linGIDs.erase( lin_it );
      globalNLGIDs.push_back( removedGID );
      std::sort( globalNLGIDs.begin(), globalNLGIDs.end() );
      finalNLGIDs.push_back( removedGID );

      // This could have been the last linear GID in the list, so erase will return lidGIDs.end().
      if (lin_it == linGIDs.end())
      {
        break;
      }
    }
  } 
  if ( numFGIDs != finalNLGIDs.size() )
  {
    std::sort( finalNLGIDs.begin(), finalNLGIDs.end() );
  }

  for (int i=0; i<numLocalRows; i++)
  {
    std::vector<int>::iterator found_it = std::find( linGIDs.begin(), linGIDs.end(), arrayGIDs[i] );

    int baseRow = solnMap->localToGlobalIndex(i);
    graph->extractGlobalRowCopy( baseRow, maxIndices, numIndices, &Indices[0] );
    std::sort( Indices.begin(), Indices.begin()+numIndices );

    // Check if this is a row associated with a linear device
    if (found_it != linGIDs.end())
    {
      lin_it = linGIDs.begin();
      std::vector<int>::iterator lin_end = linGIDs.end();
      std::vector<int>::iterator all_it = Indices.begin();
      for ( ; all_it != Indices.begin()+numIndices; all_it++ )
      {
        std::vector<int>::iterator check_it = std::find( lin_it, linGIDs.end(), *all_it );
        if ( check_it != lin_end )
        {
          linArrayCols[i].push_back( *all_it );
          lin_it = check_it;
        }
        else
        {
          linNLArrayCols[i].push_back( *all_it );
        }
      }
    }
    else
    {
      // Pick out entries that are nonlinear GIDs
      nonlin_it = globalNLGIDs.begin();
      std::vector<int>::iterator all_it = Indices.begin();
      for ( ; all_it != Indices.begin()+numIndices; all_it++ )
      {
        // Check if in nonlinear list
        std::vector<int>::iterator check_it = std::find( nonlin_it, globalNLGIDs.end(), *all_it );
        if (check_it != globalNLGIDs.end()) 
        {
          nonlinArrayCols[i].push_back( *all_it );
          nonlin_it = check_it;
        }
        else
        {
          nlLinArrayCols[i].push_back( *all_it );
        }
      }
    }
    
    // Get the number of nonzeros for these rows.
    linArrayNZs[i] = linArrayCols[i].size();
    nonlinArrayNZs[i] = nonlinArrayCols[i].size();
    linNLArrayNZs[i] = linNLArrayCols[i].size();
    nlLinArrayNZs[i] = nlLinArrayCols[i].size();

    // Record the largest nonzero values for each array.
    if (linArrayNZs[i] > maxLinArrayNZs)
      maxLinArrayNZs = linArrayNZs[i];
    if (nonlinArrayNZs[i] > maxNonlinArrayNZs)
      maxNonlinArrayNZs = nonlinArrayNZs[i];
    if (linNLArrayNZs[i] > maxLinNLArrayNZs)
      maxLinNLArrayNZs = linNLArrayNZs[i];
    if (nlLinArrayNZs[i] > maxNlLinArrayNZs)
      maxNlLinArrayNZs = nlLinArrayNZs[i];
  }

  // Get final nonlinear GID list size, just in case some were added.
  numNonlinGIDs = finalNLGIDs.size();

  // Construct maps from linear and nonlinear data.
  int tnumLinGIDs, numLinGIDs = linGIDs.size();
  pdsMgr_->getPDSComm()->sumAll( &numLinGIDs, &tnumLinGIDs, 1 );
  pdsMgr_->getPDSComm()->sumAll( &numNonlinGIDs, &tnumNonlinGIDs, 1 );

  Parallel::ParMap * linSolnMap = Parallel::createPDSParMap(tnumLinGIDs, numLinGIDs, linGIDs, 0, *(pdsMgr_->getPDSComm()));
  Parallel::ParMap * nonlinSolnMap = Parallel::createPDSParMap(tnumNonlinGIDs, numNonlinGIDs, finalNLGIDs, 0, *(pdsMgr_->getPDSComm()));
  pdsMgr_->addParallelMap(Parallel::LINEAR_SOLUTION, linSolnMap);
  pdsMgr_->addParallelMap(Parallel::NONLINEAR_SOLUTION, nonlinSolnMap);

  if (DEBUG_LINEAR)
  {
    Xyce::dout() << "Linear Solution Map:" << std::endl;
    linSolnMap->print(std::cout);

    Xyce::dout() << "Nonlinear Solution Map:" << std::endl;
    nonlinSolnMap->print(std::cout);

    // Construct graphs from linear and nonlinear data.
    Xyce::dout() << "Base solution map:" << std::endl;
    solnMap->print(std::cout);
  }

  Graph* linearGraph = Xyce::Linear::createGraph( *solnMap, linArrayNZs );
  Graph* nonlinGraph = Xyce::Linear::createGraph( *solnMap, nonlinArrayNZs );
  Graph* linearNLGraph = Xyce::Linear::createGraph( *solnMap, linNLArrayNZs );
  Graph* nlLinearGraph = Xyce::Linear::createGraph( *solnMap, nlLinArrayNZs );

  Graph* lcl_nonlinGraph = Xyce::Linear::createGraph( *nonlinSolnMap, maxNonlinArrayNZs );
  Graph* lcl_nlLinearGraph = Xyce::Linear::createGraph( *nonlinSolnMap, maxNlLinArrayNZs );
  Graph* lcl_linearGraph = Xyce::Linear::createGraph( *linSolnMap, maxLinArrayNZs );
  Graph* lcl_linearNLGraph = Xyce::Linear::createGraph( *linSolnMap, maxLinNLArrayNZs );

  for( int i = 0; i < numLocalRows; ++i )
  {
    // Add linear variable entries to graph.
    if( linArrayNZs[i] )
    {
      linearGraph->insertGlobalIndices( arrayGIDs[i], linArrayNZs[i], &(linArrayCols[i])[0] );
      lcl_linearGraph->insertGlobalIndices( arrayGIDs[i], linArrayNZs[i], &(linArrayCols[i])[0] );
    }
    // Add nonlinear variable entries to graph.
    if( nonlinArrayNZs[i] )
    {
      nonlinGraph->insertGlobalIndices( arrayGIDs[i], nonlinArrayNZs[i], &(nonlinArrayCols[i])[0] );
      lcl_nonlinGraph->insertGlobalIndices( arrayGIDs[i], nonlinArrayNZs[i], &(nonlinArrayCols[i])[0] );
    }

    // Add linear variable nonlinear column entries to graph.
    if( linNLArrayNZs[i] )
    {
      linearNLGraph->insertGlobalIndices( arrayGIDs[i], linNLArrayNZs[i], &(linNLArrayCols[i])[0] );
      lcl_linearNLGraph->insertGlobalIndices( arrayGIDs[i], linNLArrayNZs[i], &(linNLArrayCols[i])[0] );
    }
    // Add nonlinear variable entries to graph.
    if( nlLinArrayNZs[i] )
    {
      nlLinearGraph->insertGlobalIndices( arrayGIDs[i], nlLinArrayNZs[i], &(nlLinArrayCols[i])[0] );
      lcl_nlLinearGraph->insertGlobalIndices( arrayGIDs[i], nlLinArrayNZs[i], &(nlLinArrayCols[i])[0] );
    }
  }

  linearGraph->fillComplete();
  nonlinGraph->fillComplete();
  linearNLGraph->fillComplete();
  nlLinearGraph->fillComplete();

  pdsMgr_->addMatrixGraph( Parallel::GLOBAL_LINEAR_JACOBIAN, linearGraph );
  pdsMgr_->addMatrixGraph( Parallel::GLOBAL_LIN_NONLIN_JACOBIAN, linearNLGraph );
  pdsMgr_->addMatrixGraph( Parallel::GLOBAL_NONLINEAR_JACOBIAN, nonlinGraph );
  pdsMgr_->addMatrixGraph( Parallel::GLOBAL_NONLIN_LIN_JACOBIAN, nlLinearGraph );

  lcl_linearGraph->fillComplete();
  lcl_nonlinGraph->fillComplete();
  lcl_linearNLGraph->fillComplete( *nonlinSolnMap, *linSolnMap );
  lcl_nlLinearGraph->fillComplete( *linSolnMap, *nonlinSolnMap );

  pdsMgr_->addMatrixGraph( Parallel::LINEAR_JACOBIAN, lcl_linearGraph );
  pdsMgr_->addMatrixGraph( Parallel::LIN_NONLIN_JACOBIAN, lcl_linearNLGraph );
  pdsMgr_->addMatrixGraph( Parallel::NONLINEAR_JACOBIAN, lcl_nonlinGraph );
  pdsMgr_->addMatrixGraph( Parallel::NONLIN_LIN_JACOBIAN, lcl_nlLinearGraph );

  if (DEBUG_LINEAR)
  {
    std::cout << "Linear graph: " << std::endl;
    linearGraph->print(std::cout);  
    std::cout << "Local Linear graph: " << std::endl;
    lcl_linearGraph->print(std::cout);  
    std::cout << "Nonlinear graph: " << std::endl;
    nonlinGraph->print(std::cout); 
    std::cout << "Local Nonlinear graph: " << std::endl;
    lcl_nonlinGraph->print(std::cout); 
    std::cout << "Linear->Nonlinear graph: " << std::endl;
    linearNLGraph->print(std::cout);  
    std::cout << "Local Linear->Nonlinear graph: " << std::endl;
    lcl_linearNLGraph->print(std::cout);  
    std::cout << "Nonlinear->Linear graph: " << std::endl;
    nlLinearGraph->print(std::cout); 
    std::cout << "Local Nonlinear->Linear graph: " << std::endl;
    lcl_nlLinearGraph->print(std::cout); 
    const Graph * baseGraph = pdsMgr_->getMatrixGraph( Parallel::JACOBIAN );
    std::cout << "Base graph: " << std::endl;
    baseGraph->print(std::cout);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Builder::getSeparatedSolnMap
// Purpose       :
// Special Notes : This is for the separation of the solution variables related
//                 to linear and nonlinear devices.
// Scope         : Public
// Creator       : Heidi Thornquist and Ting Mei, Electrical Simulation and Modeling
// Creation Date : 04/01/15
//-----------------------------------------------------------------------------
void Builder::getSeparatedSolnMap( RCP<Parallel::ParMap>& linear_map,
                                   RCP<Parallel::ParMap>& nonlin_map
                                 ) const
{
  linear_map = Teuchos::rcp(pdsMgr_->getParallelMap(Parallel::LINEAR_SOLUTION),false);
  nonlin_map = Teuchos::rcp(pdsMgr_->getParallelMap(Parallel::NONLINEAR_SOLUTION),false);
}

//-----------------------------------------------------------------------------
// Function      : Builder::getSeparatedGraph
// Purpose       :
// Special Notes : This is for the separation of the solution variables related
//                 to linear and nonlinear devices, using linear and nonlinear map.
// Scope         : Public
// Creator       : Heidi Thornquist and Ting Mei, Electrical Simulation and Modeling
// Creation Date : 04/01/15
//-----------------------------------------------------------------------------
void Builder::getSeparatedGraph( RCP<const Graph>& linear_graph,
                                 RCP<const Graph>& linNonlin_graph,
                                 RCP<const Graph>& nonlin_graph,
                                 RCP<const Graph>& nonlinLin_graph
                               ) const
{
  linear_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::LINEAR_JACOBIAN ),false);
  linNonlin_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::LIN_NONLIN_JACOBIAN ),false);
  nonlin_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::NONLINEAR_JACOBIAN ),false);
  nonlinLin_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::NONLIN_LIN_JACOBIAN ),false);
}

//-----------------------------------------------------------------------------
// Function      : Builder::getGlobalSeparatedGraph
// Purpose       :
// Special Notes : This is for the separation of the solution variables related
//                 to linear and nonlinear devices, using global solution map.
// Scope         : Public
// Creator       : Heidi Thornquist and Ting Mei, Electrical Simulation and Modeling
// Creation Date : 04/01/15
//-----------------------------------------------------------------------------
void Builder::getGlobalSeparatedGraph( RCP<const Graph>& linear_graph,
                                       RCP<const Graph>& linNonlin_graph,
                                       RCP<const Graph>& nonlin_graph,
                                       RCP<const Graph>& nonlinLin_graph
                                     ) const
{
  linear_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::GLOBAL_LINEAR_JACOBIAN ),false);
  linNonlin_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::GLOBAL_LIN_NONLIN_JACOBIAN ),false);
  nonlin_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::GLOBAL_NONLINEAR_JACOBIAN ),false);
  nonlinLin_graph = Teuchos::rcp(pdsMgr_->getMatrixGraph( Parallel::GLOBAL_NONLIN_LIN_JACOBIAN ),false);
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

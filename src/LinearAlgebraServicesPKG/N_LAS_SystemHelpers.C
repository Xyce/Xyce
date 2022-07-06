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
// Purpose        : This is collection of non-member functions that help
//                  in the construction of linear systems.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 04/06/15
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <utility>
#include <numeric>

#include <N_LAS_SystemHelpers.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_BlockMatrix.h>

#include <N_LAS_Graph.h>

#include <N_LAS_Problem.h>
#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_LAS_QueryUtil.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Math.h>

namespace Xyce {
namespace Linear {

// Set vector values (used for .IC)
void setInitialConditions( const System& system, Vector& vector,
                           const std::map<int, double>& op )
{ 
  BlockVector* b_vector = dynamic_cast< BlockVector* >( &vector );
  if (b_vector)
  {
    std::map<int, double> new_op( op );
    system.builder().createInitialConditionOp( new_op );

    for (std::map<int, double>::const_iterator it = new_op.begin();
         it != new_op.end(); ++it)
    {
      vector[(*it).first] = (*it).second;
    }
  }
  else
  {
    for (std::map<int, double>::const_iterator it = op.begin();
         it != op.end(); ++it)
    {
      vector[(*it).first] = (*it).second;
    }
  }

  // Import overlaps for shared nodes
  vector.importOverlap();
}

void setInitialConditions( const System& system, Vector& vector,
                           const NodeNameMap& op, double value )
{
  BlockVector* b_vector = dynamic_cast< BlockVector* >( &vector );
  if (b_vector)
  {
    std::vector<int> new_op;
    new_op.reserve( vector.globalLength() );

    for (NodeNameMap::const_iterator it = op.begin(); it != op.end(); ++it)
    {
      new_op.push_back( (*it).second );
    }
    system.builder().createInitialConditionOp( new_op );

    for (std::vector<int>::iterator it = new_op.begin(); it != new_op.end(); ++it)
    {
      vector[*it] = value;
    }
  }
  else
  {
    for (NodeNameMap::const_iterator it = op.begin();
         it != op.end(); ++it)
    {
      vector[(*it).second] = value;
    }
  }

  // Import overlaps for shared nodes
  vector.importOverlap();
}

//-----------------------------------------------------------------------------
// Function      : extractValues
// Purpose       : A helper function that extracts the values into another matrix
// Special Notes : When a matrix is being separated into various components,
//               : we need to have a way to extract out the subset of the entries
//               : into different matrices.  It is assumed that the matrix has
//               : the same parallel distribution as the original, just a subset
//               : of the entries.
// Creator       : Heidi Thornquist, SNL, Electrical Systems Modeling
// Creation Date : 04/06/15
//-----------------------------------------------------------------------------
void extractValues( const Matrix& inputMatrix,
                    std::vector<Teuchos::RCP<Matrix> >& outputMatrices )
{
  int numMatrices = outputMatrices.size();

  if (numMatrices > 0)
  {
    // Get basic information from input graph.
    const Graph* inGraph = inputMatrix.getGraph();   
    int in_currNNZ, in_maxNNZ = inGraph->maxNumIndices();
    std::vector<int> inIdxs(in_maxNNZ);
    std::vector<double> inValues(in_maxNNZ);

    // Extract the values from the inputMatrix and place them in each of the 
    // outputMatrices.  This assumes that the outputMatrices have a subset of
    // the graph used to create the inputMatrix.
    for (int i=0; i<numMatrices; i++)
    {
      const Graph* outGraph = outputMatrices[i]->getGraph();   
     
      // Get basic information from this graph.
      int numMyRows = outGraph->numLocalEntities();
      int currNNZ, maxNNZ = outGraph->maxNumIndices();
      std::vector<int> outIdxs(maxNNZ);
      std::vector<double> values(maxNNZ);
      for (int j=0; j<numMyRows; j++)
      {
         int globalRow = outGraph->localToGlobalRowIndex(j);

         // Extract entries needed by the output matrix.
         outGraph->extractGlobalRowCopy( globalRow, maxNNZ, currNNZ, &outIdxs[0] );
  
         if (currNNZ > 0)
         { 
           // Extract entries in the input matrix.
           inputMatrix.getRowCopy( globalRow, in_maxNNZ, in_currNNZ, &inValues[0], &inIdxs[0] );
      
           int inIdxPtr = 0;
           for (int k=0; k<currNNZ; k++)
           {
             int searchIdx = outIdxs[k];
             bool foundIdx = false;
             for (int l=inIdxPtr; l<in_currNNZ; l++)
             {
               // If the requested index is less than the current index, it will not be found.
               if (inIdxs[l] > searchIdx)
                 break;

               // This search presumes the list of entries are sorted.
               if (searchIdx == inIdxs[l])
               {
                 values[k] = inValues[l];
                 inIdxPtr = l+1;
                 foundIdx = true;
                 break;
               }
             }
             if (!foundIdx)
               values[k] = 0.0;
           }
           // Now insert the requested values into the outputMatrices[i].
           outputMatrices[i]->putRow( globalRow, currNNZ, &values[0], &outIdxs[0] ); 
         }            
      } 
    }
  }
}

bool checkProblemForNaNs( const Linear::Problem& problem, 
                          std::vector< std::pair<int, int> >& nanEntries )
{
  Linear::Problem& prob = const_cast<Linear::Problem &>( problem );

  int numrhs = prob.getLHS()->numVectors();
  bool foundNaN = false;
  nanEntries.clear();
  std::vector< std::pair<int, int> > tmp_nanEntries;
  std::vector<double> resNorm(numrhs,0.0);

  // Check the solution vector for NaNs (using 2-norm calculation)
  prob.getLHS()->lpNorm( 2, &resNorm[0] );
  for (int i=0; i<numrhs; ++i)
  { 
    if (std::isnan(resNorm[i]) || std::isinf(resNorm[i]))
      foundNaN = true;
  }

  // If a NaN has been found, it's likely in the RHS or Matrix.
  // Otherwise, the linear solver introduced a NaN.
  if (foundNaN)
  {
    bool foundRHSNaN = false;
    int numrows = prob.getLHS()->localLength();
    prob.getRHS()->lpNorm( 2, &resNorm[0] );

    for (int i=0; i<numrhs; ++i)
    {
      if (std::isnan(resNorm[i]) || std::isinf(resNorm[i]))
        foundRHSNaN = true;
    }

    // If foundRHSNaN, look for which entries have NaNs in the RHS vector.
    std::vector< std::pair< int, int > > entriesNaN;
    if (foundRHSNaN)
    {
      const Parallel::ParMap& map = *(prob.getRHS()->pmap());

      for (int j=0; j<numrhs; ++j)
      {
        const Vector& mvCol_j = *(prob.getRHS()->getVectorView(j));
        for (int i=0; i<numrows; ++i)
        {
          const double & val = mvCol_j[ i ];
          if (std::isnan(val) || std::isinf(val))
            tmp_nanEntries.push_back( std::make_pair( map.localToGlobalIndex(i), -1 ) );
        }
      }
    }
    else // look in the Jacobian matrix entries for NaNs
    {
      if (!prob.matrixFree())
      {
        Matrix & A = *prob.getMatrix();
        const Graph & Agraph = *A.getGraph();
        int* indices;
        double* values;
        int numIndices;

        for( int i=0; i<numrows; ++i )
        {
          A.getLocalRowView( i, numIndices, values, indices );
          for ( int j=0; j<numIndices; ++j)
          {
            if (std::isnan(values[j]) || std::isinf(values[j]))
              tmp_nanEntries.push_back( std::make_pair( Agraph.localToGlobalRowIndex(i), 
                                                    Agraph.localToGlobalColIndex(indices[j]) ) );
          }
        }
      }
    }
  }

  // Now synchronize the data, if executing in parallel.
  const Parallel::Communicator& comm = *(prob.getRHS()->pdsComm());
  if (!comm.isSerial())
  {
    int numProc = comm.numProc();
    std::vector<int> totProcNaNs( numProc, 0 ), procNaNs( numProc, 0 );
    procNaNs[comm.procID()] = tmp_nanEntries.size();
    comm.sumAll( &procNaNs[0], &totProcNaNs[0], numProc );
    int totNumNaNs = std::accumulate( totProcNaNs.begin(), totProcNaNs.end(), 0 );

    if (totNumNaNs)
    {
      foundNaN = true;
      nanEntries.resize( totNumNaNs );
      std::vector<int> totNaNGIDs( 2*totNumNaNs, 0 ), NaNGIDs( 2*totNumNaNs, 0 );

      // Compute the starting pointer where the NaN pairs should be inserted
      int startPtr = 0;
      for ( int i=0; i<comm.procID(); ++i )
        startPtr += totProcNaNs[ i ];

      for ( int i=0; i<tmp_nanEntries.size(); ++i )
      {
        NaNGIDs[ 2*(startPtr+i) ] = tmp_nanEntries[i].first;
        NaNGIDs[ 2*(startPtr+i)+1 ] = tmp_nanEntries[i].second;
      }
     
      // Communicate the NaN GIDs through a global sum 
      comm.sumAll( &totNaNGIDs[0], &NaNGIDs[0], 2*totNumNaNs );

      // Reconstruct the global list of NaNs from totNaNGIDs
      for ( int i=0; i<totNumNaNs; ++i )
        nanEntries[i] = std::make_pair( totNaNGIDs[ 2*i ], totNaNGIDs[ 2*i+1 ] );
    }
  }
  else
  {
    nanEntries = tmp_nanEntries;
  } 

  return foundNaN; 
} 

bool checkVectorForNaNs( const Linear::MultiVector& vector,
                         std::vector< int >& nanEntries )
{
  int numrhs = vector.numVectors();
  bool foundNaN = false;
  nanEntries.clear();
  std::vector<int> tmp_nanEntries;
  std::vector<double> resNorm(numrhs,0.0);

  // Check the vector for NaNs (using 2-norm calculation)
  vector.lpNorm( 2, &resNorm[0] );
  for (int i=0; i<numrhs; ++i)
  {
    if (std::isnan(resNorm[i]) || std::isinf(resNorm[i]))
      foundNaN = true;
  }
  // If a NaN has been found, it's likely in the RHS or Matrix.
  // Otherwise, the linear solver introduced a NaN.
  if (foundNaN)
  {
    std::vector<int> entriesNaN;
    int numrows = vector.localLength();
    const Parallel::ParMap& map = *(vector.pmap());

    for (int j=0; j<numrhs; ++j)
    {
      const Vector& mvCol_j = *(vector.getVectorView(j));
      for (int i=0; i<numrows; ++i)
      {
        const double & val = mvCol_j[ i ];
        if (std::isnan(val) || std::isinf(val))
          tmp_nanEntries.push_back( map.localToGlobalIndex(i) );
      }
    }
    // Now synchronize the data, if executing in parallel.
    const Parallel::Communicator& comm = *(vector.pdsComm());
    if (!comm.isSerial())
    {
      int numProc = comm.numProc();
      std::vector<int> totProcNaNs( numProc, 0 ), procNaNs( numProc, 0 );
      procNaNs[comm.procID()] = tmp_nanEntries.size();
      comm.sumAll( &procNaNs[0], &totProcNaNs[0], numProc );
      int totNumNaNs = std::accumulate( totProcNaNs.begin(), totProcNaNs.end(), 0 );

      if (totNumNaNs)
      {
        foundNaN = true;
        nanEntries.resize( totNumNaNs );
        std::vector<int> totNaNGIDs( totNumNaNs, 0 );

        // Compute the starting pointer where the NaN pairs should be inserted
        int startPtr = 0;
        for ( int i=0; i<comm.procID(); ++i )
          startPtr += totProcNaNs[ i ];

        for ( int i=0; i<tmp_nanEntries.size(); ++i )
          nanEntries[ startPtr+i ] = tmp_nanEntries[i];
    
        // Communicate the NaN GIDs through a global sum
        comm.sumAll( &totNaNGIDs[0], &nanEntries[0], totNumNaNs );
      }
    }
    else
    {
      nanEntries = tmp_nanEntries;
    }
  }

  return foundNaN;
}
                   
} // namespace Linear
} // namespace Xyce

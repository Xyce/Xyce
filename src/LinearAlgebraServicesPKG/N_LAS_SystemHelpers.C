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

#include <N_LAS_SystemHelpers.h>

#include <N_PDS_ParMap.h>
#include <N_PDS_Comm.h>

#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_BlockVector.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_BlockMatrix.h>

#include <N_LAS_Graph.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_LAS_QueryUtil.h>

#include <N_ERH_ErrorMgr.h>

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
    std::cout << "Initial conditions being set on a block vector!" << std::endl;
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
                    
} // namespace Linear
} // namespace Xyce

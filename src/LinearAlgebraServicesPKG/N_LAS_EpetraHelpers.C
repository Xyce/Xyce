//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
//                  in the construction of block linear systems, like those
//                  found in AC analysis.
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL, Electrical Systems Modeling
//
// Creation Date  : 06/22/11
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>
#include <sstream>

#include <N_LAS_SystemHelpers.h>
#include <N_LAS_EpetraHelpers.h>
#include <N_LAS_EpetraMultiVector.h>
#include <N_LAS_Operator.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_EpetraProblem.h>
#include <N_LAS_EpetraImporter.h>
#include <N_LAS_EpetraGraph.h>
#include <N_LAS_EpetraBlockMultiVector.h>
#include <N_LAS_EpetraVector.h>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>

#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_BlockMapOut.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Linear {

  // Non-member creation methods.
MultiVector* createMultiVector( const Parallel::ParMap & map,
                                int numVectors )
{
  return new EpetraMultiVector( map, numVectors );
}

MultiVector* createMultiVector( const Parallel::ParMap & map,
                                const Parallel::ParMap & ol_map,
                                int numVectors )
{
  return new EpetraMultiVector( map, ol_map, numVectors );
}

Vector* createVector( const Parallel::ParMap & map ) 
{
  return new EpetraVector( map );
}

Vector* createVector( const Parallel::ParMap & map, const Parallel::ParMap & ol_map )
{
  return new EpetraVector( map, ol_map );
}

Matrix* createMatrix( const Graph* overlapGraph,
                      const Graph* baseGraph )
{
  return new EpetraMatrix( overlapGraph, baseGraph );
}

Graph* createGraph( const Parallel::ParMap & solution_overlap,
                    const Parallel::ParMap & solution_overlap_ground,
                    const std::vector<int>& numIndicesPerRow,
                    const std::vector<std::vector<int> >& rcData)
{
  return new EpetraGraph( solution_overlap, solution_overlap_ground, numIndicesPerRow, rcData );
}

Problem* createProblem( Matrix* A, MultiVector* x, MultiVector* b )
{
  return new EpetraProblem( A, x, b );
}
               
Problem* createProblem( Operator* Op, MultiVector* x, MultiVector* b )
{
  return new EpetraProblem( Op, x, b );
}

Importer* createImporter( const Parallel::ParMap & target_map, 
                          const Parallel::ParMap & source_map )
{
  return new EpetraImporter( target_map, source_map );
}

void writeToFile(const Epetra_LinearProblem& problem, std::string prefix, 
                 int file_number, bool write_map)
{
  if (write_map) {
    std::string file_name = prefix + "_BlockMap.mm";
    EpetraExt::BlockMapToMatrixMarketFile( file_name.c_str(), (problem.GetMatrix())->Map() );
  }

  std::stringstream char_file_name(""), file_name("");
  char_file_name << prefix << "_Matrix" << file_number << ".mm";
  std::string sandiaReq = "Sandia National Laboratories is a multimission laboratory managed and operated by National Technology and\n%";
  sandiaReq += " Engineering Solutions of Sandia LLC, a wholly owned subsidiary of Honeywell International Inc. for the\n%";
  sandiaReq += " U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.\n%\n% Xyce circuit matrix.\n%%";
  EpetraExt::RowMatrixToMatrixMarketFile( char_file_name.str().c_str(), *(problem.GetMatrix()), sandiaReq.c_str() );
  file_name << prefix << "_RHS" << file_number << ".mm";
  EpetraExt::MultiVectorToMatrixMarketFile( file_name.str().c_str(), *(problem.GetRHS()) );
}

void writeToFile( const Epetra_MultiVector& vector, const char * filename, 
                  bool useLIDs, bool mmFormat )
{
  int numProcs = vector.Comm().NumProc();
  int localRank = vector.Comm().MyPID();
  int masterRank = 0;

  if (!mmFormat)
  {
    for( int p = 0; p < numProcs; ++p )
    {
      //A barrier inside the loop so each processor waits its turn.
      vector.Comm().Barrier();

      if(p == localRank)
      {
        FILE *file = NULL;

        if(masterRank == localRank)
        {
          //This is the master processor, open a new file.
          file = fopen(filename,"w");

          //Write the RDP_MultiVector dimension n into the file.
          fprintf(file,"%d\n",vector.GlobalLength());
        }
        else
        {
          //This is not the master proc, open file for appending
          file = fopen(filename,"a");
        }

        //Now loop over the local portion of the RDP_MultiVector.
        int length  = vector.MyLength();
        int numVecs = vector.NumVectors();

        for (int i = 0; i < numVecs; ++i)
          for (int j = 0; j < length; ++j)
          {
            int loc = vector.Map().GID(j);
            if( useLIDs ) loc = j;
            fprintf(file,"%d %d %20.13e\n",i,loc,vector[i][j]);
          } 
        fclose(file);
      } 
    } 
  } 
  else
  {
    EpetraExt::MultiVectorToMatrixMarketFile( filename, vector );
  } 
}
               
// ///////////////////////////////////////////////////////////////////
//
// Implementation of the Xyce::Linear::EpetraTransOp class.
//
// ///////////////////////////////////////////////////////////////////

EpetraTransOp::EpetraTransOp (const Teuchos::RCP<Epetra_Operator> &Op)
  : Epetra_Op(Op)
{
}

// The version of Apply() that takes two arguments and returns int
// implements the Epetra_Operator interface.
int
EpetraTransOp::Apply (const Epetra_MultiVector &X,
                     Epetra_MultiVector &Y) const
{
  // This operation computes Y = A^{-1}*X.
  const int info = Epetra_Op->Apply( X, Y );

  return info;
}

// This implements Epetra_Operator::ApplyInverse().
int
EpetraTransOp::ApplyInverse (const Epetra_MultiVector &X,
                            Epetra_MultiVector &Y) const
{
  // This operation computes Y = A*X.
  const int info = Epetra_Op->ApplyInverse( X, Y );

  return info;
}


} // namespace Linear
} // namespace Xyce

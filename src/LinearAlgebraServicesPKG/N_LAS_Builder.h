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
// Purpose        : Builder for LAS objects, hides parallel map stuff
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

#ifndef  Xyce_LAS_Builder_h
#define  Xyce_LAS_Builder_h


#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>

#include <N_LAS_fwd.h>

// to eliminate RCP warnings, putting N_PDS_Manager header here.
#include <N_PDS_Manager.h>

class N_PDS_ParMap;

using Teuchos::RCP;
using Teuchos::rcp;

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : Builder
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/22/02
//-----------------------------------------------------------------------------
class Builder
{

public:

  // Default Constructor
  Builder()
  : pdsMgr_(0),
    lasQueryUtil_(0)
  {}

  // Destructor
  virtual ~Builder() {}

  // Registration methods for necessary utilities
  bool registerPDSManager(N_PDS_Manager * PDS_Manager)
  {
    return (pdsMgr_ = PDS_Manager);
  }

  bool registerQueryUtil(QueryUtil * LAS_QUtil)
  {
    return (lasQueryUtil_ = LAS_QUtil);
  }

  // Vector and Matrix creators which use QueryUtil and N_PDS_ParMap
  // attributes to transparently construct proper objects for this linear
  // system

  // Multivector factory with num vectors and initial value
  virtual MultiVector * createMultiVector( const int numVectors = 1, const double value = 0.0 ) const;
  // State Multivector factory with num vectors and initial value
  virtual MultiVector * createStateMultiVector( const int numVectors = 1, const double value = 0.0 ) const;
  // Store Multivector factory with num vectors and initial value
  virtual MultiVector * createStoreMultiVector( const int numVectors = 1, const double value = 0.0 ) const;
  // Vector factory with initial value
  virtual Vector * createVector( const double value = 0.0 ) const;
  // State-vector factory
  virtual Vector * createStateVector( const double value = 0.0 ) const;
  // Store-vector factory
  virtual Vector * createStoreVector( const double value = 0.0 ) const;
  // LeadCurrent-vector factory
  virtual Vector * createLeadCurrentVector( const double value = 0.0 ) const;

#if 1
  // ERK.  HACK!!! FIX THIS
  // State-vector factory
  virtual Vector * createStateQuadVector( const double value = 0.0 ) const { return createStateVector(value); }
  // Store-vector factory
  virtual Vector * createStoreQuadVector( const double value = 0.0 ) const { return createStoreVector(value); }
  // LeadCurrent-vector factory
  virtual Vector * createLeadCurrentQuadVector( const double value = 0.0 ) const { return createLeadCurrentVector(value); }
#endif

  // Matrix factory
  virtual Matrix * createMatrix( const double initialValue = 0.0 ) const;

  //Coloring Assoc with Variable Types in Solution Vector
  virtual const std::vector<int> & createSolnColoring() const;

  //Coloring needed for imposing .IC and .NODESET
  virtual const std::vector<int> & createInitialConditionColoring() const;

  virtual bool generateParMaps();

  virtual bool generateGraphs();

  // This is an on-demand capability for the builder to analyze the
  // current solution maps and Jacobian graphs and separate them based
  // on linear vs. nonlinear GID information from topology.
  virtual bool setupSeparatedLSObjects();
  virtual void getSeparatedSolnMap( RCP<N_PDS_ParMap>& linear_map, 
                                    RCP<N_PDS_ParMap>& nonlin_map
                                  ) const;

  // Return the graphs with reference to the separated solution map
  virtual void getSeparatedGraph( RCP<Graph>& linear_graph,
                                  RCP<Graph>& linNonlin_graph,
                                  RCP<Graph>& nonlin_graph,
                                  RCP<Graph>& nonlinLin_graph
                                ) const;

  // Return the graphs with reference to the global solution map
  virtual void getGlobalSeparatedGraph( RCP<Graph>& linear_graph,
                                        RCP<Graph>& linNonlin_graph,
                                        RCP<Graph>& nonlin_graph,
                                        RCP<Graph>& nonlinLin_graph
                                      ) const;

  virtual RCP<const N_PDS_ParMap> getSolutionMap() const;
  
  virtual RCP<N_PDS_ParMap> getSolutionMap();

  virtual const std::vector<int> & vnodeGIDVec() const;

  virtual N_PDS_Comm* getPDSComm() const
  {
    return pdsMgr_->getPDSComm();
  }

protected:
  mutable std::vector<int>      solnColoring_;
  mutable std::vector<int>      icColoring_;

  N_PDS_Manager *       pdsMgr_;
  QueryUtil *           lasQueryUtil_;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_LAS_Builder_h

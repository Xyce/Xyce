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

#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_PDS_Manager.h>

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
  bool registerPDSManager(Parallel::Manager * PDS_Manager)
  {
    return (pdsMgr_ = PDS_Manager);
  }

  bool registerQueryUtil(QueryUtil * LAS_QUtil)
  {
    return (lasQueryUtil_ = LAS_QUtil);
  }

  // Vector and Matrix creators which use QueryUtil and ParMap
  // attributes to transparently construct proper objects for this linear
  // system

  // Multivector factory with num vectors 
  virtual MultiVector * createMultiVector( const int numVectors = 1 ) const;
  // State Multivector factory with num vectors 
  virtual MultiVector * createStateMultiVector( const int numVectors = 1 ) const;
  // Store Multivector factory with num vectors 
  virtual MultiVector * createStoreMultiVector( const int numVectors = 1 ) const;
  // Vector factory 
  virtual Vector * createVector() const;
  // State-vector factory
  virtual Vector * createStateVector() const;
  // Store-vector factory
  virtual Vector * createStoreVector() const;
  // LeadCurrent-vector factory
  virtual Vector * createLeadCurrentVector() const;

  // Matrix factory
  virtual Matrix * createMatrix() const;

  //Coloring Assoc with Variable Types in Solution Vector
  virtual const std::vector<int> & createSolnColoring() const;

  // Convert topology op data to analysis specific op data
  virtual bool createInitialConditionOp( std::map<int,double> & op ) const
  { return false; }

  virtual bool createInitialConditionOp( std::vector<int>& op ) const
  { return false; }

  //Coloring needed for imposing .IC and .NODESET
  virtual const std::vector<int> & createInitialConditionColoring() const;

  virtual bool generateParMaps();

  virtual bool generateGraphs();

  virtual RCP<const Parallel::ParMap> getSolutionMap() const;
  
  virtual RCP<Parallel::ParMap> getSolutionMap();

  virtual const std::vector<int> & vnodeGIDVec() const;

  virtual Parallel::Communicator* getPDSComm() const
  {
    return pdsMgr_->getPDSComm();
  }

protected:
  mutable std::vector<int>      solnColoring_;
  mutable std::vector<int>      icColoring_;

  Parallel::Manager *       pdsMgr_;
  QueryUtil *               lasQueryUtil_;
};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_LAS_Builder_h

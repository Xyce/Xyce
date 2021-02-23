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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/09/00
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef N_LAS_QueryUtil_h
#define N_LAS_QueryUtil_h 1

#include <string>
#include <vector>

#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : QueryUtil
// Purpose       : Abstract class for Query Utility used by LAS to
//                 determine parameters for linear system instantiation
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
class QueryUtil
{
public:
  static void populateMetadata(IO::PkgOptionsMgr &options_manager);

  // Default constructor
  QueryUtil();

  // Destructor
  virtual ~QueryUtil() {}

  // Method to register the package options manager
  bool registerTimeOptions(const Util::OptionBlock & OB) { return true; }
  bool registerOptions(const Util::OptionBlock & OB);
     
  // Setup row/col data for linear solver including reorder.
  virtual bool setupRowCol() = 0;
  //
  // Get access to the supernode flag
  virtual bool supernodeFlag() { return supernode_; }
  //
  // Get access to the names file flag
  virtual bool namesFileFlag() { return namesFile_; }
  
  //------ accessor methods
  virtual int numGlobalRows() const = 0;
  virtual int numLocalRows() const { return numGlobalRows(); }
  virtual int numExternRows() const { return 0; }
  virtual int numGlobalExternRows() const { return 0; }

  virtual int numGlobalStateVars() const = 0;
  virtual int numLocalStateVars() const { return numGlobalStateVars(); }

  virtual int numGlobalStoreVars() const = 0;
  virtual int numLocalStoreVars() const { return numGlobalStoreVars(); }
  
  virtual int numGlobalLeadCurrentVars() const = 0;
  virtual int numLocalLeadCurrentVars() const { return numGlobalLeadCurrentVars(); }

  virtual int numGlobalNZs() const = 0;
  virtual int numLocalNZs() const { return numGlobalNZs(); }

  virtual const std::vector<int> & rowList_GID() const = 0;
  virtual const std::vector< std::pair<int,int> > & rowList_ExternGID() const
  {
    static std::vector< std::pair<int,int> > dummy;
    return dummy;
  }

  virtual const std::vector<int> & rowList_StateGID() const = 0;
  virtual const std::vector<int> & rowList_StoreGID() const = 0;
  virtual const std::vector<int> & rowList_LeadCurrentGID() const = 0;
  
  virtual const std::vector<int> & rowList_NumNZs() const = 0;

  virtual const std::vector< std::vector<int> > & rowList_ColList() const = 0;

  virtual const std::vector<int> & vnodeGIDVec() const = 0;

  virtual const std::vector<int> & vsrcGIDVec() const = 0;

  virtual const std::vector<int> & nonlinGIDVec() const = 0;

  virtual const std::vector<char> & rowList_VarType() const = 0;

  virtual bool isClean() const { return isClean_; }
  virtual void cleanRowLists() = 0;

protected:
  bool checkConnectivity_;
  bool supernode_;
  bool isClean_;
  bool namesFile_;
};

bool registerPkgOptionsMgr(QueryUtil &linear_solve_utility, IO::PkgOptionsMgr &options_manager);

} // namespace Linear
} // namespace Xyce

#endif

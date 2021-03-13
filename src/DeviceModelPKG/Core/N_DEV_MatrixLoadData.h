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
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/30/02
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_MatrixLoadData_h
#define Xyce_N_DEV_MatrixLoadData_h

#include <vector>

namespace Xyce {
namespace Device {

class colData;
class valData;

//-----------------------------------------------------------------------------
// Class         : MatrixLoadData
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/30/02
//-----------------------------------------------------------------------------
class MatrixLoadData
{
public:
  MatrixLoadData ();

  MatrixLoadData (const MatrixLoadData & right);

  ~MatrixLoadData ();

  bool initializeAll (int isizeTmp = 100);

  void resizeTestJacSolData(int size);
  void resizeTestJacQData(int size);
  void resizeTestJacStateData(int size);


  void resizeSolnSizedVectors (int size);
  void resizeStateSizedVectors (int size);

public:
  int isize;
  int isizeNumJac;

  // temporary jacobian load structures:
  std::vector<int>    cols;
  std::vector<double> vals;
  std::vector<double> Qvals;

  // temporary numerical jacobian load structures:
  std::vector<valData> val_local;
  std::vector<valData> Qval_local;
  std::vector<colData> col_local;
  std::vector<int>     row_local;
  std::vector<int>     internalFlag;

  // Structures used by the "testJacobian" function.
  std::vector<std::vector<double> > numJac;
  std::vector< std::vector<double> > saveJac;
  std::vector< std::vector<double> > devJac;
  std::vector< std::vector<double> > diffJac;
  std::vector< std::vector<double> > relJac;

  std::vector<std::vector<double> > numJacQ;
  std::vector< std::vector<double> > saveJacQ;
  std::vector< std::vector<double> > devJacQ;
  std::vector< std::vector<double> > diffJacQ;
  std::vector< std::vector<double> > relJacQ;

  std::vector<std::vector<int> > status;
  std::vector<std::vector<int> > stencil;
  std::vector<std::vector<int> > statusQ;

  std::vector<double> saveF;
  std::vector<double> pertF;
  std::vector<double> origF;
  std::vector<double> saveQ;
  std::vector<double> pertQ;
  std::vector<double> origQ;
  std::vector<double> saveB;
  std::vector<double> pertB;
  std::vector<double> origB;

  std::vector<double> saveSoln;
  std::vector<double> pertSoln;
  std::vector<double> saveCurrSoln;

  std::vector<double> saveLastState;
  std::vector<double> saveCurrState;
  std::vector<double> saveNextState;
  std::vector<double> saveStateDerivs;
};

//-----------------------------------------------------------------------------
// Class         : colData
// Purpose       : This class contains a vector of column indices.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class colData
{
public:
  colData (int isizeTmp = 100):
    isize(isizeTmp), col()
  { col.reserve(isizeTmp); }

public:
  int isize;
  std::vector<int> col;
};

//-----------------------------------------------------------------------------
// Class         : valData
// Purpose       : This class contains a vector of value indices.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/30/02
//-----------------------------------------------------------------------------
class valData
{
private:
protected:
public:
  valData (int isizeTmp = 100):
    isize(isizeTmp), val() { val.reserve(isizeTmp); }

private:
protected:
public:
  int isize;
  std::vector<double> val;
};

} // namespace Device
} // namespace Xyce

#endif


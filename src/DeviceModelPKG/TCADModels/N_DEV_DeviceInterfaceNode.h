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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/08/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceInterfaceNode_h
#define Xyce_N_DEV_DeviceInterfaceNode_h

#include <N_LAS_fwd.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceInterfaceNode
// Purpose       : This class contains information about a single circuit node
//                 which is connected to a single device electrode.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class DeviceInterfaceNode
{
public:
  DeviceInterfaceNode ()
    : eName(""),
      nName(""),
      index(-1),
      given(false),
      area(0.0),
      currentSum(0.0),
      //stateC(-1),
      chargeSum(0.0),
      numBoundaryPoints(0),
      Vckt(0.0),
      Vckt_orig(0.0),
      dIdVckt(0.0),
      dQdVckt(0.0),
      numCrossTerms(0),
      dxdvPtr(NULL),
      dxdvAllocated(false),
      Vckt_old(0.0),
      Vckt_final(0.0),
      Vckt_delta(0.0),
      Vckt_deltaC(0.0),
      Vckt_ramp(0.0),
      Vckt_ramp_old(0.0),
      neumannBCFlagV(false),
      neumannBCFlagN(false),
      neumannBCFlagP(false),
      material ("neutral"),
      materialGiven(false),
      oxideBndryFlag(false),
      oxthick(0.0),
      oxcharge(0.0)
  {
    Vcol.reserve(20);
    Ncol.reserve(20);
    Pcol.reserve(20);
    areaVector.reserve(20);
  };

public:
  std::string eName;     // device electrode name.  Should be in all upper case.
  std::string nName;     // circuit node name.

  int    index;     // index w.r.t. the order it was listed in the netlist.

  int   labelIndex; // index into the label array.

  int lid;          // local ID of the circuit node.
  int lidOffset;    // diagonal col in the ckt node matrix row. DMA only.

  std::vector<int> crossOffsets; // columns for the other ckt nodes.

  bool   given;     // did the user specify this node or not.  If not, probably gnd.

  std::vector<int> Vcol;  // column array
  std::vector<int> Ncol;  // column array
  std::vector<int> Pcol;  // column array

  // geometrical stuff:
  double area;                // total area for the edge.
  std::vector<double> areaVector;  // area for each edge node

  // state variable stuff:
  double currentSum;    // sum of currents, to be stored as a state variable.
  //int    stateC;        // state variable index, current.

  // capacitance related stuff
  double chargeSum;     // total charge on the electrode.

  //local id's (offsets)
  int li_stateC;

  // number of mesh points for this boundary.
  int numBoundaryPoints;

  // from the circuit node:
  double Vckt;

  // Vckt, before voltage limiting.
  double Vckt_orig;

  // derivative information needed for the 2-level Newton.
  double dIdVckt; // derivative of currentSum w.r.t. Vckt.
  double dQdVckt; // derivative of chargeSum w.r.t. Vckt.

  std::vector<double> dFdVckt;         // deriv. of residual w.r.t. Vckt.
  std::vector<int>    neighborNodes;   // nodes neighboring this boundary.

  std::vector<double> dQdX;     // deriv. of chargeSum w.r.t. PDE solution vars.
  std::vector<double> dIdX;     // deriv. of currentSum w.r.t. PDE solution vars.
  std::vector<int>    dIdXcols; // nodes neighboring this boundary.

  std::vector<int>    dIdXoffset; // If running with DMA, use this instead of
  // dIdXcols.

  int numCrossTerms;

  // equilibrium voltage at this electrode (if the applied voltage is
  // zero, there is still an internal voltage drop).
  std::vector<double> VequVec;

  // boundary conditions to be imposed on V,n and p.
  std::vector<double> VbcVec;
  std::vector<double> nnbcVec;
  std::vector<double> npbcVec;

  // dxdv vector:
  Linear::Vector * dxdvPtr;
  bool dxdvAllocated;

  // this map is between mesh node ID's and the indices of VbcVec
  std::map<int,int> meshGlobalToLocal;

  // These BC variables are only used for Continuation NL solves.
  // first 3 are w.r.t change in Vckt between ckt iterations.  Big change
  double Vckt_old;    // ckt value from the previous Newton solve.
  double Vckt_final;  // eventual ckt voltage value for the current solve.
  double Vckt_delta;  // total change in the Vckt. (Vckt_final-Vckt_old)

  // next 3, w.r.t change between continuation solves. small change
  double Vckt_deltaC; // incremental, intermediate change in the Vckt.
  double Vckt_ramp;   // current intermediate value of Vckt used during
  // continuation.
  double Vckt_ramp_old;   // old intermediate value of Vckt used during
  // continuation.

  // neumann BC flag;
  bool neumannBCFlagV;
  bool neumannBCFlagN;
  bool neumannBCFlagP;

  // material information:
  std::string material;
  bool   materialGiven;
  bool   oxideBndryFlag;
  double oxthick;
  double oxcharge;
};

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_DeviceInterfaceNode_h

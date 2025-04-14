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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 08/11/2010
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------


// ----------   Xyce Includes   ----------
#include <N_DEV_MembraneCS.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_DEV_SolverState.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : MembraneCS::MembraneCS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
MembraneCS::MembraneCS (const SolverState & ss1 ) : MembraneModel(ss1)
{
  // Connor Stevens has the following unknowns for the membrane
  // 1.  voltage
  // 2.  n
  // 3.  m
  // 4.  h
  // 5.  a
  // 6.  b
  // 7.  M
  // 8.  H
  // 9.  c
  // 10. Ca

  numIndependentVars_ = 10;
}


//-----------------------------------------------------------------------------
// Function      : MembraneCS::setJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneCS::setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp )
{
  int offset = numExtVars + numIndependentVars_*segmentNumber;
  int jacobianRowSize = segmentJacStamp[offset].size();

  //
  // jacobian element count by row:
  // vin:  2
  // vout: 2
  // v1:   3
  //
  // i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //
  // Vin   : g(0,1) * (V(1) - Vin ) = 0
  // Vout  : g(n,n-1) * (V(n-1) - Vout) = 0
  // Vnode : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //        plus node supporting equations (a, b, m)

  // jacobian format for full Connor Stevens Model
  //             Vin      Vout      V1   n   m   h  a  b  M  H  c  Ca  V2   n   m   h    V(nSeg)   n   m   h
  // kcl Vin    -g(0,1)            g(0,1)
  // kcl Vout           -g(n,n-1)                                      g(n,n-1)
  // kcl V1     yes                yes yes yes yes  y  y  y  y  y     yes
  // n                             yes yes
  // m                             yes     yes
  // h                             yes         yes
  // a                             yes             yes
  // b                             yes                yes
  // M                             yes                   yes
  // H                             yes                     yes
  // c                             yes                        yes yes
  // Ca                            yes                   yes yes  yes
  //
  //
  // jacobian element count by row:
  // vin:  2
  // vout: 2
  // v1:   11
  // n1:   2
  // m1:   2
  // h1:   2
  // a1:   2
  // b1:   2
  // M1:   2
  // H1:   2
  // c1:   3
  // Ca1:  4
  //

  // this resize should already have been done ?
  //segmentJacStamp[offset].resize(11);

  // This is for Vprev, Cable handles that
  /*
  if( i == 2 )
  {
    segmentJacStamp[offset][0] = 0;     // v_in
  }
  else
  {
    segmentJacStamp[offset][0] = offset-10;  // v_prev
  }
  */
  segmentJacStamp[offset][1] = offset;    // v
  segmentJacStamp[offset][2] = offset+1;  // n
  segmentJacStamp[offset][3] = offset+2;  // m
  segmentJacStamp[offset][4] = offset+3;  // h
  segmentJacStamp[offset][5] = offset+4;  // a
  segmentJacStamp[offset][6] = offset+5;  // b
  segmentJacStamp[offset][7] = offset+6;  // M
  segmentJacStamp[offset][8] = offset+7;  // H
  segmentJacStamp[offset][9] = offset+8;  // c

  // this is for Vnext, Cable handles that
  /*
  if( offset==(numVars-10) )
  {
    segmentJacStamp[offset][10] = 1;  // v_out
  }
  else
  {
    segmentJacStamp[offset][10] = offset+10; // v_next
  }
  */
  segmentJacStamp[offset+1].resize(2);   // n
  segmentJacStamp[offset+1][0] = offset;
  segmentJacStamp[offset+1][1] = offset+1;
  segmentJacStamp[offset+2].resize(2);   // m
  segmentJacStamp[offset+2][0] = offset;
  segmentJacStamp[offset+2][1] = offset+2;
  segmentJacStamp[offset+3].resize(2);   // h
  segmentJacStamp[offset+3][0] = offset;
  segmentJacStamp[offset+3][1] = offset+3;
  segmentJacStamp[offset+4].resize(2);   // a
  segmentJacStamp[offset+4][0] = offset;
  segmentJacStamp[offset+4][1] = offset+4;
  segmentJacStamp[offset+5].resize(2);   // b
  segmentJacStamp[offset+5][0] = offset;
  segmentJacStamp[offset+5][1] = offset+5;
  segmentJacStamp[offset+6].resize(2);   // M
  segmentJacStamp[offset+6][0] = offset;
  segmentJacStamp[offset+6][1] = offset+6;
  segmentJacStamp[offset+7].resize(2);   // H
  segmentJacStamp[offset+7][0] = offset;
  segmentJacStamp[offset+7][1] = offset+7;
  segmentJacStamp[offset+8].resize(3);   // c
  segmentJacStamp[offset+8][0] = offset;
  segmentJacStamp[offset+8][1] = offset+8;
  segmentJacStamp[offset+8][2] = offset+9;
  segmentJacStamp[offset+9].resize(4);   // ca
  segmentJacStamp[offset+9][0] = offset;
  segmentJacStamp[offset+9][1] = offset+6;
  segmentJacStamp[offset+9][2] = offset+7;
  segmentJacStamp[offset+9][3] = offset+9;
}

//-----------------------------------------------------------------------------
// Function      : MembraneCS::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneCS::loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeQVecPtr, double segArea)
{
}

//-----------------------------------------------------------------------------
// Function      : MembraneCS::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneCS::loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeFVecPtr, double segArea)
{
}

//-----------------------------------------------------------------------------
// Function      : MembraneCS::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneCS::loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, Linear::Vector * solnVecPtr, Linear::Matrix * dQdxMatPtr, double segArea)
{
}

//-----------------------------------------------------------------------------
// Function      : MembraneCS::loadDAEdFdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneCS::loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, Linear::Vector * solnVecPtr, Linear::Matrix * dFdxMatPtr, double segArea)
{
}

} // namespace Device
} // namespace Xyce

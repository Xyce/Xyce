//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
#include <N_DEV_Neuron_CommonEquations.h>
#include <N_DEV_MembraneHH.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_DEV_SolverState.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : MembraneHH::MembraneHH
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
MembraneHH::MembraneHH (const SolverState & ss1, double cMem, double gMem, double vRest, double eK, double gK, double eNa, double gNa):
  MembraneModel(ss1),
  cMem_(cMem),
  gMem_(gMem),
  vRest_(vRest),
  eK_(eK),
  gK_(gK),
  eNa_(eNa),
  gNa_(gNa)
{
  // Hodgkin-Huxley has the following unknows for the membrane
  // 1. voltage
  // 2. n
  // 3. m
  // 4. h
  numIndependentVars_ = 4;
}


//-----------------------------------------------------------------------------
// Function      : MembraneHH::setJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneHH::setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp )
{
  // caller sets up size of row of jac stamp to numIndependentVars_ + extra's needed for
  // its modeling.  So for a cable based device this is Vpre, Vseg, (other membrane vars), Vnext.
  // in general this is numIndependentVars_ + 2 (Vpre and Vnext).  In this routine we fill in
  // just what is needed for the membrane model

  int offset = numExtVars + numIndependentVars_*segmentNumber;
  int jacobianRowSize = segmentJacStamp[offset].size();

  // note, any external vars come before the internal vars.  Dependance on Vprev and Vnext (if
  // they are really there are handled by the caller as in a cable equation model.

  // Jacobian strcuture:
  //       Vpre    V     n     m     h   Vnext
  // kcl    yes   yes  yes   yes   yes   yes
  // n-equ        yes  yes
  // m-equ        yes        yes
  // h-equ        yes              yes

  // membrane voltage equation.
  segmentJacStamp[offset][vOffset    ] = offset;             // V   this should already have been set
  segmentJacStamp[offset][vOffset + 1] = offset+1;           // n
  segmentJacStamp[offset][vOffset + 2] = offset+2;           // m
  segmentJacStamp[offset][vOffset + 3] = offset+3;           // h

  // n equation
  segmentJacStamp[offset+1].resize(2);
  segmentJacStamp[offset+1][0] = offset;
  segmentJacStamp[offset+1][1] = offset+1;

  // m equation
  segmentJacStamp[offset+2].resize(2);
  segmentJacStamp[offset+2][0] = offset;
  segmentJacStamp[offset+2][1] = offset+2;

  // h equation
  segmentJacStamp[offset+3].resize(2);
  segmentJacStamp[offset+3][0] = offset;
  segmentJacStamp[offset+3][1] = offset+3;

}

//-----------------------------------------------------------------------------
// Function      : MembraneHH::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneHH::loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeQVecPtr, double segArea)
{
  // Each segment will have numIndependentVars_ with segment voltage being the first
  // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
  int index = segmentNumber * numIndependentVars_;
  (*daeQVecPtr)[lidIndexVector[index]] += cMem_ * segArea * (*solnVecPtr)[lidIndexVector[index]];
  // n equation
  (*daeQVecPtr)[lidIndexVector[index + 1]] += (*solnVecPtr)[lidIndexVector[index + 1]];
  // m equation
  (*daeQVecPtr)[lidIndexVector[index + 2]] += (*solnVecPtr)[lidIndexVector[index + 2]];
  // h equation
  (*daeQVecPtr)[lidIndexVector[index + 3]] += (*solnVecPtr)[lidIndexVector[index + 3]];
}

//-----------------------------------------------------------------------------
// Function      : MembraneHH::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneHH::loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, Linear::Vector * solnVecPtr, Linear::Vector * daeFVecPtr, double segArea)
{
// Each segment will have numIndependentVars_ with segment voltage being the first
  // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.

  int index = segmentNumber * numIndependentVars_;
  double vSeg = (*solnVecPtr)[lidIndexVector[index]];
  double n = (*solnVecPtr)[lidIndexVector[index + 1]];
  double m = (*solnVecPtr)[lidIndexVector[index + 2]];
  double h = (*solnVecPtr)[lidIndexVector[index + 3]];

  // membrane current equation
  (*daeFVecPtr)[lidIndexVector[index]] += Neuron::HH_Vseg_F<double>(
    vSeg, n,  m,  h, (gMem_ * segArea), vRest_, (gK_ * segArea), eK_,
    (gNa_ * segArea), eNa_ );

  // cew 2/18/11 - changed the functions used below to use voltage shifted by resting potential
  // n equation
  (*daeFVecPtr)[lidIndexVector[index + 1]] += Neuron::nEquF<double>( vSeg-vRest_, n );

  // m equation
  (*daeFVecPtr)[lidIndexVector[index + 2]] += Neuron::mEquF<double>( vSeg-vRest_, m );

  // h equation
  (*daeFVecPtr)[lidIndexVector[index + 3]] += Neuron::hEquF<double>( vSeg-vRest_, h );

}

//-----------------------------------------------------------------------------
// Function      : MembraneHH::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneHH::loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector,
                                    std::vector< std::vector< int > > & jacobianOffsets,
                                    Linear::Vector * solnVecPtr,
                                    Linear::Matrix * dQdxMatPtr,
                                    double segArea)
{
   // while lidIndexVector lists LID's for just the segment variables (just V in the case
   // of a passive cable). The jacobianOffsets includes the Vin and Vout as the first
   // two variables.  Thus, there is a constant offset of 2 for everything in jacobianOffsets

   // And, as in the Q and F load functions,  Each segment will have numIndependentVars_ with segment voltage being the first
   // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
   // in the case of the passive cable numIndependentVars_=1.

   int index = segmentNumber * numIndependentVars_;
   int row = numExternalVars_ + index;               // numExternalVars_ a contant of 2 assumed in MembraneModel base class

   // Vseg equation
   (*dQdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset]] += cMem_ * segArea;

   // for internal variables remember that order if Vseg, n, m, h

   // n equation
   (*dQdxMatPtr)[lidIndexVector[index + 1]][jacobianOffsets[row + 1][1]] += 1;

   // m equation
   (*dQdxMatPtr)[lidIndexVector[index + 2]][jacobianOffsets[row + 2][1]] += 1;

   // h equation
   (*dQdxMatPtr)[lidIndexVector[index + 3]][jacobianOffsets[row + 3][1]] += 1;
}

//-----------------------------------------------------------------------------
// Function      : MembraneHH::loadDAEdFdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembraneHH::loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector,
                                    std::vector< std::vector< int > > & jacobianOffsets,
                                    Linear::Vector * solnVecPtr,
                                    Linear::Matrix * dFdxMatPtr,
                                    double segArea)
{
   // while lidIndexVector lists LID's for just the segment variables (just V in the case
   // of a passive cable). The jacobianOffsets includes the Vin and Vout as the first
   // two variables.  Thus, there is a constant offset of 2 for everything in jacobianOffsets

   // And, as in the Q and F load functions,  Each segment will have numIndependentVars_ with segment voltage being the first
   // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
   // in the case of the passive cable numIndependentVars_=1.

   int index = segmentNumber * numIndependentVars_;
   int row = numExternalVars_ + index;       // numExternalVars_ a contant of 2 assumed in MembraneModel base class

   double vSeg = (*solnVecPtr)[lidIndexVector[index]];
   double n = (*solnVecPtr)[lidIndexVector[index + 1]];
   double m = (*solnVecPtr)[lidIndexVector[index + 2]];
   double h = (*solnVecPtr)[lidIndexVector[index + 3]];

   // Since the Sacado types for auto differentiation are templated on the number of
   // derivatives needed, use scoping brackets here so I can safely reuse some variable names

   // Vseg equation
   {
     Sacado::Fad::SFad<double,4> vVar( 4, 0, vSeg );
     Sacado::Fad::SFad<double,4> nVar( 4, 1, n );
     Sacado::Fad::SFad<double,4> mVar( 4, 2, m );
     Sacado::Fad::SFad<double,4> hVar( 4, 3, h );
     // parameters
     Sacado::Fad::SFad<double,4> gMemVar( gMem_ * segArea );
     Sacado::Fad::SFad<double,4> vRestVar( vRest_ );
     Sacado::Fad::SFad<double,4> gKVar( gK_ * segArea );
     Sacado::Fad::SFad<double,4> eKVar( eK_ );
     Sacado::Fad::SFad<double,4> gNaVar( gNa_ * segArea );
     Sacado::Fad::SFad<double,4> eNaVar( eNa_ );

     Sacado::Fad::SFad<double,4> resultFad;
     resultFad = Neuron::HH_Vseg_F( vVar, nVar, mVar, hVar, gMemVar, vRestVar, gKVar, eKVar, gNaVar, eNaVar );

     (*dFdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset]] += resultFad.dx(0);      // /dVseg
     (*dFdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset+1]] += resultFad.dx(1);    // /dn
     (*dFdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset+2]] += resultFad.dx(2);    // /dm
     (*dFdxMatPtr)[lidIndexVector[index]][jacobianOffsets[row][vOffset+3]] += resultFad.dx(3);    // /dh
   }

   // cew 2/18/11 - changed the functions used below to use voltage shifted by resting potential
   // n equation
   {
     Sacado::Fad::SFad<double,2> vVar( 2, 0, vSeg-vRest_ );
     Sacado::Fad::SFad<double,2> nVar( 2, 1, n );
     Sacado::Fad::SFad<double,2> nEquResultFad = Neuron::nEquF( vVar, nVar );

     (*dFdxMatPtr)[lidIndexVector[index + 1]][jacobianOffsets[row + 1][0]] += nEquResultFad.dx(0);    // dVseg
     (*dFdxMatPtr)[lidIndexVector[index + 1]][jacobianOffsets[row + 1][1]] += nEquResultFad.dx(1);    // dn
   }

   // m equation
   {
     Sacado::Fad::SFad<double,2> vVar( 2, 0, vSeg-vRest_ );
     Sacado::Fad::SFad<double,2> mVar( 2, 1, m );
     Sacado::Fad::SFad<double,2> mEquResultFad = Neuron::mEquF( vVar, mVar );

     (*dFdxMatPtr)[lidIndexVector[index + 2]][jacobianOffsets[row + 2][0]] += mEquResultFad.dx(0);    // dVseg
     (*dFdxMatPtr)[lidIndexVector[index + 2]][jacobianOffsets[row + 2][1]] += mEquResultFad.dx(1);    // dm
   }

   // h equation
   {
     Sacado::Fad::SFad<double,2> vVar( 2, 0, vSeg-vRest_ );
     Sacado::Fad::SFad<double,2> hVar( 2, 1, h );
     Sacado::Fad::SFad<double,2> hEquResultFad = Neuron::hEquF( vVar, hVar );

     (*dFdxMatPtr)[lidIndexVector[index + 3]][jacobianOffsets[row + 3][0]] += hEquResultFad.dx(0);    // dVseg
     (*dFdxMatPtr)[lidIndexVector[index + 3]][jacobianOffsets[row + 3][1]] += hEquResultFad.dx(1);    // dh
   }

}

} // namespace Device
} // namespace Xyce

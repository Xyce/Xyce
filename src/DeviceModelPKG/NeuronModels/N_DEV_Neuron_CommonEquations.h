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

//-----------------------------------------------------------------------------
//
// Purpose        : Common, template based, equations for neuron devices
//                  Keeping them here to avoid duplication in various
//                  device files
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/02/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron_CommonEquations_h
#define Xyce_N_DEV_Neuron_CommonEquations_h

#include <Sacado_No_Kokkos.hpp>

namespace Xyce {
namespace Device {
namespace Neuron {

// These functions represents the equations that need to be solved
// for this device.  Since Xyce loads an F and Q contribution, the
// equations are broken up into their F and Q components.  Thus there
// is a kcl1EquF() and kcl1EquQ().  Automatic differentiation will
// be used to generate all the derivatives of these equations for the
// dF/dX and dQ/dX loads


//
// from membrane patch formulation (V1 - inner membrane, v2 - outer membrane)
//

// first we list some utility functions for calculating coefficients
//
// These functions expect V to be in milli-volts and then return values that
// are in 1/ms.  Thus the extra factor's of 1000 here and there

// cew 2/18/2011
// Commenting out the Connor Stevens alphaN, betaN, alphaM, betaM, alphaH, betaH
// and replacing them with Hodgkin-Huxley version, which is different.
// Right now, the Connor Stevens membrane model isn't even
// implemented, so those equations aren't needed yet.
// Copied the HH equations from N_Dev_Neuron.h, but modified them and corresponding
// code in N_DEV_MembraneHH.C to handle subtraction of Vrest in calling code.
// At some point, we have to re-evaluate whether it's possible to have a set
// of common equations multiple devices can use.

// equations from Hodgkin-Huxley model

// potassium current, functions for activator equation
template <typename ScalarT>
static ScalarT alphaN( const ScalarT & Vn1)
{
  ScalarT vDiff = 1000.0 * Vn1;  // convert voltage to milli-volts
  ScalarT r;
  if ((vDiff > 9.99) && (vDiff < 10.01) )
  {
    r = 1.0/(10.0 * ( std::exp( (10.0 - vDiff)/10.0 )));
  }
  else
  {
    r = (10.0 - vDiff) /
      (100.0 * ( std::exp( (10.0 - vDiff)/10.0 ) - 1.0 ));
  }
  r *= 1000.0; // change from 1/ms to 1/s
  return r;
}

template <typename ScalarT>
static ScalarT betaN( const ScalarT & Vn1)
{
  ScalarT vDiff = 1000.0 * Vn1;
  ScalarT r = 0.125 * std::exp( -vDiff/80.0 );
  r *= 1000.0; // change from 1/ms to 1/s
  return r;
}

// sodium current, functions for activator equation
template <typename ScalarT>
static ScalarT alphaM( const ScalarT & Vn1)
{
  ScalarT vDiff = 1000.0 * Vn1;
  ScalarT r;
  if ((vDiff > 24.99) && (vDiff < 25.01) )
  {
    r = (1.0) /
      (( std::exp( (25.0 - vDiff)/10.0 )));
  }
  else
  {
    r = (25.0 - vDiff) /
      (10.0 * ( std::exp( (25.0 - vDiff)/10.0 ) - 1.0 ));
  }
  r *= 1000.0; // change from 1/ms to 1/s
  return r;
}

template <typename ScalarT>
static ScalarT betaM( const ScalarT & Vn1)
{
  ScalarT vDiff = 1000.0 * Vn1;
  ScalarT r = 4.0 * std::exp( -vDiff/18.0 );
  r *= 1000.0; // change from 1/ms to 1/s

  return r;
}

template <typename ScalarT>
static ScalarT alphaH( const ScalarT & Vn1)
{
  ScalarT vDiff = 1000.0 * Vn1;
  ScalarT r = 0.07 * std::exp( -vDiff/20.0 );
  r *= 1000.0; // change from 1/ms to 1/s
  return r;
}

template <typename ScalarT>
static ScalarT betaH( const ScalarT & Vn1)
{
  ScalarT vDiff = 1000.0 * Vn1;
  ScalarT r = 1.0 / ( std::exp( (30.0 - vDiff)/10.0 ) + 1.0 );
  r *= 1000.0; // change from 1/ms to 1/s
  return r;
}



// equations from the ConnorStevens Model
/*
// potassium current, functions for activator equation
template <typename ScalarT>
static ScalarT alphaN( const ScalarT Vin)
{
ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
ScalarT r = 1000.0 * (0.02 * (vScaled + 45.7)) / (1.0 - std::exp(-0.1*(vScaled+45.7)));
// result.  the 1000 factor is to change from 1/ms to 1/s
return r;
}

template <typename ScalarT>
static ScalarT betaN( const ScalarT Vin)
{
ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
ScalarT r = 1000.0 * 0.25 * std::exp( -0.0125 * (vScaled + 55.7));
// result.  the 1000 factor is to change from 1/ms to 1/s
return r;
}

// sodium current, functions for activator equation
template <typename ScalarT>
static ScalarT alphaM( const ScalarT Vin)
{
ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
ScalarT r = 1000.0 * (0.38 * (vScaled + 29.7)) / (1.0 - std::exp(-0.1*(vScaled+29.7)));
// result.  the 1000 factor is to change from 1/ms to 1/s
return r;
}

template <typename ScalarT>
static ScalarT betaM( const ScalarT Vin)
{
ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
ScalarT r = 1000.0 * 15.2 * std::exp( -0.0556 * (vScaled + 54.7));
// result.  the 1000 factor is to change from 1/ms to 1/s
return r;
}

template <typename ScalarT>
static ScalarT alphaH( const ScalarT Vin)
{
ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
ScalarT r = 1000.0 * 0.266 * std::exp( -0.05 * (vScaled + 48.0));
// result.  the 1000 factor is to change from 1/ms to 1/s
return r;
}

template <typename ScalarT>
static ScalarT betaH( const ScalarT Vin)
{
ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
ScalarT r = 1000.0 * 3.8 / (1.0 + std::exp(-0.1*(vScaled+18.0)));
// result.  the 1000 factor is to change from 1/ms to 1/s
return r;
}
*/
// a-current functions
template <typename ScalarT>
static ScalarT aInf( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = std::pow( ((0.0761 * std::exp(0.0314 * (vScaled+94.22))) / (1.0+std::exp(0.0346*(vScaled+1.17)))), 1.0/3.0);
  return r;
}

template <typename ScalarT>
static ScalarT aTau( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = (0.3632 + 1.158 / (1.0 + std::exp(0.0497 * (vScaled + 55.96)))) / 1000.0;
  // the 1000.0 factor is to change ms to s.
  return r;
}

template <typename ScalarT>
static ScalarT bInf( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = std::pow( (1.0 / (1.0 + std::exp(0.0688*(vScaled+53.3)))), 4.0);
  return r;
}

template <typename ScalarT>
static ScalarT bTau( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = (1.24 + 2.678 / (1.0 + std::exp(0.0624 * (vScaled + 50.0)))) / 1000.0;
  return r;
}

// transient calcium functions
template <typename ScalarT>
static ScalarT M_Inf( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = 1.0/(1.0 + std::exp(-(vScaled+57)/6.2));
  return r;
}

template <typename ScalarT>
static ScalarT M_Tau( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = (0.612 + 1.0/(std::exp(-(vScaled+132)/16.7) + std::exp((vScaled+16.8)/18.2)) ) / 1000.0;
  return r;
}

template <typename ScalarT>
static ScalarT H_Inf( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = 1.0 / (1.0 + std::exp((vScaled+81)/4.0));
  return r;
}

template <typename ScalarT>
static ScalarT H_Tau( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r;
  if( vScaled < -80.0 )
  {
    r = std::exp( (vScaled + 467)/66.6 ) / 1000.0;
  }
  else
  {
    r = ( 28.0 + std::exp(-(vScaled+22.0)/10.5)) / 1000.0;
  }
  return r;
}

// Calcium dependent Potassium conductances
template <typename ScalarT>
static ScalarT C_Inf( const ScalarT Vin, const ScalarT CaConc)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = (CaConc / (CaConc + 3.0)) * (1.0 / (1.0 + std::exp(-(vScaled+28.3)/12.6 )));
  return r;
}

template <typename ScalarT>
static ScalarT C_Tau( const ScalarT Vin)
{
  ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
  ScalarT r = (90.3 - 75.1/(1.0 + std::exp(-(vScaled+46)/22.7))) / 1000.0;
  return r;
}

// now the device equations
// KCL equation 1
template <typename ScalarT>
static ScalarT kcl1EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                         const ScalarT& a, const ScalarT& b, const ScalarT& MC, const ScalarT& HC, const ScalarT& CC,
                         const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE,
                         const ScalarT& Ag, const ScalarT& Ae, const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& KCaG)
{
  ScalarT powN = n * n * n * n;
  ScalarT powM = m * m * m;
  ScalarT powA = a * a * a;
  ScalarT powMC = MC * MC;
  ScalarT powCC = CC * CC * CC * CC;
  ScalarT r = memG * (Vn1 - Vn2 - restV) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE )
    + Ag * powA * b * (Vn1 - Vn2 - Ae) + CaTg * powMC * HC * (Vn1 - Vn2 - CaE) + KCaG * powCC * (Vn1 - Vn2 - Ke);
  return r;
}

template <typename ScalarT>
static ScalarT kcl1EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
{
  ScalarT r = memC * (Vn1 - Vn2);
  return r;
}

// KCL equation 2 -- -1 * equation 1 because of device symmetry
template <typename ScalarT>
static ScalarT kcl2EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                         const ScalarT& a, const ScalarT& b, const ScalarT& MC, const ScalarT& HC, const ScalarT& CC,
                         const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE,
                         const ScalarT& Ag, const ScalarT& Ae, const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& KCaG)
{
  ScalarT powN = n * n * n * n;
  ScalarT powM = m * m * m;
  ScalarT powA = a * a * a;
  ScalarT powMC = MC * MC;
  ScalarT powCC = CC * CC * CC * CC;
  ScalarT r = -1.0 * (memG * (Vn1 - Vn2 - restV) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE )
                      + Ag * powA * b * (Vn1 - Vn2 - Ae) + CaTg * powMC * HC * (Vn1 - Vn2 - CaE) + KCaG * powCC * (Vn1 - Vn2 - Ke) );
  return r;
}

template <typename ScalarT>
static ScalarT kcl2EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
{
  ScalarT r = -1.0 * memC * (Vn1 - Vn2);
  return r;
}

// Hodgkin Huxley Voltage equation for a membrane segment
// KCL equation 1
template <typename ScalarT>
static ScalarT HH_Vseg_F( const ScalarT& Vseg,const ScalarT& n, const ScalarT& m, const ScalarT& h,
                          const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE )
{
  ScalarT n4 = n * n * n * n;
  ScalarT m3 = m * m * m;
  ScalarT r = memG * (Vseg - restV) + Kg * n4 * (Vseg - Ke ) + NaG * m3 * h * (Vseg - NaE );
  return r;
}

#if 0
template <typename ScalarT>
static ScalarT HH_Vseg_Q( const ScalarT& Vseg, const ScalarT& memC )
{
  ScalarT r = memC * Vn1;
  return r;
}
#endif

// n conservation equation
template <typename ScalarT>
static ScalarT nEquF( const ScalarT& Vn, const ScalarT& n)
{
  ScalarT alpha = alphaN<ScalarT>(Vn);
  ScalarT beta = betaN<ScalarT>(Vn);
  ScalarT r = (alpha + beta) * n - alpha;
  return r;
}

template <typename ScalarT>
static ScalarT nEquQ( const ScalarT& n )
{
  ScalarT r = n;
  return r;
}

// m conservation equation
template <typename ScalarT>
static ScalarT mEquF( const ScalarT& Vn, const ScalarT& m )
{
  ScalarT alpha = alphaM<ScalarT>(Vn);
  ScalarT beta = betaM<ScalarT>(Vn);
  ScalarT r = (alpha + beta) * m - alpha;
  return r;
}

template <typename ScalarT>
static ScalarT mEquQ( const ScalarT& m )
{
  ScalarT r = m;
  return r;
}

// h conservation equation
template <typename ScalarT>
static ScalarT hEquF( const ScalarT& Vn, const ScalarT& h )
{
  ScalarT alpha = alphaH<ScalarT>(Vn);
  ScalarT beta = betaH<ScalarT>(Vn);
  ScalarT r = (alpha + beta) * h - alpha;
  return r;
}

template <typename ScalarT>
static ScalarT hEquQ( const ScalarT& h )
{
  ScalarT r = h;
  return r;
}

// a conservation equation
template <typename ScalarT>
static ScalarT aEquF( const ScalarT& Vn1, const ScalarT& a, const ScalarT& Vrest )
{
  ScalarT vDiff = Vn1; // - Vrest
  ScalarT Inf = aInf<ScalarT>(vDiff);
  ScalarT Tau = aTau<ScalarT>(vDiff);
  ScalarT r = (a - Inf)/Tau;
  return r;
}

template <typename ScalarT>
static ScalarT aEquQ( const ScalarT& a )
{
  ScalarT r = a;
  return r;
}

// b conservation equation
template <typename ScalarT>
static ScalarT bEquF( const ScalarT& Vn1, const ScalarT& b, const ScalarT& Vrest )
{
  ScalarT vDiff = Vn1; // - Vrest
  ScalarT Inf = bInf<ScalarT>(vDiff);
  ScalarT Tau = bTau<ScalarT>(vDiff);
  ScalarT r = (b - Inf)/Tau;
  return r;
}

template <typename ScalarT>
static ScalarT bEquQ( const ScalarT& b )
{
  ScalarT r = b;
  return r;
}

// M conservation equation
template <typename ScalarT>
static ScalarT M_EquF( const ScalarT& Vn1, const ScalarT& M, const ScalarT& Vrest )
{
  ScalarT vDiff = Vn1; // - Vrest
  ScalarT Inf = M_Inf<ScalarT>(vDiff);
  ScalarT Tau = M_Tau<ScalarT>(vDiff);
  ScalarT r = (M - Inf)/Tau;
  return r;
}

template <typename ScalarT>
static ScalarT M_EquQ( const ScalarT& M )
{
  ScalarT r = M;
  return r;
}

// H conservation equation
template <typename ScalarT>
static ScalarT H_EquF( const ScalarT& Vn1, const ScalarT& H, const ScalarT& Vrest )
{
  ScalarT vDiff = Vn1; // - Vrest
  ScalarT Inf = H_Inf<ScalarT>(vDiff);
  ScalarT Tau = H_Tau<ScalarT>(vDiff);
  ScalarT r = (H - Inf)/Tau;
  return r;
}

template <typename ScalarT>
static ScalarT H_EquQ( const ScalarT& H )
{
  ScalarT r = H;
  return r;
}

// C conservation equation
template <typename ScalarT>
static ScalarT C_EquF( const ScalarT& Vn1, const ScalarT& C, const ScalarT& CaConc, const ScalarT& Vrest )
{
  ScalarT vDiff = Vn1; // - Vrest
  ScalarT Inf = C_Inf<ScalarT>(vDiff, CaConc);
  ScalarT Tau = C_Tau<ScalarT>(vDiff);
  ScalarT r = (C - Inf)/Tau;
  return r;
}

template <typename ScalarT>
static ScalarT C_EquQ( const ScalarT& C )
{
  ScalarT r = C;
  return r;
}

// Calcium conservation equation
template <typename ScalarT>
static ScalarT Ca_EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& MC, const ScalarT& HC, const ScalarT& Ca,
                        const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& CaGamma, const ScalarT& CaTau )
{
  ScalarT r = CaGamma * CaTg * MC * MC * HC * (Vn1 - Vn2 - CaE) + Ca / CaTau;
  return r;
}

template <typename ScalarT>
static ScalarT Ca_EquQ( const ScalarT& Ca)
{
  ScalarT r = Ca;
  return r;
}

//
// from membrane cable formulation (Vin - inner membrane outer membrane is considered to be zero)
//
// equations from the ConnorStevens Model


// now the device equations
// KCL equation 1
template <typename ScalarT>
static ScalarT kcl1EquF( const ScalarT& VSeg, const ScalarT& VSegP, const ScalarT & VSegN, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                         const ScalarT& a, const ScalarT& b, const ScalarT& MC, const ScalarT& HC, const ScalarT& CC,
                         const ScalarT& gPrev, const ScalarT& gNext,
                         const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE,
                         const ScalarT& Ag, const ScalarT& Ae, const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& KCaG)
{
  ScalarT powN = n * n * n * n;
  ScalarT powM = m * m * m;
  ScalarT powA = a * a * a;
  ScalarT powMC = MC * MC;
  ScalarT powCC = CC * CC * CC * CC;
  ScalarT r = memG * (VSeg - restV)  - gNext * (VSegN - VSeg) - gPrev * (VSegP - VSeg);
  //ScalarT r = memG * (VSeg - restV) + Kg * powN * (VSeg - Ke ) + NaG * powM * h * (VSeg - NaE )
  //    + Ag * powA * b * (VSeg - Ae) + CaTg * powMC * HC * (VSeg - CaE) + KCaG * powCC * (VSeg - Ke) - gNext * (VSegN - VSeg) - gPrev * (VSegP - VSeg);
  return r;
}

template <typename ScalarT>
static ScalarT kcl1EquQ( const ScalarT& VSeg, const ScalarT& memC )
{
  ScalarT r = memC * VSeg;
  return r;
}


// Calcium conservation equation
template <typename ScalarT>
static ScalarT Ca_EquF( const ScalarT& Vn1, const ScalarT& MC, const ScalarT& HC, const ScalarT& Ca,
                        const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& CaGamma, const ScalarT& CaTau )
{
  ScalarT r = CaGamma * CaTg * MC * MC * HC * (Vn1 - CaE) + Ca / CaTau;
  return r;
}

} // namespace Neuron
} // namespace Device
} // namespace Xyce

#endif

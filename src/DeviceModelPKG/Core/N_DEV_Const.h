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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/17/01
//
//
//
//
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_Const_h
#define Xyce_N_DEV_Const_h

#define CONSTroot2   sqrt(2.0)
#define CONSTQ       (1.6021918e-19)     // electron charge
#define CONSTCtoK    (273.15)            // conversion of Celsius to Kelvin
#define CONSTboltz   (1.3806226e-23)     // Boltzmann's constant
#define CONSTplanck  (6.62606957e-34)    // Planck's constant  (in J-s)
#define CONSTemass   (9.10938291e-31)    // e- mass in kg.
#define CONSTREFTEMP (300.15)            // 27 degrees C

#define CONSTKoverQ  (CONSTboltz/CONSTQ) // K_b/q

#define CONSTe       exp(1.0)
#define CONSTclight  (2.99792458e8)      // speed of light m/s
#define CONSTpermFS   (4.0* M_PI * 1.0e-7 )   // permeability of free space -- usually written as mu0
                                              // not the same as the permitivity of free space defined next.
#define CONSTperm0   (8.854214871e-12)   // permitivity of free space.
#define CONSTpermOxide (3.9*CONSTperm0) 

#define CONSTvt0     (CONSTboltz * (27.0 +CONSTCtoK)/CONSTQ)

// these two constants are used to caluculate the temperature dependent
// bandgap for Si.

#define CONSTEg300   (1.1150877) // band gap for Si at T=300.15K (room temp)
#define CONSTEg0     (1.16)      // band gap for Si at T=0K. (eV)
#define CONSTalphaEg (7.02e-4)   // (eV/K)
#define CONSTbetaEg  (1108.0)    // (K)

#define CONSTNi0     (1.45e10)   // carrier concentration at room temp.

#define CONSTNMOS 1
#define CONSTPMOS -1

#define CONSTEPSOX (3.453133e-11)
#define CONSTEPSSI (1.03594e-10)

#define CONSTEPSSIL (11.7 * 8.854214871e-12)


// BSIM3 constants:
// Most of these exponential constants are taken from 3f5's BSIM3 model.
// In reality they may be machine dependent.
//
// Normally machine dependent limits which are used in 3f5 come from
// files such as limits.h, hw.h, hw_ieee.h, and/or hw_vax.h.  limits.h
// is a file that comes with the compiler, but the rest are part of 3f5.
//
#define CONSTMAX_EXPL 2.688117142e+43
#define CONSTMIN_EXPL 3.720075976e-44
#define CONSTEXPL_THRESHOLD 100.0

#define CONSTEXP_THRESHOLD 34.0
#define CONSTMAX_EXP 5.834617425e14
#define CONSTMIN_EXP 1.713908431e-15

#define CONSTDELTA_1 0.02
#define CONSTDELTA_2 0.02
#define CONSTDELTA_3 0.02
#define CONSTDELTA_4 0.02

// end of the BSIM3 constants.

// This value, 100.0 was arbitrarily chosen, mainly to
// avoid ieee errors when running Xyce without spice-style
// voltage limiters.
#define CONSTMAX_EXP_ARG 100.0

// tanh threshold constant (value for which tanh(x) - 1.0 == 0.0) - this seems
// to vary from machine to machine but the largest I found was 22.0 on g++2.96
// with -O2 (it was less with no optimization).  I've chose 20.0 as a
// comprimise.  SAH, 14 April 2003

#define CONSTTANH_THRESH 20.0

#endif


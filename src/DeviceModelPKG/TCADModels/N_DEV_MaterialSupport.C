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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/19/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_MaterialSupport.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace Device {

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getEffectiveMassN
// Purpose       : returns effective mass for electrons.
// Special Notes : Relative to free space mass.
//
//                 These are from Appendix 3 of Streetman.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getEffectiveMassN (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();

  double mass=0.0;

  if (mater == "si")
  {
    double ml = 0.98; // longitudinal mass
    double mt = 0.19; // transverse mass
    mass = pow((ml*mt*mt),1.0/3.0);
  }
  else if (mater == "ge" )
  {
    double ml = 1.64; // longitudinal mass
    double mt = 0.082; // transverse mass
    mass = pow((1.64*0.082*0.082),1.0/3.0);
  }
  else if (mater == "gaas")
  {
    //mass = 0.067;
    mass = 6.69935094e-2;
  }
  else if (mater=="inalas" || mater=="alinas")
  {
    mass = 0.074;
    //dnco=  0.020;
  }
  else if (mater=="ingaas" || mater=="gainas")
  {
    mass = 0.041;
    // dnco=  0.0083;
  }
  else if (mater == "ingap")
  {
    mass = 0.0179;
    // dnco =  0.02391;
  }
  else if (mater == "inp")
  {
    mass = 0.079;
    // dnco=  0.01734;
  }
  else
  {
    Report::UserFatal() << material << " material not recognized in getEffectiveMassN.";
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getEffectiveMassP
// Purpose       : returns effective mass for holes.
// Special Notes : Relative to free space mass.
//
//                 These are from Appendix 3 of Streetman.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getEffectiveMassP (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass=0.0;

  if (mater == "si")
  {
    double mlh = 0.16; // light hole mass
    double mhh = 0.49; // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh,1.5)),2.0/3.0);
  }
  else if (mater == "ge" )
  {
    double mlh = 0.04; // light hole mass
    double mhh = 0.28; // heavy hole mass
    mass = pow((pow(mlh, 1.5) + pow(mhh, 1.5)),2.0/3.0);
  }
  else if (mater == "gaas")
  {
    // Note: Wampler's value of dnva = 0.4318
    // Their SAND report uses a DOS effective mass as: mh*=0.571
    // Note that the electron effective mass for GaAs seems to be uncontroversial m=0.067
    // In order to match XPD, and pass historical regression tests, etc, using this number.
    //mass=0.571; 
    mass=5.71287987e-1;
  }
  else if (mater=="inalas" || mater=="alinas")
  {
    double mlh = 0.08; // light hole mass
    double mhh = 0.6;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else if (mater=="ingaas" || mater=="gainas")
  {
    double mlh = 0.05; // light hole mass
    double mhh = 0.54; // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else if (mater=="ingap")
  {
    // dnva =  0.5423;
    //mass = 0.665;
    mass = 6.65007285e-1;
  }
  else if (mater == "gan")
  {
    mass = 0.8; // http://www.ioffe.ru/SVA/NSM/Semicond/GaN/basic.html
  }
  else if (mater == "inp")
  {
    double mlh = 0.074; // light hole mass
    double mhh = 0.5;   // heavy hole mass
    mass = pow((pow(mlh,1.5) + pow(mhh, 1.5)), 2.0/3.0);
  }
  else
  {
    Report::UserFatal() << material << " material not recognized in getEffectiveMassP.";
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::get_DOS_EffectiveMassN
// Purpose       : returns effective mass for electrons for density of states
// Special Notes : See http://ecee.colorado.edu/~bart/book/effmass.htm#dosmass.
//
//   m_e_dos = Mc^(2/3) * (m_l*m_t*m_t)^(1/3)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
double MaterialSupport::get_DOS_EffectiveMassN (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();
  double mass;

  if (mater == "si")
  {
    double ml = 0.98; // longitudinal mass
    double mt = 0.19; // transverse mass
    double Mc = 6.0; // degeneracy factor (number of equivalent band minimums)
    mass = pow(Mc,2.0/3.0)*pow((ml*mt*mt),1.0/3.0);
    // Note, this should evaluate to 1.08.
  }
  else if (mater == "ge")
  {
    double ml = 1.64; // longitudinal mass
    double mt = 0.082; // transverse mass
    double Mc = 4.0; // degeneracy factor 
    mass = pow(Mc,2.0/3.0)*pow((ml*mt*mt),1.0/3.0);
    // this mass should be around 0.56
  }
  else if (mater == "gaas")
  {
    mass = 0.067; // GaAs is isotropic
  }
  else if (mater=="inalas" || mater=="alinas")
  {
    //dnco=  0.020; = mass^(1.5)
    mass = 0.074;
  }
  else if (mater=="ingaas" || mater=="gainas")
  {
    // dnco=  0.0083; = mass^(1.5)
    mass = 0.041;
  }
  else if (mater=="ingap")
  {
    // dnco =  0.02391;
    //mass = 0.0179;
    mass = 8.29952143e-2;
  }
  else if (mater == "gan")
  {
    mass = 0.2; // http://www.ioffe.ru/SVA/NSM/Semicond/GaN/basic.html
  }
  else if (mater == "inp")
  {
    mass = 0.079;
  }
  else
  {
    Report::UserFatal() << material << " material not recognized get_DOS_EffectiveMassN.";
  }

  return mass;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::get_DOS_EffectiveMassP
// Purpose       : returns effective mass for holes for density of states
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
// ----------------------------------------------------------------------------
double MaterialSupport::get_DOS_EffectiveMassP (const std::string & material)
{
  // Unlike for electrons, Mc is not applied to holes, so this is simply 
  // the same as the non-DOS version.
  return getEffectiveMassP (material);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNc
// Purpose       : density of states, conduction band
// Special Notes :
//
//    Nc = 2*(2*pi*me*kT/h^2)^(3/2)
//       = 2*(2*pi*me_DOS*m0*kT/h^2)^(3/2)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/10/2014
// ----------------------------------------------------------------------------
double MaterialSupport::getNc (const std::string & material, double temp)
{
  double h_planck(CONSTplanck); // Planck's constant  (in J-s)
  double e_mass (CONSTemass);   // e- mass in kg.
  double kb (1.3806488e-23); // boltzmann's constant (J/K)
  double meDOS = get_DOS_EffectiveMassN(material);
  double dnbnd0 = 2.0*M_PI*meDOS*e_mass*kb*temp/(h_planck*h_planck);
  double Nc = 2.0*pow(dnbnd0,1.5)/1.0e6;

  return Nc;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNv
// Purpose       : Density of states, valance band.
// Special Notes :
//
//     Nv = 2*(2*pi*mh*kT/h^2)^(3/2)
//       = 2*(2*pi*mh_DOS*m0*kT/h^2)^(3/2)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/10/2014
// ----------------------------------------------------------------------------
double MaterialSupport::getNv (const std::string & material, double temp)
{
  double h_planck(CONSTplanck); // Planck's constant  (in J-s)
  double e_mass (CONSTemass);   // e- mass in kg.
  double kb (1.3806488e-23); // boltzmann's constant (J/K)
  double mhDOS = get_DOS_EffectiveMassP(material);
  double dnbnd0 = 2.0*M_PI*mhDOS*e_mass*kb*temp/(h_planck*h_planck);
  double Nv = 2.0*pow(dnbnd0,1.5)/1.0e6;

  return Nv;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNi
// Purpose       : returns intrinsic electron concentration.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/12
// ----------------------------------------------------------------------------
double MaterialSupport::getNi (const std::string & material, double temp)
{
  ExtendedString mater = material;
  mater.toLower();
  double ni=0.0;
  double kbq = 8.6173324e-5; // boltzmann's constant  (eV/K)
  double bg = bandgap(mater,temp);
  double Nc = getNc(mater,temp);
  double Nv = getNv(mater,temp);
  ni = sqrt (Nc * Nv) * exp(-bg/(2.0 * kbq * temp));
  return ni;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getNi_old
// Purpose       : returns intrinsic electron concentration.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getNi_old (const std::string & material, double temp)
{
  ExtendedString mater = material;
  mater.toLower();
  double ni=0.0;

  if (mater == "si")
  {
    ni = 4.9e15
         * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
         * pow(6.0,0.5) * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
    // ni = 1.25e10;
  }
  else if (mater == "gaas")
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater == "ge")
  {
    ni = 4.9e15
         * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
         * 2.0 * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
   // ni = 2.5e13;
  }
  // for the next several, as they are all III-V materials, I copied the
  // gaas functions.  I *think* this is correct, as I *think* that Mc is
  // going to be 1.0 for all of these.
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else if (mater == "inp")
  {
    ni = 4.9e15
          * pow(getEffectiveMassN(mater)*getEffectiveMassP(mater),0.75)
          * pow(temp, 1.5) * exp(-bandgap(mater,temp)/
	                                   (2.0 * 8.6174e-5 * temp));
  }
  else
  {
    Report::UserError() <<  "MaterialSupport::getNi:  "
                        << material
                        <<  " material not recognized in getNi_old.";
  }

  return ni;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::getRelPerm
// Purpose       : returns relative permitivity
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::getRelPerm (const std::string & material)
{
  ExtendedString mater = material;
  mater.toLower();

  double perm;
  if (mater == "si")
  {
    perm = 11.8;
  }
  else if (mater == "sio2")
  {
    perm = 3.9;
  }
  else if (mater == "ge" )
  {
    perm = 16.0;
  }
  else if (mater == "gaas")
  {
    perm = 13.2;
  }
  else if (mater=="inalas" || mater=="alinas") // indium aluminum arsenide
  {
    perm = 12.5;
  }
  else if (mater=="ingaas" || mater=="gainas") // indium galium arsenide
  {
    perm = 14.0;
  }
  else if (mater == "gan")
  {
    perm = 11.9;
  }
  else if (mater == "inp")
  {
    perm = 12.6;
  }
  else
  {
    Report::UserFatal() << material << " material not recognized getRelPerm.";
  }

  return perm;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcRsrh
// Purpose       : Calculates schockley-read-hall recombination.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcRsrh
  (const std::string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pn = Ni*Ni;

  double A = (n*p-pn);
  double B = (tp*(n+Ni)+tn*(p+Ni));

  double arg = CONSTMAX_EXP_ARG;
  if (B >= exp(arg)) B = exp(arg);

  return (A/B);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRsrhN
// Purpose       : Calculates partial derivatives for schockley-read-hall
//                 recombination, with respect to electron density.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRsrhN
  (const std::string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pdRsrhN;
  double A1, B1, C1;
  double dAdn;
  double dBdn;

  double pn = Ni*Ni;

  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdn = (p);

  C1 = (tp*(n+Ni)+tn*(p+Ni));
  if (C1 >= exp(arg)) C1 = exp(arg);

  B1 = 1.0/C1;
  dBdn = -1.0/(C1*C1) * tp;

  pdRsrhN = dAdn * B1 + dBdn * A1;

  return pdRsrhN;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRsrhP
// Purpose       : Calculates partial derivatives for schockley-read-hall
//                 recombination, with respect to hole density.
//
// Special Notes : For this function, it shouldn't matter if the variables
//                 are scaled or not.
//
//                 The material dependence here comes indirectly, from the
//                 lifetimes, the carrier densities, and Ni, the intrinsic
//                 concentration.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRsrhP
  (const std::string & material, double ni, double n, double p, double tn, double tp)
{
  double Ni = ni;
  double pdRsrhP;
  double A1, B1, C1;
  double dAdp;
  double dBdp;

  double pn = Ni*Ni;

  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdp = (n);

  C1 = (tp*(n+Ni)+tn*(p+Ni));
  if (C1 >= exp(arg)) C1 = exp(arg);

  B1 = 1.0/C1;
  dBdp = -1.0/(C1*C1) * tn;

  pdRsrhP = dAdp * B1 + dBdp * A1;

  return pdRsrhP;
}


// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcRaug
// Purpose       : Calculates Auger recombination.
//
// Special Notes : This function MUST be called with unscaled 
//                 variables.  The coefficients Cn, Cp are in units of
//                 cm^6/s.  But they cannot be scaled correctly in 
//                 this function without knowing C0 and t0.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcRaug
  (const std::string & material, double ni, double n, double p)
{
  double Ni = ni;
  double Cn = MaterialLayer(material).augpnn;
  double Cp = MaterialLayer(material).augnpp;

  double pn = Ni*Ni;

  double A = (n*p-pn);
  double C = (Cn*n+Cp*p);

  double arg = CONSTMAX_EXP_ARG;
  if (C >= exp(arg)) C = exp(arg);

  return (A*C);
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRaugN
// Purpose       : Calculates partial derivative w.r.t. electron density
//                 for Auger recombination.
//
// Special Notes : This function MUST be called with unscaled 
//                 variables.  The coefficients Cn, Cp are in units of
//                 cm^6/s.  But they cannot be scaled correctly in 
//                 this function without knowing C0 and t0.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRaugN
  (const std::string & material, double ni, double n, double p)
{
  double Ni = ni;
  double pdRaugN;
  double A1, B1;
  double dAdn;
  double dBdn;

  double Cn = MaterialLayer(material).augpnn;
  double Cp = MaterialLayer(material).augnpp;
  double pn = Ni*Ni;
  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdn = (p);

  B1 = (Cn*n+Cp*p);
  if (B1 >= exp(arg)) B1 = exp(arg);

  dBdn = Cn;

  pdRaugN = dAdn*B1 + A1*dBdn;

  return pdRaugN;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::pdRaugP
// Purpose       : Calculates partial derivative w.r.t. hole density
//                 for Auger recombination.
//
// Special Notes : This function MUST be called with unscaled 
//                 variables.  The coefficients Cn, Cp are in units of
//                 cm^6/s.  But they cannot be scaled correctly in 
//                 this function without knowing C0 and t0.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/03
// ----------------------------------------------------------------------------
double MaterialSupport::pdRaugP
  (const std::string & material, double ni, double n, double p)
{
  double Ni = ni;
  double pdRaugP;
  double A1, B1;
  double dAdp;
  double dBdp;

  double Cn = MaterialLayer(material).augpnn;
  double Cp = MaterialLayer(material).augnpp;
  double pn = Ni*Ni;
  double arg = CONSTMAX_EXP_ARG;

  A1 = (n*p-pn);
  if (A1 >= exp(arg)) A1 = exp(arg);

  dAdp = (n);

  B1 = (Cn*n+Cp*p);
  if (B1 >= exp(arg)) B1 = exp(arg);

  dBdp = Cp;

  pdRaugP = dAdp*B1 + A1*dBdp;

  return pdRaugP;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::workfunc
// Purpose       : This function returns the workfunction
//                 of various metals
//
// Special Notes :
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/15/03
//----------------------------------------------------------------------------
double MaterialSupport::workfunc(std::string & metal)
{
  double wkfunc=0.0;

  ExtendedString metalName = metal;
  metalName.toLower ();

  if (metalName=="al")
  {
    wkfunc = 4.10;   //aluminum
  }
  else if (metalName=="ppoly")
  {
    wkfunc = 5.25;   // p+-polysilicon
  }
  else if (metalName=="npoly")
  {
    wkfunc = 4.17;   // n+-polysilicon
  }
  else if (metalName=="mo")
  {
    wkfunc = 4.53;  // molybdenum
  }
  else if (metalName=="w")
  {
    wkfunc = 4.63;  // tungsten
  }
  else if (metalName=="modi")
  {
    wkfunc = 4.80;  // molybdenum disilicide
  }
  else if (metalName=="wdi")
  {
    wkfunc = 4.80;  // tungsten disilicide
  }
  else if (metalName=="cu")
  {
    wkfunc = 4.25;   // copper
  }
  else if (metalName=="pt")
  {
    wkfunc = 5.30;   // platinum
  }
  else if (metalName=="au")
  {
    wkfunc = 4.80;   // gold
  }
  else if (metalName=="neutral")
  {
    wkfunc = 0.0;
  }
  else
  {
    Report::UserFatal() << metalName << " material not recognized.";
  }

  return wkfunc;
}
//----------------------------------------------------------------------------
// Function      : MaterialSupport::affin
// Purpose       : This function returns the electron affinity
//                 of various semiconductor materials
//
// Special Notes :
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/15/03
//---------------------------------------------------------------------------
double MaterialSupport::affin(const std::string & material)
{

  double afty=0.0;

  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")
  {
    afty = 4.17;      // silicon
  }
  else if (materialName=="ge")
  {
    afty = 4.00;     // germanium
  }
  else if (materialName=="gaas")
  {
    afty = 4.07;    // gallium arsenide
  }
  else if (materialName=="sio2")
  {
    afty = 0.97;    // silicon dioxide
  }
  else if (materialName=="nitride")
  {
    afty = 0.97;    // silicon nitride
  }
  else if (materialName=="sapphire")
  {
    afty = 0.97;     // sapphire (also known as aluminum oxide)
  }
  else
  {
    Report::UserError0() << materialName << " material not recognized.";
  }

  return afty;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::bandgap
// Purpose       : This function returns the electronic bandgap
//                 of various semiconductor materials.
//
// Special Notes : Reference for temperature-dependent semiconductor
//                 materials is "The Standard Thermodynamic Function
//                 of the Formation of Electrons and Holes in Ge, Si,
//                 GaAs, and GaP," by C. D. Thurmond, J. Electrochem. Soc.,
//                 vol. 122, p. 1133, 1975.
//
// Scope         : public
// Creator       : Deborah Fixel, SNL, Parallel Computational Sciences
// Creation Date : 7/18/03
//---------------------------------------------------------------------------
double MaterialSupport::bandgap(const std::string & material,  double temp)
{
  double gap = 0.0;
  ExtendedString materialName = material;
  materialName.toLower();

  if (materialName=="si")   // silicon
  {
    gap = 1.17 - 4.73e-4*pow(temp,2.0)/(temp + 636.0);
  }
  else if (materialName=="ge")  // germanium
  {
    gap = 0.7437 - 4.774e-4*pow(temp,2.0)/(temp + 235);
  }
  else if (materialName=="gaas")  // gallium arsenide
  {
    gap = 1.519 - 5.405e-4*pow(temp,2.0)/(temp + 204);
  }
  else if (materialName=="ingap")  
  {
    gap = 1.86098;
  }
  else if (materialName=="sio2")  // silicon dioxide
  {
    gap = 9.00;
  }
  else if (materialName=="nitride")  // silicon nitride
  {
    gap = 4.7;
  }
  else if (materialName=="sapphire")  // sapphire
  {
    gap = 4.7;
  }
  else if (materialName=="inalas" || materialName=="alinas") // indium aluminum arsenide
  {
    gap = 1.46;
  }
  else if (materialName=="ingaas" || materialName=="gainas") // indium galium arsenide
  {
    gap = 0.75;
  }
  else if (materialName=="gan")  // sapphire
  {
    gap = 3.4;
  }
  else if (materialName=="inp")  // indium phosphide
  {
    gap = 1.07;
  }
  else
  {
    Report::UserError0() << materialName << " material not recognized in bandgap.";
  }

  return gap;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::Ebgn
// Purpose       : Band-gap narrowing
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/18/2014
//---------------------------------------------------------------------------
double MaterialSupport::Ebgn (
    const std::string & material,  
    const std::string & bgnModel,
    double dope,
    bool ntype)
{
  double Ebgn=0.0;

  if (bgnModel=="slotboom") // parameters taken from Medici manual
  {
    Ebgn = slotboomEbgn( material,  dope, ntype) ;
  }
  else if (bgnModel=="bennet-wilson") // default model for Sentaurus.
  {
    Ebgn = bennetWilsonEbgn ( material,  dope, ntype) ;
  }
  else if (bgnModel=="jain") // consistent XPD code.
  {
    Ebgn = jainEbgn( material,  dope, ntype) ;
  }
  else if (bgnModel=="jain2") // consistent XPD code.
  {
    Ebgn = jain2Ebgn( material,  dope, ntype) ;
  }
  else if (bgnModel=="jain3") // consistent XPD code.
  {
    Ebgn = jain3Ebgn( material,  dope, ntype) ;
  }
  else // assume this is the default
  {
    Ebgn = jainEbgn( material,  dope, ntype) ;
  }

  return Ebgn;
}


//----------------------------------------------------------------------------
// Function      : MaterialSupport::bennetWilsonEbgn
// Purpose       : Band-gap narrowing, using the Bennet-Wilson model.
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/2/2014
//---------------------------------------------------------------------------
double MaterialSupport::bennetWilsonEbgn ( const std::string & material,  double dope, bool ntype)
{
  double Ebgn=0.0;
  if (material=="si") 
  {
    double Eref = 6.84e-3;
    double Nref = 3.162e+18;
    if (dope >= Nref) // Bennet-Wilson model
    {
      Ebgn = Eref*std::pow(log(dope/Nref),2.0);
    }
  }
  else
  {
    // bennet-wilson not implemented for this material
    Ebgn=0.0;
  }
  return Ebgn;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::slotboomEbgn
// Purpose       : Band-gap narrowing, using the Slotboom model.
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/2/2014
//---------------------------------------------------------------------------
double MaterialSupport::slotboomEbgn ( const std::string & material,  double dope, bool ntype)
{
  double Ebgn = 0.0;
  double n0_bgn=0.0;
  double v0_bgn=0.0;
  double con_bgn=0.0;

  if (material=="si")   // silicon
  {
    // Slotboom model constants:
    n0_bgn = 1.0e17; // Nref
    v0_bgn = 9.0e-3; // Eref
    con_bgn = 0.5;
  }
  else if (material=="ge")  // germanium
  {
    n0_bgn = 1.0e17;
    v0_bgn = 9.0e-3;
    con_bgn = 0.5;
  }
  else if (material=="gaas")  // gallium arsenide
  {
    n0_bgn = 1.0e17;
    v0_bgn = 0.0;
    con_bgn = 0.0;
  }
  else if (material=="ingap")
  {
    n0_bgn = 1.0e17;
    v0_bgn = 0.0;
    con_bgn = 0.0;
  }
  else if (material=="sio2")  // silicon dioxide
  {
    n0_bgn = 1.0e17; // Nref
    v0_bgn = 9.0e-3; // Eref
    con_bgn = 0.5;
  }
  else
  {
    Report::UserError() << material << " material not implemented for the Slotboom band-gap narrowing model";
    return Ebgn;
  }

  double tmp1 = log(dope/n0_bgn);
  Ebgn = (v0_bgn) * (tmp1 + sqrt(tmp1*tmp1 + con_bgn));

  return Ebgn;
}

//----------------------------------------------------------------------------
// Function      : MaterialSupport::jainEbgn
//
// Purpose       : Band-gap narrowing, using the Jain model.
//
//                 This corresponds to the BGN2 model in Davinci/Medici
//
// Special Notes : 
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/18/2014
//---------------------------------------------------------------------------
double MaterialSupport::jainEbgn ( const std::string & material,  double dope, bool ntype)
{
  double Ebgn=0.0;

  double deltaEcn=0.0;
  double deltaEvn=0.0;
  double deltaEcp=0.0;
  double deltaEvp=0.0;

  double anc_bgn = 0.0;
  double bnc_bgn = 0.0;
  double cnc_bgn = 0.0;
  double anv_bgn = 0.0;
  double bnv_bgn = 0.0;
  double cnv_bgn = 0.0;
  double apc_bgn = 0.0;
  double bpc_bgn = 0.0;
  double cpc_bgn = 0.0;
  double apv_bgn = 0.0;
  double bpv_bgn = 0.0;
  double cpv_bgn = 0.0;

  if (material=="si")   // silicon
  {  
    anc_bgn = -14.84e-3;
    bnc_bgn = 0.0;
    cnc_bgn = 0.78e-3;

    anv_bgn = 0.0;
    bnv_bgn = 15.08e-3;
    cnv_bgn = 0.74e-3;

    apc_bgn = 0.0;
    bpc_bgn = -16.27e-3;
    cpc_bgn = -0.18e-3;

    apv_bgn = 18.46e-3;
    bpv_bgn = 0.0;
    cpv_bgn = -2.63e-3;
  }
  else if (material=="ge")  // germanium
  {
    anc_bgn = -8.67e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -2.02e-3;
    anv_bgn = 0.0;
    bnv_bgn = 8.14e-3;
    cnv_bgn = 2.29e-3;
    apc_bgn = -8.21e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -2.19e-3;
    apv_bgn = 0.0;
    bpv_bgn = 9.18e-3;
    cpv_bgn = 3.58e-3; 
  }
  else if (material=="gaas")  // gallium arsenide
  {
    anc_bgn = -16.30e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (material=="ingap")
  { 
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (material=="sio2")  // silicon dioxide (use same as Si)
  {
    anc_bgn = -14.84e-3;
    bnc_bgn = 0.0;
    cnc_bgn = 0.78e-3;

    anv_bgn = 0.0;
    bnv_bgn = 15.08e-3;
    cnv_bgn = 0.74e-3;

    apc_bgn = 0.0;
    bpc_bgn = -16.27e-3;
    cpc_bgn = -0.18e-3;

    apv_bgn = 18.46e-3;
    bpv_bgn = 0.0;
    cpv_bgn = -2.63e-3; 
  }
  else if (material=="inalas" || material=="alinas") // indium aluminum arsenide
  { 
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (material=="ingaas" || material=="gainas") // indium galium arsenide
  { 
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;

    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else if (material=="inp")  // indium phosphide
  { 
    anc_bgn = -16.3e-3;
    bnc_bgn = 0.0;
    cnc_bgn = -18.13e-3;
    anv_bgn = 0.0;
    bnv_bgn = 7.47e-3;
    cnv_bgn = 72.52e-3;
    apc_bgn = -9.71e-3;
    bpc_bgn = 0.0;
    cpc_bgn = -0.47e-3;
    apv_bgn = 0.0;
    bpv_bgn = 12.19e-3;
    cpv_bgn = 3.41e-3;
  }
  else
  {
    Report::UserError() << material << " material not implemented for the Jain band-gap narrowing model";
    return fabs(Ebgn);
  }

  double dopeScaled=fabs(dope)/1.0e+18;

  if (ntype)
  {
    deltaEcn=
      +anc_bgn*std::pow(dopeScaled,(1.0/3.0))
      +bnc_bgn*std::pow(dopeScaled,(0.25))
      +cnc_bgn*std::pow(dopeScaled,(0.5));

    deltaEvn=
      +anv_bgn*std::pow(dopeScaled,(1.0/3.0))
      +bnv_bgn*std::pow(dopeScaled,(0.25))
      +cnv_bgn*std::pow(dopeScaled,(0.5));

    Ebgn = deltaEcn-deltaEvn;
  }
  else // p-type
  {
    deltaEcp=
      +apc_bgn*std::pow(dopeScaled,(1.0/3.0))
      +bpc_bgn*std::pow(dopeScaled,(0.25))
      +cpc_bgn*std::pow(dopeScaled,(0.5));

    deltaEvp=
      +apv_bgn*std::pow(dopeScaled,(1.0/3.0))
      +bpv_bgn*std::pow(dopeScaled,(0.25))
      +cpv_bgn*std::pow(dopeScaled,(0.5));

    Ebgn = deltaEcp-deltaEvp;
  }

  return fabs(Ebgn);
}


//----------------------------------------------------------------------------
// Function      : MaterialSupport::jain2Ebgn
// Purpose       : Band-gap narrowing, using the Jain model, another version.
//
//                 This one is from the original Jain paper, rather than
//                 from the DaVinci manual.  The paper doesn't have the same 
//                 materials as the Davinci manual, so the list of supported 
//                 materials is not the same
//
// Special Notes : 
// 
// Jain's paper reference:
// S. C. Jain, J. M. McGregor, D. J. Roulston, "Band‐gap narrowing in novel III‐V semiconductors"
// 1990, Journal of Applied Physics, p3747-3749, vol. 68, number 7, doi:10.1063/1.346291
// https://aip.scitation.org/doi/abs/10.1063/1.346291
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/18/2014
//---------------------------------------------------------------------------
double MaterialSupport::jain2Ebgn ( const std::string & material,  double dope, bool ntype)
{
  double deltaEbgn=0.0;

  double A_p=0.0;
  double B_p=0.0;
  double C_p=0.0;

  double A_n=0.0;
  double B_n=0.0;
  double C_n=0.0;

  if (material=="gaas")  // gallium arsenide
  {
    A_p=9.83e-9;
    B_p=3.90e-7;
    C_p=3.90e-12;

    A_n=16.5e-9;
    B_n=2.39e-7;
    C_n=91.4e-12;
  }
  else if (material=="insb")
  { 
    A_p=7.28e-9;
    B_p=2.58e-7;
    C_p=3.30e-12;

    A_n=12.2e-9;
    B_n=1.09e-7;
    C_n=604.0e-12;
  }
  else if (material=="inas")
  {
    A_p=8.34e-9;
    B_p=2.91e-7;
    C_p=4.53e-12;

    A_n=14.0e-9;
    B_n=1.97e-7;
    C_n=57.9e-12;
  }
  else if (material=="ingaas" || material=="gainas") // indium galium arsenide
  { 
    A_p=9.20e-9;
    B_p=3.57e-7;
    C_p=3.65e-12;

    A_n=15.5e-9;
    B_n=1.95e-7;
    C_n=159.0e-12;
  }
  else if (material=="gasb")
  {    
    A_p=8.07e-9;
    B_p=2.80e-7;
    C_p=4.12e-12;

    A_n=13.6e-9;
    B_n=1.66e-7;
    C_n=119.0e-12;
  }
  else if (material=="inp")  // indium phosphide
  {     
    A_p=10.3e-9;
    B_p=4.43e-7;
    C_p=3.38e-12;

    A_n=17.2e-9;
    B_n=2.62e-7;
    C_n=98.4e-12;
  }
  else if (material=="alsb")
  {
    A_p=11.5e-9;
    B_p=5.30e-7;
    C_p=3.53e-12;

    A_n=10.1e-9;
    B_n=3.09e-7;
    C_n=8.27e-7;
  }
  else if (material=="alas")
  {
    A_p=10.6e-9;
    B_p=5.47e-7;
    C_p=3.01e-12;

    A_n=9.76e-9;
    B_n=4.33e-7;
    C_n=2.93e-7;
  }
  else if (material=="gap")
  {
    A_p=12.7e-9;
    B_p=5.85e-7;
    C_p=3.90e-12;

    A_n=10.7e-9;
    B_n=3.45e-7;
    C_n=9.97e-7;
  }
  else
  {
    Report::UserError() << material << " material not implemented for jain2Ebgn function"
     <<" (which implements the Jain bandgap-narrowing model using parameters from the original paper)";
    return deltaEbgn;
  }

  double dopeScaled=fabs(dope);

  if (ntype) // n-type
  {
    deltaEbgn=
      +A_n*std::pow(dopeScaled,(1.0/3.0))
      +B_n*std::pow(dopeScaled,(0.25))
      +C_n*std::pow(dopeScaled,(0.5));
  }
  else // p-type
  {
    deltaEbgn=
      +A_p*std::pow(dopeScaled,(1.0/3.0))
      +B_p*std::pow(dopeScaled,(0.25))
      +C_p*std::pow(dopeScaled,(0.5));
  }

  return deltaEbgn;
}


//----------------------------------------------------------------------------
// Function      : MaterialSupport::jain3Ebgn
// Purpose       : Band-gap narrowing, using the Jain model, another version
// Special Notes : Ntype vs. Ptype is determined by the sign of "dope".  I 
//                 have only had time to set up a few materials.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/18/2014
//---------------------------------------------------------------------------
double MaterialSupport::jain3Ebgn ( const std::string & material,  double dope, bool ntype)
{
  double deltaEbgn=0.0;
  double A_p=0.0;
  double A_n=0.0;

  if (material=="gaas")  // gallium arsenide
  {
    A_p=2.6e-8;
    //A_n=16.5e-9;
  }
  else if (material=="ingaas" || material=="gainas") // indium galium arsenide
  { 
    A_p=2.43e-8;
    //A_n=15.5e-9;
  }
  else if (material=="gasb")  // indium phosphide
  { 
    A_p=2.22e-8;
  }
  else
  {
    Report::UserError() << material << " material not implemented for the jain3Ebn"
     <<" version of the Jain band-gap narrowing model";
    return deltaEbgn;
  }

  double dopeScaled=fabs(dope);

  if (ntype) // n-type
  {
    deltaEbgn=A_n*std::pow(dopeScaled,(1.0/3.0));
  }
  else // p-type
  {
    deltaEbgn=A_p*std::pow(dopeScaled,(1.0/3.0));
  }

  return deltaEbgn;
}

// ----------------------------------------------------------------------------
// Function      : MaterialSupport::calcLt
// Purpose       : This function calculates carrier lifetimes.
//
// Special Notes : holeFlag parameter indicates electrons or holes:
//                     holeFlag = true  -> holes
//                     holeFlag = false -> electrons
//
// This function assumes that conc is an absolute value.  
//
// This function comes from this paper:
//
//     "Analysis of High-Efficiency Silicon Solar Cells",
//     IEEE Transactions on Electron Devices, by Harry T.
//     Weaver and R. D. Nasby, vol. ED-28, no. 5, May 1981.
//
// In the paper, conc = total doping, i.e. (N_A + N_D)
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/20/03
// ----------------------------------------------------------------------------
double MaterialSupport::calcLt (bool holeFlag, double conc, std::string  material)
{
  double lt = 0.0;
  double LT0, Nref;

  conc = fabs(conc);

  if (material=="si")
  {
    if (holeFlag)
    {
      LT0   = 3.52e-5;
      Nref  = 7.1e15;
      lt = LT0 / (1.0 + conc / Nref);
    }
    else
    {
      LT0   = 3.95e-4;
      Nref  = 7.1e15;
      lt = LT0 / (1.0 + conc / Nref);
    }
  }
  else if (material=="gaas")
  {
    if (holeFlag)
    // Adding in SILVACO GaAs defaults here
    {
      LT0   = 2.0e-8;
      Nref  = 5.0e16;
      lt = LT0 / (1.0 + conc / Nref);
    }
    else
    // Adding in SILVACO GaAs defaults here
    {
      LT0   = 1.0e-9;
      Nref  = 5.0e16;
      lt = LT0 / (1.0 + conc / Nref);
    }
  }
  else
  {
    Report::UserFatal() << material << " material not implemented for carrier lifetime model.";
  }

  return lt;
}

} // namespace Device
} // namespace Xyce

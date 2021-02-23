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

//-------------------------------------------------------------------------
//
// Purpose        : This file contains similar functions to the spice3f5
//                  file, devsup.c.  It contains support routines for
//                  device models
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <algorithm>

// ----------   Xyce Includes   ----------
#include <N_UTL_Math.h>
#include <N_DEV_DeviceSupport.h>
#include <N_DEV_Const.h>
#include <N_UTL_RandomNumbers.h>

namespace Xyce {
namespace Device {

// For block gainscale homotopy
static Xyce::Util::RandomNumbers *theRandomGenerator=0;

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::lambertw
// Purpose       : provides a lambert-w function for diodes and BJT's.
// Special Notes :
//
// Purpose.  Evaluate principal branch of Lambert W function at x.
//
// w = w(x) is the value of Lambert's function.
// ierr = 0 indicates a safe return.
// ierr = 1 if x is not in the domain.
// ierr = 2 if the computer arithmetic contains a bug.
// xi may be disregarded (it is the error).
//
// Prototype: void lambertw( double, double, int, double);
//
// Reference:
// T.C. Banwell
// Bipolar transistor circuit analysis using the Lambert W-function,
// IEEE Transactions on Circuits and Systems I: Fundamental Theory
// and Applications
//
// vol. 47, pp. 1621-1633, Nov. 2000.
//
// Scope         : public
// Creator       : David Day,  SNL
// Creation Date : 04/16/02
//-----------------------------------------------------------------------------
void DeviceSupport::lambertw(double x, double &w, int &ierr, double &xi)
{
  int i=0, maxit = 10;
  const double turnpt = -exp(-1.), c1 = 1.5, c2 = .75;
  double r, r2, r3, s, mach_eps, relerr = 1., diff;
  mach_eps = 2.e-15;   // float:2e-7
  ierr = 0;

  if( x > c1)
  {
    w = c2*log(x);
    xi = log( x/ w) - w;
  }
  else
  {
    if( x >= 0.0)
    {
      w = x;
      if( x == 0. ) return;
      if( x < (1-c2) ) w = x*(1.-x + c1*x*x);
      xi = - w;
    }
    else
    {
      if( x >= turnpt)
      {
        if( x > -0.2 )
        {
          w = x*(1.0-x + c1*x*x);
          xi = log(1.0-x + c1*x*x) - w;
        }
        else
        {
          diff = x-turnpt;
          if( diff < 0.0 ) diff = -diff;
          w = -1 + sqrt(2.0*exp(1.))*sqrt(x-turnpt);
          if( diff == 0.0 ) return;
          xi = log( x/ w) - w;
        }
      }
      else
      {
        ierr = 1; // x is not in the domain.
        w = -1.0;
        return;
      }
    }
  }

  while( relerr > mach_eps  && i<maxit)
  {
     r = xi/(w+1.0);   //singularity at w=-1
     r2 = r*r;
     r3 = r2*r;
     s  = 6.*(w+1.0)*(w+1.0);
     w = w * (  1.0 + r + r2/(2.0*( w+1.0)) - (2. * w -1.0)*r3/s  );
     if( w * x < 0.0 ) w = -w;
     xi = log( x/ w) - w;

     if( x>1.0 )
     {
       relerr =  xi / w;
     }
     else
     {
       relerr =  xi;
     }
     if(relerr < 0.0 ) relerr = -relerr;
     ++i;
   }
   if( i == maxit ) ierr = 2;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::limvds
// Purpose       : limit the per-iteration change of VDS
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
double DeviceSupport::limvds ( double vnew, double vold)
{

  if(vold >= 3.5)
  {
    if(vnew > vold) vnew = std::min(vnew,(3.0 * vold) + 2.0);
    else
    {
      if (vnew < 3.5) vnew = std::max(vnew,2.0);
    }
  }
  else
  {
    if(vnew > vold) vnew = std::min(vnew, 4.0);
    else            vnew = std::max(vnew,-0.5);
  }
  return(vnew);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::pnjlim
// Purpose       : limit the per-iteration change of PN junction voltages
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
double DeviceSupport::pnjlim (
    double vnew,
    double vold,
    double vt,
    double vcrit,
    int *icheck
  )
{
  double arg;

  if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt)))
  {
    if(vold > 0)
    {
      arg = 1 + (vnew - vold) / vt;

      if(arg > 0)
      {
        vnew = vold + vt * log(arg);
      }
      else
      {
        vnew = vcrit;
      }
    }
    else
    {
      vnew = vt *log(vnew/vt);
    }

    *icheck = 1;
  }
  else
  {
    *icheck = 0;
  }

  return(vnew);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::pnjlim_new
// Purpose       : limit the per-iteration change of PN junction voltages
// Special Notes : Copied from NGSpice, which has the following comment:
//     This code has been fixed by Alan Gillespie adding limiting
//     for negative voltages
// Scope         : public
// Creator       : Tom Russo, SNL 1445, Electrical Systems Modeling
// Creation Date : 11/19/2012
//-----------------------------------------------------------------------------
double DeviceSupport::pnjlim_new (
    double vnew,
    double vold,
    double vt,
    double vcrit,
    int *icheck
  )
{
  double arg;

  if((vnew > vcrit) && (fabs(vnew - vold) > (vt + vt)))
  {
    if(vold > 0)
    {
      arg = (vnew - vold) / vt;

      if(arg > 0)
      {
        vnew = vold + vt * (2+log(arg-2));
      }
      else
      {
        vnew = vold - vt*(2+log(2-arg));
      }
    }
    else
    {
      vnew = vt *log(vnew/vt);
    }

    *icheck = 1;
  }
  else
  {
    if (vnew < 0)
    {
      if (vold > 0)
      {
        arg= -vold -1;
      }
      else
      {
        arg = 2*vold -1;
      }
      if (vnew < arg)
      {
        vnew = arg;
        *icheck=1;
      }
      else
      {
        *icheck = 0;
      }
    }
    else
    {
      *icheck = 0;
    }
  }

  return(vnew);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::fetlim
// Purpose       : limit the per-iteration change of FET voltages
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
double DeviceSupport::fetlim (
    double vnew,
    double vold,
    double vto
  )
{
  double vtsthi;
  double vtstlo;
  double vtox;
  double delv;
  double vtemp;

  vtsthi = fabs(2*(vold-vto))+2;
  vtstlo = vtsthi/2 +2;
  vtox = vto + 3.5;
  delv = vnew-vold;

  if (vold >= vto)
  {
    if(vold >= vtox)
    {
      if(delv <= 0)
      {
        // going off
        if(vnew >= vtox)
        {
          if(-delv >vtstlo) vnew =  vold - vtstlo;
        }
        else
        {
          vnew = std::max(vnew,vto+2.0);
        }
      }
      else
      {
        // staying on
        if(delv >= vtsthi) vnew = vold + vtsthi;
      }
    }
    else
    {
      // middle region
      if(delv <= 0)
      {
        // decreasing
        vnew = std::max(vnew,vto-0.5);
      }
      else
      {
        // increasing
        vnew = std::min(vnew,vto+4.0);
      }
    }
  }
  else
  {
    // off
    if(delv <= 0)
    {
      if(-delv >vtsthi) vnew = vold - vtsthi;
    }
    else
    {
      vtemp = vto + 0.5;
      if(vnew <= vtemp)
      {
        if(delv >vtstlo) vnew = vold + vtstlo;
      }
      else
      {
        vnew = vtemp;
      }
    }
  }
  return(vnew);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::cmeyer
// Purpose       : Compute the MOS overlap capacitances as functions of the
//                 device terminal voltages
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
void DeviceSupport::cmeyer (
   double vgs0,    // initial voltage gate-source
   double vgd0,    // initial voltage gate-drain
   double vgb0,    // initial voltage gate-bulk
   double von0,
   double vdsat0,
   double vgs1,    // final voltage gate-source
   double vgd1,    // final voltage gate-drain
   double vgb1,    // final voltage gate-bulk
   double covlgs,  // overlap capacitance gate-source
   double covlgd,  // overlap capacitance gate-drain
   double covlgb,  // overlap capacitance gate-bulk
   double *cgs,
   double *cgd,
   double *cgb,
   double phi,
   double cox,
   double von,
   double vdsat
 )
{
    double vdb;
    double vdbsat;
    double vddif;
    double vddif1;
    double vddif2;
    double vgbt;

    *cgs = 0;
    *cgd = 0;
    *cgb = 0;

    vgbt = vgs1-von;
    if (vgbt <= -phi)
    {
      *cgb = cox;
    }
    else if (vgbt <= -phi/2)
    {
      *cgb = -vgbt*cox/phi;
    }
    else if (vgbt <= 0)
    {
      *cgb = -vgbt*cox/phi;
      *cgs = cox/(7.5e-1*phi)*vgbt+cox/1.5;
    }
    else
    {
      vdbsat = vdsat-(vgs1-vgb1);
      vdb = vgb1-vgd1;
      if (vdbsat <= vdb)
      {
        *cgs = cox/1.5;
      }
      else
      {
        vddif = 2.0*vdbsat-vdb;
        vddif1 = vdbsat-vdb-1.0e-12;
        vddif2 = vddif*vddif;
        *cgd = cox*(1.0-vdbsat*vdbsat/vddif2)/1.5;
        *cgs = cox*(1.0-vddif1*vddif1/vddif2)/1.5;
      }
    }

    vgbt = vgs0-von0;
    if (vgbt <= -phi)
    {
      *cgb += cox;
    }
    else if (vgbt <= -phi/2)
    {
      *cgb += -vgbt*cox/phi;
    }
    else if (vgbt <= 0)
    {
      *cgb += -vgbt*cox/phi;
      *cgs += cox/(7.5e-1*phi)*vgbt+cox/1.5;
    }
    else
    {
      vdbsat = vdsat0-(vgs0-vgb0);
      vdb = vgb0-vgd0;
      if (vdbsat <= vdb)
      {
        *cgs += cox/1.5;
      }
      else
      {
        vddif = 2.0*vdbsat-vdb;
        vddif1 = vdbsat-vdb-1.0e-12;
        vddif2 = vddif*vddif;
        *cgd += cox*(1.0-vdbsat*vdbsat/vddif2)/1.5;
        *cgs += cox*(1.0-vddif1*vddif1/vddif2)/1.5;
      }
    }

    *cgs = *cgs *.5 + covlgs;
    *cgd = *cgd *.5 + covlgd;
    *cgb = *cgb *.5 + covlgb;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::qmeyer
// Purpose       : Compute the MOS overlap capacitances as functions of the
//                 device terminal voltages
//
// Special Notes : ARGSUSED  because vgb is no longer used
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/17/01
//-----------------------------------------------------------------------------
void DeviceSupport::qmeyer (
   double vgs,    // initial voltage gate-source
   double vgd,    // initial voltage gate-drain
   double vgb,    // initial voltage gate-bulk
   double von,
   double vdsat,
   double & capgs,  // non-constant portion of g-s overlap capacitance
   double & capgd,  // non-constant portion of g-d overlap capacitance
   double & capgb,  // non-constant portion of g-b overlap capacitance
   double phi,
   double cox     // oxide capactiance
 )
{
  double vds;
  double vddif;
  double vddif1;
  double vddif2;
  double vgst;

  //double vgdt;
  //double vdenom;
  //double vdenom2;


  vgst = vgs-von;
  if (vgst <= -phi)
  {
    capgb = cox/2;
    capgs = 0;
    capgd = 0;
  }
  else if (vgst <= -phi/2)
  {
    capgb = -vgst*cox/(2*phi);
    capgs = 0;
    capgd = 0;
  }
  else if (vgst <= 0)
  {
    capgb = -vgst*cox/(2*phi);
    capgs = vgst*cox/(1.5*phi)+cox/3;
    capgd = 0;
  }
  else
  {
    vds = vgs-vgd;
    if (vdsat <= vds)
    {
      capgs = cox/3;
      capgd = 0;
      capgb = 0;
    }
    else
    {

      vddif = 2.0*vdsat-vds;
      vddif1 = vdsat-vds/*-1.0e-12*/;
      vddif2 = vddif*vddif;
      capgd = cox*(1.0-vdsat*vdsat/vddif2)/3;
      capgs = cox*(1.0-vddif1*vddif1/vddif2)/3;
      capgb = 0;


      //vgdt = vgd-von;
      //vdenom=vgst + vgdt;
      //vdenom2=vdenom*vdenom;

      //capgd = cox*(1.0-vgst*vgst/vdenom2)/3.0;
      //capgs = cox*(1.0-vgdt*vgdt/vdenom2)/3.0;

    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::qmeyerderivs
// Purpose       : Computes the partial derivatives of the Meyer capacitances
//                 with respect to various voltages.
//
// Special Notes : We need this in order to make the Meyer model MPDE
//                 compatible.
//
// Scope         : public
// Creator       : Keith R. Santarelli, SNL, Electrical & Microsystems Modeling
// Creation Date : 02/08/08
//-----------------------------------------------------------------------------
void DeviceSupport::qmeyerderivs
    (
       double vgs,    // initial voltage gate-source
       double vgd,    // initial voltage gate-drain
       double vgb,    // initial voltage gate-bulk
       double von,
       double vdsat,
       double & dcapgsdvgs, //partial deriv. of capgs with respect to vgs
       double & dcapgsdvgb, //partial deriv. of capgs with respect to vgb
       double & dcapgsdvgd, //partial deriv. of capgs with respect to vgd
       double & dcapgddvgs, //partial deriv. of capgd with respect to vgs
       double & dcapgddvgb, //partial deriv. of capgd with respect to vgb
       double & dcapgddvgd, //partial deriv. of capgd with respect to vgd
       double & dcapgbdvgs, //partial deriv. of capgb with respect to vgs
       double & dcapgbdvgb, //partial deriv. of capgb with respect to vgb
       double & dcapgbdvgd, //partial deriv. of capgb with respect to vgd
       double phi,
       double cox,     // oxide capactiance
       int Dtype //transistor type
    )
{
  double vgst;
  double vds;
  double vdenom, vdenom3;
  double vgdt;

  vgst = vgs-von;
  if (vgst <= -phi)
  {
    dcapgsdvgs=0;
    dcapgsdvgb=0;
    dcapgsdvgd=0;
    dcapgddvgs=0;
    dcapgddvgb=0;
    dcapgddvgd=0;
    dcapgbdvgs=0;
    dcapgbdvgb=0;
    dcapgbdvgd=0;

  }
  else if (vgst <= -phi/2)
  {
    dcapgsdvgs=0;
    dcapgsdvgb=0;
    dcapgsdvgd=0;
    dcapgddvgs=0;
    dcapgddvgb=0;
    dcapgddvgd=0;
    dcapgbdvgs=-1.0*cox/(2.0*phi);
    dcapgbdvgb=0;
    dcapgbdvgd=0;
  }
  else if (vgst <= 0)
  {
    dcapgsdvgs=cox/(1.5*phi);
    dcapgsdvgb=0;
    dcapgsdvgd=0;
    dcapgddvgs=0;
    dcapgddvgb=0;
    dcapgddvgd=0;
    dcapgbdvgs=-1.0*cox/(2.0*phi);
    dcapgbdvgb=0;
    dcapgbdvgd=0;
  }
  else
  {
    vds = vgs-vgd;
    if (vdsat <= vds)
    {
      dcapgsdvgs=0;
      dcapgsdvgb=0;
      dcapgsdvgd=0;
      dcapgddvgs=0;
      dcapgddvgb=0;
      dcapgddvgd=0;
      dcapgbdvgs=0;
      dcapgbdvgb=0;
      dcapgbdvgd=0;
    }
    else
    {
      vgdt = vgd-von;
      vdenom=vgst + vgdt;
      vdenom3=vdenom*vdenom*vdenom;

      dcapgsdvgs=4.0/3.0*cox*vgdt*vgdt/vdenom3;
      dcapgsdvgb=0;
      dcapgsdvgd=-4.0/3.0*cox*vgst*vgdt/vdenom3;
      dcapgddvgs=-4.0/3.0*cox*vgst*vgdt/vdenom3;
      dcapgddvgb=0;
      dcapgddvgd=4.0/3.0*cox*vgst*vgst/vdenom3;
      dcapgbdvgs=0;
      dcapgbdvgb=0;
      dcapgbdvgd=0;

    }
  }

  //Now have to "type-ize" the cap. derivatives:

  //dcapgsdvgs=Dtype*dcapgsdvgs;
  //dcapgsdvgb=Dtype*dcapgsdvgb;
  //dcapgsdvgd=Dtype*dcapgsdvgd;
  //dcapgddvgs=Dtype*dcapgddvgs;
  //dcapgddvgb=Dtype*dcapgddvgb;
  //dcapgddvgd=Dtype*dcapgddvgd;
  //dcapgbdvgs=Dtype*dcapgbdvgs;
  //dcapgbdvgb=Dtype*dcapgbdvgb;
  //dcapgbdvgd=Dtype*dcapgbdvgd;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::noiseSupport
//
// Purpose       : Related to the NevalSrc function in spice, 
//                 excluding the multiplication by gain.
//
// Special Notes : The reason this function excludes the gain multiplication is
//                 that the devices in Xyce don't typically have access to the 
//                 complex version of the solution vector.  
//
//                 Noise sources are assumed to be on a branch connecting two
//                 circuit nodes.  The gain for that branch is computed as
//
//                 realVal = realSolution(node1) - realSolution(node2)
//                 imagVal = imagSolution(node1) - imagSolution(node2)
//                 gain = (realVal*realVal) + (imagVal*imagVal)
//
//                 AC and NOISE analysis in Xyce were implemented in a manner 
//                 such that the complex linear system is 100% constructed 
//                 in the relevant analysis classes in a real equivalent form.
//                 The devices don't know about it.
//
//                 Most of the Xyce code was written to support DC and transient,
//                 and the complex versions of objects such as the solution vector
//                 are not known in the device package.  Analysis types that
//                 rely on complex numbers were added much later, and we had a 
//                 desire to not hack the devices any more than necessary.
//
//                 So, the multiplication by gain happens up in the N_ANP_NOISE
//                 class, where it has enough information to do it.
//
//                 This is also the reason that this function doesn't bother 
//                 supporting the "GAIN" type, as it should just be returning 
//                 "noise=1" in that case.  However, up in the device model 
//                 it is important to understand what is going on in that case.  
//                 Any time you see a call in the original Spice noise model 
//                 to NevalSrc, but with the "GAIN" type, you should exclude 
//                 that call (or its equivalent, to this function) from the 
//                 Xyce version.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/17/2014
//-----------------------------------------------------------------------------
void DeviceSupport::noiseSupport (
  double & noise, double & lnNoise, 
  const int type, const double param, const double temp)
{
  switch (type) 
  {
    case SHOTNOISE: // param is the dc current in a semiconductor 
      noise = 2.0 * CONSTQ * fabs(param);          
      lnNoise = log( std::max(noise,N_MINLOG) ); 
      break;

    case THERMNOISE: // param is the conductance of a resistor 
      noise = 4.0 * CONSTboltz * temp * param;         
      lnNoise = log( std::max(noise,N_MINLOG) );
      break;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::contVds
// Purpose       : continuation adjustment for MOSFET drain-source voltage.
//
// Special Notes : min*vds is the value returned for vds when alpha=0.
//
//                 This idea is based, loosely, on a paper by Jaijeet
//                 Rosychowdhury.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/04/03
//-----------------------------------------------------------------------------
double DeviceSupport::contVds (double vds, double alpha, double min)
{
  if (min <= 0.0) min = 0.3;
  return ( vds * (alpha*(1.0-min) + min) );
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::contVgst
// Purpose       : continuation adjustment for MOSFET drain-source voltage.
//
// Special Notes : vgstConst is the value returned for vgst when alpha=0.
//
//                 This idea is based, loosely, on a paper by Jaijeet
//                 Rosychowdhury.
//
//                 The alpha=0 condition, for which vgs, or vgst is
//                 considered constant, essentially makes the device a
//                 single-state device.   No matter what actual voltage is
//                 applied across Vg-Vs, the device will act as though
//                 there is a fixed applied voltage.  If vgstconst is set
//                 to a high voltage (like the default, 3.0) then the
//                 devices are all in an "on" state for alpha=0.  If
//                 vgstconst is set to zero, then the devices are all in an
//                 off state for alpha=0.
//
//                 As the alpha parameter is swept from zero to 1, at some
//                 point the MOSFET will start behaving, functionally like
//                 a MOSFET.  At alpha=0, it doesn't - it kind of acts like
//                 a capacitor.  At what point in the sweep this
//                 functionality kicks in, depends on the point at which
//                 the Ids vs Vgs curve has an intersection with Ids=0, for
//                 reasonable values of Vgs.  (reasonable values probably
//                 being -1 to +5 volts)
//
//                 If Vgstconst is set to zero, this intersection happens
//                 immediately, on the first step where alpha is nonzero.
//                 If Vgstconst is fairly high (say 4.0 volts), this
//                 intersection happens late.  If it is too high, it never
//                 happens.
//
//                 Unfortunately, during the sweep, once it kicks in, it
//                 kicks in for every MOSFET almost all at once, and there
//                 is a sharp transition in the continuation curve.  The
//                 trick to making this homotopy work well, is to make this
//                 transition as gentle as possible, I think.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/04/03
//-----------------------------------------------------------------------------
double DeviceSupport::contVgst
  (double vgst, double alpha, double vgstConst)
{
  return ((alpha)*vgst + (1.0-alpha)*vgstConst);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::getGainScaleBlockID
// Purpose       : determines which block to associate the mosfet device
//                 This is used for a multi-block gainscale continuation.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 01/28/05
//-----------------------------------------------------------------------------
int DeviceSupport::getGainScaleBlockID(int numBlocks)
{
  if (theRandomGenerator==0)
    theRandomGenerator = new Xyce::Util::RandomNumbers(0,false);

  double val = theRandomGenerator->uniformRandom();

  // get rid of negatives
  val = val * val;
  val = sqrt(val);

  double interval = 1.0 / numBlocks;

  for (int i = 0; i < numBlocks; ++i) {
    double lowBound = ((double) i) * interval;
    double highBound = ((double) (i+1)) * interval;
    if ((val >= lowBound) && (val < highBound)) {
      return i;
    }
  }

  return (numBlocks-1);
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::getGainScaleBlockID
// Purpose       : computes random perturbation for stretched homotopy.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 01/28/05
//-----------------------------------------------------------------------------
double DeviceSupport::getRandomPerturbation()
{
  if (theRandomGenerator==0)
    theRandomGenerator = new Xyce::Util::RandomNumbers(0,false);

  double val = theRandomGenerator->uniformRandom();
  val = val * val;
  val = sqrt(val);
  return val;
}

//-----------------------------------------------------------------------------
// Function      : DeviceSupport::getGainScaleBlockID
// Purpose       : sets seed for random number generator used in getRandomPerturbation().
// Scope         : public
// Creator       : Richard Schie, Electrical Systems Modeling
// Creation Date : 10/01/12
//-----------------------------------------------------------------------------
int DeviceSupport::SetSeed(long seedIn)
{
  if (theRandomGenerator==0)
    theRandomGenerator = new Xyce::Util::RandomNumbers(seedIn,false);
  theRandomGenerator->seedRandom( seedIn, false );
  return 0;
}

} // namespace Device
} // namespace Xyce

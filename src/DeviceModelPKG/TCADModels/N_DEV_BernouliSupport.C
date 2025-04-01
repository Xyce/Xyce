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
// Creation Date  : 08/01/04
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_BernouliSupport.h>
#include <N_ERH_ErrorMgr.h>

// These are the old Bernouli breakpoint numbers.  They were generated
// using the SGF program brkpnts.c . They are appropriate for linux.
// The functions in this file are designed to generate these.  
// In practice, these seem to work fine.
#define BP0_BERN    -3.742945958815751e+01
#define BP1_BERN    -1.948848145621305e-02
#define BP2_BERN    1.230611609815494e-02
#define BP3_BERN    3.742945958815751e+01
#define BP4_BERN    7.451332191019408e+02
#define BP0_DBERN   -4.117119704160766e+01
#define BP1_DBERN   -3.742945958815751e+01
#define BP2_DBERN   -1.848271746976161e-02
#define BP3_DBERN   8.806697697210611e-03
#define BP4_DBERN   3.742945958815751e+01
#define BP5_DBERN   7.451332191019408e+02
#define BP0_AUX1    -8.301056680276218e-03
#define BP1_AUX1    8.301056680276218e-03
#define BP0_DAUX1   -4.826242066078996e-03
#define BP1_DAUX1   4.826242066078996e-03
#define BP0_AUX2    -4.436141955583643e+01
#define BP1_AUX2    3.680808162809191e+01
#define BP2_AUX2    7.451332191019419e+02
#define BP0_DAUX2   -7.451332191019419e+02
#define BP1_DAUX2   -4.436141955583643e+01
#define BP2_DAUX2   3.680808162809191e+01
#define BP3_DAUX2   7.451332191019419e+02
#define BP0_MISC    7.097827128183643e+02

#define PRECISION 1.0e-15
#define MAX_ITERATIONS 100

namespace Xyce {
namespace Device {

namespace {

// This code comes from the brkpnts.c program of the SGF framework

double Bbp0a(double x)  { return exp(x) - 1.0; }
double Bbp0b(double x)  { return - 1.0; }

double Bbp1a(double x)  { return x / (exp(x) - 1.0); }
double Bbp1b(double x)
{
  return 1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0));
}

double Bbp2a(double x)
{
  return 1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0));
}

double Bbp2b(double x)
{
  return x * exp(-x) / (1.0 - exp(-x));
}

double Bbp3a(double x)  { return 1.0 - exp(-x); }
double Bbp3b(double x)  { return 1.0; }

double Bbp4a(double x)  { return x * exp(-x); }
double Bbp4b(double x)  { return 0.0; }

double dBbp0a(double x) { return (1.0 - x) * exp(x) - 1.0; }
double dBbp0b(double x) { return -1.0; }

double dBbp2a(double x)
{
  return ((1.0 - x) * exp(x) - 1.0) / ((exp(x) - 1.0) * (exp(x) - 1.0));
}

double dBbp2b(double x) { return -0.5 + x/6.0 * (1.0 - x*x/30.0); }

double dBbp3a(double x) { return -0.5 + x/6.0 * (1.0 - x*x/30.0); }
double dBbp3b(double x)
{
  return (exp(-x)*(1.0 - x) - exp(-2.0 * x))/((1.0 - exp(-x))*(1.0 - exp(-x)));
}

double dBbp5a(double x) { return exp(-x) * (1.0 - x) - exp(-2.0 * x); }
double dBbp5b(double x) { return 0.0; }

double AUX1bp0a(double x) { return x / sinh(x); }
double AUX1bp0b(double x) { return 1.0 - x*x/6.0 * (1.0 - 7.0*x*x/60.0); }

double dAUX1bp0a(double x)
{
  return (sinh(x) - x*cosh(x)) / (sinh(x) * sinh(x));
}

double dAUX1bp0b(double x) { return -x/3.0 * (1.0 - 7*x*x/30.0); }

double AUX2bp0a(double x) { return 1.0; }
double AUX2bp0b(double x) { return 1.0 + exp(x); }

double AUX2bp1a(double x) { return 1.0 + exp(x); }
double AUX2bp1b(double x) { return exp(x); }

double AUX2bp2a(double x) { return exp(-x); }
double AUX2bp2b(double x) { return 0.0; }

double dAUX2bp0a(double x) { return exp(x); }
double dAUX2bp0b(double x) { return 0.0; }

} // namespace <empty>

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::BernouliSupport
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
BernouliSupport::BernouliSupport (bool regenerate):

  // aux1 function breakpoints:
  bp0_AUX1   ( BP0_AUX1 ),
  bp1_AUX1   ( BP1_AUX1 ),

  // aux1 derivative function breakpoints:
  bp0_DAUX1  ( BP0_DAUX1 ),
  bp1_DAUX1  ( BP1_DAUX1 ),

  // aux2 function breakpoints:
  bp0_AUX2   ( BP0_AUX2 ),
  bp1_AUX2   ( BP1_AUX2 ),
  bp2_AUX2   ( BP2_AUX2 ),

  // aux2 derivative function breakpoints:
  bp0_DAUX2  ( BP0_DAUX2 ),
  bp1_DAUX2  ( BP1_DAUX2 ),
  bp2_DAUX2  ( BP2_DAUX2 ),
  bp3_DAUX2  ( BP3_DAUX2 ),

  // This is the log of Xyce::Util::MachineDependentParams::DoubleMax()
  bp0_MISC   ( BP0_MISC )
{

  // generate constants internally, to override the hardired defaults.
  // years ago doing this caused problems on some platforms so they were 
  // left optional.
  if (regenerate)
  {
    // Aux1 breakpoints:
    bp0_AUX1 = Secant(AUX1bp0a, AUX1bp0b, -1.0e+00);
    bp1_AUX1 = Secant(AUX1bp0a, AUX1bp0b,  1.0e+00);

    // Aux1 derivative breakpoints:
    bp0_DAUX1 = Secant(dAUX1bp0a, dAUX1bp0b, -1.0e+00);
    bp1_DAUX1 = Secant(dAUX1bp0a, dAUX1bp0b,  1.0e+00);

    // Aux2 breakpoints:
    bp0_AUX2 = Asymptotic(AUX2bp0a, AUX2bp0b,  0.0e+00, -1.0e+02);
    bp1_AUX2 = Asymptotic(AUX2bp1a, AUX2bp1b,  0.0e+00,  1.0e+02);
    bp2_AUX2 = Asymptotic(AUX2bp2a, AUX2bp2b,  0.0e+00,  1.0e+02);

    // Aux2 derivative breakpoints:
    bp0_DAUX2 = Asymptotic(dAUX2bp0a, dAUX2bp0b,  0.0e+00, -1.0e+02);
    bp1_DAUX2 = bp0_AUX2;
    bp2_DAUX2 = bp1_AUX2;
    bp3_DAUX2 = bp2_AUX2;

    // miscellaneous
    bp0_MISC = log(Xyce::Util::MachineDependentParams::DoubleMax());
  }

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout().width(21); Xyce::dout().precision(13); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << std::endl;
    Xyce::dout() << "Bernouli function breakpoints: " <<std::endl;
    Xyce::dout() << "bp0_AUX1  = " << bp0_AUX1  << std::endl;
    Xyce::dout() << "bp1_AUX1  = " << bp1_AUX1  << std::endl;
    Xyce::dout() << "bp0_DAUX1 = " << bp0_DAUX1 << std::endl;
    Xyce::dout() << "bp1_DAUX1 = " << bp1_DAUX1 << std::endl;

    Xyce::dout() << "bp0_AUX2  = " << bp0_AUX2  << std::endl;
    Xyce::dout() << "bp1_AUX2  = " << bp1_AUX2  << std::endl;
    Xyce::dout() << "bp2_AUX2  = " << bp2_AUX2  << std::endl;
    Xyce::dout() << "bp0_DAUX2 = " << bp0_DAUX2 << std::endl;
    Xyce::dout() << "bp1_DAUX2 = " << bp1_DAUX2 << std::endl;
    Xyce::dout() << "bp2_DAUX2 = " << bp2_DAUX2 << std::endl;
    Xyce::dout() << "bp3_DAUX2 = " << bp3_DAUX2 << std::endl;
    Xyce::dout() << "bp0_MISC  = " << bp0_MISC << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::~BernouliSupport
// Purpose       : constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
BernouliSupport::~BernouliSupport ()
{

}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::BernouliSupport
// Purpose       : copy constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
BernouliSupport::BernouliSupport
  (const BernouliSupport & right)
{

}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::BernouliSupport
// Purpose       : copy constructor.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
int BernouliSupport::sign(double x)
{
  if      (x < 0.0) return(-1);
  else if (x > 0.0) return(+1);
  else              return(0);
}

//-----------------------------------------------------------------------------
// Function      : double BernouliSupport::Bisection
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
double BernouliSupport::Bisection
  (FUNC func1, FUNC func2, double Xpos, double Xneg)
{
  double Fpos = func1(Xpos) - func2(Xpos);
  double Fneg = func1(Xneg) - func2(Xneg);
  double Xmid, Fmid, Xlast;

  if    (Fpos == 0.0) return(Xpos);
  else if (Fneg == 0.0) return(Xneg);
  else if ((Fpos > 0.0) && (Fneg < 0.0)) ;
  else if ((Fpos < 0.0) && (Fneg > 0.0))
  {
    Xmid = Xpos;
    Xpos = Xneg;
    Xneg = Xmid;
  }
  else
  {
    Report::DevelFatal() << "BernouliSupport::Bisection: "
                         << " Initial interval may not contain a root";
  }

  Xlast = 0.0;
  do
  {
    Xmid = 0.5 * (Xpos + Xneg);
    Fmid = func1(Xmid) - func2(Xmid);
    if      (Fmid < 0.0)  Xneg = Xmid;
    else if (Fmid > 0.0)  Xpos = Xmid;
    if (Xlast == Xmid) return(Xmid);
    else Xlast = Xmid;
  } while (Xneg != Xpos);

  return(Xmid);

}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::Secant
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
double BernouliSupport::Secant(FUNC func1, FUNC func2, double x1)
{
  double slope, dx, x3, f3;
  int    s3, iteration;

  double x2 = 0.9 * x1;
  double f1 = func1(x1) - func2(x1);
  double f2 = func1(x2) - func2(x2);
  int    s2 = sign(x2);

  for(;;)
  {
    iteration = 0;
    slope = (f2 - f1) / (x2 - x1);
    dx = f2 / slope;
    x3 = x2 - dx;
    f3 = func1(x3) - func2(x3);
    s3 = sign(x3);

    while ((fabs(f3) >= fabs(f2)) || (s3 != s2))
    {
      dx /= 2.0;
      x3 += dx;
      f3  = func1(x3) - func2(x3);
      s3  = sign(x3);
      if (++iteration > MAX_ITERATIONS)
      {
        if (fabs(f2) <= 100.0 * PRECISION) return(x2);
        Report::DevelFatal() <<  "BernouliSupport::Secant: "
                             << " method not converging.";
      }
    }

    x1 = x2;
    x2 = x3;
    f1 = f2;
    f2 = f3;

    if ((fabs(dx / x2) <= PRECISION) || (fabs(f2) <= PRECISION)) break;
  }

  return(x3);
}

//-----------------------------------------------------------------------------
// Function      : BernouliSupport::Asymptotic
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/01/04
//-----------------------------------------------------------------------------
double BernouliSupport::Asymptotic(FUNC func1, FUNC func2, double x, double dx)
{
  double test = 1.0;
  while (1)
  {

    if (x==0.0)  test = 1.0;
    else         test = fabs(dx/x);
    if (test <= PRECISION) return(x);

    while (func1(x) != func2(x))
      x += dx;
    dx *= -0.1;

    if (x==0.0)  test = 1.0;
    else         test = fabs(dx/x);
    if (test <= PRECISION) return(x);

    while (func1(x) == func2(x))
      x += dx;
    dx *= -0.1;
  }
}

} // namespace Device
} // namespace Xyce

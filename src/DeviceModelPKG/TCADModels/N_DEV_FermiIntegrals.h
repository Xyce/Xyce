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
// Purpose        : This file contains the fermi integral functors
//
// Special Notes  :
//
// Creator        : Lawrence C Musson, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/05/2014
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_FermiIntegrals_h
#define Xyce_N_DEV_FermiIntegrals_h

// ---------- Standard Includes ----------

#include <cmath>
#include <time.h>

#include <Sacado_No_Kokkos.hpp>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_MaterialSupport.h>
#include <N_DEV_BernouliSupport.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_ScalingVars.h>

// ---------- Forward Declarations ----------

typedef Sacado::Fad::SFad<double,10> pdeFadType;

namespace Xyce 
{
  namespace Device 
  {


//-----------------------------------------------------------------------------
// Class         : inverse_fermi_one_half_N
// Purpose       : inverse fermi-dirac integral.  Implemented as a functor.
// Special Notes :
// Creator       : Eric Keiter, SNL
// Creation Date : 07/01/11
//-----------------------------------------------------------------------------
class inverse_fermi_one_half_N
{
  private:
    double d__1, d__2, d__3;
    double a1, a2, a3, a4, a5, a6, a7, a8, x10, y10, yp10, x20, y20, yp20, c1, c2;
    double pi, delx, dely;

  public:
    inverse_fermi_one_half_N () // this stuff gets called 1x.
    {
      double c_b2 = 4.0/3.0;
      pi = 2.0*asin(1.0);

      a1 = sqrt(2.0) / 4.0;
      a2 = 0.1875 - sqrt(3.0) / 9.0;
      a3 = sqrt(2.0) * 5.0 / 48.0 + 0.125 - sqrt(6.0) / 9.0;
      a4 = sqrt(2.0) * 5.0 / 32.0 + 1585.0/6912.0 - sqrt(3.0) * 5.0/24.0 - sqrt(5.0) / 25.0;
      d__1 = sqrt(pi) * 3.0/4.0;
      a5 = pow(d__1, c_b2);
      a6 = 4.0/3.0;
      a7 = pi   * pi   / 6.0;
      a8 = 1.0/3.0;
      x10 = 7.5;
      d__3 = x10*x10;
      y10 = log(x10) + a1 * x10 + a2 * (x10*x10) + a3*(x10*x10*x10) + a4*(d__3 * d__3);
      yp10 = 1.0/x10 + a1 + a2 * 2.0 * x10 + a3 * 3.0 * (x10 * x10) + a4 * 4.0 * (x10 * x10 * x10);
      x20 = 8.5;
      y20 = sqrt(a5 * pow(x20, a6) - a7);
      yp20 = 0.5 / sqrt(a5 * pow(x20, a6) - a7) * a6 * a5 * pow(x20, a8);
      delx = 0.5;
      dely = y20 - y10;
      c1 = dely * 0.5 / (delx * delx) - yp10 * 0.75 / delx - yp20 * 0.25 / delx;
      c2 = dely * 0.5 / (delx * delx) - yp20 * 0.75 / delx - yp10 * 0.25 / delx;
    }

    template <typename ScalarT>
      ScalarT operator()(const ScalarT & ratio)
      {
        ScalarT ret_val = 0.0;
        ScalarT tempVal = 0.0;

        // Joyce-Dixon expressions as used in Medici
        if (ratio > 0.0 && ratio <= 7.5)
        {
          tempVal = ratio*ratio;
          ret_val = log(ratio) + a1 * ratio + a2*(ratio*ratio) + a3*(ratio*ratio*ratio) + a4*(tempVal*tempVal);
        }

        // These next two clauses from Sam Myers
        if (ratio > 7.5 && ratio <= 8.0)
        {
          ScalarT diff = ratio - 7.5;
          ret_val = y10 + yp10*diff + c1*(diff*diff);
        }
        if (ratio > 8. && ratio < 8.5)
        {
          ScalarT diff = 8.5-ratio;
          ret_val = y20 - yp20*diff - c2*(diff*diff);
        }
        if (ratio >= 8.5)
        {
          ret_val = sqrt(a5 * pow(ratio, a6) - a7);
        }
        return ret_val;
      }
};


//-----------------------------------------------------------------------------
// Class         : fermi_one_half
// Purpose       : fermi-dirac integral.  Implemented as a functor.
// Special Notes :
// Creator       : Lawrence C Musson
// Creation Date : 06/04/14
//-----------------------------------------------------------------------------
class fermi_one_half
{
  private:

  double pi;

  public:
    fermi_one_half () 
    {
      pi = 4.0*atan(1.0);
    }

    template <typename ScalarT>
      ScalarT operator()(const ScalarT & arg)
      {

        // Reference:  "The Approximation of the Fermi-Dirac Integral F1/2(eta)"
        // by D. Bednarczyk and J. Bednarczyk, Physics Letters, Vol. 64A, No. 4,
        // 9 January 1978, pp. 409-410.
        ScalarT nu_eta = pow(arg, 4.0) + 50.0 +
          33.6*arg*(1.0 - 0.68*exp(-0.17*pow(arg+1.0,2)));

        ScalarT xi = 3.0*sqrt(pi)/(4.0*pow(nu_eta,0.375));

        ScalarT neta_uarg = 1.0/(exp(-arg)+xi);

        return neta_uarg;
      }

};

//-----------------------------------------------------------------------------
// Class         : fermi_minus_one_half
// Purpose       : fermi-dirac integral.  Implemented as a functor.
// Special Notes :
// Creator       : Lawrence C Musson
// Creation Date : 06/04/14
//-----------------------------------------------------------------------------
//template <typename ScalarT> 
  class fermi_minus_one_half
  {

    //Halen Pulfrey approximation
  private:

    double a1[7],a2[9],a3[9],n1[7];

  public:
    fermi_minus_one_half () 
    {
      a1[0] =  1.12837;
      a1[1] = -0.470698;
      a1[2] = -0.453108;
      a1[3] = -228.975;
      a1[4] =  8303.5;
      a1[5] = -118124.0;
      a1[6] =  632895.0;
      
      a2[0] =  0.604856;
      a2[1] =  0.380080;
      a2[2] =  0.059320;
      a2[3] = -0.014526;
      a2[4] = -0.004222;
      a2[5] =  0.001335;
      a2[6] =  0.000291;
      a2[7] = -0.000159;
      a2[8] =  0.000018;

      a3[0] =  0.638086;
      a3[1] =  0.292266;
      a3[2] =  0.159486;
      a3[3] = -0.077691;
      a3[4] =  0.018650;
      a3[5] = -0.002736;
      a3[6] =  0.000249;
      a3[7] = -0.000013;
      a3[8] =  2.9814e-07;

      // For negative arguments
      n1[0] = 0.999909;
      n1[1] = 0.706781;
      n1[2] = 0.572752;
      n1[3] = 0.466318;
      n1[4] = 0.324511;
      n1[5] = 0.152889;
      n1[6] = 0.033673;
    }

      template <typename ScalarT>
      ScalarT operator()(const ScalarT & arg)
      {

      ScalarT ret_val = 0.0;

      if (arg <= 0.0)
        {
          double sign = 1.0;
          double index;
          for (int i=0; i < 7; ++i)
            {
              index = i;
              ret_val += sign*n1[i]*std::exp((index+1.0)*arg);
              sign *= -1.0;
            }
        }
      else if (arg >= 5.0)
        {
          ScalarT mult = std::sqrt(arg);
          for (int i=0; i < 7; ++i)
            {
              ret_val += a1[i] / std::pow(arg, 2.0*i);
            }
          ret_val = mult * ret_val;
        }
      else if (arg < 2.5)
        {
          for (int i=0; i < 9; ++i)
            {
              ret_val += a2[i] * std::pow(arg, i);
            }
        }
      else // 2.5 <= arg < 5.0
        {
          for (int i=0; i < 9; ++i)
            {
              ret_val += a3[i] * std::pow(arg, i);
            }
        }
       
        return ret_val;
      }

};




  } //Device namespace
}  //Xyce namespace

#endif  //Xyce_N_DEV_FermiIntegrals_h


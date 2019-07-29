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
// Purpose        : This file contains the PDE device instance base class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DevicePDEInstance_h
#define Xyce_N_DEV_DevicePDEInstance_h

// ---------- Standard Includes ----------
#include <N_UTL_Math.h>
#include <time.h>

#include <Sacado_No_Kokkos.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_MaterialSupport.h>
#include <N_DEV_BernouliSupport.h>
#include <N_DEV_Const.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_DopeInfo.h>
#include <N_DEV_ScalingVars.h>
#include <N_DEV_FermiIntegrals.h>

// ---------- Forward Declarations ----------

typedef Sacado::Fad::SFad<double,10> pdeFadType;

namespace Xyce {
namespace Device {


//-----------------------------------------------------------------------------
// Class         : DevicePDEInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DevicePDEInstance : public DeviceInstance
{
  public:
    DevicePDEInstance(
        const InstanceBlock &       IB,
        ParametricData<void> &      parametric_data,
        const FactoryBlock &        factory_block);

    virtual ~DevicePDEInstance  () {};

    bool isPDEDevice() const { return true; }

  private:
    DevicePDEInstance(const DevicePDEInstance & right);
    DevicePDEInstance &operator=(const DevicePDEInstance & right);

  public:
    // Fermi-Dirac integral
    double fermi_one_half_B(double arg)
    {
      // Reference:  "The Approximation of the Fermi-Dirac Integral F1/2(eta)"
      // by D. Bednarczyk and J. Bednarczyk, Physics Letters, Vol. 64A, No. 4,
      // 9 January 1978, pp. 409-410.
      double pi = 4.0*atan(1.0);

      double nu_eta = pow(arg, 4.0) + 50.0 +
        33.6*arg*(1.0 - 0.68*exp(-0.17*pow(arg+1.0,2)));

      double xi = 3.0*sqrt(pi)/(4.0*pow(nu_eta,0.375));

      return 1.0/(exp(-arg)+xi);
    }

    double getVoltDepHoleDens (double Vmin, double V, double Na)
    {
      return  Na * exp ( std::min(CONSTMAX_EXP_ARG, ((Vmin-V)/Ut)) );
    }

    double getVoltDepElecDens (double Vmax, double V, double Nd)
    {
      return  Nd * exp ( std::min(CONSTMAX_EXP_ARG, ((V-Vmax)/Ut)) );
    }

    double aux1 (double x);
    double aux2 (double x);


    double pd1aux1(double x);
    double pd1aux2(double x);

    double Jn (double n1, double n2, double E, double u, double h);

    double dJndV1 (double n1, double n2, double E, double u, double h);
    double dJndV2 (double n1, double n2, double E, double u, double h);
    double dJndn1 (double n1, double n2, double E, double u, double h);
    double dJndn2 (double n1, double n2, double E, double u, double h);

    double Jp (double p1, double p2, double E, double u, double h);

    double dJpdV1 (double p1, double p2, double E, double u, double h);
    double dJpdV2 (double p1, double p2, double E, double u, double h);
    double dJpdn1 (double p1, double p2, double E, double u, double h);
    double dJpdn2 (double p1, double p2, double E, double u, double h);

    // charge dependent current density calculations
    double J_qdep (double n1, double n2, double E, double u, double h, int z);

    pdeFadType aux1 (pdeFadType & x)
    {
      pdeFadType retVal=0.0;
      if      (x < -bernSupport.bp0_MISC) x = -bernSupport.bp0_MISC;
      else if (x >  bernSupport.bp0_MISC) x =  bernSupport.bp0_MISC;

      if (x <= bernSupport.bp0_AUX1) retVal=(x / sinh(x));
      else if (x <= bernSupport.bp1_AUX1) retVal=(1 - x*x/6.0*(1.0 - 7.0*x*x/60.0));
      else retVal=(x / sinh(x));

      return retVal;
    }

    pdeFadType aux2 (pdeFadType & x)
    {
      pdeFadType retVal=0.0;

      if (x <= bernSupport.bp0_AUX2) retVal=(1.0);
      else if (x <= bernSupport.bp1_AUX2) retVal=(1.0 / (1.0 + exp(x)));
      else if (x <= bernSupport.bp2_AUX2) retVal=(exp(-x));
      else retVal=(0.0);

      return retVal;
    }


    pdeFadType nMidpoint(pdeFadType & n1, pdeFadType & n2, pdeFadType & E, double h, int z);

    double J_qdep (double n1, double n2, double E, pdeFadType & u, double h, int z)
    { return J_qdep (n1, n2, E, u.val(), h, z); }

    double dJdV1_qdep (double n1, double n2, double E, double u, double h, int z);
    double dJdV2_qdep (double n1, double n2, double E, double u, double h, int z);
    double dJdn1_qdep (double n1, double n2, double E, double u, double h, int z);
    double dJdn2_qdep (double n1, double n2, double E, double u, double h, int z);

    double dJdV1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdV2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdn1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdn2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdp1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdp2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdbm1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdbm2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdpp1_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    double dJdpp2_qdep (double n1, double n2, double E, const pdeFadType & u, double h, int z);
    //

    double erf(double x);
    double pd1erf(double x);

    const std::string timeDateStamp();
    const std::string tecplotTimeDateStamp();

    // np0 calculation
    template <typename ScalarT>
    ScalarT np0_calculation( ScalarT const& elec_dens,
                             ScalarT const& hole_dens,
                             ScalarT const& Ni,
                             ScalarT const& cond_band,
                             ScalarT const& vale_band,
                             ScalarT const& eff_dens_cond,
                             ScalarT const& eff_dens_vale,
                             ScalarT const& temp);

    // equilbrium concentrations: (using FD)
    template <typename ScalarT>
    void n0_and_p0 ( ScalarT const& elec_dens,
                     ScalarT const& hole_dens,
                     ScalarT const& Ni,
                     ScalarT const& cond_band,
                     ScalarT const& vale_band,
                     ScalarT const& eff_dens_cond,
                     ScalarT const& eff_dens_vale,
                     ScalarT const& temp,
                     ScalarT & n0,
                     ScalarT & p0);

    //Fermi-Dirac fluxes
    template <typename ScalarT>
      ScalarT FDCarrierFlux( ScalarT n1,
                             ScalarT n2,
                             ScalarT V1,
                             ScalarT V2,
                             ScalarT mu,
                             ScalarT temp,
                             double h,
                             double z,
                             double DOS);
    
  public:
    // physical constants:
    double Temp;     // operating temperature           (K)
    double charge;   // electron charge                 (C)
    double kb;       // boltzmann's constant            (J/K)
    double Vt;       // thermal voltage
    double Ut;       // thermal voltage, scaled.
    double e0;       // permittivity of vacuum          (F/cm)
    double eSi;      // relative permittivity of Si
    double eSiO2;    // relative permittivity of SiO2
    double eps;      // permittivity of Si              (F/cm)
    double Ni;       // intrinsic concentration of Si   (cm^-3)
    double h_planck; // planck's constant
    double e_mass;   // electron mass

    // scaling variables:
    double x0_user;  // distance scaling, as set by the user (cm)
    double C0_user;  // concentration scaling, as set by the user (cm^-3)
    double t0_user;  // time scaling, as set by the user (sec)

    ScalingVars scalingVars;

    std::map<std::string, DopeInfo *> dopeInfoMap;

    // continuation parameters:
    double maxVoltDelta;
    bool enableContinuationCalled;
    double continuationAlpha;

    bool sensOn;
    bool sensProcess;
    bool meshSensMod;
    bool dopingSensMod;
    bool photogenSensMod;

    std::string mobModelName;
    bool fieldDependentMobility;
    bool fieldDependentMobilityGiven;
    std::string bulkMaterial;
    bool variablesScaled;

    //MaterialSupport matSupport;
    BernouliSupport bernSupport;

    std::string outputName;  // added to remove the Y%PDE% prefix.

    // inverse fermi integral function functor.
    inverse_fermi_one_half_N fdinvObj;
    fermi_one_half fonehalfObj;
    fermi_minus_one_half fminusonehalfObj;
  protected:

  private:
    template <typename T> int sgn(T val)
    {
      return (val > T(0)) - (val < T(0));
    }

};

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::timeDateStamp_
// Purpose       : get current date and time and format for .PRINT output
// Special Notes : inline
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
inline const std::string DevicePDEInstance::timeDateStamp()
{
  const time_t now = time( NULL );
  char timeDate[ 80 ];

  // format for output
  strftime( timeDate, 80, "TIME='%I:%M:%S %p' DATE='%b %d, %Y' ",
      localtime( &now ) );

  return std::string( timeDate );
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::tecplotTimeDateStamp_
// Purpose       : Get current date and time and format for .PRINT output
// Special Notes : tecplot version of timeDateStamp_.
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 9/6/04
//-----------------------------------------------------------------------------
inline const std::string DevicePDEInstance::tecplotTimeDateStamp()
{
  const time_t now = time( NULL );
  char timeDate[ 80 ];

  // format for output
  strftime( timeDate, 80, "TIME= \" %I:%M:%S %p %b %d, %Y \" ",
      localtime( &now ) );

  return std::string( timeDate );
}



//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::np0_calculation
// Purpose       : np0 is calculated using Fermi-Dirac statistics
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT DevicePDEInstance::np0_calculation
( ScalarT const& elec_dens,
  ScalarT const& hole_dens,
  ScalarT const& Ni,
  ScalarT const& cond_band,
  ScalarT const& vale_band,
  ScalarT const& eff_dens_cond,
  ScalarT const& eff_dens_vale,
  ScalarT const& temp)
{
  ScalarT product = 0.0;
  ScalarT a_ = 0.0;
  ScalarT eta_ = 0.0;
  ScalarT n0 = 0.0;
  ScalarT p0 = 0.0;
  ScalarT kbq = 8.6173324e-5; // boltzmann's constant in eV K^-1

  if (elec_dens > hole_dens)
  {
    a_ = elec_dens - hole_dens;

    // Solve for equilibrium electron concentration
    n0 = (a_ + std::sqrt(a_*a_ + 4.0*Ni*Ni))/2.0;
    ScalarT n_ratio = n0/eff_dens_cond;
    eta_ = fdinvObj(n_ratio);
    ScalarT fermi_lev = cond_band + (kbq *temp)*eta_;
    p0 = eff_dens_vale * exp((vale_band - fermi_lev)/(kbq*temp));
    n0 = a_ + p0;
    product = n0*p0;
  }
  else
  {
    a_ = hole_dens - elec_dens;

    // Solve for equilibrium hole concentration
    p0 = (a_ + std::sqrt(a_*a_ + 4.0*Ni*Ni))/2.0;
    ScalarT p_ratio = p0/eff_dens_vale;
    eta_ = fdinvObj(p_ratio);
    ScalarT fermi_lev = vale_band - (kbq*temp)*eta_;
    n0 = eff_dens_cond * exp((fermi_lev - cond_band)/(kbq*temp));
    p0 = a_ + n0;
    product = n0 * p0;
  }

  return product;
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::n0_and_p0
//
// Purpose       : computes equilibrium concentrations of n and p using
//                 Fermi-Dirac statistics.
//
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 7/25/11
//-----------------------------------------------------------------------------
template <typename ScalarT>
void DevicePDEInstance::n0_and_p0
( ScalarT const& elec_dens,
  ScalarT const& hole_dens,
  ScalarT const& Ni,
  ScalarT const& cond_band,
  ScalarT const& vale_band,
  ScalarT const& eff_dens_cond,
  ScalarT const& eff_dens_vale,
  ScalarT const& temp,
  ScalarT & n0,
  ScalarT & p0)
{
  ScalarT a_ = 0.0;
  ScalarT eta_ = 0.0;
  ScalarT kbq = 8.6173324e-5; // boltzmann's constant in eV K^-1

  if (elec_dens > hole_dens)
  {
    a_ = elec_dens - hole_dens;

    // Solve for equilibrium electron concentration
    n0 = (a_ + std::sqrt(a_*a_ + 4.0*Ni*Ni))/2.0;
    ScalarT n_ratio = n0/eff_dens_cond;
    eta_ = fdinvObj(n_ratio);
    ScalarT fermi_lev = cond_band + (kbq *temp)*eta_;
    p0 = eff_dens_vale * exp((vale_band - fermi_lev)/(kbq*temp));
    n0 = a_ + p0;
  }
  else
  {
    a_ = hole_dens - elec_dens;

    // Solve for equilibrium hole concentration
    p0 = (a_ + std::sqrt(a_*a_ + 4.0*Ni*Ni))/2.0;
    ScalarT p_ratio = p0/eff_dens_vale;
    eta_ = fdinvObj(p_ratio);
    ScalarT fermi_lev = vale_band - (kbq*temp)*eta_;
    n0 = eff_dens_cond * exp((fermi_lev - cond_band)/(kbq*temp));
    p0 = a_ + n0;
  }
}

//-----------------------------------------------------------------------------
// Function      : DevicePDEInstance::FDCarrierFlux
//
// Purpose       : computes carrier fluxes with Fermi-Dirac statistics
//
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 9/30/14
//-----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT DevicePDEInstance::FDCarrierFlux(ScalarT n1,
                                         ScalarT n2,
                                         ScalarT V1,
                                         ScalarT V2,
                                         ScalarT mu,
                                         ScalarT temp,
                                         double h,
                                         double z,
                                         double DOS)
{

  ScalarT kbq = 8.6173324e-5; // boltzmann's constant in eV K^-1

  ScalarT kT = kbq*temp;

  ScalarT fermiArg1 = n1/DOS;

  ScalarT fermiArg2 = n2/DOS;

  ScalarT fdFactor = -0.5*z*(fdinvObj(fermiArg2) - fdinvObj(fermiArg1));

  ScalarT commonFactor = (V1-V2)/(2.0*kT);

  //If there's no field, it could present a problem

  ScalarT buffer = kT*1.e-3;

  if(-z*commonFactor < 0.0)buffer *= -1.0;

  ScalarT fluxDenominator = exp(-z*commonFactor+buffer/(2.0*kT)) - exp(z*commonFactor-buffer/(2.0*kT));

  ScalarT concProduct = n1*n2;

  if(n1*n2 < 0.0)concProduct = abs(buffer);

  ScalarT fluxPrefactor = -mu*(-z*(V1-V2)+buffer)/h*sqrt(concProduct);

  ScalarT fluxNumerator = exp(-z*(commonFactor + fdFactor)) - 
                          exp( z*(commonFactor + fdFactor));

  //If commonFactor AND commonFactor relative to fdFactor becomes large enough, the formulas can produce
  //overflow.  However, really, fluxNumerator/fluxDenominator -> 1.0.  Capture this event.

  if(abs(commonFactor) > 50.0 || abs(fdFactor) > 50.00)
    {
      fluxNumerator = 1.0;
      if(fdFactor < 0.0)
        fluxDenominator = -1.0;
      else
        fluxDenominator = 1.0;
    }

  ScalarT flux = fluxPrefactor*fluxNumerator/fluxDenominator;

  return flux;

}


} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_DevicePDEInstance_h

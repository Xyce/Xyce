//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Purpose        : strategy pattern for reactions with different temperature
//                  dependent rate constant styles
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 07/27/2006
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <string>
#include <vector>

#include <N_DEV_Const.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Specie.h>
#include <N_DEV_RateConstantCalculators.h>
#include <N_DEV_MaterialLayer.h>

#include<iostream>
#include <N_UTL_Math.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::SimpleRateCalculator
// Purpose       : constructor for "simple" reaction rate calculator
//                 This is one that has no temperature dependence, and scales
//                 as concentration*time
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  SimpleRateCalculator::SimpleRateCalculator(double k, double C0, double t0,
                                             double x0)
    : K(k)
  {
    setScaleFactors(C0,t0,x0);
    myReactionName = "SimpleRateCalculator";
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::SimpleRateCalculator
// Purpose       : Copy constructor for "simple" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  SimpleRateCalculator::SimpleRateCalculator(SimpleRateCalculator &right)
    : myReactionName(right.myReactionName),
      K(right.K),
      rk0(right.rk0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::Clone
// Purpose       : Copy self operation for "simple" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  SimpleRateCalculator *SimpleRateCalculator::Clone()
  {
    return new SimpleRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::computeRateConstant
// Purpose       : returns rate constant for simple rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double SimpleRateCalculator::computeRateConstant(double T)
  {
    return(K);
  }


//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::computeRateConstant
// Purpose       : returns rate constant for simple rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/14
//-----------------------------------------------------------------------------
  double SimpleRateCalculator::computeRateConstant(double T,
                                                   std::vector<double> &concs,
                                                   std::vector<double> &constant_vec)
  {
    return(K);
  }

//-----------------------------------------------------------------------------
// Function      : SimpleRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for simple rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double SimpleRateCalculator::rateConstantScaleFactor()
  {
    return (rk0);
  }


//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::CaptureRateCalculator
// Purpose       : constructor for capture reaction rate calculator
//                 This is one that models electron or hole capture reactions.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  CaptureRateCalculator::CaptureRateCalculator(double sigma, double v,
                                               double C0, double t0,
                                               double x0)
  {
    //K=sigma*v;  //velocity is incorporated later
    K=sigma;
    setScaleFactors(C0,t0,x0);
    myReactionName = "CaptureRateCalculator";
  }

//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::CaptureRateCalculator
// Purpose       : Copy constructor for "capture" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  CaptureRateCalculator::CaptureRateCalculator(CaptureRateCalculator &right)
    :myReactionName(right.myReactionName),
     K(right.K),
     rk0(right.rk0)
  {
  }
//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::Clone
// Purpose       : Copy self operation for "capture" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  CaptureRateCalculator *CaptureRateCalculator::Clone()
  {
    return new CaptureRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::computeRateConstant
// Purpose       : returns rate constant for capture rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double CaptureRateCalculator::computeRateConstant(double T)
  {
    return(K);
  }


//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::computeRateConstant
// Purpose       : returns rate constant for capture rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double CaptureRateCalculator::computeRateConstant(double T,
                                                    std::vector<double> &concs,
                                                    std::vector<double> &constant_vec)
  {
    return(K);
  }

//-----------------------------------------------------------------------------
// Function      : CaptureRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for capture rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double CaptureRateCalculator::rateConstantScaleFactor()
  {
    return (rk0);
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::EmissionRateCalculator
// Purpose       : constructor for emission reaction rate calculator
//                 This is one that models electron or hole emission reactions.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  EmissionRateCalculator::EmissionRateCalculator(double sigma, double v,
                                               double N, double Energy,
                                               double C0, double t0,
                                               double x0)
    : E(Energy)
  {
    //K_f=sigma*v*N;
    K_f=sigma;  //velocity and DOS are incorporated later

    setScaleFactors(C0,t0,x0);
    myReactionName = "Emission Rate Calculator";
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::EmissionRateCalculator
// Purpose       : Copy constructor for "emission" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  EmissionRateCalculator::EmissionRateCalculator(EmissionRateCalculator &right)
    :myReactionName(right.myReactionName),
     K_f(right.K_f),
     E(right.E),
     T0(right.T0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::Clone
// Purpose       : Copy self operation for "emission" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  EmissionRateCalculator *EmissionRateCalculator::Clone()
  {
    return new EmissionRateCalculator(*this);
  }
//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::computeRateConstant
// Purpose       : returns rate constant for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double EmissionRateCalculator::computeRateConstant(double T)
  {
    double KbT=CONSTboltz*T/CONSTQ;
    return(K_f*exp(-E/KbT));
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::computeRateConstant
// Purpose       : returns rate constant for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double EmissionRateCalculator::computeRateConstant(double T,
                                                     std::vector<double> &concs,
                                                     std::vector<double> &constant_vec)
  {
    double KbT=CONSTboltz*T/CONSTQ;
    return(K_f*exp(-E/KbT));
  }

//-----------------------------------------------------------------------------
// Function      : EmissionRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double EmissionRateCalculator::rateConstantScaleFactor()
  {
    return (T0);
  }

//-----------------------------------------------------------------------------
// Function      : FDEmissionRateCalculator::FDEmissionRateCalculator
// Purpose       : constructor for emission reaction rate calculator
//                 This is one that models electron or hole emission reactions.
//                 Fermi-Dirac stats are included in this.
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  FDEmissionRateCalculator::FDEmissionRateCalculator(double sigma,
                                                     double Energy,
                                                     double v,
                                                     double C0, double t0,
                                                     double x0)
    : E(Energy)
  {
    K_f=sigma;  //velocity and other things are accounted later
    setScaleFactors(C0,t0,x0);
    myReactionName = "Emission Rate Calculator";
  }

//-----------------------------------------------------------------------------
// Function      : FDEmissionRateCalculator::FDEmissionRateCalculator
// Purpose       : Copy constructor for "emission" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  FDEmissionRateCalculator::FDEmissionRateCalculator(FDEmissionRateCalculator &right)
    :myReactionName(right.myReactionName),
     K_f(right.K_f),
     E(right.E),
     T0(right.T0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : FDEmissionRateCalculator::Clone
// Purpose       : Copy self operation for "emission" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  FDEmissionRateCalculator *FDEmissionRateCalculator::Clone()
  {
    return new FDEmissionRateCalculator(*this);
  }
//-----------------------------------------------------------------------------
// Function      : FDEmissionRateCalculator::computeRateConstant
// Purpose       : returns rate constant for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double FDEmissionRateCalculator::computeRateConstant(double T)
  {
    return(K_f);
  }

//-----------------------------------------------------------------------------
// Function      : FDEmissionRateCalculator::computeRateConstant
// Purpose       : returns rate constant for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double FDEmissionRateCalculator::computeRateConstant(double T,
                                                     std::vector<double> &concs,
                                                     std::vector<double> &constant_vec)
  {
    return(K_f);
  }

//-----------------------------------------------------------------------------
// Function      : FDEmissionRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for emission rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double FDEmissionRateCalculator::rateConstantScaleFactor()
  {
    return (T0);
  }


//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::ComplexRateCalculator
// Purpose       : constructor for "complex" reaction rate calculator
//                 This is one that models two discrete species forming a
//                 complex, e.g.:
//                       V0 + VM -> VVM
//
// Special Notes : For this to work, there must either be two species in the
//                 reactants list, or one species with 2.0 as the stochiometric
//                 coefficient.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexRateCalculator::ComplexRateCalculator(
    std::vector<Specie> &VariableSpecies, std::vector<Specie> &ConstantSpecies,
     std::vector< std::pair<int,double> > &Reactants,
    double iC0, double t0, double x0)
    :coulombAttraction(false),
     classC0(iC0)
      
  {
    int ij;

    myReactionName = "Complex Rate Calculator";

    // Check assumptions:
    if ( ! ((Reactants.size() == 1 && Reactants[0].second == 2.0) ||
            (Reactants.size() == 2 && Reactants[0].second == 1.0 &&
             Reactants[1].second == 1.0)))
    {
      Report::UserError msg;
      msg << "Invalid attempt to use complex rate method.  This method is only valid for binary complexing reactions:\n";
      if (Reactants.size() == 1)
      {
        msg << "   Only one reactant specified, but its stoichimetric coefficient is not 2.";
      }
      else if (Reactants.size() == 2)
      {
        msg << "   Two reactants specified, but both stoichimetric coefficient are not 1.";
      }
      else
      {
        msg << "   More than two reactants specified.";
      }
    }

    if (Reactants[0].first >= 0)
      Specie1 = &(VariableSpecies[Reactants[0].first]);
    else
      Specie1 = &(ConstantSpecies[-(Reactants[0].first+1)]);

    // Handle case where there's only one species with a coefficient of 2.0
    if (Reactants.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
    else
    {
      if (Reactants[1].first >= 0)
        Specie2 = &(VariableSpecies[Reactants[1].first]);
      else
        Specie2 = &(ConstantSpecies[-(Reactants[1].first+1)]);
    }
    ij    =Specie1->getChargeState();
    ij   *=Specie2->getChargeState();


    // Only divide reaction_distance_factor by T in one special case

    //NOTE:  The reaction_distance_factor that is being calculated here no longer
    //has a distance in it.  This is because this coefficient is being calculated
    //before there is a known bulk material.  And the reaction distance is material
    //dependent and also a function of concentration in gallium arsenide. -- LCM
    Tdep=false;
    if (ij>0)
      reaction_distance_factor=0.0;
    else if (ij == 0)
      {
        coulombAttraction = false;
        reaction_distance_factor = 4*M_PI;
        chargeNumberProduct = 0.0;
      }
    else
    {
      coulombAttraction = true;
      reaction_distance_factor = 4*M_PI;
      chargeNumberProduct = -(double)ij;
      Tdep=true;
    }

    setScaleFactors(classC0,t0,x0);
  }


//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::ComplexRateCalculator
// Purpose       : Copy constructor for "complex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexRateCalculator::ComplexRateCalculator(ComplexRateCalculator &right)
    :myReactionName(right.myReactionName),
     Specie1(right.Specie1),
     Specie2(right.Specie2),
     reaction_distance_factor(right.reaction_distance_factor),
     Tdep(right.Tdep),
     rk0(right.rk0),
     coulombAttraction(right.coulombAttraction),
     classC0(right.classC0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::Clone
// Purpose       : Copy self operation for "complex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexRateCalculator *ComplexRateCalculator::Clone()
  {
    return new ComplexRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::computeRateConstant
// Purpose       : returns rate constant for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double ComplexRateCalculator::computeRateConstant(double T)
  {

    if (Tdep)
      return(reaction_distance_factor/T);
    else
      return(reaction_distance_factor);


    /*  This is always done in reaction->getRate so that Bourgoin Corbett effects can 
        be properly incorporated. --LCM
    if (Tdep)
      return(reaction_distance_factor/T*(Specie1->getDiffusionCoefficient(T)
                                        +Specie2->getDiffusionCoefficient(T)));
    else
      return(reaction_distance_factor*(Specie1->getDiffusionCoefficient(T)
                                       +Specie2->getDiffusionCoefficient(T)));
    */

  }


//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::computeRateConstant
// Purpose       : returns rate constant for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double ComplexRateCalculator::computeRateConstant(double T,
                                                    std::vector<double> &concs,
                                                    std::vector<double> &constant_vec)
  {
    if (Tdep)
      return(reaction_distance_factor/T);
    else
      return(reaction_distance_factor);

    /*  This is always done in reaction->getRate so that Bourgoin Corbett effects can 
        be properly incorporated. --LCM
    if (Tdep)
      return(reaction_distance_factor/T*(Specie1->getDiffusionCoefficient(T,concs,constant_vec)
                                        +Specie2->getDiffusionCoefficient(T)));
    else
      return(reaction_distance_factor*(Specie1->getDiffusionCoefficient(T,concs,constant_vec)
                                       +Specie2->getDiffusionCoefficient(T,concs,constant_vec)));
    */

  }

//-----------------------------------------------------------------------------
// Function      : ComplexRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double ComplexRateCalculator::rateConstantScaleFactor()
  {
   return (rk0);
  }



//-----------------------------------------------------------------------------
// Function      : ComplexMultiplierRateCalculator::ComplexMultiplierRateCalculator
// Purpose       : constructor for "complex" reaction rate calculator
//                 This is one that models two discrete species forming a
//                 complex, e.g.:
//                       V0 + VM -> VVM
//                 This method multiplies the rate by a user specified constant
// Special Notes : For this to work, there must either be two species in the
//                 reactants list, or one species with 2.0 as the stochiometric
//                 coefficient.
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 06/26/2014
//-----------------------------------------------------------------------------
  ComplexMultiplierRateCalculator::ComplexMultiplierRateCalculator(
    std::vector<Specie> &VariableSpecies, std::vector<Specie> &ConstantSpecies,
     std::vector< std::pair<int,double> > &Reactants,
    double C0, double t0, double x0, double mult)
  {
    int ij;

    multiplier = mult;

    myReactionName = "Complex Multiplier Rate Calculator";
    // Check assumptions:
    if ( ! ((Reactants.size() == 1 && Reactants[0].second == 2.0) ||
            (Reactants.size() == 2 && Reactants[0].second == 1.0 &&
             Reactants[1].second == 1.0)))
    {
      Report::UserError msg;
      msg <<"Invalid attempt to use complex rate method.  This method is only valid for binary complexing reactions:\n";
      if (Reactants.size() == 1)
      {
        msg << "   Only one reactant specified, but its stoichimetric coefficient is not 2.";
      }
      else if (Reactants.size() == 2)
      {
        msg << "   Two reactants specified, but both stoichimetric coefficient are not 1.";
      }
      else
      {
        msg << "   More than two reactants specified.";
      }
    }

    if (Reactants[0].first >= 0)
      Specie1 = &(VariableSpecies[Reactants[0].first]);
    else
      Specie1 = &(ConstantSpecies[-(Reactants[0].first+1)]);

    // Handle case where there's only one species with a coefficient of 2.0
    if (Reactants.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
    else
    {
      if (Reactants[1].first >= 0)
        Specie2 = &(VariableSpecies[Reactants[1].first]);
      else
        Specie2 = &(ConstantSpecies[-(Reactants[1].first+1)]);
    }
    ij    =Specie1->getChargeState();
    ij   *=Specie2->getChargeState();



    // Only divide reaction_distance_factor by T in one special case

    //NOTE:  The reaction_distance_factor that is being calculated here no longer
    //has a distance in it.  This is because this coefficient is being calculated
    //before there is a known bulk material.  And the reaction distance is material
    //dependent and also a function of concentration in gallium arsenide. -- LCM
    Tdep=false;
    if (ij>0)
      reaction_distance_factor=0.0;
    else if (ij == 0)
      {
        coulombAttraction = false;
        //reaction_distance_factor = 4*M_PI*5.0e-8;
        reaction_distance_factor = 4*M_PI*multiplier;
        chargeNumberProduct = 0.0;
      }
    else
    {
      coulombAttraction = true;
      //reaction_distance_factor = 4*M_PI*1.4e-4*(-ij);
      reaction_distance_factor = 4*M_PI*multiplier;
      chargeNumberProduct = -(double)ij;
      Tdep=true;
    }

    setScaleFactors(C0,t0,x0);
  }


//-----------------------------------------------------------------------------
// Function      : ComplexMultiplierRateCalculator::ComplexMultiplierRateCalculator
// Purpose       : Copy constructor for "complex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexMultiplierRateCalculator::ComplexMultiplierRateCalculator(ComplexMultiplierRateCalculator &right)
    :myReactionName(right.myReactionName),
     Specie1(right.Specie1),
     Specie2(right.Specie2),
     reaction_distance_factor(right.reaction_distance_factor),
     Tdep(right.Tdep),
     rk0(right.rk0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : ComplexMultiplierRateCalculator::Clone
// Purpose       : Copy self operation for "complex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  ComplexMultiplierRateCalculator *ComplexMultiplierRateCalculator::Clone()
  {
    return new ComplexMultiplierRateCalculator(*this);
  }

//-----------------------------------------------------------------------------
// Function      : ComplexMultiplierRateCalculator::computeRateConstant
// Purpose       : returns rate constant for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double ComplexMultiplierRateCalculator::computeRateConstant(double T)
  {
    if (Tdep)
      return(multiplier*reaction_distance_factor/T*(Specie1->getDiffusionCoefficient(T)
                                                    +Specie2->getDiffusionCoefficient(T)));
    else
      return(multiplier*reaction_distance_factor*(Specie1->getDiffusionCoefficient(T)
                                                  +Specie2->getDiffusionCoefficient(T)));

  }


//-----------------------------------------------------------------------------
// Function      : ComplexMultiplierRateCalculator::computeRateConstant
// Purpose       : returns rate constant for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double ComplexMultiplierRateCalculator::computeRateConstant(double T,
                                                    std::vector<double> &concs,
                                                    std::vector<double> &constant_vec)
  {
    if (Tdep)
      return(multiplier*reaction_distance_factor/T*(Specie1->getDiffusionCoefficient(T,concs,constant_vec)
                                                    +Specie2->getDiffusionCoefficient(T)));
    else
      return(multiplier*reaction_distance_factor*(Specie1->getDiffusionCoefficient(T,concs,constant_vec)
                                                  +Specie2->getDiffusionCoefficient(T,concs,constant_vec)));

  }

//-----------------------------------------------------------------------------
// Function      : ComplexMultiplierRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for complex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double ComplexMultiplierRateCalculator::rateConstantScaleFactor()
  {
   return (rk0);
  }





//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::DecomplexRateCalculator
// Purpose       : constructor for "complex" reaction rate calculator
//                 This is one that models two discrete species decomposing from
//                 a complex, e.g.:
//                       VMM->V0 + VM
//
// Special Notes : For this to work, there must either be two species in the
//                 products list, or one species with 2.0 as the stochiometric
//                 coefficient.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 5/04/09
//-----------------------------------------------------------------------------
  DecomplexRateCalculator::DecomplexRateCalculator(
     std::vector<Specie> &VariableSpecies, std::vector<Specie> &ConstantSpecies,
     std::vector< std::pair<int,double> > &Reactants,
     std::vector< std::pair<int,double> > &Products,
     double bindingEnergy, double degenAB, double degenA, double degenB,
     double siliconConcentration,
     double C0, double t0, double x0)
    :deltaE(bindingEnergy),
     gammaA(degenA),
     gammaB(degenB),
     gammaAB(degenAB),
     concSi(siliconConcentration),
     c0(C0)
  {
    int ij;

    myReactionName = "Decomplex Rate Calculator";

    // Check assumptions:
    if ( ! ((Products.size() == 1 && Products[0].second == 2.0) ||
            (Products.size() == 2 && Products[0].second == 1.0 &&
             Products[1].second == 1.0)))
    {
      Report::UserError msg;
      msg << "Invalid attempt to use decomplex rate method.  This method is only valid for decomplexing reactions with two products:\n";
      if (Products.size() == 1)
      {
        msg << "   Only one product specified, but its stoichimetric coefficient is not 2.";
      }
      else if (Products.size() == 2)
      {
        msg << "   Two products specified, but both stoichimetric coefficient are not 1.";
      }
      else
      {
        msg << "   More than two products specified.";
      }
    }

    if (Products[0].first >= 0)
      Specie1 = &(VariableSpecies[Products[0].first]);
    else
      Specie1 = &(ConstantSpecies[-(Products[0].first+1)]);

    // Handle case where there's only one species with a coefficient of 2.0
    if (Products.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
    else
    {
      if (Products[1].first >= 0)
        Specie2 = &(VariableSpecies[Products[1].first]);
      else
        Specie2 = &(ConstantSpecies[-(Products[1].first+1)]);
    }
    ij    =Specie1->getChargeState();
    ij   *=Specie2->getChargeState();


    // Only divide reaction_distance_factor by T in one special case

    //NOTE:  The reaction_distance_factor that is being calculated here no longer
    //has a distance in it.  This is because this coefficient is being calculated
    //before there is a known bulk material.  And the reaction distance is material
    //dependent and also a function of concentration in gallium arsenide. -- LCM
    Tdep=false;
    if (ij>0)
      reaction_distance_factor=0.0;
    else if (ij == 0)
      {
        coulombAttraction = false;
        reaction_distance_factor = 4*M_PI;
        chargeNumberProduct = 0.0;
      }
    else
    {
      coulombAttraction = true;
      reaction_distance_factor = 4*M_PI;
      chargeNumberProduct = -(double)ij;
      Tdep=true;
    }

    setScaleFactors(C0,t0,x0);
  }


//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::DecomplexRateCalculator
// Purpose       : Copy constructor for "decomplex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 5/04/09
//-----------------------------------------------------------------------------
  DecomplexRateCalculator::DecomplexRateCalculator(DecomplexRateCalculator &right)
    :myReactionName(right.myReactionName),
     Specie1(right.Specie1),
     Specie2(right.Specie2),
     reaction_distance_factor(right.reaction_distance_factor),
     Tdep(right.Tdep),
     deltaE(right.deltaE),
     gammaA(right.gammaA),
     gammaB(right.gammaB),
     gammaAB(right.gammaAB),
     concSi(right.concSi),
     rk0(right.rk0),
     c0(right.c0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::Clone
// Purpose       : Copy self operation for "decomplex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  DecomplexRateCalculator *DecomplexRateCalculator::Clone()
  {
    return new DecomplexRateCalculator(*this);
  }


//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::computeRateConstant
// Purpose       : returns rate constant for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double DecomplexRateCalculator::computeRateConstant(double T)
  {
    double R;
    double KbT=CONSTboltz*T/CONSTQ;
    double k;
    double D1=Specie1->getDiffusionCoefficient(T);
    double D2=Specie2->getDiffusionCoefficient(T);
    if (Tdep)
      R=(reaction_distance_factor/T);
    else
      R=reaction_distance_factor;

    //k=(R*(D1 +D2)*(concSi)*((gammaA*gammaB)/gammaAB)*exp(-deltaE/KbT));

    k=(R*(concSi)*((gammaA*gammaB)/gammaAB)*exp(-deltaE/KbT));

    return k;
  }


//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::computeRateConstant
// Purpose       : returns rate constant for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
  double DecomplexRateCalculator::computeRateConstant(double T,
                                                      std::vector<double> &concs,
                                                      std::vector<double> &constant_vec)
  {
    double R;
    double KbT=CONSTboltz*T/CONSTQ;
    double k;
    double D1=Specie1->getDiffusionCoefficient(T,concs,constant_vec);
    double D2=Specie2->getDiffusionCoefficient(T,concs,constant_vec);
    if (Tdep)
      R=(reaction_distance_factor/T);
    else
      R=reaction_distance_factor;

    //k=(R*(D1 +D2)*(concSi)*((gammaA*gammaB)/gammaAB)*exp(-deltaE/KbT));

    k=(R*(concSi)*((gammaA*gammaB)/gammaAB)*exp(-deltaE/KbT));

    return k;
  }

//-----------------------------------------------------------------------------
// Function      : DecomplexRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 7/31/06
//-----------------------------------------------------------------------------
  double DecomplexRateCalculator::rateConstantScaleFactor()
  {
   return (rk0);
  }




//-----------------------------------------------------------------------------
// Function      : BourgoinCorbettRateCalculator::BourgoinCorbettRateCalculator
// Purpose       : constructor for "bourgoin corbett" reaction rate calculator
//                 This is one that models two discrete species decomposing from
//                 a complex, e.g.:
//                       SII0 + B0 + H -> BI0
//
// Special Notes : 
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
  BourgoinCorbettHoleRateCalculator::BourgoinCorbettHoleRateCalculator(
     std::vector<Specie> &VariableSpecies, std::vector<Specie> &ConstantSpecies,
     std::vector< std::pair<int,double> > &Reactants,
     std::vector< std::pair<int,double> > &Products,
     double sigma,
     double C0, double t0, double x0)
    :c0(C0)
  {
    int ij;

    // Check assumptions:
    if ( ! (Reactants.size() == 3 && Reactants[0].second == 1.0  &&
            Reactants[1].second == 1.0 && Reactants[2].second == 1.0) )
    {
      Report::UserError() << "Invalid attempt to use rate method.  This method is only valid for ternary reactions";
    }

    if (Reactants[0].first >= 0)
      Specie1 = &(VariableSpecies[Reactants[0].first]);
    else
      Specie1 = &(ConstantSpecies[-(Reactants[0].first+1)]);

    // Handle case where there's only one species with a coefficient of 2.0
    if (Reactants.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
    else
    {
      if (Reactants[1].first >= 0)
        Specie2 = &(VariableSpecies[Reactants[1].first]);
      else
        Specie2 = &(ConstantSpecies[-(Reactants[1].first+1)]);
    }
    ij    =Specie1->getChargeState();
    ij   *=Specie2->getChargeState();

    double carrierThermalVelocity;

    if(Reactants[2].first == -1)
      carrierThermalVelocity = 20335471.413078606;
    else
      carrierThermalVelocity = 16805108.930336751;


    // Only divide reaction_distance_factor by T in one special case
    Tdep=false;
    if (ij>0)
      reaction_distance_factor=0.0;
    else if (ij == 0)
      {
      reaction_distance_factor = 4*M_PI*5.0e-8*sigma* carrierThermalVelocity * 5.0e-8 * 5.0e-8/6.0;
      //reaction_distance_factor = 4*M_PI*5.43e-8;
      //reaction_distance_factor = 4*3.1416*5.43e-8;
      }
    else
    {
      reaction_distance_factor = 4*M_PI*1.4e-4*(-ij)*sigma*carrierThermalVelocity*5.0e-8*5.0e-8/6.0;
      //reaction_distance_factor = 4*M_PI*1.40419676681412915e-4*(-ij)*sigma*carrierThermalVelocity*5.0e-8*5.0e-8/6.0;
      Tdep=true;
    }

    setScaleFactors(C0,t0,x0);
  }


//-----------------------------------------------------------------------------
// Function      : BourgoinCorbettHoleRateCalculator::BourgoinCorbettHoleRateCalculator
// Purpose       : Copy constructor
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
  BourgoinCorbettHoleRateCalculator::BourgoinCorbettHoleRateCalculator(BourgoinCorbettHoleRateCalculator &right)
    :Specie1(right.Specie1),
     Specie2(right.Specie2),
     reaction_distance_factor(right.reaction_distance_factor),
     Tdep(right.Tdep),
     sigma(right.sigma),
     rk0(right.rk0),
     c0(right.c0)
  {
  }

//-----------------------------------------------------------------------------
// Function      : BourgoinCorbettHoleRateCalculator::Clone
// Purpose       : Copy self operation for "decomplex" reaction rate calculator
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
  BourgoinCorbettHoleRateCalculator *BourgoinCorbettHoleRateCalculator::Clone()
  {
    return new BourgoinCorbettHoleRateCalculator(*this);
  }


//-----------------------------------------------------------------------------
// Function      : BourgoinCorbettHoleRateCalculator::computeRateConstant
// Purpose       : returns rate constant for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
  double BourgoinCorbettHoleRateCalculator::computeRateConstant(double T)
  {
    double R;
    double KbT=CONSTboltz*T/CONSTQ;
    double k;
    double D1=Specie1->getDiffusionCoefficient(T);
    double D2=Specie2->getDiffusionCoefficient(T);
    if (Tdep)
      R=(reaction_distance_factor/T);
    else
      R=reaction_distance_factor;

    k=R;

    return k;
  }


//-----------------------------------------------------------------------------
// Function      : BourgoinCorbettHoleRateCalculator::computeRateConstant
// Purpose       : returns rate constant for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
  double BourgoinCorbettHoleRateCalculator::computeRateConstant(double T,
                                                                std::vector<double> &concs,
                                                                std::vector<double> &constant_vec)
  {
    double R;
    double KbT=CONSTboltz*T/CONSTQ;
    double k;
    double D1=Specie1->getDiffusionCoefficient(T);
    double D2=Specie2->getDiffusionCoefficient(T);
    if (Tdep)
      R=(reaction_distance_factor/T);
    else
      R=reaction_distance_factor;

    k=R;

    return k;
  }

//-----------------------------------------------------------------------------
// Function      : BourgoinCorbettHoleRateCalculator::rateConstantScaleFactor
// Purpose       : returns rate scaling factor for decomplex rate style
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
  double BourgoinCorbettHoleRateCalculator::rateConstantScaleFactor()
  {
   return (rk0);
  }



 
} // namespace Device
} // namespace Xyce

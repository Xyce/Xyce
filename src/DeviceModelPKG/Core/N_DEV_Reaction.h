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
// Purpose        : Provide a generic class for reactions using simple
//                  law-of-mass-action kinetics, i.e.:
//                       A+2B+C -> D + 3E + F
//                  implies that the reaction happens at a rate:
//                       reactionRate = k*[A][B]^2[C]
//                  where [X] means "concentration of species X"
//
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 03/20/06
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_DEV_REACTION_H
#define N_DEV_REACTION_H
#include <iosfwd>
#include <vector>

#include <N_DEV_Specie.h>
#include <N_DEV_RateConstantCalculators.h>
#include <N_DEV_FermiIntegrals.h>
#include <N_DEV_MaterialLayer.h>

// ---------- Forward Declarations ----------

typedef Sacado::ELRFad::DFad<double> FDFadType;

namespace Xyce {
namespace Device {

class Reaction
{
public:
  Reaction();
  Reaction(std::vector< std::pair<int,double> > & ,
           std::vector< std::pair<int,double> > &,
           double);
  Reaction(const Reaction &right);
  ~Reaction();
  void setRateConstant(double);
  void setRateConstantFromCalculator(double T);
  void setRateConstantFromCalculator(double T, 
                                     std::vector<double> &concs,
                                     std::vector<double> &constant_vec);
  void scaleRateConstant(double);
  void scaleRateConstantFromCalculator();
  void unscaleRateConstantFromCalculator();
  void setReactants(std::vector< std::pair<int,double> > &products );
  void addReactant(int species,double stoich);
  void setProducts(std::vector< std::pair<int,double> > &products );
  void addProduct(int species,double stoich);
  double getRate(std::vector<double> &concentrations,
                 std::vector<double> &constants);
  void getDdt(std::vector<double> &concentrations,
              std::vector<double> &constants,
              std::vector<double> &ddt);
  void getDRateDC(std::vector<double> &concentrations,
                  std::vector<double> &constants,
                  std::vector<double> &dratedc);
  void getDRateDConst(int constNum,
                      std::vector<double> &concentrations,
                      std::vector<double> &constants,
                      double &dratedc);
  void getJac(std::vector<double> &concentrations,
              std::vector<double> &constants,
              std::vector<std::vector<double> > &jac);
  void getDFdConst(int constantNumber,
                   std::vector<double> &concentrations,
                   std::vector<double> &constants,
                   std::vector<double> &dFdConst);

  void getJacobianVC(std::vector<double> &concentrations,
                     std::vector<double> &constants,
                     std::vector<std::vector<double> > &jac,
                     std::vector<double> &constVec);

  void output ( const std::vector<Specie> & species,
                std::ostream & os ) const;

  void setSimpleRateCalculator(double k, double C0, double t0, double x0);
  void setCaptureRateCalculator(double sigma, double v, double C0, double t0,
                                double x0);
  void setEmissionRateCalculator(double sigma, double v, double N,
                                 double Energy, double C0, double t0,
                                 double x0);
  void setFDEmissionRateCalculator(int carrierIndex, double sigma, double Energy, 
                                   double v, double C0, double t0, double x0);
  void setComplexRateCalculator(std::vector<Specie> &VariableSpecies,
                                std::vector<Specie> &ConstantSpecies,
                                double C0, double t0, double x0);
  void setComplexMultiplierRateCalculator(std::vector<Specie> &VariableSpecies,
                                          std::vector<Specie> &ConstantSpecies,
                                          double C0, double t0, double x0, 
                                          double multiplier);
  void setDecomplexRateCalculator(std::vector<Specie> &VariableSpecies,
                                  std::vector<Specie> &ConstantSpecies,
                                  double bindingEnergy,
                                  double gammaAB, double gammaA, double gammaB,
                                  double concSi,
                                  double C0, double t0, double x0);

  void setBourgoinCorbettHoleRateCalculator(std::vector<Specie> &VariableSpecies,
                                            std::vector<Specie> &ConstantSpecies,
                                            double sigma,double C0, double t0, double x0);

  int getCarrierEmissionIndex();

  void setTemperature(double T);

  void setMaterial(MaterialLayer *material, double Temp);
  void setCoefficient(double Temp);
  void setRxnVariableCoeffs(bool variableCoeff){variableRateCoefficient = variableCoeff;}

  template <class ScalarT>
  ScalarT getFDEmissionRate(std::vector<ScalarT> &concentrations,
                            std::vector<ScalarT> &constants);

  template <class ScalarT>
  ScalarT getRateVC(std::vector<ScalarT> &concentrations,
                    std::vector<ScalarT> &constants);

  inline void setScaleFactors(double C0i, double t0, double x0)
  {
    C0 = C0i;
    if (myRateCalc)
      myRateCalc->setScaleFactors(C0i,t0,x0);
  };

  Reaction & operator=(const Reaction & right);

private:
  void setDependency(int cSize);
  void setConstDependency(int cSize);
  std::vector< std::pair<int,double> > theReactants;
  std::vector< std::pair<int,double> > theProducts;
  double theRateConstant;
  int numconcs; // size of vector of concentrations
  int numconsts; // size of vector of constants
  int carrierEmissionIndex;
  std::vector<int> concDependency;
  std::vector<int> constDependency;
  RateCalculator *myRateCalc;
  double C0;
  double temperature,energy;
  inverse_fermi_one_half_N fdinvObj;
  MaterialLayer *material;
  bool variableRateCoefficient;
  std::string myReactionName;
  bool FADVectorsAllocated;
  std::vector<FDFadType> defects;
  std::vector<FDFadType> carriers;
  bool variableCoefficient;

  //Coefficient data
  int coefficientType;
  double tolerance;
  double charge;
  double peq;
  double constCoeff;
  double unshieldedLength;
  double latticeConstant;
  template <class ScalarT>
  ScalarT rxnCoefficient(std::vector<ScalarT> &concentrations,
                         std::vector<ScalarT> &constants);

  template <class ScalarT>
  ScalarT complexCoefficient(std::vector<ScalarT> &concentrations,
                             std::vector<ScalarT> &constants);

  Specie *Specie1,*Specie2;
  int chargeProduct;
  double diffusionCoefficient1,diffusionCoefficient2;
  //The following is needed for bourgoin corbett 
  int carrierBCIndex;
  double hopLength,thermalVelocity,sigmaBC;
  int carrierCharge;

  enum materialEnum_t {si,gaas,gan};
  materialEnum_t materialEnum;


};

//-----------------------------------------------------------------------------
// Function      : Reaction::setRateConstant
// Purpose       : Accessor function to set rate constant
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void Reaction::setRateConstant(double rateConst)
{
  theRateConstant=rateConst;
}

//-----------------------------------------------------------------------------
// Function      : Reaction::scaleRateConstant
// Purpose       : Accessor function to scale rate constant
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void Reaction::scaleRateConstant(double scalar)
{
  theRateConstant *= scalar;
}


//-----------------------------------------------------------------------------
// Function      : Reaction::getEmissionCarrierIndex
// Purpose       : Accessor function to emission carrier index
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 07/14/2014
//-----------------------------------------------------------------------------
inline int Reaction::getCarrierEmissionIndex()
{
  return carrierEmissionIndex;
}


//-----------------------------------------------------------------------------
// Function      : Reaction::setTemperature
// Purpose       : Accessor function to set the temperature
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 07/14/2014
//-----------------------------------------------------------------------------
inline void Reaction::setTemperature(double T)
{
  temperature = T;
}

//-----------------------------------------------------------------------------
// Function      : Reaction::getFDEmissionRate
// Purpose       : compute and return the carrier emission rate with FD stats.  
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 07/11/2014
//-----------------------------------------------------------------------------

template <class ScalarT>
ScalarT Reaction::getFDEmissionRate(std::vector<ScalarT> &concentrations,
                                    std::vector<ScalarT> &constants)
{

  ScalarT rRate=0.0;

  int rSize=theReactants.size();
  int pSize=theProducts.size();

  int i;
  int species;
  double stoich;
  ScalarT c;

  // Reaction rate determined by law of mass action

  ScalarT localCoeff = rxnCoefficient(concentrations, constants);

  rRate = theRateConstant*localCoeff;
  // product of concentrations of reactants raised to the power of
  // stoichiometric coefficient

  for (i=0;i<rSize;++i)
  {
    species=theReactants[i].first;
    stoich=theReactants[i].second;

    if (species>=0)
    {
      c=concentrations[species];
    }
    else
    {
      c=constants[-(species+1)];
    }
       
    if (stoich != 1.0)
    {
      rRate *= pow(c,stoich);
    }
    else
    {
      rRate *= c;
    }
  }

  double KbT=CONSTboltz*temperature/CONSTQ;

  double DOS;

  if(carrierEmissionIndex==0)
    DOS = material->Nc;
  else
    DOS = material->Nv;

  if(C0*constants[carrierEmissionIndex] <1.e12)
  {
    ScalarT exponential = exp(-energy/KbT);
    rRate *= DOS*exponential;
  }
  else
  {
    DOS /= C0;

    ScalarT ratio = constants[carrierEmissionIndex]/DOS;

    ScalarT iFI = KbT*fdinvObj(ratio);

    ScalarT exponential = exp(-(iFI + energy)/KbT);
   
    rRate*=C0*constants[carrierEmissionIndex]*exponential;
  }

  return rRate;

}


//-----------------------------------------------------------------------------
// Function      : Reaction::getRateVC
// Purpose       : compute and return the reaction rate.  This is normally the
//                 first step  of getDdt, but we might need this for
//                 other purposes (e.g. a kludge for getting electron/hole
//                 recombination/regeneration rates even when electrons and
//                 holes are constant species).  In this version the rate coefficients
//                 are CONCENTRAION AND temperature dependent.  Moreover, the method
//                 is templated so that it can be FADed cleanly.
// Special Notes : It is tacitly assumed that the concentrations and
//                 constants in the provided vector are in an order
//                 consistent with the labeling in the reactants and
//                 products vectors.  The contributions to the species
//                 time derivatives are summed into the provided output
//                 vector because it is assumed that there will be a network
//                 of reactions, not just one.
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 08/26/2014
//-----------------------------------------------------------------------------
template <class ScalarT>
ScalarT Reaction::getRateVC(std::vector<ScalarT> &concentrations,
                            std::vector<ScalarT> &constants)
{
  int rSize=theReactants.size();
  int pSize=theProducts.size();
  ScalarT reactionRate;
  int i;
  int species;
  double stoich;
  ScalarT c;

  // Reaction rate determined by law of mass action

  ScalarT localCoeff = rxnCoefficient(concentrations, constants);

  reactionRate=theRateConstant*localCoeff;

  // product of concentrations of reactants raised to the power of
  // stoichiometric coefficient

  for (i=0;i<rSize;++i)
  {
    species=theReactants[i].first;
    stoich=theReactants[i].second;

    if (species>=0)
    {
      c=concentrations[species];
    }
    else
    {
      c=constants[-(species+1)];
    }

    if (stoich != 1.0)
    {
      reactionRate *= pow(c,stoich);
    }
    else
    {
      reactionRate *= c;
    }

  }

  return reactionRate;
}


//---------------------------------------------------------
//
// What follows is a set of templated functions that provide 
// concentration dependent pieces of the rate coefficients.
//
//---------------------------------------------------------

template <class ScalarT>
ScalarT Reaction::rxnCoefficient(std::vector<ScalarT> &concentrations,
                                 std::vector<ScalarT> &constants)
{
  
  ScalarT rxnCoeff;

  switch(coefficientType)
  {
  case 0 :
    rxnCoeff = constCoeff;
    break;

  case 1:

    rxnCoeff = complexCoefficient(concentrations, constants);
    break;

  default :
    //unrecognized reaction name
    Report::UserError() << "Unrecognized reaction coefficient type in rection coefficient evaluator";
  }

  return rxnCoeff;
};

template <class ScalarT>
ScalarT Reaction::complexCoefficient(std::vector<ScalarT> &concentrations,
                                     std::vector<ScalarT> &constants)
{

  //The following is done to jibe with Myers and Wampler for GaAs

  ScalarT localCoeff;

  ScalarT ADTEMP = temperature;

  if(chargeProduct == 0)
    if(materialEnum == si || materialEnum == gan)
      localCoeff = latticeConstant;
    else if(materialEnum == gaas)
    {
      ScalarT concsDensity=0.0;
      for(int i=0 ; i<concentrations.size() ; ++i)
        concsDensity += concentrations[i];

      ScalarT rdefTilde = pow(tolerance+C0*concsDensity,1.0/3.0); 
       
      ScalarT rdef = 0.5/rdefTilde;

      localCoeff =  1.0/(1.0/rdef + 1.0/latticeConstant);
    }
    else
    {
      //unrecognized material in complexCoefficient in 
      Report::UserError() << "Unrecognized material in complex reaction coefficient evaluator for neutral charge product";
    }
  else
  {
    //if(material->material=="si")
    if(materialEnum == si || materialEnum == gan)
      localCoeff = -(ScalarT)chargeProduct*unshieldedLength;
    //else if(material->material == "gaas")
    else if(materialEnum == gaas)
    {

      // Defect / Doping Concentration Dependence
      ScalarT concsDensity=0.0;

      charge=1.6021918e-19; //C

      peq = charge/(13.1*8.854214871e-14);  //From Myers & Wampler

      for(int i=0 ; i<concentrations.size() ; ++i)
        concsDensity += concentrations[i];

      ScalarT rdefTilde = pow(tolerance+C0*concsDensity,1.0/3.0); 

      ScalarT rdef = 0.5/rdefTilde;
 
      //Carrier Dependence

      ScalarT TeV = ADTEMP/11604.0;  //temperature in eV

      ScalarT deby0 = sqrt(TeV/peq);

      ScalarT carrierConcentration = tolerance + C0*constants[0];

      if(tolerance + C0*constants[0] < 0.0) carrierConcentration = 1.e-16;

      ScalarT carrier1Debye=deby0/sqrt(carrierConcentration);

      carrierConcentration = tolerance + C0*constants[1];

      if(carrierConcentration < 0.0)carrierConcentration = 1.0;
       
      ScalarT carrier2Debye=deby0/sqrt(carrierConcentration);
       
      ScalarT debye= 1.0/(1.0/carrier1Debye+1.0/carrier2Debye);

      // scrt(z)= 0.65d0*log(1.0d0+z/0.65d0) ! approximates root of screened coulomb potential
         
      ScalarT arg = - (ScalarT)chargeProduct*unshieldedLength/debye;
       
      ScalarT reactionRadiusTilde = 0.65*log(1.0 + arg/0.65)*debye;
       
      //if carriers < 0, default to using rdef

      if(constants[0] < 0.0 || constants[1] < 0)
        reactionRadiusTilde = 1.0e12;
       
      localCoeff =  1.0/(1.0/rdef + 1.0/reactionRadiusTilde);
    }
    else
    {
      //unrecognized material in complexCoefficient in 
      Report::UserError() << "Unrecognized material in complex reaction coefficient evaluator for non-neutral charge product";
    }
  }
 

  ScalarT diffusion = diffusionCoefficient1 + diffusionCoefficient2;

  if(carrierBCIndex >= 0)
    diffusion += sigmaBC*thermalVelocity*C0*constants[carrierBCIndex]*hopLength*hopLength/6.0;

  return localCoeff*diffusion;

}


} // namespace Device
} // namespace Xyce

#endif

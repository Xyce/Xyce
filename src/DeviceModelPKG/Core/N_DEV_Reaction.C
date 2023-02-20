//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

#include <Xyce_config.h>

// Standard includes
#include <vector>
#include <N_UTL_Math.h>
#include <iostream>

// Xyce includes
#include <N_UTL_fwd.h>
#include <N_UTL_FeatureTest.h>
#include <N_DEV_Reaction.h>
#include <N_ERH_ErrorMgr.h>

#include <iostream>


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : Reaction::Reaction
// Purpose       : Default constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
Reaction::Reaction()
  :theRateConstant(0.0),
   numconcs(0),
   numconsts(0),
   carrierEmissionIndex(-1),
   myRateCalc(0),
   material(0),
   FADVectorsAllocated(false),
   carrierBCIndex(-1)
{
  theReactants.resize(0);
  theProducts.resize(0);
  concDependency.resize(0);
  constDependency.resize(0);
}

//-----------------------------------------------------------------------------
// Function      : Reaction::Reaction
// Purpose       : constructor
// Special Notes : Used when vector of reactants, products, and rate constant
//                 are available at time of construction
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
Reaction::Reaction(std::vector< std::pair<int,double> > & reactants,
                   std::vector< std::pair<int,double> > & products,
                   double rateConstant)
  :
  theReactants(reactants),
  theProducts(products),
  theRateConstant(rateConstant),
  numconcs(0),
  numconsts(0),
  carrierEmissionIndex(-1),
  myRateCalc(0),
  C0(1.0),
  energy(0.0),
  material(0),
  myReactionName("UNDEFINED"),
  FADVectorsAllocated(false),
  tolerance(1.e12),
  charge(1.6021918e-19),
  peq(charge/(13.1*8.854214871e-14)),
  constCoeff(1.0),
  unshieldedLength(0.0),
  latticeConstant(0.0),
  Specie1(0),
  Specie2(0),
  chargeProduct(0.0),
  carrierBCIndex(-1),
  carrierCharge(0)
{
  concDependency.resize(0);
  constDependency.resize(0);
}

//-----------------------------------------------------------------------------
// Function      : Reaction::Reaction
// Purpose       : copy constructor
// Special Notes : Needed because we keep STL vectors of these objects, and
//                 we keep a pointer to the rate constant calculator.  Kills
//                 us during all the push_backs without this copy constructor,
//                 as soon as the push starts moving existing objects around.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
Reaction::Reaction(const Reaction &right)
  :theReactants(right.theReactants),
   theProducts(right.theProducts),
   theRateConstant(right.theRateConstant),
   numconcs(right.numconcs),
   numconsts(right.numconsts),
   carrierEmissionIndex(right.carrierEmissionIndex),
   concDependency(right.concDependency),
   constDependency(right.constDependency),
   C0(right.C0),
   energy(right.energy),
   material(right.material),
   myReactionName(right.myReactionName),
   FADVectorsAllocated(false),
   tolerance(right.tolerance),
   charge(right.charge),
   peq(right.peq),
   constCoeff(right.constCoeff),
   unshieldedLength(right.unshieldedLength),
   latticeConstant(right.latticeConstant),
   Specie1(right.Specie1),
   Specie2(right.Specie2),
   chargeProduct(right.chargeProduct),
   carrierBCIndex(right.carrierBCIndex),
   carrierCharge(right.carrierCharge)
{
  if (right.myRateCalc) // Never, ever copy the pointer to the rate calculator!
  {
    myRateCalc = right.myRateCalc->Clone();
  }
  else
  {
    myRateCalc=0;
  }
}


//-----------------------------------------------------------------------------
// Function      : Reaction::operator=
// Purpose       : assignment operator
// Special Notes : Needed because we keep STL vectors of these objects, and
//                 we keep a pointer to the rate constant calculator.  Kills
//                 us during all the push_backs without this copy constructor,
//                 as soon as the push starts moving existing objects around.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
Reaction &
Reaction::operator=(Reaction const & right)
{

  if (this == &right) return *this;    // assignment to self

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "We're doing an assignment of reaction! " << std::endl;
  }

  // otherwise:
  theReactants = right.theReactants;
  theProducts = right.theProducts;
  concDependency = right.concDependency;
  constDependency = right.constDependency;
  theRateConstant = right.theRateConstant;
  numconcs = right.numconcs;
  numconsts = right.numconsts;

  if (right.myRateCalc) // Never, ever copy the pointer to the rate calculator!
  {

    if (myRateCalc) delete myRateCalc; // we already have one, kill it

    myRateCalc=right.myRateCalc->Clone();
  }
  else
  {
    myRateCalc=0;
  }

  return *this;

}

//-----------------------------------------------------------------------------
// Function      : Reaction::Reaction
// Purpose       : constructor
// Special Notes : Used when vector of reactants, products, and rate constant
//                 are available at time of construction
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
Reaction::~Reaction()
{
  if (myRateCalc)
  {
    delete myRateCalc;
    myRateCalc=0;
  }
}

//-----------------------------------------------------------------------------
// Function      : Reaction::setReactants
// Purpose       : copy reactant/stoichiometric coefficients into object
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::setReactants(std::vector< std::pair<int,double> > & reactants)
{
  theReactants=reactants;
}

//-----------------------------------------------------------------------------
// Function      : Reaction::addReactant
// Purpose       : add reactant/stoichiometric coefficients to reactants vector
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::addReactant(int species,double stoich)
{
  std::vector<std::pair<int,double> >::iterator iter;
  std::vector<std::pair<int,double> >::iterator iter_end=theReactants.end();
  bool foundSpecies=false;

  // Make sure we only have each that appears on the LHS of a reaction
  // listed only once.  If the user has specified one twice, combine them
  // into a single term with an augmented stoichiometric coefficient.
  for (iter=theReactants.begin(); iter != iter_end; iter++)
  {
    if (iter->first == species)
    {
      iter->second += stoich;
      foundSpecies=true;
      break;
    }
  }
  // Only if we didn't find the species in our existing list should we
  // add a new entry to the vector.
  if (!foundSpecies)
    theReactants.push_back(std::pair<int,double>(species,stoich));
}


//-----------------------------------------------------------------------------
// Function      : Reaction::addProduct
// Purpose       : add reactant/stoichiometric coefficients to vector
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::addProduct(int species,double stoich)
{
  std::vector<std::pair<int,double> >::iterator iter;
  std::vector<std::pair<int,double> >::iterator iter_end=theProducts.end();
  bool foundSpecies=false;

  // Make sure we only have each that appears on the RHS of a reaction
  // listed only once.  If the user has specified one twice, combine them
  // into a single term with an augmented stoichiometric coefficient.
  for (iter=theProducts.begin(); iter != iter_end; iter++)
  {
    if (iter->first == species)
    {
      iter->second += stoich;
      foundSpecies=true;
      break;
    }
  }
  // Only if we didn't find the species in our existing list should we
  // add a new entry to the vector.
  if (!foundSpecies)
    theProducts.push_back(std::pair<int,double>(species,stoich));
}

//-----------------------------------------------------------------------------
// Function      : Reaction::setProducts
// Purpose       : copy product/stoichiometric coefficients into object
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::setProducts(std::vector< std::pair<int,double> > & products)
{
  theProducts=products;
}

//-----------------------------------------------------------------------------
// Function      : Reaction::getRate
// Purpose       : compute and return the reaction rate.  This is normally the
//                 first step  of getDdt, but we might need this for
//                 other purposes (e.g. a kludge for getting electron/hole
//                 recombination/regeneration rates even when electrons and
//                 holes are constant species)
// Special Notes : It is tacitly assumed that the concentrations and
//                 constants in the provided vector are in an order
//                 consistent with the labeling in the reactants and
//                 products vectors.  The contributions to the species
//                 time derivatives are summed into the provided output
//                 vector because it is assumed that there will be a network
//                 of reactions, not just one.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
double Reaction::getRate(std::vector<double> &concentrations,
                         std::vector<double> &constants)
{
  int rSize=theReactants.size();
  int pSize=theProducts.size();
  double reactionRate;
  int i;
  int species;
  double stoich;
  double c;


  //*************************************************************
  // The following block is an LCM modification
  // I'm calling the veriable coefficient version of getRate 
  // to accommodate variable and material dependent coefficients
  //*************************************************************

  reactionRate=getRateVC(concentrations,constants); 

  return reactionRate; 

  //*************************************************************
  // The preceding block is an LCM modification  -- Code from 
  // here down is dead and will be removed.
  //*************************************************************

  // Reaction rate determined by law of mass action
  reactionRate=theRateConstant;
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

//-----------------------------------------------------------------------------
// Function      : Reaction::getDRateDC
// Purpose       : return a vector of the derivatives of this reaction's
//                 rate with respect to each species
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/13/06
//-----------------------------------------------------------------------------
void Reaction::getDRateDC(std::vector<double> &concentrations,
                          std::vector<double> &constants,
                          std::vector<double> &dratedc)
{
  int cSize=concentrations.size();
  int rSize=theReactants.size();
  int species;
  double stoich,c;
  int i,j;

  if (numconcs != cSize)
  {
    setDependency(cSize);
  }

  // dratedc is the vector of derivatives of the reaction rate w.r.t.
  // concentrations.  Since the reaction rate is just a product of factors
  // involving individual concentrations, we will just be multiplying in
  // terms, so we initialize only those that are necessary to theRateConstant.
  // Only those elements of dratedc that represent concentrations that the
  // rate depends on will be initialized to non-zero.
  for (i=0;i<cSize;++i)
  {
    if (concDependency[i]==1)    // the rxn rate depends on this species
    {
      dratedc[i]=theRateConstant;
    }
  }


  // The reaction rate is K*(product as i=0,rSize){[Xi]^{stoich[i]}}
  // compute dR/dXj
  for (i=0;i<rSize;++i)  // loop over all species on left side of reaction
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


    // for each species in the network, if we depend on that species,
    // update the product in the dR/dXj
    if (stoich != 1.0)
    {
      for (j=0;j<cSize;++j)
      {
        if (concDependency[j] != 0)
        {
          if (j==species)
          {
            dratedc[j] *= stoich*pow(c,stoich-1.0);
          }
          else
          {
            dratedc[j] *= pow(c,stoich);
          }
        }
      }
    }
    else
    {
      for (j=0;j<cSize;++j)
      {
        // when j==species we'd just multiply by 1.0
        if (j!=species && concDependency[j] != 0)
        {
          dratedc[j] *= c;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Reaction::getDRateDConst
// Purpose       : return a vector of the derivatives of this reaction's
//                 rate with respect to each constant species
//
// Special Notes : Very much similar to getDRateDC, but computes
//                 derivative w.r.t. constants.  Like getDRateDC, this
//                 will ultimately be used to compute jacobian elements,
//                 but unlike that routine, it will be used to pick out
//                 specific constants for special treatment, and unlike
//                 getDRateDC, it only computes the derivative with respect
//                 to a single constant, not all constants..
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/17/06
//-----------------------------------------------------------------------------
void Reaction::getDRateDConst(int constNum,
                              std::vector<double> &concentrations,
                              std::vector<double> &constants,
                              double &dratedc)
{
  int cSize=constants.size();
  int rSize=theReactants.size();
  int species;
  double stoich,c;
  int i,j;

  if (numconsts != cSize)
  {
    setConstDependency(cSize);
  }

  if (constDependency[constNum]==0)    // the rxn rate does not depend on this
                                       //species
  {
    dratedc=0.0;
  }
  else
  {
    dratedc=theRateConstant;

    // The reaction rate is K*(product as i=0,rSize){[Xi]^{stoich[i]}}
    // compute dR/dXj
    for (i=0;i<rSize;++i)  // loop over all species on left side of reaction
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


      // for each species in the network, if we depend on that species,
      // update the product in the dR/dXj
      if (stoich != 1.0)
      {
        if (species < 0 && constNum ==-(species+1))  // this specie is the one
          // we're differntiating wrt
        {
          dratedc *= stoich*pow(c,stoich-1.0);
        }
        else
        {
          dratedc *= pow(c,stoich);
        }
      }
      else
      {
        // when constNum is the specie under consideration we'd  multiply by 1.0
        if (constNum !=-(species+1))
        {
          dratedc *= c;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Reaction::getDdt
// Purpose       : sum time derivatives of species into provided vector
// Special Notes : It is tacitly assumed that the concentrations and
//                 constants in the provided vector are in an order
//                 consistent with the labeling in the reactants and
//                 products vectors.  The contributions to the species
//                 time derivatives are summed into the provided output
//                 vector because it is assumed that there will be a network
//                 of reactions, not just one.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::getDdt(std::vector<double> &concentrations,
                      std::vector<double> &constants,
                      std::vector<double> &ddt)
{
  int rSize=theReactants.size();
  int pSize=theProducts.size();
  double reactionRate;
  int i;
  int species;
  double stoich;

  // Reaction rate determined by law of mass action
  if(carrierEmissionIndex < 0)
    {
      //reactionRate=getRate(concentrations,constants);
      reactionRate=getRateVC(concentrations,constants); 
    }
  else
    {
      reactionRate=getFDEmissionRate(concentrations,constants);
    }

  // Now update time derivatives:

  for (i=0;i<rSize;++i)
  {
    species=theReactants[i].first;
    stoich=theReactants[i].second;
    if (species>=0)
    {
      ddt[species] -= stoich*reactionRate;
    }
  }

  for (i=0;i<pSize;++i)
  {
    species=theProducts[i].first;
    stoich=theProducts[i].second;
    if (species>=0)
    {
      ddt[species] += stoich*reactionRate;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Reaction::getJac
// Purpose       : sum derivatives of time derivatives of species with
//                 respect to concentration into provided matrix
// Special Notes : It is tacitly assumed that the concentrations and
//                 constants in the provided vector are in an order
//                 consistent with the labeling in the reactants and
//                 products vectors.  The contributions to the species
//                 second derivatives are summed into the provided output
//                 vector because it is assumed that there will be a network
//                 of reactions, not just one.
//
//                 An implicit assumption here is that the species numbering
//                 does not change through the run --- that is, once we
//                 call this routine with a vector of concentrations, it
//                 is assumed from that point forward that element i of the
//                 vector always refers to the same species in subsequent
//                 calls.  If it is necessary to change that over time, the
//                 reaction network will probably have to be reconstructed
//                 and re-initialized.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::getJac(std::vector<double> &concentrations,
                      std::vector<double> &constants,
                      std::vector<std::vector<double> > &jac)
{
  int cSize=concentrations.size();
  int constSize=constants.size();
  int rSize=theReactants.size();
  int pSize=theProducts.size();
  int species;
  double stoich,c;
  std::vector<double> dratedc(cSize,0.0);
  int i,j;
  int FADSize=cSize+constSize;

  // here is where the assumption of static species numbering is coded.
  // Should it be necessary to change that assumption, here's where the work
  // would be needed.
  //LCM --  This must be modified to use the FADed versions with Fermi-Dirac statistics

  //allocate FAD arrays once
  if(!FADVectorsAllocated)
    {
      defects.resize(cSize);
      carriers.resize(constSize);
      FADVectorsAllocated = true;
    }

  setDependency(cSize);
  if(carrierEmissionIndex <0)
    {
      //getDRateDC(concentrations, constants, dratedc); 
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      FDFadType tempdRdCFD;
      for( i=0 ; i<2 ; ++i)
        {
          carriers[i] = constants[i];
          //carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          defects[i].diff(i+2,FADSize);
        }

      tempdRdCFD = getRateVC(defects, carriers);
      for(int j=0 ; j<cSize ; ++j)
        {
          dratedc[j] = tempdRdCFD.dx(j+2);
        }
    }
  else
    {
      if (numconcs != cSize)
        {
          setDependency(cSize);
        }
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      FDFadType tempdRdCFD;
      for( i=0 ; i<2 ; ++i)
        {
          carriers[i] = constants[i];
          //carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          defects[i].diff(i+2,FADSize);
        }

      tempdRdCFD = getFDEmissionRate(defects, carriers);
      for(int j=0 ; j<cSize ; ++j)
        {
          dratedc[j] = tempdRdCFD.dx(j+2);
        }
    }


  /*
  if(carrierEmissionIndex <0)
    {
      getDRateDC(concentrations, constants, dratedc); 
    }
  else
    {
      if (numconcs != cSize)
        {
          setDependency(cSize);
        }
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      FDFadType tempdRdCFD;
      std::vector<FDFadType> defects(cSize);
      std::vector<FDFadType> carriers(2);
      for( i=0 ; i<2 ; ++i)
        {
          carriers[i] = constants[i];
          carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          defects[i].diff(i+2,FADSize);
        }

      tempdRdCFD = getFDEmissionRate(defects, carriers);
      for(int j=0 ; j<cSize ; ++j)
        {
          dratedc[j] = tempdRdCFD.dx(j+2);
        }
    }
  */
  // we now know how the reaction rate depends on the various concentrations.
  // Now we can assemble the jacobian, which is d^2[Xi]/(dt d[Xj])


  // species on the left of the reaction
  for (i=0;i<rSize;++i)
  {
    species=theReactants[i].first;
    stoich=theReactants[i].second;
    if (species>=0)
    {
      for (j=0;j<cSize;++j)
      {
        if (concDependency[j]!=0)
        {
          jac[species][j] -= stoich*dratedc[j];
        }
      }
    }
  }
  // species on the right
  for (i=0;i<pSize;++i)
  {
    species=theProducts[i].first;
    stoich=theProducts[i].second;
    if (species>=0)
    {
      for (j=0;j<cSize;++j)
      {
        if (concDependency[j]!=0)
        {
          jac[species][j] += stoich*dratedc[j];
        }
      }
    }
  }

}


//-----------------------------------------------------------------------------
// Function      : Reaction::getJac
// Purpose       : sum derivatives of time derivatives of species with
//                 respect to concentration into provided matrix
//
// Scope         : public
// Creator       : Lawrence C Musson, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/3/14
//-----------------------------------------------------------------------------
  void Reaction::getJacobianVC(std::vector<double> &concentrations,
                               std::vector<double> &constants,
                               std::vector<std::vector<double> > &jac,
                               std::vector<double> &constVec)
{
  int constSize=constants.size();
  int cSize=concentrations.size();
  int rSize=theReactants.size();
  int pSize=theProducts.size();
  int species;
  double stoich,c;
  std::vector<double> dratedc(cSize,0.0);
  int i,j;
  std::vector<double> dratedConst(constSize,0.0);
  int FADSize=cSize+constSize;


  // here is where the assumption of static species numbering is coded.
  // Should it be necessary to change that assumption, here's where the work
  // would be needed.

  //allocate FAD arrays once
  if(!FADVectorsAllocated)
    {
      defects.resize(cSize);
      carriers.resize(constSize);
      FADVectorsAllocated = true;
    }

  setDependency(cSize);
  if (numconsts != constSize)
    {
      setConstDependency(constSize);
    }


  if(carrierEmissionIndex <0)
    {
      //getDRateDC(concentrations, constants, dratedc); 
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      FDFadType tempdRdCFD;
      for( i=0 ; i<constSize ; ++i)
        {
          carriers[i] = constants[i];
          carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          defects[i].diff(i+constSize,FADSize);
        }

      tempdRdCFD = getRateVC(defects, carriers);
      for(int j=0 ; j<cSize ; ++j)
        {
          dratedc[j] = tempdRdCFD.dx(j+constSize);
        }
      for(int j=0 ; j<constSize ; ++j)
        dratedConst[j] = tempdRdCFD.dx(j);

    }
  else
    {
      if (numconcs != cSize)
        {
          setDependency(cSize);
        }
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      FDFadType tempdRdCFD;
      for( i=0 ; i<constSize ; ++i)
        {
          carriers[i] = constants[i];
          carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          defects[i].diff(i+constSize,FADSize);
        }

      tempdRdCFD = getFDEmissionRate(defects, carriers);
      for(int j=0 ; j<cSize ; ++j)
        {
          dratedc[j] = tempdRdCFD.dx(j+constSize);
        }
      for(int j=0 ; j<constSize ; ++j)
        dratedConst[j] = tempdRdCFD.dx(j);
    }


  // we now know how the reaction rate depends on the various concentrations.
  // Now we can assemble the jacobian, which is d^2[Xi]/(dt d[Xj])


  // species on the left of the reaction
  for (i=0;i<rSize;++i)
  {
    species=theReactants[i].first;
    stoich=theReactants[i].second;
    if (species>=0)
    {
      for (j=0;j<cSize;++j)
      {
          jac[species][j] -= stoich*dratedc[j];
      }
    }
  }
  // species on the right
  for (i=0;i<pSize;++i)
  {
    species=theProducts[i].first;
    stoich=theProducts[i].second;
    if (species>=0)
    {
      for (j=0;j<cSize;++j)
      {
          jac[species][j] += stoich*dratedc[j];
      }
    }
  }

    // we now know how the reaction rate depends on the various constants
    // we can assemble the derivatives we need.

  for(int constIdx=0 ; constIdx<constSize ; ++constIdx)
    {
      int offset = constIdx*cSize;
      // species on the left of the reaction
      for (i=0;i<rSize;++i)
        {
          species=theReactants[i].first;
          stoich=theReactants[i].second;
          if (species>=0)
            constVec[species+offset] -= stoich*dratedConst[constIdx];
        }
      // species on the right
      for (i=0;i<pSize;++i)
        {
          species=theProducts[i].first;
          stoich=theProducts[i].second;
          if (species>=0)
            constVec[species+offset] += stoich*dratedConst[constIdx];
        }
    }

}

//-----------------------------------------------------------------------------
// Function      : Reaction::getDFdConst
// Purpose       : sum derivatives of time derivatives of species with
//                 respect to concentration into provided vector
// Special Notes : It is tacitly assumed that the concentrations and
//                 constants in the provided vector are in an order
//                 consistent with the labeling in the reactants and
//                 products vectors.  The contributions to the species
//                 second derivatives are summed into the provided output
//                 vector because it is assumed that there will be a network
//                 of reactions, not just one.
//
//                 An implicit assumption here is that the species numbering
//                 does not change through the run --- that is, once we
//                 call this routine with a vector of concentrations, it
//                 is assumed from that point forward that element i of the
//                 vector always refers to the same species in subsequent
//                 calls.  If it is necessary to change that over time, the
//                 reaction network will probably have to be reconstructed
//                 and re-initialized.
//
//                 This method sums into a vector containing the derivatives
//                 of variable species with respect to the specified constant
//                 specie.  It'll be used as a fragment of a larger jacobian.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::getDFdConst(int constantNumber,
                                      std::vector<double> &concentrations,
                                      std::vector<double> &constants,
                                      std::vector<double> &dFdConst)
{
  int constSize=constants.size();
  int cSize=concentrations.size();
  int rSize=theReactants.size();
  int pSize=theProducts.size();
  int species;
  double stoich,c;
  double dratedc;
  int i;
  int FADSize=cSize+constSize;

  // here is where the assumption of static species numbering is coded.
  // Should it be necessary to change that assumption, here's where the work
  // would be needed.

  //allocate FAD arrays once
  if(!FADVectorsAllocated)
    {
      defects.resize(cSize);
      carriers.resize(constSize);
      FADVectorsAllocated = true;
    }

  if(carrierEmissionIndex <0)
    {
      //getDRateDConst(constantNumber,concentrations, constants, dratedc);
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns

      int cSize = concentrations.size();
      if (numconsts != constSize)
        {
          setConstDependency(constSize);
        }

      FDFadType tempdRdCFD;
      for( i=0 ; i<2 ; ++i)
        {
          carriers[i] = constants[i];
          carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<concentrations.size() ; ++i)
        {
          defects[i] = concentrations[i];
          //defects[i].diff(i+2,FADSize);
        }

      tempdRdCFD = getRateVC(defects, carriers);
      dratedc = tempdRdCFD.dx(constantNumber);
    }
  else
    {
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      int cSize = concentrations.size();
      if (numconsts != constSize)
        {
          setConstDependency(constSize);
        }

      FDFadType tempdRdCFD;
      for( i=0 ; i<2 ; ++i)
        {
          carriers[i] = constants[i];
          carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          //defects[i].diff(i+2,FADSize);
        }
      tempdRdCFD = getFDEmissionRate(defects, carriers);
      dratedc = tempdRdCFD.dx(constantNumber);
    }

    /*
  if(carrierEmissionIndex <0)
    {
      getDRateDConst(constantNumber,concentrations, constants, dratedc);
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
    }
  else
    {
      //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
      int cSize = concentrations.size();
      if (numconsts != constSize)
        {
          setConstDependency(constSize);
        }

      FDFadType tempdRdCFD;
      std::vector<FDFadType> defects(cSize);
      std::vector<FDFadType> carriers(2);
      for( i=0 ; i<2 ; ++i)
        {
          carriers[i] = constants[i];
          carriers[i].diff(i,FADSize);
        }
      for( i=0 ; i<cSize ; ++i)
        {
          defects[i] = concentrations[i];
          defects[i].diff(i+2,FADSize);
        }
      tempdRdCFD = getFDEmissionRate(defects, carriers);
      dratedc = tempdRdCFD.dx(constantNumber);
    }
  */

  // no point going through these loops if this reaction doesn't depend
  // on this constant!
  if (constDependency[constantNumber] != 0)
  {
    // we now know how the reaction rate depends on the various constants
    // we can assemble the derivatives we need.

    // species on the left of the reaction
    for (i=0;i<rSize;++i)
    {
      species=theReactants[i].first;
      stoich=theReactants[i].second;
      if (species>=0)
      {
        dFdConst[species] -= stoich*dratedc;
      }
    }
    // species on the right
    for (i=0;i<pSize;++i)
    {
      species=theProducts[i].first;
      stoich=theProducts[i].second;
      if (species>=0)
      {
        dFdConst[species] += stoich*dratedc;
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Reaction::setSimpleRateCalculator
// Purpose       : set the rate calculator for this reaction to a simple one
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::setSimpleRateCalculator(double k_in, double C0, double t0,
                                             double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "simplerate";

  myRateCalc = dynamic_cast<RateCalculator *> (new SimpleRateCalculator(k_in, C0, t0, x0));
}

//-----------------------------------------------------------------------------
// Function      : Reaction::setCaptureRateCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 capture rate type (A0+E->AM, for example)
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::setCaptureRateCalculator(double sigma, double v,
                                              double C0, double t0, double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "capturerate";

  //N.B., it's labeled v (velocity), but really it's the charge of the carrier

  if(v < 0.0)
    carrierCharge = -1;
  else if(v > 0.0)
    carrierCharge =  1;
  else
    {
      Report::UserError() <<"Illegal carrier charge";
    }

  myRateCalc = dynamic_cast<RateCalculator *>(new CaptureRateCalculator(sigma, v, C0, t0, x0));
}


//-----------------------------------------------------------------------------
// Function      : Reaction::setEmissionRateCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 emission rate type (AM->A0+E, for example)
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::setEmissionRateCalculator(double sigma, double v,
                                               double N, double Energy,
                                               double C0, double t0, double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "emissionrate";

  //N.B., it's labeled v (velocity), but really it's the charge of the carrier

  if(v < 0.0)
    carrierCharge = -1;
  else if(v > 0.0)
    carrierCharge =  1;
  else
    {
      Report::UserError() << "Illegal carrier charge";
    }


  myRateCalc = dynamic_cast<RateCalculator *> (
    new EmissionRateCalculator(sigma, v, N, Energy, C0, t0, x0));
}

//-----------------------------------------------------------------------------
// Function      : Reaction::setFDEmissionRateCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 emission rate type (AM->A0+E, for example).  This
//                 include Fermi-Dirac stats
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 06/30/2014
//-----------------------------------------------------------------------------
  void Reaction::setFDEmissionRateCalculator(int carrierIndex, double sigma, double Energy, 
                                             double v, double C0i, double t0, double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  //N.B., it's labeled v (velocity), but really it's the charge of the carrier

  if(v < 0.0)
    carrierCharge = -1;
  else if(v > 0.0)
    carrierCharge =  1;
  else
    {
      Report::UserError() << "Illegal carrier charge";
    }

  carrierEmissionIndex = carrierIndex;  //emitted carrier
  C0 = C0i;  //concentration
  energy = Energy; //emission energy

  myReactionName = "fdemissionrate";

  myRateCalc = dynamic_cast<RateCalculator *> (
     new FDEmissionRateCalculator(sigma, Energy, v, C0, t0, x0));
}


//-----------------------------------------------------------------------------
// Function      : Reaction::setComplexRateCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 complex rate type (A+B->AB, for example)
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::setComplexRateCalculator(std::vector<Specie> &V,
                                              std::vector<Specie> &C,
                                              double C0, double t0, double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "complexrate";

  //We need to store the reactant list for convenient access

  myRateCalc = dynamic_cast<RateCalculator *> (
    new ComplexRateCalculator(V,C,theReactants,C0, t0, x0));

  //The sanity check on this was already handled in the constructor for myRateCalc

  if (theReactants[0].first >= 0)
    Specie1 = &(V[theReactants[0].first]);
  else
    Specie1 = &(C[-(theReactants[0].first+1)]);

  // Handle case where there's only one species with a coefficient of 2.0
  if (theReactants.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
  else
    {
      if (theReactants[1].first >= 0)
        Specie2 = &(V[theReactants[1].first]);
      else
        Specie2 = &(C[-(theReactants[1].first+1)]);
    }

  chargeProduct    =Specie1->getChargeState();
  chargeProduct   *=Specie2->getChargeState();

}

//-----------------------------------------------------------------------------
// Function      : Reaction::setComplexRateMultiplierCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 complex rate type (A+B->AB, for example) times
//               : the specified multiplier
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 06/26/2014
//-----------------------------------------------------------------------------
void Reaction::setComplexMultiplierRateCalculator(std::vector<Specie> &V,
                                                  std::vector<Specie> &C,
                                                  double C0, double t0, double x0,
                                                  double multiplier)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "complexmultiplierrate";

  myRateCalc = dynamic_cast<RateCalculator *> (
    new ComplexRateCalculator(V,C,theReactants,C0, t0, x0));


  //The sanity check on this was already handled in the constructor for myRateCalc

  if (theReactants[0].first >= 0)
    Specie1 = &(V[theReactants[0].first]);
  else
    Specie1 = &(C[-(theReactants[0].first+1)]);

  // Handle case where there's only one species with a coefficient of 2.0
  if (theReactants.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
  else
    {
      if (theReactants[1].first >= 0)
        Specie2 = &(V[theReactants[1].first]);
      else
        Specie2 = &(C[-(theReactants[1].first+1)]);
    }

  chargeProduct    =Specie1->getChargeState();
  chargeProduct   *=Specie2->getChargeState();
}

//-----------------------------------------------------------------------------
// Function      : Reaction::setDecomplexRateCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 decomplex rate type (AB->A+B, for example)
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::setDecomplexRateCalculator(std::vector<Specie> &V,
                                              std::vector<Specie> &C,
                                              double bindingEnergy,
                                              double gammaAB, double gammaA,
                                              double gammaB, double concSi,
                                              double C0, double t0, double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "decomplexrate";

  myRateCalc = dynamic_cast<RateCalculator *> (
     new DecomplexRateCalculator(V,C,theReactants,theProducts,
                                              bindingEnergy, gammaAB,gammaA,
                                              gammaB, concSi,
                                              C0, t0, x0));

  //The sanity check on this was already handled in the constructor for myRateCalc

  if (theProducts[0].first >= 0)
    Specie1 = &(V[theProducts[0].first]);
  else
    Specie1 = &(C[-(theProducts[0].first+1)]);

  // Handle case where there's only one species with a coefficient of 2.0
  if (theProducts.size() == 1)
    {
      Specie2 = Specie1;   // that way we can just treat as A+A instead of 2A
    }
  else
    {
      if (theProducts[1].first >= 0)
        Specie2 = &(V[theProducts[1].first]);
      else
        Specie2 = &(C[-(theProducts[1].first+1)]);
    }

  chargeProduct    =Specie1->getChargeState();
  chargeProduct   *=Specie2->getChargeState();

}


//-----------------------------------------------------------------------------
// Function      : Reaction::setDecomplexRateCalculator
// Purpose       : set the rate calculator for this reaction to the
//                 decomplex rate type (AB->A+B, for example)
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 03/24/2014
//-----------------------------------------------------------------------------
void Reaction::setBourgoinCorbettHoleRateCalculator(std::vector<Specie> &V,
                                                    std::vector<Specie> &C, double sigma,
                                                    double C0, double t0, double x0)
{
  if (myRateCalc) // we already have a calculator set
  {
    delete myRateCalc; // kill it
    myRateCalc=0;
  }

  myReactionName = "bouroincorbettholerate";

  myRateCalc = dynamic_cast<RateCalculator *> (
     new BourgoinCorbettHoleRateCalculator(V,C,theReactants,theProducts, sigma,
                                           C0, t0, x0));
}

//-----------------------------------------------------------------------------
// Function      : Reaction::setRateConstantFromCalculator
// Purpose       : set the rate constant from the temperature using our rate
//                 calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::setRateConstantFromCalculator(double T)
{
  if (myRateCalc) // if we *have* a rate calculator
  {
    if(material->material=="UNDEFINED" || material->material=="")
      {
        Report::UserError() << "The bulk material " << material->material << " is unknown";
      }

    if(!myRateCalc->isBulkMaterialSet())
      myRateCalc->setBulkMaterial(material->material);

    setRateConstant(myRateCalc->computeRateConstant(T));
  }
  // otherwise just leave it where it is (which might be 0)
}


//-----------------------------------------------------------------------------
// Function      : Reaction::setRateConstantFromCalculator
// Purpose       : set the rate constant from the temperature using our rate
//                 calculator
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
void Reaction::setRateConstantFromCalculator(double T, 
                                     std::vector<double> &concs,
                                     std::vector<double> &constant_vec)
{
  if (myRateCalc) // if we *have* a rate calculator
    {
      if(material->material=="UNDEFINED" || material->material=="")
        {
          Report::UserError() << "The bulk material " << material->material << " is unknown";
        }
    setRateConstant(myRateCalc->computeRateConstant(T,concs,constant_vec));

  }
  // otherwise just leave it where it is (which might be 0)
}

//-----------------------------------------------------------------------------
// Function      : Reaction::scaleRateConstantFromCalculator()
// Purpose       : scale the rate constant using values stored in our rate
//                 calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::scaleRateConstantFromCalculator()
{
  if (myRateCalc) // if we *have* a rate calculator
  {
    theRateConstant *= myRateCalc->rateConstantScaleFactor();
  }
  // otherwise just leave it where it is
}

//-----------------------------------------------------------------------------
// Function      : Reaction::unscaleRateConstantFromCalculator()
// Purpose       : unscale the rate constant using values stored in our rate
//                 calculator
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/01/06
//-----------------------------------------------------------------------------
void Reaction::unscaleRateConstantFromCalculator()
{
  if (myRateCalc) // if we *have* a rate calculator
  {
    theRateConstant /= myRateCalc->rateConstantScaleFactor();
  }
  // otherwise just leave it where it is
}


//-----------------------------------------------------------------------------
// Function      : Reaction::setDependency
// Purpose       : Set the vector that denotes which of the cSize species
//                 that the network tracks are actually used in this one
//                 reaction as reactants, that is, the species that are
//                 used in computing the reaction rate.
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
void Reaction::setDependency(int cSize)
{
  int rSize=theReactants.size();
  int i;
  numconcs=cSize;
  concDependency.resize(numconcs,0);

  for (i=0;i<rSize;++i)
  {
    if (theReactants[i].first >= 0)
    {
      concDependency[theReactants[i].first]=1;
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Reaction::setConstDependency
// Purpose       : Set the vector that denotes which of the cSize constant
//                 species that the network contains are actually used in this
//                 one reaction as reactants, that is, the species that are
//                 used in computing the reaction rate.
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/17/06
//-----------------------------------------------------------------------------
void Reaction::setConstDependency(int cSize)
{
  int rSize=theReactants.size();
  int i;
  numconsts=cSize;
  constDependency.resize(numconsts,0);

  for (i=0;i<rSize;++i)
  {
    if (theReactants[i].first < 0)
    {
      constDependency[-(theReactants[i].first+1)]=1;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Reaction::output
// Purpose       : Accessor function to set rate constant
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 5/27/06
//-----------------------------------------------------------------------------
void Reaction::output
  ( const std::vector<Specie> & species, std::ostream & os ) const
{
  int i=0;
  int isize = theReactants.size();
  os << "   Rxn: ";
  bool firstPrintDone=false;
  for (i=0;i<isize;++i)
  {
    int speciesIndex = theReactants[i].first;
    if (speciesIndex >= 0)
    {
      if (firstPrintDone) os << " + ";
      double tmp = theReactants[i].second;
      if (tmp > 1.0)
        os << " "<<tmp<< " * ";
      os.setf(std::ios::right); os.width(3);
      os << species[speciesIndex].getName();

      firstPrintDone = true;
    }
  }

  os << " = ";
  isize = theProducts.size();
  firstPrintDone=false;
  for (i=0;i<isize;++i)
  {
    int speciesIndex = theProducts[i].first;
    if (speciesIndex >= 0)
    {
      if (firstPrintDone) os << " + ";
      double tmp = theProducts[i].second;
      if (tmp > 1.0)
        os << " " << tmp << " * ";
      os.setf(std::ios::right); os.width(3);

      os << species[speciesIndex].getName();
      firstPrintDone = true;
    }
  }

  os << "    Rate Constant: ";
  os.precision(8); os.setf(std::ios::scientific);
  os << theRateConstant;

  os << std::endl;

}


//-----------------------------------------------------------------------------
// Function      : Reaction::setMaterial
// Purpose       : Set the bulk material in the reaction object
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 5/27/06
//-----------------------------------------------------------------------------

  void Reaction::setMaterial(MaterialLayer *mat, double Temp)
  {
    material = mat;

    ExtendedString mN = material->material;
    mN.toLower();
    material->material = mN;

    if(material->material == "si")
      materialEnum = si;
    else if(material->material == "gaas")
      materialEnum = gaas;
    else if(material->material == "gan")
      materialEnum = gan;

    //now with the material known, we can set some 
    //material dependent quantities and arrange for 
    //other, concentration dependent coefficients
    //to be set up.

    if (myRateCalc)
      myRateCalc->computeRateConstant(Temp);
  }



//-----------------------------------------------------------------------------
// Function      : Reaction::setCoefficient
// Purpose       : Set the coefficient for the reactions
// Special Notes : This is different than the rate constant.  It's most often
//                 just 1.0.  But for some complexing reactions, there are
//                 concentration dependent values to be used.
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 5/27/06
//-----------------------------------------------------------------------------

  void Reaction::setCoefficient(double Temp)
  {

    //now with the material and reaction type known, 
    //we can arrange for other, concentration 
    //dependent coefficients to be set up.

    if(myReactionName == "simplerate")
      {
        constCoeff = 1.0;
        coefficientType = 0;
      }
    else if(myReactionName == "capturerate")
      {
        //For this, the coefficient is going to be carrier velocity
        if(carrierCharge == 1)
          constCoeff = material->holeThermalV;
        else
          constCoeff = material->electronThermalV;
        coefficientType = 0;
      }
    else if(myReactionName == "emissionrate")
      {
        //For this, the coefficient is going to be the product
        //of carrier velocity and DOS
        if(carrierCharge == 1)
          {
            constCoeff = material->holeThermalV;
            constCoeff *= material->Nv;
         }
        else
          {
            constCoeff = material->electronThermalV;
            constCoeff *= material->Nc;
          }
        coefficientType = 0;
      }
    else if(myReactionName == "complexrate" || myReactionName == "complexmultiplierrate"
            || myReactionName == "decomplexrate") 
      {
        tolerance = 1.e12;
        unshieldedLength = material->defectReactionRadius;
        latticeConstant = material->latticeConstant;
        coefficientType = 1;

        diffusionCoefficient1 = Specie1->getDiffusionCoefficient(Temp);
        diffusionCoefficient2 = Specie2->getDiffusionCoefficient(Temp);

        //Need to check if Bourogin Corbett enhanced diffusion is turned on
        if(Specie1->getEnhancedDiffusion())
          {
            carrierBCIndex = Specie1->getBCCarrierIndex();
            hopLength = Specie1->getBCHopLength();
            sigmaBC = Specie1->getBCSigma();
            if(carrierBCIndex == 0)
              thermalVelocity = material->electronThermalV;
            else
              thermalVelocity = material->holeThermalV;
          }
        else if(Specie2->getEnhancedDiffusion())
          {
            carrierBCIndex = Specie2->getBCCarrierIndex();
            hopLength = Specie2->getBCHopLength();
            sigmaBC = Specie2->getBCSigma();
            if(carrierBCIndex == 0)
              thermalVelocity = material->electronThermalV;
            else
              thermalVelocity = material->holeThermalV;
          }
        else
          {
            carrierBCIndex = -1;
            hopLength = 0.0;
            thermalVelocity = 0.0;
          }
      }
    else if(myReactionName == "fdemissionrate")
      {
        //For this, the coefficient is going to be the product
        //of carrier velocity and DOS
        if(carrierCharge == 1)
          {
            constCoeff = material->holeThermalV;
          }
        else
          {
            constCoeff = material->electronThermalV;
          }
        coefficientType = 0;
      }
    else
      {
        //unrecognized reaction name
        Report::UserError() << "Unrecognized reaction name in setting rate coefficient creation " << myReactionName;
      }

    return;

  }

} // namespace Device
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 03/20/2006
//
//
//
//
//-----------------------------------------------------------------------------
#ifndef N_DEV_ReactionNetwork_H
#define N_DEV_ReactionNetwork_H

#include <N_UTL_fwd.h>
#include <iosfwd>
#include <vector>
#include <map>
#include <string>

#include <N_DEV_Reaction.h>
#include <N_DEV_Specie.h>
#include <N_DEV_MaterialLayer.h>

#include <N_ERH_ErrorMgr.h>
// ---------- Forward Declarations ----------

// Note: sstream and strstream are only needed here because of all the
// inlined functions.
#include <sstream>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_ReactionNetwork
// Purpose       :
// Special Notes :
// Creator       : Tom Russo, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
class ReactionNetwork
{
public:

  ReactionNetwork(
      const SolverState & solver_state,
      const std::string &name = "NoName");

  ReactionNetwork(const ReactionNetwork & right);
  virtual ~ReactionNetwork();

  void setReactionNetworkFromFile(const NetlistLocation &netlist_location, const std::string &fileName);

  void addReaction(const std::string &name);
  void addReaction(const std::string &name, Reaction &reaction);
  void addReactant(const std::string &name, const std::string &reactant, double stoich);
  void addProduct(const std::string &name, const std::string &reactant, double stoich);

  // set rate constant calculator for each type of reaction
  void setSimpleCalc(const std::string &name, double k);
  void setCaptureCalc(const std::string &name, double sigma, double v);
  void setEmissionCalc(const std::string &name, double sigma, double v, double N, double E);
  void setElectronCaptureCalc(const std::string &name, double sigma);
  void setElectronEmissionCalc(const std::string &name, double sigma, double E);
  void setFDElectronEmissionCalc(const std::string &name, double sigma, double E);
  void setHoleCaptureCalc(const std::string &name, double sigma);
  void setHoleEmissionCalc(const std::string &name, double sigma, double E);
  void setFDHoleEmissionCalc(const std::string &name, double sigma, double E);
  void setComplexCalc(const std::string &name);
  void setComplexMultiplierCalc(const std::string &name, double multiplier);
  void setDecomplexCalc(const std::string &name, double bindingEnergy,
                        double gammaAB, double gammaA, double gammaB, double concSi);
  void setBourgoinCorbettHoleCalc(const std::string &name,double sigma);
  void setBourgoinCorbettCalc(const std::string &speciesName,const std::string &carrierName,
                              double sigma, double hopDistance);

  void setRateConstant(const std::string &name, double k);
  void scaleRateConstant(const std::string &name, double kscale);

  void setScaleParams(double c, double t, double x);
  void setMaterial(MaterialLayer *material, double Temp);
  void setCoefficients(double Temp);
  void setRxnVariableCoeffs(bool variableCoeffs);

  // Use rate constant calculators
  void scaleRateConstantFromCalculator(const std::string &name);
  void unscaleRateConstantFromCalculator(const std::string &name);
  void setRateConstantFromCalculator(const std::string &name, double T);
  void setRateConstantsFromCalc(double T);
  void setRateConstantFromCalculator(const std::string &name, double T, 
                                     std::vector<double> &concs,
                                     std::vector<double> &constant_vec);
  void setRateConstantsFromCalc(double T, 
                                std::vector<double> &concs,
                                std::vector<double> &constant_vec);
  void scaleRateConstantsFromCalc();
  void unscaleRateConstantsFromCalc();

  void setSpecies(std::vector<Specie> &theSpeciesVect);
  void addSpecie(const Specie &aSpecie);
  void setConstants(std::vector<Specie> &theConstantsVect);
  void addConstant(const Specie &aConstant);

  int getReactionNum(const std::string name);

  void addSourceTerm(const std::string &speciesName, const std::string &expressionStr);
  void addSourceTerm(const std::string &speciesName,Util::Expression *expression);
  void addMasterSourceTerm(const std::string &speciesName);

  void addInitialCondition(const std::string &speciesName, double value);
  std::pair<std::string,double> getInitialCondition(int i);
  int getNumInitialConditions();

#if 0
  void setSimTime(double time);
  void setSimDT(double step);
  double getBreakpointTime();
#else
  bool getBreakPoints(std::vector<Util::BreakPoint> & breakPointTimes);
#endif
  inline void setSourceScaleFac(double scf) {sourceScaleFac=scf;};
  inline void setMasterSourceValue(double msv) {masterSourceValue=msv;};

  void getDdt(std::vector<double> &concs,std::vector<double> &constants,
              std::vector<double> &ddt);
  void getJac(std::vector<double> &concs, std::vector<double> &constants,
              std::vector<std::vector<double> >&jac);
  void getDFdConst(const std::string &constantName,
                   std::vector<double> &concs, std::vector<double> &constants,
                   std::vector<double> &dFdConst);

  void getJacobianVC(std::vector<double> &concs, std::vector<double> &constants,
                     std::vector<std::vector<double> >&jac, std::vector<double> &constVec);

  double getRate(std::vector<double> &concs,std::vector<double> &constants,
                 std::vector<int> &captureVect, std::vector<int> &emissionVect);
  void getDRateDC(std::vector<double> &concs,std::vector<double> &constants,
                  std::vector<int> &captureVect, std::vector<int> &emissionVect,
                  std::vector<double>&dratedc);

  void getDRateDConst(std::vector<double> &concs,std::vector<double> &constants,
                      std::vector<int> &captureVect, std::vector<int> &emissionVect,
                      std::vector<double>&dratedc);

  double getCaptureLifetime(std::vector<double> &concs,std::vector<double> &constants,
                            std::vector<int> &captureVect,double &concentration);
  void getCaptureLifetimes(std::vector<double> &concs,std::vector<double> &constants,
                           std::vector<int> &captureVect,double &concentration,
                           std::vector<double> &lifetimes);
  double getELifetime(std::vector<double>&concs,std::vector<double>&constants);
  double getHLifetime(std::vector<double>&concs,std::vector<double>&constants);
  void getELifetimes(std::vector<double>&concs,std::vector<double>&constants,
                     std::vector<double> &lifetimes);
  void getHLifetimes(std::vector<double>&concs,std::vector<double>&constants,
                     std::vector<double> &lifetimes);

  double getERate(std::vector<double> &concs,std::vector<double> &constants);
  double getHRate(std::vector<double> &concs,std::vector<double> &constants);

  void getDERateDC(std::vector<double> &concs,std::vector<double> &constants,
                   std::vector<double>&dratedc);
  void getDHRateDC(std::vector<double> &concs,std::vector<double> &constants,
                   std::vector<double>&dratedc);

  void getDERateDConst(std::vector<double> &concs,std::vector<double> &constants,
                       std::vector<double>&dratedConst);
  void getDHRateDConst(std::vector<double> &concs,std::vector<double> &constants,
                       std::vector<double>&dratedConst);

  int getNumSpecies();
  int getNumConstants();
  const std::string & getSpeciesName(int i);
  const std::string & getConstantsName(int i);
  int getSpeciesNum(const std::string &name);
  int getConstantNum(const std::string &name);
  int getReactantNum(const std::string &name);
  bool reactantExist(const std::string &name);
  bool constantExist(const std::string &name);

  void setName(const std::string &name);
  void clear();
  void output(std::ostream & os) const;

  double getDiffusionCoefficient (const std::string & name, const double temp);
  double getDiffusionCoefficient (int specie, const double temp);
  double getDiffusionCoefficient (int specie, const double temp, 
                                  std::vector<double> concs, 
                                  std::vector<double> carrierConcs);

  int getChargeState (const std::string & name);
  int getChargeState (int specie);

  void setApplySources(bool flag);

private:
  Reaction &getReaction(const std::string &name);
  Reaction &getReaction(int i);

  std::map<std::string,int> speciesMap;
  std::vector<Specie> species;
  std::map<std::string,int> constantsMap;
  std::vector<Specie> constants;
  std::vector<std::pair<std::string,double> > initialConditions;
  std::vector<Reaction> theReactions;
  std::map<std::string,int> reactionNamesMap;
  std::vector<std::string> reactionNames;
  std::string myName;
  std::vector<int> electronCaptureReactions;
  std::vector<int> holeCaptureReactions;
  std::vector<int> electronEmissionReactions;
  std::vector<int> holeEmissionReactions;
  std::vector< std::pair<int,Util::Expression *> > theSourceTerms;
  std::vector<int> masterSourceSpecies;
  double masterSourceValue;
  double sourceScaleFac;
  double C0;
  double t0;
  double x0;
  MaterialLayer *material;
  bool variableRateCoeffs;

  bool applySources;

  const SolverState & solState_;
};

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::clear
// Purpose       : make the reaction network an empty one
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::clear()
{
  speciesMap.clear();
  species.clear();
  constantsMap.clear();
  constants.clear();
  theReactions.clear();
  reactionNamesMap.clear();
  setName("NoName");
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setName
// Purpose       : Accessor function to set network name
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setName(const std::string &name)
{
  myName=name;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getNumSpecies
// Purpose       : Accessor function to returning number of species recorded
// Special Notes : Only solution species returned, does not include constants
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getNumSpecies()
{
  return (species.size());
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getNumConstants
// Purpose       : Accessor function to returning number of constant species
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/10/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getNumConstants()
{
  return (constants.size());
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getSpeciesName
// Purpose       : Accessor function to returning name of indicated species
// Special Notes : Only returns solution species names, not constants
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline const std::string & ReactionNetwork::getSpeciesName(int i)
{
  return species[i].getName();
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getConstantsName
// Purpose       : Accessor function to returning name of indicated constant
//                 species
// Special Notes : Only returns solution constants names, not variable species
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline const std::string & ReactionNetwork::getConstantsName(int i)
{
  return constants[i].getName();
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getSpeciesNum
// Purpose       : Accessor function to return number of named variable specie
// Special Notes :  returns -1 if the specified specie is not a variable
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getSpeciesNum(const std::string &name)
{
  std::map<std::string,int>::iterator n_i;
  n_i = speciesMap.find(name);

  if (n_i == speciesMap.end())
  {
    return -1;
  }
  else
  {
    return n_i->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::reactantExist
// Purpose       : Returns false if the named specie doesn't exist, or if
//                 it is a constant.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/24/07
//-----------------------------------------------------------------------------
inline bool ReactionNetwork::reactantExist(const std::string &name)
{
  bool retFlag(true);
  int i=getSpeciesNum(name);
  if (i == -1)
  {
    retFlag = false;  // this reactant does not exist.
  }

  return retFlag;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::constantExist
// Purpose       : Returns false if the named specie doesn't exist, or if
//                 it is a constant.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 4/24/07
//-----------------------------------------------------------------------------
inline bool ReactionNetwork::constantExist(const std::string &name)
{
  bool retFlag(true);
  int i=getConstantNum(name);
  if (i == -1)
  {
    retFlag = false;  // this constant does not exist.
  }

  return retFlag;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getReactantNum
// Purpose       : Accessor function to return number of named specie
// Special Notes :  Returns negative numbers for constants, positive for
//                  variables
//                 The species in question better exist, coz there's no way
//                 to return an error condition other than bombing.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getReactantNum(const std::string &name)
{
  int i=getSpeciesNum(name);
  if (i == -1)
  {
    i=getConstantNum(name);
    if (i == -1)
    {
      Report::UserFatal() << "Invalid species name specified: " << name;
    }
    i = -(i+1);
  }
  return i;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getConstantNum
// Purpose       : Accessor function to return number of named constant specie
// Special Notes :  returns -1 if the specified specie is not a constant
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getConstantNum(const std::string &name)
{
  std::map<std::string,int>::iterator n_i;
  n_i = constantsMap.find(name);

  if (n_i == constantsMap.end())
  {
    return -1;
  }
  else
  {
    return n_i->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addInitialCondition
// Purpose       : add an initial condition to the list
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::addInitialCondition(const std::string &speciesName, double value)
{
  initialConditions.push_back(std::pair<std::string,double>(speciesName,value));
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getNumInitialConditions
// Purpose       : get an initial condition from the list
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline int
ReactionNetwork::getNumInitialConditions()
{
  return(initialConditions.size());
}
//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getInitialCondition
// Purpose       : get an initial condition from the list
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline std::pair<std::string,double>
ReactionNetwork::getInitialCondition(int i)
{
  return(initialConditions[i]);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getReactionNum
// Purpose       : Accessor function to return number of named reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getReactionNum(std::string name)
{
  std::map<std::string,int>::iterator n_i;
  n_i = reactionNamesMap.find(name);

  if (n_i == reactionNamesMap.end())
  {
    return -1;
  }
  else
  {
    return n_i->second;
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getReaction
// Purpose       : Accessor function to returning reference to indexed reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------
inline Reaction & ReactionNetwork::getReaction(int i)
{
  return theReactions[i];
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getReaction
// Purpose       : Accessor function to returning reference to named reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/20/06
//-----------------------------------------------------------------------------

inline Reaction &ReactionNetwork::getReaction(const std::string &name)
{
  int ni;
  ni=getReactionNum(name);
  if (ni == -1)
  {
    Report::UserFatal() << "Attempt to access non-existant reaction " << name;
  }

  return theReactions[ni];
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setSimpleCalc
// Purpose       : set the named reaction's rate calculator to type Simple
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setSimpleCalc(const std::string &name, double k)
{
  getReaction(name).setSimpleRateCalculator(k,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setCaptureCalc
// Purpose       : set the named reaction's rate calculator to type Capture
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setCaptureCalc(const std::string &name, double sigma, double v)
{
  getReaction(name).setCaptureRateCalculator(sigma,v,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setElectronCaptureCalc
// Purpose       : set the named reaction's rate calculator to type Capture
// Special Notes : Specifically uses v_n=2.3e7 cm/s
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setElectronCaptureCalc(const std::string &name, double sigma)
{
  //carrier velocities are a function of the bulk device material
  //The reaction network is set up prior to the device and so material
  //properties are unknown when this method is called.  Instead of passing
  //through a velocity, I'm passing through a carrier charge so that
  //when the mulk material is known, this reaction can identify the carrier
  //that is being captured and the velocity can be set appropriately for
  //the device material. --LCM
  getReaction(name).setCaptureRateCalculator(sigma,-1.0,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setHoleCaptureCalc
// Purpose       : set the named reaction's rate calculator to type Capture
// Special Notes : uses v_p=1.9e7 cm/s
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setHoleCaptureCalc(const std::string &name, double sigma)
{
  //carrier velocities are a function of the bulk device material
  //The reaction network is set up prior to the device and so material
  //properties are unknown when this method is called.  Instead of passing
  //through a velocity, I'm passing through a carrier charge so that
  //when the mulk material is known, this reaction can identify the carrier
  //that is being captured and the velocity can be set appropriately for
  //the device material. --LCM
  getReaction(name).setCaptureRateCalculator(sigma,+1.0,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setEmissionCalc(const std::string &name, double sigma, double v, double N, double E)
{
  getReaction(name).setEmissionRateCalculator(sigma,v,N,E,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setElectronEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :Uses v_n=2.3e7 cm/s and N_c=2.86e19 cm^{-3}
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setElectronEmissionCalc(const std::string &name, double sigma, double E)
{
  //carrier velocities and DOS are a function of the bulk device material
  //The reaction network is set up prior to the device and so material
  //properties are unknown when this method is called.  Instead of passing
  //through a velocity, I'm passing through a carrier charge so that
  //when the mulk material is known, this reaction can identify the carrier
  //that is being emitted and the velocity can be set appropriately for
  //the device material. --LCM
  getReaction(name).setEmissionRateCalculator(sigma,-1.0,2.86e19,E,C0,t0,x0);
}
//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setHoleEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :  uses v_p=1.9e7 cm/s and N_v=2.66e19 cm^{-3}
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setHoleEmissionCalc(const std::string &name, double sigma, double E)
{
  //carrier velocities and DOS are a function of the bulk device material
  //The reaction network is set up prior to the device and so material
  //properties are unknown when this method is called.  Instead of passing
  //through a velocity, I'm passing through a carrier charge so that
  //when the mulk material is known, this reaction can identify the carrier
  //that is being emitted and the velocity can be set appropriately for
  //the device material. --LCM
  getReaction(name).setEmissionRateCalculator(sigma,+1.0,2.66e19,E,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setComplexCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setComplexCalc(const std::string &name)
{
  getReaction(name).setComplexRateCalculator(species,constants,C0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setComplexMultiplierCalc
// Purpose       : set the named reaction's rate calculator to type complex with
//               : a multiplier on the rate
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 06/26/2014
//-----------------------------------------------------------------------------
 inline void ReactionNetwork::setComplexMultiplierCalc(const std::string &name, double multiplier)
{
  getReaction(name).setComplexMultiplierRateCalculator(species,constants,C0,t0,x0,multiplier);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setDecomplexCalc
// Purpose       : set the named reaction's rate calculator to type emission
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setDecomplexCalc(const std::string &name,
                                              double bindingEnergy,
                                              double gammaAB,
                                              double gammaA,
                                              double gammaB,
                                              double concSi)
{

  getReaction(name).setDecomplexRateCalculator(species,constants,
                                               bindingEnergy,gammaAB,gammaA,
                                               gammaB,concSi,C0,t0,x0);
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setBourgoinCorbettHoleCalc
// Purpose       : set the named reaction's rate calculator to type Bourgoin
//                 Corbett hole-enhanced diffusion
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/25/14
//-----------------------------------------------------------------------------
 inline void ReactionNetwork::setBourgoinCorbettHoleCalc(const std::string &name,double sigma)
{

  getReaction(name).setBourgoinCorbettHoleRateCalculator(species,constants,
                                                     sigma,
                                                     C0,t0,x0);
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 05/27/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::output(std::ostream & os) const
{
  int i;

  for (i=0;i< species.size();++i)
  {
    os << "species["<<i<<"] = " << species[i].getName() << std::endl;
  }
  os << std::endl;


  for (i=0;i<theReactions.size();++i)
  {
    os << reactionNames[i];
    theReactions[i].output( species, os);
  }

  if (electronCaptureReactions.size() != 0)
  {
    os << "Electron Capture Reactions: " << std::endl;
    for ( i=0; i<electronCaptureReactions.size(); ++i)
    {
      os << "  Reaction number " << electronCaptureReactions[i] << "("
         << reactionNames[electronCaptureReactions[i]] << ")" << std::endl;
    }
  }
  if (holeCaptureReactions.size() != 0)
  {
    os << "Hole Capture Reactions: " << std::endl;
    for ( i=0; i<holeCaptureReactions.size(); ++i)
    {
      os << "  Reaction number " << holeCaptureReactions[i] << "("
         << reactionNames[holeCaptureReactions[i]] << ")" << std::endl;
    }
  }
  if (electronEmissionReactions.size() != 0)
  {
    os << "Electron Emission Reactions: " << std::endl;
    for ( i=0; i<electronEmissionReactions.size(); ++i)
    {
      os << "  Reaction number " << electronEmissionReactions[i] << "("
         << reactionNames[electronEmissionReactions[i]] << ")" << std::endl;
    }
  }
  if (holeEmissionReactions.size() != 0)
  {
    os << "Hole Emission Reactions: " << std::endl;
    for ( i=0; i<holeEmissionReactions.size(); ++i)
    {
      os << "  Reaction number " << holeEmissionReactions[i] << "("
         << reactionNames[holeEmissionReactions[i]] << ")" << std::endl;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 05/27/06
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const ReactionNetwork & rn)
{
  os << "Reaction Network: " << std::endl;
  rn.output(os);

  return os;
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getERate
// Purpose       : Compute the total rate at which electrons are "consumed" or
//                 "produced" by all the capture and emission reactions,
//                 if there are any.  This can be used even if the electron
//                 concentration is held fixed (it'll just be the sum of all
//                 reaction rates involving electrons).
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+E or A+E=>B and will be incorrect if
//                 any reaction not of this form is given a name with
//                 "_ELECTRON_EMISSION" or "_ELECTRON_CAPTURE" in it.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getERate(std::vector<double> &concs,
                                        std::vector<double> &constant_vec)
{
  return getRate(concs,constant_vec,electronCaptureReactions,
                 electronEmissionReactions);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getELifetime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getELifetime(std::vector<double> &concs,
                                            std::vector<double> &constant_vec)
{

  double eConc;
  int concNum=getReactantNum("E");
  // Bleah
  if (concNum < 0)
    eConc = constant_vec[-(concNum+1)];
  else
    eConc = concs[concNum];


  return getCaptureLifetime(concs,constant_vec,electronCaptureReactions,
                            eConc);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getELifetimes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::getELifetimes(std::vector<double> &concs,
                                           std::vector<double> &constant_vec,
                                           std::vector<double> &lifetimes)
{

  double eConc;
  int concNum=getReactantNum("E");
  // Bleah
  if (concNum < 0)
    eConc = constant_vec[-(concNum+1)];
  else
    eConc = concs[concNum];


  getCaptureLifetimes(concs,constant_vec,electronCaptureReactions,
                      eConc,lifetimes);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDERateDC
// Purpose       : compute vector of derivatives of net electron emission
//                 rate (e.g. result of getERate) with respect to concentration
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::getDERateDC(std::vector<double> &concs,
                                         std::vector<double> &constant_vec,
                                         std::vector<double> &dRatedC)
{
  getDRateDC(concs,constant_vec,electronCaptureReactions,
             electronEmissionReactions,dRatedC);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDERateDConst
// Purpose       : compute vector of derivatives of net electron emission
//                 rate (e.g. result of getERate) with respect to concentration
//
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 11/15/08
//-----------------------------------------------------------------------------
inline void ReactionNetwork::getDERateDConst(std::vector<double> &concs,
                                             std::vector<double> &constant_vec,
                                             std::vector<double> &dRatedConst)
{
  getDRateDConst(concs,constant_vec,electronCaptureReactions,
                 electronEmissionReactions,dRatedConst);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getHRate
// Purpose       : Compute the total rate at which holes are "consumed" or
//                 "produced" by all the capture and emission reactions,
//                 if there are any.  This can be used even if the electron
//                 concentration is held fixed (it'll just be the sum of all
//                 reaction rates involving electrons).
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+H or A+H=>B and will be incorrect if
//                 any reaction not of this form is given a name with
//                 "_HOLE_EMISSION" or "_HOLE_CAPTURE" in it.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getHRate(std::vector<double> &concs,
                                        std::vector<double> &constant_vec)
{
  return getRate(concs,constant_vec,holeCaptureReactions,
                 holeEmissionReactions);;
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getHLifetime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getHLifetime(std::vector<double> &concs,
                                            std::vector<double> &constant_vec)
{

  double hConc;
  int concNum=getReactantNum("H");
  // Bleah
  if (concNum < 0)
    hConc = constant_vec[-(concNum+1)];
  else
    hConc = concs[concNum];

  return getCaptureLifetime(concs,constant_vec,holeCaptureReactions,
                            hConc);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getHLifetimes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::getHLifetimes(std::vector<double> &concs,
                                           std::vector<double> &constant_vec,
                                           std::vector<double> &lifetimes)
{

  double hConc;
  int concNum=getReactantNum("H");
  // Bleah
  if (concNum < 0)
    hConc = constant_vec[-(concNum+1)];
  else
    hConc = concs[concNum];


  getCaptureLifetimes(concs,constant_vec,holeCaptureReactions,
                      hConc,lifetimes);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDHRateDC
// Purpose       : compute vector of derivatives of net hole emission
//                 rate (e.g. result of getHRate) with respect to concentration
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::getDHRateDC(std::vector<double> &concs,
                                         std::vector<double> &constant_vec,
                                         std::vector<double> &dRatedC)
{
  getDRateDC(concs,constant_vec,holeCaptureReactions,
             holeEmissionReactions,dRatedC);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDHRateDConst
// Purpose       : compute vector of derivatives of net hole emission
//                 rate (e.g. result of getHRate) with respect to concentration
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline void ReactionNetwork::getDHRateDConst(std::vector<double> &concs,
                                             std::vector<double> &constant_vec,
                                             std::vector<double> &dRatedConst)
{
  getDRateDConst(concs,constant_vec,holeCaptureReactions,
                 holeEmissionReactions,dRatedConst);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDiffusionCoefficient
// Purpose       :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getDiffusionCoefficient
(const std::string & name, const double temp)
{
  int num = getSpeciesNum(name);
  double D = 0.0;

  if (num < 0)
    D = 0.0;
  else
    D = species[num].getDiffusionCoefficient(temp);

  return D;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDiffusionCoefficient
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getDiffusionCoefficient
(int specie, const double temp)
{
  double D = 0.0;

  if (specie < 0)
    D = 0.0;
  else
    D = species[specie].getDiffusionCoefficient(temp);

  return D;
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDiffusionCoefficient
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
inline double ReactionNetwork::getDiffusionCoefficient
  (int specie, const double temp, 
                                  std::vector<double> concs, 
                                  std::vector<double> carrierConcs)
{
  double D = 0.0;

  if (specie < 0)
    D = 0.0;
  else
    D = species[specie].getDiffusionCoefficient(temp, concs, carrierConcs);

  return D;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getChargeState
// Purpose       :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 02/22/2012
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getChargeState
(const std::string & name)
{
  int num = getSpeciesNum(name);
  int z = 0;

  if (num < 0)
    z = 0;
  else
    z = species[num].getChargeState();

  return z;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getChargeState
// Purpose       :
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 10/31/06
//-----------------------------------------------------------------------------
inline int ReactionNetwork::getChargeState
(int specie)
{
  int z = 0;

  if (specie < 0)
    z = 0;
  else
    z = species[specie].getChargeState();

  return z;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setApplySources
// Purpose       :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 04/19/09
//-----------------------------------------------------------------------------
inline void ReactionNetwork::setApplySources(bool flag)
{
  applySources = flag;
}

} // namespace Device
} // namespace Xyce

#endif

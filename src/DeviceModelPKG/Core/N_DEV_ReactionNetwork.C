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
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 03/20/2006
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


//Standard includes

#include <iostream>
#include <sstream>

// Xyce Includes
#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_CheckIfValidFile.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_DEV_ReactionNetwork.h>
#include <N_DEV_MaterialLayer.h>

#include <N_DEV_SolverState.h>

#include <N_UTL_BreakPoint.h>

#include <N_DEV_ExpressionGroupWrapper.h>

#ifdef Xyce_REACTION_PARSER
// Grrrr.  Stupid bison 2.4 stopped putting the pre-prologue into the header.
// need this forward declaration
namespace Xyce {
namespace Device {
class ReactionLexer;
}}

#include "N_DEV_ReactionParser.hxx"
// BLEAH!   This is here DUPLICATED from N_DEV_ReactionParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
#undef yyFlexLexer
  /* CAREFUL watch continuations! */
#define YY_DECL \
int ReactionLexer::getToken(ReactionParser::semantic_type *lvalp,  \
                            location *llocp, \
                            std::map<std::string, int> &theSpeciesIDs)

  // YECH!  Work around very stupid way that multiple parsers/lexers are 
  // handled.
  // Bison's "%name-prefix" is implemented as a #define yylex "prefix"lex
  // which BREAKS flex's C++ lexer: it contains a method "yylex" in the
  // yyFlexLexer class.  Unless we do this kludge, that method gets renamed
  // with the define as well, and the result is a broken set of classes
#undef yylex
#include <FlexLexer.h>
#include <N_DEV_ReactionLexer.h>
  // if we actually *used* yylex anywhere here it would be necessary to 
  // undo that kludge.  Note that because of this stupidity, if the 
  // "%name-prefix" is changed, this line needs to be changed, too.
  // BUT we don't actually use yylex anywhere in this file, so let's
  // leave yylex undefined.  If later it turns out that this *becomes*
  // necessary, uncomment the next line.
  //  #define yylex XyceDevicelex

#endif

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::ReactionNetwork
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
ReactionNetwork::ReactionNetwork(
    const SolverState & solver_state,
    const std::string &name)
  : myName(name),
    sourceScaleFac(1.0),
    C0(1.0),
    t0(1.0),
    x0(1.0),
    applySources(true), 
    solState_(solver_state)
{
  theReactions.reserve(10); // try to cut down on copies
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::ReactionNetwork
// Purpose       : Copy constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
ReactionNetwork::ReactionNetwork(
   const ReactionNetwork &right)
  : speciesMap(right.speciesMap),
    species(right.species),
    constantsMap(right.constantsMap),
    constants(right.constants),
    initialConditions(right.initialConditions),
    theReactions(right.theReactions),
    reactionNamesMap(right.reactionNamesMap),
    reactionNames(right.reactionNames),
    myName(right.myName),
    electronCaptureReactions(right.electronCaptureReactions),
    holeCaptureReactions(right.holeCaptureReactions),
    electronEmissionReactions(right.electronEmissionReactions),
    holeEmissionReactions(right.holeEmissionReactions),
    sourceScaleFac(right.sourceScaleFac),
    C0(right.C0),
    t0(right.t0),
    x0(right.x0),
    material(right.material),
    applySources(right.applySources),
    solState_(right.solState_)
{

  // Can't just copy the vector of source terms, coz those are pointers.
  // Need to do deep copy:
  int i,sourceSize;
  sourceSize=right.theSourceTerms.size();
  theSourceTerms.reserve(sourceSize);
  for (i=0;i<sourceSize;++i)
  {
    Util::Expression *newSource;
    // make a new one as a copy of the right one
    newSource = new Util::Expression(*((right.theSourceTerms[i]).second));
    theSourceTerms.push_back(std::pair<int,Util::Expression *>(
                                (right.theSourceTerms[i]).first,
                                newSource) );
  }

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::~ReactionNetwork
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
ReactionNetwork::~ReactionNetwork()
{

  // It is no longer the case that we store only copies of expression
  // pointers allocated from other classes.  We now *must* destroy expressions
  // we have created in parsing.  Since it is difficult to keep track of
  // which expression pointers were passed in to us and which we created
  // ourselves, we no longer store those passed-in pointers, but rather
  // pointers to copies of the pointed-to expression.

  // Walk backwards through the sourceTerms vector, deleting the expression
  // and removing the pair from the vector.
  while (!theSourceTerms.empty())
  {
    std::pair<int,Util::Expression *> tempPair = theSourceTerms.back();
    theSourceTerms.pop_back();
    delete tempPair.second;
  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setReactionNetworkFromFile
// Purpose       : Parse a reaction network specification file and set the
//                 network up.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/24/06
//-----------------------------------------------------------------------------
void
ReactionNetwork::setReactionNetworkFromFile(const NetlistLocation &netlist_location, const std::string &fileName)
{
#ifdef Xyce_REACTION_PARSER
  if (fileName != "") // do nothing if no file specified
  {
    // Get the reaction set from a file
    std::map<std::string,int> theSpeciesIDs;

    // Error out if the user-specified file does not exist, cannot
    // be opened, or is a directory name rather than a file name.  
    // See SON Bug 785 and SRN Bug 2100 for more details.
    if ( !(Util::checkIfValidFile(fileName)) )
    {
      Report::UserFatal() << "Cannot find reaction specification file " << fileName;
    }

    // temporary hard-coded filename
    std::ifstream reactionFile(fileName.c_str(),std::ios::in);
    if (reactionFile.is_open())
    {
      ReactionLexer theLexer(netlist_location, fileName, &reactionFile);
      XyceDevice::ReactionParser theParser(&theLexer, theSpeciesIDs, *this);

      if (theParser.parse() != 0)
      {
        if (DEBUG_DEVICE)
        {
          Xyce::dout() << *this << std::endl;
        }
      }
    }
    else
    {
      Report::UserFatal() << "Cannot open reaction specification file " << fileName;
    }
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setSpecies
// Purpose       : register a list of species names used in reactions.
//
// Special Notes : Once registered this way, it is assumed that all vectors
//                 of species concentrations are ordered in the same way as
//                 the provided vector of species names.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setSpecies(std::vector<Specie> &theSpeciesVect)
{

  int i;
  int nspec=theSpeciesVect.size();

  speciesMap.clear();
  species=theSpeciesVect;

  for (i=0;i<nspec;++i)
  {
    speciesMap[theSpeciesVect[i].getName()]=i;
  }

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setConstants
// Purpose       : register a list of constant species names used in reactions.
//
// Special Notes : Once registered this way, it is assumed that all
//                 vectors of constant species concentrations are ordered in
//                 the same way as the provided vector of species names.
//
//                 Constant species are those whose concentrations remain
//                 unaltered by the reaction network, e.g. those that are
//                 fed into the system at fixed rate as in a CSTR.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setConstants(std::vector<Specie> &theConstantsVect)
{

  int i;
  int nspec=theConstantsVect.size();

  constantsMap.clear();
  constants = theConstantsVect;

  for (i=0;i<nspec;++i)
  {
    constantsMap[theConstantsVect[i].getName()]=i;
  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addSpecie
// Purpose       : add one to the list of species names used in reactions.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addSpecie(const Specie &aSpecie)
{

  int i=species.size(); // get index of next one we add
  species.push_back(aSpecie);
  speciesMap[aSpecie.getName()]=i;

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addConstant
// Purpose       : add one to the list of constant names used in reactions.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addConstant(const Specie &aConstant)
{

  int i=constants.size(); // get index of next one we add
  constants.push_back(aConstant);
  constantsMap[aConstant.getName()]=i;

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addReaction
// Purpose       : Create a new named reaction in the network
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addReaction(const std::string &name)
{
  Xyce::Device::Reaction dummy;
  //Reaction dummy;

  addReaction(name, dummy);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addReaction
// Purpose       : Add a reaction object into the network and give it a name
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
  void ReactionNetwork::addReaction(const std::string &name, Xyce::Device::Reaction &reaction)
{
  // first check if there is already a reaction of that name in the
  // reactions map:

  if (reactionNamesMap.find(name) != reactionNamesMap.end())
  {
    Report::DevelFatal() << "  Attempt to add reaction duplicate name " << name;
  }

  //try to cut down on calls to copy constructor
  if (theReactions.capacity()<=theReactions.size())
  {
    theReactions.reserve(theReactions.size()+10);
  }
  theReactions.push_back(reaction);
  reactionNames.push_back(name);
  reactionNamesMap[name]=theReactions.size()-1;

  // Now check if the reaction has "_ELECTRON_CAPTURE",
  // "_ELECTRON_EMISSION", "_HOLE_CAPTURE" or "_HOLE_EMISSION" in it.
  // If so, keep track of it in a separate vector for such things.
  // we do this so we can later get at them quickly without doing name
  // searches.
  if (name.find("_ELECTRON_CAPTURE",0) != std::string::npos)
  {
    electronCaptureReactions.push_back(theReactions.size()-1);
  }
  else if (name.find("_HOLE_CAPTURE",0) != std::string::npos)
  {
    holeCaptureReactions.push_back(theReactions.size()-1);
  }
  else if (name.find("_ELECTRON_EMISSION",0) != std::string::npos)
  {
    electronEmissionReactions.push_back(theReactions.size()-1);
  }
  else if (name.find("_HOLE_EMISSION",0) != std::string::npos)
  {
    holeEmissionReactions.push_back(theReactions.size()-1);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addReactant
// Purpose       : add a species/stoichiometric coefficient pair to the
//                 reactant list of a specified named reaction
//                 The species is specified using its name, and can be either
//                 one of those in the list provided by setSpecies or
//                 setConstants.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addReactant(const std::string &name, const std::string &reactant, double stoich)
{

  std::map<std::string,int>::iterator n_i;
  // check that name exists in map
  int reactionNum;
  reactionNum=getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to add reactant " << reactant
                         << " to non-existant reaction " << name;
  }
  else
  {
    int speciesNum;
    // now check if we know this reactant

    n_i = speciesMap.find(reactant);
    if (n_i == speciesMap.end())
    {
      // not a solution species --- is it a constant?
      n_i = constantsMap.find(reactant);
      if (n_i == constantsMap.end())
      { // nope
        Report::DevelFatal() << "attempt to add unknown reactant " << reactant
                             << " to reaction number " << reactionNum
                             << "("<<name<<")";
      }
      else
      {
        speciesNum = -(n_i->second+1);
      }
    }
    else
    {
      speciesNum = n_i->second;
    }


    getReaction(reactionNum).addReactant(speciesNum,stoich);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addProduct
// Purpose       : add a species/stoichiometric coefficient pair to the
//                 product list of a specified named reaction
//                 The species is specified using its name, and must be
//                 one of those in the list provided by setSpecies
//                 (must not be one from the Constants list).
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addProduct(const std::string &name, const std::string &product, double stoich)
{

  std::map<std::string,int>::iterator n_i;

  // check that reaction name exists in map
  // check that name exists in map
  int reactionNum;
  reactionNum=getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to add product " << product
                         << " to non-existant reaction " << name;
  }
  else
  {
    int speciesNum;
    // now check if we know this product

    n_i = speciesMap.find(product);
    if (n_i == speciesMap.end())
    {
      // not a solution species --- is it a constant?
      n_i = constantsMap.find(product);
      if (n_i == constantsMap.end())
      { // nope
        Report::DevelFatal()  << "attempt to add unknown product " << product
                              << " to reaction number " << reactionNum
                              << "("<<name<<")";
      }
      else
      {
        std::ostringstream ost;

#ifdef Xyce_RXN_WARNINGS
        Report::UserWarning() << " Specified constant species " << product
                              << " as product of reaction number " << reactionNum
                              << "(" << name << ")"<<std::endl
                              << " IGNORING that product";
#endif
      }
    }
    else
    {
      speciesNum = n_i->second;
      getReaction(reactionNum).addProduct(speciesNum,stoich);
    }


  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setRateConstant
// Purpose       : set the rate constant for the named reaction
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setRateConstant(const std::string &name, double k)
{
  int reactionNum;
  // check that reaction name exists in map
  reactionNum=getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to set rate constant of non-existant reaction "
                         << name;
  }
  else
  {
    getReaction(reactionNum).setRateConstant(k);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::scaleRateConstant
// Purpose       : set the rate constant for the named reaction
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 06/01/06
//-----------------------------------------------------------------------------
void ReactionNetwork::scaleRateConstant(const std::string &name, double kscale)
{

  int reactionNum;
  // check that reaction name exists in map
  reactionNum = getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to scale rate constant of non-existant reaction "
                         << name;
  }
  else
  {
    getReaction(reactionNum).scaleRateConstant(kscale);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setRateConstantFromCalculator
// Purpose       : set the rate constant for the named reaction according
//                 to its saved method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setRateConstantFromCalculator (const std::string &name, double T)
{

  int reactionNum;
  // check that reaction name exists in map
  reactionNum = getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to scale rate constant of non-existant reaction "
                         << name;
  }
  else
  {
    getReaction(reactionNum).setRateConstantFromCalculator(T);
  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setRateConstantFromCalculator
// Purpose       : set the rate constant for the named reaction according
//                 to its saved method
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL, Electrical and Microsystems Modeling
// Creation Date : 04/17/2014
//-----------------------------------------------------------------------------
void ReactionNetwork::setRateConstantFromCalculator (const std::string &name, double T,
                                     std::vector<double> &concs,
                                     std::vector<double> &constant_vec)
{

  int reactionNum;
  // check that reaction name exists in map
  reactionNum = getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to scale rate constant of non-existant reaction "
                         << name;
  }
  else
  {
    getReaction(reactionNum).setRateConstantFromCalculator(T,concs,constant_vec);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::scaleRateConstantFromCalculator
// Purpose       : scale the rate constant for the named reaction according
//                 to its saved method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
void ReactionNetwork::scaleRateConstantFromCalculator(const std::string &name)
{

  int reactionNum;
  // check that reaction name exists in map
  reactionNum = getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to scale rate constant of non-existant reaction "
                         << name;
  }
  else
  {
    getReaction(reactionNum).scaleRateConstantFromCalculator();
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::unscaleRateConstantFromCalculator
// Purpose       : unscale the rate constant for the named reaction according
//                 to its saved method
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/01/06
//-----------------------------------------------------------------------------
void ReactionNetwork::unscaleRateConstantFromCalculator (const std::string &name)
{

  int reactionNum;
  // check that reaction name exists in map
  reactionNum = getReactionNum(name);
  if (reactionNum == -1)
  {
    Report::DevelFatal() << " Attempt to scale rate constant of non-existant reaction "
                         << name;
  }
  else
  {
    getReaction(reactionNum).unscaleRateConstantFromCalculator();
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setRateConstantsFromCalc
// Purpose       : set the rate constant for the all reaction according
//                 to their saved methods
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/02/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setRateConstantsFromCalc(double T)
{
  int numReacts=theReactions.size();
  int i;

  for (i=0; i<numReacts; ++i)
  {
    theReactions[i].setRateConstantFromCalculator(T);
    theReactions[i].setTemperature(T);
  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setRateConstantsFromCalc
// Purpose       : set the rate constant for the all reaction according
//                 to their saved methods
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/02/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setRateConstantsFromCalc(double T,
                                std::vector<double> &concs,
                                std::vector<double> &constant_vec)
{
  int numReacts=theReactions.size();
  int i;
  for (i=0; i<numReacts; ++i)
  {
    theReactions[i].setRateConstantFromCalculator(T,concs,constant_vec);
    theReactions[i].setTemperature(T);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::scaleRateConstantsFromCalc
// Purpose       : scale the rate constant for the all reaction according
//                 to their saved methods
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/02/06
//-----------------------------------------------------------------------------
void ReactionNetwork::scaleRateConstantsFromCalc()
{
  int numReacts=theReactions.size();
  int i;
  for (i=0; i<numReacts; ++i)
  {
    theReactions[i].scaleRateConstantFromCalculator();
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::unscaleRateConstantsFromCalc
// Purpose       : unscale the rate constant for the all reaction according
//                 to their saved methods
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/02/06
//-----------------------------------------------------------------------------
void ReactionNetwork::unscaleRateConstantsFromCalc()
{
  int numReacts=theReactions.size();
  int i;
  for (i=0; i<numReacts; ++i)
  {
    theReactions[i].unscaleRateConstantFromCalculator();
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setScaleParams
// Purpose       : set the scale parameters and pass them down to any existing
//                 reactions.  Done this way, we make sure that all reactions
//                 in the network are reset to the same set of scale
//                 parameters.
// Special Notes : If called before reactions defined, it's OK, because as
//                 each reaction's rate calculator is created, the stored
//                 C0, t0, and x0 are passed to it.  This function does the
//                 calls down to existing reactions just in case they were
//                 previously defined.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/24/06
//-----------------------------------------------------------------------------
void ReactionNetwork::setScaleParams(double c,double t, double x)
{
  std::vector<Xyce::Device::Reaction>::iterator reactionIter=theReactions.begin();
  std::vector<Xyce::Device::Reaction>::iterator reactionIterEnd=theReactions.end();
  C0=c;
  t0=t;
  x0=x;

  for (; reactionIter != reactionIterEnd; ++reactionIter)
  {
    reactionIter->setScaleFactors(C0,t0,x0);
  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setmaterial
// Purpose       : Set the bulk semiconductor material
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 08/19/14
//-----------------------------------------------------------------------------
  void ReactionNetwork::setMaterial(MaterialLayer *mat, double Temp)
{
  material = mat;

  std::vector<Reaction>::iterator reactionIter=theReactions.begin();
  std::vector<Reaction>::iterator reactionIterEnd=theReactions.end();
  for (; reactionIter != reactionIterEnd; ++reactionIter)
    reactionIter->setMaterial(material, Temp);

  //This is a good opportunity to check the species for bourgoin-corbett
  //enhanced diffusion and then to set some parameters
  for(int i=0 ; i<species.size() ; ++i)
    {
      if(species[i].getBCCarrierCharge() == 1)
        {
          species[i].setBCThermalVelocity(material->holeThermalV);
        }
      else if(species[i].getBCCarrierCharge() == -1)
        {
          species[i].setBCThermalVelocity(material->electronThermalV);
        }
    }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setmaterial
// Purpose       : Set the bulk semiconductor material
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 08/19/14
//-----------------------------------------------------------------------------
  void ReactionNetwork::setCoefficients(double Temp)
{
  std::vector<Reaction>::iterator reactionIter=theReactions.begin();
  std::vector<Reaction>::iterator reactionIterEnd=theReactions.end();
  for (; reactionIter != reactionIterEnd; ++reactionIter)
    reactionIter->setCoefficient(Temp);
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setmaterial
// Purpose       : Set the bulk semiconductor material
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 08/19/14
//-----------------------------------------------------------------------------
  void ReactionNetwork::setRxnVariableCoeffs(bool variableCoeffs)
{
  variableRateCoeffs=variableCoeffs;

  std::vector<Reaction>::iterator reactionIter=theReactions.begin();
  std::vector<Reaction>::iterator reactionIterEnd=theReactions.end();
  for (; reactionIter != reactionIterEnd; ++reactionIter)
    reactionIter->setRxnVariableCoeffs(variableCoeffs);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addSourceTerm
// Purpose       : add an expression for a source term
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/04/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addSourceTerm(const std::string &speciesName, const std::string &expressionStr)
{
  if (applySources)
  {
    int speciesNum=getSpeciesNum(speciesName);
    if (speciesNum >= 0) // the species exists
    {
      Util::Expression *foo= new Util::Expression(solState_.getGroupWrapper()->expressionGroup_,expressionStr);
      theSourceTerms.push_back( std::pair<int,Util::Expression *>(speciesNum, foo));
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addSourceTerm
// Purpose       : add an expression for a source term given
//                 already-constructed expression object.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/04/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addSourceTerm(const std::string &speciesName, Util::Expression *expression)
{
  int speciesNum=getSpeciesNum(speciesName);
  Util::Expression * ExpressionCopy=new Util::Expression(*expression);

  if (speciesNum >= 0) // the species exists
  {
    theSourceTerms.push_back(
       std::pair<int,Util::Expression *>(speciesNum,ExpressionCopy));
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::addMasterSourceTerm
// Purpose       : add an expression for a source term that depends on a
//                 "master" source value.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/04/06
//-----------------------------------------------------------------------------
void ReactionNetwork::addMasterSourceTerm(const std::string &speciesName)
{
  int speciesNum=getSpeciesNum(speciesName);

  if (speciesNum >= 0)   // the species exists
  {
    masterSourceSpecies.push_back(speciesNum);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool ReactionNetwork::getBreakPoints(std::vector<Util::BreakPoint> & breakPointTimes)
{
  std::vector< std::pair<int,Util::Expression *> >::iterator iterSource=
    theSourceTerms.begin();
  std::vector< std::pair<int,Util::Expression *> >::iterator source_end=
    theSourceTerms.end();

  for (;iterSource != source_end; ++iterSource)
  {
    (iterSource->second)->getBreakPoints(breakPointTimes);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDdt
// Purpose       : compute the time derivative of all species concentrations
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::getDdt(std::vector<double> &concs,std::vector<double> &constant_vec,
              std::vector<double> &ddt)
{
  int rSize=theReactions.size();
  int i;
  for (i=0;i<rSize;++i)
  {
    theReactions[i].getDdt(concs,constant_vec,ddt);
  }

  // add in the source terms
  std::vector< std::pair<int,Util::Expression *> >::iterator iterSource=
    theSourceTerms.begin();
  std::vector< std::pair<int,Util::Expression *> >::iterator source_end=
    theSourceTerms.end();

  for (;iterSource != source_end; ++iterSource)
  {
    double return_val;
    (iterSource->second)->evaluateFunction(return_val); // ignore return val?  That's
                                               // what's done everywhere else
    ddt[(iterSource->first)] +=    sourceScaleFac*return_val;
  }

  // Add in the master source terms:
  std::vector<int>::iterator iterMasterSource=masterSourceSpecies.begin();
  std::vector<int>::iterator masterSource_end=masterSourceSpecies.end();
  for (; iterMasterSource!= masterSource_end; ++iterMasterSource)
  {
    ddt[*iterMasterSource] += masterSourceValue*sourceScaleFac;
  }

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getJac
// Purpose       : compute the jacobian of the reaction network, the
//                 derivatives of the time derivatives with respect to species
//                 concentration.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::getJac(std::vector<double> &concs, std::vector<double> &constant_vec,
              std::vector<std::vector<double> >&jac)
{
  int rSize=theReactions.size();
  int i;

  for (i=0;i<rSize;++i)
  {
    theReactions[i].getJac(concs,constant_vec,jac);
  }
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDFdConst
// Purpose       : compute the derivatives of the RHS with respect to
//                 a specified constant concentration.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/17/06
//-----------------------------------------------------------------------------
void ReactionNetwork::getDFdConst( const std::string & constantName,
                                              std::vector<double> &concs,
                                              std::vector<double> &constant_vec,
                                              std::vector<double> &dFdConst)
{
  int rSize=theReactions.size();
  int cSize=concs.size();
  int i;
  int constNum=getConstantNum(constantName);

 dFdConst.resize(cSize);
  for (i=0;i<cSize;++i)
    dFdConst[i]=0.0;

  for (i=0;i<rSize;++i)
  {
    theReactions[i].getDFdConst(constNum,concs,constant_vec,dFdConst);
  }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getJac
// Purpose       : compute the jacobian of the reaction network, the
//                 derivatives of the time derivatives with respect to species
//                 concentration.
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/2014
//-----------------------------------------------------------------------------
void ReactionNetwork::getJacobianVC(std::vector<double> &concs, std::vector<double> &constant_vec,
                                    std::vector<std::vector<double> >&jac, std::vector<double> &constVec)
{
  int rSize=theReactions.size();
  int cSize=concs.size();
  int i;

  constVec.resize(2*cSize);
  for (i=0 ; i<constVec.size() ; ++i)
    constVec[i]=0.0;

  for (i=0;i<rSize;++i)
  {
    theReactions[i].getJacobianVC(concs,constant_vec,jac,constVec);
  }

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getRate
// Purpose       : Compute the total rate at which species are "consumed" or
//                 "produced" by all the capture and emission reactions,
//                 if there are any.  This can be used even if the species
//                 concentration is held fixed (it'll just be the sum of all
//                 reaction rates involving electrons).
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+E or A+E=>B and will be incorrect if
//                 the number of any reaction not of this form is listed in
//                 the capture or emission vectors.
//
//                 Generic version in which the list of capture and emission
//                 reactions are passed in
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/13/06
//-----------------------------------------------------------------------------
double ReactionNetwork::getRate(std::vector<double> &concs,
                                           std::vector<double> &constant_vec,
                                           std::vector<int> &captureVect,
                                           std::vector<int> &emissionVect)
{
  int i;
  double rate=0;

  for (i=0;i<captureVect.size();++i)
  {
    rate -= theReactions[captureVect[i]].getRateVC(concs, constant_vec);
  }

  for (i=0;i<emissionVect.size();++i)
  {
    if(theReactions[emissionVect[i]].getCarrierEmissionIndex() < 0)
      {
        rate += theReactions[emissionVect[i]].getRateVC(concs,constant_vec);
      }
    else
      {
        rate += theReactions[emissionVect[i]].getFDEmissionRate(concs,constant_vec);
      }
  }


  return rate;
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getCaptureLifetime
// Purpose       : Given a list of capture reaction numbers and the
//                 species number of the concentration consumed by those
//                 reactions, return the lifetime.
// Special Notes : Returns -1 if lifetime would be infinite (zero capture rate)
//                 This is a kludge, but one that is easily tested by anyone
//                 who uses this.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
double ReactionNetwork::getCaptureLifetime(std::vector<double> &concs,
                                           std::vector<double> &constant_vec,
                                           std::vector<int> &captureVect,
                                           double & concentration)
{
  int i;
  double rate=0;

  for (i=0;i<captureVect.size();++i)
  {
    rate += theReactions[captureVect[i]].getRateVC(concs, constant_vec);
  }

  // inverse of lifetime = rate/concentration, so...
  if (rate > 0) // don't try to compute infinite lifetime!
    return concentration/rate;
  else
    return (-1.0); // bollocks, but let's not quibble.   This is just as
                   // invalid a lifetime as infinity would be!

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getCaptureLifetimes
// Purpose       : Given a list of capture reaction numbers and the
//                 species number of the concentration consumed by those
//                 reactions, return a vector of lifetimes due to each
//                 reaction.
// Special Notes : Returns -1 if lifetime would be infinite (zero capture rate)
//                 This is a kludge, but one that is easily tested by anyone
//                 who uses this.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/20/06
//-----------------------------------------------------------------------------
void ReactionNetwork::getCaptureLifetimes(std::vector<double> &concs,
                                           std::vector<double> &constant_vec,
                                           std::vector<int> &captureVect,
                                           double & concentration,
                                           std::vector<double> &lifetimes)
{
  int i;
  lifetimes.resize(captureVect.size());

  for (i=0;i<captureVect.size();++i)
  {
    lifetimes[i] = theReactions[captureVect[i]].getRateVC(concs, constant_vec);
    // inverse of lifetime = rate/concentration, so...
    if (lifetimes[i] > 0) // don't try to compute infinite lifetime!
      lifetimes[i] = concentration/lifetimes[i];
    else
      lifetimes[i]=-1.0; // bollocks, but let's not quibble.   This is just as
                         // invalid a lifetime as infinity would be!
  }

}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDRateDC
// Purpose       : Compute the vector of derivatives of the Rate (as returned
//                 by getRate) with respect to concentrations
//
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+E or A+E=>B and will be incorrect if
//                 the number of any reaction not of this form is listed in
//                 the capture or emission vectors.
//
//                 Generic version in which the list of capture and emission
//                 reactions are passed in
//
//                 Note that dRatedC is zeroed, so if you need to sum
//                 terms from multiple species, you need to pass a temporary
//                 to this routine and sum terms in to a different vector.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/13/06
//-----------------------------------------------------------------------------
void ReactionNetwork::getDRateDC(std::vector<double> &concs,
                                            std::vector<double> &constant_vec,
                                            std::vector<int> &captureVect,
                                            std::vector<int> &emissionVect,
                                            std::vector<double> &dRatedC)
{
  int i,j;
  int cSize=dRatedC.size();
  int const_size=constant_vec.size();
  std::vector<double> tempdRdC(cSize);
  int FADSize = cSize+const_size;

  for (i=0;i<cSize;++i)
  {
    dRatedC[i]=0;
  }

  //The following sets up the FAD types in order to get sensitvities of the capture Rxns
  FDFadType tempdRdCFD;
  std::vector<FDFadType> defects(cSize);
  std::vector<FDFadType> carriers(const_size);
  for( i=0 ; i<const_size ; ++i)
    {
      carriers[i] = constant_vec[i];
      //carriers[i].diff(i,FADSize);
    }
  for( i=0 ; i<cSize ; ++i)
    {
      defects[i] = concs[i];
      defects[i].diff(i+const_size,FADSize);
    }

  for (i=0;i<captureVect.size();++i)
    {
      tempdRdCFD = theReactions[captureVect[i]].getRateVC(defects, carriers);

      for (j=0;j<cSize;++j)
        {
          dRatedC[j] -= tempdRdCFD.dx(j+const_size);
        }
    }

  for (i=0;i<emissionVect.size();++i)
    {
      tempdRdC.assign(cSize,0.0);
      if(theReactions[emissionVect[i]].getCarrierEmissionIndex()<0)
        {
          tempdRdCFD = theReactions[emissionVect[i]].getRateVC(defects, carriers);
          for(int j=0 ; j<cSize ; ++j)
            {
              tempdRdC[j] = tempdRdCFD.dx(j+const_size);
            }
          //theReactions[emissionVect[i]].getDRateDC(concs, constant_vec, tempdRdC);
        }
      else
        {
          tempdRdCFD = theReactions[emissionVect[i]].getFDEmissionRate(defects, carriers);
          for(int j=0 ; j<cSize ; ++j)
            {
              tempdRdC[j] = tempdRdCFD.dx(j+const_size);
            }
        }


      for (j=0;j<cSize;++j)
        {
          dRatedC[j] += tempdRdC[j];
        }
    }
}


//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::getDRateDConst
//
// Purpose       : Compute the vector of derivatives of the Rate (as returned
//                 by getRate) with respect to constants, like E and H.
//
// Special Notes : Assumes that all emission and capture reactions are
//                 of the form B=>A+E or A+E=>B and will be incorrect if
//                 the number of any reaction not of this form is listed in
//                 the capture or emission vectors.
//
//                 Generic version in which the list of capture and emission
//                 reactions are passed in
//
//                 Note that dRatedConst is zeroed, so if you need to sum
//                 terms from multiple species, you need to pass a temporary
//                 to this routine and sum terms in to a different vector.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 11/15/08
//-----------------------------------------------------------------------------
void ReactionNetwork::getDRateDConst(std::vector<double> &concs,
                                                std::vector<double> &constant_vec,
                                                std::vector<int> &captureVect,
                                                std::vector<int> &emissionVect,
                                                std::vector<double> &dRatedConst)
{
  int i,j;
  int cSize=dRatedConst.size();
  int concSize=concs.size();
  int constSize=constant_vec.size();
  std::vector<double> tempdRdConst(cSize);
  int FADSize=concSize+constSize;

  for (i=0;i<cSize;++i)
  {
    dRatedConst[i]=0;
  }

  //The following sets up the FAD types in order to get sensitvities of the FD emission Rxns
  FDFadType tempdRdConstFD;
  std::vector<FDFadType> defects(concs.size());
  std::vector<FDFadType> carriers(2);
  for( i=0 ; i<cSize ; ++i)
    {
      carriers[i] = constant_vec[i];
      carriers[i].diff(i,FADSize);
    }
  for( i=0 ; i<concs.size() ; ++i)
    {
      defects[i] = concs[i];
      //defects[i].diff(i+2,FADSize);
    }

  for (i=0;i<captureVect.size();++i)
  {
    tempdRdConst.assign(cSize,0.0);

    tempdRdConstFD = theReactions[captureVect[i]].getRateVC(defects, carriers);

    for (j=0;j<cSize;++j)
    {
      dRatedConst[j] -= tempdRdConstFD.dx(j);
    }
  }


  for (i=0;i<emissionVect.size();++i)
  {
    tempdRdConst.assign(cSize,0.0);

      if(theReactions[emissionVect[i]].getCarrierEmissionIndex()<0)
        {
          tempdRdConstFD = theReactions[emissionVect[i]].getRateVC(defects, carriers);
          for(int j=0 ; j<cSize ; ++j)
            {
              tempdRdConst[j] = tempdRdConstFD.dx(j);
            }
        }
      else
        {
          tempdRdConstFD = theReactions[emissionVect[i]].getFDEmissionRate(defects, carriers);
          for(int j=0 ; j<cSize ; ++j)
            {
              tempdRdConst[j] = tempdRdConstFD.dx(j);
            }
        }

      for (j=0;j<cSize;++j)
        dRatedConst[j] += tempdRdConst[j];
  }

}



//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setBourgoinCorbettCalc
// Purpose       : set the named Bourgoin Corbett parameters in the species objects
//                 Corbett enhanced diffusion
// Special Notes :
// Scope         : public
// Creator       : Lawrence C Musson, SNL
// Creation Date : 03/25/14
//-----------------------------------------------------------------------------
  void ReactionNetwork::setBourgoinCorbettCalc(const std::string &speciesName,const std::string &carrierName,
                                               double sigma, double hopDistance)
  {



    //int i=species.size(); // get index of next one we add
    //species.push_back(aSpecie);
    //speciesMap[aSpecie.getName()]=i;

    int speciesIndex = speciesMap[speciesName];

    int carrierIndex=-1;

    double thermalVelocity=0.0;


    // for (int i=0 ; i<constantsMap.size() ; ++i)
    //if(carrierName == constants
    carrierIndex = constantsMap[carrierName];

    int BCCarrierCharge=0;

    if(carrierName == "E")
      {
        BCCarrierCharge = -1;
      }
    else if(carrierName == "H")
      {
        BCCarrierCharge = 1;
      }
    else
    {
      Report::DevelFatal() << "ReactionNetwork::setBourgoinCorbettCalc: Illegal carrier for BC enhancement: " << carrierName;
    }

    species[speciesIndex].setBCEnhancedDiffusion(carrierIndex,sigma,BCCarrierCharge,hopDistance);

  }

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setFDElectronEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission with Fermi-Dirac Stats
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 06/30/2014
//-----------------------------------------------------------------------------
void ReactionNetwork::setFDElectronEmissionCalc(const std::string &name, double sigma, double E)
{

  int carrierIndex = 0;

  carrierIndex = constantsMap["E"];

  //carrier velocities and DOS are a function of the bulk device material
  //The reaction network is set up prior to the device and so material
  //properties are unknown when this method is called.  Instead of passing
  //through a velocity, I'm passing through a carrier charge so that
  //when the mulk material is known, this reaction can identify the carrier
  //that is being emitted and the velocity can be set appropriately for
  //the device material. --LCM
  double v = -1.0;
  getReaction(name).setFDEmissionRateCalculator(carrierIndex,sigma,E,v,1.0,t0,x0);
}

//-----------------------------------------------------------------------------
// Function      : ReactionNetwork::setFDHoleEmissionCalc
// Purpose       : set the named reaction's rate calculator to type emission with Fermi-Dirac Stats
// Scope         : public
// Creator       : Lawrence C Musson
// Creation Date : 06/30/2014
//-----------------------------------------------------------------------------
void ReactionNetwork::setFDHoleEmissionCalc(const std::string &name, double sigma, double E)
{

  int carrierIndex = 1;

  carrierIndex = constantsMap["H"];

  //carrier velocities and DOS are a function of the bulk device material
  //The reaction network is set up prior to the device and so material
  //properties are unknown when this method is called.  Instead of passing
  //through a velocity, I'm passing through a carrier charge so that
  //when the mulk material is known, this reaction can identify the carrier
  //that is being emitted and the velocity can be set appropriately for
  //the device material. --LCM
  double v = +1.0;
  getReaction(name).setFDEmissionRateCalculator(carrierIndex,sigma,E,v,1.0,t0,x0);
}


} // namespace Device
} // namespace Xyce

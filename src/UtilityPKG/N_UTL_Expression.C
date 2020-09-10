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

//-------------------------------------------------------------------------
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Eric R. Keiter, SNL
//
// Creation Date : 04/17/08
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_DEV_Const.h>

// ---------- Standard Includes ----------

#include <iterator>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <sstream> 
// ----------   Xyce Includes   ----------
#include <N_UTL_Expression.h>

#include <newExpression.h>
#include <xyceExpressionGroup.h>
#include <mainXyceExpressionGroup.h>
#include <N_UTL_ExpressionInternals.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : Expression::Expression
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
Expression::Expression( 
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & baseGrp_,
    const std::string & exp, 
    const std::vector<std::string> & functionArgStringVec)
  :
   newExpPtr_(NULL),
   grp_(baseGrp_)
{
  newExpPtr_ = Teuchos::rcp(new Xyce::Util::newExpression(exp,grp_) );

  if (!(functionArgStringVec.empty()))
  {
    newExpPtr_->setFunctionArgStringVec(functionArgStringVec);
  }

  newExpPtr_->lexAndParseExpression();
}

//-----------------------------------------------------------------------------
// Function      : Expression::Expression
// Purpose       : Copy Constructor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
Expression::Expression( const Expression & right)
  :
   newExpPtr_(right.newExpPtr_),
   grp_(right.grp_)
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::operator=
// Purpose       : assignment operator
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 1/7/2019
//-----------------------------------------------------------------------------
Expression& Expression::operator=(const Expression& right) 
{
  grp_ = right.grp_;
  newExpPtr_ = right.newExpPtr_;
  //if (*(newExpPtr_->parsed())) newExpPtr_->lexAndParseExpression();
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : Expression::~Expression
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
Expression::~Expression ()
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::parsed
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::parsed() const 
{
  return newExpPtr_->parsed();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getFuncSize
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
int Expression::getFuncSize()
{
  return newExpPtr_->getFuncOpVec().size();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getFuncNames
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::getFuncNames (std::vector<std::string> & funcNames)
{
  funcNames = newExpPtr_-> getFuncNameVec ();
}


//-----------------------------------------------------------------------------
// Function      : Expression::getFuncPrototypeArgStrings
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::getFuncPrototypeArgStrings(std::vector<std::string> & arguments)
{
  newExpPtr_->getFuncPrototypeArgStrings(arguments);
} 

//-----------------------------------------------------------------------------
// Function      : Expression::attachFunctionNode
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::attachFunctionNode (const std::string & funcName, const Expression & exp)
{
  newExpPtr_->attachFunctionNode(funcName,exp.newExpPtr_);
}

//-----------------------------------------------------------------------------
// Function      : Expression::attachParameterNode
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::attachParameterNode (const std::string & paramName, const Expression & exp, enumParamType type)
{
#if 0
  Xyce::dout() << "attachParameterNode name = " << paramName << " which is ";
  if (type==DOT_PARAM)               { Xyce::dout() << "a .param" <<std::endl; }
  else if (type==DOT_GLOBAL_PARAM)   { Xyce::dout() << "a .global_param" <<std::endl; }
  else if ( type== SUBCKT_ARG_PARAM} { Xyce::dout() << "a subcircuit parameter argument" <<std::endl; }
#endif
  newExpPtr_->attachParameterNode(paramName,exp.newExpPtr_, type);
}


//-----------------------------------------------------------------------------
// Function      : Expression::getFunctionArgStringVec
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getFunctionArgStringVec ()
{
  return newExpPtr_->getFunctionArgStringVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_type
// Purpose       : Finds the type of an input quantity name
//
// Special Notes : This is only called from N_TOP_CktNode_Dev.C, in the 
//                 function: CktNode_Dev::getDepSolnVars.  It only needs to 
//                 detect 3 types:  XEXP_NODE, XEXP_INSTANCE and  XEXP_LEAD
//
//                 At the moment it does not detect XEXP_LEAD.  So, not done yet.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::get_type ( const std::string & var )
{
  newExpPtr_->setupVariousAstArrays();

  int retVal=0; 

  std::string tmpName = var;
  Xyce::Util::toUpper(tmpName);

  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & voltMap = newExpPtr_->getVoltOpMap ();
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & currMap = newExpPtr_->getCurrentOpMap ();
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & leadMap = newExpPtr_->getLeadCurrentOpMap ();

  if ( voltMap.find(tmpName) != voltMap.end() )
  {
    retVal = XEXP_NODE;
  }
  else if ( currMap.find(tmpName) != currMap.end() )
  {
    retVal = XEXP_INSTANCE;
  }
  else if ( leadMap.find(tmpName) != leadMap.end() )
  {
    retVal = XEXP_LEAD;
  }
  else
  {
    newExpPtr_->dumpParseTree(Xyce::dout());
    Xyce::dout() << "Error. Xyce::Util::Expression::get_type.  Cannot find type for " << var 
      << " in expression: " << newExpPtr_->getExpressionString()  <<std::endl;
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::make_constant
// Purpose       : Convert a 'string' placeholder into a constant
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::make_constant (const std::string & var, const double & val, enumParamType type)
{
  newExpPtr_->setupVariousAstArrays();

#if 0
  Xyce::dout() << "make_constant name = " << var << " which is ";
  if (type==DOT_PARAM)               { Xyce::dout() << "a .param" <<std::endl; }
  else if (type==DOT_GLOBAL_PARAM)   { Xyce::dout() << "a .global_param" <<std::endl; }
  else if ( type== SUBCKT_ARG_PARAM} { Xyce::dout() << "a subcircuit parameter argument" <<std::endl; }
#endif
  bool retVal=false; // ERK.  check this.
  retVal = newExpPtr_->make_constant (var,val, type);
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::make_var
//
// Purpose       : Convert a 'string' placeholder into a variable
//
// Special Notes : I set up the new expression library to compute derivatives 
//                 with respect to any voltage nodes and also any source 
//                 currents present in the expression.  I did not set up the new
//                 expression library to compute derivatives w.r.t. parameters.
//                 This seemed like a fraught exercise, as params and 
//                 global_params were often just proxies for other, complicated 
//                 expressions.  Why differentiate them? why not differentiation 
//                 something more atomistic?
//
//                 But there was a solution.  The old expresison library also 
//                 only assumed you wanted voltages and currents differentiated.  
//                 There is a "make_var" function, which allows you to tag 
//                 certain variables as needing differentiation.  This is useful
//                 for Xyce::Util::newExpression as well, for obvious reasons.
//
//                 For the reasons given above, make_var can only operate on
//                 paramOp classes.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::make_var (std::string const & var, enumParamType type)
{ 
#if 0
  Xyce::dout() << "mak_var name = " << var << " which is ";
  if (type==DOT_PARAM)               { Xyce::dout() << "a .param" <<std::endl; }
  else if (type==DOT_GLOBAL_PARAM)   { Xyce::dout() << "a .global_param" <<std::endl; }
  else if ( type== SUBCKT_ARG_PARAM} { Xyce::dout() << "a subcircuit parameter argument" <<std::endl; }
#endif
  return newExpPtr_->make_var(var, type);
}

//-----------------------------------------------------------------------------
// Function      : Expression::setGroup
// Purpose       : 
// Special Notes :
// Scope         : 
// Creator       : Eric Keiter, SNL
// Creation Date : 8/22/2020
//-----------------------------------------------------------------------------
void Expression::setGroup( Teuchos::RCP<baseExpressionGroup> & grp )
{
  newExpPtr_->setGroup(grp);
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedParams
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/20/2020
//-----------------------------------------------------------------------------
void Expression::getUnresolvedParams (std::vector<std::string> & params) const
{
  newExpPtr_->setupVariousAstArrays();

#if 0
  newExpPtr_->dumpParseTree(Xyce::dout());
#endif

  params.clear();
#if 0
  std::vector<Teuchos::RCP<astNode<usedType> > > & paramOpVec = newExpPtr_->getParamOpVec();
  for (int ii=0;ii<paramOpVec.size();ii++)
  {
    Teuchos::RCP<paramOp<usedType> > parPtr = Teuchos::rcp_dynamic_cast<paramOp<usedType> > (paramOpVec[ii]);

    if( !(parPtr->getIsConstant())  && !(parPtr->getIsVar())  && !(parPtr->getIsAttached()) ) 
    {
      std::string tmpName = paramOpVec[ii]->getName();
      std::vector<std::string>::iterator it = std::find(params.begin(), params.end(), tmpName);
      if (it == params.end())
      {
        params.push_back( tmpName );
#if 0
        Xyce::dout() << "newExpression::getUnresolvedParams for " << newExpPtr_->getExpressionString() 
          << " pushing back " << tmpName << std::endl;
#endif
      }
    }
  }
#else
  std::vector<std::string> & unresolvedParamNameVec = newExpPtr_->getUnresolvedParamNameVec();
  if (!(unresolvedParamNameVec.empty())) { params.insert(params.end(),unresolvedParamNameVec.begin(), unresolvedParamNameVec.end()); }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedParams
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/9/2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getUnresolvedParams () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getUnresolvedParamNameVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getParams
// Purpose       : 
// Special Notes : ERK: Fix this.  
//                 It is figuring out a unique list every single time, by 
//                 using "find"
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getParams (std::vector<std::string> & params) const
{
  newExpPtr_->setupVariousAstArrays();

  params.clear();
  for (int ii=0;ii<newExpPtr_->getParamOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getParamOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(params.begin(), params.end(), tmpName);
    if (it == params.end())
    {
      params.push_back( tmpName );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getVoltageNodes   (std::vector<std::string> & nodes) const
{
  newExpPtr_->setupVariousAstArrays();

  nodes.clear();
  std::vector<std::string> & voltNames = newExpPtr_->getVoltNameVec ();
  if (!(voltNames.empty())) { nodes.insert(nodes.end(),voltNames.begin(), voltNames.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getVoltageNodes () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getVoltNameVec ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDeviceCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getDeviceCurrents (std::vector<std::string> & devices) const
{
  newExpPtr_->setupVariousAstArrays();

  devices.clear();
  std::vector<std::string> & currentNames = newExpPtr_->getCurrentNameVec ();
  if (!(currentNames.empty())) { devices.insert(devices.end(),currentNames.begin(), currentNames.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDeviceCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getDeviceCurrents () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getCurrentNameVec ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getLeadCurrents (std::vector<std::string> & leads) const
{
  newExpPtr_->setupVariousAstArrays();

  leads.clear();
  for (int ii=0;ii<newExpPtr_->getLeadCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getLeadCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end())
    {
      leads.push_back( tmpName );
    }
  }

  // experiment:
  for (int ii=0;ii<newExpPtr_->getBsrcCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getBsrcCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end())
    {
      leads.push_back( tmpName );
    }
  }
  // experiment:   In at least some cases, what is really being requested is branch calculations, which can be either lead currents or power.
  for (int ii=0;ii<newExpPtr_->getPowerOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getPowerOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end())
    {
      leads.push_back( tmpName );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrentsExcludeBsrc
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getLeadCurrentsExcludeBsrc (std::vector<std::string> & leads) const
{
  newExpPtr_->setupVariousAstArrays(); 

  leads.clear();
  for (int ii=0;ii<newExpPtr_->getLeadCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getLeadCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end())
    {
      leads.push_back( tmpName );
    }
  }

  // In at least some cases, what is really being requested is branch calculations, which can be either lead currents or power.
  for (int ii=0;ii<newExpPtr_->getPowerOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getPowerOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end())
    {
      leads.push_back( tmpName );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getFunctions
// Purpose       : 
//
// Special Notes : ERK: Fix this.  
//                 It is figuring out a unique list every single time, by 
//                 using "find"
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getFunctions (std::vector<std::string> & funcs) const
{
  newExpPtr_->setupVariousAstArrays();

  funcs.clear();
  for (int ii=0;ii<newExpPtr_->getFuncOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getFuncOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(funcs.begin(), funcs.end(), tmpName);
    if (it == funcs.end())
    {
      funcs.push_back( tmpName );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedFunctions
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/23/2020
//-----------------------------------------------------------------------------
void Expression::getUnresolvedFunctions (std::vector<std::string> & funcs) const
{
  newExpPtr_->setupVariousAstArrays();

  funcs.clear();
  std::vector<std::string> & unresolvedFuncNameVec = newExpPtr_->getUnresolvedFuncNameVec();
  if (!(unresolvedFuncNameVec.empty())) { funcs.insert(funcs.end(),unresolvedFuncNameVec.begin(), unresolvedFuncNameVec.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedFunctions
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/23/2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getUnresolvedFunctions () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getUnresolvedFuncNameVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getSpecials
// Purpose       : 
// Special Notes : does this need to catch GMIN as well?
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getSpecials (std::vector<std::string> & specials) const
{
  newExpPtr_->setupVariousAstArrays();

  specials.clear();
  if (newExpPtr_->getTimeDependent()) { specials.push_back(std::string("TIME")); }
  if (newExpPtr_->getTempDependent()) { specials.push_back(std::string("TEMP")); }
  if (newExpPtr_->getVTDependent()) { specials.push_back(std::string("VT")); }
  if (newExpPtr_->getFreqDependent()) { specials.push_back(std::string("FREQ")); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVariables
//
// Purpose       : This function returns the names of parameters in the expression 
//                 which have an ongoing external dependency.    This basically 
//                 means .global_params, and not .params in Xyce. 
//
//                 The param class has several booleans in it, that pertain to 
//                 different aspects of external resolution.  But, most of them 
//                 are unrelated to whether or not they were originally a .param 
//                 or a .global_param.  So, I created another boolean to indicate 
//                 when or not a parameter is a .param.  If it is, then it can't 
//                 be a "variable".  If it is not, then it is considered a "variable".
//
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getVariables (std::vector<std::string> & variables) const
{
  newExpPtr_->setupVariousAstArrays();

  variables.clear();
#if 0
  std::vector<Teuchos::RCP<astNode<usedType> > > & paramOpVec = newExpPtr_->getParamOpVec();
  for (int ii=0;ii<paramOpVec.size();ii++)
  {
    Teuchos::RCP<paramOp<usedType> > parPtr = Teuchos::rcp_dynamic_cast<paramOp<usedType> > (paramOpVec[ii]);

    if (  parPtr->getParamType() == DOT_GLOBAL_PARAM ) 
    {
      std::string tmpName = paramOpVec[ii]->getName();
      std::vector<std::string>::iterator it = std::find(variables.begin(), variables.end(), tmpName);
      if (it == variables.end())
      {
        variables.push_back( tmpName );
      }
    }
  }

#if 1
  std::vector<std::string> & globalParamNameVec = newExpPtr_->getGlobalParamNameVec();
  if (globalParamNameVec.size() != variables.size())
  {
    std::cout << "Expression::getVariables problem!" <<std::endl;
    std::cout << "variables.size = " << variables.size() <<  " globalParamNameVec.size() = " << globalParamNameVec.size() <<std::endl;

    for (int ii=0;ii<variables.size();ii++)
    {
      std::cout << "variables["<<ii<<"] = " << variables[ii] <<std::endl;
    }
    for (int ii=0;ii<globalParamNameVec.size();ii++)
    {
      std::cout << "globalParamNameVec["<<ii<<"] = " << globalParamNameVec[ii] <<std::endl;
    }

    exit(0);
  }
#endif

#else
  std::vector<std::string> & globalParamNameVec = newExpPtr_->getGlobalParamNameVec();
  if (!(globalParamNameVec.empty())) { variables.insert(variables.end(),globalParamNameVec.begin(), globalParamNameVec.end()); }
#endif

#if 0
  if ( !(variables.empty()) )
  {
    Xyce::dout() << "Expression::getVariables call for " << newExpPtr_->getExpressionString() << std::endl;
    for (int ii=0;ii<variables.size();ii++) { Xyce::dout() << variables[ii] << std::endl; }
    newExpPtr_->dumpParseTree(Xyce::dout());
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVariables
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getVariables() const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getGlobalParamNameVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getPowerCalcs
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getPowerCalcs       (std::vector<std::string> & powerCalcs) const
{
  newExpPtr_->setupVariousAstArrays();

  powerCalcs.clear();
  for (int ii=0;ii<newExpPtr_->getPowerOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getPowerOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(powerCalcs.begin(), powerCalcs.end(), tmpName);
    if (it == powerCalcs.end())
    {
      powerCalcs.push_back( tmpName );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getEverything
// Purpose       : 
// Special Notes : ERK: Fix this.  
//                 It is figuring out a unique list every single time, by 
//                 using "find"
//
//                 This is an attempt to squeeze out some extra performance.
//                 In the IO package, there are several places in which 
//                 these 5 vectors of strings are requested.  They were requested 
//                 via 5 separate function calls, but that seemed wasteful.  This
//                 function combines them.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getEverything (
    std::vector<std::string> & nodes,
    std::vector<std::string> & devices,
    std::vector<std::string> & leads,
    std::vector<std::string> & variables,
    std::vector<std::string> & specials 
    ) const
{
  newExpPtr_->setupVariousAstArrays();

  nodes.clear(); devices.clear(); leads.clear(); variables.clear(); specials.clear();

  // voltage nodes
  std::vector<std::string> & voltNames = newExpPtr_->getVoltNameVec ();
  if (!(voltNames.empty())) { nodes.insert(nodes.end(),voltNames.begin(), voltNames.end()); }

  // solution current
  std::vector<std::string> & currentNames = newExpPtr_->getCurrentNameVec ();
  if (!(currentNames.empty())) { devices.insert(devices.end(),currentNames.begin(), currentNames.end()); }

  // lead currents
  for (int ii=0;ii<newExpPtr_->getLeadCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getLeadCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end()) { leads.push_back( tmpName ); }
  }

  // more lead currents
  for (int ii=0;ii<newExpPtr_->getBsrcCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getBsrcCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end()) { leads.push_back( tmpName ); }
  }

  // more, more lead currents:   In at least some cases, what is really being requested is 
  // branch calculations, which can be either lead currents or power.
  for (int ii=0;ii<newExpPtr_->getPowerOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getPowerOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(leads.begin(), leads.end(), tmpName);
    if (it == leads.end()) { leads.push_back( tmpName ); }
  }

  //variables (global params):
  for (int ii=0;ii<newExpPtr_->getParamOpVec().size();ii++)
  {
    Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (newExpPtr_->getParamOpVec()[ii]);
    if (  parOp->getParamType() == DOT_GLOBAL_PARAM ) 
    {
      std::string tmpName = newExpPtr_->getParamOpVec()[ii]->getName();
      std::vector<std::string>::iterator it = std::find(variables.begin(), variables.end(), tmpName);
      if (it == variables.end()) { variables.push_back( tmpName ); }
    }
  }

  // specials:
  if (newExpPtr_->getTimeDependent()) { specials.push_back(std::string("TIME")); }
  if (newExpPtr_->getTempDependent()) { specials.push_back(std::string("TEMP")); }
  if (newExpPtr_->getVTDependent()) { specials.push_back(std::string("VT")); }
  if (newExpPtr_->getFreqDependent()) { specials.push_back(std::string("FREQ")); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getIsConstant
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getIsConstant ()
{
  return newExpPtr_->getIsConstant();
}

//-----------------------------------------------------------------------------
// Function      : Expression::setTemperature
// Purpose       :
// Special Notes : This is ONLY called when you want to override the
//                 'circuit' temperature.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/19/2020
//-----------------------------------------------------------------------------
bool Expression::setTemperature   (const double & temp)
{
  return newExpPtr_->setTemperature(temp);
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_expression
// Purpose       : Returns a string of the expression, post replacements
//
//                 This is primarily used as a diagnostic.  
//
// Special Notes : In the old expression library, the returned string was 
//                 potentially modified, as part of its process for resolving
//                 functions and parameters.  So it wasn't necessarily the same
//                 string as was originally passed into the expression library.
//
//                 In the new expression library that is not the case.  In the
//                 new expression library functions and parameters are not 
//                 resolved via string substitution.  They are resolved by 
//                 attaching nodes to the AST.  So nowadays, this function
//                 returns 100% the same thing as the function "get_input".
//
//                 There maybe a few use cases where this would be useful 
//                 in the newExpression library.  Adding a feature to it where
//                 it created a new expression string from the AST would not 
//                 be too hard; it would be akin to what is already done for 
//                 dumping the expression tree but in a more compact form.
//
//                 One use case where it would still be useful would be for 
//                 params that are *not* attached (like const params), 
//                 and also for node aliases. 
//
//                 But setting this up is not a high priority at the moment.
//
//                 The function should probably have a more descriptive name.
//                 From the name alone, it was hard for me to understand the
//                 difference between this and get_input. (which is part of why
//                 they are for now equivalent)
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
std::string Expression::get_expression () const
{
  return newExpPtr_->getExpressionString();
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluate
// Purpose       : Evaluate expression and derivatives using stored input values
// Special Notes : 
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluate ( std::complex<double> & exp_r, std::vector< std::complex<double> > & deriv_r)
{
  bool retVal=true;
#ifdef USE_TYPE_DOUBLE
  double result;
  std::vector<double> derivs;
  retVal = newExpPtr_->evaluate( result, derivs );
  exp_r = std::complex<double>(result,0.0);
  deriv_r.resize(derivs.size(),0.0);
  for(int ii=0;ii<derivs.size();ii++) { deriv_r[ii] = std::complex<double>(derivs[ii],0.0); } // could use a lambda here
#else
  retVal = newExpPtr_->evaluate( exp_r, deriv_r );
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluateFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluateFunction ( std::complex<double> & exp_r, bool efficiencyOn )
{
  bool retVal=true; 
#ifdef USE_TYPE_DOUBLE
  double result;
  retVal = newExpPtr_->evaluateFunction( result, efficiencyOn );
  exp_r = std::complex<double>(result,0.0);
#else
  retVal = newExpPtr_->evaluateFunction ( exp_r, efficiencyOn );
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluate
// Purpose       : Evaluate expression and derivatives using stored input values
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluate ( double & exp_r, std::vector<double> & deriv_r)
{
  bool retVal=true;
#ifdef USE_TYPE_DOUBLE
  retVal = newExpPtr_->evaluate( exp_r, deriv_r );
#else
  std::complex<double> result;
  std::vector<std::complex<double> > derivs;
  retVal = newExpPtr_->evaluate( result, derivs );

  exp_r = std::real(result);
  deriv_r.resize(derivs.size(),0.0);
  for(int ii=0;ii<derivs.size();ii++) {  deriv_r[ii] = std::real(derivs[ii]); } // could use a lambda here
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluateFunction
// Purpose       : Evaluate expression using stored input values
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluateFunction ( double & exp_r, bool efficiencyOn )
{
  bool retVal=true; 
#ifdef USE_TYPE_DOUBLE
  retVal = newExpPtr_->evaluateFunction ( exp_r, efficiencyOn );
#else
  std::complex<double> result;
  retVal = newExpPtr_->evaluateFunction ( result, efficiencyOn );
  exp_r = std::real(result);
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::clearOldResult
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 09/04/2020
//-----------------------------------------------------------------------------
void Expression::clearOldResult ()
{
  return newExpPtr_->clearOldResult();
}

//-----------------------------------------------------------------------------
// Function      : Expression::setupBreakPoints
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 08/25/2020
//-----------------------------------------------------------------------------
void Expression::setupBreakPoints()
{
  return newExpPtr_->setupBreakPoints();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getBreakPoints
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::getBreakPoints(std::vector<Util::BreakPoint> &breakPointTimes)
{
  return newExpPtr_->getBreakPoints(breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_input
// Purpose       : Return expression input string
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
const std::string & Expression::get_input () const
{
  return newExpPtr_->getExpressionString();
}

//-----------------------------------------------------------------------------
// Function      : Expression::replace_name
// Purpose       : Change the name of an input quantity
//
// Special Notes : ERK.  This called from the N_IO_DistToolBase.C file/class
//                 in the function DistToolBase::instantiateDevice.
//
//                 "Input quantity" in this case means voltage nodes (XEXP_NODE), 
//                 device instances (XEXP_INSTANCE), and lead currents (XEXP_LEAD).
//
//                 Sometimes, they are specified in expressions without their 
//                 names being fully resolved.  ie, the expression is inside of 
//                 a subcircuit, and thus implicitly assumes the full prefix.
//
//                 So, this function adds the full prefix to these names, 
//                 so they can be fully resolved.
//
//                 The function DistToolBase::instantiateDevice calls this 
//                 function twice for some reason that I don't (yet) understand.
//                 It seems to require 2 passes to properly update the name. ie,
//                 "name" -> "; name" -> "prefix:name".
//  
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::replace_name ( const std::string & old_name,
                                const std::string & new_name)
{
  return newExpPtr_->replaceName( old_name, new_name );
}

//-----------------------------------------------------------------------------
// Function      : Expression::isTimeDependent
// Purpose       : Return true if expression is either explicitly OR implicitly
//                 time dependent
// Special Notes : 
// Scope         :
// Creator       : 
// Creation Date : 10/07/2013
//-----------------------------------------------------------------------------
bool Expression::isTimeDependent() const
{
  // ERK. This is important so that (for example) capacitors with expression 
  // dependent capacitance aren't calling updateTemperature over and over 
  // again.   If they think an expression dependent parameter is time 
  // dependent, then they are obligated to keep re-evaluating it.
  // So, this needs to give the right answer.
  bool explicitTimeDep = newExpPtr_->getTimeDependent();
  bool sdtDep = !(newExpPtr_->getSdtOpVec().empty());
  bool ddtDep = !(newExpPtr_->getDdtOpVec().empty());
  return (explicitTimeDep || sdtDep || ddtDep);
}

//-----------------------------------------------------------------------------
// Function      : Expression::isFreqDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isFreqDependent() const
{
  return newExpPtr_->getFreqDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::isSolutionDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isSolutionDependent() const
{
  return ( !(newExpPtr_->getVoltOpVec().empty()) || 
           !(newExpPtr_->getCurrentOpVec().empty()) );
}

//-----------------------------------------------------------------------------
// Function      : Expression::isRandomDependent
// Purpose       : Return true if expression dependent on GAUSS, AGAUSS or RAND
// Special Notes : This is only based on local dependence, from parsing.
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isRandomDependent() const
{
  if ( !(newExpPtr_->getLocalAgaussOpVec().empty())  ) { return true; }
  if ( !(newExpPtr_->getLocalGaussOpVec().empty()) ) { return true; }
  if ( !(newExpPtr_->getLocalAunifOpVec().empty()) ) { return true; }
  if ( !(newExpPtr_->getLocalUnifOpVec().empty())  ) { return true; }
  if ( !(newExpPtr_->getLocalRandOpVec().empty())  ) { return true; }
  if ( !(newExpPtr_->getLocalTwoArgLimitOpVec().empty()) ) { return true; }

  return false;
}

//-----------------------------------------------------------------------------
// Function      : Expression::dumpParseTree
// Purpose       : Dump out the parse tree for an expression
// Special Notes : 
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 9/9/10
//-----------------------------------------------------------------------------
void Expression::dumpParseTree()
{
  newExpPtr_->dumpParseTree(Xyce::dout());
}

//-----------------------------------------------------------------------------
// Function      : Expression::seedRandom
// Purpose       : Public method to initialize random number generator
//                 used by rand(), gauss() and agauss() functions.
// Special Notes : 
// Scope         :
// Creator       : Tom Russo
// Creation Date : 03/28/17
//-----------------------------------------------------------------------------
///
/// Public method to initialize random number generator
///
/// @author Tom Russo
/// @date 03/28/17
///
void Expression::seedRandom(long seed)
{
  //if(useNewExpressionLibrary_)
  //{
  //}
  //else
  {
    //ExpressionInternals::seedRandom(seed);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::treatAsTempAndConvert
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::treatAsTempAndConvert()
{
  newExpPtr_->treatAsTempAndConvert();
}

//-----------------------------------------------------------------------------
// Function      : Expression::clearProcessSuccessfulTimeStepMap
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::clearProcessSuccessfulTimeStepMap ()
{
  newExpression::clearProcessSuccessfulTimeStepMap ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::processSuccessfulTimeStep
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::processSuccessfulTimeStep ()
{
  newExpPtr_->processSuccessfulTimeStep ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getNumDdt
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
int Expression::getNumDdt()
{
  return newExpPtr_->getNumDdt();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDdtVals
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getDdtVals (std::vector<double> & vals)
{
  return newExpPtr_->getDdtVals(vals);
}

//-----------------------------------------------------------------------------
// Function      : Expression::setDdtDerivs
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setDdtDerivs (std::vector<double> & vals)
{
  return newExpPtr_->setDdtDerivs(vals);
}

// random operator information
//-----------------------------------------------------------------------------
// Function      : Expression::getAgaussData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getAgaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAgauss = newExpPtr_->getLocalAgaussOpVec();
  for (int ii=0;ii<localAgauss.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<agaussOp<usedType> > agOp = Teuchos::rcp_static_cast<agaussOp<usedType> > (localAgauss[ii]);

    usedType mu    = agOp->getMu();
    usedType alpha = agOp->getAlpha();
    usedType n     = agOp->getN();

    sampling_param.opName     = "AGAUSS";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_AGAUSS;
    sampling_param.type       = "NORMAL";
    sampling_param.mean       = std::real(mu);
    sampling_param.stdDev     = std::real(alpha)/std::real(n);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getGaussData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getGaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localGauss = newExpPtr_->getLocalGaussOpVec();
  for (int ii=0;ii<localGauss.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<gaussOp<usedType> > gaOp = Teuchos::rcp_static_cast<gaussOp<usedType> > (localGauss[ii]);

    usedType mu        = gaOp->getMu();
    usedType rel_alpha = gaOp->getAlpha();
    usedType n         = gaOp->getN();
    usedType alpha     = rel_alpha*mu;

    sampling_param.opName     = "GAUSS";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_GAUSS;
    sampling_param.type       = "NORMAL";
    sampling_param.mean       = std::real(mu);
    sampling_param.stdDev     = std::real(alpha)/std::real(n);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getAunifData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getAunifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAunif = newExpPtr_->getLocalAunifOpVec();
  for (int ii=0;ii<localAunif.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<aunifOp<usedType> > auOp = Teuchos::rcp_static_cast<aunifOp<usedType> > (localAunif[ii]);

    usedType mu    = auOp->getMu();
    usedType alpha = auOp->getAlpha();

    sampling_param.opName     = "AUNIF";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_AUNIF;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = std::real(mu)-std::real(alpha);
    sampling_param.stopVal  = std::real(mu)+std::real(alpha);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnifData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getUnifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localUnif = newExpPtr_->getLocalUnifOpVec();
  for (int ii=0;ii<localUnif.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<unifOp<usedType> > unOp = Teuchos::rcp_static_cast<unifOp<usedType> > (localUnif[ii]);

    usedType mu        = unOp->getMu();
    usedType rel_alpha = unOp->getAlpha();
    usedType alpha     = rel_alpha*mu;

    sampling_param.opName     = "UNIF";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_UNIF;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = std::real(mu)-std::real(alpha);
    sampling_param.stopVal  = std::real(mu)+std::real(alpha);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getRandData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getRandData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localRand = newExpPtr_->getLocalRandOpVec();
  for (int ii=0;ii<localRand.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<randOp<usedType> > auOp = Teuchos::rcp_static_cast<randOp<usedType> > (localRand[ii]);

    sampling_param.opName     = "RAND";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_RAND;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = 0.0;
    sampling_param.stopVal  = 1.0;

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLimitData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getLimitData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localTwoArgLimit = newExpPtr_->getLocalTwoArgLimitOpVec();
  for (int ii=0;ii<localTwoArgLimit.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<twoArgLimitOp<usedType> > auOp = Teuchos::rcp_static_cast<twoArgLimitOp<usedType> > (localTwoArgLimit[ii]);

    sampling_param.opName     = "LIMIT";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_LIMIT;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = 0.0;
    sampling_param.stopVal  = 1.0;

    sampleVec.push_back(sampling_param);
  }
}


} // namespace Util
} // namespace Xyce


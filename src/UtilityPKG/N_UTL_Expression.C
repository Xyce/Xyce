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
  // ERK; removing the beginning and ending brace should really be handled by flex/bison, 
  // but I was in a hurry today.
  std::string expCopy = exp;
  if ( !(expCopy.empty()))
  {
    if (expCopy[0]== '{' && expCopy[expCopy.size()-1]=='}')
    {
      expCopy.erase(0,1);// lop off open curly brace
      expCopy.erase(expCopy.length()-1); // lop off close curly brace
    }
  }

  newExpPtr_ = Teuchos::rcp(new Xyce::Util::newExpression(expCopy,grp_) );

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
void Expression::attachParameterNode (const std::string & paramName, const Expression & exp, bool isDotParam)
{
#if 0
  std::cout << "attachParameterNode name = " << paramName << " which is ";
  if (isDotParam) { std::cout << "a .param" <<std::endl; }
    else          { std::cout << "a not .param" <<std::endl; }
#endif
  newExpPtr_->attachParameterNode(paramName,exp.newExpPtr_, isDotParam);
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
// Function      : Expression::set
// Purpose       : Set the value of the expression to a string
// Special Notes : ERK: is this needed?  Does anyone call it?
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::set ( const std::string & exp )
{
  bool retVal = false; 

  std::string expCopy = exp;

  if ( !(expCopy.empty()))
  {
    if (expCopy[0]== '{' && expCopy[expCopy.size()-1]=='}')
    {
      expCopy.erase(0,1);// lop off open curly brace
      expCopy.erase(expCopy.length()-1); // lop off close curly brace
    }
  }

  if ( !(Teuchos::is_null(newExpPtr_)) )
  {
    newExpPtr_->clear();
    newExpPtr_->setExpressionString (expCopy);
  }
  else
  {
    newExpPtr_ = Teuchos::rcp(new Xyce::Util::newExpression(expCopy, grp_) );
  }
  newExpPtr_->lexAndParseExpression();

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::getSymbolTable
// Purpose       : Returns the symbol table
// Special Notes : ERK.  The need for this may go away by the time I'm done implementing newExpression
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 08/19/2016
//-----------------------------------------------------------------------------
void Expression::getSymbolTable(std::vector< ExpressionSymbolTableEntry > & theSymbolTable ) const
{ 
  //ERK.  Like with everything else, stuff that is old-expression specific I am 
  //handling here, rather than in Xyce::Util::newExpression.

  //  local versions of these
  //  set them up, and then go thru a similar loop as the getSymbolTable function in N_UTL_ExpressionInternals
  std::vector<int> varTypes;            ///< array of types of variables
  std::vector<std::string> varValues;   ///< array of values of variables
  std::string leadDesignator;           ///< lead designator for current variables

  theSymbolTable.clear();

  int var_type = XEXP_NODE;
  for (int ii=0;ii<newExpPtr_->getVoltOpVec().size();ii++)
  {
    int size = newExpPtr_->getVoltOpVec()[ii]->getNodeNames().size();

    for (int jj=0;jj<size;jj++)
    {
      std::string tmpName = newExpPtr_->getVoltOpVec()[ii]->getNodeNames()[jj] ;
      std::vector<std::string>::iterator it = std::find(varValues.begin(), varValues.end(), tmpName);
      if (it == varValues.end())
      {
        varValues.push_back( tmpName );
        varTypes.push_back(var_type);
        leadDesignator.push_back(' ');
      }
    }
  }

  var_type = XEXP_INSTANCE;
  for (int ii=0;ii<newExpPtr_->getCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(varValues.begin(), varValues.end(), tmpName);
    if (it == varValues.end())
    {
      varValues.push_back( tmpName );
      varTypes.push_back(var_type);
      leadDesignator.push_back(' ');
    }
  }

  var_type = XEXP_LEAD; // ERK.  I haven't figured this out yet, but need to.


  var_type = XEXP_STRING; // for some mysterious reason, this means params and global_params
  for (int ii=0;ii<newExpPtr_->getParamOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getParamOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(varValues.begin(), varValues.end(), tmpName);
    if (it == varValues.end())
    {
      varValues.push_back( tmpName );
      varTypes.push_back(var_type);
      leadDesignator.push_back(' ');
    }
  }

  var_type = XEXP_SPECIAL; // ERK.   Does this need to track gmin?
  if (newExpPtr_->getTimeDependent()) { varValues.push_back(std::string("TIME")); varTypes.push_back(var_type); leadDesignator.push_back(' '); }
  if (newExpPtr_->getTempDependent()) { varValues.push_back(std::string("TEMP")); varTypes.push_back(var_type); leadDesignator.push_back(' '); }
  if (newExpPtr_->getVTDependent()) { varValues.push_back(std::string("VT")); varTypes.push_back(var_type); leadDesignator.push_back(' '); }
  if (newExpPtr_->getFreqDependent()) { varValues.push_back(std::string("FREQ")); varTypes.push_back(var_type); leadDesignator.push_back(' '); }

  // VARIABLE is a global parameter that needs to be updated at each call.
  var_type = XEXP_VARIABLE;
#if 0
  Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
  const std::vector<std::string> & names = xyceGroup->getNames();
  int nameSize = names.size();
  for (int ii=0;ii<nameSize;++ii)
  {
    varValues.push_back(names[ii]);
    varTypes.push_back(var_type);
    leadDesignator.push_back(' ');
  }
#endif

  var_type = XEXP_FUNCTION;
  for (int ii=0;ii<newExpPtr_->getFuncOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getFuncOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(varValues.begin(), varValues.end(), tmpName);
    if (it == varValues.end())
    {
      varValues.push_back( tmpName );
      varTypes.push_back(var_type);
      leadDesignator.push_back(' ');
    }
  }

  var_type = XEXP_NODAL_COMPUTATION; // haven't figure this out yet

  int size = varValues.size();
  for (int ii=0;ii<size;++ii)
  {
    theSymbolTable.push_back(ExpressionSymbolTableEntry(varValues[ii],varTypes[ii],leadDesignator[ii]));
  }
  
  return;
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
  int retVal=0; 

  std::string tmpName = var;
  Xyce::Util::toUpper(tmpName);

  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & voltMap = newExpPtr_->getVoltOpNames ();
  const std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & currMap = newExpPtr_->getCurrentOpNames ();

  if ( voltMap.find(tmpName) != voltMap.end() )
  {
    retVal = XEXP_NODE;
  }
  else if ( currMap.find(tmpName) != currMap.end() )
  {
    retVal = XEXP_INSTANCE;
  }
  else
  {
    newExpPtr_->dumpParseTree(std::cout);
    std::cout << "Error. Xyce::Util::Expression::get_type.  Cannot find type for " << var << std::endl;
  }

  //std::cout << "Expression::get_type ( const std::string & var ) " << std::endl;
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
bool Expression::make_constant (const std::string & var, const double & val, bool isDotParam)
{
#if 0
  std::cout << "make_constant name = " << var << " which is ";
  if (isDotParam) { std::cout << "a .param" <<std::endl; }
    else          { std::cout << "a not .param" <<std::endl; }
#endif
  bool retVal=false; // ERK.  check this.
  retVal = newExpPtr_->make_constant (var,val, isDotParam);
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
bool Expression::make_var (std::string const & var, bool isDotParam)
{ 
#if 0
  std::cout << "mak_var name = " << var << " which is ";
  if (isDotParam) { std::cout << "a .param" <<std::endl; }
    else          { std::cout << "a not .param" <<std::endl; }
#endif
  bool retVal=false; 
  retVal = newExpPtr_->make_var(var, isDotParam);
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedParams
// Purpose       : 
// Special Notes : ERK: Fix this.  
//                 It is figuring out a unique list every single time, by 
//                 using "find"
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/20/2020
//-----------------------------------------------------------------------------
void Expression::getUnresolvedParams (std::vector<std::string> & params) const
{
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
      }
    }
  }
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
// Special Notes : ERK: Fix this.  
//                 It is figuring out a unique list every single time, by 
//                 using "find"
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getVoltageNodes   (std::vector<std::string> & nodes) const
{
  for (int ii=0;ii<newExpPtr_->getVoltOpVec().size();ii++)
  {
    int size = newExpPtr_->getVoltOpVec()[ii]->getNodeNames().size();

    for (int jj=0;jj<size;jj++)
    {
      std::string tmpName = newExpPtr_->getVoltOpVec()[ii]->getNodeNames()[jj] ;
      std::vector<std::string>::iterator it = std::find(nodes.begin(), nodes.end(), tmpName);
      if (it == nodes.end())
      {
        nodes.push_back( tmpName );
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDeviceCurrents
// Purpose       : 
// Special Notes : ERK: Fix this.  
//                 It is figuring out a unique list every single time, by 
//                 using "find"
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getDeviceCurrents (std::vector<std::string> & devices) const
{
  for (int ii=0;ii<newExpPtr_->getCurrentOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getCurrentOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(devices.begin(), devices.end(), tmpName);
    if (it == devices.end())
    {
      devices.push_back( tmpName );
    }
  }
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
  //params = newExpPtr_->
  //std::cout << "Error. Xyce::Util::Expression::getLeadCurrents not yet implemented." <<std::endl;
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
// Function      : Expression::getSpecials
// Purpose       : 
// Special Notes : does this need to catch GMIN as well?
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getSpecials (std::vector<std::string> & specials) const
{
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
  for (int ii=0;ii<newExpPtr_->getParamOpVec().size();ii++)
  {
    if (  !(newExpPtr_->getParamOpVec()[ii]->getIsDotParam ())  )
    {
      std::string tmpName = newExpPtr_->getParamOpVec()[ii]->getName();
      std::vector<std::string>::iterator it = std::find(variables.begin(), variables.end(), tmpName);
      if (it == variables.end())
      {
        variables.push_back( tmpName );
      }
    }
  }
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
// Function      : Expression::getNodalCalcs
// Purpose       : 
// Special Notes : this may not be necessary for new expression.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getNodalComputation (std::vector<std::string> & nodalCalcs) const
{

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
int Expression::evaluate ( std::complex<double> & exp_r, std::vector< std::complex<double> > & deriv_r)
{
  int retVal=0;
#ifdef USE_TYPE_COMPLEX
  retVal = newExpPtr_->evaluate( exp_r, deriv_r );
#else
  double result;
  std::vector<double> derivs;
  retVal = newExpPtr_->evaluate( result, derivs );
  exp_r = std::complex<double>(result,0.0);
  deriv_r.resize(derivs.size(),0.0);
  for(int ii=0;ii<derivs.size();ii++) { deriv_r[ii] = std::complex<double>(derivs[ii],0.0); } // could use a lambda here
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
int Expression::evaluateFunction ( std::complex<double> & exp_r )
{
  int retVal=0; 
#ifdef USE_TYPE_COMPLEX
  retVal = newExpPtr_->evaluateFunction ( exp_r );
#else
  double result;
  retVal = newExpPtr_->evaluateFunction( result );
  exp_r = std::complex<double>(result,0.0);
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
int Expression::evaluate ( double & exp_r, std::vector<double> & deriv_r)
{
  int retVal=0;
#ifdef USE_TYPE_COMPLEX
  std::complex<double> result;
  std::vector<std::complex<double> > derivs;
  retVal = newExpPtr_->evaluate( result, derivs );

  exp_r = std::real(result);
  deriv_r.resize(derivs.size(),0.0);
  for(int ii=0;ii<derivs.size();ii++) {  deriv_r[ii] = std::real(derivs[ii]); } // could use a lambda here
#else
  retVal = newExpPtr_->evaluate( exp_r, deriv_r );
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
int Expression::evaluateFunction ( double & exp_r )
{
  int retVal=0; 
#ifdef USE_TYPE_COMPLEX
  std::complex<double> result;
  retVal = newExpPtr_->evaluateFunction ( result );
  exp_r = std::real(result);
#else
  retVal = newExpPtr_->evaluateFunction ( exp_r );
#endif
  return retVal;
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
  bool retVal=false; 
  //std::cout << "NOTE:  replace_name just got called on " << old_name << " to now be " << new_name <<std::endl;

  bool found=false;
  {
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & voltMap = newExpPtr_->getVoltOpNames ();

    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > >::iterator iter = voltMap.find(old_name);

    if (iter != voltMap.end())
    {
      std::vector<Teuchos::RCP<astNode<usedType> > > & astVec = iter->second;

      for(int ii=0;ii<astVec.size();++ii)
      {
        Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (astVec[ii]);
        std::vector<std::string> & nodes = voltOp->getVoltageNodes();
        for(int jj=0;jj<nodes.size();++jj)
        {
          if(nodes[jj]==old_name)
          {
            nodes[jj] = new_name;
          }
        }
      }
      voltMap[new_name] = astVec;
      voltMap.erase(old_name);
      found=true;
    }
  }

  if(!found)
  {
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & currMap = newExpPtr_->getCurrentOpNames ();
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > >::iterator iter = currMap.find(old_name);

    if (iter != currMap.end())
    {
      std::vector<Teuchos::RCP<astNode<usedType> > > & astVec = iter->second;

      for(int ii=0;ii<astVec.size();++ii)
      {
        Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (astVec[ii]);
        currOp->setCurrentDevice(new_name);
      }

      currMap[new_name] = astVec;
      currMap.erase(old_name);
      found=true;
    }
  }
  
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::isTimeDependent
// Purpose       : Return true if expression is either explicitly or implicitly
//                 time dependent
// Special Notes : The ExpressionInternals::isTimeDependent method only returns
//                 true if the expression is implicitly time dependent.
//
//                 It is impossible at this time for indirect time dependence
//                 through global_params to be detected through this method.
//
// Scope         :
// Creator       : Richard Schiek, SNL
// Creation Date : 10/07/2013
//-----------------------------------------------------------------------------
bool Expression::isTimeDependent() const
{
  // ERK. Not done.  Probably do need this.
  //return false;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Expression::isRandomDependent
// Purpose       : Return true if expression dependent on GAUSS, AGAUSS or RAND
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isRandomDependent() const
{
  // ERK. Not done.  Probably do need this.
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
  newExpPtr_->dumpParseTree(std::cout);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::seedRandom
// Purpose       : Public method to initialize random number generator
//                 used by rand(), gauss() and agauss() functions.
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

} // namespace Util
} // namespace Xyce


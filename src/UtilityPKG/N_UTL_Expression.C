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
void Expression::attachParameterNode (const std::string & paramName, const Expression & exp)
{
  newExpPtr_->attachParameterNode(paramName,exp.newExpPtr_);
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
#if 0
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp(new xyceExpressionGroup() );
    grp_ = xyceGroup;
#endif
    //newExpPtr_ = new Xyce::Util::newExpression(expCopy, grp_);
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

  var_type = XEXP_SPECIAL; // ERK.  This doesn't yet track external specials dependencies
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
// Function      : Expression::get_names
//
// Purpose       : This function returns the names of various entities present 
//                 in a parsed expression by type.
//
// Special Notes : These notes pertain to the Xyce::Util::newExpression implementation.
//
//                 This function, which is a holdover of the old expression 
//                 library, is vaguely named.  Its name, "get_names" begs the 
//                 question, "names of what?"
//
//                 I believe that this function returns a vector of strings, 
//                 and the strings are all names of different types of entities 
//                 which may or may not be present in the expression.  These are
//                 generally entities which will need some kind of external 
//                 information, from outside the expression class, 
//                 to be evaluated.
//
//                 So, to use the old expression library, the calling code would 
//                 first call this function, to get the list of named entities.  
//                 It would then see if it could resolve them all.  For example, 
//                 a named entity of type XEXP_NODE refers to a voltage node,
//                 and so the next logical thing to do for those entities 
//                 would be to see if that node actually existed in the 
//                 circuit, and what machinery was needed to obtain its value.
//
//                 Similar thing with XEXP_FUNCTION, for .FUNC references, etc.
//
//                 The full list of types is as follows.
//
//                 enum XEXP_TYPES
//                 {
//                   XEXP_ALL,            // 0    everything
//                   XEXP_NODE,           // 1    nodal variables, does this include currents? (no)
//                   XEXP_INSTANCE,       // 2    current variables, from voltage sources, current sources
//                   XEXP_LEAD,           // 3    current variables, but not from solution vec. ( like I(R1) )
//                   XEXP_STRING,         // 4    as-yet unresolved strings.  Once preliminary resolution is done, this is .func arguments.  
//                   XEXP_SPECIAL,        // 5    special returns vars like TIME
//                   XEXP_VARIABLE,       // 6    also global params, apparently.  Also, it set_var targets?
//                   XEXP_FUNCTION,       // 7    these are .funcs
//                   XEXP_NODAL_COMPUTATION, // 8  these are things like power P()
//                   XEXP_COUNT     // total
//                 };
//
//                 Note, the implied order of "variables" is the order of the types in the enum.
//
//                 The old expression library expects that *all* of these entities will 
//                 be passed into it when the function "evaluate" or "evaluateFunction" 
//                 is called, and they will be passed in as a std::vector<double>.   This
//                 vector will include all of the above things, if they exist.  So, they are
//                 not passed in via separate calls, or in separate containers.  That is, 
//                 the global parameter values will be passed in along with the voltage 
//                 and current values, etc, all in the same vector.  It is up to the user to
//                 establish the order, using the order_names function.
//
//                 I don't like this design.  So I am attempting to have as much of 
//                 it as possible be handled directly in the interface class 
//                 (N_UTL_Expression), rather than in this class.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void Expression::get_names(int const & type, std::vector<std::string> & names ) const
{
  switch (type)
  {
    case XEXP_ALL:  // ERK.  I don't think this one gets called by anyone
      break;

    case XEXP_NODE:
      getVoltageNodes(names);
      break;

    case XEXP_INSTANCE:
      getDeviceCurrents(names);
      break;

    case XEXP_LEAD: // ERK.  I haven't figured this out yet, but need to.
      getLeadCurrents(names);
      break;

    case XEXP_STRING: // unresolved strings.  
      // This is called in a few use cases:
      //
      // (1) to obtain the function arguments specified in the function definition.  
      // ie, if you have:
      //
      // .func abc(x,y)
      //
      // the old Xyce code creates an expression (probably called "functionPrototype") from 
      // the string "abc(x,y)" and then requests the "strings" back, which will be x,y.  
      // For it to work properly, the string vector needs to be in the same order as 
      // they were specified in the prototype.
      //
      // (2) to obtain what are probably function arguments in a function body.  
      // ie, if you have:
      //
      // .param a=2.0
      // .func abc(x,y) {x+y+5*a}
      //
      // Then the function body is {x+y+5*a}.  A "resolution" will figure out that "a" 
      // is a .param, and mark it accordingly, but it will NOT find x and y, as they are 
      // not .params or .global_params.  As a result, x and y will still be in the list 
      // of "strings". (I think.  check this).  This seems a bit backward - the code should
      // already know that x,y are the arguments, if it has already executed use case (1), 
      // above.
      //
      // Note, for this to work, it cannot return any param names (strings) that have 
      // previously had "set_constant" or "set_var" called on them.
      getUnresolvedParams(names);
      break;

    case XEXP_SPECIAL: // ERK.  This doesn't yet track external specials dependencies
      if (newExpPtr_->getTimeDependent()) { names.push_back(std::string("TIME")); }
      if (newExpPtr_->getTempDependent()) { names.push_back(std::string("TEMP")); }
      if (newExpPtr_->getVTDependent()) { names.push_back(std::string("VT")); }
      if (newExpPtr_->getFreqDependent()) { names.push_back(std::string("FREQ")); }
      break;

    case XEXP_VARIABLE:
      //names.insert(names.end(),(xyceGroup->getNames()).begin(), (xyceGroup->getNames()).end());
      break;

    case XEXP_FUNCTION:
      getFunctions(names);
      break;

    case XEXP_NODAL_COMPUTATION:
      break;

    case XEXP_COUNT:
      break;

    default:
      break;

  }

#if 0
  for(int ii=0;ii<names.size();++ii)
  {
    std::cout << "N_UTL_Expression::get_names:  names["<<ii<<"] = " << names[ii] << std::endl;
  }
#endif

  //std::cout << "Expression::get_names(int const & type, std::vector<std::string> & names ) const " << std::endl;
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_type
// Purpose       : Finds the type of an input quantity name
// Special Notes :
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

  std::cout << "Expression::get_type ( const std::string & var ) " << std::endl;
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
bool Expression::make_constant (const std::string & var, const double & val)
{
  bool retVal=false; // ERK.  check this.
  retVal = newExpPtr_->make_constant (var,val);
  //std::cout << "Expression::make_constant (const std::string & var, const double & val).   var = " << var << "  val = " << val << std::endl;
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
bool Expression::make_var (std::string const & var)
{
  bool retVal=false; 
  retVal = newExpPtr_->make_var(var);
  //std::cout << "Expression::make_var (std::string const & var)  var = " << var << std::endl;
  return retVal;
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
// Special Notes :
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
// Special Notes :
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
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes :
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
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes : ERK.  I haven't figured this out yet, but need to.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getLeadCurrents   (std::vector<std::string> & leads) const
{
  //params = newExpPtr_->
  std::cout << "Error. Xyce::Util::Expression::getLeadCurrents not yet implemented." <<std::endl;
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes : 
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
// Function      : Expression::get_num
// Purpose       : Returns the number of input quantities of a requested type
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::get_num(int const & type)
{
  int retVal=0; 

  std::vector<std::string> tmpNames;

  switch (type)
  {
    case XEXP_ALL:
      retVal = newExpPtr_->getVoltOpVec().size() + newExpPtr_->getCurrentOpVec().size() + newExpPtr_->getParamOpVec().size();
      break;

    case XEXP_NODE:
      retVal = newExpPtr_->getVoltOpVec().size();
      break;

    case XEXP_INSTANCE:
      retVal = newExpPtr_->getCurrentOpVec().size();
      break;

    case XEXP_LEAD:
      break;

    case XEXP_STRING: 
      getUnresolvedParams(tmpNames);
      retVal = tmpNames.size();
      break;

    case XEXP_SPECIAL:
      break;

    case XEXP_VARIABLE:
      break;

    case XEXP_FUNCTION:
      retVal = newExpPtr_->getFuncOpVec().size();
      break;

    case XEXP_NODAL_COMPUTATION:
      break;

    case XEXP_COUNT:
      break;

    default:
      break;
  }

  {
    std::map<int, std::string>  typeMap;
    typeMap[0] = std::string( "ALL");
    typeMap[1] = std::string( "NODE");
    typeMap[2] = std::string( "INSTANCE");
    typeMap[3] = std::string( "LEAD");
    typeMap[4] = std::string( "STRING");
    typeMap[5] = std::string( "SPECIAL");
    typeMap[6] = std::string( "VARIABLE");
    typeMap[7] = std::string( "FUNCTION");
    typeMap[8] = std::string( "NODAL_COMPUTATION");

    std::cout << "Expression::get_num(int const & type) for " << newExpPtr_->getExpressionString() << " type[" << type << "]="<<typeMap[type] << " num = " << retVal << std::endl;
  }
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
// Function      : Expression::set_sim_time
// Purpose       : Set 'time' special variable in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::set_sim_time(double time)
{
  // ERK.  this can go.
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Expression::set_sim_dt
// Purpose       : Set time step special variable (dt) in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 12/18/2017
//-----------------------------------------------------------------------------
bool Expression::set_sim_dt(double dt)
{
  // ERK.  this can go.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Expression::set_temp
// Purpose       : Set 'temp' special variable in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::set_temp(double const & tempIn)
{
  // ERK.  this can go.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Expression::set_sim_freq
// Purpose       : Set time step special variable (freq) in expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 12/18/2017
//-----------------------------------------------------------------------------
bool Expression::set_sim_freq(double freq)
{
  // ERK.  this can go.
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Expression::set_accepted_time
// Purpose       : Set accepted time for converged soltion
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void Expression::set_accepted_time(double const time)
{ 
  // ERK.  this can go.
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_break_time
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
double Expression::get_break_time()
{
  double retVal=0.0; 
  //newExpPtr_->evaluate
  // ERK. Note, I shouldn't have to process this list of BP at all, if the API was any good.
  // The API should simply request this vector of breakpoints, and I should return it.
  std::vector<Xyce::Util::BreakPoint> breakPointTimes;
  newExpPtr_->getBreakPoints ( breakPointTimes );
#if 0
  Xyce::Util::BreakPointLess breakPointLess_ = Xyce::Util::BreakPoint::defaultTolerance_;
  std::sort ( breakPointTimes.begin(), breakPointTimes.end(), breakPointLess_ );
  std::vector<Xyce::Util::BreakPoint>::iterator it = std::unique ( breakPointTimes.begin(), breakPointTimes.end());
  breakPointTimes.resize( std::distance (breakPointTimes.begin(), it ));
#endif

  // ERK. This is a total kludge.  This really should just be 
  // replaced by a "getBreakPoints(std::vector<breakpoint> & bpVec)" 
  // call that follows the same patterns as all the other getBreakPoints calls 
  // throughout Xyce (especially in device package).
  //
  // Having logic here to pull out a single BP is silly.  There is better logic 
  // for that sort of thing in the time integrator.
  //
  // Part of the reason for this (bad) structure is that the old expression library
  // doesn't setup breakpoints in a smart way.  There are no sources in the old
  // library that have precomputed formula for breakpoints. (unlike the spice 
  // sources in the device package).  Instead, it "solves" for when the next 
  // breakpoint should be, using a Newton-ish loop.  It does this for *all* 
  // time-dependent expressions, no matter what the nature of their 
  // time dependence.  It does have the benefit of identifying discontinuities 
  // in functions like "STP" (the step function) which do not have set breakpoint times.
  // But it is silly to apply it to things like the PWL or PULSE source, which are sources
  // with known, fixed breakpoints.
  //
  // I should probably attempt to apply the Newton-ish loop to functions like STP, 
  // however.  I had not considered that.
  //
  // Another issue; for time-dependent expressions that do *not* have any 
  // solution variable dependence, Bsrc's have some hidden (weird) behavior.
  // In the Bsrc, if the numExtVars==0, then the evaluate function is not called
  // on the primary expression.  At all.    But, somehow, mysteriously, it gets 
  // updated.  
  //
  // Follow up: I think I just figured this out.  The Bsrc has 2 expression-dependent parameters: V and I.
  // If the core expression does NOT depend on solution variables, then the expression value is simply 
  // set to V or I (depending our src type).   V and/or I are updated during the more 
  // global "updateDependentParams" call, using an "evaluateFunction" call.  If there are no solution vars,
  // then this is sufficient.  However, if there *are* solution vars, then derivatives are needed,
  // so, then device takes over and calls "evaluate" prior to and/or during the load functions.
  //
  double simTime = newExpPtr_->getTime();
  int size = breakPointTimes.size();
  double min = 1.0e+99;
  for (int ii=0;ii<size;++ii)
  {
    double bpTime = breakPointTimes[ii].value();
    double delta = bpTime-simTime;
    if (delta > 0.0 && delta < min)
    {
      min = delta;
      retVal = bpTime;
    }
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_break_time_i
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
double Expression::get_break_time_i()
{
  double retVal=0.0; 
  // ERK.  I don't understand this yet.  Do we need it?
  return retVal;
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
// Function      : Expression::order_names
// Purpose       : Put input quantity names in a particular order (used for
//                 replace_func which requires identical ordering for expression
//                 and user defined function
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::order_names(std::vector<std::string> const & new_names)
{
  // ERK. this can go.
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Expression::replace_func
//
// Purpose       : Replace user defined function with its definition in expression
// 
// Special Notes : ERK.  For the "new" expression library, this function will 
//                 eventually be obsolete.  For now, it is being used to 
//                 integrate the new library in while using the old API.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::replace_func (std::string const & func_name,
                                   Expression & func_def,
                                   int numArgs)
{

  // ERK. This can go.
  int retVal=0; 
#if 0
  Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
  if (!(func_def.newExpPtr_->parsed()))
  {
    func_def.newExpPtr_->lexAndParseExpression();
  }
  xyceGroup->addFunction(func_name, func_def.newExpPtr_);
#endif
  return numArgs;
}

//-----------------------------------------------------------------------------
// Function      : Expression::replace_var
// Purpose       : Replace a variable usage with a parsed sub-expression
//
// Special Notes : This is used for subcircuit parameters that cannot be
//                 fully resolved to a constant because they have global
//                 parameter usage.
//
//                 ERK: See bug 1801.
//                 http://charleston.sandia.gov/bugzilla/show_bug.cgi?id=1801
//
//                 And also tests:
//                 Xyce_Regression/Netlists/Certification_Tests/BUG_1801
//
//                 ERK.  Do I go thru the "group" for this one?  or not bother.
//                 Note that parameters,etc that need this step will probably 
//                 fail to be resolved in my early "resolveExpression" functions.
//
//                 ERK. 4/21/2020.  Addendum.  This is about more than what 
//                 I describe above.  This is also for dealing with parameters 
//                 that are not simply numerical values, but are expressions. 
//                 If the RHS of a global_param is an expression, it isn't 
//                 really correct to call "make_constant" on it,and it also 
//                 isn't correct to call "make_var".  So, instead, this function
//                 is called in this case.
//
//                 This is a perfect candidate for using the "attachParameterNode" 
//                 function.
//
//
// Scope         :
// Creator       : Thomas Russo, SNL
// Creation Date : 08/10/2010
//-----------------------------------------------------------------------------
int Expression::replace_var(
  const std::string &   var_name,
  const Expression &    subexpr)
{
  int retVal=0; 
  //std::cout << "NOTE:  replace_var (expr version) just got called on " << var_name <<std::endl;
  attachParameterNode (var_name, subexpr);
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::replace_var
// Purpose       : Replace a variable usage with a parsed sub-expression
//
// Special Notes : This is used for subcircuit parameters that cannot be
//                 fully resolved to a constant because they have global
//                 parameter usage.
// Scope         :
// Creator       : Thomas Russo, SNL
// Creation Date : 08/10/2010
//-----------------------------------------------------------------------------
int Expression::replace_var (std::string const & var_name,
                             Op::Operator *op )
{
  int retVal=0; 
  {
  std::cout << "NOTE:  replace_var (op version) just got called on " << var_name <<std::endl;
  std::cout << "replace_var (op version) is not implemented yet for newExpression" <<std::endl;
  exit(0);
  }
  return retVal;
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
// Function      : Expression::getNumDdt
// Purpose       : Return the number of ddt() calls in the expression
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::getNumDdt ()
{
  // ERK. Not done.  Might not need it.
  int retVal=0; 
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDdtVals
// Purpose       : Return the most recent arguments of ddt() in the expression
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void Expression::getDdtVals ( std::vector<double> & vals )
{ 
  // ERK. Not done.  Might not need it.
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::setDdtDerivs
// Purpose       : Set the evaluated value of the ddt functions
// Special Notes : This is normally done with derivative values from the
//                 time integration package
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void Expression::setDdtDerivs ( std::vector<double> & vals )
{ 
  // ERK. Not done.  Might not need it.
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::num_vars
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::num_vars() const
{
  // ERK. Not done.  Might not need it.
  int retVal=0; 
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


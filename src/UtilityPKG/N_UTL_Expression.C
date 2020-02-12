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

#include <sstream> 
// ----------   Xyce Includes   ----------
#include <N_UTL_Expression.h>

#include <newExpression.h>
#include <xyceExpressionGroup.h>
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
Expression::Expression( const std::string & exp, bool useNew )
  :
   useNewExpressionLibrary_(useNew),
   newExpPtr_(NULL),
   expPtr_(NULL)
{

  if(useNewExpressionLibrary_)
  {
    if (exp!=std::string(""))
    {
      Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp(new xyceExpressionGroup() );
      grp_ = xyceGroup;
      newExpPtr_ = new newExpression(exp, grp_);
      newExpPtr_->lexAndParseExpression();
    }
  }
  else
  {
    expPtr_ = new ExpressionInternals(exp);
  }
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
   newExpPtr_(NULL),
   expPtr_(NULL)
{
  if(useNewExpressionLibrary_)
  {
    newExpPtr_ = new newExpression( *(right.newExpPtr_));
    newExpPtr_->lexAndParseExpression();
  }
  else
  {
    expPtr_ = new ExpressionInternals( *(right.expPtr_));
  }
  return;
}

#ifdef NEW_EXPRESSION
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
  expPtr_ = new ExpressionInternals( *(right.expPtr_));
  expPtr_->lexAndParseExpression();
  return *this;
}
#endif

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
  if(useNewExpressionLibrary_)
  {
    delete newExpPtr_;
  }
  else
  {
    delete expPtr_;
  }
  return;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
bool Expression::parsed() const 
{
  if(useNewExpressionLibrary_)
    return newExpPtr_->parsed();
  else
    return expPtr_->parsed();
}


//-----------------------------------------------------------------------------
// Function      : Expression::set
// Purpose       : Set the value of the expression to a string
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::set ( const std::string & exp )
{
  bool retVal = false; 

  if(useNewExpressionLibrary_)
  {
    retVal = newExpPtr_->set (exp);
  }
  else
  {
    retVal = expPtr_->set (exp);
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::getSymbolTable
// Purpose       : Returns the symbol table
// Special Notes :
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 08/19/2016
//-----------------------------------------------------------------------------
void Expression::getSymbolTable(std::vector< ExpressionSymbolTableEntry > & theSymbolTable ) const
{ 
  {
#ifdef NEW_EXPRESSION
#else
    expPtr_->getSymbolTable(theSymbolTable);
#endif
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_names
// Purpose       : Returns the names of input quantities by type
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
void Expression::get_names(int const & type, std::vector<std::string> & names ) const
{ 
  {
#ifdef NEW_EXPRESSION
#else
    expPtr_->get_names(type,names);
#endif
  }
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
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->get_type (var);
#endif
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
bool Expression::make_constant (const std::string & var, const double & val)
{
  bool retVal=false; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->make_constant (var,val);
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::make_var
// Purpose       : Convert a 'string' placeholder into a variable
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::make_var (std::string const & var)
{
  bool retVal=false; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->make_var(var);
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::differentiate
// Purpose       : Form the analytic derivative trees for all variables
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::differentiate ()
{
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->differentiate ();
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::set_var
// Purpose       : Sets the value of an input quantity
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::set_var ( const std::string & var, const double & val)
{
#ifdef NEW_EXPRESSION
  Teuchos::RCP<xyceExpressionGroup> testGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
  testGroup->setSoln( var, val );
  return true;
#else
  bool retVal=false; 
  {
    retVal = expPtr_->set_var (var, val);
  }
  return retVal;
#endif
}

//-----------------------------------------------------------------------------
// Function      : Expression::set_vars
// Purpose       : Sets the values of all input quantities
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::set_vars ( const std::vector<double> & vals )
{
  bool retVal=false; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->set_vars ( vals );
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_expression
// Purpose       : Returns a string of the expression
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
std::string Expression::get_expression () const
{
  std::string retVal; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->get_expression ();
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_derivative
// Purpose       : Returns a string of a derivative
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
std::string Expression::get_derivative ( std::string const & var )
{
  std::string retVal; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->get_derivative ( var );
#endif
  }
  return retVal;
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
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->get_num(type);
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluate
// Purpose       : Evaluate expression and derivatives using provided input values
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::evaluate ( double & exp_r,
                                 std::vector<double> & deriv_r,
                                 std::vector<double> & vals )
{
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->evaluate ( exp_r, deriv_r, vals );
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluateFunction
// Purpose       : Evaluate expression using provided input values.  
// Special Notes : This is for cases in which the user does not need 
//                 the derivatives.
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 04/14/08
//-----------------------------------------------------------------------------
int Expression::evaluateFunction ( double & exp_r, std::vector<double> & vals )
{
  int retVal=0;
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->evaluateFunction ( exp_r, vals );
#endif
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
int Expression::evaluate ( double & exp_r,
                                 std::vector<double> & deriv_r)
{
  int retVal=0;
  {
    retVal = expPtr_->evaluate ( exp_r, deriv_r);
  }
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
  {
    retVal = expPtr_->evaluateFunction ( exp_r );
  }
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
#ifdef NEW_EXPRESSION
  Teuchos::RCP<xyceExpressionGroup> testGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
  testGroup->setTime(time);
  return true;
#else
  return expPtr_->set_sim_time(time);
#endif
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
#ifdef NEW_EXPRESSION
  return false;
#else
  return expPtr_->set_sim_dt(dt);
#endif
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
  bool retVal=false;
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->set_temp(tempIn);
#endif
  }
  return retVal;
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
#ifdef NEW_EXPRESSION
  return false;
#else
  return expPtr_->set_sim_freq(freq);
#endif
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
  {
#ifdef NEW_EXPRESSION
#else
    expPtr_->set_accepted_time(time);
#endif
  }
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
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->get_break_time();
#endif
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
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->get_break_time_i();
#endif
  }
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
const std::string & Expression::get_input ()
{
#ifdef NEW_EXPRESSION
  return expPtr_->getExpressionString();
#else
  return expPtr_->get_input ();
#endif
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
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->order_names(new_names);
#endif
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::replace_func
//
// Purpose       : Replace user defined function with its definition in expression
// 
// Special Notes : ERK.  I think for the "new" expression library, this function 
//                 should (probably) be obsolete, unless there is a use case where
//                 it is absolutely necessary.  Currently I don't see a use case
//                 for it.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
int Expression::replace_func (std::string const & func_name,
                                   Expression & func_def,
                                   int numArgs)
{
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->replace_func (func_name, *(func_def.expPtr_), numArgs);
#endif
  }
  return retVal;
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
// Scope         :
// Creator       : Thomas Russo, SNL
// Creation Date : 08/10/2010
//-----------------------------------------------------------------------------
int Expression::replace_var(
  const std::string &   var_name,
  const Expression &    subexpr)
{
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->replace_var (var_name, *(subexpr.expPtr_));
#endif
  }
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
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->replace_var (var_name, op);
#endif
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
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->replace_name ( old_name, new_name);
#endif
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
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_-> getNumDdt ();
#endif
  }
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
  {
#ifdef NEW_EXPRESSION
#else
    expPtr_-> getDdtVals ( vals );
#endif
  }
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
  {
#ifdef NEW_EXPRESSION
#else
    expPtr_->setDdtDerivs ( vals );
#endif
  }
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
  int retVal=0; 
  {
#ifdef NEW_EXPRESSION
#else
    retVal = expPtr_->num_vars();
#endif
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
#ifdef NEW_EXPRESSION
  return false;
#else
  bool implicitTimeDep = expPtr_->isImplicitTimeDepedent();
  bool explicitTimeDep = false;
  std::vector<std::string> specials;
  expPtr_->get_names(XEXP_SPECIAL, specials);
  if (!specials.empty())
  {
    explicitTimeDep=(std::find(specials.begin(), specials.end(), "TIME") != specials.end());
  }
  return (implicitTimeDep || explicitTimeDep);
#endif
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
#ifdef NEW_EXPRESSION
#else
  return ( expPtr_->isRandomDepedent() );
#endif
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
#ifdef NEW_EXPRESSION
  expPtr_->dumpParseTree(std::cout);
#else
  expPtr_->dumpParseTree();
#endif
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
#ifdef NEW_EXPRESSION
#else
  ExpressionInternals::seedRandom(seed);
#endif
}

} // namespace Util
} // namespace Xyce


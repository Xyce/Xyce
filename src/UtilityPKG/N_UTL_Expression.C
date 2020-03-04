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
   namesSet_(false),
   newExpPtr_(NULL),
   expPtr_(NULL)
{

  if(useNewExpressionLibrary_)
  {
    if (exp!=std::string(""))
    {
      Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp(new xyceExpressionGroup() );
      grp_ = xyceGroup;

      // ERK; removing the beginning and ending brace should really be handled by flex/bison, 
      // but I was in a hurry today.
      std::string expCopy = exp;
      if (expCopy[0]== '{' && expCopy[expCopy.size()-1]=='}')
      {
        expCopy.erase(0,1);// lop off open curly brace
        expCopy.erase(expCopy.length()-1); // lop off close curly brace
      }

      newExpPtr_ = new newExpression(expCopy, grp_);
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
   useNewExpressionLibrary_(right.useNewExpressionLibrary_),
   namesSet_(right.namesSet_),
   newExpPtr_(NULL),
   grp_(right.grp_),
   expPtr_(NULL)
{
  if(useNewExpressionLibrary_)
  {
    newExpPtr_ = new newExpression( *(right.newExpPtr_));
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
  useNewExpressionLibrary_ = right.useNewExpressionLibrary_;
  namesSet_ = right.namesSet_;
  grp_ = right.grp_;

  if(useNewExpressionLibrary_)
  {
    newExpPtr_ = new newExpression( *(right.newExpPtr_) );
    newExpPtr_->lexAndParseExpression();
  }
  else
  {
    expPtr_ = new ExpressionInternals( *(right.expPtr_));
  }

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
  {
    return newExpPtr_->parsed();
  }
  else
  {
    return expPtr_->parsed();
  }
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
  if(useNewExpressionLibrary_)
  {
    //newExpPtr_->getSymbolTable(theSymbolTable);
  }
  else
  {
    expPtr_->getSymbolTable(theSymbolTable);
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_names
//
// Purpose       : This function returns the names of various entities present 
//                 in a parsed expression by type.
//
// Special Notes : These notes pertain to the newExpression implementation.
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
//                   XEXP_STRING,         // 4    for some mysterious reason, this means params and global_params
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
  if(useNewExpressionLibrary_)
  {
    switch (type)
    {
      case XEXP_ALL:
        break;

      case XEXP_NODE:
        for (int ii=0;ii<newExpPtr_->getVoltOpVec().size();ii++)
        {
          int size = newExpPtr_->getVoltOpVec()[ii]->getNodeNames().size();

          for (int jj=0;jj<size;jj++)
          {
            std::string tmpName = newExpPtr_->getVoltOpVec()[ii]->getNodeNames()[jj] ;
            std::vector<std::string>::iterator it = std::find(names.begin(), names.end(), tmpName);
            if (it == names.end())
            {
              names.push_back( tmpName );
            }
          }
        }
        break;

      case XEXP_INSTANCE:
        for (int ii=0;ii<newExpPtr_->getCurrentOpVec().size();ii++)
        {
          std::string tmpName = newExpPtr_->getCurrentOpVec()[ii]->getName();
          std::vector<std::string>::iterator it = std::find(names.begin(), names.end(), tmpName);
          if (it == names.end())
          {
            names.push_back( tmpName );
          }
        }
        break;

      case XEXP_LEAD:
        break;

      case XEXP_STRING: // for some mysterious reason, this means params and global_params
        for (int ii=0;ii<newExpPtr_->getParamOpVec().size();ii++)
        {
          std::string tmpName = newExpPtr_->getParamOpVec()[ii]->getName();
          std::vector<std::string>::iterator it = std::find(names.begin(), names.end(), tmpName);
          if (it == names.end())
          {
            names.push_back( tmpName );
          }
        }
        break;

      case XEXP_SPECIAL:
        break;

      case XEXP_VARIABLE:
        break;

      case XEXP_FUNCTION:
        for (int ii=0;ii<newExpPtr_->getFuncOpVec().size();ii++)
        {
          std::string tmpName = newExpPtr_->getFuncOpVec()[ii]->getName();
          std::vector<std::string>::iterator it = std::find(names.begin(), names.end(), tmpName);
          if (it == names.end())
          {
            names.push_back( tmpName );
          }
        }
        break;

      case XEXP_NODAL_COMPUTATION:
        break;

      case XEXP_COUNT:
        break;

      default:
        break;

    }
  }
  else
  {
    expPtr_->get_names(type,names);
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
  if(useNewExpressionLibrary_)
  {
    std::string tmpName = var;
    Xyce::Util::toUpper(tmpName);

    const std::unordered_map<std::string,int> & voltMap = newExpPtr_->getVoltOpNames ();
    const std::unordered_map<std::string,int> & currMap = newExpPtr_->getCurrentOpNames ();

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
      std::cout << "Error. newExpression::get_type.  Cannot find type for " << var << std::endl;
    }
  }
  else
  {
    retVal = expPtr_->get_type (var);
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
  if(useNewExpressionLibrary_)
  {
    retVal = newExpPtr_->make_constant (var,val);
  }
  else
  {
    retVal = expPtr_->make_constant (var,val);
  }
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
//                 for newExpression as well, for obvious reasons.
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
  if(useNewExpressionLibrary_)
  {
#if 0
    std::string tmpParName = var;
    Xyce::Util::toUpper(tmpParName);
    //std::unordered_map<std::string,Teuchos::RCP<astNode<usedType> > > & paramOpMap = 

   //   newExpPtr_->getParamOpMap () ;

    if (newExpPtr_->getParamOpMap().find(tmpParName) != newExpPtr_->getParamOpMap().end())
    {
      newExpPtr_->getParamOpMap()[tmpParName].setVar();
    }
#endif

    newExpPtr_->setVar(var);
  }
  else
  {
    retVal = expPtr_->make_var(var);
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
  if(useNewExpressionLibrary_)
  {
    retVal = newExpPtr_->differentiate ();
  }
  else
  {
    retVal = expPtr_->differentiate ();
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
    xyceGroup->setSolutionVal( var, val );
    return true;
  }
  else
  {
    bool retVal=false; 
    {
      retVal = expPtr_->set_var (var, val);
    }
    return retVal;
  }
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);

    std::cout << "Expression::set_vars" << std::endl;

    if (!namesSet_) // kludge
    {
      std::vector<std::string> names;
      get_names( XEXP_NODE, names); // for now just nodes. make XEXP_ALL later
      get_names( XEXP_INSTANCE, names); // for now just nodes. make XEXP_ALL later
      xyceGroup->setNames ( names );
      namesSet_ = true;

      for (int ii=0;ii<names.size();++ii) { std::cout << "names["<<ii<<"] = " << names[ii] << std::endl; }
    }

    for (int ii=0;ii<vals.size();++ii) { std::cout << "vals["<<ii<<"] = " << vals[ii] << std::endl; }

    xyceGroup->setVals ( vals );
  }
  else
  {
    retVal = expPtr_->set_vars ( vals );
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
  if(useNewExpressionLibrary_)
  {
    //std::cout << "Expression::get_expression not implemented for newExpression library yet" <<std::endl;
    //exit(0);
    retVal = newExpPtr_->getExpressionString(); // note, for new expression, this is not a reconstruction
  }
  else
  {
    retVal = expPtr_->get_expression ();
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
  if(useNewExpressionLibrary_)
  {
    std::cout << "Expression::get_derivative not implemented for newExpression library yet" <<std::endl;
    exit(0);
  }
  else
  {
    retVal = expPtr_->get_derivative ( var );
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);

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
        retVal = newExpPtr_->getParamOpVec().size();
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
  }
  else
  {
    retVal = expPtr_->get_num(type);
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);

    if (!namesSet_) // kludge
    {
      std::vector<std::string> names;
      get_names( XEXP_NODE, names); // for now just nodes. make XEXP_ALL later
      get_names( XEXP_INSTANCE, names); // for now just nodes. make XEXP_ALL later

      // get the global param names.
      newExpPtr_->getGlobalParamNames ( names );

      xyceGroup->setNames ( names );
      namesSet_ = true;
    }

    xyceGroup->setVals ( vals );
    retVal = newExpPtr_->evaluate( exp_r, deriv_r);
  }
  else
  {
    retVal = expPtr_->evaluate ( exp_r, deriv_r, vals );
  }
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluateFunction
// Purpose       : Evaluate expression using provided input values.  
// Special Notes : This is for cases in which the user does not need 
//                 the derivatives.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/14/08
//-----------------------------------------------------------------------------
int Expression::evaluateFunction ( double & exp_r, std::vector<double> & vals )
{
  int retVal=0;
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);

    if (!namesSet_) // kludge
    {
      std::vector<std::string> names;
      get_names( XEXP_NODE, names); // for now just nodes. make XEXP_ALL later
      get_names( XEXP_INSTANCE, names); // for now just nodes. make XEXP_ALL later
      xyceGroup->setNames ( names );
      namesSet_ = true;
    }

    xyceGroup->setVals ( vals );
    retVal = newExpPtr_->evaluateFunction ( exp_r );
  }
  else
  {
    retVal = expPtr_->evaluateFunction ( exp_r, vals );
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
int Expression::evaluate ( double & exp_r, std::vector<double> & deriv_r)
{
  int retVal=0;
  if(useNewExpressionLibrary_)
  {
    retVal = newExpPtr_->evaluate( exp_r, deriv_r );
  }
  else
  {
    retVal = expPtr_->evaluate( exp_r, deriv_r );
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
  if(useNewExpressionLibrary_)
  {
    retVal = newExpPtr_->evaluateFunction ( exp_r );
  }
  else
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
    xyceGroup->setTime(time);
    return true;
  }
  else
  {
    return expPtr_->set_sim_time(time);
  }
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
    return xyceGroup->setTimeStep(dt);
  }
  else
  {
    return expPtr_->set_sim_dt(dt);
  }
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
    return xyceGroup->setTemp(tempIn);
  }
  else
  {
    retVal = expPtr_->set_temp(tempIn);
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
  if(useNewExpressionLibrary_)
  {
    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
    return xyceGroup->setFreq(freq);
  }
  else
  {
    return expPtr_->set_sim_freq(freq);
  }
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    expPtr_->set_accepted_time(time);
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_->get_break_time();
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_->get_break_time_i();
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
  if(useNewExpressionLibrary_)
  {
    return newExpPtr_->getExpressionString();
  }
  else
  {
    return expPtr_->get_input ();
  }
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
  if(useNewExpressionLibrary_)
  {
    newExpPtr_->setFunctionArgStringVec ( new_names );
  }
  else
  {
    retVal = expPtr_->order_names(new_names);
  }
  return retVal;
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
  int retVal=0; 
  if(useNewExpressionLibrary_)
  {
    Xyce::Util::newExpression funcExpr = *(func_def.newExpPtr_) ; // copy construction
    funcExpr.lexAndParseExpression();

    Teuchos::RCP<xyceExpressionGroup> xyceGroup = Teuchos::rcp_static_cast<xyceExpressionGroup>(grp_);
    xyceGroup->addFunction(func_name, funcExpr);
    newExpPtr_->resolveExpression();

    return numArgs;
  }
  else
  {
    retVal = expPtr_->replace_func (func_name, *(func_def.expPtr_), numArgs);
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_->replace_var (var_name, *(subexpr.expPtr_));
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_->replace_var (var_name, op);
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_->replace_name ( old_name, new_name);
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_-> getNumDdt ();
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    expPtr_-> getDdtVals ( vals );
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    expPtr_->setDdtDerivs ( vals );
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
  if(useNewExpressionLibrary_)
  {
  }
  else
  {
    retVal = expPtr_->num_vars();
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
  if(useNewExpressionLibrary_)
  {
    return false;
  }
  else
  {
    bool implicitTimeDep = expPtr_->isImplicitTimeDepedent();
    bool explicitTimeDep = false;
    std::vector<std::string> specials;
    expPtr_->get_names(XEXP_SPECIAL, specials);
    if (!specials.empty())
    {
      explicitTimeDep=(std::find(specials.begin(), specials.end(), "TIME") != specials.end());
    }
    return (implicitTimeDep || explicitTimeDep);
  }
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
  if(useNewExpressionLibrary_)
  {
    return false;
  }
  else
  {
    return ( expPtr_->isRandomDepedent() );
  }
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
  if(useNewExpressionLibrary_)
  {
    newExpPtr_->dumpParseTree(std::cout);
  }
  else
  {
    expPtr_->dumpParseTree();
  }
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
    ExpressionInternals::seedRandom(seed);
  }
}

} // namespace Util
} // namespace Xyce


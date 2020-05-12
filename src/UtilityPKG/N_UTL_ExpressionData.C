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
// Purpose        : Handle data for one expression object
//
// Special Notes  :
//
// Creator        : Richard Schiek
//
// Creation Date  : 8/24/2009
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_Op.h>
#include <N_IO_OptionBlock.h>
#include <N_LAS_Vector.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Param.h>
#include <N_UTL_ExpressionSymbolTable.h>

#include <mainXyceExpressionGroup.h>
#include <outputsXyceExpressionGroup.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : ExpressionData::ExpressionData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::ExpressionData (
      const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group,
      const std::string &   expression)
  : expression_(0),
    expressionString_(expression),
    state_(NOT_SETUP),
    sensitivitiesPossible_(true)
    //,
    //expressionGroup_(group) // don't want to use this
{

  Teuchos::RCP<mainXyceExpressionGroup> mainGroup = Teuchos::rcp_dynamic_cast<mainXyceExpressionGroup>(group);

  Teuchos::RCP<outputsXyceExpressionGroup> outputsGroup = 
    Teuchos::rcp(new outputsXyceExpressionGroup(
      mainGroup->comm_,
      mainGroup->top_,
      mainGroup->analysisManager_,
      mainGroup->deviceManager_,
      mainGroup->outputManager_
        ) );

  expressionGroup_ = outputsGroup;
}


//-----------------------------------------------------------------------------
// Function      : ExpressionData::~ExpressionData
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
ExpressionData::~ExpressionData ()
{
  delete expression_;

  for (Util::Op::OpList::iterator it = expressionOps_.begin(); it != expressionOps_.end(); ++it)
    delete *it;
}


//-----------------------------------------------------------------------------
// Function      : ExpressionData::parsed
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool ExpressionData::parsed() const
{
  return expression_ ? expression_->parsed() : false;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------
double ExpressionData::evaluate(
  Parallel::Machine     comm,
  double                current_circuit_time,
  double                current_circuit_dt,
  const Linear::Vector *  solnVecPtr,
  const Linear::Vector *  stateVecPtr,
  const Linear::Vector *  stoVecPtr,
  const Linear::Vector *  solnVecImagPtr) const
{
  if (state_ == NOT_SETUP)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Must call setup() prior to evaluate()";
  }
  else if (state_ == PARSE_FAILED)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Expression parse failed";
  }
  else if (state_ == UNRESOLVED_SYMBOL)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Unresolved symbols in expression";
  }

  double value = 0.0;

  //if (solnVecPtr)
  {

#if 0
    // ERK.  This doesn't belong here any more, but keeping for reference.  
    // This kind of things belongs in the group.

    // loop over expressionOps_ to get all the values.
    variableValues_.clear();
    for (Util::Op::OpList::const_iterator it = expressionOps_.begin(); it != expressionOps_.end(); ++it)
    {
      Util::Op::OpData opDataTmp(0, solnVecPtr, solnVecImagPtr, stateVecPtr, stoVecPtr, 0);
      variableValues_.push_back( Util::Op::getValue(comm, *(*it), opDataTmp).real());
    }
#endif

    if (expression_)
    {
      expression_->evaluateFunction(value);
    }
  }

  return value;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/16/2015
//-----------------------------------------------------------------------------
double ExpressionData::evaluate(
  Parallel::Machine             comm,
  double                        current_circuit_time,
  double                        current_circuit_dt,
  const Util::Op::OpData &      op_data) const
{
  if (state_ == NOT_SETUP)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Must call setup() prior to evaluate()";
  }
  else if (state_ == PARSE_FAILED)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Expression parse failed";
  }
  else if (state_ == UNRESOLVED_SYMBOL)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Unresolved symbols in expression";
  }

  double value = 0.0;

#if 0
  // ERK.  This doesn't belong here any more, but keeping for reference.  
  // This kind of things belongs in the group.

  // loop over expressionOps_ to get all the values.
  variableValues_.clear();
  for (Util::Op::OpList::const_iterator it = expressionOps_.begin(); it != expressionOps_.end(); ++it)
  {
    variableValues_.push_back( Util::Op::getValue(comm, *(*it), op_data).real());
  }
#endif

  if (expression_)
  {
    expression_->evaluateFunction(value);
  }

  return value;
}


//-----------------------------------------------------------------------------
// Function      : ExpressionData::evaluate
// Purpose       : Evaluate result and derivatives
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/25/2018
//-----------------------------------------------------------------------------
void ExpressionData::evaluate(
  Parallel::Machine             comm,
  double                        current_circuit_time,
  double                        current_circuit_dt,
  const Util::Op::OpData &      op_data,
  double &result, 
  std::vector< double > &derivs 
  ) const
{
  if (state_ == NOT_SETUP)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Must call setup() prior to evaluate()";
  }
  else if (state_ == PARSE_FAILED)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Expression parse failed";
  }
  else if (state_ == UNRESOLVED_SYMBOL)
  {
    Report::DevelFatal().in("ExpressionData::evaluate") << "Unresolved symbols in expression";
  }

#if 0
  // ERK.  This doesn't belong here any more, but keeping for reference.  
  // This kind of things belongs in the group.

  // loop over expressionOps_ to get all the values.
  variableValues_.clear();
  for (Util::Op::OpList::const_iterator it = expressionOps_.begin(); it != expressionOps_.end(); ++it)
  {
    variableValues_.push_back( Util::Op::getValue(comm, *(*it), op_data).real());
  }
#endif

  if (expression_)
  {
    expression_->evaluate( result, derivs);
  }

  return;
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes : just the declaration, definition at end of file
// Creator       : Tom Russo, SNL
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
namespace {
void convertNodalComputation(std::string &nodalComputation, ParamList &paramList);
} // namespace (unnamed)

//-----------------------------------------------------------------------------
// Function      : ExpressionData::setup
//
// Purpose       : Manages all the basic setup for this class.
//
// Special Notes : Originally, all this stuff was in the constructor,
//                 but it needs to happen after topology is
//                 completely done setting itself up, and the
//                 constructor call was too early for that.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/10/04
//-----------------------------------------------------------------------------
ExpressionData::State
ExpressionData::setup(
  Parallel::Machine                     comm,
  const Util::Op::BuilderManager &      op_builder_manager,
  const Util::ParamMap &                context_function_map,
  const Util::ParamMap &                context_param_map,
  const Util::ParamMap &                context_global_param_map)
{
  int unresolved_symbol_count = 0;

  // allocate expression pointer if we need to
  if( expression_ == NULL )
  {
    expression_ = new Expression(expressionGroup_, expressionString_);
  }

  if (!expression_->parsed())
  {
    state_ = PARSE_FAILED;
    return state_;
  }

  // resolve user-defined functions
  {
    std::vector<std::string> global_function_names;
    expression_->getFuncNames(global_function_names);
 
    std::vector<std::string>::iterator it = global_function_names.begin();
    std::vector<std::string>::iterator end = global_function_names.end();
    //const Util::ParamMap & context_function_map = output_manager.getMainContextFunctionMap();
    //expression_->setFunctionMap(context_function_map);
    //
    for ( ; it != end; ++it)
    {
      Util::ParamMap::const_iterator paramMapIter = context_function_map.find(*it);

      if (paramMapIter == context_function_map.end())
      {
        Report::UserError0() << "Cannot find global function definition for " << *it 
          << " in expression " << expression_->get_expression();
        break;
      }

      const Util::Param &functionParameter = (*paramMapIter).second;

      std::string functionPrototype(functionParameter.tag());
      std::string functionBody(functionParameter.stringValue());

      // The function prototype is defined here as a string whose
      // value is the  function name together with its parenthese
      // enclosed comma separated argument list. To resolve a
      // function, create an expression from the function prototype
      // and get its ordered list of arguments via get_names, then
      // create an expression from the function definition and
      // order its names from that list. Finally, replace the
      // function in the expression to be resolved.
      Util::Expression prototypeExression(expressionGroup_, functionPrototype);
      std::vector<std::string> arguments = prototypeExression.getFunctionArgStringVec ();

      // in the parameter we found, pull out the RHS expression and attach
      if(functionParameter.getType() == Xyce::Util::EXPR)
      {
        Util::Expression & expToBeAttached 
          = const_cast<Util::Expression &> (functionParameter.getValue<Util::Expression>());

        // attach the node
        expression_->attachFunctionNode(*it, expToBeAttached);
      }
      else
      {
        Xyce::dout()  << "functionParameter is not EXPR type!!!" <<std::endl;

        switch (functionParameter.getType()) 
        {
          case Xyce::Util::STR:
            Xyce::dout()  <<"It is STR type: " <<  functionParameter.stringValue();
            break;
          case Xyce::Util::DBLE:
            Xyce::dout()  <<"It is DBLE type: " <<  functionParameter.getImmutableValue<double>();
            break;
          case Xyce::Util::EXPR:
            Xyce::dout()  <<"It is EXPR type: " << functionParameter.getValue<Util::Expression>().get_expression();
            break;
          default:
            Xyce::dout()  <<"It is default type (whatever that is): " << functionParameter.stringValue();
        }
      }
    }
  }

  // this varNames vec is a list of string representations of all of
  // the vars in the expression.  
  ParamList param_list;

  // query the expression object for all of its dependent vars.
  std::vector<Xyce::Util::ExpressionSymbolTableEntry> theSymbolTable;
  
  expression_->getSymbolTable(theSymbolTable);
  std::vector<Xyce::Util::ExpressionSymbolTableEntry>::iterator it = theSymbolTable.begin();
  std::vector<Xyce::Util::ExpressionSymbolTableEntry>::iterator end = theSymbolTable.end(); 

  for ( ; it != end; ++it)
  {
    std::string varName = it->name;
    int varType = it->type;
    char varLeadDesignator = it->leadDesignator;
    
    // based on the type of variable, create the needed Param
    // objects for setParmContextType_ to work.

    switch (varType)
    {
      case XEXP_NODAL_COMPUTATION: // this covers things like P(Res)
        // deconstruct the string and turn it into params, push back into
        // param_list
        sensitivitiesPossible_=false; 
        convertNodalComputation(varName, param_list);
        break;

      case XEXP_NODE: // traditional voltage nodes  Even drops (like V(1,2)) handled here
        param_list.push_back(Param( "V" , 1 ));
        param_list.push_back(Param( varName , 0.0 ));
        break;

      case XEXP_INSTANCE:
      {
        std::string currentName("I");
        if( (varLeadDesignator != 0) && (varLeadDesignator!=' ') )
        {
          currentName = currentName + varLeadDesignator;
        }
        param_list.push_back( Param( currentName , 1 ) );
        param_list.push_back( Param( varName , 0.0 ) );
      } 
      break;

      case XEXP_LEAD:
      {
        sensitivitiesPossible_=false;// only difference with XEXP_INSTANCE

        std::string currentName("I");
        if( (varLeadDesignator != 0) && (varLeadDesignator!=' ') )
        {
          currentName = currentName + varLeadDesignator;
        }
        param_list.push_back( Param( currentName , 1 ) );
        param_list.push_back( Param( varName , 0.0 ) );
      }
      break;

      case XEXP_SPECIAL:
        sensitivitiesPossible_=false;
        param_list.push_back( Param( varName , 0.0 ) );
        break;

      case XEXP_STRING:
      {
        Util::ParamMap::const_iterator param_it = context_param_map.find(varName);

        if (param_it != context_param_map.end())
        {
          const Util::Param &replacement_param = (*param_it).second;

          if ( replacement_param.getType() == Xyce::Util::STR ||
               replacement_param.getType() == Xyce::Util::DBLE )
          {
            bool isDotParam=true;
            if (!expression_->make_constant(varName, replacement_param.getImmutableValue<double>(),isDotParam)  )
            {
              Report::UserWarning0() << "Problem converting parameter " << varName << " to its value.";
            }
          }
          else if (replacement_param.getType() == Xyce::Util::EXPR)
          {
            // ERK.  Need to add error messaging to the "attachParameterNode" function to (if need be) output this:
            //  Report::UserWarning0() << "Problem inserting expression " << replacement_param.getValue<Util::Expression>().get_expression()
            //                       << " as substitute for " << varName << " in expression " << expressionString;
            //
            bool isDotParam=true;
            expression_->attachParameterNode (varName, replacement_param.getValue<Util::Expression>(), isDotParam);
          }
        }
        else
        {
          // this block of code will check if the current string is in the global parameter map.
          // If it is, then it will call the "make_var" function for this string.
          // In the old expression library, this marks the string as being something that the 
          // calling code will need to set.  It does not do anything else.
          // Later, the string will be added to the globalParams container, and also the 
          // expVarNames vector.
          param_it = context_global_param_map.find(varName);
          if (param_it != context_global_param_map.end())
          {
            //globalParams.push_back(varName);
            param_list.push_back( Param( "GLOBAL_PARAMETER" , varName ) ); // ERK check this

            if(param_it->second.getType() == Xyce::Util::EXPR)
            {
              Util::Expression & expToBeAttached = const_cast<Util::Expression &> (param_it->second.getValue<Util::Expression>());
              expression_->attachParameterNode(varName, expToBeAttached);
            }
            else
            {
              if (!expression_->make_var(varName))
              {
                Report::UserWarning0() << "Problem setting global parameter " << varName;
              }
            }
          }
          else
          {
            Report::UserWarning0() << "This field: " << varName 
              << " from the expression " << expression_->get_expression() << " is not resolvable";
          }
        }
      }
      break;

      case XEXP_VARIABLE:
      {
        sensitivitiesPossible_=false;
        // this case is a global param that must be resolved at each use because
        // it can change during a simulation
        param_list.push_back( Param( "GLOBAL_PARAMETER" , varName ) );
      }
      break;

      default:
      {
        Report::UserError0() << "Can't find context for expression variable " << varName << " in expression "
                               << expressionString_ << std::endl
                               << "Please check to ensure this parameter is correct and set in your netlist.";
        ++unresolved_symbol_count;
      }
    }
  }

  Parallel::AllReduce(comm, MPI_MIN, &unresolved_symbol_count, 1);

  if (unresolved_symbol_count > 0 )
  {
    state_ = UNRESOLVED_SYMBOL;
    return state_;
  }

  state_ = READY;

#if 0
  {
  ParamList::iterator first = param_list.begin();
  ParamList::iterator theEnd = param_list.end();
  ParamList::iterator iter;
  int i=0;
  for (iter=first;iter!=theEnd;iter++,i++)
  {
    std::cout << "param_list " << i << " " << *iter;
  }
  }
#endif

  Util::Op::makeOps(comm, op_builder_manager, NetlistLocation(), param_list.begin(), param_list.end(), std::back_inserter(expressionOps_));

#if 0
  {
  // recall that: typedef std::vector<Operator *> OpList;  So NOT a LIST!!!  arg
  Util::Op::OpList::iterator first = expressionOps_.begin();
  Util::Op::OpList::iterator theEnd = expressionOps_.end();
  Util::Op::OpList::iterator iter;
  int i=0;
  for (iter=first;iter!=theEnd;iter++,i++)
  {
    Op::Operator * opPtr = (*iter);
    Op::Operator & opRef = *opPtr;
    const std::vector<std::string> & args = opRef.getArgs();
    std::cout << "expressionOps_ " << i << "  " << opRef.getName();
    for (int j=0;j<args.size();j++)
    {
      std::cout << "   " << args[j];
    }
    std::cout << std::endl;
  }
  }
#endif

  return state_;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void ExpressionData::getExpressionArgs(std::vector<std::string> & args) 
{
  args.clear();
  Util::Op::OpList::iterator first = expressionOps_.begin();
  Util::Op::OpList::iterator theEnd = expressionOps_.end();
  Util::Op::OpList::iterator iter;
  for (iter=first;iter!=theEnd;iter++)
  {
    Op::Operator * opPtr = (*iter);
    Op::Operator & opRef = *opPtr;
    const std::vector<std::string> & argsLocal = opRef.getArgs();
    args.insert(args.end(), argsLocal.begin(), argsLocal.end());
  }
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : Tom Russo, SNL
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : convertNodalComputation
// Purpose       : given a nodal expression string (e.g. "VM(A,B)"),
//                 construct the set of Params that makeOps would expect for
//                 it
// Special Notes :
//
// Scope         : file-local
// Creator       : Tom Russo
// Creation Date : 11/27/2013
//-----------------------------------------------------------------------------
void convertNodalComputation(
  std::string &         nodalComputation,
  ParamList &           paramList)
{
  ParamList tempParamList;

  std::size_t firstParen = nodalComputation.find_first_of("(");
  std::size_t lastParen = nodalComputation.find_first_of("(");
  // the length of the name of the param is actually equal to the position
  // of the first paren
  std::string compName=nodalComputation.substr(0,firstParen);
  std::string args=nodalComputation.substr(firstParen+1,nodalComputation.length()-compName.length()-2);

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Processing nodalComputation : " << nodalComputation << Util::push<< std::endl
                 << "name of computation: " << compName << std::endl
                 << "args: " << args << Util::push << std::endl;

  std::size_t firstComma=args.find_first_of(",");
  while (firstComma != std::string::npos)
  {
    std::string arg = args.substr(0,firstComma);
    std::size_t argsLength = args.length();
    args = args.substr(firstComma+1,argsLength-arg.length()-1);
    firstComma = args.find_first_of(",");
    tempParamList.push_back(Param(arg,0.0));
    if (DEBUG_EXPRESSION)
      Xyce::dout() << "arg " << arg << std::endl;
  }

  tempParamList.push_back(Param(args, 0.0));

  if (DEBUG_EXPRESSION)
    Xyce::dout() << "Remaining arg " << args << std::endl
                 << "There were " << tempParamList.size() << " args." << Util::pop << std::endl
                 << Util::pop << std::endl;

  paramList.push_back(Param(compName, static_cast<int>(tempParamList.size())));
  std::copy (tempParamList.begin(), tempParamList.end(), std::back_inserter(paramList));
}

} // namespace (unnammed)
} // namespace Util
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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
// Creation Date : 02/16/2015
//-----------------------------------------------------------------------------
void ExpressionData::evaluate(
  Parallel::Machine             comm,
  double                        current_circuit_time,
  double                        current_circuit_dt,
  const Util::Op::OpData &      op_data,
  double                        &result
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

  if (expression_)
  {
    Teuchos::RCP<outputsXyceExpressionGroup> outputsGroup = Teuchos::rcp_dynamic_cast<outputsXyceExpressionGroup>(expressionGroup_);
    outputsGroup->setOpData(op_data);
    expression_->processSuccessfulTimeStep();
    expression_->evaluateFunction(result);
    expression_->clearOldResult();
  }

  return;
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

  if (expression_)
  {
    Teuchos::RCP<outputsXyceExpressionGroup> outputsGroup = Teuchos::rcp_dynamic_cast<outputsXyceExpressionGroup>(expressionGroup_);
    outputsGroup->setOpData(op_data);

    expression_->processSuccessfulTimeStep();
    expression_->evaluate( result, derivs);
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : ExpressionData::evaluate
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/16/2015
//-----------------------------------------------------------------------------
void ExpressionData::evaluate(
  Parallel::Machine             comm,
  double                        current_circuit_time,
  double                        current_circuit_dt,
  const Util::Op::OpData &      op_data,
  std::complex<double>          &result
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

  if (expression_)
  {
    Teuchos::RCP<outputsXyceExpressionGroup> outputsGroup = Teuchos::rcp_dynamic_cast<outputsXyceExpressionGroup>(expressionGroup_);
    outputsGroup->setOpData(op_data);
    expression_->processSuccessfulTimeStep();
    expression_->evaluateFunction(result);
    expression_->clearOldResult();
  }

  return;
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
  std::complex<double> &result, 
  std::vector< std::complex<double> > &derivs 
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

  if (expression_)
  {
    Teuchos::RCP<outputsXyceExpressionGroup> outputsGroup = Teuchos::rcp_dynamic_cast<outputsXyceExpressionGroup>(expressionGroup_);
    outputsGroup->setOpData(op_data);

    expression_->processSuccessfulTimeStep();
    expression_->evaluate( result, derivs);
  }

  return;
}

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
        Xyce::dout()  << "ExpressionData::setup.  functionParameter is not EXPR type!!!" <<std::endl;

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

  // resolve .param and .global_params
  const std::vector<std::string> params = expression_->getUnresolvedParams(); // params has to be a copy, not a reference!
  for (int ii=0;ii<params.size();ii++)
  {
    std::string varName = params[ii];
    Util::ParamMap::const_iterator param_it = context_param_map.find(varName);

    if (param_it != context_param_map.end())
    {
      const Util::Param &replacement_param = (*param_it).second;

      if ( replacement_param.getType() == Xyce::Util::STR ||
           replacement_param.getType() == Xyce::Util::DBLE )
      {
        enumParamType paramType=DOT_PARAM;
        if (!expression_->make_constant(varName, replacement_param.getImmutableValue<double>(),paramType)  )
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
        enumParamType paramType=DOT_PARAM;
        expression_->attachParameterNode (varName, replacement_param.getValue<Util::Expression>(), paramType);
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
        const Util::Param &replacement_param = param_it->second;

        if(replacement_param.getType() == Xyce::Util::EXPR)
        {
          Util::Expression & expToBeAttached = const_cast<Util::Expression &> (replacement_param.getValue<Util::Expression>());
          expression_->attachParameterNode(varName, expToBeAttached);
        }
        else
        {
          // this should never happen now.
        }
      }
      else
      {
        // disabling this warning because it is sometimes wrong.  If I have {R1:R} on the .print line, and R1 
        // does exist in the circuit, then it *is* resolvable.  But, this function will not correctly figure
        // that out, because in this function it only checks if it is a .param or a .global_param.   The 
        // expression library doesn't know enough yet to exclude R1:R from the "unresolved params" vector.
        // But it will print just fine, once the simulation is up and running, because the outputs group
        // is able to find it correctly.
        //
        // Also, if R1:R isn't legitimately resolvable (i.e. R1 isn't present 
        // in the circuit), then the ouptuts group will complain about it later and 
        // issue a fatal error.  So, this warning isn't necessary, and is sometimes wrong.
        //
        //Report::UserWarning0() << "This field: " << varName 
          //<< " from the expression " << expression_->get_expression() << " is not resolvable";
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

  // Check ops now.
  //
  // ERK.  The old expression library allocated all the necessary ops in this function, and 
  // those allocations provided early error checking to see if the expression included any
  // invalid ops.
  //
  // The new expression library provides values to the expression using the group object.
  // So, all Ops needed by the expression are allocated and evaluated by the group.
  //
  // Currently, the outputs group allocates those Ops as it needs them, so the easiest 
  // way to force them to be allocated at this stage is to simply evaluate the expression.
  // The result of that evaluation is discarded.
  //
  if (expression_)
  {
    double result;
    Teuchos::RCP<outputsXyceExpressionGroup> outputsGroup = Teuchos::rcp_dynamic_cast<outputsXyceExpressionGroup>(expressionGroup_);
    expression_->evaluateFunction(result);
    expression_->clearOldResult();
  }

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

} // namespace Util
} // namespace Xyce

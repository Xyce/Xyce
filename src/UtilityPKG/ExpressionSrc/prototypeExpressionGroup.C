
#include <iostream>
#include <unordered_map>
#include <string>

#include <prototypeExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : prototypeExpressionGroup::prototypeExpressionGroup
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/11/2019
//-------------------------------------------------------------------------------
prototypeExpressionGroup::prototypeExpressionGroup 
    (Xyce::Util::netlistData & nd) 
    : nd_(nd) 
{

  std::vector<double> & solution = nd_.getSolution();
  std::vector<double> & solutionDeriv = nd_.getSolutionDeriv();
  std::vector<double> & solutionInteg = nd_.getSolutionIntegral ();
  std::unordered_map<std::string,Xyce::Util::nodeData> & nodes = nd_.getNodes ();

  solution.resize(nodes.size());
  solutionDeriv.resize(nodes.size());
  solutionInteg.resize(nodes.size());

  if (!(nodes.empty()))
  {
    for ( auto it = nodes.begin(); it != nodes.end(); ++it )
    {
      int solIndex = it->second.index_;
      double value = Xyce::Util::Value(it->second.value_);

      if (solIndex>-1 && solIndex < solution.size())
      {
        solution[solIndex] = value;
        solutionDeriv[solIndex] = value; //ERK just doing this b/c I don't feel like computing actual derivatives in a prototype
        solutionInteg[solIndex] = value; //ERK just doing this b/c I don't feel like computing actual derivatives in a prototype
      }
      else
      {
        std::cout << "Solution out of range!!!" << std::endl;
      }
    }
  }
}

//-------------------------------------------------------------------------------
// Function      : prototypeExpressionGroup::resolveExpression
//
// Purpose       : This function attempts to resolve as many of the 
//                 symbols in the expression as possible.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/4/2019
//-------------------------------------------------------------------------------
bool prototypeExpressionGroup::resolveExpression (Xyce::Util::newExpression & exp)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & paramOpVector = exp.getParamOpVec();

  //---------------------------------------------------------------------------
  // get the unresolved functions
  std::vector<Teuchos::RCP<astNode<usedType> > > & funcOpVector = exp.getFuncOpVec ();
  std::vector<Teuchos::RCP<astNode<usedType> > > & unresolvedFuncOpVector = exp.getUnresolvedFuncOpVec();

  //---------------------------------------------------------------------------
  // Attempt to resolve the unresolved functions (i.e., find them in the functions container, 
  // and then assign the node pointer to the symbol)
  int funcOpSize = funcOpVector.size();
  for (int ii=0;ii<funcOpSize;++ii)
  {
    Xyce::Util::newExpression externalExp;
    if ( getFunction(funcOpVector[ii]->getName(),externalExp) ) // found it
    {
      funcOpVector[ii]->setNode(externalExp.getAst());

      Teuchos::RCP<funcOp<usedType> > tmpPtr = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (funcOpVector[ii]);
      tmpPtr->setFuncArgs(  externalExp.getFunctionArgOpVec() );
    }
    else
    {
      unresolvedFuncOpVector.push_back(funcOpVector[ii]);
    }
  }

  //---------------------------------------------------------------------------
  // Attempt to resolve the unresolved params and global parameters 
  // (i.e., find them in the params/global_params containers, and then 
  // assign the node pointer to the symbol)
  std::vector<Teuchos::RCP<astNode<usedType> > > & unresolvedParamOpVector = exp.getUnresolvedParamOpVector();
  int paramOpVectorSize = paramOpVector.size();
  for (int ii=0;ii<paramOpVectorSize;++ii)
  {
    Xyce::Util::newExpression externalExp;
    if ( getParam(paramOpVector[ii]->getName(),externalExp) ) // found it
    {
      paramOpVector[ii]->setNode(externalExp.getAst());
    }
    else
    {
      if (getGlobalParam(paramOpVector[ii]->getName(),externalExp)) // found it
      {
        paramOpVector[ii]->setNode(externalExp.getAst());
      }
      else
      {
        unresolvedParamOpVector.push_back(paramOpVector[ii]);
      }
    }
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : prototypeExpressionGroup::getFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool prototypeExpressionGroup::getFunction
(const std::string & name, Xyce::Util::newExpression & exp)
{
  std::unordered_map<std::string, Xyce::Util::functionData > & functions = nd_.getFunctions();
  bool retval=true;

  if (functions.find(name) != functions.end()) 
  { 
    exp = *(functions[name].exp_.get()); // this makes a copy of the expression
  }
  else { retval = false; }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : prototypeExpressionGroup::getParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool prototypeExpressionGroup::getParam
(const std::string & name, Xyce::Util::newExpression & exp)
{
  std::unordered_map<std::string, Xyce::Util::paramData > & params = nd_.getParams() ;
  bool retval=true;

  if (params.find(name) != params.end()) 
  { 
    exp = *(params[name].exp_.get()); // copy the underlying expression
  }
  else { retval = false; }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : prototypeExpressionGroup::getGlobalParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool prototypeExpressionGroup::getGlobalParam
(const std::string & name, Xyce::Util::newExpression & exp)
{
  std::unordered_map<std::string, Xyce::Util::globalParamData > & global_params = nd_.getGlobalParams() ;

  bool retval=true;

  if (global_params.find(name) != global_params.end()) 
  { 
    exp = *(global_params[name].exp_.get()); // copy the underlying expression
  }
  else { retval = false; }

  return retval;
}

}
}


#include <iostream>
#include <unordered_map>
#include <string>

#include <xyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::xyceExpressionGroup
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/11/2019
//-------------------------------------------------------------------------------
xyceExpressionGroup::xyceExpressionGroup () :
  time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), dt_(0.0), alpha_(0.0)
{

}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::resolveExpression
//
// Purpose       : This function attempts to resolve as many of the 
//                 symbols in the expression as possible.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/4/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::resolveExpression (Xyce::Util::newExpression & exp)
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


bool xyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  bool success=true;
  retval = 0.0;
  std::string tmp = nodeName;
  Xyce::Util::toUpper(tmp);

  std::vector<std::string>::iterator it = std::find(names_.begin(), names_.end(), tmp);
  if (it != names_.end())
  {
    int index = it - names_.begin();
    retval = dvals_[index];
    std::cout << "Solution variable " << nodeName << " found by the xyceExpresionGroup! value = " << retval << std::endl;
  }
  else // not found
  {
    std::cout << "ERROR.  Solution variable " << nodeName << " not found by the xyceExpresionGroup!" << std::endl;
  }

  return success; // FIX THIS
}

bool xyceExpressionGroup::getGlobalParameterVal (const std::string & nodeName, double & retval )
{
  bool success=true;
  retval = 0.0;
  std::string tmp = nodeName;
  Xyce::Util::toUpper(tmp);

  std::vector<std::string>::iterator it = std::find(names_.begin(), names_.end(), tmp);
  if (it != names_.end())
  {
    int index = it - names_.begin();
    retval = dvals_[index];
  }

  return success; // FIX THIS
}


//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getFunction
(const std::string & name, Xyce::Util::newExpression & exp)
{
  bool retval=true;

  std::string upperName = name;
  Xyce::Util::toUpper(upperName);

  if (functions_.find(upperName) != functions_.end()) { exp = functions_[upperName]; }
  else { retval = false; }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::addFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
void xyceExpressionGroup::addFunction (const std::string & name, Xyce::Util::newExpression & exp)
{
  std::string upperName = name;
  Xyce::Util::toUpper(upperName);
  functions_[upperName] = exp;
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getParam
(const std::string & name, Xyce::Util::newExpression & exp)
{
  bool retval=true;
  return retval; // FIX THIS
}

//-------------------------------------------------------------------------------
// Function      : xyceExpressionGroup::getGlobalParam
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/28/2019
//-------------------------------------------------------------------------------
bool xyceExpressionGroup::getGlobalParam
(const std::string & name, Xyce::Util::newExpression & exp)
{
  bool retval=true;
  return retval; // FIX THIS
}

}
}

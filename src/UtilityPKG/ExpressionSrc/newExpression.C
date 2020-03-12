
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <cmath>

#include "newExpression.h"

#if( defined HAVE__ISNAN_AND__FINITE_SUPPORT )
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#else
#define isnan(x) std::isnan(x)
#define isinf(x) std::isinf(x)
#endif

//-------------------------------------------------------------------------------
// Expression Lexer/Parser header stuff
//
// Grrrr.  Stupid bison 2.4 stopped putting the pre-prologue into the header.
// need this forward declaration
namespace Xyce {
namespace Util {
class ExpressionLexer;
}}

#include "ExpressionParser.hxx"
// BLEAH!   This is here DUPLICATED from ExpressionParser.yxx
// because of STUPID choice in Bison 2.3 to put the post-prologue into the
// .cxx file instead of the .hxx file that Bison 2.1 used to put it in.
#undef yyFlexLexer
  /* CAREFUL watch continuations! */
#define YY_DECL \
int ExpressionLexer::getToken(ExpressionParser::semantic_type *lvalp,  \
                            location *llocp)

  // YECH!  Work around very stupid way that multiple parsers/lexers are
  // handled.
  // Bison's "%name-prefix" is implemented as a #define yylex "prefix"lex
  // which BREAKS flex's C++ lexer: it contains a method "yylex" in the
  // yyFlexLexer class.  Unless we do this kludge, that method gets renamed
  // with the define as well, and the result is a broken set of classes
#undef yylex
#undef yyFlexLexer
#define yyFlexLexer expFlexLexer
#include <FlexLexer.h>
#include "ExpressionLexer.h"
  // if we actually *used* yylex anywhere here it would be necessary to
  // undo that kludge.  Note that because of this stupidity, if the
  // "%name-prefix" is changed, this line needs to be changed, too.
  // BUT we don't actually use yylex anywhere in this file, so let's
  // leave yylex undefined.  If later it turns out that this *becomes*
  // necessary, uncomment the next line.
  //  #define yylex XyceDevicelex
//-------------------------------------------------------------------------------

namespace Xyce {
namespace Util {


//-------------------------------------------------------------------------------
// Function      : newExpression::lexAndParseExpression
// Purpose       : Lexes and Parses the expression string
//
// Special Notes : The lexer and parser objects are local to this function.
//                 If successful, It will have created an Abstract Syntax Tree,
//                 and the top level pointer to that tree will be set.
//
//                 If this function hasn't been called yet, the astNodePtr will be
//                 NULL.  After it has been called, it should contain a valid
//                 pointer value.
//
//                 Also, in addition to allocating the AST, the parsing process
//                 will populate the paramOpMap object.  Prior to calling this function
//                 that map will be empty.  After calling this function that map
//                 may have some stuff in it.
//
//                 This is not called automatically.  After the expression is
//                 allocated, a second call must be made to parse it.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 10/2?/2019
//-------------------------------------------------------------------------------
bool newExpression::lexAndParseExpression()
{
  if (traditionalParse_)
  {
    std::string fileName("bogusTestFile");
    std::stringstream expressionStringStream ( expressionString_ );
    Xyce::Util::ExpressionLexer expLexer(fileName,&expressionStringStream);
    XyceExpression::ExpressionParser expParser(&expLexer,*this);
    int retCode = expParser.parse();
    parsed_ = (retCode == 0);
  }
  else
  {
    parsed_=true;
  }

  // pull the function arguments (if they are present) out of the parameter vector
  {
    // The paramOpMap_ was set up during parsing (see the ExpressionParser.yxx file and the
    // ExpressionParser.cxx file.  All the code for setting up the paramOpMap_ is there)
    //
    // The functionArgStringVec was hopefully set right after the expression was allocated.
    // (if not, code below won't work).  The functionArgStringVec contains the "prototype"
    // ie. for .func f(x,y) {2*x+3*y}, "x" and "y" would be found in the
    // functionArgStringVec object.  Since they are passed into the function, they should
    // not be considered as params or global_params.
    //
    // Note, if using the flex/bison NetlistParser,
    // when each .func statement is processed, a functionData object is created, and different
    // fields are populated, including "args_", which is where functionArgStringVec comes from
    // (it gets copied over in ExpressionTest).  So, anyway, paramOpMap comes from expression
    // parsing, but functionArgStringVec comes from Netlist parsing.
    //
    // the paramOpMap_, immediately after parsing will contain
    // both regular params and function arguments.  It cannot tell the difference yet.  This
    // next bit of code is designed to pull them apart.  When finished, the elements of
    // the param container should contain ONLY parameters and no function arguments.
    //
    // If the functionArgStringVec object is empty, then it is assumed that this expression
    // is not a function, and thus it has no arguments to resolve.
    //
    // Also: this search is performed before any attempt is made to resolve
    // parameters.  So, if the user has used the same symbol for both a function
    // argument and a parameter, the function argument interpretation will take
    // priority, since it is being checked first.  This seems reasonable, as the
    // function argument will be more of a "local" variable.

    int stringArgsSize=functionArgStringVec_.size();
    std::vector<std::string>::iterator paramIter;
    if (stringArgsSize>0)
    {
      functionArgOpVec_.clear();
      functionArgOpVec_.resize(stringArgsSize,(getGarbageParam()));
      for (int ii=0;ii<stringArgsSize;++ii)
      {
        paramIter = std::find(paramNameVec_.begin(),paramNameVec_.end(), functionArgStringVec_[ii]);
        if (paramIter != paramNameVec_.end()) // found it
        {
          int index = std::distance(paramNameVec_.begin(),paramIter);
          functionArgOpVec_[ii] = paramOpVec_[index];
          paramNameVec_.erase(paramNameVec_.begin()+index);
          paramOpVec_.erase(paramOpVec_.begin()+index);
          functionArgOpVec_[ii]->setFunctionArgType();
        }
      }
    }
  }

  {
    for (int ii=0;ii<voltOpVec_.size();++ii)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (voltOpVec_[ii]);
      std::vector<std::string> & tmp = voltOp->getVoltageNodes();
      for (int jj=0;jj<tmp.size();++jj)
      {
        voltOpNames_[tmp[jj]].push_back(voltOpVec_[ii]);
      }
    }

    for (int ii=0;ii<currentOpVec_.size();++ii)
    {
      Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);
      std::string tmp = currOp->getCurrentDevice();
      currentOpNames_[tmp].push_back(currentOpVec_[ii]);
    }
  }

#if 0
  dumpParseTree(std::cout);
#endif

  return parsed_;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::resolveExpression
// Purpose       : resolve parameters, and external functions.
//
// This must be a separate phase from the setup in lexAndParseExpression function, as
// *all* the relevant expressions must be allocated and have gone thru their initial
// set up before the group can perform the next phase.
//
// The role of the "group_" is to provide a lookup for external .func,
// .param and .global_param expressions, that must be used to resolve them in
// this expression.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 11/5/2019
//-------------------------------------------------------------------------------
void newExpression::resolveExpression ()
{
  //---------------------------------------------------------------------------
  // Attempt to resolve the unresolved functions. Get them from the group,
  // which is the newExpression class connection to other expressions
  // and then assign the node pointer to the symbol.
  int funcOpSize = funcOpVec_.size();
  for (int ii=0;ii<funcOpSize;++ii)
  {
    Xyce::Util::newExpression exp;
    if ( group_->getFunction(funcOpVec_[ii]->getName(),exp) ) // found it
    {
      funcOpVec_[ii]->setNode(exp.getAst());

      Teuchos::RCP<funcOp<usedType> > tmpPtr = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (funcOpVec_[ii]);
      if ( !(Teuchos::is_null(tmpPtr)) )
      {
        tmpPtr->setFuncArgs(  exp.getFunctionArgOpVec() );
        externalDependencies_ = true;
      }
    }
    else
    {
      unresolvedFuncOpVec_.push_back(funcOpVec_[ii]);
    }
  }

  //---------------------------------------------------------------------------
  // Attempt to resolve the unresolved params and global parameters, using the
  // same approach.
  int paramOpVec_Size = paramOpVec_.size();
  for (int ii=0;ii<paramOpVec_Size;++ii)
  {
    Xyce::Util::newExpression exp;
    if ( group_->getParam(paramOpVec_[ii]->getName(),exp) ) // found it
    {
      paramOpVec_[ii]->setNode(exp.getAst());
      externalDependencies_ = true;
    }
    else
    {
      if (group_->getGlobalParam(paramOpVec_[ii]->getName(),exp)) // found it
      {
        paramOpVec_[ii]->setNode(exp.getAst());
        externalDependencies_ = true;
      }
      else
      {
        unresolvedParamOpVec_.push_back(paramOpVec_[ii]);
      }
    }
  }
  expressionResolved_ = true;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void newExpression::clear ()
{
  // copied from destructor
  if (astNodePtrPtr_)
  {
    delete astNodePtrPtr_;
  }
  if ( tableNodePtrPtr_ )
  {
    delete tableNodePtrPtr_;
  }

  // vectors of pointers to RCPs of AST nodes; use constructors
  for (int ii=0;ii<masterAstNodeVec_.size();ii++)
  {
    delete masterAstNodeVec_[ii];
  }
  masterAstNodeVec_.clear();

  for (int ii=0;ii<srcAstNodeVec_.size();ii++)
  {
    delete srcAstNodeVec_[ii];
  }
  masterAstNodeVec_.clear();


  expressionString_ = std::string("");
  parsed_ = false;
  derivsSetup_ = false;
  astArraysSetup_ = false;
  astNodePtrPtr_ = NULL;
  tableNodePtrPtr_ = NULL;
  numDerivs_ = 0;
  traditionalParse_ = true;
  externalDependencies_ = false;

  functionArgStringVec_.clear();
  functionArgOpVec_.clear();

  paramNameVec_.clear();
  paramOpVec_.clear();
  unresolvedParamOpVec_.clear();

  funcOpVec_.clear();
  unresolvedFuncOpVec_.clear();

  voltOpVec_.clear();
  unresolvedVoltOpVec_.clear();
  voltOpNames_.clear();

  currentOpVec_.clear();
  unresolvedCurrentOpVec_.clear();
  currentOpNames_.clear();

  bpTol_ = 0.0;
  timeStep_ = 0.0;
  timeStepAlpha_ = 0.0;
  timeStepPrefac_ = 0.0;

  derivIndexVec_.clear();

  return;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
bool newExpression::make_constant (std::string const & var, usedType const & val)
{
  std::string tmpParName = var;
  Xyce::Util::toUpper(tmpParName);
  bool retval=false;

  std::vector<std::string>::iterator paramIter;
  paramIter = std::find(paramNameVec_.begin(),paramNameVec_.end(), tmpParName);
  if (paramIter != paramNameVec_.end()) // found it
  {
    int index = std::distance(paramNameVec_.begin(),paramIter);
    Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[index]);
    parOp->setValue(val);
    retval=true;
  }

  return retval;
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
bool newExpression::make_var (std::string const & var)
{
  std::string tmpParName = var;
  Xyce::Util::toUpper(tmpParName);
  bool retval=false;

  std::vector<std::string>::iterator paramIter;
  paramIter = std::find(paramNameVec_.begin(),paramNameVec_.end(), tmpParName);
  if (paramIter != paramNameVec_.end()) // found it
  {
    int index = std::distance(paramNameVec_.begin(),paramIter);
    Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[index]);
    parOp->unsetValue(); // just to be safe "unset" the value
    parOp->setVar();
    retval = true; // just means we found it
  }

  return retval;
}

//-------------------------------------------------------------------------------
// ERK.  do we really need this function and the one above (make_var)?
//-------------------------------------------------------------------------------
void newExpression::setVar(const std::string & var)
{
  std::string tmpParName = var;
  Xyce::Util::toUpper(tmpParName);
  std::vector<std::string>::iterator paramIter;
  paramIter = std::find(paramNameVec_.begin(),paramNameVec_.end(), tmpParName);
  if (paramIter != paramNameVec_.end()) // found it
  {
    int index = std::distance(paramNameVec_.begin(),paramIter);
    Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[index]);
    parOp->setVar();
  }
}

//-------------------------------------------------------------------------------
// ERK.  This probably isn't needed.  In the old expression library, this
// function`appears to be equivalent (kind of) to the
// newExpression::setupDerivatives_  function.
//-------------------------------------------------------------------------------
int newExpression::differentiate ()
{
  std::cout << "newExpression::differentiate not set up yet." <<std::endl;
  exit(0);
  return -1;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setupDerivatives_
//
// Purpose       : this is yet another phase in setup, which must be called after
//                 the "resolveExpression" calls (above) for all relevant expressions
//                 have been completed first.  If they are not, then there
//                 will be some NULL pointers in the AST tree, which will cause
//                 seg faults.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter
// Creation Date : 11/5/2019
//-------------------------------------------------------------------------------
void newExpression::setupDerivatives_ ()
{
  // figure out the derivative indices
  // This assumes we need derivatives with respect to all
  // voltages and currents.
  //
  // I was considering automatically doing parameters as well, but
  // that isn't what the old library does.  Params are only differentiated
  // if "make_var" is called on them.
  {
    numDerivs_=0;
    derivIndexVec_.clear();

    for (int ii=0;ii<voltOpVec_.size();ii++)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (voltOpVec_[ii]);
      std::vector<std::string> & nodes = voltOp->getVoltageNodes();

      if (nodes.size() == 1) // putting this here b/c I now want to handle the V(A,B) case differently than I planned for
      {
        std::string tmp = nodes[0]; Xyce::Util::toUpper(tmp);
        std::unordered_map<std::string, int>::iterator mapIter;
        mapIter = derivNodeIndexMap_.find(tmp);
        if (mapIter == derivNodeIndexMap_.end()) { derivNodeIndexMap_[tmp] = numDerivs_; numDerivs_++; }
        derivIndexVec_.push_back(derivIndexPair_(voltOpVec_[ii],derivNodeIndexMap_[tmp]));
      }
      else
      {
        std::cout << "ERROR. derivatives not correct for 2-node V(A,B) specification" <<std::endl;
      }
    }

    for (int ii=0;ii<currentOpVec_.size();ii++)
    {
      Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);
      std::string tmp = currOp->getCurrentDevice();
      Xyce::Util::toUpper(tmp);
      std::unordered_map<std::string, int>::iterator mapIter;
      mapIter = derivNodeIndexMap_.find(tmp);
      if (mapIter == derivNodeIndexMap_.end())
      {
        derivNodeIndexMap_[tmp] = numDerivs_; numDerivs_++;
      }
      derivIndexPair_ currentNodePair(currentOpVec_[ii],derivNodeIndexMap_[tmp]);
      derivIndexVec_.push_back(currentNodePair);
    }

    // Complication for params:
    // Unlike voltages and currents, params can be assigned to each other, etc, and they can be *anything*.
    // So, for example if I have
    //
    //  .param a = {10*c+b}
    //  .param b = {5*x}
    //  .param c = {sqrt(y)}
    //  .param x = 10.0
    //  .param y = 20.0
    //
    // the most likely params I would want differentiated are x and y.
    // Differentiating w.r.t. a,b,c might not make sense, as they are
    // placeholders for other expressions.
    //
    // So, possibly I need to do this only w.r.t. terminal parameters,
    // that are set to simple numerical values.  That would mean
    // that for the above set of expressions, I only would do x and y.
    //
    // Solution: only differentiate params that have their "setVar" boolean set to true.
    //

    for (int ii=0;ii<paramOpVec_.size();ii++)
    {
      Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
      if (parOp->getVar())
      {
        std::string tmp = parOp->getName();
        Xyce::Util::toUpper(tmp);
        std::unordered_map<std::string, int>::iterator mapIter;
        mapIter = derivNodeIndexMap_.find(tmp);
        if (mapIter == derivNodeIndexMap_.end())
        {
          derivNodeIndexMap_[tmp] = numDerivs_; numDerivs_++;
        }
        derivIndexPair_ parNodePair(paramOpVec_[ii],derivNodeIndexMap_[tmp]);
        derivIndexVec_.push_back(parNodePair);
      }
    }

    if (group_->isOption (std::string("verbose")))
    {
      std::cout << "---------------------------------------------------------------"<<std::endl;
      std::cout << "For expression " << expressionString_ << " there are ";
      std::cout << numDerivs_ << " derivative index pairs, ";
      std::cout << voltOpVec_.size() << " voltage nodes, ";
      std::cout << currentOpVec_.size() << " current nodes, ";
      std::cout << paramOpVec_.size() << " param nodes ";
      std::cout << functionArgOpVec_.size() << " function argument nodes ";
      std::cout <<std::endl;

      if (numDerivs_>0)
      {
        std::cout << "Derivative indices:"<<std::endl;
        for (int ii=0;ii<numDerivs_;ii++)
        {
          std::cout << ii << ": index = " << derivIndexVec_[ii].second;
          derivIndexVec_[ii].first->output(std::cout);
          //std::cout << std::endl;
        }
      }
      if (!functionArgOpVec_.empty())
      {
        std::cout << "Function Arguments:" <<std::endl;
        for (int ii=0;ii<functionArgOpVec_.size();ii++)
        {
          functionArgOpVec_[ii]->output(std::cout);
        }
        std::cout <<std::endl;
      }
      std::cout << "---------------------------------------------------------------"<<std::endl;
    }
  }

  derivsSetup_ = true;
}

// ERK.  12/26/2019. This may be refactored away later.
// This is not the best way to do this.
// But, I needed it for several parserUnitTests to pass,
// post-RCP refactor.
//
// Notes from 2/7/2020
//
// (1) This function had a bug in it that was not revealed until I made calling
// this function mandatory when "evaluate" or "evaluateFunction" is called.
//
// (2) The flaw was in handling function arguments.  The "getInterestingOps"
// function was using dummy args instead of actual args.  This has been fixed now.
//
// (3) another flaw, which is not fixed, is that the code treats function arguments
// are parameters.  And even after their paramNodes are "set" by the actual args,
// they will still get included in the param array.
//
// (4) In principal, they don't need to be in this array, as they are function args
// NOT regular params such as global params, etc.  So possibly the right thing to
// do is create a different Op for function args.   This is not done yet.
//
// (5) this second issue does not break the code, but adds an inefficiency and
// memory bloat.
//
// private function
void newExpression::setupVariousAstArrays_()
{
#if 0
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Various arrays for expression: " << expressionString_ <<std::endl;
  if (externalDependencies_) std::cout << "externalDependencies_ = true" <<std::endl;
  else std::cout << "externalDependencies_ = false" <<std::endl;

  std::cout << "Array sizes BEFORE update:" <<std::endl;
  std::cout << "Size of paramOpVec_ = " << paramOpVec_.size() << std::endl;
  std::cout << "Size of funcOpVec_ = " << funcOpVec_.size() << std::endl;
  std::cout << "Size of voltOpVec_ = " << voltOpVec_.size() << std::endl;
  std::cout << "Size of currentOpVec_ = " << currentOpVec_.size() << std::endl;

  std::cout << "paramOpVec_:" << std::endl;
  for (int ii=0;ii<paramOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    paramOpVec_[ii]->output(std::cout,0);
  }

  std::cout << "funcOpVec_:" << std::endl;
  for (int ii=0;ii<funcOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    funcOpVec_[ii]->output(std::cout,0);
  }

  std::cout << "voltOpVec_:" << std::endl;
  for (int ii=0;ii<voltOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    voltOpVec_[ii]->output(std::cout,0);
  }

  std::cout << "currentOpVec_:" << std::endl;
  for (int ii=0;ii<currentOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    currentOpVec_[ii]->output(std::cout,0);
  }
#endif

  if (externalDependencies_)
  {
    paramOpVec_.clear();
    funcOpVec_.clear();
    voltOpVec_.clear();
    currentOpVec_.clear();

    if( !(Teuchos::is_null(astNodePtr_)) )
    {
      if (astNodePtr_->paramType())
      {
        if ( !(astNodePtr_->getFunctionArgType()) )  // parameters are occasionally function arguments.  Don't include those
        {
          paramOpVec_.push_back(astNodePtr_);
        }
      }
      if (astNodePtr_->funcType())
      {
        funcOpVec_.push_back(astNodePtr_);
      }
      if (astNodePtr_->voltageType())
      {
        voltOpVec_.push_back(astNodePtr_);
      }
      if (astNodePtr_->currentType())
      {
        currentOpVec_.push_back(astNodePtr_);
      }

      astNodePtr_->getInterestingOps(paramOpVec_,funcOpVec_, voltOpVec_,currentOpVec_);
      // this won't work with RCP
      //std::sort (paramOpVec_.begin(), paramOpVec_.end());
      //std::unique (paramOpVec_.begin(), paramOpVec_.end());
    }

    for (int ii=0;ii<voltOpVec_.size();++ii)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (voltOpVec_[ii]);
      std::vector<std::string> & tmp = voltOp->getVoltageNodes();
      for (int jj=0;jj<tmp.size();++jj)
      {
        voltOpNames_[tmp[jj]].push_back(voltOpVec_[ii]);
      }
    }

    for (int ii=0;ii<currentOpVec_.size();++ii)
    {
      Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);
      std::string tmp = currOp->getCurrentDevice();
      currentOpNames_[tmp].push_back(currentOpVec_[ii]);
    }
  }


  // need to more thoroughly set up the specials booleans.  Check external dependencies for these nodes.



#if 0
  std::cout << "Array sizes AFTER update:" <<std::endl;

  std::cout << "Size of paramOpVec_ = " << paramOpVec_.size() << std::endl;
  std::cout << "Size of funcOpVec_ = " << funcOpVec_.size() << std::endl;
  std::cout << "Size of voltOpVec_ = " << voltOpVec_.size() << std::endl;
  std::cout << "Size of currentOpVec_ = " << currentOpVec_.size() << std::endl;

  std::cout << "paramOpVec_:" << std::endl;
  for (int ii=0;ii<paramOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    paramOpVec_[ii]->output(std::cout,0);
  }

  std::cout << "funcOpVec_:" << std::endl;
  for (int ii=0;ii<funcOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    funcOpVec_[ii]->output(std::cout,0);
  }

  std::cout << "voltOpVec_:" << std::endl;
  for (int ii=0;ii<voltOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    voltOpVec_[ii]->output(std::cout,0);
  }

  std::cout << "currentOpVec_:" << std::endl;
  for (int ii=0;ii<currentOpVec_.size();ii++)
  {
    std::cout << ii << " ";
    currentOpVec_[ii]->output(std::cout,0);
  }
#endif

  astArraysSetup_ = true;
};

//-------------------------------------------------------------------------------
// these two functions return int error codes in the original expression library
//-------------------------------------------------------------------------------
int newExpression::evaluate (usedType &result, std::vector< usedType > &derivs)
{
  int retVal=0;
  if (parsed_)
  {
    if (!expressionResolved_) { resolveExpression(); }

    if (!astArraysSetup_)
    {
      setupVariousAstArrays_ ();
    }

    if (!derivsSetup_)
    {
      setupDerivatives_ ();
    }

#if 0
    dumpParseTree(std::cout);
#endif

    int err1 = evaluateFunction (result);

    if (derivs.size() != numDerivs_) {derivs.clear(); derivs.resize(numDerivs_);}

    for (int ii=0;ii<derivIndexVec_.size();ii++) { derivIndexVec_[ii].first->setDerivIndex(derivIndexVec_[ii].second); }

    for (int ii=0;ii<numDerivs_;++ii) { derivs[ii] = astNodePtr_->dx(ii); }

    for (int ii=0;ii<derivIndexVec_.size();ii++) { derivIndexVec_[ii].first->unsetDerivIndex(); }
  }
  else
  {
    std::cout << "Error.  Expression " << expressionString_ << " is not parsed yet" << std::endl;
  }

  // fix these properly for std::complex later.
  for(int ii=0;ii<derivs.size();++ii)
  {
    if ( isnan(std::real(derivs[ii])) )
    {
      derivs[ii] = 0.0;
    }
    if ( isinf(std::real(derivs[ii])) ) 
    {
      derivs[ii] = 1.0e+10; // fix this
    }
  }
  // old expression library returns EXPRerrno, which is a static variable.
  // If it is zero, everything is cool.
  return retVal;
};

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
int newExpression::evaluateFunction (usedType &result)
{
  int retVal=0;
  if (parsed_)
  {
    if (!expressionResolved_) { resolveExpression(); }

    if (!astArraysSetup_) { setupVariousAstArrays_ (); }

    // get solution values we need from the group
    {
      for (int ii=0;ii<voltOpVec_.size();ii++)
      {
        Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (voltOpVec_[ii]);
        std::vector<std::string> & nodes = voltOp->getVoltageNodes();
        std::vector<usedType> & vals = voltOp->getVoltageVals();

        for (int jj=0;jj<nodes.size();jj++)
        {
          group_->getSolutionVal(nodes[jj], vals[jj]);
        }
      }

      for (int ii=0;ii<currentOpVec_.size();ii++)
      {
        Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);

        usedType val;
        group_->getSolutionVal(currOp->getCurrentDevice(),val);
        currOp->setCurrentVal ( val );
      }

     // ERK. The global parameter setting here should eventually go away.
     // It is only here b/c I am trying to maintain the API to the old expression library for now.
      for (int ii=0;ii<paramOpVec_.size();++ii)
      {
        Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);

        // the "isVar" boolean currently serves two purposes.
        //
        // (1) if we want derivatives w.r.t. it.
        //    (this meaning will probably persist after old API is gone, albeit with a different name)
        //
        // (2) if it should be considered a dynamic variable that gets its values from the vals array.
        //    (this meaning will be obsolete later, once old API is gone)
        //
        // The old API is set up to set *everything* thru a single std::vector of values,
        // that is passed into the expression as a function argument to the "evaluate" function.
        // This vector includes both
        //
        //
        if ( parOp->getVar() )
        {
          usedType val;
          group_->getGlobalParameterVal(parOp->getName(),val);
          parOp->setValue(val);
        }
      }
    }

    timeNodePtr_->setValue(group_->getTime());
    tempNodePtr_->setValue(group_->getTemp()); // assuming correct units.  Conversion happens in the group
    vtNodePtr_->setValue(group_->getVT());
    freqNodePtr_->setValue(group_->getFreq());

    bpTol_ = group_->getBpTol();
    timeStep_ = group_->getTimeStep ();
    timeStepAlpha_ = group_->getTimeStepAlpha ();
    timeStepPrefac_ = group_->getTimeStepPrefac ();

    result = astNodePtr_->val();

    // ERK: fix this failsafe properly for std::complex 
    if (isnan(std::real(result)))
    {
      result = 0.0;
    }
    if (isinf(std::real(result)))
    {
      result = 1.0e+20; // ERK fix this
    }
  }
  else
  {
    std::cout << "Error.  Expression " << expressionString_ << " is not parsed yet" << std::endl;
  }

  return retVal;
};


// some of the parameter and function objects are stored in multiple containers.
void newExpression::setFunctionArgStringVec (const std::vector<std::string> & args)
{
  functionArgStringVec_ = args;
  int size = functionArgStringVec_.size();
  for (int ii=0;ii<size;ii++)
  {
    Xyce::Util::toUpper(functionArgStringVec_[ii]);
  }
};
}
}

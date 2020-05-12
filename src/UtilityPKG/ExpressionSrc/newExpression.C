
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
  
#if 0
  std::cout << "lexAndParseExpression for " << expressionString_ <<std::endl;
#endif

  if (traditionalParse_)
  {
    std::string fileName("bogusTestFile");
    std::stringstream expressionStringStream ( expressionString_ );

    // if the expressionString_ is empty, bison will throw an error.  
    // Plus, no point in parsing an empty string.
    if(expressionString_.empty())
    {
      parsed_ = true;  // ERK.   If this is false, and someone calls "evaluate", triggers a fatal error. So let it be true.
    }
    else
    {
      Xyce::Util::ExpressionLexer expLexer(fileName,&expressionStringStream);
      XyceExpression::ExpressionParser expParser(&expLexer,*this);
      int retCode = expParser.parse();
      parsed_ = (retCode == 0);
    }
  }
  else
  {
    parsed_=true;
  }

  // pull the function arguments (if they are present) out of the parameter vector
  {
    // The following comments are pertainent to the RHS of .func expressions.  For example,
    // if the netlist has the following function: .func f(x,y) {2*x+3*y}, I am calling:
    //
    // LHS ->  f(x,y) 
    // RHS ->  {2*x+3*y}
    //
    // Xyce processes the LHS and RHS of a .func declaration separately, but the 
    // RHS needs to know information from the LHS to do it right.  The following 
    // comments and code below are all about how to correctly handle "x" and "y" 
    // in the RHS expression {2*x+3*y}.  They need to be treated differently 
    // than params.
    //
    // The paramOpVec_ was set up during parsing of the RHS. See the ExpressionParser.yxx 
    // file and the Bison-produced ExpressionParser.cxx file.  All the code for 
    // setting up the paramOpVec_ is there.
    //
    // The functionArgStringVec is not set up during RHS parsing, as just parsing the 
    // RHS of an expression doesn't have enough information to know which things 
    // are function arguments and which things are parameters.  So, functionArgStringVec
    // should be set up and passed into newExpression BEFORE ::lexAndParseExpression 
    // is called, probably right after the expression was allocated.  If not, then
    // the code below will not work.  The functionArgStringVec contains the "prototype"
    // arguments for a function.  It must come from parsing of the LHS of the 
    // .func declaration.  This code (newExpression) doesn't care how this parsing is 
    // done, of course.  But note that the Xyce IO package actually uses the 
    // expression library to determine functionArgStringVec, but it does so by 
    // allocating a completely different expression object of the LHS "f(x,y)" 
    // and parsing it.
    //
    // In the example .func f(x,y) {2*x+3*y}, the functionArgStringVec object would contain
    // "x" and "y".  Since they are passed into the function, they should
    // not be considered as params or global_params, and excluded from any operations
    // that are specific to params/global_params.  As noted, however, at parse time 
    // for the RHS it is impossible to tell the difference.
    //
    // Anyway, the paramOpVec_, immediately after parsing the RHS will contain
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

  // set up names vectors for voltage and current nodes.
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

  // if dependent on a special, add relevant specials node to the relevant specials vector
  if(isTimeDependent_) { timeOpVec_.push_back(timeNodePtr_); }
  if(isTempDependent_) { tempOpVec_.push_back(tempNodePtr_); }
  if(isVTDependent_) { vtOpVec_.push_back(vtNodePtr_); }
  if(isFreqDependent_) { freqOpVec_.push_back(freqNodePtr_); }
  if(isGminDependent_) { gminOpVec_.push_back(gminNodePtr_); }

#if 0
  dumpParseTree(std::cout);
#endif

  return parsed_;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::attachFunctionNode
//
// Purpose       : resolve external functions.  This function is to provide a 
//                 different method for function resolution, which doesn't rely on 
//                 the group.  Instead, the calling code is responsible for finding
//                 the expression that must be associated with this function name, 
//                 and then passing both in.  This function then performs the attachment.
//
// This must be a separate phase from the setup in lexAndParseExpression function, as
// *all* the relevant expressions must be allocated and have gone thru their initial
// set up before the group can perform the next phase.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/10/2020
//-------------------------------------------------------------------------------
bool newExpression::attachFunctionNode(const std::string & funcName, const Teuchos::RCP<Xyce::Util::newExpression> expPtr)
{
  bool retval=true;

  externalExpressions_.push_back(expPtr);

  // should I use upper?
  std::string funcNameUpper=funcName;
  Xyce::Util::toUpper(funcNameUpper);

  if (funcOpMap_.find(funcNameUpper) != funcOpMap_.end())
  {
    std::vector<Teuchos::RCP<astNode<usedType> > > & tmpVec = funcOpMap_[funcNameUpper];
    for (int ii=0;ii<tmpVec.size();++ii) 
    { 
      if ( !(Teuchos::is_null(tmpVec[ii])) ) 
      {
#if 0
        std::cout << "newExpression::attachFunctionNode. Attaching " << funcName << " to " << expressionString_ <<std::endl;
#endif
        tmpVec[ii]->setNode(expPtr->getAst()); 
      }
      else { retval=false; }

      Teuchos::RCP<funcOp<usedType> > castedFuncPtr = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (tmpVec[ii]);
      if ( !(Teuchos::is_null(castedFuncPtr)) )
      {
        castedFuncPtr->setFuncArgs( expPtr->getFunctionArgOpVec() ); 
      }
      else { retval=false; }
    }
    externalDependencies_ = true;
  }
  else { retval=false; }
  return retval;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::attachParameterNode
//
// Purpose       : resolve external parameters.  This function is to provide a 
//                 different method for parameter resolution, which doesn't rely on 
//                 the group.  Instead, the calling code is responsible for finding
//                 the expression that must be associated with this parameter name, 
//                 and then passing both in.  This function then performs the attachment.
//
// This must be a separate phase from the setup in lexAndParseExpression function, as
// *all* the relevant expressions must be allocated and have gone thru their initial
// set up before the group can perform the next phase.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/10/2020
//-------------------------------------------------------------------------------
bool newExpression::attachParameterNode(const std::string & paramName, const Teuchos::RCP<Xyce::Util::newExpression> expPtr, bool isDotParam)
{
  bool retval=true;

  externalExpressions_.push_back(expPtr);

  std::string paramNameUpper=paramName;
  Xyce::Util::toUpper(paramNameUpper);
  std::vector<std::string>::iterator nameIter = std::find(paramNameVec_.begin(),paramNameVec_.end(), paramNameUpper);
  if ( nameIter != paramNameVec_.end() )
  {
    int index = std::distance(paramNameVec_.begin(),nameIter);
    Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[index]);
    parOp->setNode(expPtr->getAst());
    parOp->setIsAttached();
    if(isDotParam){ parOp->setIsDotParam(); }
    externalDependencies_ = true;
  }
  else { retval=false; }
  return retval;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::clear
//
// Purpose       : Empties/resets out everything in the newExpression class, 
//                 and puts the class in the state that it should be in
//                 prior to calling the function lexAndParseExpression
//
// Special Notes : I am not sure this function needs to exist and it is hard to maintain reliably.
// Scope         :
// Creator       : Eric Keiter
// Creation Date : ??
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
  srcAstNodeVec_.clear();

  expressionString_ = std::string("");
  parsed_ = false;
  derivsSetup_ = false;
  astArraysSetup_ = false;
  bpTol_ = 0.0;
  timeStep_ = 0.0;
  timeStepAlpha_ = 0.0;
  timeStepPrefac_ = 0.0;
  numDerivs_ = 0;
  traditionalParse_ = true;
  externalDependencies_ = false;
  isTimeDependent_ = false;
  isTempDependent_ = false;
  isVTDependent_ = false;
  isFreqDependent_ = false;
  isGminDependent_ = false;

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

  derivIndexVec_.clear();

  return;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::make_constant
//
// Purpose       : This was originally needed for the old API, but it is useful 
//                 for the new expression library as well.   It applies a 
//                 specified value to a specified parameter.
//
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/??/2020
//-------------------------------------------------------------------------------
bool newExpression::make_constant (std::string const & var, usedType const & val, bool isDotParam)
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
    parOp->setIsConstant();
    if(isDotParam){ parOp->setIsDotParam(); }
    retval=true;
  }
  else
  {
    std::cout << "newExpression::make_constant  ERROR.  Could not find parameter " << tmpParName <<std::endl;
  }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::make_var
//
// Purpose       : Needed for the old API.   This sets the "var" boolean flag on 
//                 a specified parameter. In the old API, this means two things:
//                   (1) it is a global parameter rather than a regular parameter.
//                   (2) it should have derivatives computed.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : ??
//-------------------------------------------------------------------------------
bool newExpression::make_var (std::string const & var, bool isDotParam)
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
    parOp->setIsVar();
    if(isDotParam){ parOp->setIsDotParam(); }
    retval = true; // just means we found it
  }
  else
  {
    std::cout << "newExpression::make_var  ERROR.  Could not find parameter " << tmpParName <<std::endl;
  }

  return retval;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setupDerivatives_
//
// Purpose       : this is yet another phase in setup, which must be called after
//                 external entities (if they exist) have been resolved.
//
//                 If they are not, then there
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
    // Solution: only differentiate params that have their "setIsVar" boolean set to true.
    //

    for (int ii=0;ii<paramOpVec_.size();ii++)
    {
      Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
      if (parOp->getIsVar())
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
  }

  derivsSetup_ = true;
}

#define NEW_EXP_OUTPUT_ARRAY(VECTOR) \
  if ( !(VECTOR.empty()) )  { \
  os << #VECTOR << " (size="<<VECTOR.size()<<"):" << std::endl; \
  for (int ii=0;ii<VECTOR.size();ii++) \
  { \
    os << ii << " "; \
    VECTOR[ii]->output(os,0); \
  } }

//-------------------------------------------------------------------------------
// Function      : newExpression::outputVariousAstArrays_
// Purpose       : debug output
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/1/2020
//-------------------------------------------------------------------------------
void newExpression::outputVariousAstArrays( std::ostream & os )
{
  os << "Various arrays for expression: " << expressionString_ <<std::endl;
  if (externalDependencies_) os << "externalDependencies_ = true" <<std::endl;
  else os << "externalDependencies_ = false" <<std::endl;

NEW_EXP_OUTPUT_ARRAY(paramOpVec_)
NEW_EXP_OUTPUT_ARRAY(unresolvedParamOpVec_)
NEW_EXP_OUTPUT_ARRAY(funcOpVec_)
NEW_EXP_OUTPUT_ARRAY(voltOpVec_)
NEW_EXP_OUTPUT_ARRAY(currentOpVec_)
NEW_EXP_OUTPUT_ARRAY(powerOpVec_)
NEW_EXP_OUTPUT_ARRAY(internalDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY(dnoNoiseDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY(dniNoiseDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY(oNoiseOpVec_)
NEW_EXP_OUTPUT_ARRAY(iNoiseOpVec_)
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setupVariousAstArrays_
// Purpose       : see below
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 12/26/2019
//-------------------------------------------------------------------------------
void newExpression::setupVariousAstArrays_()
{
#if 0
  std::cout << "Array sizes BEFORE update:" <<std::endl;
  outputVariousAstArrays(std::cout);
#endif

  if (externalDependencies_)
  {

    // 1. setup arrays that require full AST traversal:
    paramOpVec_.clear();
    funcOpVec_.clear();
    voltOpVec_.clear();
    currentOpVec_.clear();
    internalDevVarOpVec_.clear();

    if( !(Teuchos::is_null(astNodePtr_)) )
    {
      if (astNodePtr_->paramType())
      {
        if ( !(astNodePtr_->getFunctionArgType()) )  // parameters are occasionally function arguments.  Don't include those
        {
          paramOpVec_.push_back(astNodePtr_);
        }
      }
      if (astNodePtr_->funcType())    { funcOpVec_.push_back(astNodePtr_); }
      if (astNodePtr_->voltageType()) { voltOpVec_.push_back(astNodePtr_); }
      if (astNodePtr_->currentType()) { currentOpVec_.push_back(astNodePtr_); }
      if (astNodePtr_->internalDeviceVarType()) { internalDevVarOpVec_.push_back(astNodePtr_); }

      if (astNodePtr_->dnoNoiseVarType()) { dnoNoiseDevVarOpVec_.push_back(astNodePtr_); }
      if (astNodePtr_->dniNoiseVarType()) { dniNoiseDevVarOpVec_.push_back(astNodePtr_); }
      if (astNodePtr_->oNoiseType())      { oNoiseOpVec_.push_back(astNodePtr_); }
      if (astNodePtr_->iNoiseType())      { iNoiseOpVec_.push_back(astNodePtr_); }

      opVectors_.isTimeDependent = isTimeDependent_;
      opVectors_.isTempDependent = isTempDependent_;
      opVectors_.isVTDependent = isVTDependent_;
      opVectors_.isFreqDependent = isFreqDependent_;
      opVectors_.isGminDependent = isGminDependent_;

      astNodePtr_->getInterestingOps( opVectors_  );

      if (opVectors_.isTimeDependent) isTimeDependent_ = true;
      if (opVectors_.isTempDependent) isTempDependent_ = true;
      if (opVectors_.isVTDependent  ) isVTDependent_   = true;
      if (opVectors_.isGminDependent) isGminDependent_ = true;
    }

#if 0
    funcNameVec_.clear();
    //funcOpMap_.clear();
    for (int ii=0;ii<funcOpVec_.size();++ii)
    {
      Teuchos::RCP<funcOp<usedType> > functionOp = Teuchos::rcp_static_cast<funcOp<usedType> > (funcOpVec_[ii]);
      std::string tmp = functionOp->getName();
      Xyce::Util::toUpper(tmp);
      //funcOpMap_[tmp].push_back(funcOpVec_[ii]);
      funcNameVec_.push_back(tmp);
    }
    
    //paramNameVec_.clear();
    //paramOpMap_.clear();
    for (int ii=0;ii<paramOpVec_.size();++ii)
    {
      Teuchos::RCP<paramOp<usedType> > parameterOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
      std::string tmp = parameterOp->getName();
      Xyce::Util::toUpper(tmp);
      //paramOpMap_[tmp].push_back(paramOpVec_[ii]);
      //paramNameVec_.push_back(tmp);
    }
#endif

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

    // 2. setup arrays that require traversal of expression objects (rather than AST nodes)
    //    For specials, this is more appropriate the AST traversal, as there will be at
    //    most one "time" node per expression object.   I think this should be a faster 
    //    traversal, generally.

    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getTimeNodes(timeOpVec_); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getTempNodes(tempOpVec_); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getVtNodes(vtOpVec_); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getFreqNodes(freqOpVec_); }
    for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getGminNodes(gminOpVec_); }
 
    isTimeDependent_ = !(timeOpVec_.empty());
    isTempDependent_ = !(tempOpVec_.empty());
    isVTDependent_   = !(vtOpVec_.empty());
    isFreqDependent_ = !(freqOpVec_.empty());
    isGminDependent_ = !(gminOpVec_.empty());

#if 0
    std::cout << "Array sizes AFTER update:" <<std::endl;
    outputVariousAstArrays(std::cout);
#endif
  }

  // check if this is expression is a constant
 
  bool noVariableParams = paramOpVec_.empty();
  if (!noVariableParams)
  {
    noVariableParams = true;
    for (int ii=0;ii<paramOpVec_.size();ii++)
    {
      Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
      if ( (parOp->getIsVar()) ) { noVariableParams = false; }
    }
  }

  if (
    !isTimeDependent_  &&
    !isTempDependent_  &&
    !isVTDependent_    &&
    !isFreqDependent_  &&
    !isGminDependent_  &&  // not relevant
    noVariableParams &&
    (unresolvedParamOpVec_.empty()) &&
    //(funcOpVec_.empty()) &&  // not relevant
    (voltOpVec_.empty()) &&
    (currentOpVec_.empty()) &&
    (powerOpVec_.empty()) &&
    (internalDevVarOpVec_.empty()) &&
    (dnoNoiseDevVarOpVec_.empty()) &&
    (dniNoiseDevVarOpVec_.empty()) &&
    (oNoiseOpVec_.empty()) &&
    (iNoiseOpVec_.empty())
      )
    {
      isConstant_ = true;
    }

#if 0
  if (isConstant_)
  {
    std::cout << "expression: " << expressionString_ << " is constant" << std::endl;
  }
  else
  {
    std::cout << "expression: " << expressionString_ << " is constant" << std::endl;
  }
#endif

  astArraysSetup_ = true;
};

//-------------------------------------------------------------------------------
// Function      : newExpression::getValuesFromGroup
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
void newExpression::getValuesFromGroup()
{
  // get solution values we need from the group
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
    if ( !(group_->getSolutionVal(currOp->getCurrentDevice(),val) ) ) // ERK.  reconsider the logic
    {
      group_->getCurrentVal(currOp->getCurrentDevice(),val);
    }
    currOp->setCurrentVal ( val );
  }

#if 0
  // don't remember why I wrote this 2x.
  for (int ii=0;ii<currentOpVec_.size();ii++)
  {
    Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);

    usedType val;
    if ( !(group_->getSolutionVal(currOp->getCurrentDevice(),val) ) ) // ERK.  reconsider the logic
    {
      group_->getCurrentVal(currOp->getCurrentDevice(),val);
    }
    currOp->setCurrentVal ( val );
  }
#endif

  for (int ii=0;ii<internalDevVarOpVec_.size();ii++)
  {
    Teuchos::RCP<internalDevVarOp<usedType> > intVarOp = Teuchos::rcp_static_cast<internalDevVarOp<usedType> > (internalDevVarOpVec_[ii]);

    usedType val;
    group_->getInternalDeviceVar(intVarOp->getInternalVarDevice(),val);
    intVarOp->setInternalDeviceVar ( val );
  }

  for (int ii=0;ii<paramOpVec_.size();++ii)
  {
    Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);

    // the "isVar" boolean currently serves two purposes.  If it is true that means:
    //
    // (1) we want derivatives w.r.t. it.
    // (2) it should be considered a dynamic variable that gets its values externally
    //
    if ( parOp->getIsVar() )
    {
      usedType val;
      group_->getGlobalParameterVal(parOp->getName(),val);
      parOp->setValue(val);
    }
  }

  for (int ii=0;ii<dnoNoiseDevVarOpVec_.size();ii++)
  {
    Teuchos::RCP<dnoNoiseVarOp<usedType> > dnoOp = Teuchos::rcp_static_cast<dnoNoiseVarOp<usedType> > (dnoNoiseDevVarOpVec_[ii]);

    usedType val;
    group_->getDnoNoiseDeviceVar(dnoOp->getNoiseDevice(),val);
    dnoOp->setNoiseVar ( val );
  }

  for (int ii=0;ii<dniNoiseDevVarOpVec_.size();ii++)
  {
    Teuchos::RCP<dniNoiseVarOp<usedType> > dniOp = Teuchos::rcp_static_cast<dniNoiseVarOp<usedType> > (dniNoiseDevVarOpVec_[ii]);
    usedType val;
    group_->getDniNoiseDeviceVar(dniOp->getNoiseDevice(),val);
    dniOp->setNoiseVar ( val );
  }

  for (int ii=0;ii<oNoiseOpVec_.size();ii++)
  {
    Teuchos::RCP<oNoiseOp<usedType> > onoiseOp = Teuchos::rcp_static_cast<oNoiseOp<usedType> > (oNoiseOpVec_[ii]);
    usedType val;
    group_->getONoise(val);
    onoiseOp->setNoiseVar ( val );
  }

  for (int ii=0;ii<iNoiseOpVec_.size();ii++)
  {
    Teuchos::RCP<iNoiseOp<usedType> > inoiseOp = Teuchos::rcp_static_cast<iNoiseOp<usedType> > (iNoiseOpVec_[ii]);
    usedType val;
    group_->getINoise(val);
    inoiseOp->setNoiseVar ( val );
  }

  for (int ii=0;ii<powerOpVec_.size();ii++)
  {
    Teuchos::RCP<powerOp<usedType> > pwrOp = Teuchos::rcp_static_cast<powerOp<usedType> > (powerOpVec_[ii]);
    usedType val;
    group_->getPower (pwrOp->getPowerDevice(), val);
    pwrOp->setPowerVal ( val );
  }

  for (int ii=0;ii<timeOpVec_.size();ii++) { timeOpVec_[ii]->setValue(group_->getTime()); }
  for (int ii=0;ii<tempOpVec_.size();ii++) { tempOpVec_[ii]->setValue(group_->getTemp()); } // Conversion to correct units in group 
  for (int ii=0;ii<vtOpVec_.size();ii++)   { vtOpVec_[ii]->setValue(group_->getVT()); }
  for (int ii=0;ii<freqOpVec_.size();ii++) { freqOpVec_[ii]->setValue(group_->getFreq()); }
  for (int ii=0;ii<gminOpVec_.size();ii++) { gminOpVec_[ii]->setValue(group_->getGmin()); }

  bpTol_ = group_->getBpTol();
  startingTimeStep_ = group_->getStartingTimeStep();
  finalTime_ = group_->getFinalTime();

  int srcSize = srcAstNodeVec_.size();
  for (int ii=0;ii< srcSize; ii++) 
  { 
    (*(srcAstNodeVec_[ii]))->setBreakPointTol(bpTol_); 
    (*(srcAstNodeVec_[ii]))->setStartingTimeStep(startingTimeStep_);
    (*(srcAstNodeVec_[ii]))->setFinalTime(finalTime_);
  }

  timeStep_ = group_->getTimeStep ();
  timeStepAlpha_ = group_->getTimeStepAlpha ();
  timeStepPrefac_ = group_->getTimeStepPrefac ();
}

//-------------------------------------------------------------------------------
// these two functions return int error codes in the original expression library
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Function      : newExpression::evaluate
// Purpose       : evaluates the expression, including derivatives
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
int newExpression::evaluate (usedType &result, std::vector< usedType > &derivs)
{
  int retVal=0;
  if (parsed_)
  {
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
    exit(0);
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
// Function      : newExpression::evaluateFunction
// Purpose       : evaluates the expression without derivatives
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
int newExpression::evaluateFunction (usedType &result)
{
  int retVal=0;
  if (parsed_)
  {
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }

    if ( !(unresolvedFuncOpVec_.empty()) )
    {
      std::cout << "ERROR.  Unresolved functions in expression " << expressionString_ <<std::endl;
      for(int ii=0;ii<unresolvedFuncOpVec_.size();++ii)
      {
        std::cout << "unresolvedFuncOpVec_[" << ii << "] = " << unresolvedFuncOpVec_[ii]->getName() <<std::endl;
      }
    }

    if (!isConstant_) 
    {
      getValuesFromGroup();
    }

    if (!isConstant_ || !evaluateFunctionCalledBefore_)
    {
      result = astNodePtr_->val();
      evaluateFunctionCalledBefore_ = true;

      // ERK: fix this failsafe properly for std::complex 
      if (isnan(std::real(result))) { result = 0.0; }
      if (isinf(std::real(result))) { result = 1.0e+20; }

      savedResult_ = result;
    }
    else
    {
      result = savedResult_;
    }
  }
  else
  {
    std::cout << "Error.  Expression " << expressionString_ << " is not parsed yet" << std::endl;
    exit(0);
  }

  return retVal;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::getBreakPoints
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/23/2020
//-------------------------------------------------------------------------------
bool newExpression::getBreakPoints (std::vector<Xyce::Util::BreakPoint> & breakPointTimes )
{
  int srcSize = srcAstNodeVec_.size();
  for (int ii=0;ii< srcSize; ii++) 
  { 
    (*(srcAstNodeVec_[ii]))->getBreakPoints(breakPointTimes); 
  }

#if 0
  if (srcSize>0)
  {
    std::cout << "newExpression::getBreakPoints. Expression " << expressionString_ << "  Number of breakpoints = " << breakPointTimes.size() <<std::endl;
    for (int ii=0;ii<breakPointTimes.size();ii++)
    {
      std::cout << "bp["<<ii<<"] = " << breakPointTimes[ii].value() <<std::endl;
    }
  }
#endif

  return true;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setFunctionArgStringVec
// Purpose       : 
// Special Notes : this must be set prior to lexAndParseExpression.
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
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

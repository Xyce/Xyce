//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-------------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <cmath>

#include <newExpression.h>
#include <N_ERH_Message.h>

//-------------------------------------------------------------------------------
// Expression Lexer/Parser header stuff
//
// Grrrr.  Stupid bison 2.4 stopped putting the pre-prologue into the header.
// need this forward declaration
namespace Xyce {
namespace Util {

template <typename ScalarT>
class tableArgs
{
  public:
    std::string * keywordPtr;
    Teuchos::RCP<astNode<ScalarT> > node;
};

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

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to true then logical AND is &&, logical
// OR is || and ^ is a synonym for exponentiation.  The default is false.
bool useHspiceMath;

// This value is derived from the -hspice-ext command line option.  It is
// set, based on that command line option, in the constructor for the
// IO::ParsingMgr class.  If set to true, then the separator is a period 
// character.  If false, the it is a colon.
bool useHspiceSeparator;
char separator;

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
  Xyce::dout() << "lexAndParseExpression for " << expressionString_ <<std::endl;
#endif
  if (traditionalParse_)
  {
    std::stringstream expressionStringStream ( expressionString_ );

    // if the expressionString_ is empty, bison will throw an error.
    // Plus, no point in parsing an empty string.
    if(expressionString_.empty())
    {
      parsed_ = true;  // ERK.   If this is false, and someone calls "evaluate", triggers a fatal error. So let it be true.
    }
    else
    {
      Xyce::Util::ExpressionLexer expLexer(expressionString_, &expressionStringStream);
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
      Teuchos::RCP<paramOp<usedType> > garbageParamOpPtr = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

      functionArgOpVec_.clear();
      functionArgOpVec_.resize(stringArgsSize,(garbageParamOpPtr));
      for (int ii=0;ii<stringArgsSize;++ii)
      {
        paramIter = std::find(paramNameVec_.begin(),paramNameVec_.end(), functionArgStringVec_[ii]);
        if (paramIter != paramNameVec_.end()) // found it
        {
          int index = std::distance(paramNameVec_.begin(),paramIter);
          functionArgOpVec_[ii] = paramOpVec_[index];
          paramNameVec_.erase(paramNameVec_.begin()+index);
          paramOpVec_.erase(paramOpVec_.begin()+index);
          paramOpMap_.erase( functionArgStringVec_[ii] );
          functionArgOpVec_[ii]->setFunctionArgType();
        }
      }
    }
  }

  // set up shallow dependence flags.
  {
    isShallowTimeDependent_ = isTimeDependent_;
    isShallowTempDependent_ = isTempDependent_;
    isShallowVTDependent_ = isVTDependent_;
    isShallowFreqDependent_ = isFreqDependent_;
    isShallowGminDependent_ = isGminDependent_;
  }

  isVariableDependent_ = !(globalParamNameVec_.empty()); // this is not reliable this stage, because the globalParamNameVec is always empty at this point.
  isVoltageNodeDependent_ = !(voltNameVec_.empty());
  isDeviceCurrentDependent_ = !(currentNameVec_.empty());
  isLeadCurrentDependent_ = !(leadCurrentNameVec_.empty());
  isLeadCurrentDependentExcludeBsrc_ = !(leadCurrentExcludeBsrcNameVec_.empty());

  // if dependent on a special, add relevant specials node to the relevant specials vector
  if(isTimeDependent_) // ERK:  should there be a separate boolean for dtDependent?
  {
    timeOpVec_.push_back(timeNodePtr_);
    dtOpVec_.push_back(dtNodePtr_);
  }
  if(isTempDependent_) { tempOpVec_.push_back(tempNodePtr_); }
  if(isVTDependent_) { vtOpVec_.push_back(vtNodePtr_); }
  if(isFreqDependent_) { freqOpVec_.push_back(freqNodePtr_); }
  if(isGminDependent_) { gminOpVec_.push_back(gminNodePtr_); }

#if 0
  dumpParseTree(Xyce::dout());
#endif

  // let AC analysis know if it needs to produce RF param output
  if ( !(Teuchos::is_null(group_)) )
  {
    if  ( !(sparamOpVec_.empty()) ) { group_->setRFParamsRequested(std::string("S")); }
    if  ( !(yparamOpVec_.empty()) ) { group_->setRFParamsRequested(std::string("Y")); }
    if  ( !(zparamOpVec_.empty()) ) { group_->setRFParamsRequested(std::string("Z")); }
  }

  unresolvedParamNameVec_ = paramNameVec_;
  unresolvedFuncNameVec_ = funcNameVec_;

  checkIsConstant_();
  astArraysSetup_ = true;

#if 0
  // this is a test to ensure that adding a "global" AST layer works for all use cases.  
  // I turn this on for unit tests, optionally.  But it should only be "on" for unit 
  // testing and not otherwise.
  if(parsed_)
  {
    setAsGlobal();
  }
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
bool newExpression::attachFunctionNode(
    const std::string & funcName,
    const Teuchos::RCP<Xyce::Util::newExpression> expPtr)
{
  bool retval=false;

  if ( !(Teuchos::is_null(expPtr)) )
  {
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
          if ( !(Teuchos::is_null( expPtr->getAst() )))
          {
            tmpVec[ii]->setNode(expPtr->getAst());
            Teuchos::RCP<funcOp<usedType> > castedFuncPtr = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (tmpVec[ii]);
            if ( !(Teuchos::is_null(castedFuncPtr)) )
            {
              castedFuncPtr->setFuncArgs( expPtr->getFunctionArgOpVec() );

              int size1 = castedFuncPtr->getFuncArgs().size();
              int size2 = expPtr->getFunctionArgOpVec().size();
              // getFunctionArgStringVec
              if (size1 != size2)
              {
                std::string errMsg = "Wrong number of arguments for user defined function " + castedFuncPtr->getName() + "(";
                for (int ii=0; ii<  expPtr->getFunctionArgStringVec().size();ii++)
                {
                  errMsg += expPtr->getFunctionArgStringVec()[ii];
                  if (size2 > 1 && ii < size2-1) { errMsg += ","; }
                }
                errMsg += ") in expression " + expressionString_;
                Xyce::Report::UserError() << errMsg;
              }

              castedFuncPtr->setSdtArgs( expPtr->getLocalSdtOpVec() );
              castedFuncPtr->setDdtArgs( expPtr->getLocalDdtOpVec() );
              retval=true;

              // remove from the unresolvedFunction container.  Currently this container isn't maintained or used.
            }
          }
        }
      }
      externalDependencies_ = true;
      astArraysSetup_ = false;

      isVariableDependent_      = isVariableDependent_ || expPtr->getVariableDependent();
      isVoltageNodeDependent_   = isVoltageNodeDependent_ || expPtr->getVoltageNodeDependent();
      isDeviceCurrentDependent_ = isDeviceCurrentDependent_ || expPtr->getDeviceCurrentDependent();
      isLeadCurrentDependent_   = isLeadCurrentDependent_ || expPtr->getLeadCurrentDependent();
      isLeadCurrentDependentExcludeBsrc_ = isLeadCurrentDependentExcludeBsrc_ || expPtr->getLeadCurrentDependentExcludeBsrc();

      isTimeDependent_ = isTimeDependent_ || expPtr->getTimeDependent();
      isTempDependent_ = isTempDependent_ || expPtr->getTempDependent();
      isVTDependent_ = isVTDependent_ || expPtr->getVTDependent();
      isFreqDependent_ = isFreqDependent_ || expPtr->getFreqDependent();
    }
    else { retval=false; }
  }

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
bool newExpression::attachParameterNode(
    const std::string & paramName,
    const Teuchos::RCP<Xyce::Util::newExpression> expPtr,
    enumParamType type)
{
  bool retval=false;

  if ( !(Teuchos::is_null(expPtr)) )
  {
    externalExpressions_.push_back(expPtr);
    std::string paramNameUpper=paramName;
    Xyce::Util::toUpper(paramNameUpper);
    if ( paramOpMap_.find( paramNameUpper ) != paramOpMap_.end() )
    {
      std::vector<Teuchos::RCP<astNode<usedType> > > & nodeVec = paramOpMap_[paramNameUpper];
      int size = nodeVec.size();
      for (int ii=0;ii<size;ii++)
      {
        Teuchos::RCP<astNode<usedType> > & node  = nodeVec[ii];
        Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (node);
        if ( !(Teuchos::is_null( expPtr->getAst() )))
        {
          parOp->setNode(expPtr->getAst());
          parOp->setIsAttached();
          parOp->setParamType(type);
          externalDependencies_ = true;
          astArraysSetup_ = false;
          retval=true;

#if 1
          std::vector<std::string>::iterator it = std::find(unresolvedParamNameVec_.begin(), unresolvedParamNameVec_.end(), paramNameUpper);
          if (it != unresolvedParamNameVec_.end())
          {
            int index = std::distance(unresolvedParamNameVec_.begin(),it);
            unresolvedParamNameVec_.erase(unresolvedParamNameVec_.begin()+index);
          }

          if (type == DOT_GLOBAL_PARAM)
          {
            it = std::find(globalParamNameVec_.begin(), globalParamNameVec_.end(), paramNameUpper);
            if (it == globalParamNameVec_.end())
            {
              globalParamNameVec_.push_back(paramNameUpper);
              isVariableDependent_ = !(globalParamNameVec_.empty());
            }
          }
          else
          {
            it = std::find(globalParamNameVec_.begin(), globalParamNameVec_.end(), paramNameUpper);
            if (it != globalParamNameVec_.end())
            {
              int index = std::distance(globalParamNameVec_.begin(),it);
              globalParamNameVec_.erase(globalParamNameVec_.begin()+index);
              isVariableDependent_ = !(globalParamNameVec_.empty());
            }
          }
#endif
          isVariableDependent_      = isVariableDependent_ || expPtr->getVariableDependent();
          isVoltageNodeDependent_   = isVoltageNodeDependent_ || expPtr->getVoltageNodeDependent();
          isDeviceCurrentDependent_ = isDeviceCurrentDependent_ || expPtr->getDeviceCurrentDependent();
          isLeadCurrentDependent_   = isLeadCurrentDependent_ || expPtr->getLeadCurrentDependent();
          isLeadCurrentDependentExcludeBsrc_ = isLeadCurrentDependentExcludeBsrc_ || expPtr->getLeadCurrentDependentExcludeBsrc();

          isTimeDependent_ = isTimeDependent_ || expPtr->getTimeDependent();
          isTempDependent_ = isTempDependent_ || expPtr->getTempDependent();
          isVTDependent_ = isVTDependent_ || expPtr->getVTDependent();
          isFreqDependent_ = isFreqDependent_ || expPtr->getFreqDependent();
        }
      }
    }
  }
  return retval;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::clear
//
// Purpose       : Empties/resets out everything in the newExpression class,
//                 and puts the class in the state that it should be in
//                 prior to calling the function lexAndParseExpression
//
// Special Notes : This function is hard to maintain reliably. Unfortunately, I
//                 can't (yet?) remove it
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : ??
//-------------------------------------------------------------------------------
void newExpression::clear ()
{
  // copied from destructor

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
  unresolvedParamNameVec_.clear();
  paramOpVec_.clear();
  paramOpMap_.clear();

  funcNameVec_.clear();
  unresolvedFuncNameVec_.clear();
  funcOpVec_.clear();
  funcOpMap_.clear();

  voltNameVec_.clear();
  voltOpVec_.clear();
  voltOpMap_.clear();

  currentNameVec_.clear();
  currentOpVec_.clear();
  currentOpMap_.clear();

  leadCurrentNameVec_.clear();
  leadCurrentExcludeBsrcNameVec_.clear();
  leadCurrentOpVec_.clear();

  bsrcCurrentOpVec_.clear();

  powerOpVec_.clear();

  internalDevVarOpVec_.clear();

  dnoNoiseDevVarOpVec_.clear();

  dniNoiseDevVarOpVec_.clear();

  oNoiseOpVec_.clear();

  iNoiseOpVec_.clear();

  sdtOpVec_.clear();
  ddtOpVec_.clear();
  phaseOpVec_.clear();
  sparamOpVec_.clear();
  yparamOpVec_.clear();
  zparamOpVec_.clear();
  agaussOpVec_.clear();
  localAgaussOpVec_.clear();
  gaussOpVec_.clear();
  localGaussOpVec_.clear();
  aunifOpVec_.clear();
  localAunifOpVec_.clear();
  unifOpVec_.clear();
  localUnifOpVec_.clear();
  randOpVec_.clear();
  localRandOpVec_.clear();
  twoArgLimitOpVec_.clear();
  localTwoArgLimitOpVec_.clear();

  srcAstNodeVec_.clear();
  stpAstNodeVec_.clear();
  compAstNodeVec_.clear();
  limitAstNodeVec_.clear();

  timeOpVec_.clear();
  dtOpVec_.clear();
  tempOpVec_.clear();
  vtOpVec_.clear();
  freqOpVec_.clear();
  gminOpVec_.clear();

  externalExpressions_.clear();

  derivIndexVec_.clear();
  derivNodeIndexMap_.clear();

  return;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::make_constant
//
// Purpose       : This function applies a specified value to a specified
//                 parameter, and also sets the "constant" boolean to true.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/??/2020
//-------------------------------------------------------------------------------
bool newExpression::make_constant (
    std::string const & var,
    usedType const & val,
    enumParamType type)
{
  std::string paramNameUpper = var;
  Xyce::Util::toUpper(paramNameUpper);
  bool retval=false;

  if ( paramOpMap_.find( paramNameUpper ) != paramOpMap_.end() )
  {
    std::vector<Teuchos::RCP<astNode<usedType> > > & nodeVec = paramOpMap_[paramNameUpper];
    int size = nodeVec.size();
    for (int ii=0;ii<size;ii++)
    {
      Teuchos::RCP<astNode<usedType> > & node  = nodeVec[ii];
      Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (node);
      parOp->setValue(val);
      parOp->setIsConstant();
      parOp->setParamType(type);
    }
    retval=true;

    std::vector<std::string>::iterator it = std::find(unresolvedParamNameVec_.begin(), unresolvedParamNameVec_.end(), paramNameUpper);
    if (it != unresolvedParamNameVec_.end())
    {
      int index = std::distance(unresolvedParamNameVec_.begin(),it);
      unresolvedParamNameVec_.erase(unresolvedParamNameVec_.begin()+index);
    }

    if (type == DOT_GLOBAL_PARAM)
    {
      it = std::find(globalParamNameVec_.begin(), globalParamNameVec_.end(), paramNameUpper);
      if (it == globalParamNameVec_.end())
      {
        globalParamNameVec_.push_back(paramNameUpper);
        isVariableDependent_ = !(globalParamNameVec_.empty());
      }
    }
    else
    {
      it = std::find(globalParamNameVec_.begin(), globalParamNameVec_.end(), paramNameUpper);
      if (it != globalParamNameVec_.end())
      {
        int index = std::distance(globalParamNameVec_.begin(),it);
        globalParamNameVec_.erase(globalParamNameVec_.begin()+index);
        isVariableDependent_ = !(globalParamNameVec_.empty());
      }
    }
    checkIsConstant_();
  }
  else
  {
    Xyce::Report::UserError()
      << "newExpression::make_constant  ERROR.  Could not find parameter "
      << paramNameUpper
      << " in expression: " << expressionString_ <<std::endl;
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
  if (!derivsSetup_)
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
        Teuchos::RCP<voltageOp<usedType> > voltOp
          = Teuchos::rcp_static_cast<voltageOp<usedType> > (voltOpVec_[ii]);

        const std::string & node = voltOp->getVoltageNode();

        std::string tmp = node; Xyce::Util::toUpper(tmp);
        std::unordered_map<std::string, int>::iterator mapIter;
        mapIter = derivNodeIndexMap_.find(tmp);
        if (mapIter == derivNodeIndexMap_.end()) { derivNodeIndexMap_[tmp] = numDerivs_; numDerivs_++; }
        derivIndexVec_.push_back(derivIndexPair_(voltOpVec_[ii],derivNodeIndexMap_[tmp]));
      }

      for (int ii=0;ii<currentOpVec_.size();ii++)
      {
        Teuchos::RCP<currentOp<usedType> > currOp
          = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);
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
      // Unlike voltages and currents, params can be assigned to each
      // other, etc, and they can be *anything*.
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
      // Solution: only differentiate params that have their "setIsVar"
      // boolean set to true.
      //
      for (int ii=0;ii<paramOpVec_.size();ii++)
      {
        Teuchos::RCP<paramOp<usedType> > parOp
          = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
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
// Function      : newExpression::outputVariousAstArrays
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
NEW_EXP_OUTPUT_ARRAY(funcOpVec_)
NEW_EXP_OUTPUT_ARRAY(voltOpVec_)
NEW_EXP_OUTPUT_ARRAY(currentOpVec_)
NEW_EXP_OUTPUT_ARRAY(leadCurrentOpVec_)
NEW_EXP_OUTPUT_ARRAY(powerOpVec_)
NEW_EXP_OUTPUT_ARRAY(internalDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY(dnoNoiseDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY(dniNoiseDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY(oNoiseOpVec_)
NEW_EXP_OUTPUT_ARRAY(iNoiseOpVec_)
NEW_EXP_OUTPUT_ARRAY(sdtOpVec_)
NEW_EXP_OUTPUT_ARRAY(ddtOpVec_)
NEW_EXP_OUTPUT_ARRAY(stpAstNodeVec_)
NEW_EXP_OUTPUT_ARRAY(compAstNodeVec_)
NEW_EXP_OUTPUT_ARRAY(limitAstNodeVec_)
NEW_EXP_OUTPUT_ARRAY(phaseOpVec_)
NEW_EXP_OUTPUT_ARRAY(timeOpVec_)
NEW_EXP_OUTPUT_ARRAY(dtOpVec_)
}

#define NEW_EXP_OUTPUT_ARRAY_SIZE(VECTOR) \
  if ( !(VECTOR.empty()) )  { \
  os << #VECTOR << " (size="<<VECTOR.size()<<"):" << std::endl; \
  }

//-------------------------------------------------------------------------------
// Function      : newExpression::outputVariousAstArraySizes
// Purpose       : debug output
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/1/2020
//-------------------------------------------------------------------------------
void newExpression::outputVariousAstArraySizes( std::ostream & os )
{
  os << "Various arrays for expression: " << expressionString_ <<std::endl;
  if (externalDependencies_) os << "externalDependencies_ = true" <<std::endl;
  else os << "externalDependencies_ = false" <<std::endl;

NEW_EXP_OUTPUT_ARRAY_SIZE(paramOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(funcOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(voltOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(currentOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(leadCurrentOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(powerOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(internalDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(dnoNoiseDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(dniNoiseDevVarOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(oNoiseOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(iNoiseOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(sdtOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(ddtOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(stpAstNodeVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(compAstNodeVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(limitAstNodeVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(phaseOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(timeOpVec_)
NEW_EXP_OUTPUT_ARRAY_SIZE(dtOpVec_)
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setupVariousAstArrays
//
// Purpose       : This function builds up lists of all the relevant data
//                 structures required to understand the external depedencies of
//                 this expression.
//
//                 It is only necessary to call it when the expression has external
//                 dependencies.  The reason is that during parsing, all the
//                 information that can be determined during parsing is added to
//                 these data structures.  External dependencies, however,
//                 cannot be fully evaluated during the parse phase.
//
//                 Ideally, however, the data structutes below should be added to
//                 during the two "attach" functions:
//
//                   newExpression::attachFunctionNode
//                   newExpression::attachParameterNode
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 12/26/2019
//-------------------------------------------------------------------------------
void newExpression::setupVariousAstArrays()
{
  if (!astArraysSetup_)
  {
#if 0
    Xyce::dout() << "Array sizes BEFORE update:" <<std::endl;
    outputVariousAstArrays(Xyce::dout());
#endif
    if (externalDependencies_)
    {
      // setup arrays that require full AST traversal:
      paramOpVec_.clear();
      funcOpVec_.clear();
      voltOpVec_.clear();
      currentOpVec_.clear();
      leadCurrentOpVec_.clear();
      bsrcCurrentOpVec_.clear();
      powerOpVec_.clear();
      internalDevVarOpVec_.clear();
      dnoNoiseDevVarOpVec_.clear();
      dniNoiseDevVarOpVec_.clear();
      oNoiseOpVec_.clear();
      iNoiseOpVec_.clear();
      sdtOpVec_.clear();
      ddtOpVec_.clear();
      srcAstNodeVec_.clear();
      stpAstNodeVec_.clear();
      compAstNodeVec_.clear();
      limitAstNodeVec_.clear();
      phaseOpVec_.clear();
      sparamOpVec_.clear();
      yparamOpVec_.clear();
      zparamOpVec_.clear();

      agaussOpVec_.clear();
      gaussOpVec_.clear();
      aunifOpVec_.clear();
      unifOpVec_.clear();
      randOpVec_.clear();
      twoArgLimitOpVec_.clear();

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
        if (astNodePtr_->leadCurrentType()) { leadCurrentOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->bsrcCurrentType()) { bsrcCurrentOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->powerType()) { powerOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->internalDeviceVarType()) { internalDevVarOpVec_.push_back(astNodePtr_); }

        if (astNodePtr_->dnoNoiseVarType()) { dnoNoiseDevVarOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->dniNoiseVarType()) { dniNoiseDevVarOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->oNoiseType())      { oNoiseOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->iNoiseType())      { iNoiseOpVec_.push_back(astNodePtr_); }

        if (astNodePtr_->sdtType())      { sdtOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->ddtType())      { ddtOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->srcType())      { srcAstNodeVec_.push_back(astNodePtr_); }
        if (astNodePtr_->stpType())      { stpAstNodeVec_.push_back(astNodePtr_); }
        if (astNodePtr_->compType())      { compAstNodeVec_.push_back(astNodePtr_); }
        if (astNodePtr_->limitType())      { limitAstNodeVec_.push_back(astNodePtr_); }
        if (astNodePtr_->phaseType())    { phaseOpVec_.push_back(astNodePtr_); }

        if (astNodePtr_->sparamType())    { sparamOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->yparamType())    { yparamOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->zparamType())    { zparamOpVec_.push_back(astNodePtr_); }

        if (astNodePtr_->agaussType())    { agaussOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->gaussType())    { gaussOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->aunifType())    { aunifOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->unifType())    { unifOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->randType())    { randOpVec_.push_back(astNodePtr_); }
        if (astNodePtr_->twoArgLimitType())    { twoArgLimitOpVec_.push_back(astNodePtr_); }

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

      // populate some vectors that depend on the above AST traversal,
      // but that were not updated during it. (these could easily be included in the
      // traversal, so perhaps do that later)
      voltOpMap_.clear();
      voltNameVec_.clear(); // needed by "getVoltageNodes"
      for (int ii=0;ii<voltOpVec_.size();++ii)
      {
        Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (voltOpVec_[ii]);
        const std::string & node = voltOp->getVoltageNode();

        if ( voltOpMap_.find(node) == voltOpMap_.end() )
        {
          std::vector<Teuchos::RCP<astNode<usedType> > > vec;
          vec.push_back(voltOpVec_[ii]);
          voltOpMap_[node] = vec;
          voltNameVec_.push_back( node );
        }
        else
        {
          voltOpMap_[node].push_back(voltOpVec_[ii]);
        }
      }

      currentOpMap_.clear();
      currentNameVec_.clear(); // needed by "getDeviceCurrents"
      for (int ii=0;ii<currentOpVec_.size();++ii)
      {
        Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (currentOpVec_[ii]);
        std::string tmp = currOp->getCurrentDevice();

        if ( currentOpMap_.find(tmp) == currentOpMap_.end() )
        {
          std::vector<Teuchos::RCP<astNode<usedType> > > vec;
          vec.push_back(currentOpVec_[ii]);
          currentOpMap_[tmp] = vec;
          currentNameVec_.push_back( tmp );
        }
        else
        {
          currentOpMap_[tmp].push_back(currentOpVec_[ii]);
        }
      }

      leadCurrentNameVec_.clear();
      leadCurrentExcludeBsrcNameVec_.clear();
      for (int ii=0;ii<leadCurrentOpVec_.size();++ii)
      {
        Teuchos::RCP<leadCurrentOp<usedType> > leadCurrOp = Teuchos::rcp_static_cast<leadCurrentOp<usedType> > (leadCurrentOpVec_[ii]);
        std::string tmp = leadCurrOp->getLeadCurrentDevice();

        std::string tmpName = leadCurrentOpVec_[ii]->getName();

        std::vector<std::string>::iterator it = std::find(leadCurrentNameVec_.begin(), leadCurrentNameVec_.end(), tmpName);
        if (it == leadCurrentNameVec_.end())
        {
          leadCurrentNameVec_.push_back( tmpName );
          leadCurrentExcludeBsrcNameVec_.push_back( tmpName );
        }
      }

      // Bsrc's are special.  Depending on where in Xyce we are,
      // they should sometimes be included in lead currents, and sometimes not.
      for (int ii=0;ii<bsrcCurrentOpVec_.size();ii++)
      {
        std::string tmpName = bsrcCurrentOpVec_[ii]->getName();
        std::vector<std::string>::iterator it = std::find(leadCurrentNameVec_.begin(), leadCurrentNameVec_.end(), tmpName);
        if (it == leadCurrentNameVec_.end())
        {
          leadCurrentNameVec_.push_back( tmpName );
        }
      }

      // In at least some cases, what is really being requested is branch calculations, which can be either lead currents or power.
      for (int ii=0;ii<powerOpVec_.size();ii++)
      {
        std::string tmpName = powerOpVec_[ii]->getName();
        std::vector<std::string>::iterator it = std::find(leadCurrentNameVec_.begin(), leadCurrentNameVec_.end(), tmpName);
        if (it == leadCurrentNameVec_.end())
        {
          leadCurrentNameVec_.push_back( tmpName );
          leadCurrentExcludeBsrcNameVec_.push_back( tmpName );
        }
      }

      // 8/6/2020.  ERK.  I originally thought that paramNameVec_ wasn't used after parsing,
      // but I was wrong.  It needs updating here, because it is used by functions such as
      // make_const and make_var.  (those functions *should* use the map, instead)
      //
      // A vector is needed during parsing, b/c I need to keep the params in their
      // original order for some purposes. (I tried a map during parsing and it broke
      // some things, b/c maps change the order)
      //
      // But, once those issues are over with, a map should be better.  So, make_const and make_var
      // really should use the map instead.  Then once that happens, it won't be necessary to
      // update paramNameVec here.  But for now, to make those functions work properly, it is
      // what it is.
      //
      // Note, it just occured to me(8/7/2020):  the paramNameVec usage in make_var and make_const will
      // find the first occurance of a particular parameter.  If there are duplicates, it won't
      // find the others.  With the map object, it wouldn't have this problem.
      //
      // During parsing, there is code in the parser to ensure that there aren't any
      // duplicates in the paramOpVec_ or in the parmaNameVec_.  But once external nodes are attached,
      // it isn't possible to enforce this anymore.  So, at that stage, the paramNameVec_
      // could have duplicates in it, and make_var and make_const will break.
      //
      // 9/10/2020: Update.  The paramNameVec can indeed have duplicates in it.  But the
      // corresponding entries in the paramOpVec_ are usually duplicates as well.
      // So, just finding the first works OK, but it is wasteful.
      //
      // I believe this happening because the AST traversal
      // (in the function call astNodePtr_->getInterestingOps( opVectors_  ); above )
      // makes no attempt to  avoid duplicates.  It just pushes them all back onto the paramOpVec
      // as it traverses the tree.
      paramNameVec_.clear();
      paramOpMap_.clear();
      unresolvedParamNameVec_.clear();
      globalParamNameVec_.clear();

      std::unordered_map<std::string,std::vector<unsigned long int> > paramIDMap;

      for (int ii=0;ii<paramOpVec_.size();++ii)
      {
        Teuchos::RCP<paramOp<usedType> > parPtr = Teuchos::rcp_dynamic_cast<paramOp<usedType> > (paramOpVec_[ii]);
        std::string tmp = parPtr->getName();

        if (paramOpMap_.find(tmp) == paramOpMap_.end())
        {
          std::vector<Teuchos::RCP<astNode<usedType> > > vec;
          vec.push_back(paramOpVec_[ii]);
          paramOpMap_[tmp] = vec;

          std::vector<unsigned long int> IDvec;
          IDvec.push_back(parPtr->getId());
          paramIDMap[tmp] = IDvec;
        }
        else
        {
          unsigned long int id = parPtr->getId();
          std::vector<unsigned long int> & IDvec = paramIDMap[tmp];
          std::vector<unsigned long int>::iterator idIter = std::find(IDvec.begin(),IDvec.end(),id) ;
          if ( idIter == IDvec.end() )
          {
            paramOpMap_[tmp].push_back(paramOpVec_[ii]);
            IDvec.push_back(id);
          }
        }

        paramNameVec_.push_back(tmp);

        if( !(parPtr->getIsConstant())  && !(parPtr->getIsAttached()) )
        {
          std::string tmpName = paramOpVec_[ii]->getName();
          std::vector<std::string>::iterator it = std::find(unresolvedParamNameVec_.begin(), unresolvedParamNameVec_.end(), tmpName);
          if (it == unresolvedParamNameVec_.end())
          {
            unresolvedParamNameVec_.push_back( tmpName );
          }
        }

        if (  parPtr->getParamType() == DOT_GLOBAL_PARAM )
        {
          std::string tmpName = paramOpVec_[ii]->getName();
          std::vector<std::string>::iterator it = std::find(globalParamNameVec_.begin(), globalParamNameVec_.end(), tmpName);
          if (it == globalParamNameVec_.end())
          {
            globalParamNameVec_.push_back( tmpName );
          }
        }
      }

      // setup unresolvedFuncNameVec_
      unresolvedFuncNameVec_.clear();
      for (int ii=0;ii<funcOpVec_.size();ii++)
      {
        Teuchos::RCP<funcOp<usedType> > funPtr = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (funcOpVec_[ii]);

        if( !(funPtr->getNodeResolved()) || !(funPtr->getArgsResolved()) )
        {
          std::string tmpName = funcOpVec_[ii]->getName();
          std::vector<std::string>::iterator it = std::find(unresolvedFuncNameVec_.begin(), unresolvedFuncNameVec_.end(), tmpName);
          if (it == unresolvedFuncNameVec_.end())
          {
            unresolvedFuncNameVec_.push_back( tmpName );
          }
        }
      }

      // setup arrays that require traversal of expression objects (rather than AST nodes)
      // For specials, this is more appropriate the AST traversal, as there will be at
      // most one "time" node per expression object.   I think this should be a faster
      // traversal, generally.

      for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getTimeNodes(timeOpVec_); }
      for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getDtNodes(dtOpVec_); }
      for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getTempNodes(tempOpVec_); }
      for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getVtNodes(vtOpVec_); }
      for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getFreqNodes(freqOpVec_); }
      for (int ii=0;ii<externalExpressions_.size();ii++) { externalExpressions_[ii]->getGminNodes(gminOpVec_); }

      isTimeDependent_ = !( timeOpVec_.empty() && dtOpVec_.empty() );
      isTempDependent_ = !(tempOpVec_.empty());
      isVTDependent_   = !(vtOpVec_.empty());
      isFreqDependent_ = !(freqOpVec_.empty());
      isGminDependent_ = !(gminOpVec_.empty());

#if 0
      Xyce::dout() << "Array sizes AFTER update:" <<std::endl;
      outputVariousAstArrays(Xyce::dout());
#endif

      if ( !(Teuchos::is_null(group_)) )
      {
        if  ( !(sparamOpVec_.empty()) ) { group_->setRFParamsRequested(std::string("S")); }
        if  ( !(yparamOpVec_.empty()) ) { group_->setRFParamsRequested(std::string("Y")); }
        if  ( !(zparamOpVec_.empty()) ) { group_->setRFParamsRequested(std::string("Z")); }
      }
    }

    checkIsConstant_();
    astArraysSetup_ = true;
    groupSetup_ = false;
  }
};

#define NEW_EXP_ADD_TO_VEC1(VEC1,VEC2) \
  if ( !(VEC2.empty()) )  { VEC1.insert( VEC1.end(), VEC2.begin(), VEC2.end() ); }

#define NEW_EXP_ADD_TO_VEC2(VEC1,VEC2) \
  if ( !(VEC2.empty()) ) { for (int ii=0;ii<VEC2.size();ii++) { VEC1.push_back(VEC2[ii]); } }

//-------------------------------------------------------------------------------
// Function      : newExpression::checkIsConstant_
//
// Purpose       : This function is used by the N_UTL_Param function
//                 "isExpressionConstant", which is used to determine if an
//                 expression can be turned into an immutable data type such as
//                 DBLE or INT.
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 2020
//-------------------------------------------------------------------------------
void newExpression::checkIsConstant_()
{
  // check if this is expression is a constant
  if ( !(unresolvedParamNameVec_.empty()) )
  {
    isConstant_=false;
    return;
  }

  bool noVariableParams = paramOpVec_.empty();
  if (!noVariableParams)
  {
    if ( (unresolvedParamNameVec_.empty()) ) // only do this check if everything is resolved. Otherwise it is a pointless check
    {
      noVariableParams = true;
      for (int ii=0;ii<paramOpVec_.size();ii++)
      {
        Teuchos::RCP<paramOp<usedType> > parOp
          = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
        if (  parOp->getParamType() == DOT_GLOBAL_PARAM ) { noVariableParams = false;  break;}
      }
    }
  }

  if (
    !isTimeDependent_  &&
    !isTempDependent_  &&
    !isVTDependent_    &&
    !isFreqDependent_  &&
    !isGminDependent_  &&  // not relevant
    noVariableParams &&
    //(funcOpVec_.empty()) &&  // not relevant but based on comments above consistent with num_vars in the old library.  Leaving it out is the right thing to do and doesn't appear to break anything.
    (voltOpVec_.empty()) &&
    (currentOpVec_.empty()) &&
    (leadCurrentOpVec_.empty()) &&
    (powerOpVec_.empty()) &&
    (internalDevVarOpVec_.empty()) &&
    (dnoNoiseDevVarOpVec_.empty()) &&
    (dniNoiseDevVarOpVec_.empty()) &&
    (oNoiseOpVec_.empty()) &&
    (iNoiseOpVec_.empty()) &&
    ( (localAgaussOpVec_.empty())  )  &&
    ( (localGaussOpVec_.empty()) )  &&
    ( (localAunifOpVec_.empty()) )  &&
    ( (localUnifOpVec_.empty())  )  &&
    ( (localRandOpVec_.empty())  )  &&
    ( (localTwoArgLimitOpVec_.empty()) )
      )
    {
      isConstant_ = true;
    }
  else
  {
      isConstant_ = false;
  }

#if 0
  if (isConstant_)
  {
    Xyce::dout() << "expression: " << expressionString_
      << " is constant" << std::endl;
  }
  else
  {
    Xyce::dout() << "expression: " << expressionString_
      << " is constant" << std::endl;
  }
#endif
}

//-------------------------------------------------------------------------------
// Function      : newExpression::getValuesFromGroup_
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-------------------------------------------------------------------------------
bool newExpression::getValuesFromGroup_()
{
  bool noChange=true;
#if 0
  std::cout << "newExpression::getValuesFromGroup_() expression = "
    << expressionString_ << std::endl;
#endif

  // this function will get nearly everything, including voltages, currents, parameters, etc.
  noChange = group_->putValues(*this);

  // specials:
  if ( !(timeOpVec_.empty()) )
  {
    usedType oldtime = timeOpVec_[0]->getValue();
    usedType newtime = group_->getTime();
    if(oldtime != newtime)
    {
      for (int ii=0;ii<timeOpVec_.size();ii++) { timeOpVec_[ii]->setValue(newtime); }
      noChange = false;
    }
  }

  if ( !(dtOpVec_.empty()) )
  {
    usedType oldtimestep = dtOpVec_[0]->getValue();
    usedType newtimestep = group_->getTimeStep();
    if(oldtimestep != newtimestep)
    {
      for (int ii=0;ii<dtOpVec_.size();ii++) { dtOpVec_[ii]->setValue(newtimestep); }
      noChange = false;
    }
  }

  if ( !(tempOpVec_.empty()) )
  {
    usedType oldtemp = tempOpVec_[0]->getValue();
    usedType newtemp = group_->getTemp();

    if (overrideGroupTemperature_) { newtemp = overrideTemp_; }

    if(oldtemp != newtemp)
    {
      for (int ii=0;ii<tempOpVec_.size();ii++) { tempOpVec_[ii]->setValue(newtemp); }
      noChange = false;
    }
  }

  if ( !(vtOpVec_.empty()) )
  {
    usedType oldvt = vtOpVec_[0]->getValue();
    usedType newvt = group_->getVT ();
    if(oldvt != newvt)
    {
      for (int ii=0;ii<vtOpVec_.size();ii++)   { vtOpVec_[ii]->setValue(newvt); }
      noChange = false;
    }
  }

  if ( !(freqOpVec_.empty()) )
  {
    usedType oldfreq = freqOpVec_[0]->getValue();
    usedType newfreq = group_->getFreq();
    if(oldfreq != newfreq)
    {
      for (int ii=0;ii<freqOpVec_.size();ii++) { freqOpVec_[ii]->setValue(newfreq); }
      noChange = false;
    }
  }

  if ( !(gminOpVec_.empty()) )
  {
    usedType oldgmin = gminOpVec_[0]->getValue();
    usedType newgmin = group_->getGmin();
    if(oldgmin != newgmin)
    {
      for (int ii=0;ii<gminOpVec_.size();ii++) { gminOpVec_[ii]->setValue(newgmin); }
      noChange = false;
    }
  }

  // ERK: should the rest of these affect "noChange" boolean?
  if ( !(srcAstNodeVec_.empty())  || 
       !(stpAstNodeVec_.empty())  ||  
       !(compAstNodeVec_.empty()) ||  
       !(limitAstNodeVec_.empty()) )
  {
    bpTol_ = group_->getBpTol();
  }

  if ( !(srcAstNodeVec_.empty()) )
  {
    startingTimeStep_ = group_->getStartingTimeStep();
    finalTime_ = group_->getFinalTime();
    int srcSize = srcAstNodeVec_.size();
    for (int ii=0;ii< srcSize; ii++)
    {
      (srcAstNodeVec_[ii])->setBreakPointTol(bpTol_);
      (srcAstNodeVec_[ii])->setStartingTimeStep(startingTimeStep_);
      (srcAstNodeVec_[ii])->setFinalTime(finalTime_);
    }
  }

  if ( !(stpAstNodeVec_.empty()) )
  {
    int stpSize = stpAstNodeVec_.size();
    for (int ii=0;ii< stpSize; ii++)
    {
      (stpAstNodeVec_[ii])->setBreakPointTol(bpTol_);
    }
  }

  if ( !(compAstNodeVec_.empty()) )
  {
    int compSize = compAstNodeVec_.size();
    for (int ii=0;ii< compSize; ii++)
    {
      (compAstNodeVec_[ii])->setBreakPointTol(bpTol_);
    }
  }

  if ( !(limitAstNodeVec_.empty()) )
  {
    int limitSize = limitAstNodeVec_.size();
    for (int ii=0;ii< limitSize; ii++)
    {
      (limitAstNodeVec_[ii])->setBreakPointTol(bpTol_);
    }
  }

#if 0
  // these seem like they aren't used, commenting out.

  double oldTime_ = time_;
  time_ = group_->getTime();
  timeStep_ = group_->getTimeStep ();
  timeStepAlpha_ = group_->getTimeStepAlpha ();
  timeStepPrefac_ = group_->getTimeStepPrefac ();

  unsigned int oldStepNumber_ = stepNumber_;
  stepNumber_ = group_->getStepNumber ();
#endif

  // this probably never changes ... only set it 1x
  if ( !(phaseOpVec_.empty()) )
  {
    phaseOutputUsesRadians_ = group_->getPhaseOutputUsesRadians();
    for (int ii=0;ii<phaseOpVec_.size();ii++)
    {
      Teuchos::RCP<phaseOp<usedType> > phOp = Teuchos::rcp_static_cast<phaseOp<usedType> > (phaseOpVec_[ii]);
      phOp->setPhaseOutputUsesRadians( phaseOutputUsesRadians_ );
    }
  }

#if 0
  for (int ii=0;ii<twoArgLimitOpVec_.size();ii++)
  {
    Teuchos::RCP<twoArgLimitOp<usedType> > talOp = Teuchos::rcp_static_cast<twoArgLimitOp<usedType> > (twoArgLimitOpVec_[ii]);
  }
#endif

  return noChange;
}


//-------------------------------------------------------------------------------
// Function      : newExpression::setAsGlobal
// Purpose       : Adds an extra layer to the AST.  
//
// Special Notes : It is mostly harmless to call this on *any* expression.  
//                 However, it does add an extra layer of indirection to the AST, 
//                 so it is best to only apply this when it is really needed.  
//                 The main use case is when the expression is the RHS of a 
//                 .global_param statement.  Global params are variables, and can 
//                 be "set" by various operations in Xyce.  When they are set, the 
//                 new value needs to override whatever the RHS expression originally 
//                 was.  There is no guarentee that the RHS was a simple number; 
//                 it might have originally been a whole other expression.  
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/20/2021
//-------------------------------------------------------------------------------
void newExpression::setAsGlobal()
{
  Teuchos::RCP<globalParamLayerOp<usedType> > paramLayerPtr = Teuchos::RCP<globalParamLayerOp<usedType> >(new globalParamLayerOp<usedType>());
  paramLayerPtr->setNode(astNodePtr_);
  Teuchos::RCP<astNode<usedType> > newAstNodePtr = paramLayerPtr;
  setAstPtr(newAstNodePtr);
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void newExpression::setValue(usedType val)
{
  Teuchos::RCP<globalParamLayerOp<usedType> > paramLayerPtr = Teuchos::rcp_static_cast<globalParamLayerOp<usedType> > (astNodePtr_);
  paramLayerPtr->setValue(val);
}

//-------------------------------------------------------------------------------
// these two functions return int error codes in the original expression library
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Function      : newExpression::evaluate
// Purpose       : evaluates the expression, including derivatives
// Special Notes : returns "true" if expression has changed, false if not.
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-------------------------------------------------------------------------------
bool newExpression::evaluate (usedType &result, std::vector< usedType > &derivs)
{
  bool retVal=true;
  if (parsed_)
  {
    setupVariousAstArrays ();
    setupDerivatives_ ();

    if (!groupSetup_)
    {
      groupSetup_=group_->setupGroup(*this);
    }

#if 0
    Xyce::dout() << "Parse Tree for " << expressionString_ << std::endl;
    dumpParseTree(Xyce::dout());
#endif

#if 0
    // old, slower way
    retVal = evaluateFunction (result); // for now don't check anything beyond what evaluateFunction checks
    if (derivs.size() != numDerivs_) {derivs.clear(); derivs.resize(numDerivs_);}

    if ( !(Teuchos::is_null(astNodePtr_)) )
    {
      for (int ii=0;ii<derivIndexVec_.size();ii++) { derivIndexVec_[ii].first->setDerivIndex(derivIndexVec_[ii].second); }
      for (int ii=0;ii<numDerivs_;++ii) { derivs[ii] = astNodePtr_->dx(ii); }
      for (int ii=0;ii<derivIndexVec_.size();ii++) { derivIndexVec_[ii].first->unsetDerivIndex(); }
    }
#else
    //  new faster way, using dx2().  This allows all the calculations to be performed in a single AST traversal.
    bool noChange = getValuesFromGroup_(); // this was inside of evaluateFunction
    
    if (derivs.size() != numDerivs_) {derivs.clear(); derivs.resize(numDerivs_);}
    if ( !(Teuchos::is_null(astNodePtr_)) )
    {
      for (int ii=0;ii<derivIndexVec_.size();ii++) { derivIndexVec_[ii].first->setDerivIndex(derivIndexVec_[ii].second); }

      astNodePtr_->dx2(result,derivs);

      // this block was in evaluateFunction
      Util::fixNan(result);
      Util::fixInf(result);
      retVal = (result != savedResult_);
      savedResult_ = result;

      for (int ii=0;ii<derivIndexVec_.size();ii++) { derivIndexVec_[ii].first->unsetDerivIndex(); }
    }
#endif
  }
  else
  {
    Xyce::Report::UserError() << "Error.  Expression "
      << expressionString_
      << " was not successfully parsed." << std::endl;
  }

  for(int ii=0;ii<derivs.size();++ii)
  {
    Util::fixNan(derivs[ii]);// this was previously set to return 0.0, but the function returns 1.0e+50.  check this
    Util::fixInf(derivs[ii]);// this was previously for 1.0e+10, but this function call uses 1.0e+50.  check this.
  }

  // old expression library returns EXPRerrno, which is a static variable.
  // If it is zero, everything is cool.
  return retVal;
};

//-------------------------------------------------------------------------------
// Function      : newExpression::evaluateFunction
// Purpose       : evaluates the expression without derivatives
// Special Notes : returns "true" if expression has changed, false if not.
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-------------------------------------------------------------------------------
bool newExpression::evaluateFunction (usedType &result, bool efficiencyOn)
{
  bool retVal=true;
  if (parsed_)
  {
    setupVariousAstArrays ();

    if (!groupSetup_)
    {
      groupSetup_=group_->setupGroup(*this);
    }

#if 0
    Xyce::dout() << "newExpression::evaluateFunction. about to evaluate expression tree for " << expressionString_ << std::endl;
    dumpParseTree(Xyce::dout());
#endif

    bool noChange = getValuesFromGroup_();
    bool doTheEvaluation = true;
    if (efficiencyOn && noChange)
    {
      if ( evaluateFunctionCalledBefore_ )
      {
        doTheEvaluation = false;
      }
      else
      {
        doTheEvaluation = true;
        evaluateFunctionCalledBefore_  = true;
      }
    }

    if (doTheEvaluation)
    {
      if ( !(Teuchos::is_null(astNodePtr_)) )
      {
#if 1
        result = astNodePtr_->val();
#else
        // this is a test, to use with unit tests to ensure that the "result" evaluation in dx2 matches that of val().
        std::vector<usedType> derivs; // if empty, then dx2 should run OK.
        astNodePtr_->dx2(result, derivs);
#endif
        Util::fixNan(result);
        Util::fixInf(result);
      }
      retVal = (result != savedResult_);

#if 0
      Xyce::dout().width(20); Xyce::dout().precision(13);
      Xyce::dout().setf(std::ios::scientific);

      if (retVal)
      {
        Xyce::dout() << "newExpression::evaluateFunction. just evaluated expression tree for " << expressionString_ << " result = " << result << " retVal(changed) = true" << std::endl;
      }
      else
      {
        Xyce::dout() << "newExpression::evaluateFunction. just evaluated expression tree for " << expressionString_ << " result = " << result << " retVal(changed) = false" << std::endl;
      }

      dumpParseTree(Xyce::dout());
#endif



      savedResult_ = result;
    }
    else
    {
      retVal = false;
      result = savedResult_;
#if 0
      Xyce::dout()
        << "newExpression::evaluateFunction. just skipped evaluating the expression tree (b/c constant) for "
        << expressionString_ << " result = " << result << std::endl;
#endif
    }
  }
  else
  {
    Xyce::dout()
      << "Error.  Expression "
      << expressionString_
      << " is not parsed yet" << std::endl;
    exit(0);
  }

#if 0
  Xyce::dout()
    << "newExpression::evaluateFunction. just evaluated expression tree for "
    << expressionString_ << " result = " << result << std::endl;
#endif

  return retVal;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setupBreakPoints
//
// Purpose       : This is called when doing something
//                 like a .step around a transient simulation, at the beginning
//                 of each .step iteration.
//
// Special Notes : The tableOp class has a meaningful setupBreakPoints function.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 8/25/2020
//-------------------------------------------------------------------------------
void newExpression::setupBreakPoints ()
{
  if(isTimeDependent_)
  {
    int srcSize = srcAstNodeVec_.size();
    for (int ii=0;ii< srcSize; ii++)
    { (srcAstNodeVec_[ii])->setupBreakPoints(); }

    //int stpSize = stpAstNodeVec_.size();
    //for (int ii=0;ii< stpSize; ii++)
    //{ (stpAstNodeVec_[ii])->setupBreakPoints(); }

    //int compSize = compAstNodeVec_.size();
    //for (int ii=0;ii< compSize; ii++)
    //{ (compAstNodeVec_[ii])->setupBreakPoints(); }
  }

  return;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::getBreakPoints
// Purpose       :
// Special Notes : do not need to be sorted; other parts of Xyce will sort them
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/23/2020
//-------------------------------------------------------------------------------
bool newExpression::getBreakPoints (
    std::vector<Xyce::Util::BreakPoint> & breakPointTimes )
{
  if(isTimeDependent_)
  {
    int srcSize = srcAstNodeVec_.size();
    for (int ii=0;ii< srcSize; ii++)
    { (srcAstNodeVec_[ii])->getBreakPoints(breakPointTimes); }

    int stpSize = stpAstNodeVec_.size();
    for (int ii=0;ii< stpSize; ii++)
    { (stpAstNodeVec_[ii])->getBreakPoints(breakPointTimes); }

    int compSize = compAstNodeVec_.size();
    for (int ii=0;ii< compSize; ii++)
    { (compAstNodeVec_[ii])->getBreakPoints(breakPointTimes); }

    int limitSize = limitAstNodeVec_.size();
    for (int ii=0;ii< limitSize; ii++)
    { (limitAstNodeVec_[ii])->getBreakPoints(breakPointTimes); }

#if 0
    {
      Xyce::dout() << "newExpression::getBreakPoints. Expression "
        << expressionString_ << "  Number of breakpoints = "
        << breakPointTimes.size() <<std::endl;
      for (int ii=0;ii<breakPointTimes.size();ii++)
      {
        Xyce::dout() << "bp["<<ii<<"] = "
          << breakPointTimes[ii].value() <<std::endl;
      }
    }
#endif
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::replaceName
// Purpose       :
// Special Notes : ERK.  This called (via N_UTL_Expression) from the
//                 N_IO_DistToolBase.C file/class
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
//                 function to handle subcircuit argument nodes.  This
//                 requires a 2-pass procedure to avoid conflicts.  A test
//                 case where this is important is the bug 1806 test.  In
//                 the circuit bug_1806_2.cir, you have the following:
//
//                 .param upgefukt=1.5
//                 .subckt foobar  a b m n o p
//                 Bfoobar a b V={upgefukt*sqrt((V(m)-V(n))**2+(v(o)-v(p))**2)}
//                 .ends
//
//                 Xomgwtf 10 0 a 0 c 0 foobar
//
//
//                 The parsed expression has the nodes m,n,o and p.
//
//                 But, these are names that are local to the subcircuit.  The
//                 names  passed into the subcircuit are different.  More than
//                 one of them (which originally had different names) become
//                 ground.  Also, the "internal" "node a" and the "external"
//                 "node a" are different.
//
//                 To avoid conflicts, the node names are initially changed
//                 with a semicolon prefix.  This is normally an invalid character,
//                 so a node named ";0" would never have come from the parser.  It
//                 can only have come from this procedure.  This is the first
//                 pass.  The second pass basically removes all the semicolons.
//
//
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 2020
//-------------------------------------------------------------------------------
bool newExpression::replaceName (
    const std::string & old_name,
    const std::string & new_name)
{
  setupVariousAstArrays ();

  bool found=false;
  {
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > >::iterator iter = voltOpMap_.find(old_name);
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > >::iterator checkNewName = voltOpMap_.find(new_name);

    if (iter != voltOpMap_.end())
    {
      std::vector<Teuchos::RCP<astNode<usedType> > > & astVec = iter->second;
      for(int ii=0;ii<astVec.size();++ii)
      {
        Teuchos::RCP<voltageOp<usedType> > voltOp = Teuchos::rcp_static_cast<voltageOp<usedType> > (astVec[ii]);
        const std::string & node = voltOp->getVoltageNode();
        if(node == old_name) { voltOp->setVoltageNode(new_name); } 
      }

      if (checkNewName != voltOpMap_.end()) // "new" name already exists, so combine the vectors
      {
        std::vector<Teuchos::RCP<astNode<usedType> > > & astVec2 = checkNewName->second;
        astVec2.insert( astVec2.end(),  astVec.begin(), astVec.end() );
      }
      else // "new" name does not exist
      {
        voltOpMap_[new_name] = astVec;
      }
      voltOpMap_.erase(old_name);

      found=true;
    }

    std::vector<std::string>::iterator voltNameIter = std::find(voltNameVec_.begin(), voltNameVec_.end(), old_name);
    if (voltNameIter != voltNameVec_.end()) { *voltNameIter = new_name; }
  }

  if(!found)
  {
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > >::iterator iter = currentOpMap_.find(old_name);
    std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > >::iterator checkNewName = currentOpMap_.find(new_name);

    if (iter != currentOpMap_.end())
    {
      std::vector<Teuchos::RCP<astNode<usedType> > > & astVec = iter->second;

      for(int ii=0;ii<astVec.size();++ii)
      {
        Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (astVec[ii]);
        currOp->setCurrentDevice(new_name);
      }

      if (checkNewName != currentOpMap_.end()) // "new" name already exists, so combine the vectors
      {
        std::vector<Teuchos::RCP<astNode<usedType> > > & astVec2 = checkNewName->second;
        astVec2.insert( astVec2.end(),  astVec.begin(), astVec.end() );
      }
      else
      {
        currentOpMap_[new_name] = astVec;
      }
      currentOpMap_.erase(old_name);
      found=true;
    }

    std::vector<std::string>::iterator currentNameIter = std::find(currentNameVec_.begin(), currentNameVec_.end(), old_name);
    if (currentNameIter != currentNameVec_.end()) { *currentNameIter = new_name; }
  }

#if 0
  dumpParseTree(Xyce::dout());
#endif
  return found;
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

//-------------------------------------------------------------------------------
// Function      : newExpression::treatAsTempAndConvert()
// Purpose       : This function is used when both of the following are true:
//
//   (1) the expression is the RHS of a device "temp" parameter
//   (2) when the units need to be "auto converted" from C to K units.
//
//   ie, when this if statement is true: if (isTempParam(par.tag()) && param.getAutoConvertTemperature())
//
// Special Notes : For the old expression library, this was accomplished by
// getting the string of the original expression and then modifying it to
// include "+ CONSTCtoK", and then using that new string to create a new expression.
//
// This approach doesn't work with the new expression library, b/c the new expression
// library doesn't handle external dependencies via string replacements.
//
// The CONSTCtoK modification to the string was happening in the device entity, long
// after the IO package was done with its work on parameters.  And, thus, long after
// parameter and function resolutions.
//
// With the new expression library, doing a string modification for CONSTCtoK could
// not include those IO-based param and func resolutions.  They were lost, because they were
// not in the string.
//
// The right way to fix this is not by re-parsing with a new string.  instead, it is better
// to just modify the AST, which is what this function does.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/17/2020
//-------------------------------------------------------------------------------
void newExpression::treatAsTempAndConvert()
{
  Teuchos::RCP<astNode<usedType> > CtoK = getCtoKNode ();
  Teuchos::RCP<astNode<usedType> > newTopPtr = Teuchos::RCP<astNode<usedType> >(new binaryAddOp<usedType>  (astNodePtr_, CtoK));
  setAstPtr(newTopPtr);
}

//-------------------------------------------------------------------------------
// Function      : newExpression::setTemperatuure
// Purpose       :
// Special Notes : this is only called to overide the getTemeprature call to the
//                 group.  This will not happen very much; only when a device
//                 model has an internal self-heating model, and thus has a
//                 local temperature.  As of this writing, the only device that
//                 I know of that needs this is the thermal resistor.
//
//                 ERK.  Correction, it isn't quite true.    Currently (7/25/2020)
//                 this function gets called thru at least two pathways.
//
//                 One is via the DeviceMgr::updateTemperature function.   The
//                 other is via the DeviceEntity::updateDependentParameters
//                 function, but only the one that has "temp" as a function argument.
//                 That special version of "DeviceEntity::updateDependentParameters"
//                 seems to only be called from the thermal resistor.
//
//                 But the "updateTemperature" function is called from other places,
//                 probably whenever .STEP is invoked on TEMP.  I think.  I haven't
//                 tracked this down yet.
//
//                 Both of these pathways are used by the
//                 "MULTIPLICITY/thermal_resistor.cir" test case, which uses
//                 a thermal resistor, and also performs a .STEP over temperature.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/19/2020
//-------------------------------------------------------------------------------
bool newExpression::setTemperature (const double & temp)
{
  bool changed = false; // the value of changed is now irrelevant
  overrideGroupTemperature_ = true;
  overrideTemp_ = temp-CONSTCtoK;
  return changed;
}

//-------------------------------------------------------------------------------
// Function      : newExpression::processSuccessfulTimeStep
//
// Purpose       : Tells relevant AST nodes to update their time arrays, etc.
//                 Some AST nodes, such as DDT and SDT (time derivative and
//                 time integral, respectively) maintain arrays of data, where the
//                 array index refers to time.  These arrays have to be rotated
//                 when the time step advances, but it is difficult to know inside
//                 of an AST node when this advance should happen.  This function
//                 call is the method to let the SDT and DDT operators know that
//                 they need to update.
//
// Special Notes : I have been debating exactly how to handle this issue, and I
//                 am not yet sure if this is the best way.  For now, this function
//                 is an experiment.
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/7/2020
//-------------------------------------------------------------------------------
void newExpression::processSuccessfulTimeStep ()
{
  // this is a kludge, to test an idea.
  // Two things that I don't like about it, but were convenient for now:
  // (1) using static data in the ast classes.
  // (2) relying on "evaluateFunction" to process succesful time steps.
  //     (i.e. rotate appropriate state data for SDT and DDT, but just once)
  //
  // The nice thing about both of these is that they allowed me to not
  // have to hack every AST class.  The evaluateFunction call will traverse
  // the AST tree via the val() call.  So, I didn't have to add another
  // set of traversal functions.  I just had to hack in some conditional
  // code into the sdtOp and ddtOp val() functions.  The drawback is that
  // it involves unncessary computational cost, etc.
  //
  // The static data is for similar reasons.  I didn't feel like
  // doing the work (and all the typing) to do the same thing in a
  // non-static way.
  //
  // Here is what I would like to do, but haven't done yet:
  // (1) set up a new traversal function that will create a container of AST
  // nodes, which only contains nodes requiring state management.
  // (2) then, after calling that new traversal function once, to set up
  // that container, use that container to do the state management required
  // at the end of each successful timem step.
  //
  if( !(sdtOpVec_.empty())  ||  !(ddtOpVec_.empty()) )
  {
    staticsContainer::processSuccessfulStepFlag = true;
    usedType result;
    evaluateFunction (result);
    staticsContainer::processSuccessfulStepFlag = false;
  }
}

//-----------------------------------------------------------------------------
//
// ERK.  Note: NONE of what follows here will work if ddt is
// inside of a .func that is called more than once!!!!
//
// fix later.  This needs better book-keeping.
//
// For ddts that are inside of .funcs, there are the following complications.
// (1) arg to ddt is different once the passed values are sub'd in.
// (2) the value of ddt passed in by setDdtVals has to happen at the correct time.
//
// Another way of thinking about it;  the ddt evaluation depends on a history
// (or state) and that history goes with the *call* to ddt, not the allocation
// of ddt.
//
// processing of sdt will have the same problem.
//
//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::getDdtVals
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void newExpression::getDdtVals   (std::vector<double> & vals)
{
  std::vector<std::complex<double> > cmplxVals;
  getDdtVals(cmplxVals);
  vals.clear();
  vals.resize(cmplxVals.size());
  for (int ii=0;ii<cmplxVals.size();ii++) { vals[ii] = std::real(cmplxVals[ii]); }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::getDdtVals
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void newExpression::getDdtVals   (std::vector<std::complex<double> > & vals)
{
  vals.clear();
  vals.resize(ddtOpVec_.size());
  for (int ii=0;ii<ddtOpVec_.size();ii++)
  {
    Teuchos::RCP<ddtOp<usedType> > ddt = Teuchos::rcp_static_cast<ddtOp<usedType> > (ddtOpVec_[ii]);
    vals[ii] = ddt->getDdtArg();
  }
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::setDdtDerivs
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void newExpression::setDdtDerivs (std::vector<double> & vals)
{
  for (int ii=0;ii<ddtOpVec_.size();ii++)
  {
    Teuchos::RCP<ddtOp<usedType> > ddt = Teuchos::rcp_static_cast<ddtOp<usedType> > (ddtOpVec_[ii]);
    usedType tmp(vals[ii]);
    ddt->setDdtDeriv(tmp);
  }
}

//-----------------------------------------------------------------------------
// Function      : setDdt
// Purpose       : Helper function for setDdtDerivs
// Special Notes : double version
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void setDdt(Teuchos::RCP<ddtOp<double> > & ddt, std::complex<double> & tmpVal)
{
  ddt->setDdtDeriv(std::real(tmpVal));
}

//-----------------------------------------------------------------------------
// Function      : setDdt
// Purpose       : Helper function for setDdtDerivs
// Special Notes : std::complex<double> version
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void setDdt(Teuchos::RCP<ddtOp<std::complex<double> > > & ddt, std::complex<double> & tmpVal)
{
  ddt->setDdtDeriv(tmpVal);
}

//-----------------------------------------------------------------------------
// Function      : ExpressionInternals::setDdtDerivs
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void newExpression::setDdtDerivs (std::vector<std::complex<double> > & vals)
{
  for (int ii=0;ii<ddtOpVec_.size();ii++)
  {
    Teuchos::RCP<ddtOp<usedType> > ddt = Teuchos::rcp_static_cast<ddtOp<usedType> > (ddtOpVec_[ii]);
    setDdt(ddt,vals[ii]);
  }
}

}
}

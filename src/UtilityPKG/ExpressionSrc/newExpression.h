

#ifndef newExpression_H
#define newExpression_H

// std includes
#include<string>
#include<complex>
#include<utility>
#include<vector>
#include<list>

// trilinos includes
#include <Teuchos_RCP.hpp>

// Xyce includes.
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Interface_Enum_Types.h>

// new code includes:
#include <ast.h>
#include <ExpressionType.h>
#include <expressionGroup.h>

namespace Xyce {
namespace Util {


#if 0
static struct constant {
    const char *name;
    double value;
} constants[] = {
    { "EXP",  exp(1.0) },
    { "PI", M_PI }
} ;
static std::string specSigs[] = {"TIME", "TEMP", "VT", "FREQ"};

bool ExpressionInternals::set_temp(double const & tempIn)
{
  set_var(std::string("VT"), tempIn*CONSTKoverQ);

  return set_var(std::string("TEMP"), tempIn-CONSTCtoK);
}

#endif

//-----------------------------------------------------------------------------
// Class         : newExpression
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 10/28/2019
//-----------------------------------------------------------------------------
class newExpression
{
public:
  newExpression () :
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    expressionResolved_(false),
    expressionFunctionsResolved_(false),
    expressionParametersResolved_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    bpTol_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    numDerivs_(0),
    traditionalParse_(true),
    externalDependencies_(false),
    isTimeDepdendent_(false),
    isTempDepdendent_(false),
    isVTDepdendent_(false),
    isFreqDepdendent_(false)
  {};

  // primary constructor
  newExpression ( std::string const & exp, Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_(exp),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    expressionResolved_(false),
    expressionFunctionsResolved_(false),
    expressionParametersResolved_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    bpTol_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    numDerivs_(0),
    traditionalParse_(true),
    externalDependencies_(false),
    isTimeDepdendent_(false),
    isTempDepdendent_(false),
    isVTDepdendent_(false),
    isFreqDepdendent_(false)
  {
    // The bison file is officially case-insensitive.  So converting the
    // input string to all upper case is not necessary for it to work.
    //
    // However:
    //
    // Xyce mostly deals with netlist strings by converting to upper case.
    // So, when the code outside of bison code interacts with it, Bison will not
    // have converted it to upper or lower case.  It will be the original case,
    // whatever that was in the netlist.
    //
    // The simplest, easiest way to make all of this work is to simply
    // upcase the whole string.
    //
    Xyce::Util::toUpper(expressionString_);

    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
  };


  // special "big table" constructor
  newExpression (const std::vector<usedType> & xvals, const std::vector<usedType> & yvals,
      Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_("TIME"),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    expressionResolved_(false),
    expressionFunctionsResolved_(false),
    expressionParametersResolved_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    bpTol_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    numDerivs_(0),
    traditionalParse_(false),
    externalDependencies_(false),
    isTimeDepdendent_(false),
    isTempDepdendent_(false),
    isVTDepdendent_(false),
    isFreqDepdendent_(false)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());

    Teuchos::RCP<astNode<usedType> > time_base = timeNodePtr_;
    tableNodePtrPtr_ = new Teuchos::RCP<tableOp<usedType> >(new tableOp<usedType> (time_base, xvals, yvals));
    astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*tableNodePtrPtr_);
    setAstPtr(*astNodePtrPtr_);
  };

  // another "big table" constructor
  newExpression (Teuchos::RCP<astNode<usedType> > & left,
      const std::vector<usedType> & xvals, const std::vector<usedType> & yvals,
      Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_("TIME"),
    parsed_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    expressionResolved_(false),
    expressionFunctionsResolved_(false),
    expressionParametersResolved_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    bpTol_(0.0),
    timeStep_(0.0),
    timeStepAlpha_(0.0),
    timeStepPrefac_(0.0),
    numDerivs_(0),
    traditionalParse_(false),
    externalDependencies_(false),
    isTimeDepdendent_(false),
    isTempDepdendent_(false),
    isVTDepdendent_(false),
    isFreqDepdendent_(false)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("GARBAGE")));

    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TIME")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("TEMP")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("VT")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("FREQ")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());

    tableNodePtrPtr_ = new Teuchos::RCP<tableOp<usedType> >(new tableOp<usedType> (left, xvals, yvals));
    astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*tableNodePtrPtr_);
    setAstPtr(*astNodePtrPtr_);

    if( !(Teuchos::is_null(left)) )
    {
      if (left->paramType())   { paramOpVec_.push_back(left); }
      if (left->funcType())    { funcOpVec_.push_back(left); }
      if (left->voltageType()) { voltOpVec_.push_back(left); }
      if (left->currentType()) { currentOpVec_.push_back(left); }
      left->getInterestingOps(paramOpVec_,funcOpVec_, voltOpVec_,currentOpVec_);

#if 0
      paramNameVec_.clear();
      std::vector<Teuchos::RCP<astNode<usedType> > >::iterator iterParamOp;
      for (iterParamOp=paramOpVec_.begin();iterParamOp!=paramOpVec_.end();)
      {
        std::string name = iterParamOp->getName();
        std::vector<std::string>::iterator nameIter = std::find(paramNameVec_.begin(),paramNameVec_.end(),name);

        if (nameIter != paramNameVec_.end())
        {
          paramNameVec_.push_back(name);
          ++iterParamOp;
        }
        else
        {
          iterParamOp = paramOpVec_.erase(iterParamOp);
        }
      }
#endif
    }
// setup voltOpNames here?
// setup currentOpNames here?

  };

  // copy constructor - this may need work
  // This is necessary b/c of things like the masterAstNodeVec_ object, which uses raw pointers
  newExpression (const newExpression & right) :
    group_(right.group_),
    expressionString_(right.expressionString_),
    parsed_(right.parsed_),
    derivsSetup_(right.derivsSetup_),
    astArraysSetup_(right.astArraysSetup_),
    expressionResolved_(right.expressionResolved_),
    expressionFunctionsResolved_(right.expressionFunctionsResolved_),
    expressionParametersResolved_(right.expressionParametersResolved_),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    functionArgStringVec_(right.functionArgStringVec_),
    functionArgOpVec_ (right.functionArgOpVec_),
    paramOpVec_(right.paramOpVec_),
    unresolvedParamOpVec_(right.unresolvedParamOpVec_),
    funcOpVec_(right.funcOpVec_),
    unresolvedFuncOpVec_(right.unresolvedFuncOpVec_),
    voltOpVec_(right.voltOpVec_),
    unresolvedVoltOpVec_(right.unresolvedVoltOpVec_),
    voltOpNames_(right.voltOpNames_),
    currentOpVec_(right.currentOpVec_),
    unresolvedCurrentOpVec_(right.unresolvedCurrentOpVec_),
    currentOpNames_(right.currentOpNames_),
    bpTol_(right.bpTol_),
    timeStep_(right.timeStep_),
    timeStepAlpha_(right.timeStepAlpha_),
    timeStepPrefac_(right.timeStepPrefac_),
    derivIndexVec_ (right.derivIndexVec_),
    derivNodeIndexMap_(right.derivNodeIndexMap_),
    numDerivs_(right.numDerivs_),
    traditionalParse_(right.traditionalParse_),
    externalDependencies_(right.externalDependencies_),
    isTimeDepdendent_(right.isTimeDepdendent_),
    isTempDepdendent_(right.isTempDepdendent_),
    isVTDepdendent_(right.isVTDepdendent_),
    isFreqDepdendent_(right.isFreqDepdendent_)
  {
    garbageParamOpPtr_ = right.garbageParamOpPtr_;
    timeNodePtr_ = right.timeNodePtr_;
    tempNodePtr_ = right.tempNodePtr_;
    vtNodePtr_   = right.vtNodePtr_;
    freqNodePtr_ = right.freqNodePtr_;
    piNodePtr_   = right.piNodePtr_;
    astNodePtr_ = right.astNodePtr_; // copy over the whole tree

    // pointers to RCP's. use constructors
    if (right.astNodePtrPtr_)
    {
      astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*(right.astNodePtrPtr_));
    }
    if (right.tableNodePtrPtr_)
    {
      tableNodePtrPtr_ = new Teuchos::RCP< tableOp<usedType> >(*(right.tableNodePtrPtr_));
    }

    // vectors of pointers to RCPs of AST nodes; use constructors
    for (int ii=0;ii<right.masterAstNodeVec_.size();ii++)
    {
      masterAstNodeVec_.push_back( new Teuchos::RCP<astNode<usedType> >(*(right.masterAstNodeVec_[ii])) );
    }

    for (int ii=0;ii<right.srcAstNodeVec_.size();ii++)
    {
      srcAstNodeVec_.push_back( new Teuchos::RCP<astNode<usedType> >(*(right.srcAstNodeVec_[ii])) );
    }
  };

  // assignment operator
  // This is necessary b/c of things like the masterAstNodeVec_ object, which uses raw pointers
 	newExpression & operator =(const newExpression & right)
  {
    group_ = right.group_;
    expressionString_ = right.expressionString_;
    parsed_ = right.parsed_;
    derivsSetup_ = right.derivsSetup_;
    astArraysSetup_ = right.astArraysSetup_;
    expressionResolved_ = right.expressionResolved_;
    expressionFunctionsResolved_ = right.expressionFunctionsResolved_;
    expressionParametersResolved_ = right.expressionParametersResolved_;
    astNodePtrPtr_ = NULL;
    tableNodePtrPtr_ = NULL;
    functionArgStringVec_ = right.functionArgStringVec_;
    functionArgOpVec_  = right.functionArgOpVec_;
    paramOpVec_ = right.paramOpVec_;
    unresolvedParamOpVec_ = right.unresolvedParamOpVec_;
    funcOpVec_ = right.funcOpVec_;
    unresolvedFuncOpVec_ = right.unresolvedFuncOpVec_;
    voltOpVec_ = right.voltOpVec_;
    voltOpNames_ = right.voltOpNames_;
    unresolvedVoltOpVec_ = right.unresolvedVoltOpVec_;
    currentOpVec_ = right.currentOpVec_;
    currentOpNames_ = right.currentOpNames_;
    unresolvedCurrentOpVec_ = right.unresolvedCurrentOpVec_;

    bpTol_ = right.bpTol_;
    timeStep_ = right.timeStep_;
    timeStepAlpha_ = right.timeStepAlpha_;
    timeStepPrefac_ = right.timeStepPrefac_;

    derivIndexVec_ = right.derivIndexVec_;
    derivNodeIndexMap_ = right.derivNodeIndexMap_;
    numDerivs_ = right.numDerivs_;
    traditionalParse_ = right.traditionalParse_;
    externalDependencies_ = right.externalDependencies_;

    isTimeDepdendent_ = right.isTimeDepdendent_;
    isTempDepdendent_ = right.isTempDepdendent_;
    isVTDepdendent_ = right.isVTDepdendent_;
    isFreqDepdendent_ = right.isFreqDepdendent_;

    garbageParamOpPtr_ = right.garbageParamOpPtr_;
    timeNodePtr_ = right.timeNodePtr_;
    tempNodePtr_ = right.tempNodePtr_;
    vtNodePtr_   = right.vtNodePtr_;
    freqNodePtr_ = right.freqNodePtr_;
    piNodePtr_   = right.piNodePtr_;
    astNodePtr_ = right.astNodePtr_; // copy over the whole tree

    // pointers to RCP's. use constructors
    if (right.astNodePtrPtr_)
    {
      astNodePtrPtr_ = new Teuchos::RCP<astNode<usedType> >(*(right.astNodePtrPtr_));
    }
    if (right.tableNodePtrPtr_)
    {
      tableNodePtrPtr_ = new Teuchos::RCP< tableOp<usedType> >(*(right.tableNodePtrPtr_));
    }

    // vectors of pointers to RCPs of AST nodes; use constructors
    for (int ii=0;ii<right.masterAstNodeVec_.size();ii++)
    {
      masterAstNodeVec_.push_back( new Teuchos::RCP<astNode<usedType> >(*(right.masterAstNodeVec_[ii])) );
    }

    for (int ii=0;ii<right.srcAstNodeVec_.size();ii++)
    {
      srcAstNodeVec_.push_back( new Teuchos::RCP<astNode<usedType> >(*(right.srcAstNodeVec_[ii])) );
    }

    return *this;
  };

  // destructor
  ~newExpression ()
  {
    // pointers to RCP's. use constructors
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
  };

  bool lexAndParseExpression();

  //bool resolveExpression();
  bool attachFunctionNode(const std::string & funcName, Teuchos::RCP<Xyce::Util::newExpression> expPtr);
  bool attachParameterNode(const std::string & paramName, Teuchos::RCP<Xyce::Util::newExpression> expPtr);

  void clear(); // reset expression to the state it should be before lexAndParseExpression

  bool parsed() const { return parsed_; };
  bool derivsSetup () const { return derivsSetup_; };
  bool astArraysSetup () const { return astArraysSetup_; }
  bool expressionResolved() const {return expressionResolved_; }

  void setExpressionString (std::string const & exp)  { expressionString_ = exp; }

  bool make_constant (std::string const & var, usedType const & val);
  bool make_var (std::string const & var);

  void getGlobalParamNames ( std::vector<std::string> & names )
  {
    for (int ii=0;ii<paramOpVec_.size();++ii)
    {
      Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (paramOpVec_[ii]);
      if ( parOp->getVar() )
      {
        names.push_back( parOp->getName() );
      }
    }
  };

  // added this accessor to allow the possibility of external resolution, rather than via the group.
  void setExpressionResolved (bool res) { expressionResolved_ = res; }

  void setAstPtr(Teuchos::RCP<astNode<usedType> > & astNodePtr) { astNodePtr_ = astNodePtr; };

  // these two functions return int error codes in the original expression library
  int evaluate (usedType &result, std::vector< usedType > &derivs);
  int evaluateFunction (usedType &result);

  void dumpParseTree(std::ostream & os) { if ( !(Teuchos::is_null(astNodePtr_)) ){astNodePtr_->output(os); }}

  bool getBreakPoints (std::vector<Xyce::Util::BreakPoint> & breakPointTimes )
  {
    for (int ii=0;ii< srcAstNodeVec_.size(); ii++) { (*(srcAstNodeVec_[ii]))->getBreakPoints(breakPointTimes); }
    return true;
  }

  Teuchos::RCP<astNode<usedType> > & getAst() {return astNodePtr_;}

  Teuchos::RCP<paramOp<usedType> > & getGarbageParam() {return garbageParamOpPtr_;}

  Teuchos::RCP<astNode<usedType> > getTimeNode () { return timeNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getTempNode () { return tempNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getVtNode () { return vtNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getFreqNode () { return freqNodePtr_; }
  Teuchos::RCP<astNode<usedType> > getPiNode () { return piNodePtr_; }

  // some of the parameter and function objects are stored in multiple containers.
  void setFunctionArgStringVec (const std::vector<std::string> & args);

  std::vector<std::string> & getFunctionArgStringVec () { return functionArgStringVec_; };

  std::vector< Teuchos::RCP<astNode<usedType> > > & getFunctionArgOpVec() { return functionArgOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getParamOpVec () { return paramOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedParamOpVector() {  return unresolvedParamOpVec_; };
  //std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getParamOpMap () { return paramOpMap_; };
  std::vector<std::string> & getParamNameVec () { return paramNameVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getFuncOpVec () { return funcOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedFuncOpVec() { return unresolvedFuncOpVec_; };
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getFuncOpMap () { return funcOpMap_; };
  std::vector< std::string > & getFuncNameVec () { return funcNameVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getVoltOpVec () { return voltOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedVoltOpVec() { return unresolvedVoltOpVec_; };
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getVoltOpNames ()
  {
    //if (!expressionResolved_) { resolveExpression(); }
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return voltOpNames_;
  };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getCurrentOpVec () { return currentOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedCurrentOpVec() { return unresolvedCurrentOpVec_; };
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > & getCurrentOpNames ()
  {
    //if (!expressionResolved_) { resolveExpression(); }
    if (!astArraysSetup_) { setupVariousAstArrays_ (); }
    return currentOpNames_;
  };

  void codeGen( std::ostream & os ) { astNodePtr_->codeGen(os); os << ";" <<std::endl; };

  std::vector< Teuchos::RCP<astNode<usedType> > * > & getMasterNodeVec() { return masterAstNodeVec_; }
  std::vector< Teuchos::RCP<astNode<usedType> > * > & getSrcNodeVec() { return srcAstNodeVec_;}

  const std::string & getExpressionString() { return expressionString_; };

  double getTime() { return std::real(timeNodePtr_->val()); };

  bool getTimeDependent() { return isTimeDepdendent_; }
  void setTimeDependent(bool val) { isTimeDepdendent_ = val; }

  bool getTempDependent() { return isTempDepdendent_; }
  void setTempDependent(bool val) { isTempDepdendent_ = val; }

  bool getVTDependent() { return isVTDepdendent_; }
  void setVTDependent(bool val) { isVTDepdendent_ = val; }

  bool getFreqDependent() { return isFreqDepdendent_; }
  void setFreqDependent(bool val) { isFreqDepdendent_ = val; }

  // note: I don't particularly like these next 2 functions, but they are needed
  // to support the old expression API.
  //
  // this function is only used to determine function arguments of a function prototype
  // So if we have .func abc(x,y) {x+y+10}
  // At a certain point the prototype abc(x,y) will get parsed, and x,y are the prototype args.
  //
  // The old Xyce source code uses an expression to parse abc(x,y) so I'm attempting to 
  // do the same here. When doing this, it is necessary to keep the args in order, 
  // which is why I can't rely on the getParamOpVec function (as the paramOpVec is in 
  // the same order as the paramOpMap)
  void getFuncPrototypeArgStrings ( std::vector< std::string > & funcArgStrings )
  {
    //if (!expressionResolved_) { resolveExpression(); }

    funcArgStrings.clear();
    if(!(funcOpVec_.empty()))
    {
      Teuchos::RCP<funcOp<usedType> > funcOperator = Teuchos::rcp_dynamic_cast<funcOp<usedType> > (funcOpVec_[0]) ;
      if ( !(Teuchos::is_null(funcOperator)) )
      {
        std::vector< Teuchos::RCP<astNode<usedType> > > & funcArgs = funcOperator->getFuncArgs();
        for(int ii=0;ii<funcArgs.size();++ii) { funcArgStrings.push_back(funcArgs[ii]->getName()); }
      }
    }
    return;
  }

  // when parsing the function prototype, the function name is needed as well.
  // Note; this was needed by the OLD API.  May not be needed now.
  void getFuncPrototypeName ( std::string & prototypeName) 
  {
    //if (!expressionResolved_) { resolveExpression(); }

    if(!(funcOpVec_.empty()))
    {
      prototypeName = funcOpVec_[0]->getName();
    }
  }

private:
  void setupDerivatives_ ();
  void setupVariousAstArrays_ () ;

  Teuchos::RCP<baseExpressionGroup> group_;
  std::string expressionString_;
  bool parsed_;
  bool derivsSetup_;
  bool astArraysSetup_;
  bool expressionResolved_;
  bool expressionFunctionsResolved_;
  bool expressionParametersResolved_;

  Teuchos::RCP<astNode<usedType> > astNodePtr_;
  Teuchos::RCP<astNode<usedType> > * astNodePtrPtr_;
  Teuchos::RCP<tableOp<usedType> > * tableNodePtrPtr_;

  // function argument vectors.  Both strings and operators
  std::vector< std::string > functionArgStringVec_;

  Teuchos::RCP<paramOp<usedType> > garbageParamOpPtr_; // convenient for unused params/args

  std::vector< Teuchos::RCP<astNode<usedType> > > functionArgOpVec_;

  std::vector<std::string> paramNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > paramOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedParamOpVec_;
  //std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > paramOpMap_; 

  std::vector<std::string> funcNameVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > funcOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedFuncOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > funcOpMap_; 

  std::vector<Teuchos::RCP<astNode<usedType> > > voltOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedVoltOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > voltOpNames_;

  std::vector<Teuchos::RCP<astNode<usedType> > > currentOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedCurrentOpVec_;
  std::unordered_map<std::string,std::vector<Teuchos::RCP<astNode<usedType> > > > currentOpNames_;

  // master vector of nodes.  This is only used for deleting the ast tree in
  // the destructor.  The tree should be deleted by marching down the
  // branches of the tree, as some of the nodes use the same pointer.
  std::vector< Teuchos::RCP<astNode<usedType> > * > masterAstNodeVec_;

  // time integration related variables
  double bpTol_;
  double timeStep_;
  double timeStepAlpha_;
  double timeStepPrefac_;

  // vector of independent sources, but only those with breakpoints.  This
  // vector is ONLY used for obtaining breakpoints.
  std::vector< Teuchos::RCP<astNode<usedType> > * > srcAstNodeVec_;

  // const and specials nodes:
  Teuchos::RCP<specialsOp<usedType> > timeNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > tempNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > vtNodePtr_;
  Teuchos::RCP<specialsOp<usedType> > freqNodePtr_;
  Teuchos::RCP<piConstOp<usedType> > piNodePtr_;

  // derivative related stuff
  typedef  std::pair< Teuchos::RCP<astNode<usedType> > , int > derivIndexPair_;
  std::vector< derivIndexPair_  > derivIndexVec_;

  std::unordered_map<std::string,int> derivNodeIndexMap_;
  int numDerivs_;

  bool traditionalParse_;
  bool externalDependencies_; // true if expression includes a call to a .func, .param or .global_param

  bool isTimeDepdendent_;
  bool isTempDepdendent_;
  bool isVTDepdendent_;
  bool isFreqDepdendent_;
};

}
}

#endif



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
  newExpression () {};

  // primary constructor
  newExpression ( std::string const & exp, Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_(exp),
    parsed_(false),
    resolved_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    voltOpVecSetup_(false),
    currentOpVecSetup_(false),
    numDerivs_(0),
    traditionalParse_(true),
    externalDependencies_(false)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("garbage")));

    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("time")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("temp")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("vt")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("freq")));
    piNodePtr_   = Teuchos::rcp(new piConstOp<usedType>  ());
  };


  // special "big table" constructor
  newExpression (const std::vector<usedType> & xvals, const std::vector<usedType> & yvals,
      Teuchos::RCP<baseExpressionGroup> & group ) :
    group_(group),
    expressionString_("time"),
    parsed_(false),
    resolved_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    voltOpVecSetup_(false),
    currentOpVecSetup_(false),
    numDerivs_(0),
    traditionalParse_(false),
    externalDependencies_(false)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("garbage")));

    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("time")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("temp")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("vt")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("freq")));
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
    expressionString_("time"),
    parsed_(false),
    resolved_(false),
    derivsSetup_(false),
    astArraysSetup_(false),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    voltOpVecSetup_(false),
    currentOpVecSetup_(false),
    numDerivs_(0),
    traditionalParse_(false),
    externalDependencies_(false)
  {
    garbageParamOpPtr_ = Teuchos::rcp(new paramOp<usedType> (std::string("garbage")));

    timeNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("time")));
    tempNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("temp")));
    vtNodePtr_   = Teuchos::rcp(new specialsOp<usedType> (std::string("vt")));
    freqNodePtr_ = Teuchos::rcp(new specialsOp<usedType> (std::string("freq")));
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
    }
  };

  // copy constructor - this may need work
  // This is necessary b/c of things like the masterAstNodeVec_ object, which uses raw pointers
  newExpression (const newExpression & right) :
    group_(right.group_),
    expressionString_(right.expressionString_),
    parsed_(right.parsed_),
    resolved_(right.resolved_),
    derivsSetup_(right.derivsSetup_),
    astArraysSetup_(right.astArraysSetup_),
    astNodePtrPtr_(NULL),
    tableNodePtrPtr_(NULL),
    functionArgStringVec_(right.functionArgStringVec_),
    paramOpMap_(right.paramOpMap_),
    functionArgOpVec_ (right.functionArgOpVec_),
    paramOpVec_(right.paramOpVec_),
    unresolvedParamOpVec_(right.unresolvedParamOpVec_),
    funcOpVec_(right.funcOpVec_),
    unresolvedFuncOpVec_(right.unresolvedFuncOpVec_),
    voltOpVecSetup_(right.voltOpVecSetup_),
    voltOpVec_(right.voltOpVec_),
    unresolvedVoltOpVec_(right.unresolvedVoltOpVec_),
    currentOpVecSetup_(right.currentOpVecSetup_),
    currentOpVec_(right.currentOpVec_),
    unresolvedCurrentOpVec_(right.unresolvedCurrentOpVec_),
    derivIndexVec_ (right.derivIndexVec_),
    derivNodeIndexMap_(right.derivNodeIndexMap_),
    numDerivs_(right.numDerivs_),
    traditionalParse_(right.traditionalParse_),
    externalDependencies_(right.externalDependencies_)
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
    resolved_ = right.resolved_;
    derivsSetup_ = right.derivsSetup_;
    astArraysSetup_ = right.astArraysSetup_;
    astNodePtrPtr_ = NULL;
    tableNodePtrPtr_ = NULL;
    functionArgStringVec_ = right.functionArgStringVec_;
    paramOpMap_ = right.paramOpMap_;
    functionArgOpVec_  = right.functionArgOpVec_;
    paramOpVec_ = right.paramOpVec_;
    unresolvedParamOpVec_ = right.unresolvedParamOpVec_;
    funcOpVec_ = right.funcOpVec_;
    unresolvedFuncOpVec_ = right.unresolvedFuncOpVec_;
    voltOpVecSetup_ = right.voltOpVecSetup_;
    voltOpVec_ = right.voltOpVec_;
    unresolvedVoltOpVec_ = right.unresolvedVoltOpVec_;
    currentOpVecSetup_ = right.currentOpVecSetup_;
    currentOpVec_ = right.currentOpVec_;
    unresolvedCurrentOpVec_ = right.unresolvedCurrentOpVec_;
    derivIndexVec_ = right.derivIndexVec_;
    derivNodeIndexMap_ = right.derivNodeIndexMap_;
    numDerivs_ = right.numDerivs_;
    traditionalParse_ = right.traditionalParse_;
    externalDependencies_ = right.externalDependencies_;

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
    masterAstNodeVec_.clear();
  };

  bool lexAndParseExpression();

  void resolveExpression();

  bool parsed() const { return parsed_; };
  bool derivsSetup () const { return derivsSetup_; };
  bool astArraysSetup () const { return astArraysSetup_; }

  bool set(std::string const & exp); // is this needed?  

  //void getSymbolTable (std::vector< ExpressionSymbolTableEntry > & names) const;

#if 0
  void get_names (int const & type, std::vector< std::string > & names) const;
#else

#endif

  int get_type (std::string const & var);
  bool make_constant (std::string const & var, double const & val);
  bool make_var (std::string const & var);


  void setAstPtr(Teuchos::RCP<astNode<usedType> > & astNodePtr) { astNodePtr_ = astNodePtr; };

  // these two functions return int error codes in the original expression library
  int evaluate (usedType &result, std::vector< usedType > &derivs);
  int evaluateFunction (usedType &result);

  void dumpParseTree(std::ostream & os) { astNodePtr_->output(os); }

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
  void setFunctionArgStringVec (const std::vector<std::string> & args) { functionArgStringVec_ = args; };
  std::vector<std::string> & getFunctionArgStringVec () { return functionArgStringVec_; };

  std::unordered_map<std::string,Teuchos::RCP<astNode<usedType> > > & getParamOpMap () { return paramOpMap_; };

  std::vector< Teuchos::RCP<astNode<usedType> > > & getFunctionArgOpVec() { return functionArgOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getParamOpVec () { return paramOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedParamOpVector() {  return unresolvedParamOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getFuncOpVec () { return funcOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedFuncOpVec() { return unresolvedFuncOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getVoltOpVec () { return voltOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedVoltOpVec() { return unresolvedVoltOpVec_; };

  std::vector<Teuchos::RCP<astNode<usedType> > > & getCurrentOpVec () { return currentOpVec_; };
  std::vector<Teuchos::RCP<astNode<usedType> > > & getUnresolvedCurrentOpVec() { return unresolvedCurrentOpVec_; };

  void codeGen( std::ostream & os ) { astNodePtr_->codeGen(os); os << ";" <<std::endl; };

  std::vector< Teuchos::RCP<astNode<usedType> > * > & getMasterNodeVec() { return masterAstNodeVec_; }
  std::vector< Teuchos::RCP<astNode<usedType> > * > & getSrcNodeVec() { return srcAstNodeVec_;}

  const std::string & getExpressionString() { return expressionString_; };

private:
  void setupDerivatives_ ();
  void setupVariousAstArrays_ () ;

  Teuchos::RCP<baseExpressionGroup> group_;
  std::string expressionString_;
  bool parsed_;
  bool resolved_;
  bool derivsSetup_;
  bool astArraysSetup_;
  Teuchos::RCP<astNode<usedType> > astNodePtr_;
  Teuchos::RCP<astNode<usedType> > * astNodePtrPtr_;
  Teuchos::RCP<tableOp<usedType> > * tableNodePtrPtr_;

  // function argument vectors.  Both strings and operators
  std::vector< std::string > functionArgStringVec_;

  Teuchos::RCP<paramOp<usedType> > garbageParamOpPtr_; // convenient for unused params/args

  std::unordered_map<std::string,Teuchos::RCP<astNode<usedType> > > paramOpMap_;
  std::vector< Teuchos::RCP<astNode<usedType> > > functionArgOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > paramOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedParamOpVec_;

  std::vector<Teuchos::RCP<astNode<usedType> > > funcOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedFuncOpVec_;

  bool voltOpVecSetup_;
  std::vector<Teuchos::RCP<astNode<usedType> > > voltOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedVoltOpVec_;

  bool currentOpVecSetup_;
  std::vector<Teuchos::RCP<astNode<usedType> > > currentOpVec_;
  std::vector<Teuchos::RCP<astNode<usedType> > > unresolvedCurrentOpVec_;

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
};

}
}

#endif

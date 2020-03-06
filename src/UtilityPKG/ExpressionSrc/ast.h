
#ifndef ast_H
#define ast_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <cmath>
#include <complex>

#include <Teuchos_RCP.hpp>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_Interpolators.h>
#include <N_UTL_ExtendedString.h>

template <typename ScalarT>
class astNode;

template <typename ScalarT>
class funcOp;

template <typename ScalarT>
class paramOp;

template <typename ScalarT>
class voltageOp;

template <typename ScalarT>
class currentOp;

inline void yyerror(std::vector<std::string> & s);

//-------------------------------------------------------------------------------
// base node class
template <typename ScalarT>
class astNode
{
  public:
    astNode() {};
    virtual ~astNode() {};

    astNode( Teuchos::RCP<astNode<ScalarT> > &left ):
      leftAst_(left) {};

    astNode(Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      leftAst_(left),rightAst_(right) {};

    virtual ScalarT val() = 0;
    virtual ScalarT dx(int i) = 0;
    virtual void output(std::ostream & os, int indent=0) = 0;
    virtual void codeGen (std::ostream & os ) = 0;

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) {};
    virtual void unsetNode() {};

    virtual void setFuncArgs(std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpArgVec ) {};

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes) { return true;}

    virtual void setDerivIndex(int i) {};
    virtual void unsetDerivIndex() {};

    virtual void setValue(ScalarT val) {}; // supports specialsOp, and paramOp. otherwise no-op
    virtual void unsetValue() {};          // supports specialsOp, and paramOp. otherwise no-op

    // base class no-ops.  Derived functions only in paramOp, only called from ddx.
    virtual void setVar() {};
    virtual void unsetVar() {};
    virtual bool getVar() { return false; }

    virtual bool numvalType()      { return false; };
    virtual bool paramType()       { return false; };
    virtual bool funcType()        { return false; };
    virtual bool voltageType()     { return false; };
    virtual bool currentType()     { return false; };

    virtual bool getFunctionArgType() { return false; };
    virtual void setFunctionArgType() {};
    virtual void unsetFunctionArgType() {};

    virtual std::string getName () { return std::string(""); };
    virtual std::vector<std::string> getNodeNames() { std::vector<std::string> tmp; return tmp; }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(leftAst_)) )
      {
        if (leftAst_->paramType())   { paramOpVector.push_back(leftAst_); }
        if (leftAst_->funcType())    { funcOpVector.push_back(leftAst_); }
        if (leftAst_->voltageType()) { voltOpVector.push_back(leftAst_); }
        if (leftAst_->currentType()) { currentOpVector.push_back(leftAst_); }
        leftAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(rightAst_)) )
      {
        if (rightAst_->paramType())   { paramOpVector.push_back(rightAst_); }
        if (rightAst_->funcType())    { funcOpVector.push_back(rightAst_); }
        if (rightAst_->voltageType()) { voltOpVector.push_back(rightAst_); }
        if (rightAst_->currentType()) { currentOpVector.push_back(rightAst_); }
        rightAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(leftAst_)) )
      {
        if (leftAst_->paramType()) { paramOpVector.push_back(leftAst_); }
        leftAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(rightAst_)) )
      {
        if (rightAst_->paramType()) { paramOpVector.push_back(rightAst_); }
        rightAst_->getParamOps(paramOpVector);
      }
    }

    // func arg ops are of class paramOp, but have been identified as being
    // function arguments.  (thus being excluded from the getParamOps function, above)
    // This is needed by the ddx function
    //
    // ERK. This needs to be propagated down to other AST classes.  2/8/2020. Not done yet.
    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(leftAst_)) )
      {
        if (leftAst_->getFunctionArgType()) { funcArgOpVector.push_back(leftAst_); }
        leftAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(rightAst_)) )
      {
        if (rightAst_->getFunctionArgType()) { funcArgOpVector.push_back(rightAst_); }
        rightAst_->getFuncArgOps(funcArgOpVector);
      }
    }

    // func Ops are operators of class funcOp.
    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector) 
    {
      if( !(Teuchos::is_null(leftAst_)) )
      {
        if (leftAst_->funcType()) { funcOpVector.push_back(leftAst_); }
        leftAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(rightAst_)) )
      {
        if (rightAst_->funcType()) { funcOpVector.push_back(rightAst_); }
        rightAst_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(leftAst_)) )
      {
        if (leftAst_->voltageType()) { voltOpVector.push_back(leftAst_); }
        leftAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(rightAst_)) )
      {
        if (rightAst_->voltageType()) { voltOpVector.push_back(rightAst_); }
        rightAst_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(leftAst_)) )
      {
        if (leftAst_->currentType()) { currentOpVector.push_back(leftAst_); }
        leftAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(rightAst_)) )
      {
        if (rightAst_->currentType()) { currentOpVector.push_back(rightAst_); }
        rightAst_->getCurrentOps(currentOpVector);
      }
    }

  protected:
    Teuchos::RCP<astNode<ScalarT> > leftAst_;
    Teuchos::RCP<astNode<ScalarT> > rightAst_;
};

#include "astbinary.h"
#include "astfuncs.h"
#include "astcomp.h"
#include "ast_spice_src.h"

//-------------------------------------------------------------------------------
// pow "^"  operator
template <typename ScalarT>
class powOp : public astNode<ScalarT>
{
  public:
    powOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    };

    virtual ScalarT val() { return std::pow(this->leftAst_->val(), this->rightAst_->val());}

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;
      ScalarT retVal = 0.0;

      if (rightConst_ && !leftConst_) {if (lef->val() != 0.0) { retVal = rig->val()*lef->dx(i)/lef->val()*std::pow(lef->val(),rig->val()) ; }}
      else if (!rightConst_ && leftConst_) {if (lef->val() != 0.0) { retVal = std::log(lef->val())*std::pow(lef->val(),rig->val())*rig->dx(i); }}
      else {if (lef->val() != 0.0) { retVal = (rig->dx(i)*std::log(lef->val())+rig->val()*lef->dx(i)/lef->val())*std::pow(lef->val(),rig->val());}}
      return  retVal;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "power operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::pow(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

  private:
    bool rightConst_;
    bool leftConst_;

};

//-------------------------------------------------------------------------------
// atan2(y,x) operator
template <typename ScalarT>
class atan2Op : public astNode<ScalarT>
{
  public:
    atan2Op (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right), rightConst_(true),leftConst_(false)
    {
      rightConst_ = this->rightAst_->numvalType();
      leftConst_ = this->leftAst_->numvalType();
    };

    virtual ScalarT val() { return std::atan2(this->leftAst_->val(), this->rightAst_->val());}

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & lef = this->leftAst_;
      Teuchos::RCP<astNode<ScalarT> > & rig = this->rightAst_;
      ScalarT retVal = 0.0;

      if (rightConst_ && !leftConst_) { retVal = (rig->val()*lef->dx(i))/ (lef->val()*lef->val() + rig->val()*rig->val()); }
      else if (!rightConst_ && leftConst_) { retVal = (-lef->val()*rig->dx(i)) / (lef->val()*lef->val() + rig->val()*rig->val()); }
      else { retVal = (rig->val()*lef->dx(i) - lef->val()*rig->dx(i))/ (lef->val()*lef->val() + rig->val()*rig->val()) ; }
      return  retVal;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "atan2 operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::atan2(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

  private:
    bool rightConst_;
    bool leftConst_;

};

//-------------------------------------------------------------------------------
// phase operator
template <typename ScalarT>
class phaseOp : public astNode<ScalarT>
{
  public:
    phaseOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::arg(this->leftAst_->val()); }

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (phase) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "phase operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::arg(";
      this->leftAst_->codeGen(os);
      os << ")";
    }
};

//-------------------------------------------------------------------------------
// real part  operator
template <typename ScalarT>
class realOp : public astNode<ScalarT>
{
  public:
    realOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::real(this->leftAst_->val()); }

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (real) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "real operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::real(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// imag part  operator
template <typename ScalarT>
class imagOp : public astNode<ScalarT>
{
  public:
    imagOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::imag(this->leftAst_->val()); }

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (imag) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "imag operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::imag(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// max operator
template <typename ScalarT>
class maxOp : public astNode<ScalarT>
{
  public:
    maxOp ( Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right ): astNode<ScalarT>(left,right) {};

    virtual ScalarT val()
    { return std::max( std::real(this->leftAst_->val()), std::real(this->rightAst_->val()) ); }

    virtual ScalarT dx(int i)
    {
      bool cmp = std::real(this->leftAst_->val()) < std::real(this->rightAst_->val());
      return cmp?(std::real(this->rightAst_->dx(i))):(std::real(this->leftAst_->dx(i)));
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "max operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::max(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }
};

//-------------------------------------------------------------------------------
// min operator
template <typename ScalarT>
class minOp : public astNode<ScalarT>
{
  public:
    minOp ( Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right ): astNode<ScalarT>(left,right) {};

    virtual ScalarT val() { return std::min( std::real(this->leftAst_->val()), std::real(this->rightAst_->val()) ); }

    virtual ScalarT dx(int i)
    {
      bool cmp = std::real(this->rightAst_->val()) < std::real(this->leftAst_->val()) ;
      return (!cmp)?(std::real(this->leftAst_->dx(i))):(std::real(this->rightAst_->dx(i)));
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "min operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::min(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// unary logical NOT sign
template <typename ScalarT>
class unaryNotOp : public astNode<ScalarT>
{
  public:
    unaryNotOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      return ((std::real(this->leftAst_->val())==0)?1:0);
    }

    virtual ScalarT dx(int i) {return 0.0;}

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unary NOT operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(!";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// unary minus sign
template <typename ScalarT>
class unaryMinusOp : public astNode<ScalarT>
{
  public:
    unaryMinusOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return (-(this->leftAst_->val())); }

    virtual ScalarT dx(int i) { return (-(this->leftAst_->dx(i))); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unary minus operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(-";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// unary plus sign
template <typename ScalarT>
class unaryPlusOp : public astNode<ScalarT>
{
  public:
    unaryPlusOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return (+(this->leftAst_->val())); }

    virtual ScalarT dx(int i) { return (+(this->leftAst_->dx(i))); }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unary plus operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(+";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};


//-------------------------------------------------------------------------------
template <typename ScalarT>
class numval : public astNode<ScalarT>
{
  public:
    numval (ScalarT d): astNode<ScalarT>(),number(d) {};
    numval (std::complex<ScalarT> d): astNode<ScalarT>(),number(std::real(d)) {};

    virtual ScalarT val() { return number; }
    virtual ScalarT dx(int i) {return 0.0;}
    ScalarT number;

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "numval number = " << number << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this later for formatting
      os << number;
    }

    virtual bool numvalType() { return true; };
};

//-------------------------------------------------------------------------------
template <>
class numval<std::complex<double>> : public astNode<std::complex<double>>
{
  public:

    numval (std::complex<double> d): astNode<std::complex<double> >(),number(d) {};

    virtual std::complex<double> val() {return number;}
    virtual std::complex<double> dx(int i) {return (std::complex<double>(0.0,0.0));}
    std::complex<double> number;

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "numval number = " << number << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this later for formatting
      os << "std::complex<double>" << number;
    }

    virtual bool numvalType() { return true; };
};

//-------------------------------------------------------------------------------
//
// This is the parameter Op class.
//
// This is still a work in progress.  I originally wrote it to primarily behave 
// like a global_param.  ie, something that could dynamically change, and was 
// attached to an external expression via the "paramNode" object, which points 
// to the top of the ast tree of another expression.
//
// I have been modifying this slowly to make it so that it can also (under 
// some circumstances) behave like a Xyce .param.  ie, a hardwired constant.
// This aspect is not 100% fleshed out yet.  
//
// For a simple evaluation of the "val" function, this class will:
//
//    (1) call the "val" function of the underlying external syntax tree 
//    (2) return the scalar quantity "number_" 
//
// For derivative calculation, there is at least one use case that mixes these two 
// modes of operation together.  So I still need to think about that.  Possibly 
// there should be two different kinds of classes for this, to avoid confusion.
//
template <typename ScalarT>
class paramOp: public astNode<ScalarT>
{
  public:
    paramOp (std::string par):
      astNode<ScalarT>(),
      paramName_(par), 
      thisIsAFunctionArgument_(false),
      isVar(false),
      derivIndex_(-1)
  {
    numvalNode_ = Teuchos::rcp(new numval<ScalarT> (0.0));
    paramNode_ = numvalNode_;
    savedParamNode_ = numvalNode_;
  };

    paramOp (std::string par, Teuchos::RCP<astNode<ScalarT> > & tmpNode):
      astNode<ScalarT>(),
      paramName_(par), 
      paramNode_(tmpNode),
      savedParamNode_(tmpNode),
      thisIsAFunctionArgument_(false),
      isVar(false),
      derivIndex_(-1)
  {
    numvalNode_ = Teuchos::rcp(new numval<ScalarT> (0.0));
  };

    virtual ScalarT val() { return paramNode_->val(); }
    virtual ScalarT dx(int i) 
    { 
      ScalarT retval=0.0;
      if (isVar) { retval = (derivIndex_==i)?1.0:0.0; }
      else       { retval = paramNode_->dx(i); }
      return retval;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "parameter : " << paramName_ << " = " << val() << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << paramName_;
    }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) { paramNode_ = tmpNode; savedParamNode_ = tmpNode; };
    virtual void unsetNode() { paramNode_ = numvalNode_; };

    virtual void setValue(ScalarT val) { numvalNode_->number = val; paramNode_ = numvalNode_; };
    virtual void unsetValue() { paramNode_ = savedParamNode_; };

    virtual void setDerivIndex(int i) { derivIndex_=i; };
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    virtual std::string getName() { return paramName_; }

    virtual bool paramType() { return !thisIsAFunctionArgument_ ; };

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(paramNode_)) )
      {
        if (paramNode_->paramType())   { paramOpVector.push_back(paramNode_); }
        if (paramNode_->funcType())    { funcOpVector.push_back(paramNode_); }
        if (paramNode_->voltageType()) { voltOpVector.push_back(paramNode_); }
        if (paramNode_->currentType()) { currentOpVector.push_back(paramNode_); }
        paramNode_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(paramNode_)) )
      {
        if (paramNode_->paramType()) { paramOpVector.push_back(paramNode_); }
        paramNode_->getParamOps(paramOpVector);
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(paramNode_)) )
      {
        if (paramNode_->getFunctionArgType()) { funcArgOpVector.push_back(paramNode_); }
        paramNode_->getFuncArgOps(funcArgOpVector);
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(paramNode_)) )
      {
        if (paramNode_->funcType()) { funcOpVector.push_back(paramNode_); }
        paramNode_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(paramNode_)) )
      {
        if (paramNode_->voltageType()) { voltOpVector.push_back(paramNode_); }
        paramNode_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(paramNode_)) )
      {
        if (paramNode_->currentType()) { currentOpVector.push_back(paramNode_); }
        paramNode_->getCurrentOps(currentOpVector);
      }
    }

    virtual bool getFunctionArgType() { return thisIsAFunctionArgument_; };
    virtual void setFunctionArgType() { thisIsAFunctionArgument_ = true;};
    virtual void unsetFunctionArgType() { thisIsAFunctionArgument_ = true;};

    // the variable "isVar" is to support the old expression library API.
    // If true, it means that this parameter is one of the variables included 
    // in the "vars" array that is passed into the functions expression::evalauate 
    // and expression::evaluateFunction.
    void setVar() { isVar = true; }
    void unsetVar() { isVar = false; }
    bool getVar() { return isVar; }

  private:
    // data:
    std::string paramName_;

    Teuchos::RCP<astNode<ScalarT> > paramNode_;
    Teuchos::RCP<astNode<ScalarT> > savedParamNode_;
    Teuchos::RCP<numval<ScalarT> > numvalNode_;

    bool thisIsAFunctionArgument_;

    bool isVar;

    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class voltageOp: public astNode<ScalarT>
{
  public:
    voltageOp (std::vector<std::string> voltageNodes):
      astNode<ScalarT>(),
      voltageNodes_(voltageNodes),
      voltageVals_(voltageNodes.size(),0.0),
      number_(0.0),
      derivIndex_(-1)
    {
      int voltNodeSize = voltageNodes_.size();
      for (int ii=0;ii<voltNodeSize;++ii)
      {
        Xyce::Util::toUpper(voltageNodes_[ii]);
      }
    };

    virtual ScalarT val()
    {
      number_ = (voltageNodes_.size() == 2)?(voltageVals_[0]-voltageVals_[1]):voltageVals_[0];
      return number_;
    }

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Voltage node:" << std::endl;
      int voltNodeSize = voltageNodes_.size();
      if (voltNodeSize > 0)
      {
        for (int ii=0;ii<voltNodeSize;++ii)
        {
          os << std::setw(indent) << " ";
          os << "node " << ii << ":  V(" << voltageNodes_[ii] << ") = " << voltageVals_[ii] <<std::endl;
        }
      }
      os << std::setw(indent) << " ";
      os << "value = " << val() << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      if (voltageNodes_.size() == 1)
      {
        os << "V_";
        os << voltageNodes_[0];
      }
      else if (voltageNodes_.size() == 2)
      {
        os << "(V_";
        os << voltageNodes_[0];
        os << "-V_";
        os << voltageNodes_[1];
        os << ")";
      }
    }

    virtual void setVals(const std::vector<ScalarT> & vals)
    {
      int size = voltageVals_.size();
      if (vals.size() != size)
      {
        std::string tmp = "Voltage Args size don't match for V(";
        tmp+= (size > 0)?(voltageNodes_[0]):"";
        tmp+= (size > 1)?(","+voltageNodes_[1]):"";
        tmp += "). ";
        //tmp += "Size = " + size;
        std::vector<std::string> errStr(1,tmp);
        yyerror(errStr);
      }

      for (int ii=0;ii<voltageVals_.size();ii++) { voltageVals_[ii] = vals[ii]; }
    }

    virtual void setDerivIndex(int i) { derivIndex_=i; };
    virtual void unsetDerivIndex() { derivIndex_=-1; };

    std::vector<std::string> & getVoltageNodes() { return voltageNodes_; }
    std::vector<ScalarT> & getVoltageVals() { return voltageVals_; }

    virtual bool voltageType() { return true; };

    virtual std::vector<std::string> getNodeNames() { return voltageNodes_; }

  private:
    // data:
    std::vector<std::string> voltageNodes_;
    std::vector<ScalarT> voltageVals_;
    ScalarT number_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
template <typename ScalarT>
class currentOp: public astNode<ScalarT>
{
  public:
    currentOp (std::string currentDevice):
      astNode<ScalarT>(),
      number_(0.0),
      currentDevice_(currentDevice),
      derivIndex_(-1)
    {
      Xyce::Util::toUpper(currentDevice_);
    };

    virtual ScalarT val() {return number_;}

    virtual ScalarT dx(int i) { return (derivIndex_==i)?1.0:0.0; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Isrc : device = " << currentDevice_ <<std::endl;
      os << std::setw(indent) << " " << "value = " << val() <<std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "I_";
      os << currentDevice_;
    }

    virtual void setDerivIndex(int i) {derivIndex_=i;};
    virtual void unsetDerivIndex() {derivIndex_=-1;};

    void setCurrentDevice(const std::string & devName) { currentDevice_ = devName; }
    std::string getCurrentDevice() { return currentDevice_; }
    ScalarT getCurrentVal () { return number_; }
    void setCurrentVal (ScalarT n) { number_ = n; }

    virtual bool currentType() { return true; };

    virtual std::string getName () { return currentDevice_; }

  private:
// data:
    ScalarT number_;
    std::string currentDevice_;
    int derivIndex_;
};

//-------------------------------------------------------------------------------
// This class represents the invocation of a .FUNC (function).  It represents a
// specific usage of a function that was declared elsewhere.
//
// This needs a variable number of arguments, and those args that are passed
// in are the unique ones for this usage.
//
// The functionNode_ pointer is the "top" node of an AST tree for the function
// that is being invoked by this class.  It has to be resolved externally, as
// it was defined somewhere else.  It is a separate expression.
//
// The "dummy" args are also defined externally.  They are the generic arguments
// that were part of the .FUNC declaration.
//
// So, for example:
//
// .func f(x,y) { x+2*y } <- function declaration, set up elsewhere.
// { 10*f(2,3) }          <- expression using the function
//
//external to this AST:
// {x+2*y} = pointed to by functionNode_
// x,y     = dummyFuncArgs_  - parameter nodes, used by functionNode_ AST. Their setNode function must be called to set them to specific values
//
//internal to this funcOp class:
// 2,3     = funcArgs_ - generic AST nodes.  Could be anything.  These are the nodes that are "set" onto the dummyArgs. 
// f(2,3)  = function call handled by this class
//
//internal to this AST:
// 10      = handled elsewhere in this AST tree as a numval node
// *       = handled elsewhere in this AST tree as a binyarMultOp node
//
template <typename ScalarT>
class funcOp: public astNode<ScalarT>
{
  public:
    // functions:
    funcOp ( std::string name, std::vector<Teuchos::RCP<astNode<ScalarT> > > * args):
      astNode<ScalarT>(),
      funcName_(name),
      funcArgs_(*args),
      number_(0.0),
      nodeResolved_(false),
      argsResolved_(false)
    {};

    virtual ScalarT val()
    {
      if (nodeResolved_ && argsResolved_)
      {
        if (funcArgs_.size() != dummyFuncArgs_.size())
        {
          std::vector<std::string> errStr;
          errStr.push_back(std::string("FuncOp Function Args sizes don't match for: "));
          errStr.push_back(funcName_); 
          errStr.push_back(std::string("funcArgs size = ") + std::to_string(funcArgs_.size()) );
          errStr.push_back(std::string("dummyFuncArgs size = ") + std::to_string(funcArgs_.size()));
          yyerror(errStr);
        }
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
        number_ = functionNode_->val();
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
      }
      else
      {
        std::string tmp = "Function " + funcName_ + " is not resolved";
        std::vector<std::string> errStr(1,tmp);
        yyerror(errStr);
      }
      return number_;
    }

    virtual ScalarT dx(int i)
    {
      if (nodeResolved_ && argsResolved_)
      {
        if (funcArgs_.size() != dummyFuncArgs_.size())
        {
          std::string tmp = "FuncOp Function Args sizes don't match for " + funcName_;
          std::vector<std::string> errStr(1,tmp);
          errStr.push_back(std::string("funcArgs size = ") + std::to_string(funcArgs_.size()) );
          errStr.push_back(std::string("dummyFuncArgs size = ") + std::to_string(funcArgs_.size()));
          yyerror(errStr);
        }

        // Two phases, do do a complete chain rule calculation: 
        //
        // all the "d" symbols should be partials:
        // chain rule :  F′(x) = F'(x) + f′(g(x)) * g′(x)  ->  df/dx = df/dx + df/dp * dp/dx 
        //
        // phase 1:  F'(x) = df/dx.
        //
        // For this phase, the funcArgs are in the "full" form - 
        //   ie, if they represent an AST tree, we use the whole tree to evaluate.
        number_ = 0.0;
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
        number_ = functionNode_->dx(i);
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore

        //std::cout << "number_ = " << number_ << std::endl;

        // phase 2:  f′(g(x)) * g′(x) = df/dp * dp/dx
        //
        // g(x) = funcArg->val().  This should be evaluated inside of dx call.
        // g'(x) = funcArg->dx(i)
        // f'(g(x)) = functionNode_->dx(ii);
        //
        // For this phase, the funcArgs are simple params.  
        //   ie, they don't have an AST tree, just a number.
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) 
        { 
          // the index is intentionally negative, starting at -1.  This is so it 
          // doesn't conflict with derivative indices that were already set at 
          // the top of the tree in the newExpression::evaluate function.
          int index=-ii-1;
          dummyFuncArgs_[ii]->setValue ( funcArgs_[ii]->val() ); 
          dummyFuncArgs_[ii]->setDerivIndex ( index );
        }

        for (int ii=0;ii<dummyFuncArgs_.size();++ii)  // loop over args (p).  ii = p index, i = x index
        {
          int index=-ii-1;
          ScalarT delta = functionNode_->dx(index) *  funcArgs_[ii]->dx(i);
          number_ += delta; // functionNode_->dx(index) *  funcArgs_[ii]->dx(i);
          //std::cout << "ii="<< ii << "  functionNode_->dx(-ii-1) = " << functionNode_->dx(index) << "  funcArgs_[ii]->dx(i) = " << funcArgs_[ii]->dx(i) << "  delta = " << delta << " number_ = " << number_ << std::endl;
        }

        for (int ii=0;ii<dummyFuncArgs_.size();++ii) 
        { 
          dummyFuncArgs_[ii]->unsetValue ();
          dummyFuncArgs_[ii]->unsetDerivIndex ();
        } // restore

      }
      return number_;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " " << "function: " << funcName_ << ":" <<std::endl;;
      os << std::setw(indent) << " " << "function args: " << std::endl; indent++;
      for (int ii=0;ii<funcArgs_.size();++ii) { funcArgs_[ii]->output(os,indent+1); }

      if( !(Teuchos::is_null(functionNode_)) )
      {
        os << std::setw(indent) << " " << "functionNode_ ("<<funcName_<<") details: " << std::endl;


        if (funcArgs_.size() != dummyFuncArgs_.size())
        {
          std::vector<std::string> errStr;
          errStr.push_back(std::string("FuncOp Function Args sizes don't match for: "));
          errStr.push_back(funcName_); 
          errStr.push_back(std::string("funcArgs size = ") + std::to_string(funcArgs_.size()) );
          errStr.push_back(std::string("dummyFuncArgs size = ") + std::to_string(funcArgs_.size()));
          yyerror(errStr);
        }

        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
        functionNode_->output(os,indent+2);
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
      }
      else
      {
        os << std::setw(indent) << " " << "functionNode_ is not resolved " << std::endl;
      }
      os << std::setw(indent) << " " << "val = " << val() << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << funcName_;
      os << "(";
      int size = funcArgs_.size();
      for (int ii=0;ii<size;ii++)
      {
        funcArgs_[ii]->codeGen(os);
        if (ii < size-1) {os << ",";}
      }
      os << ")";
    }

    virtual void setNode(Teuchos::RCP<astNode<ScalarT> > & tmpNode) { functionNode_ = tmpNode; nodeResolved_ = true; };
    virtual void unsetNode()
    {
      //functionNode_ = NULL;
      nodeResolved_ = false;
      number_ = 0.0;
    };

    virtual void setFuncArgs(std::vector< Teuchos::RCP<paramOp<ScalarT> > > & tmpParamVec )
    {
      dummyFuncArgs_.clear(); dummyFuncArgs_.resize(tmpParamVec.size());
      for (int ii=0;ii<tmpParamVec.size();++ii)
      {
        dummyFuncArgs_[ii] = tmpParamVec[ii];
      }
      argsResolved_ = true;
    };

    virtual void setFuncArgs(std::vector< Teuchos::RCP<astNode<ScalarT> > > & tmpArgVec )
    {
      dummyFuncArgs_.clear(); dummyFuncArgs_.resize(tmpArgVec.size());
      for (int ii=0;ii<tmpArgVec.size();++ii)
      {
        dummyFuncArgs_[ii] = tmpArgVec[ii];
      }
      argsResolved_ = true;
    };

    virtual bool funcType()    { return true; };

    virtual std::string getName() { return funcName_; }

    // still experimenting with these:
    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(functionNode_)) )
      {
        if (functionNode_->paramType())   { paramOpVector.push_back(functionNode_); }
        if (functionNode_->funcType())    { funcOpVector.push_back(functionNode_); }
        if (functionNode_->voltageType()) { voltOpVector.push_back(functionNode_); }
        if (functionNode_->currentType()) { currentOpVector.push_back(functionNode_); }

        if (funcArgs_.size() != dummyFuncArgs_.size())
        {
          std::vector<std::string> errStr;
          errStr.push_back(std::string("FuncOp Function Args sizes don't match for: "));
          errStr.push_back(funcName_); 
          errStr.push_back(std::string("funcArgs size = ") + std::to_string(funcArgs_.size()) );
          errStr.push_back(std::string("dummyFuncArgs size = ") + std::to_string(funcArgs_.size()));
          yyerror(errStr);
        }

        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->setNode( funcArgs_[ii] ); }
        functionNode_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
        for (int ii=0;ii<dummyFuncArgs_.size();++ii) { dummyFuncArgs_[ii]->unsetNode(); } // restore
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(functionNode_)) )
      {
        if (functionNode_->paramType()) { paramOpVector.push_back(functionNode_); }
        functionNode_->getParamOps(paramOpVector);
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(functionNode_)) )
      {
        if (functionNode_->getFunctionArgType()) { funcArgOpVector.push_back(functionNode_); }
        functionNode_->getFuncArgOps(funcArgOpVector);
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(functionNode_)) )
      {
        if (functionNode_->funcType()) { funcOpVector.push_back(functionNode_); }
        functionNode_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(functionNode_)) )
      {
        if (functionNode_->voltageType()) { voltOpVector.push_back(functionNode_); }
        functionNode_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(functionNode_)) )
      {
        if (functionNode_->currentType()) { currentOpVector.push_back(functionNode_); }
        functionNode_->getCurrentOps(currentOpVector);
      }
    }

  private:
// data:
    std::string funcName_;
    std::vector<Teuchos::RCP<astNode<ScalarT> > > funcArgs_;  // the unique args that are passed in to this instance of a function
    std::vector<Teuchos::RCP<astNode<ScalarT> > > dummyFuncArgs_;  // generic args that the functionNode_ owns; ie, that are used to evaluate it.  They have to be temporarily replaced whenever the function is called.
    ScalarT number_;
    Teuchos::RCP<astNode<ScalarT> > functionNode_;
    bool nodeResolved_;
    bool argsResolved_;
};

//-------------------------------------------------------------------------------
// sign corrected x raised to y power
//
//  pwrs(x,y) = x^y     if x > 0
//  pwrs(x,y) = 0       if x = 0
//  pwrs(x,y) = -(-x)^y if x < 0
//
template <typename ScalarT>
class pwrsOp : public astNode<ScalarT>
{
  public:
    pwrsOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right) {};

    virtual ScalarT val()
    {
      ScalarT ret=0.0;
      if (std::real(this->leftAst_->val()) >= 0)
      {
        ret = std::pow(this->leftAst_->val(), this->rightAst_->val());
      }
      else if (std::real(this->leftAst_->val()) < 0)
      {
        ret = std::pow(-(this->leftAst_->val()), this->rightAst_->val());
      }
      return ret;
    }

    virtual ScalarT dx (int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (pwrs) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
// fix this:
#if 0
      astNode<ScalarT> & lef = *(this->leftAst_);
      astNode<ScalarT> & rig = *(this->rightAst_);
      return  ((rig.dx(i)*std::log(lef.val())+rig.val()*lef.dx(i)/lef.val())*std::pow(lef.val(),rig.val()));
#endif
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "pwrs operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      if (std::real(this->leftAst_->val()) < 0) { os << "-"; }
      os << "std::pow(";
      this->leftAst_->codeGen(os);
      os << ",";
      this->rightAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
//
// +1 if x > 0
//  0 if x = 0
// -1 if x < 0
//
template <typename ScalarT>
class sgnOp : public astNode<ScalarT>
{
  public:
    sgnOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      ScalarT ret = 0.0;
      ret = ( (std::real(this->leftAst_->val())>0)?+1:ret );
      ret = ( (std::real(this->leftAst_->val())<0)?-1:ret );
      return ret;
    }

    virtual ScalarT dx(int i)
    {
#if 0
      return (-(this->leftAst_->dx(i)));
#else
      return ScalarT(0.0);
#endif
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "sgn operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "signbit(";
      this->leftAst_->codeGen(os);
      os << ")";
    }
};

//-------------------------------------------------------------------------------
//  sign(x,y) = sgn(y)|x|  sign of y times absolute value of x
template <typename ScalarT>
class signOp : public astNode<ScalarT>
{
  public:
    signOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right): astNode<ScalarT>(left,right) {};

    virtual ScalarT val()
    {
      ScalarT y = 0.0;
      y = ( (std::real(this->rightAst_->val())>0)?+1:y );
      y = ( (std::real(this->rightAst_->val())<0)?-1:y );

      ScalarT x =  (std::abs(this->leftAst_->val()));
      return (y*x);
    }

    virtual ScalarT dx (int i)
    {
      ScalarT y = 0.0;
      y = ( (std::real(this->rightAst_->val())>0)?+1:y );
      y = ( (std::real(this->rightAst_->val())<0)?-1:y );
      ScalarT dx = (std::real(this->leftAst_->val()) >= 0 ? this->leftAst_->dx(i) : ScalarT(-this->leftAst_->dx(i)));
      return (y*dx);
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "sign(x,y) = (sgn(y)|x|) op " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "signbit(";
      this->rightAst_->codeGen(os);
      os << ")*std::abs(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// least integer greater or equal to variable x
template <typename ScalarT>
class ceilOp : public astNode<ScalarT>
{
  public:
    ceilOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::ceil(std::real(this->leftAst_->val())); }

    virtual ScalarT dx(int i)
    {
      // derivative is undefined at integers and 0.0 elsewhere
      return  0.0;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ceil operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }
   
    virtual void codeGen (std::ostream & os )
    {
      os << "std::ceil(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// greatest integer less than or equal to variable x
template <typename ScalarT>
class floorOp : public astNode<ScalarT>
{
  public:
    floorOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val() { return std::floor(std::real(this->leftAst_->val())); }

    virtual ScalarT dx(int i)
    {
      // derivative is undefined at integers and 0.0 elsewhere
      return  0.0;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "floor operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "std::floor(";
      this->leftAst_->codeGen(os);
      os << ")";
    }

};

//-------------------------------------------------------------------------------
// integer part of the real variable x
template <typename ScalarT>
class intOp : public astNode<ScalarT>
{
  public:
    intOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      int tmp = std::real(this->leftAst_->val()) ;
      ScalarT ret = tmp;
      return ret;
    }

    virtual ScalarT dx(int i) { return  0.0; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "int operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "static_cast<int>( std::real(";
      this->leftAst_->codeGen(os);
      os << "))";
    }

};

//-------------------------------------------------------------------------------
// if statement op (or ternary op)
template <typename ScalarT>
class ifStatementOp : public astNode<ScalarT>
{
  public:
    ifStatementOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &zAst):
      astNode<ScalarT>(xAst,yAst),
      zAst_(zAst) {};

    virtual ScalarT val()
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);
      return ((std::real(x->val()))?(y->val()):(z->val()));
    };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);
      return ((std::real(x->val()))?(y->dx(i)):(z->dx(i)));
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "if statement operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      zAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "((";
      this->leftAst_->codeGen(os);
      os << ")?(";
      this->rightAst_->codeGen(os);
      os << "):(";
      zAst_->codeGen(os);
      os << "))";
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType())   { paramOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->funcType())    { funcOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType())   { paramOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->funcType())    { funcOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->paramType())   { paramOpVector.push_back(zAst_); }
        if (zAst_->funcType())    { funcOpVector.push_back(zAst_); }
        if (zAst_->voltageType()) { voltOpVector.push_back(zAst_); }
        if (zAst_->currentType()) { currentOpVector.push_back(zAst_); }
        zAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType()) { paramOpVector.push_back(this->leftAst_); }
        this->leftAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType()) { paramOpVector.push_back(this->rightAst_); }
        this->rightAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->paramType()) { paramOpVector.push_back(zAst_); }
        zAst_->getParamOps(paramOpVector);
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->getFunctionArgType()) { funcArgOpVector.push_back(zAst_); }
        zAst_->getFuncArgOps(funcArgOpVector);
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->funcType()) { funcOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->funcType()) { funcOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->funcType()) { funcOpVector.push_back(zAst_); }
        zAst_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        this->leftAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        this->rightAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->voltageType()) { voltOpVector.push_back(zAst_); }
        zAst_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->currentType()) { currentOpVector.push_back(zAst_); }
        zAst_->getCurrentOps(currentOpVector);
      }
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > zAst_;
};

//-------------------------------------------------------------------------------
// x limited to range y to z
// y if x < y
// x if y < x < z
// z if x > z
template <typename ScalarT>
class limitOp : public astNode<ScalarT>
{
  public:
    limitOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &zAst):
      astNode<ScalarT>(xAst,yAst),
      zAst_(zAst) {};

    virtual ScalarT val()
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      return ((std::real(x->val())<std::real(y->val()))?
          std::real(y->val()):
          ((std::real(x->val())>std::real(z->val()))?std::real(z->val()):std::real(x->val())));
    };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & z = (zAst_);

      return ((std::real(x->val())<std::real(y->val()))?0.0:((std::real(x->val())>std::real(z->val()))?0.0:1.0));

    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "limit operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      zAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      /// fix this
      os << "LIMIT";
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType())   { paramOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->funcType())    { funcOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType())   { paramOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->funcType())    { funcOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->paramType())   { paramOpVector.push_back(zAst_); }
        if (zAst_->funcType())    { funcOpVector.push_back(zAst_); }
        if (zAst_->voltageType()) { voltOpVector.push_back(zAst_); }
        if (zAst_->currentType()) { currentOpVector.push_back(zAst_); }
        zAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType()) { paramOpVector.push_back(this->leftAst_); }
        this->leftAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType()) { paramOpVector.push_back(this->rightAst_); }
        this->rightAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->paramType()) { paramOpVector.push_back(zAst_); }
        zAst_->getParamOps(paramOpVector);
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->getFunctionArgType()) { funcArgOpVector.push_back(zAst_); }
        zAst_->getFuncArgOps(funcArgOpVector);
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->funcType()) { funcOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->funcType()) { funcOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->funcType()) { funcOpVector.push_back(zAst_); }
        zAst_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        this->leftAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        this->rightAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->voltageType()) { voltOpVector.push_back(zAst_); }
        zAst_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(zAst_)) )
      {
        if (zAst_->currentType()) { currentOpVector.push_back(zAst_); }
        zAst_->getCurrentOps(currentOpVector);
      }
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > zAst_;
};


//-------------------------------------------------------------------------------
// step function :  1 if x > 0
template <typename ScalarT>
class stpOp : public astNode<ScalarT>
{
  public:
    stpOp (Teuchos::RCP<astNode<ScalarT> > &left):
      astNode<ScalarT>(left) {};

    virtual ScalarT val() { return ((std::real(this->leftAst_->val()))>0)?1.0:0.0; }

    virtual ScalarT dx (int i) { return 0.0; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "step function operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      /// fix this
      os << "STP";
    }

};

//-------------------------------------------------------------------------------
// ramp function x if x > 0; 0 otherwise
template <typename ScalarT>
class urampOp : public astNode<ScalarT>
{
  public:
    urampOp (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      return ((std::real(this->leftAst_->val()))>0)?(std::real(this->leftAst_->val())):0.0;
    }

    virtual ScalarT dx (int i)
    {
      return ((std::real(this->leftAst_->val()))>0)?1.0:0.0;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "uramp operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "(((std::real(";
      this->leftAst_->codeGen(os);
      os << "))>0)?(std::real(";
      this->leftAst_->codeGen(os);
      os << ")):0.0)";
    }

};

//-------------------------------------------------------------------------------
// TABLE(x,y,z,*)
// f(x) where f(y) = z
// piecewise linear interpolation, multiple (y,z) pairs can be specified
template <typename ScalarT>
class tableOp : public astNode<ScalarT>
{
  public:
    // functions:
    tableOp (Teuchos::RCP<astNode<ScalarT> > &input, std::vector<Teuchos::RCP<astNode<ScalarT> > > * args):
      astNode<ScalarT>(), tableArgs_(*args), allNumVal_(true), input_(input)
      {
        int size = tableArgs_.size();
        if (size % 2)
        {
          std::vector<std::string> errStr(1,std::string("AST node (table) needs an even number of arguments")); yyerror(errStr);
        }
        else
        {
          allNumVal_=true; ta_.resize(size/2); ya_.resize(size/2); dya_.resize(size/2,0.0);
          for (int ii=0,jj=0;ii<size;ii+=2,jj++)
          {
            ta_[jj] = (tableArgs_)[ii]->val();
            ya_[jj] = (tableArgs_)[ii+1]->val();
            if (!( (tableArgs_)[ii]->numvalType() && (tableArgs_)[ii+1]->numvalType() ) ) { allNumVal_ = false; }
          }
        }
      };

    // special constructor for values read in from a file, that are now stored in std::vector objects
    tableOp (Teuchos::RCP<astNode<ScalarT> > & input, const std::vector<ScalarT> & xvals, const std::vector<ScalarT> & yvals):
      astNode<ScalarT>(), allNumVal_(true), input_(input)
      {
        allNumVal_=true; int size = xvals.size(); int size2=yvals.size();
        if (size != size2)
        {
          std::vector<std::string> errStr(1,std::string("AST node (table) needs x and y vectors to be the same size.")); yyerror(errStr);
        }
        ta_.resize(size); ya_.resize(size); dya_.resize(size,0.0);

        for (int ii=0;ii<size;ii++)
        {
          ta_[ii] = xvals[ii];
          ya_[ii] = yvals[ii];
        }
      };

    virtual ScalarT val()
    {
      ScalarT y = 0.0;
      if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again
      {
        int size = tableArgs_.size();
        for (int ii=0,jj=0;ii<size;ii+=2,jj++)
        {
          ta_[jj] = (tableArgs_)[ii]->val();
          ya_[jj] = (tableArgs_)[ii+1]->val();
        }
        yInterpolator_.init(ta_,ya_); // for linear, this isn't necessary, but for others it is
      }
      ScalarT input = std::real(this->input_->val());
      yInterpolator_.eval(ta_,ya_, input, y);
      return y;
    };

    virtual ScalarT dx(int i)
    {
      ScalarT dydx = 0.0;
      if (!allNumVal_)  // if not all pure numbers, then initialize the arrays again.  Otherwise, dydx=0
      {
        int size = tableArgs_.size();
        for (int ii=0,jj=0;ii<size;ii+=2,jj++)
        {
          ta_[jj] = (tableArgs_)[ii]->val();
          dya_[jj] = (tableArgs_)[ii+1]->dx(i);
        }
        dyInterpolator_.init(ta_,dya_); // for linear, this isn't necessary, but for others it is
        ScalarT input = std::real(this->input_->val());
        dyInterpolator_.eval(ta_,dya_, input, dydx);
      }
      return dydx;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "table operator " << std::endl;
      //++indent;
      //this->input_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "TABLE";
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {

      if( !(Teuchos::is_null(input_)) )
      {
        if (input_->paramType())   { paramOpVector.push_back(input_); }
        if (input_->funcType())    { funcOpVector.push_back(input_); }
        if (input_->voltageType()) { voltOpVector.push_back(input_); }
        if (input_->currentType()) { currentOpVector.push_back(input_); }
        input_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          if( !(Teuchos::is_null(tableArgs_[ii])) )
          {
            if (tableArgs_[ii]->paramType())   { paramOpVector.push_back(tableArgs_[ii]); }
            if (tableArgs_[ii]->funcType())    { funcOpVector.push_back(tableArgs_[ii]); }
            if (tableArgs_[ii]->voltageType()) { voltOpVector.push_back(tableArgs_[ii]); }
            if (tableArgs_[ii]->currentType()) { currentOpVector.push_back(tableArgs_[ii]); }
            tableArgs_[ii]->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
          }
        }
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(input_)) )
      {
        if (input_->paramType()) { paramOpVector.push_back(input_); }
        input_->getParamOps(paramOpVector);
      }

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          if( !(Teuchos::is_null(tableArgs_[ii])) )
          {
            if (tableArgs_[ii]->paramType()) { paramOpVector.push_back(tableArgs_[ii]); }
            tableArgs_[ii]->getParamOps(paramOpVector);
          }
        }
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(input_)) )
      {
        if (input_->getFunctionArgType()) { funcArgOpVector.push_back(input_); }
        input_->getFuncArgOps(funcArgOpVector);
      }

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          if( !(Teuchos::is_null(tableArgs_[ii])) )
          {
            if (tableArgs_[ii]->getFunctionArgType()) { funcArgOpVector.push_back(tableArgs_[ii]); }
            tableArgs_[ii]->getFuncArgOps(funcArgOpVector);
          }
        }
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(input_)) )
      {
        if (input_->funcType()) { funcOpVector.push_back(input_); }
        input_->getFuncOps(funcOpVector);
      }

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          if( !(Teuchos::is_null(tableArgs_[ii])) )
          {
            if (tableArgs_[ii]->funcType()) { funcOpVector.push_back(tableArgs_[ii]); }
            tableArgs_[ii]->getFuncOps(funcOpVector);
          }
        }
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(input_)) )
      {
        if (input_->voltageType()) { voltOpVector.push_back(input_); }
        input_->getVoltageOps(voltOpVector);
      }

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          if( !(Teuchos::is_null(tableArgs_[ii])) )
          {
            if (tableArgs_[ii]->voltageType()) { voltOpVector.push_back(tableArgs_[ii]); }
            tableArgs_[ii]->getVoltageOps(voltOpVector);
          }
        }
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(input_)) )
      {
        if (input_->currentType()) { currentOpVector.push_back(input_); }
        input_->getCurrentOps(currentOpVector);
      }

      if (!allNumVal_)
      {
        int size=tableArgs_.size();
        for(int ii=0;ii<size;ii++)
        {
          if( !(Teuchos::is_null(tableArgs_[ii])) )
          {
            if (tableArgs_[ii]->currentType()) { currentOpVector.push_back(tableArgs_[ii]); }
            tableArgs_[ii]->getCurrentOps(currentOpVector);
          }
        }
      }
    }

  private:
// data:
    std::vector<Teuchos::RCP<astNode<ScalarT> > > tableArgs_;
    bool allNumVal_;
    std::vector<ScalarT> ta_; // using ta for name instead of xa so as not to confuse meaning of dx function
    std::vector<ScalarT> ya_;
    std::vector<ScalarT> dya_;
    Xyce::Util::linear<ScalarT> yInterpolator_; // possibly make this a user choice
    Xyce::Util::linear<ScalarT> dyInterpolator_; // possibly make this a user choice
    Teuchos::RCP<astNode<ScalarT> > input_;
};

//-------------------------------------------------------------------------------
// time integral of x
template <typename ScalarT>
class sdtOp : public astNode<ScalarT>
{
  public:
    sdtOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      std::vector<std::string> errStr(1,std::string("AST node (sdt) without an val function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (sdt) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "sdt (time integral) operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "SDT";
    }

};

//-------------------------------------------------------------------------------
// time derivative of x
template <typename ScalarT>
class ddtOp : public astNode<ScalarT>
{
  public:
    ddtOp(Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};

    virtual ScalarT val()
    {
      std::vector<std::string> errStr(1,std::string("AST node (ddt) without an val function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (ddt) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ddt (time integral) operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "DDT";
    }

};

//-------------------------------------------------------------------------------
// ddx(f(x),x)
// partial derivative of f (x) with respect to x
template <typename ScalarT>
class ddxOp : public astNode<ScalarT>
{
  public:
    ddxOp (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):
      astNode<ScalarT>(left,right) , foundX_(false)
    {
      resolveArg_();
    }

    void resolveArg_()
    {
      // check for params
      //
      // For most uses of "paramType" it is important to exclude function 
      // arguments, which are also allocated as parameters, but are usually 
      // not used in the same way as parameters.  This is the one use case where
      // function arguments must be included.
      //
      // ERK:
      // Note, as currently written, these are very simple linear searches
      // to match the differential variable with a parameter (or voltage or current)
      //
      // This is fine for expressions with very small numbers of parameters, etc.
      // Not fine for really large containers.  Fix later.
      //
      if ( this->rightAst_->paramType() || this->rightAst_->getFunctionArgType() )
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > paramOpVector;
        if ( this->leftAst_->paramType() || this->leftAst_->getFunctionArgType() ) 
        { paramOpVector.push_back(this->leftAst_); }

        this->leftAst_->getParamOps(paramOpVector); 
        this->leftAst_->getFuncArgOps(paramOpVector);

        // now match the user-specified right hand argument with a parameter inside
        // the function specified by the left-hand argument.  This is a klunky,
        // inefficient search and should be replaced with something better.
        // (same is true below for voltages and currents)
        std::string tmp = this->rightAst_->getName();
        if (!(tmp.empty()))
        {
          Xyce::Util::toUpper(tmp);
          for (int ii=0;ii<paramOpVector.size();ii++)
          {
            std::string tmp2 = paramOpVector[ii]->getName();
            if (!(tmp2.empty()))
            {
              Xyce::Util::toUpper(tmp2);
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = paramOpVector[ii]; }
            }
          }
        }
      }
      else if (this->rightAst_->voltageType())
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > voltOpVector;
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }

        this->leftAst_->getVoltageOps(voltOpVector);

        std::vector<std::string> tmp = this->rightAst_->getNodeNames();
        if (!(tmp.empty()))
        {
          for (int jj=0;jj<tmp.size();jj++) Xyce::Util::toUpper(tmp[jj]);

          for (int ii=0;ii<voltOpVector.size();ii++)
          {
            std::vector<std::string> tmp2 = this->rightAst_->getNodeNames();
            if (!(tmp2.empty()))
            {
              for (int jj=0;jj<tmp2.size();jj++) Xyce::Util::toUpper(tmp2[jj]);

              if (tmp.size() == tmp2.size())
              {
                if (tmp==tmp2) { foundX_ = true; astNodeX_ = voltOpVector[ii]; }
              }
            }
          }
        }
      }
      else if (this->rightAst_->currentType())
      {
        std::vector<Teuchos::RCP<astNode<ScalarT> > > currentOpVector;

        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }

        this->leftAst_->getCurrentOps(currentOpVector);

        std::string tmp = this->rightAst_->getName();
        if (!(tmp.empty()))
        {
          Xyce::Util::toUpper(tmp);
          for (int ii=0;ii<currentOpVector.size();ii++)
          {
            std::string tmp2 = currentOpVector[ii]->getName();
            if (!(tmp2.empty()))
            {
              Xyce::Util::toUpper(tmp2);
              if (tmp==tmp2) { foundX_ = true; astNodeX_ = currentOpVector[ii]; }
            }
          }
        }
      }
      else // unsupported type
      {
        std::vector<std::string> errStr(1,std::string("DDX unsupported type"));
        yyerror(errStr);
      }

      if (!foundX_)
      {
        std::vector<std::string> errStr(1,std::string("DDX argument not resolved"));
        yyerror(errStr);
      }
    };

    virtual ScalarT val()
    {
      ScalarT ret = 0.0;
      if( !(Teuchos::is_null( astNodeX_)))
      {
        astNodeX_->setDerivIndex(0);
        astNodeX_->setVar();
        ret = this->leftAst_->dx(0);
        astNodeX_->unsetDerivIndex();
        astNodeX_->unsetVar();
      }
      return ret;
    };

    virtual ScalarT dx(int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (ddx) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "ddx (time integral) operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "DDX";
    }

  private:
    bool foundX_;
    Teuchos::RCP<astNode<ScalarT> > astNodeX_;
};

//-------------------------------------------------------------------------------
// Random number sampled from normal distribution with
// mean μ and standard deviation (α)/n
template <typename ScalarT>
class agaussOp : public astNode<ScalarT>
{
  public:
    agaussOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &nAst):
      astNode<ScalarT>(xAst,yAst),
      nAst_(nAst) {};

    virtual ScalarT val()
    {
      Teuchos::RCP<astNode<ScalarT> > & mu    = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & alpha = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & n     = (nAst_);
      std::vector<std::string> errStr(1,std::string("AST node (agauss) without an val function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & mu    = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & alpha = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & n     = (nAst_);
      std::vector<std::string> errStr(1,std::string("AST node (agauss) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "agauss operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      nAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "AGAUSS";
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType())   { paramOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->funcType())    { funcOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType())   { paramOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->funcType())    { funcOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->paramType())   { paramOpVector.push_back(nAst_); }
        if (nAst_->funcType())    { funcOpVector.push_back(nAst_); }
        if (nAst_->voltageType()) { voltOpVector.push_back(nAst_); }
        if (nAst_->currentType()) { currentOpVector.push_back(nAst_); }
        nAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType()) { paramOpVector.push_back(this->leftAst_); }
        this->leftAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType()) { paramOpVector.push_back(this->rightAst_); }
        this->rightAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->paramType()) { paramOpVector.push_back(nAst_); }
        nAst_->getParamOps(paramOpVector);
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->getFunctionArgType()) { funcArgOpVector.push_back(nAst_); }
        nAst_->getFuncArgOps(funcArgOpVector);
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->funcType()) { funcOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->funcType()) { funcOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->funcType()) { funcOpVector.push_back(nAst_); }
        nAst_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        this->leftAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        this->rightAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->voltageType()) { voltOpVector.push_back(nAst_); }
        nAst_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->currentType()) { currentOpVector.push_back(nAst_); }
        nAst_->getCurrentOps(currentOpVector);
      }
    }


  private:
    Teuchos::RCP<astNode<ScalarT> > nAst_;
};

//-------------------------------------------------------------------------------
// Random number sampled from normal distribution with
// mean μ and standard deviation (α ∗ μ )/n
template <typename ScalarT>
class gaussOp : public astNode<ScalarT>
{
  public:
    gaussOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &nAst):
      astNode<ScalarT>(xAst,yAst),
      nAst_(nAst) {};

    virtual ScalarT val()
    {
      Teuchos::RCP<astNode<ScalarT> > & mu    = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & alpha = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & n     = (nAst_);
      std::vector<std::string> errStr(1,std::string("AST node (gauss) without an val function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & mu    = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & alpha = (this->rightAst_);
      Teuchos::RCP<astNode<ScalarT> > & n     = (nAst_);
      std::vector<std::string> errStr(1,std::string("AST node (gauss) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "gauss operator " << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      nAst_->output(os,indent+1);
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "GAUSS";
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType())   { paramOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->funcType())    { funcOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType())   { paramOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->funcType())    { funcOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->paramType())   { paramOpVector.push_back(nAst_); }
        if (nAst_->funcType())    { funcOpVector.push_back(nAst_); }
        if (nAst_->voltageType()) { voltOpVector.push_back(nAst_); }
        if (nAst_->currentType()) { currentOpVector.push_back(nAst_); }
        nAst_->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector);
      }
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->paramType()) { paramOpVector.push_back(this->leftAst_); }
        this->leftAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->paramType()) { paramOpVector.push_back(this->rightAst_); }
        this->rightAst_->getParamOps(paramOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->paramType()) { paramOpVector.push_back(nAst_); }
        nAst_->getParamOps(paramOpVector);
      }
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->getFunctionArgType()) { funcArgOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncArgOps(funcArgOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->getFunctionArgType()) { funcArgOpVector.push_back(nAst_); }
        nAst_->getFuncArgOps(funcArgOpVector);
      }
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->funcType()) { funcOpVector.push_back(this->leftAst_); }
        this->leftAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->funcType()) { funcOpVector.push_back(this->rightAst_); }
        this->rightAst_->getFuncOps(funcOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->funcType()) { funcOpVector.push_back(nAst_); }
        nAst_->getFuncOps(funcOpVector);
      }
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->voltageType()) { voltOpVector.push_back(this->leftAst_); }
        this->leftAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->voltageType()) { voltOpVector.push_back(this->rightAst_); }
        this->rightAst_->getVoltageOps(voltOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->voltageType()) { voltOpVector.push_back(nAst_); }
        nAst_->getVoltageOps(voltOpVector);
      }
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
      if( !(Teuchos::is_null(this->leftAst_)) )
      {
        if (this->leftAst_->currentType()) { currentOpVector.push_back(this->leftAst_); }
        this->leftAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(this->rightAst_)) )
      {
        if (this->rightAst_->currentType()) { currentOpVector.push_back(this->rightAst_); }
        this->rightAst_->getCurrentOps(currentOpVector);
      }

      if( !(Teuchos::is_null(nAst_)) )
      {
        if (nAst_->currentType()) { currentOpVector.push_back(nAst_); }
        nAst_->getCurrentOps(currentOpVector);
      }
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > nAst_;
};

//-------------------------------------------------------------------------------
// random number between 0 and 1 sampled from a uniform distribution
template <typename ScalarT>
class randOp : public astNode<ScalarT>
{
  public:
    randOp (): astNode<ScalarT>() {};

    virtual ScalarT val()
    {
      std::vector<std::string> errStr(1,std::string("AST node (rand) without an val function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual ScalarT dx (int i)
    {
      std::vector<std::string> errStr(1,std::string("AST node (rand) without a dx function"));
      yyerror(errStr);
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "rand operator " << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "RAND";
    }

  private:
};

//-------------------------------------------------------------------------------
// specials
//  {"TIME", "TEMP", "VT", "FREQ"}
template <typename ScalarT>
class specialsOp : public astNode<ScalarT>
{
  public:
    specialsOp (std::string typeName) : astNode<ScalarT>(), type_(typeName), value_(0.0) {};

    virtual ScalarT val() { return value_; };
    virtual ScalarT dx (int i) { return ScalarT(0.0); };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << type_ << " operator.  val = " << value_ << std::endl;
    }

    virtual void codeGen (std::ostream & os ) { os << value_; }

    std::string getType() { return type_; }
    void getType(std::string t) { type_ = t; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; }

  private:
    std::string type_;
    ScalarT value_;
};

//-------------------------------------------------------------------------------
// constants
template <typename ScalarT>
class piConstOp : public astNode<ScalarT>
{
  public:
    piConstOp (): astNode<ScalarT>() {};

    virtual ScalarT val() { return ScalarT(M_PI); };
    virtual ScalarT dx (int i) { return ScalarT(0.0); };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "pi const operator.  val = " << ScalarT(M_PI) << std::endl;
    }

    virtual void codeGen (std::ostream & os ) { os << ScalarT(M_PI); }

  private:
};

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
inline void yyerror(std::vector<std::string> & s)
{
  for (int i=0;i<s.size();++i)
  {
    std::cerr << "\t" << s[i] << std::endl;
  }
}

#endif

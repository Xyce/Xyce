//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 08/02/2020
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef ast_random_H
#define ast_random_H

//-------------------------------------------------------------------------------
// Random number sampled from normal distribution with
// mean μ and standard deviation (α)/n
// This uses absolute variation
// The "n" argument, AKA num_sigmas is optional
// The 4th argument, "multiplier" is not supported
// 
// Here is what the multiplier does in codes that support it:
// multiplier  - An optional argument. If you do not specify a multiplier, the default is 1. 
// If specified, code recalculates many times and saves the largest deviation. The resulting 
// parameter value might be greater than or less than nominal_val. The resulting distribution is bimodal.
//
template <typename ScalarT>
class agaussOp : public astNode<ScalarT>
{
  public:
    agaussOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args):
      astNode<ScalarT>(args), 
      value_(0.0),
      setValueCalledBefore_(false)
    {
      if (args.size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (agauss) needs at least 2 arguments.")); yyerror(errStr);
      }

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child = this->childrenAstNodes_;
      if (child.size() < 4) { child.resize(4); }

      if (args.size() > 2) { } else { child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(1.0)); }

      if (args.size() > 3)
      {
        std::vector<std::string> errStr(1,std::string("AST node (agauss) accepts at most 3 arguments.  Multiplier argument is not supported")); yyerror(errStr);
      }
      else { child[3] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(1.0)); }

      // should check to make sure that mu, alpha and n are simple constant numbers
      Teuchos::RCP<astNode<ScalarT> > & mu_ = this->childrenAstNodes_[0];
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, the random operators are all applied to real-valued .params.
    // update this if it changes.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "agauss operator id = " << this->id_ << std::endl;
      ++indent;
      this->childrenAstNodes_[0]->output(os,indent+1);
      this->childrenAstNodes_[1]->output(os,indent+1);
      this->childrenAstNodes_[2]->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "agauss operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "AGAUSS";
    }

    virtual bool agaussType()     { return true; };

    ScalarT getMu ()    { return (this->childrenAstNodes_[0]->val()); };
    ScalarT getAlpha () { return (this->childrenAstNodes_[1]->val()); };
    ScalarT getN ()     { return (this->childrenAstNodes_[2]->val()); };
    ScalarT getMult ()  { return (this->childrenAstNodes_[3]->val()); };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }

    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<agaussOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<agaussOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 

      this->childrenAstNodes_[0]->accept(visitor, this->childrenAstNodes_[0]);  // mu
      this->childrenAstNodes_[1]->accept(visitor, this->childrenAstNodes_[1]);  // alpha
      this->childrenAstNodes_[2]->accept(visitor, this->childrenAstNodes_[2]);  // n
      this->childrenAstNodes_[3]->accept(visitor, this->childrenAstNodes_[3]);  // mult
    }

  private:
    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Random number sampled from normal distribution with
// mean μ and standard deviation (α ∗ μ )/n
// This uses relative variation
// The "n" argument, AKA num_sigmas is optional
// The 4th argument, "multiplier" is not supported
template <typename ScalarT>
class gaussOp : public astNode<ScalarT>
{
  public:
    gaussOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args):
      astNode<ScalarT>(args), 
      value_(0.0),
      setValueCalledBefore_(false)
    {
      if (args.size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (gauss) needs at least 2 arguments.")); yyerror(errStr);
      }

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child = this->childrenAstNodes_;
      if (child.size() < 4) { child.resize(4); }

      if (args.size() > 2) { } else { child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(1.0)); }
      if (args.size() > 3)
      {
        std::vector<std::string> errStr(1,std::string("AST node (gauss) accepts at most 3 arguments.  Multiplier argument is not supported")); yyerror(errStr);
      }
      else { child[3] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(1.0)); }

      // should check to make sure that mu, alpha and n are simple constant numbers
      Teuchos::RCP<astNode<ScalarT> > & mu_ = this->childrenAstNodes_[0];
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, the random operators are all applied to real-valued .params.
    // update this if it changes.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "gauss operator id = " << this->id_ << std::endl;
      ++indent;
      this->childrenAstNodes_[0]->output(os,indent+1);
      this->childrenAstNodes_[1]->output(os,indent+1);
      this->childrenAstNodes_[2]->output(os,indent+2);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "gauss operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "GAUSS";
    }

    virtual bool gaussType()      { return true; };

    ScalarT getMu ()    { return (this->childrenAstNodes_[0]->val()); };
    ScalarT getAlpha () { return (this->childrenAstNodes_[1]->val()); };
    ScalarT getN ()     { return (this->childrenAstNodes_[2]->val()); };
    ScalarT getMult ()  { return (this->childrenAstNodes_[3]->val()); };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }

    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<gaussOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<gaussOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 

      this->childrenAstNodes_[0]->accept(visitor, this->childrenAstNodes_[0]);  // mu
      this->childrenAstNodes_[1]->accept(visitor, this->childrenAstNodes_[1]);  // alpha
      this->childrenAstNodes_[2]->accept(visitor, this->childrenAstNodes_[2]);  // n
      this->childrenAstNodes_[3]->accept(visitor, this->childrenAstNodes_[3]);  // mult
    }

  private:
    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Random number sampled from uniform distribution with
// mean μ and standard deviation (α)/n
// This uses absolute variation
// The "n" argument, AKA num_sigmas which is used in AGAUSS is not used for AUNIF
// The 3rd argument, "multiplier" is not supported
template <typename ScalarT>
class aunifOp : public astNode<ScalarT>
{
  public:
    aunifOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args):
      astNode<ScalarT>(args), 
      value_(0.0),
      setValueCalledBefore_(false)
    {
      if (args.size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (aunif) needs at least 2 argument.")); yyerror(errStr);
      }

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child = this->childrenAstNodes_;
      if (child.size() < 3) { child.resize(3); }

      if (args.size() > 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (aunif) accepts at most 2 arguments.  Multiplier argument is not supported")); yyerror(errStr);
      }
      else
      {
        child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(1.0));
      }

      // should check to make sure that mu, alpha are simple constant numbers
      Teuchos::RCP<astNode<ScalarT> > & mu_ = this->childrenAstNodes_[0];
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, the random operators are all applied to real-valued .params.
    // update this if it changes.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "aunif operator id = " << this->id_ << std::endl;
      ++indent;
      this->childrenAstNodes_[0]->output(os,indent+1);
      this->childrenAstNodes_[1]->output(os,indent+1);
      this->childrenAstNodes_[2]->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "aunif operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "AUNIF";
    }

    virtual bool aunifType()      { return true; };

    ScalarT getMu ()    { return (this->childrenAstNodes_[0]->val()); };
    ScalarT getAlpha () { return (this->childrenAstNodes_[1]->val()); };
    ScalarT getMult ()  { return (this->childrenAstNodes_[2]->val()); };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<aunifOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<aunifOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->childrenAstNodes_[0]->accept(visitor, this->childrenAstNodes_[0]); 
      this->childrenAstNodes_[1]->accept(visitor, this->childrenAstNodes_[1]); 
      this->childrenAstNodes_[2]->accept(visitor, this->childrenAstNodes_[2]); 
    }

  private:

    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Random number sampled from uniform distribution with
// mean μ and standard deviation (α ∗ μ )/n
// This uses relative variation
// The "n" argument, AKA num_sigmas which is used in AGAUSS is not used for UNIF
template <typename ScalarT>
class unifOp : public astNode<ScalarT>
{
  public:
    unifOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args):
      astNode<ScalarT>(args), 
      value_(0.0),
      setValueCalledBefore_(false)
    {
      if (args.size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (unif) needs at least 2 argument.")); yyerror(errStr);
      }

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child = this->childrenAstNodes_;
      if (child.size() < 3) { child.resize(3); }

      if (args.size() > 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (unif) accepts at most 2 arguments.  Multiplier argument is not supported")); yyerror(errStr);
      }
      else { child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(1.0)); }

      // should check to make sure that mu, alpha and n are simple constant numbers
      Teuchos::RCP<astNode<ScalarT> > & mu_ = this->childrenAstNodes_[0];
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, the random operators are all applied to real-valued .params.
    // update this if it changes.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unif operator id = " << this->id_ << std::endl;
      ++indent;
      this->childrenAstNodes_[0]->output(os,indent+1);
      this->childrenAstNodes_[1]->output(os,indent+1);
      this->childrenAstNodes_[2]->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "unif operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "UNIF";
    }

    virtual bool unifType()       { return true; };

    ScalarT getMu ()    { return (this->childrenAstNodes_[0]->val()); };
    ScalarT getAlpha () { return (this->childrenAstNodes_[1]->val()); };
    ScalarT getMult ()  { return (this->childrenAstNodes_[2]->val()); };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<unifOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<unifOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->childrenAstNodes_[0]->accept(visitor, this->childrenAstNodes_[0]); 
      this->childrenAstNodes_[1]->accept(visitor, this->childrenAstNodes_[1]); 
      this->childrenAstNodes_[2]->accept(visitor, this->childrenAstNodes_[2]); 
    }

  private:
    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// random number between 0 and 1 sampled from a uniform distribution
template <typename ScalarT>
class randOp : public astNode<ScalarT>
{
  public:
    randOp (): astNode<ScalarT>(),
      value_(0.0),
      setValueCalledBefore_(false)
    {
      value_ = 0.5;
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, the random operators are all applied to real-valued .params.
    // update this if it changes.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "rand operator id = " << this->id_ << std::endl;
    }

    virtual void compactOutput(std::ostream & os){ output(os,0);}

    virtual void codeGen (std::ostream & os )
    {
      // fix this
      os << "RAND";
    }

    virtual bool randType()       { return true; };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<randOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<randOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
    } // 2nd dispatch

  private:
    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Hspice version of limit
//
// This is specified in an expression as: LIMIT(nominal_val, abs_variation)
//
// If running sampling, it will return either
// (nominal_val + abs_variation) or (nominal_val - abs_variation).
//
// Which of those two values are returned depends on a random number
// between -1 and +1.  If that random number is > 0 then it returns
// (nominal_val + abs_variation).
// Otherwise it returns (nominal_val - abs_variation).
//
// If not running sampling, it returns the nominal value.
//
// In Hspice, limit has an optional 3rd argument, which specifies a
// multiplier.  This multiplier serves the same function as it does in
// every other random operator.  Currently, Xyce doesn't support it.
// Also, in Xyce a 3-argument limit will invoke the Pspice version of
// limit. (see limitOp in the ast.h file).
template <typename ScalarT>
class twoArgLimitOp : public astNode<ScalarT>
{
  public:
    twoArgLimitOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst):
      astNode<ScalarT>(xAst,yAst),
      value_(0.0)
    {
      // should check to make sure that mu, alpha and n are simple constant numbers
      Teuchos::RCP<astNode<ScalarT> > & nominal = (this->childrenAstNodes_[0]);
      Teuchos::RCP<astNode<ScalarT> > & abs_variation = (this->childrenAstNodes_[1]);
      value_ = (nominal->val());
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & nominal = (this->childrenAstNodes_[0]);
      Teuchos::RCP<astNode<ScalarT> > & abs_variation = (this->childrenAstNodes_[1]);
      return (nominal->dx (i));
    };

    // ERK check this
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs)
    {
      result = value_;
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, the random operators are all applied to real-valued .params.
    // update this if it changes.
    virtual bool getIsComplex () { return false; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "Hpsice limit operator id = " << this->id_ << std::endl;
      ++indent;
      this->childrenAstNodes_[0]->output(os,indent+1);
      this->childrenAstNodes_[1]->output(os,indent+1);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "twoArgLimit operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      /// fix this
      os << "TWO_ARG_LIMIT";
    }

    virtual bool twoArgLimitType() { return true; };

    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; }

    ScalarT getNominal()   { return (this->childrenAstNodes_[0]->val());  }
    ScalarT getVariation() { return (this->childrenAstNodes_[1]->val()); }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      Teuchos::RCP<twoArgLimitOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<twoArgLimitOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); 
      this->childrenAstNodes_[0]->accept(visitor, this->childrenAstNodes_[0]); 
      this->childrenAstNodes_[1]->accept(visitor, this->childrenAstNodes_[1]); 
    }

  private:
    ScalarT value_;
};

#endif

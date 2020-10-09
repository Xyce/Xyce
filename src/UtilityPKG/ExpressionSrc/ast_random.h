//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
//
// ERK.  Need to make the "n" parameter optional!  
template <typename ScalarT>
class agaussOp : public astNode<ScalarT>
{
  public:
    agaussOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &nAst):
      astNode<ScalarT>(xAst,yAst),
      nAst_(nAst),
      mu_(xAst),
      alpha_(yAst),
      n_(nAst_),
      value_(0.0),
      setValueCalledBefore_(false)
    {
      // should check to make sure that mu, alpha and n are simple constant numbers
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "agauss operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      nAst_->output(os,indent+1);
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

    ScalarT getMu ()    { return (this->leftAst_->val()); };
    ScalarT getAlpha () { return (this->rightAst_->val()); };
    ScalarT getN ()     { return nAst_->val() ; };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }

    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_) AST_GET_INTERESTING_OPS(nAst_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_) AST_GET_STATE_OPS(nAst_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) AST_GET_PARAM_OPS(nAst_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) AST_GET_FUNC_ARG_OPS(nAst_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) AST_GET_FUNC_OPS(nAst_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) AST_GET_VOLT_OPS(nAst_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) AST_GET_CURRENT_OPS(nAst_)
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) AST_GET_TIME_OPS(nAst_)
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > nAst_;

    Teuchos::RCP<astNode<ScalarT> > & mu_;
    Teuchos::RCP<astNode<ScalarT> > & alpha_;
    Teuchos::RCP<astNode<ScalarT> > & n_;

    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Random number sampled from normal distribution with
// mean μ and standard deviation (α ∗ μ )/n
// This uses relative variation
//
// ERK.  Need to make the "n" parameter optional!  
template <typename ScalarT>
class gaussOp : public astNode<ScalarT>
{
  public:
    gaussOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst, Teuchos::RCP<astNode<ScalarT> > &nAst):
      astNode<ScalarT>(xAst,yAst),
      nAst_(nAst),
      mu_(xAst),
      alpha_(yAst),
      n_(nAst_),
      value_(0.0),
      setValueCalledBefore_(false)
    {
      // should check to make sure that mu, alpha and n are simple constant numbers
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "gauss operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
      nAst_->output(os,indent+1);
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

    ScalarT getMu ()    { return (this->leftAst_->val()); };
    ScalarT getAlpha () { return (this->rightAst_->val()); };
    ScalarT getN ()     { return nAst_->val() ; };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }

    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_) AST_GET_INTERESTING_OPS(nAst_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_) AST_GET_STATE_OPS(nAst_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) AST_GET_PARAM_OPS(nAst_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) AST_GET_FUNC_ARG_OPS(nAst_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) AST_GET_FUNC_OPS(nAst_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) AST_GET_VOLT_OPS(nAst_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) AST_GET_CURRENT_OPS(nAst_)
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) AST_GET_TIME_OPS(nAst_)
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > nAst_;

    Teuchos::RCP<astNode<ScalarT> > & mu_;
    Teuchos::RCP<astNode<ScalarT> > & alpha_;
    Teuchos::RCP<astNode<ScalarT> > & n_;

    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Random number sampled from uniform distribution with
// mean μ and standard deviation (α)/n
// This uses absolute variation
//
// ERK.  Need to make the "n" parameter optional!    (currently excluded)
template <typename ScalarT>
class aunifOp : public astNode<ScalarT>
{
  public:
    aunifOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst):
      astNode<ScalarT>(xAst,yAst),
      mu_(xAst),
      alpha_(yAst),
      value_(0.0),
      setValueCalledBefore_(false)
    {
      // should check to make sure that mu, alpha and n are simple constant numbers
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "aunif operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
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

    ScalarT getMu ()    { return (this->leftAst_->val()); };
    ScalarT getAlpha () { return (this->rightAst_->val()); };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_) 
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_) 
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) 
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) 
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) 
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) 
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) 
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) 
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > & mu_;
    Teuchos::RCP<astNode<ScalarT> > & alpha_;
    //Teuchos::RCP<astNode<ScalarT> > & n_;

    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Random number sampled from uniform distribution with
// mean μ and standard deviation (α ∗ μ )/n
// This uses relative variation
//
// ERK.  Need to make the "n" parameter optional!    (currently excluded)
template <typename ScalarT>
class unifOp : public astNode<ScalarT>
{
  public:
    unifOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst):
      astNode<ScalarT>(xAst,yAst),
      mu_(xAst),
      alpha_(yAst),
      value_(0.0),
      setValueCalledBefore_(false)
    {
      // should check to make sure that mu, alpha and n are simple constant numbers
      value_ = mu_->val();
    };

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "unif operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
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

    ScalarT getMu ()    { return (this->leftAst_->val()); };
    ScalarT getAlpha () { return (this->rightAst_->val()); };

    bool getSetValueCalledBefore() { return setValueCalledBefore_; }
    ScalarT getValue() { return value_; }
    void setValue(ScalarT val) { value_ = val; setValueCalledBefore_=true; }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_) 
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_) 
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_) 
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_) 
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_) 
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_) 
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_) 
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) 
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > & mu_;
    Teuchos::RCP<astNode<ScalarT> > & alpha_;
    //Teuchos::RCP<astNode<ScalarT> > & n_;

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
      // should check to make sure that mu, alpha and n are simple constant numbers
      value_ = 0.5;
    };    
    
    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      ScalarT ret = 0.0;
      return ret;
    };

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

  private:
    ScalarT value_;
    bool setValueCalledBefore_;
};

//-------------------------------------------------------------------------------
// Hspice version of limit, which is used to specify a probability distribution.
// The old expression library didn't support that, but returned limit(x,y)=x+y.
// That is what is implemented here.
template <typename ScalarT>
class twoArgLimitOp : public astNode<ScalarT>
{
  public:
    twoArgLimitOp (Teuchos::RCP<astNode<ScalarT> > &xAst, Teuchos::RCP<astNode<ScalarT> > &yAst):
      astNode<ScalarT>(xAst,yAst),
      value_(0.0)
    {
      // should check to make sure that mu, alpha and n are simple constant numbers
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);
      value_ = (x->val()+y->val());
    };    

    virtual ScalarT val() { return value_; };

    virtual ScalarT dx (int i)
    {
      Teuchos::RCP<astNode<ScalarT> > & x = (this->leftAst_);
      Teuchos::RCP<astNode<ScalarT> > & y = (this->rightAst_);

      return (x->dx (i) + y->dx (i));
    };

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "twoArgLimit operator id = " << this->id_ << std::endl;
      ++indent;
      this->leftAst_->output(os,indent+1);
      this->rightAst_->output(os,indent+1);
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

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS2(leftAst_) AST_GET_INTERESTING_OPS2(rightAst_)
    }

    virtual void getStateOps(stateOpVectorContainers<ScalarT> & ovc)
    {
AST_GET_STATE_OPS2(leftAst_) AST_GET_STATE_OPS2(rightAst_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(leftAst_) AST_GET_PARAM_OPS(rightAst_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(leftAst_) AST_GET_FUNC_ARG_OPS(rightAst_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(leftAst_) AST_GET_FUNC_OPS(rightAst_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(leftAst_) AST_GET_VOLT_OPS(rightAst_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(leftAst_) AST_GET_CURRENT_OPS(rightAst_)
    }

    virtual void getTimeOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & timeOpVector)
    {
AST_GET_TIME_OPS(leftAst_) AST_GET_TIME_OPS(rightAst_) 
    }

  private:
    ScalarT value_;
};

#endif

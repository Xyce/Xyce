#ifndef _astfuncs_h_
#define _astfuncs_h_

#define AST_OP_MACRO(NAME,DX)                                                          \
template <typename ScalarT>                                                            \
class NAME ## Op : public astNode<ScalarT>                                             \
{                                                                                      \
  public:                                                                              \
    NAME ## Op (Teuchos::RCP<astNode<ScalarT> > &left): astNode<ScalarT>(left) {};     \
                                                                                       \
    virtual ScalarT val() { return std::NAME(this->leftAst_->val()); }                 \
                                                                                       \
    virtual ScalarT dx(int i) { return DX; }                                           \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << #NAME << " operator " << std::endl;                                        \
      ++indent;                                                                        \
      this->leftAst_->output(os,indent+1);                                             \
    }                                                                                  \
                                                                                       \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "std::" #NAME << "(";                                                      \
      this->leftAst_->codeGen(os);                                                     \
      os << ")";                                                                       \
    }                                                                                  \
};

AST_OP_MACRO( sqrt, (this->leftAst_->dx(i)/(2.*std::sqrt(this->leftAst_->val()))))
AST_OP_MACRO( exp, (this->leftAst_->dx(i)*std::exp(this->leftAst_->val())))
AST_OP_MACRO( abs, (std::real(this->leftAst_->val()) >= 0 ? this->leftAst_->dx(i) : ScalarT(-this->leftAst_->dx(i))))
AST_OP_MACRO( sin, (this->leftAst_->dx(i)*std::cos(this->leftAst_->val())))
AST_OP_MACRO( cos, (-this->leftAst_->dx(i)*std::sin( this->leftAst_->val())))
AST_OP_MACRO( acos, ( -this->leftAst_->dx(i)/std::sqrt(1.-this->leftAst_->val()*this->leftAst_->val()) ))
AST_OP_MACRO( acosh, (this->leftAst_->dx(i)/std::sqrt((this->leftAst_->val()-1.)*(this->leftAst_->val()+1.))))
AST_OP_MACRO( asin, (this->leftAst_->dx(i)/std::sqrt(1.-this->leftAst_->val()*this->leftAst_->val())))
AST_OP_MACRO( asinh, (this->leftAst_->dx(i)/std::sqrt(1.+this->leftAst_->val()*this->leftAst_->val())))
AST_OP_MACRO( atan, (this->leftAst_->dx(i)/(1.+this->leftAst_->val()*this->leftAst_->val())))
AST_OP_MACRO( atanh, (this->leftAst_->dx(i)/(1.-this->leftAst_->val()*this->leftAst_->val())))
AST_OP_MACRO( cosh, (this->leftAst_->dx(i)*std::sinh(this->leftAst_->val())))
AST_OP_MACRO( log, (this->leftAst_->dx(i)/this->leftAst_->val()))
AST_OP_MACRO( log10, (this->leftAst_->dx(i)/(std::log(ScalarT(10))*this->leftAst_->val())))
AST_OP_MACRO( sinh, (this->leftAst_->dx(i)*std::cosh(this->leftAst_->val())))
AST_OP_MACRO( tan, (this->leftAst_->dx(i)*(1.+std::tan(this->leftAst_->val())*std::tan(this->leftAst_->val()))))
AST_OP_MACRO( tanh, (this->leftAst_->dx(i)/(std::cosh(this->leftAst_->val())*std::cosh(this->leftAst_->val()))))

#endif

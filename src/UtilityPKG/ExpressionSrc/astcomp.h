
#ifndef _astcomp_h_
#define _astcomp_h_

#define AST_CMP_OP_MACRO(NAME,FCTQUOTE,VAL,DX)                                \
template <typename ScalarT>                                                            \
class NAME : public astNode<ScalarT>                                                   \
{                                                                                      \
  public:                                                                              \
    NAME (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):                            \
        astNode<ScalarT>(left,right) {};                                               \
                                                                                       \
    virtual ScalarT val(){return VAL; }                                                \
                                                                                       \
    virtual ScalarT dx(int i) { return DX; }                                           \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << FCTQUOTE " operator " << std::endl;                                                     \
      ++indent;                                                                        \
      this->leftAst_->output(os,indent+1);                                             \
      this->rightAst_->output(os,indent+1);                                            \
    }                                                                                  \
                                                                                       \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "(";                                                                       \
      this->leftAst_->codeGen(os);                                                     \
      os << FCTQUOTE;                                                                   \
      this->rightAst_->codeGen(os);                                                    \
      os << ")";                                                                       \
    }                                                                                  \
};                          

AST_CMP_OP_MACRO(  gtOp,  ">", ((std::real(this->leftAst_->val()) > std::real(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  ltOp,  "<", ((std::real(this->leftAst_->val()) < std::real(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  neOp, "!=", (((this->leftAst_->val()) != (this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  eqOp, "==", (((this->leftAst_->val()) == (this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  geOp, ">=", ((std::real(this->leftAst_->val()) >= std::real(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  leOp, "<=", ((std::real(this->leftAst_->val()) <= std::real(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO(  orOp, "||", ((std::real(this->leftAst_->val()) || std::real(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO( andOp, "&&", ((std::real(this->leftAst_->val()) && std::real(this->rightAst_->val()))? 1 : 0), (0.0))
AST_CMP_OP_MACRO( xorOp,  "^", (((std::real(this->leftAst_->val()) > 0 && std::real(this->rightAst_->val()) <= 0) || (std::real(this->leftAst_->val()) <= 0 && std::real(this->rightAst_->val()) > 0))?1.0:0.0), (0.0))
#endif


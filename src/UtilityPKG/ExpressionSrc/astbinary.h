
#ifndef _astbinary_h_
#define _astbinary_h_

#define AST_BIN_OP_MACRO(NAME,FCTQUOTE,FCTCODE,VAL,DX)                                \
template <typename ScalarT>                                                            \
class NAME : public astNode<ScalarT>                                                   \
{                                                                                      \
  public:                                                                              \
    NAME (Teuchos::RCP<astNode<ScalarT> > &left, Teuchos::RCP<astNode<ScalarT> > &right):                \
        astNode<ScalarT>(left,right) {};                                               \
                                                                                       \
    virtual ScalarT val(){return VAL; }                                                \
                                                                                       \
    virtual ScalarT dx(int i) { return DX; }                                           \
                                                                                       \
    virtual void output(std::ostream & os, int indent=0)                               \
    {                                                                                  \
      os << std::setw(indent) << " ";                                                  \
      os << FCTQUOTE << std::endl;                                                     \
      ++indent;                                                                        \
      this->leftAst_->output(os,indent+1);                                             \
      this->rightAst_->output(os,indent+1);                                            \
    }                                                                                  \
                                                                                       \
    virtual void codeGen (std::ostream & os )                                          \
    {                                                                                  \
      os << "(";                                                                       \
      this->leftAst_->codeGen(os);                                                     \
      os << FCTCODE;                                                                   \
      this->rightAst_->codeGen(os);                                                    \
      os << ")";                                                                       \
    }                                                                                  \
                                                                                       \
};                          

AST_BIN_OP_MACRO(
    binaryAddOp,
    "binary add ",
	  "+",
	  (this->leftAst_->val() + this->rightAst_->val()),
	  (this->leftAst_->dx (i) + this->rightAst_->dx (i))
    )
AST_BIN_OP_MACRO(
    binaryMinusOp,
    "binary minus ",
	  "-",
	  (this->leftAst_->val() - this->rightAst_->val()),
	  (this->leftAst_->dx (i) - this->rightAst_->dx (i))
    )
AST_BIN_OP_MACRO(
    binaryMulOp,
    "binary multiply ",
	  "*",
	  (this->leftAst_->val() * this->rightAst_->val()),
	  (this->leftAst_->dx(i) * this->rightAst_->val() + this->rightAst_->dx(i) * this->leftAst_->val())
    )
AST_BIN_OP_MACRO(
    binaryDivOp,
    "binary division ",
	  "/",
	  (this->leftAst_->val() / this->rightAst_->val()),
	  ((this->leftAst_->dx(i) * this->rightAst_->val() - this->rightAst_->dx(i) * this->leftAst_->val()) / (this->rightAst_->val() * this->rightAst_->val()))
    )
AST_BIN_OP_MACRO(
    binaryModOp,
    "modulus operator ",
	  "%",
	  (static_cast<int>(std::real(this->leftAst_->val())) % static_cast<int>(std::real(this->rightAst_->val()))),
	  (0.0)
    )
#endif


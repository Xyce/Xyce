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
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

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


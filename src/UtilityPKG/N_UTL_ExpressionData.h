//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Class         : N_UTL_ExpressionData
//
// Purpose       : This class manages a single expression.
//
// Special Notes : This class owns a single expression class pointer,
//                 and all the auxilliary data needed to manage its
//                 usage.
//
//                 I felt that this would be a good way to avoid
//                 bloat in functions like outputPRINT.
//
//                 I made all the data public.  Right now, the
//                 evaluate function just returns the value of the
//                 expression.
//
//                 The derivative values are all in
//                 the valDerivs vector after each evaluation.  If
//                 it becomes necessary later to access derivatives
//                 for some reason, a developer only has to access
//                 the public valDerivs vector, right after an
//                 evaluate function call.
//
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/09/04
//-----------------------------------------------------------------------------

#ifndef N_UTL_ExpressionData_H
#define N_UTL_ExpressionData_H

#include <string>
#include <vector>

// trilinos includes
#include <Teuchos_RCP.hpp>

#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

namespace Xyce {
namespace Util {

class baseExpressionGroup;

class ExpressionData
{
public:
  enum State {NOT_SETUP, PARSE_FAILED, UNRESOLVED_SYMBOL, READY};

  ExpressionData (
      const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group,
      const std::string &expression);

  ~ExpressionData();

  const std::string &getExpression() const 
  {
    return expressionString_;
  }

  void getExpressionArgs(std::vector<std::string> & args);

  bool parsed() const;

  bool getSensitivitiesPossible () const { return sensitivitiesPossible_; }

  State setup(
    Parallel::Machine                   comm,
    const Util::Op::BuilderManager &    op_builder_manager,
    const Util::ParamMap &              context_function_map,
    const Util::ParamMap &              context_param_map,
    const Util::ParamMap &              context_global_param_map);

  void evaluate(
      Parallel::Machine comm, 
      double current_circuit_time, 
      double current_circuit_dt, 
      const Util::Op::OpData & opData,
      double &result
      ) const;

  void evaluate(
      Parallel::Machine comm, 
      double current_circuit_time, 
      double current_circuit_dt, 
      const Util::Op::OpData & opData,
      double &result, 
      std::vector< double > &derivs 
      ) const;

  // complex versions
  void evaluate(
      Parallel::Machine comm, 
      double current_circuit_time, 
      double current_circuit_dt, 
      const Util::Op::OpData & opData,
      std::complex<double>  &result
      ) const;

  void evaluate(
      Parallel::Machine comm, 
      double current_circuit_time, 
      double current_circuit_dt, 
      const Util::Op::OpData & opData,
      std::complex<double> &result, 
      std::vector< std::complex<double> > &derivs 
      ) const;

   bool getIsComplex () const;

private:
  Expression *                  expression_;                    ///< Compiled Expression
  std::string                   expressionString_;              ///< Original expression string
  State                         state_;                         ///< Parse state
  Op::OpList                    expressionOps_;                 ///< Ops to compute variables
  mutable std::vector<double>   variableValues_;                ///< Cache of computed variables
  bool                          sensitivitiesPossible_;

  Teuchos::RCP<Xyce::Util::baseExpressionGroup> expressionGroup_; ///< required for setting up expressions
};

} // namespace Util
} // namespace Xyce

#endif  // N_UTL_ExpressionData_H

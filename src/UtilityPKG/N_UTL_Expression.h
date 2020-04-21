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
// Creator        : Eric Keiter
//
// Creation Date  : 04/17/08
//
//-----------------------------------------------------------------------------

#ifndef N_UTL_Expression_H
#define N_UTL_Expression_H

// ---------- Standard Includes ----------
#include <vector>
#include <string>
#include <list>
#include <iosfwd>

#include <N_UTL_fwd.h>
#include <N_UTL_Interface_Enum_Types.h>
#include <N_UTL_ExpressionSymbolTable.h>

#include <N_UTL_NoCase.h>

#include <Teuchos_RCP.hpp>


namespace Xyce {
namespace Util {

class Param;
typedef unordered_map<std::string, Param, Xyce::HashNoCase, Xyce::EqualNoCase> ParamMap;

class mainXyceExpressionGroup;
class newExpression;
class ExpressionInternals;
class baseExpressionGroup;

//-----------------------------------------------------------------------------
// Class         : Expression
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
class Expression
{
  friend mainXyceExpressionGroup;

public:
  Expression (
      const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & baseGrp_,
      std::string const & exp = std::string(), 
      const std::vector<std::string> & functionArgStringVec = std::vector<std::string>(), 
      bool useNew=true);

  Expression (const Expression &);
#ifdef NEW_EXPRESSION
  Expression& operator=(const Expression& right) ; 
#endif
  ~Expression (void);

  bool parsed() const;

  // ERK new expression stuff
  int getFuncSize();
  void getFuncNames (std::vector<std::string> & funcNames);
  void getFuncPrototypeArgStrings(std::vector<std::string> & arguments);
  void attachFunctionNode (std::string & funcName, Expression & exp); 
  void attachParameterNode (std::string & paramName, Expression & exp); 

  const std::vector<std::string> & getFunctionArgStringVec ();

  // ERK old expressionstuff.  Many of these need to be removed.
  // The entire "var" approach, in which the stuff being set comes 
  // thru a single double-precision vector is bad.
  //
  // It makes the code needlessly complex and confusing, for what is fundamentally
  // a very simple operation.
  //
  bool set (std::string const & exp);
  void getSymbolTable (std::vector< ExpressionSymbolTableEntry > & names) const;
  void get_names (int const & type, std::vector< std::string > & names) const;
  int get_type (std::string const & var);
  bool make_constant (std::string const & var, double const & val);
  bool make_var (std::string const & var);


  // ERK new expression stuff:
  void getParams         (std::vector<std::string> & params) const;
  void getVoltageNodes   (std::vector<std::string> & nodes) const;
  void getDeviceCurrents (std::vector<std::string> & devices) const;
  void getLeadCurrents   (std::vector<std::string> & leads) const;
  void getFunctions      (std::vector<std::string> & funcs) const;

  int differentiate();

  // ERK new expression stuff, where external dependencies are set in separate 
  // calls.  This is better than set_var/s, but not best.  Best would be to
  // have an all-powerful group that can "pull" these values, rather than
  // have them be set.  These are appropriate only if I continue to use the 
  // lightweight group, and that is not my long term plan.
  //
  bool setParam         (const std::string &, const double &);
  bool setVoltageNode   (const std::string &, const double &);
  bool setDeviceCurrent (const std::string &, const double &);
  bool setLeadCurrent   (const std::string &, const double &);

  std::string get_expression (void) const;
  std::string get_derivative(std::string const & var);
  int get_num(int const & type);

  int evaluate (double &result, std::vector< double > &derivs);
  int evaluateFunction (double &result);

  bool set_sim_time (double time);
  bool set_sim_dt (double dt);
  bool set_temp (double const & temp);
  bool set_sim_freq (double dt);
  void set_accepted_time (double const time);

  // ERK these two functions should also be removed and replaced with
  // a traditional "getBreakPoints" function call, which is used every
  // where else in Xyce.
  double get_break_time (void);
  double get_break_time_i (void);

  const std::string & get_input (void);

  // ERK.  Many of the functions below, here need to go.  Probably not all, but many.
  // order_names and replace_func will definitely be gone.
  int order_names (std::vector< std::string > const & new_names);
  int replace_func (std::string const & func_name, Expression & func_def, int numArgs);
  bool replace_name (const std::string & old_name, const std::string & new_name);
  int replace_var(std::string const & var_name, const Expression & subexpr);
  int replace_var(const std::string & var_name, Op::Operator *op);
  int getNumDdt();
  void getDdtVals (std::vector<double> &);
  void setDdtDerivs (std::vector<double> &);
  int num_vars() const;


  // ERK.  these are fine
  bool isTimeDependent() const;
  bool isRandomDependent() const;
  void dumpParseTree();

  static void seedRandom(long seed);

  //Teuchos::RCP<Xyce::Util::newExpression> getNewExpRCP() { return newExpPtr_; }

private:

  bool useNewExpressionLibrary_;
  bool namesSet_;
  Teuchos::RCP<Xyce::Util::newExpression> newExpPtr_;
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp_;
  ExpressionInternals *expPtr_;

};

} // namespace Util
} // namespace Xyce

#endif // N_UTL_EXPRESSION_H


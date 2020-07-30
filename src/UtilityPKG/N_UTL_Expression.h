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
#include <expressionParamTypes.h>

#include <N_UTL_NoCase.h>
#include <N_UTL_BreakPoint.h>

#include <Teuchos_RCP.hpp>

#include <N_ANP_SweepParam.h>

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
      const std::vector<std::string> & functionArgStringVec = std::vector<std::string>());

  Expression (const Expression &);
  Expression& operator=(const Expression& right) ; 
  ~Expression (void);

  bool parsed() const;

  int getFuncSize();
  void getFuncNames (std::vector<std::string> & funcNames);
  void getFuncPrototypeArgStrings(std::vector<std::string> & arguments);
  void attachFunctionNode (const std::string & funcName, const Expression & exp); 
  void attachParameterNode (const std::string & paramName, const Expression & exp, enumParamType type=DOT_GLOBAL_PARAM); 

  const std::vector<std::string> & getFunctionArgStringVec ();

  // ERK some old expressionstuff.  Many of these need to be removed.
  int get_type (std::string const & var);
  bool make_constant (std::string const & var, double const & val, enumParamType type=DOT_GLOBAL_PARAM);
  bool make_var (std::string const & var, enumParamType type=DOT_GLOBAL_PARAM);

  // ERK new expression stuff.  These kind of replace "get_names"
  void getParams              (std::vector<std::string> & params) const;
  void getUnresolvedParams    (std::vector<std::string> & params) const;
  void getVoltageNodes        (std::vector<std::string> & nodes) const;
  void getDeviceCurrents      (std::vector<std::string> & devices) const;
  void getLeadCurrents        (std::vector<std::string> & leads) const;
  void getLeadCurrentsExcludeBsrc (std::vector<std::string> & leads) const;
  void getFunctions           (std::vector<std::string> & funcs) const;
  void getUnresolvedFunctions (std::vector<std::string> & funcs) const;
  void getSpecials            (std::vector<std::string> & specials) const;
  void getVariables           (std::vector<std::string> & variables) const;
  void getPowerCalcs          (std::vector<std::string> & powerCalcs) const;

  bool getIsConstant ();

  bool setTemperature   (const double & temp);

  std::string get_expression (void) const;

  int evaluate (std::complex<double> &result, std::vector< std::complex<double> > &derivs);
  int evaluateFunction (std::complex<double> &result, bool efficiencyOn=false);

  int evaluate (double &result, std::vector< double > &derivs);
  int evaluateFunction (double &result, bool efficiencyOn=false);

  bool getBreakPoints(std::vector<Util::BreakPoint> &breakPointTimes);

  const std::string & get_input (void) const;

  bool replace_name (const std::string & old_name, const std::string & new_name);  // this is for voltage names

  bool isTimeDependent() const;
  bool isRandomDependent() const;
  void dumpParseTree();

  static void seedRandom(long seed);

  void treatAsTempAndConvert();

  void processSuccessfulTimeStep ();

  // ddt information.  This is for Bsrc support of ddt.
  int getNumDdt();
  void getDdtVals (std::vector<double> & vals);
  void setDdtDerivs (std::vector<double> & vals);

  // random operator information
  void getAgaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void getGaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void getAunifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void getUnifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void getRandData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void getLimitData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);

private:

  Teuchos::RCP<Xyce::Util::newExpression> newExpPtr_;
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp_;
};

} // namespace Util
} // namespace Xyce

#endif // N_UTL_EXPRESSION_H


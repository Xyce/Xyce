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
  void replaceParameterNode (const std::string & paramName, const Expression & exp);
  void multiplyByExternalExpression (const Expression & exp); 

  const std::vector<std::string> & getFunctionArgStringVec ();

  bool make_constant (std::string const & var, double const & val, enumParamType type=DOT_GLOBAL_PARAM);

  bool make_constant (std::string const & var, std::complex<double> const & val, enumParamType type=DOT_GLOBAL_PARAM);

  void setAsGlobal();
  void setValue(double val);
  void setValue(std::complex<double> val);

  void setGroup( Teuchos::RCP<baseExpressionGroup> & grp );
  Teuchos::RCP<baseExpressionGroup> getGroup();

  // ERK new expression stuff.  These kind of replace "get_names"
  void getUnresolvedParams        (std::vector<std::string> & params) const;
  void getVoltageNodes            (std::vector<std::string> & nodes) const;
  void getDeviceCurrents          (std::vector<std::string> & devices) const;
  void getLeadCurrents            (std::vector<std::string> & leads) const;
  void getLeadCurrentsExcludeBsrc (std::vector<std::string> & leads) const;
  void getUnresolvedFunctions     (std::vector<std::string> & funcs) const;
  void getShallowSpecials         (std::vector<std::string> & specials) const;
  void getSpecials                (std::vector<std::string> & specials) const;
  void getVariables               (std::vector<std::string> & variables) const;
  void getPowerCalcs              (std::vector<std::string> & powerCalcs) const;

  const std::vector<std::string> & getVoltageNodes () const;
  const std::vector<std::string> & getDeviceCurrents () const;

  const std::vector<std::string> & getUnresolvedParams () const;
  const std::vector<std::string> & getUnresolvedFunctions () const;
  const std::vector<std::string> & getVariables() const;
  const std::vector<std::string> & getLeadCurrents        () const;
  const std::vector<std::string> & getLeadCurrentsExcludeBsrc () const;

  bool getVariableDependent() const;
  void setVariableDependent(bool val);
  bool getVoltageNodeDependent() const;
  bool getDeviceCurrentDependent() const;
  bool getLeadCurrentDependent() const;
  bool getLeadCurrentDependentExcludeBsrc() const;
  bool getSpecialsDependent() const;

  bool getScheduleDependent() const;

  bool getIsConstant ();

  bool getIsComplex ();

  bool setTemperature   (const double & temp);

  std::string get_expression (void) const;
  std::string get_original_expression (void) const;

  void generateExpressionString (std::string & expStr);

  bool updateForStep ();

  bool evaluate (std::complex<double> &result, std::vector< std::complex<double> > &derivs);
  bool evaluateFunction (std::complex<double> &result, bool efficiencyOn=false);

  bool evaluate (double &result, std::vector< double > &derivs);
  bool evaluateFunction (double &result, bool efficiencyOn=false);

  void clearOldResult ();

  bool getBreakPoints(std::vector<Util::BreakPoint> &breakPointTimes);

  const std::string & get_input (void) const;

  bool replace_name (const std::string & old_name, const std::string & new_name);  // this is for voltage names
  bool replaceParameterName (const std::string & old_name, const std::string & new_name);  // this is for parameter names  

  bool isTimeDependent() const;
  bool isFreqDependent() const;
  bool isSolutionDependent() const;
  bool isRandomDependent() const;
  bool isOriginalRandomDependent() const;
  void dumpParseTree();

  static void seedRandom(long seed);

  void treatAsTempAndConvert();

  static void clearProcessSuccessfulTimeStepMap ();
  void processSuccessfulTimeStep ();

  // ddt information.  This is for Bsrc support of ddt.
  int getNumDdt();
  void getDdtVals (std::vector<double> & vals);
  void setDdtDerivs (std::vector<double> & vals);

  // random operator information
  void getAgaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void setAgaussValue(int index, double value);
  void setAgaussValue(int index, std::complex<double> value);

  void getGaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void setGaussValue(int index, double value);
  void setGaussValue(int index, std::complex<double> value);

  void getAunifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void setAunifValue(int index, double value);
  void setAunifValue(int index, std::complex<double> value);

  void getUnifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void setUnifValue(int index, double value);
  void setUnifValue(int index, std::complex<double> value);

  void getRandData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void setRandValue(int index, double value);
  void setRandValue(int index, std::complex<double> value);

  void getLimitData(std::vector<Xyce::Analysis::SweepParam> & sampleVec);
  void setLimitValue(int index, double value);
  void setLimitValue(int index, std::complex<double> value);

private:

  Teuchos::RCP<Xyce::Util::newExpression> newExpPtr_;
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> grp_;

#ifdef USE_TYPE_DOUBLE
  double result_;
  std::vector<double> derivs_;
#else
  std::complex<double> result_;
  std::vector<std::complex<double> > derivs_;
#endif
};

} // namespace Util
} // namespace Xyce

#endif // N_UTL_EXPRESSION_H


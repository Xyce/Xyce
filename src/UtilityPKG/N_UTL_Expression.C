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

//-------------------------------------------------------------------------
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Eric R. Keiter, SNL
//
// Creation Date : 04/17/08
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_DEV_Const.h>

// ---------- Standard Includes ----------

#include <iterator>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <sstream> 
// ----------   Xyce Includes   ----------
#include <N_UTL_Expression.h>

#include <newExpression.h>
#include <mainXyceExpressionGroup.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : Expression::Expression
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
Expression::Expression( 
    const Teuchos::RCP<Xyce::Util::baseExpressionGroup> & baseGrp_,
    const std::string & exp, 
    const std::vector<std::string> & functionArgStringVec)
  :
   newExpPtr_(NULL),
   grp_(baseGrp_),
#ifdef USE_TYPE_DOUBLE
  result_(0.0)
#else
  result_(std::complex<double>(0.0,0.0))
#endif
{
  newExpPtr_ = Teuchos::rcp(new Xyce::Util::newExpression(exp,grp_) );

  if (!(functionArgStringVec.empty()))
  {
    newExpPtr_->setFunctionArgStringVec(functionArgStringVec);
  }

  newExpPtr_->lexAndParseExpression();
}

//-----------------------------------------------------------------------------
// Function      : Expression::Expression
// Purpose       : Copy Constructor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
Expression::Expression( const Expression & right)
  :
   newExpPtr_(right.newExpPtr_),
   grp_(right.grp_)
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::operator=
// Purpose       : assignment operator
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 1/7/2019
//-----------------------------------------------------------------------------
Expression& Expression::operator=(const Expression& right) 
{
  grp_ = right.grp_;
  newExpPtr_ = right.newExpPtr_;
  //if (*(newExpPtr_->parsed())) newExpPtr_->lexAndParseExpression();
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : Expression::~Expression
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
Expression::~Expression ()
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : Expression::parsed
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::parsed() const 
{
  return newExpPtr_->parsed();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getFuncSize
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
int Expression::getFuncSize()
{
  return newExpPtr_->getFuncOpVec().size();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getFuncNames
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::getFuncNames (std::vector<std::string> & funcNames)
{
  funcNames = newExpPtr_-> getFuncNameVec ();
}


//-----------------------------------------------------------------------------
// Function      : Expression::getFuncPrototypeArgStrings
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::getFuncPrototypeArgStrings(std::vector<std::string> & arguments)
{
  newExpPtr_->getFuncPrototypeArgStrings(arguments);
} 

//-----------------------------------------------------------------------------
// Function      : Expression::attachFunctionNode
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::attachFunctionNode (const std::string & funcName, const Expression & exp)
{
  newExpPtr_->attachFunctionNode(funcName,exp.newExpPtr_);
}

//-----------------------------------------------------------------------------
// Function      : Expression::attachParameterNode
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/9/2020
//-----------------------------------------------------------------------------
void Expression::attachParameterNode (const std::string & paramName, const Expression & exp, enumParamType type)
{
  newExpPtr_->attachParameterNode(paramName,exp.newExpPtr_, type);
}

//-----------------------------------------------------------------------------
// Function      : Expression::replaceParameterNode
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 3/15/2023
//-----------------------------------------------------------------------------
void Expression::replaceParameterNode (const std::string & paramName, const Expression & exp)
{
  newExpPtr_->replaceParameterNode(paramName,exp.newExpPtr_);
}

//-----------------------------------------------------------------------------
// Function      : Expression::multiplyByExternalExpression
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/18/2022
//-----------------------------------------------------------------------------
void Expression::multiplyByExternalExpression(const Expression & exp)
{
  newExpPtr_->multiplyByExternalExpression(exp.newExpPtr_);
}

//-----------------------------------------------------------------------------
// Function      : Expression::getFunctionArgStringVec
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getFunctionArgStringVec ()
{
  return newExpPtr_->getFunctionArgStringVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::make_constant
// Purpose       : Convert a 'string' placeholder into a constant
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::make_constant (const std::string & var, const double & val, enumParamType type)
{
  return newExpPtr_->make_constant (var,val, type);
}


//-----------------------------------------------------------------------------
// Function      : Expression::make_constant
// Purpose       : Convert a 'string' placeholder into a complex constant
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/15/2023
//-----------------------------------------------------------------------------
bool Expression::make_constant (const std::string & var, const std::complex<double> & val, enumParamType type)
{
  return newExpPtr_->make_constant (var,val, type);
}


//-----------------------------------------------------------------------------
// Function      : Expression::setAsGlobal
//
// Purpose       : Add extra layer to the AST to make it easier to handle 
//                 as a global param.
//
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/20/2021
//-----------------------------------------------------------------------------
void Expression::setAsGlobal ()
{ 
  return newExpPtr_->setAsGlobal();
}

//-----------------------------------------------------------------------------
// Function      : Expression::setValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/20/2021
//-----------------------------------------------------------------------------
void Expression::setValue(double val)
{ 
  return newExpPtr_->setValue(val);
}

//-----------------------------------------------------------------------------
// Function      : Expression::setValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/20/2021
//-----------------------------------------------------------------------------
void Expression::setValue(std::complex<double> val)
{ 
  return newExpPtr_->setValue(val);
}

//-----------------------------------------------------------------------------
// Function      : Expression::setGroup
// Purpose       : 
// Special Notes :
// Scope         : 
// Creator       : Eric Keiter, SNL
// Creation Date : 8/22/2020
//-----------------------------------------------------------------------------
void Expression::setGroup( Teuchos::RCP<baseExpressionGroup> & grp )
{
  newExpPtr_->setGroup(grp);
}

//-----------------------------------------------------------------------------
// Function      : Expression::getGroup
// Purpose       : 
// Special Notes :
// Scope         : 
// Creator       : Eric Keiter, SNL
// Creation Date : 3/2/2021
//-----------------------------------------------------------------------------
Teuchos::RCP<baseExpressionGroup> Expression::getGroup()
{
  return newExpPtr_->getGroup();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedParams
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 4/20/2020
//-----------------------------------------------------------------------------
void Expression::getUnresolvedParams (std::vector<std::string> & params) const
{
  newExpPtr_->setupVariousAstArrays();
  params.clear();
  std::vector<std::string> & unresolvedParamNameVec = newExpPtr_->getUnresolvedParamNameVec();
  if (!(unresolvedParamNameVec.empty())) { params.insert(params.end(),unresolvedParamNameVec.begin(), unresolvedParamNameVec.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedParams
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 9/9/2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getUnresolvedParams () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getUnresolvedParamNameVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getVoltageNodes   (std::vector<std::string> & nodes) const
{
  newExpPtr_->setupVariousAstArrays();

  nodes.clear();
  std::vector<std::string> & voltNames = newExpPtr_->getVoltNameVec ();
  if (!(voltNames.empty())) { nodes.insert(nodes.end(),voltNames.begin(), voltNames.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVoltageNodes
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getVoltageNodes () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getVoltNameVec ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDeviceCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getDeviceCurrents (std::vector<std::string> & devices) const
{
  newExpPtr_->setupVariousAstArrays();

  devices.clear();
  std::vector<std::string> & currentNames = newExpPtr_->getCurrentNameVec ();
  if (!(currentNames.empty())) { devices.insert(devices.end(),currentNames.begin(), currentNames.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDeviceCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getDeviceCurrents () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getCurrentNameVec ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getLeadCurrents (std::vector<std::string> & leads) const
{
  newExpPtr_->setupVariousAstArrays();
  leads.clear();
  std::vector<std::string> & leadCurrentNameVec = newExpPtr_->getLeadCurrentNameVec ();
  if (!(leadCurrentNameVec.empty())) { leads.insert(leads.end(),leadCurrentNameVec.begin(), leadCurrentNameVec.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrents
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getLeadCurrents () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getLeadCurrentNameVec ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrentsExcludeBsrc
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getLeadCurrentsExcludeBsrc (std::vector<std::string> & leads) const
{
  newExpPtr_->setupVariousAstArrays(); 
  leads.clear();
  std::vector<std::string> & leadCurrentExcludeBsrcNameVec = newExpPtr_->getLeadCurrentExcludeBsrcNameVec ();
  if (!(leadCurrentExcludeBsrcNameVec.empty())) { leads.insert(leads.end(),leadCurrentExcludeBsrcNameVec.begin(), leadCurrentExcludeBsrcNameVec.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrentsExcludeBsrc
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getLeadCurrentsExcludeBsrc () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getLeadCurrentExcludeBsrcNameVec ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedFunctions
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/23/2020
//-----------------------------------------------------------------------------
void Expression::getUnresolvedFunctions (std::vector<std::string> & funcs) const
{
  newExpPtr_->setupVariousAstArrays();

  funcs.clear();
  std::vector<std::string> & unresolvedFuncNameVec = newExpPtr_->getUnresolvedFuncNameVec();
  if (!(unresolvedFuncNameVec.empty())) { funcs.insert(funcs.end(),unresolvedFuncNameVec.begin(), unresolvedFuncNameVec.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnresolvedFunctions
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 6/23/2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getUnresolvedFunctions () const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getUnresolvedFuncNameVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getSpecials
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getSpecials (std::vector<std::string> & specials) const
{
  newExpPtr_->setupVariousAstArrays();

  specials.clear();
  if (newExpPtr_->getTimeDependent()) { specials.push_back(std::string("TIME")); }
  if (newExpPtr_->getTempDependent()) { specials.push_back(std::string("TEMP")); }
  if (newExpPtr_->getVTDependent()) { specials.push_back(std::string("VT")); }
  if (newExpPtr_->getFreqDependent()) { specials.push_back(std::string("FREQ")); }
  if (newExpPtr_->getGminDependent()) { specials.push_back(std::string("GMIN")); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getShallowSpecials
// Purpose       : 
// Special Notes : does this need to catch GMIN as well?
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2021
//-----------------------------------------------------------------------------
void Expression::getShallowSpecials (std::vector<std::string> & specials) const
{
  specials.clear();
  if (newExpPtr_->getShallowTimeDependent()) { specials.push_back(std::string("TIME")); }
  if (newExpPtr_->getShallowTempDependent()) { specials.push_back(std::string("TEMP")); }
  if (newExpPtr_->getShallowVTDependent()) { specials.push_back(std::string("VT")); }
  if (newExpPtr_->getShallowFreqDependent()) { specials.push_back(std::string("FREQ")); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVariables
//
// Purpose       : This function returns the names of parameters in the expression 
//                 which have an ongoing external dependency.    This basically 
//                 means .global_params, and not .params in Xyce. 
//
//                 The param class has several booleans in it, that pertain to 
//                 different aspects of external resolution.  But, most of them 
//                 are unrelated to whether or not they were originally a .param 
//                 or a .global_param.  So, I created another boolean to indicate 
//                 when or not a parameter is a .param.  If it is, then it can't 
//                 be a "variable".  If it is not, then it is considered a "variable".
//
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getVariables (std::vector<std::string> & variables) const
{
  newExpPtr_->setupVariousAstArrays();

  variables.clear();
  std::vector<std::string> & globalParamNameVec = newExpPtr_->getGlobalParamNameVec();
  if (!(globalParamNameVec.empty())) { variables.insert(variables.end(),globalParamNameVec.begin(), globalParamNameVec.end()); }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVariables
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
const std::vector<std::string> & Expression::getVariables() const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getGlobalParamNameVec();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getPowerCalcs
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
void Expression::getPowerCalcs       (std::vector<std::string> & powerCalcs) const
{
  newExpPtr_->setupVariousAstArrays();

  powerCalcs.clear();
  for (int ii=0;ii<newExpPtr_->getPowerOpVec().size();ii++)
  {
    std::string tmpName = newExpPtr_->getPowerOpVec()[ii]->getName();
    std::vector<std::string>::iterator it = std::find(powerCalcs.begin(), powerCalcs.end(), tmpName);
    if (it == powerCalcs.end())
    {
      powerCalcs.push_back( tmpName );
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVariableDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getVariableDependent() 
{
  return newExpPtr_->getVariableDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getVoltageNodeDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getVoltageNodeDependent() 
{
  return newExpPtr_->getVoltageNodeDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDeviceCurrentDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getDeviceCurrentDependent() 
{
  return newExpPtr_->getDeviceCurrentDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrentDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getLeadCurrentDependent() 
{
  return newExpPtr_->getLeadCurrentDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLeadCurrentDependentExcludeBsrc
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getLeadCurrentDependentExcludeBsrc() 
{
  return newExpPtr_->getLeadCurrentDependentExcludeBsrc();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getSpecialsDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getSpecialsDependent() 
{
  bool retval =   
    (newExpPtr_->getTimeDependent()) ||
    (newExpPtr_->getTempDependent()) ||
    (newExpPtr_->getVTDependent()) ||
    (newExpPtr_->getFreqDependent()) ||
    (newExpPtr_->getGminDependent());

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : Expression::getScheduleDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2022
//-----------------------------------------------------------------------------
bool Expression::getScheduleDependent() const
{
  newExpPtr_->setupVariousAstArrays();
  return newExpPtr_->getScheduleDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getIsConstant
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2020
//-----------------------------------------------------------------------------
bool Expression::getIsConstant ()
{
  return newExpPtr_->getIsConstant();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getIsComplex
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/14/2023
//-----------------------------------------------------------------------------
bool Expression::getIsComplex ()
{
  return newExpPtr_->getIsComplex ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::setTemperature
// Purpose       :
// Special Notes : This is ONLY called when you want to override the
//                 'circuit' temperature.
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/19/2020
//-----------------------------------------------------------------------------
bool Expression::setTemperature   (const double & temp)
{
  return newExpPtr_->setTemperature(temp);
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_expression
// Purpose       : Returns the most up-to-date string of the expression,
//                 including modifications.
//
// Special Notes : This was originally for debugging, but it now is used 
//                 in various non-debugging ways in the netlist parser.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
std::string Expression::get_expression () const
{
  return newExpPtr_->getExpressionString();
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_original_expression
// Purpose       : Returns the original string of the expression, without modifications.
//
// Special Notes : mostly for debugging.
//
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
std::string Expression::get_original_expression () const
{
  return newExpPtr_->getOriginalExpressionString();
}

//-----------------------------------------------------------------------------
// Function      : Expression::update
// Purpose       : Update expression for .STEP, etc.
//
// Special Notes : This is for efficiency, so that some aspects of updating
//                 an expression don't have to happen during "evaluate" or
//                 "evaluateFunction" calls, which have to happen every Newton step.
//                 Some updates only happen at the beginning of .STEP
//                 iterations.
//
//                 This returns a "true" if anything was meaningfully
//                 updated, otherwise false.
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/10/2023
//-----------------------------------------------------------------------------
bool Expression::updateForStep ()
{
  return newExpPtr_->updateForStep();
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluate
// Purpose       : Evaluate expression and derivatives using stored input values
// Special Notes : 
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluate ( std::complex<double> & exp_r, std::vector< std::complex<double> > & deriv_r)
{
  bool retVal=true;
#ifdef USE_TYPE_DOUBLE
  retVal = newExpPtr_->evaluate( result_, derivs_ );
  exp_r = std::complex<double>(result_,0.0);
  if (derivs_.size() != deriv_r.size()) {deriv_r.clear(); deriv_r.resize(derivs_.size());}
  for(int ii=0;ii<derivs_.size();ii++) { deriv_r[ii] = std::complex<double>(derivs_[ii],0.0); } 
#else
  retVal = newExpPtr_->evaluate( exp_r, deriv_r );
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluateFunction
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluateFunction ( std::complex<double> & exp_r, bool efficiencyOn )
{
  bool retVal=true; 
#ifdef USE_TYPE_DOUBLE
  retVal = newExpPtr_->evaluateFunction( result_, efficiencyOn );
  exp_r = std::complex<double>(result_,0.0);
#else
  retVal = newExpPtr_->evaluateFunction ( exp_r, efficiencyOn );
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluate
// Purpose       : Evaluate expression and derivatives using stored input values
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluate ( double & exp_r, std::vector<double> & deriv_r)
{
  bool retVal=true;
#ifdef USE_TYPE_DOUBLE
  retVal = newExpPtr_->evaluate( exp_r, deriv_r );
#else
  retVal = newExpPtr_->evaluate( result_, derivs_ );
  exp_r = std::real(result_);
  if (derivs_.size() != deriv_r.size()) {deriv_r.clear(); deriv_r.resize(derivs_.size());}
  for(int ii=0;ii<derivs_.size();ii++) {  deriv_r[ii] = std::real(derivs_[ii]); }
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::evaluateFunction
// Purpose       : Evaluate expression using stored input values
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::evaluateFunction ( double & exp_r, bool efficiencyOn )
{
  bool retVal=true; 
#ifdef USE_TYPE_DOUBLE
  retVal = newExpPtr_->evaluateFunction ( exp_r, efficiencyOn );
#else
  retVal = newExpPtr_->evaluateFunction ( result_, efficiencyOn );
  exp_r = std::real(result_);
#endif
  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Expression::clearOldResult
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 09/04/2020
//-----------------------------------------------------------------------------
void Expression::clearOldResult ()
{
  return newExpPtr_->clearOldResult();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getBreakPoints
// Purpose       : Returns next breakpoint time
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::getBreakPoints(std::vector<Util::BreakPoint> &breakPointTimes)
{
  return newExpPtr_->getBreakPoints(breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : Expression::get_input
// Purpose       : Return expression input string
// Special Notes :
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
const std::string & Expression::get_input () const
{
  return newExpPtr_->getExpressionString();
}

//-----------------------------------------------------------------------------
// Function      : Expression::replace_name
// Purpose       : Change the name of an input quantity
//
// Special Notes : ERK.  This called from the N_IO_DistToolBase.C file/class
//                 in the function DistToolBase::instantiateDevice.
//
//                 "Input quantity" in this case means voltage nodes (XEXP_NODE), 
//                 device instances (XEXP_INSTANCE), and lead currents (XEXP_LEAD).
//
//                 Sometimes, they are specified in expressions without their 
//                 names being fully resolved.  ie, the expression is inside of 
//                 a subcircuit, and thus implicitly assumes the full prefix.
//
//                 So, this function adds the full prefix to these names, 
//                 so they can be fully resolved.
//
//                 The function DistToolBase::instantiateDevice calls this 
//                 function twice for some reason that I don't (yet) understand.
//                 It seems to require 2 passes to properly update the name. ie,
//                 "name" -> "; name" -> "prefix:name".
//  
// Scope         :
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/17/08
//-----------------------------------------------------------------------------
bool Expression::replace_name ( const std::string & old_name,
                                const std::string & new_name)
{
  return newExpPtr_->replaceName( old_name, new_name );
}

bool Expression::replaceParameterName ( const std::string & old_name,
                                        const std::string & new_name)
{
  return newExpPtr_->replaceParameterName( old_name, new_name );
}

//-----------------------------------------------------------------------------
// Function      : Expression::isTimeDependent
// Purpose       : Return true if expression is either explicitly OR implicitly
//                 time dependent
// Special Notes : 
// Scope         :
// Creator       : 
// Creation Date : 10/07/2013
//-----------------------------------------------------------------------------
bool Expression::isTimeDependent() const
{
  // ERK. This is important so that (for example) capacitors with expression 
  // dependent capacitance aren't calling updateTemperature over and over 
  // again.   If they think an expression dependent parameter is time 
  // dependent, then they are obligated to keep re-evaluating it.
  // So, this needs to give the right answer.
  bool explicitTimeDep = newExpPtr_->getTimeDependent();
  bool sdtDep = !(newExpPtr_->getSdtOpVec().empty());
  bool ddtDep = !(newExpPtr_->getDdtOpVec().empty());
  return (explicitTimeDep || sdtDep || ddtDep);
}

//-----------------------------------------------------------------------------
// Function      : Expression::isFreqDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isFreqDependent() const
{
  return newExpPtr_->getFreqDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::isSolutionDependent
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isSolutionDependent() const
{
  return ( !(newExpPtr_->getVoltOpVec().empty()) || 
           !(newExpPtr_->getCurrentOpVec().empty()) );
}

//-----------------------------------------------------------------------------
// Function      : Expression::isRandomDependent
// Purpose       : Return true if expression dependent on GAUSS, AGAUSS or RAND
//
// Special Notes : This is only based on local dependence, from parsing as well 
//                 as nodes that have been "replaced".  This excludes parts of 
//                 the tree that have been "attached".  
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isRandomDependent() const
{
  return newExpPtr_->getIsShallowRandomDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::isOriginalRandomDependent
//
// Purpose       : Return true if expression dependent on GAUSS, AGAUSS or RAND
//
// Special Notes : This is only based on local dependence, from parsing. It 
//                 excludes nodes that have been replaced, and also excludes 
//                 nodes that have been attached.
// Scope         :
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Expression::isOriginalRandomDependent() const
{
  return newExpPtr_->getIsOriginalShallowRandomDependent();
}

//-----------------------------------------------------------------------------
// Function      : Expression::dumpParseTree
// Purpose       : Dump out the parse tree for an expression
// Special Notes : 
// Scope         :
// Creator       : Tom Russo, SNL
// Creation Date : 9/9/10
//-----------------------------------------------------------------------------
void Expression::dumpParseTree()
{
  newExpPtr_->dumpParseTree(Xyce::dout());
}

//-----------------------------------------------------------------------------
// Function      : Expression::seedRandom
// Purpose       : Public method to initialize random number generator
//                 used by rand(), gauss() and agauss() functions.
// Special Notes : 
// Scope         :
// Creator       : Tom Russo
// Creation Date : 03/28/17
//-----------------------------------------------------------------------------
///
/// Public method to initialize random number generator
///
/// @author Tom Russo
/// @date 03/28/17
///
void Expression::seedRandom(long seed)
{
}

//-----------------------------------------------------------------------------
// Function      : Expression::treatAsTempAndConvert
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::treatAsTempAndConvert()
{
  newExpPtr_->treatAsTempAndConvert();
}

//-----------------------------------------------------------------------------
// Function      : Expression::clearProcessSuccessfulTimeStepMap
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::clearProcessSuccessfulTimeStepMap ()
{
  newExpression::clearProcessSuccessfulTimeStepMap ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::processSuccessfulTimeStep
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::processSuccessfulTimeStep ()
{
  newExpPtr_->processSuccessfulTimeStep ();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getNumDdt
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
int Expression::getNumDdt()
{
  return newExpPtr_->getNumDdt();
}

//-----------------------------------------------------------------------------
// Function      : Expression::getDdtVals
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getDdtVals (std::vector<double> & vals)
{
  return newExpPtr_->getDdtVals(vals);
}

//-----------------------------------------------------------------------------
// Function      : Expression::setDdtDerivs
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setDdtDerivs (std::vector<double> & vals)
{
  return newExpPtr_->setDdtDerivs(vals);
}

// random operator information
//-----------------------------------------------------------------------------
// Function      : Expression::getAgaussData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getAgaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAgauss = newExpPtr_->getLocalAgaussOpVec();
  for (int ii=0;ii<localAgauss.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<agaussOp<usedType> > agOp = Teuchos::rcp_static_cast<agaussOp<usedType> > (localAgauss[ii]);

    usedType mu    = agOp->getMu();
    usedType alpha = agOp->getAlpha();
    usedType n     = agOp->getN();

    sampling_param.opName     = "AGAUSS";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_AGAUSS;
    sampling_param.type       = "NORMAL";
    sampling_param.mean       = std::real(mu);
    sampling_param.stdDev     = std::real(alpha)/std::real(n);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setAgaussValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setAgaussValue(int index, double value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAgauss = newExpPtr_->getLocalAgaussOpVec();
  if (index < localAgauss.size() && index >= 0)
  {
    localAgauss[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setAgaussValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setAgaussValue(int index, std::complex<double> value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAgauss = newExpPtr_->getLocalAgaussOpVec();
  if (index < localAgauss.size() && index >= 0)
  {
    localAgauss[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getGaussData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getGaussData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localGauss = newExpPtr_->getLocalGaussOpVec();
  for (int ii=0;ii<localGauss.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<gaussOp<usedType> > gaOp = Teuchos::rcp_static_cast<gaussOp<usedType> > (localGauss[ii]);

    usedType mu        = gaOp->getMu();
    usedType rel_alpha = gaOp->getAlpha();
    usedType n         = gaOp->getN();
    usedType alpha     = rel_alpha*mu;

    sampling_param.opName     = "GAUSS";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_GAUSS;
    sampling_param.type       = "NORMAL";
    sampling_param.mean       = std::real(mu);
    sampling_param.stdDev     = std::real(alpha)/std::real(n);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setGaussValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setGaussValue(int index, double value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localGauss = newExpPtr_->getLocalGaussOpVec();
  if (index < localGauss.size() && index >= 0)
  {
    localGauss[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setGaussValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setGaussValue(int index, std::complex<double> value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localGauss = newExpPtr_->getLocalGaussOpVec();
  if (index < localGauss.size() && index >= 0)
  {
    localGauss[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getAunifData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getAunifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAunif = newExpPtr_->getLocalAunifOpVec();
  for (int ii=0;ii<localAunif.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<aunifOp<usedType> > auOp = Teuchos::rcp_static_cast<aunifOp<usedType> > (localAunif[ii]);

    usedType mu    = auOp->getMu();
    usedType alpha = auOp->getAlpha();

    sampling_param.opName     = "AUNIF";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_AUNIF;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = std::real(mu)-std::real(alpha);
    sampling_param.stopVal  = std::real(mu)+std::real(alpha);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setAunifValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setAunifValue(int index, double value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAunif = newExpPtr_->getLocalAunifOpVec();
  if (index < localAunif.size() && index >= 0)
  {
    localAunif[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setAunifValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setAunifValue(int index, std::complex<double> value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localAunif = newExpPtr_->getLocalAunifOpVec();
  if (index < localAunif.size() && index >= 0)
  {
    localAunif[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getUnifData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getUnifData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localUnif = newExpPtr_->getLocalUnifOpVec();
  for (int ii=0;ii<localUnif.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<unifOp<usedType> > unOp = Teuchos::rcp_static_cast<unifOp<usedType> > (localUnif[ii]);

    usedType mu        = unOp->getMu();
    usedType rel_alpha = unOp->getAlpha();
    usedType alpha     = rel_alpha*mu;

    sampling_param.opName     = "UNIF";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_UNIF;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = std::real(mu)-std::real(alpha);
    sampling_param.stopVal  = std::real(mu)+std::real(alpha);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setUnifValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setUnifValue(int index, double value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localUnif = newExpPtr_->getLocalUnifOpVec();
  if (index < localUnif.size() && index >= 0)
  {
    localUnif[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setUnifValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setUnifValue(int index, std::complex<double> value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localUnif = newExpPtr_->getLocalUnifOpVec();
  if (index < localUnif.size() && index >= 0)
  {
    localUnif[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getRandData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getRandData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localRand = newExpPtr_->getLocalRandOpVec();
  for (int ii=0;ii<localRand.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<randOp<usedType> > auOp = Teuchos::rcp_static_cast<randOp<usedType> > (localRand[ii]);

    sampling_param.opName     = "RAND";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_RAND;
    sampling_param.type       = "UNIFORM";

    sampling_param.startVal = 0.0;
    sampling_param.stopVal  = 1.0;

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setRandValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setRandValue(int index, double value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localRand = newExpPtr_->getLocalRandOpVec();
  if (index < localRand.size() && index >= 0)
  {
    localRand[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setRandValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setRandValue(int index, std::complex<double> value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localRand = newExpPtr_->getLocalRandOpVec();
  if (index < localRand.size() && index >= 0)
  {
    localRand[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::getLimitData
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::getLimitData(std::vector<Xyce::Analysis::SweepParam> & sampleVec)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localTwoArgLimit = newExpPtr_->getLocalTwoArgLimitOpVec();
  for (int ii=0;ii<localTwoArgLimit.size();ii++)
  {
    Xyce::Analysis::SweepParam sampling_param;

    Teuchos::RCP<twoArgLimitOp<usedType> > auOp = Teuchos::rcp_static_cast<twoArgLimitOp<usedType> > (localTwoArgLimit[ii]);

    sampling_param.opName     = "LIMIT";
    sampling_param.astOpIndex = ii;
    sampling_param.astType    = Util::AST_LIMIT;
    sampling_param.type       = "LIMIT";

    usedType nominal        = auOp->getNominal();
    usedType abs_variation  = auOp->getVariation();

    sampling_param.startVal = std::real(nominal)-std::real(abs_variation);
    sampling_param.stopVal  = std::real(nominal)+std::real(abs_variation);

    sampleVec.push_back(sampling_param);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setLimitValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setLimitValue(int index, double value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localTwoArgLimit = newExpPtr_->getLocalTwoArgLimitOpVec();
  if (index < localTwoArgLimit.size() && index >= 0)
  {
    localTwoArgLimit[index]->setValue(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : Expression::setLimitValue
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
void Expression::setLimitValue(int index, std::complex<double> value)
{
  std::vector<Teuchos::RCP<astNode<usedType> > > & localTwoArgLimit = newExpPtr_->getLocalTwoArgLimitOpVec();
  if (index < localTwoArgLimit.size() && index >= 0)
  {
    localTwoArgLimit[index]->setValue(value);
  }
}

} // namespace Util
} // namespace Xyce


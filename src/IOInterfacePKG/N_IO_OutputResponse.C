//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ERH_Message.h>
#include <N_IO_Op.h>
#include <N_IO_OutputResponse.h>
#include <N_LAS_Vector.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace IO {

OutputResponse::OutputResponse()
  : responseFileName_("response.out"),
    responseVarList_(),
    responseVarPtr_(0),
    responseNames_(),
    numResponseVars_(0)
{}

OutputResponse::~OutputResponse()
{
  for (Util::Op::OpList::iterator it = responseVarList_.begin(); it != responseVarList_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : OutputResponse::saveResponseVarValues
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 08/11/04
//-----------------------------------------------------------------------------
void
OutputResponse::saveResponseVarValues(
  Parallel::Machine     comm,
  double                time,
  const Linear::Vector &  solnVecPtr)
{
  // save the independant variable
  int varNumber = 0;

  responseVarPtr_->at(varNumber) = time;

  // if (outputState_.dcSweepVector_.empty())  // DNS: is this a reliable way to determine if transient?
  // {
  //   responseVarPtr_->at(varNumber) = outputState_.circuitTime_;
  // }
  // else
  // {
  //   responseVarPtr_->at(varNumber) = outputState_.dcSweepVector_[ dcLoopNumber_ ].currentVal;
  // }
  varNumber++;

  std::vector<complex> result_list;
  getValues(comm, responseVarList_, Util::Op::OpData(0, &solnVecPtr, 0, 0, 0, 0), result_list);

  for (int i = 0; i < result_list.size(); ++i)
  {
    responseVarPtr_->at(varNumber) = result_list[i].real();
    varNumber++;
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputResponse::finalizeResponseVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------
void
OutputResponse::finalizeResponseVars(
  double        time)
{
  // save the resuts of any end point variables if there are any
  if (numResponseVars_ != 0)
  {
    // save the independant variable
    int varNumber = 0;
    responseVarPtr_->at(varNumber) = time;

    // if (outputState_.dcSweepVector_.empty())  // DNS: is this a reliable way to determine if transient?
    // {
    //   responseVarPtr_->at(v  arNumber) = outputState_.circuitTime_;
    // }
    // else
    // {
    //   responseVarPtr_->at(varNumber) = outputState_.dcSweepVector_[ dcLoopNumber_ ].currentVal;
    // }
    varNumber++;
  }
}


// routines to tell the output manager which variables external programs will
// need as output.  By default we'll only remember the last timepoint
// or dc step unless asked to track all history.

//-----------------------------------------------------------------------------
// Function      : OutputResponse::registerResponseVars
// Purpose       : Create an objective from string submitted by external program
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------

bool
OutputResponse::registerResponseVars(
  Parallel::Machine             comm,
  Util::Op::BuilderManager &    op_builder_manager,
  const std::string &           objString,
  std::vector<double> *         varVectorPtr)
{
  bool result = true;
  responseVarPtr_ = varVectorPtr;

  ExtendedString sVal(objString);

  sVal.toUpper();
  if (sVal.size() < 3 || sVal[1] != '(' || sVal[sVal.size()-1] != ')') {
    Report::DevelFatal0() << "OutputResponse::registerResponseVars: response var not of format V() or I(): '" << objString << "'";
  }
  numResponseVars_++;

  Util::ParamList pList;

  Util::Param parameter;
  parameter.setTag(sVal.substr(0, 1));
  parameter.setVal(1.0);
  pList.push_back(parameter);

  parameter.setTag(sVal.substr(2, sVal.size()-3));
  parameter.setVal(0.0);
  pList.push_back(parameter);

  makeOps(comm, op_builder_manager, NetlistLocation(), pList.begin(), pList.end(), std::back_inserter(responseVarList_));

  return result;
}

void
OpCallback<void>::send(
  Parallel::Machine             comm,
  const Linear::Vector *          real_solution_vector,
  const Linear::Vector *          imaginary_solution_vector,
  const Linear::Vector *          state_vector,
  const Linear::Vector *          store_vector,
  const std::vector<double> *   sens_function,
  const std::vector<double> *   sens_dOdPDirect,
  const std::vector<double> *   sens_dOdPDirectScaled,
  const std::vector<double> *   sens_dOdPAdjoint,
  const std::vector<double> *   sens_dOdPAdjointScaled) const
{
  execute(Util::Op::getValue(comm, op_, Util::Op::OpData(0, real_solution_vector, imaginary_solution_vector,
                                                         state_vector, store_vector, 0, 0, 0, 0, 0, sens_function,
                                                         sens_dOdPDirect, sens_dOdPDirectScaled, sens_dOdPAdjoint, sens_dOdPAdjointScaled)));
                                                         // 4 additional null's passed in for lead current vects.  Need to fix.  RLS 10/2014
}

void
OutputResponse::send(
  Parallel::Machine             comm,
  const Linear::Vector *          real_solution_vector,
  const Linear::Vector *          imaginary_solution_vector,
  const Linear::Vector *          state_vector,
  const Linear::Vector *          store_vector,
  const std::vector<double> *   sens_function,
  const std::vector<double> *   sens_dOdPDirect,
  const std::vector<double> *   sens_dOdPDirectScaled,
  const std::vector<double> *   sens_dOdPAdjoint,
  const std::vector<double> *   sens_dOdPAdjointScaled)
{
  for (OpCallbackVector::const_iterator it = callbacks_.begin(), end = callbacks_.end(); it != end; ++it)
  {
    const OpCallback<void> &callback = *(*it);
    callback.send(comm, real_solution_vector, imaginary_solution_vector,
                  state_vector, store_vector, sens_function,
                  sens_dOdPDirect, sens_dOdPDirectScaled, sens_dOdPAdjoint, sens_dOdPAdjointScaled);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputResponse::setExternalNetlistParams
// Purpose       : Called from owning class to set any external parameters
//                 set by the user.  Note, these parameter lists may also
//                 request response functions that the output manger must report.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 08/07/2012
//-----------------------------------------------------------------------------
void OutputResponse::setExternalNetlistParams(
  const StringPairVector &      externalNetlistParams)
{
  // externalParams can contain section names from dakota like
  // variables 2, responses 4, derivatives 4.  If we don't find any sections tags
  // then we can assume that all the parameters are just variable names to set.
  // if we find tags, then use the "responses" section to record what response functions
  // need to be reported.

  StringPairVector *section = 0;
  for (StringPairVector::const_iterator it = externalNetlistParams.begin(), end = externalNetlistParams.end(); it != end; it++)
  {
    if ((*it).first == "variables")
      section = &variablesUsedInSimulation_;
    else if ((*it).first == "functions")
      section = &responseFunctionsRequested_;
    else if ((*it).first == "derivative_variables")
      section = &derivativeVariablesRequested_;
    else if ((*it).first == "analysis_components")
      section = &analysisComponentsRequested_;
    else if (section)
      section->push_back(*it);
  }
}

} // namespace IO
} // namespace Xyce

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
// Purpose        : Output Manager
//
// Special Notes  :
//
// Creator        : David G. Baur, Raytheon
//
// Creation Date  : 11/11/2014
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <vector>

#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_ANP_Op.h>
#include <N_ANP_AnalysisManager.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_OpBuilder.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Analysis {

struct AnalysisInitialTimeOpBuilder : public Util::Op::Builder
{
  AnalysisInitialTimeOpBuilder(const AnalysisManager &analysis_manager)
    : analysisManager_(analysis_manager)
  {}

  virtual ~AnalysisInitialTimeOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<AnalysisInitialTimeOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (param_tag == "ANALYSIS_INITIAL_TIME") {
      new_op  = new AnalysisInitialTimeOp(param_tag, analysisManager_);
      new_op->addArg(param_string);
    }

    return new_op;
  }

private:
  const AnalysisManager &       analysisManager_;
};

struct AnalysisFinalTimeOpBuilder : public Util::Op::Builder
{
  AnalysisFinalTimeOpBuilder(const AnalysisManager &analysis_manager)
    : analysisManager_(analysis_manager)
  {}

  virtual ~AnalysisFinalTimeOpBuilder()
  {}

  virtual void registerCreateFunctions(Util::Op::BuilderManager &builder_manager) const
  {
    builder_manager.addCreateFunction<AnalysisFinalTimeOp>();
  }

  virtual Util::Op::Operator *makeOp(Util::ParamList::const_iterator &it) const
  {
    Util::Op::Operator *new_op = 0;
    const std::string &param_tag = (*it).tag();
    const std::string &param_string = (*it).stringValue();

    if (param_tag == "ANALYSIS_FINAL_TIME") {
      new_op  = new AnalysisFinalTimeOp(param_tag, analysisManager_);
      new_op->addArg(param_string);
    }

    return new_op;
  }

private:
  const AnalysisManager &       analysisManager_;
};

//-----------------------------------------------------------------------------
// Function      : registerOpBuilders
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void registerOpBuilders(Util::Op::BuilderManager &builder_manager, Parallel::Machine comm, AnalysisManager &analysis_manager)
{
  builder_manager.addBuilder(new AnalysisInitialTimeOpBuilder(analysis_manager));
  builder_manager.addBuilder(new AnalysisFinalTimeOpBuilder(analysis_manager));
}

} // namespace Analysis
} // namespace Xyce

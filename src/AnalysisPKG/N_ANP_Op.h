//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#ifndef Xyce_N_ANP_Op_h
#define Xyce_N_ANP_Op_h

#include <iterator>

#include <N_ANP_fwd.h>
#include <N_UTL_Op.h>

namespace Xyce {
namespace Analysis {

class AnalysisInitialTimeOp : public Util::Op::Op<AnalysisInitialTimeOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  AnalysisInitialTimeOp(const std::string &name, const AnalysisManager &analysis_manager)
    : Base(name),
      analysisManager_(analysis_manager)
  {}

  virtual ~AnalysisInitialTimeOp()
  {}

  static complex get(const AnalysisInitialTimeOp &op, const Util::Op::OpData &op_data);

  const AnalysisManager &       analysisManager_;
};

class AnalysisFinalTimeOp : public Util::Op::Op<AnalysisFinalTimeOp, Util::Op::ReduceNone, Util::Op::EvalNoop>
{
public:
  AnalysisFinalTimeOp(const std::string &name, const AnalysisManager &analysis_manager)
    : Base(name),
      analysisManager_(analysis_manager)
  {}

  virtual ~AnalysisFinalTimeOp()
  {}

  static complex get(const AnalysisFinalTimeOp &op, const Util::Op::OpData &op_data);

  const AnalysisManager &       analysisManager_;
};


} // namespace Device
} // namespace Xyce

#endif // Xyce_N_ANP_Op_h

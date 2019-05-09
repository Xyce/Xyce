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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputResults_h
#define Xyce_N_IO_OutputResults_h

#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>

#include <N_ANP_StepEvent.h>
#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>

namespace Xyce {
namespace IO {

typedef Util::ListenerAutoSubscribe<Analysis::StepEvent> StepEventListener;
typedef std::vector<Util::ExpressionData *> ResultVector;

class OutputResults : public Analysis::ProcessorBase,
                      public StepEventListener
{
public:
  OutputResults(
    Parallel::Machine           comm,
    Analysis::AnalysisManager & analysis_manager,
    OutputMgr &                 output_manager);

  virtual ~OutputResults();

private:
  OutputResults(const OutputResults &);
  OutputResults &operator=(const OutputResults &);

public:
  bool addResultParams(const Util::OptionBlock &option_block);

  void setup(Parallel::Machine comm, OutputMgr &output_manager);

  void output(
    Parallel::Machine                   comm,
    double                              circuit_time,
    double                              circuit_dt,
    const Analysis::SweepVector &       step_sweep_vector,
    int                                 step_loop_number,
    const Linear::Vector &              solution_vector,
    const Linear::Vector &              state_vector,
    const Linear::Vector &              store_vector,
    const Linear::Vector &              lead_current_vector,
    const Linear::Vector &              junction_voltage_vector);

private:
  void notify(const Analysis::StepEvent &step_event);
  void steppingComplete();

private:
  Parallel::Machine             comm_;
  Analysis::AnalysisManager &   analysisManager_;
  OutputMgr &                   outputManager_;
  std::ostream *                os_;                            ///< ostream for .result output
  ResultVector                  resultVector_;                  ///< Expressions from .RESULT
  bool                          noIndexResult_;
};

bool
registerOutputResultsFactory(
  Analysis::FactoryBlock &      factory_block,
  Parallel::Machine             comm)
;

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_OutputResults_h

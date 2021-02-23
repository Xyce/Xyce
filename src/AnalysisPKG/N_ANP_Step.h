//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Step analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Step_h
#define Xyce_N_ANP_Step_h

#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : Step
// Purpose       : Step analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Step : public AnalysisBase
{
public:
  Step(AnalysisManager &analysis_manager, Loader::Loader &loader, AnalysisBase &child_analysis)
    : AnalysisBase(analysis_manager, "Step"),
      analysisManager_(analysis_manager),
      loader_(loader),
      outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
      childAnalysis_(child_analysis),
      stepSweepVector_(),
      stepLoopSize_(0)
  {}

  virtual ~Step()
  {}

  bool setAnalysisParams(const Util::OptionBlock & paramsBlock);
  bool setDataStatements(const Util::OptionBlock & paramsBlock);

  bool convertDataToSweepParams();

  const TimeIntg::TIAParams &getTIAParams() const; // override
  TimeIntg::TIAParams &getTIAParams(); // override

  virtual bool getDCOPFlag() const;

protected:
  virtual bool doRun();
  virtual bool doInit();
  virtual bool doLoopProcess();
  virtual bool doProcessSuccessfulStep();
  virtual bool doProcessFailedStep();
  virtual bool doFinish();
  virtual bool doHandlePredictor() { return true; }

private:
  AnalysisManager &     analysisManager_;
  Loader::Loader &      loader_;
  OutputMgrAdapter &    outputManagerAdapter_;
  AnalysisBase &        childAnalysis_;
  SweepVector           stepSweepVector_;
  int                   stepLoopSize_;

  std::map< std::string, std::vector<std::string> > dataNamesMap_;
  std::map< std::string, std::vector< std::vector<double> > > dataTablesMap_;
};

bool registerStepFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_Step_h

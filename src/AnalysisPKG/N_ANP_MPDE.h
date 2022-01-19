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
// Purpose        : MPDE analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Todd Coffey, 1414, Ting Mei 1437
//
// Creation Date  : 07/23/08
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_MPDE_h
#define Xyce_N_ANP_MPDE_h

#include <N_ANP_fwd.h>
#include <N_TOP_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_ANP_StepEvent.h>
#include <N_UTL_Listener.h>
#include <N_UTL_OptionBlock.h>

class N_MPDE_Manager;

namespace Xyce {
namespace Analysis {

typedef Util::ListenerAutoSubscribe<StepEvent> StepEventListener;

//-------------------------------------------------------------------------
// Class         : MPDE
// Purpose       : MPDE analysis class
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class MPDE : public AnalysisBase,
             public StepEventListener
{
public:
  MPDE(
    AnalysisManager &                   analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager);

  virtual ~MPDE();
  

  void notify(const StepEvent &event);

  const TimeIntg::TIAParams &getTIAParams() const; // override
  TimeIntg::TIAParams &getTIAParams(); // override

  bool getDCOPFlag() const;

protected:
  void finalExpressionBasedSetup();
  bool doRun();
  bool doInit();
  bool doLoopProcess();
  bool processSuccessfulDCOP();
  bool processFailedDCOP();
  bool doProcessSuccessfulStep();
  bool doProcessFailedStep();
  bool doFinish();
  bool doHandlePredictor();

public:
  bool finalVerboseOutput();

  N_MPDE_Manager &getMPDEManager()
  {
    return *mpdeManager_;
  }

private:
  AnalysisManager &     analysisManager_;
  Loader::Loader &      loader_;
  Linear::System &      linearSystem_;
  Nonlinear::Manager &  nonlinearManager_;
  Topo::Topology &      topology_;
  N_MPDE_Manager *      mpdeManager_;
};

bool registerMPDEFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_MPDE_h


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
// Purpose       : 
// Special Notes :
// Creator       : 
// Creation Date : 
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_ROL_TRAN_Optimizationh
#define Xyce_N_ANP_ROL_TRAN_Optimizationh

#include <Xyce_config.h>

#ifndef FD_HESSIAN
#define FD_HESSIAN 0 // 0 to set all hessvec to zero (Gauss-Newton Hessian); else FD Hessian
#endif

#ifndef IDENTITY_PRECONDITIONER
#define IDENTITY_PRECONDITIONER 1
#endif

#ifdef Xyce_ROL

namespace Xyce {
namespace Nonlinear {
template <typename ScalarT>
class objectiveFunctionData;
}
}

class ROL_Objective;

// ---------- ROL Includes ------------//

#include "ROL_StdVector.hpp"
#include "ROL_Vector.hpp"
#include "ROL_XyceVector.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <map>

// ------------ Xyce Includes -----------//

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_ROL.h>
#include <N_ANP_Transient.h>
#include <N_ANP_OutputMgrAdapter.h>

#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Solver.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Builder.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_NLS_Sensitivity.h>
#include <N_NLS_fwd.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : ROL_TRAN
// Purpose       : Thin inheritance layer of Transient class for ROL
// Special Notes :
// Creator       : Heidi Thornquist, SNL
// Creation Date : 12/13/22
//-------------------------------------------------------------------------
class ROL_TRAN: public Transient
{
public:

  ROL_TRAN(
      AnalysisManager & analysis_manager, 
      Nonlinear::Manager & nonlinear_manager,
      Loader::Loader & loader, 
      Linear::System & linear_system,
      Topo::Topology & topology,
      IO::InitialConditionsManager & initial_conditions_manager,
      IO::RestartMgr & restart_manager)
  : Transient( analysis_manager, &linear_system, nonlinear_manager, loader, topology, initial_conditions_manager, restart_manager ),
    analysisManager_( analysis_manager ),
    nonlinearManager_( nonlinear_manager ),
    loader_( loader ),
    topology_( topology ),
    initialConditionsManager_( initial_conditions_manager ),
    linearSystem_( linear_system ),
    restartManager_( restart_manager ),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    stepLoopSize_(0),
    numParams_(0)
  {}

  virtual ~ROL_TRAN() { doFree(); }
/*
  bool doAllocations(int nc, int nz);
  std::map< std::string, std::vector<std::string> >& getDataNamesMap() { return dataNamesMap_; }
  std::map< std::string, std::vector< std::vector<double> > >& getDataTablesMap() { return dataTablesMap_; }
  bool createObjectives(const std::vector<ROL_Objective>& objVec);
*/

  using Transient::setTimeIntegratorOptions;

  using Transient::doFinish;
  using Transient::doLoopProcess;
  using Transient::doProcessFailedStep;
  using Transient::doHandlePredictor;

  bool doInit();
  bool doRun();
  bool doProcessSuccessfulStep();

  Teuchos::RCP<::ROL::Objective_SimOpt<RealT> > obj_;

  std::vector<Linear::Vector *>         solutionPtrVector_;
  std::vector<Linear::Vector *>         statePtrVector_;
  std::vector<Linear::Vector *>         constraintPtrVector_;
  std::vector<Linear::Vector *>         mydfdpPtrVector_;
  std::vector<Linear::Vector *>         mydqdpPtrVector_;
  std::vector<Linear::Vector *>         mydbdpPtrVector_;
  std::vector<Linear::Vector *>         mysensRHSPtrVector_;

private:

  bool doFree() { return true; }

  std::vector<ROL_Objective>            objVec_;

  AnalysisManager &                     analysisManager_;
  Nonlinear::Manager &                  nonlinearManager_; // TT
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  Linear::System &                      linearSystem_;
  IO::RestartMgr &                      restartManager_; 
  OutputMgrAdapter &                    outputManagerAdapter_;
  int                                   stepLoopSize_;
  int                                   numParams_;
};


} // namespace Analysis
} // namespace Xyce`

#endif // Xyce_ROL

#endif // Xyce_N_ANP_ROL_TRAN_Optimization

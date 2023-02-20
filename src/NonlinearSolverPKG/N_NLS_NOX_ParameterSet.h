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
// Purpose        : Specification file which declares an interface common to
//                  all supported nonlinear solver algorithms.  The Manager
//                  class uses this interface to call a concrete algorithm.
//
// Special Notes  : This is the "Strategy" class in the Strategy design
//                  pattern.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_ParameterSet_h
#define Xyce_N_NLS_NOX_ParameterSet_h

#include <vector>

#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_NoCase.h>
#include <N_TIA_DataStore.h>

#include <N_IO_InitialConditions.h>

#include "NOX.H"
#include "Teuchos_RCP.hpp"

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS_NonLinearSolver
// Purpose       : Nonlinear Solver Abstract Class
// Creator       : Tammy Kolda, SNL, 8950
// Creation Date : 2/5/02
//-----------------------------------------------------------------------------
class ParameterSet {

public:

  ParameterSet(Xyce::Nonlinear::AnalysisMode mode);
  ~ParameterSet();
  bool setOptions(const Xyce::Util::OptionBlock& OB);
  bool setLocaOptions(const Xyce::Util::OptionBlock& OB, bool saveCopy=true);
  bool applySavedLocaOptions()
  {
    bool ret=true;
    if (savedLocaOptions_)
    {
      ret=setLocaOptions(savedLocaOB_, false);
    }
    return ret;
  }

  bool setOutputOptions(int myPID, int outputProcess);
  bool createStatusTests(
    Parallel::Machine                           comm,
    Xyce::TimeIntg::DataStore*                  data_store,
    Xyce::Loader::NonlinearEquationLoader&      nonlinear_equation_loader,
    Xyce::Linear::Solver &                      linear_solver, 
    Xyce::Linear::Vector *                      maskVectorPtr = 0);             

  Teuchos::RCP<NOX::StatusTest::Generic> getStatusTests();
  bool getVectorParam (const std::string &, int, double &);
  bool getVectorParam (const std::string &, int, std::string &);
  int getVectorParamSize(const std::string& vectorName);
  int getStatusTestReturnCode() const;
  void setStatusTestReturnCodes (const Xyce::Nonlinear::ReturnCodes & retCodesTmp);

  Teuchos::RCP<Teuchos::ParameterList> getAllParams();
  Teuchos::RCP<Teuchos::ParameterList> getNoxParams();
  Teuchos::RCP<Teuchos::ParameterList> getLocaParams();
  Teuchos::RCP<Teuchos::ParameterList> getDebugParams();
  int getNoxSolverType() const;
  void setNoxSolverType(int type);

  bool getContinuationSpecifiedFlag () const
  {
    return continuationSpecified_;
  }

  inline int  getDebugLevel() const
  {
    return debugLevel_;
  };

  inline int getDebugMinTimeStep() const
  {
    return debugMinTimeStep_;
  };

  inline int getDebugMaxTimeStep() const
  {
    return debugMaxTimeStep_;
  };

  inline double getDebugMinTime() const
  {
    return debugMinTime_;
  };

  inline double getDebugMaxTime() const
  {
    return debugMaxTime_;
  };

  inline bool getScreenOutputFlag () const
  {
    return screenOutputFlag_;
  };

  double getMaxNormF() const;
  int getMaxNormFindex () const;

  inline bool isParamsSet() const
  {
    return isParamsSet_;
  };

  inline void set_gstepping_min_value (double val)
  {
    gstepping_min_value_=val;
  }

  inline void set_gstepping_minimum_conductance (double val)
  {
    gstepping_minimum_conductance_ = val;
  }

  Teuchos::RCP<AugmentLinSys>
    createAugmentLinearSystem(Xyce::Linear::System* ls) const;

  // Create augmented linear system, IC version
    Teuchos::RCP<AugmentLinSys>
    createAugmentLinearSystem(Xyce::Linear::System* ls,
                              Xyce::IO::InitialConditionsData::NodeLidValueMap & op,
                              bool gminStepping=false) const;

private:

  void unsupportedOption_(const std::string& tag);
  bool parseOptionBlock_(const Xyce::Util::OptionBlock& OB);

private:

  Teuchos::RCP<Teuchos::ParameterList> allParams_;
  Teuchos::ParameterList& noxParams_;
  Teuchos::ParameterList& locaParams_;
  Teuchos::ParameterList& debugParams_;
  Teuchos::ParameterList statusTestParams_;

  std::map<std::string, std::vector<Xyce::Util::Param> > vectorParams;

  // Combo of all tests
  Teuchos::RCP<NOX::StatusTest::Combo> comboPtr_;

  // Vector containing all the tests we create
  std::vector< Teuchos::RCP<NOX::StatusTest::Generic> > tests_;

  bool isParamsSet_;
  bool isStatusTestsSet_;

  bool continuationSpecified_;

  Xyce::Nonlinear::AnalysisMode mode_;

  int noxSolver;

  // This object acts as a factory for the augment system strategy.
  enum VoltageListType {
    VLT_DOFS,
    VLT_Node,
    VLT_None
  };

  VoltageListType voltageListType_;

  double voltageScaleFactor_;

  // Minimum parameter value used in gstepping
  double gstepping_min_value_;

  // Evil, evil gmin stepping hack, a residual conductance that should almost
  // always be zero.
  double gstepping_minimum_conductance_;

  bool savedLocaOptions_;
  Xyce::Util::OptionBlock savedLocaOB_;

  // debug info:
  int  debugLevel_;
  int debugMinTimeStep_;
  int debugMaxTimeStep_;
  double debugMinTime_;
  double debugMaxTime_;
  bool screenOutputFlag_;
  bool maskingFlag_;
};

}}} // namespace N_NLS_NOX

#endif


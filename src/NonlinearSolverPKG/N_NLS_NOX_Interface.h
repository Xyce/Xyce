//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Special Notes  : 
//
// Creator        : Tammy Kolda
//
// Creation Date  : 01/31/2002
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_Interface_h
#define Xyce_N_NLS_NOX_Interface_h

#include <N_IO_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>

#include <N_IO_InitialConditions.h>
#include <N_NLS_Manager.h>	// defines AnalysisMode
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_NOX_ParameterSet.h>

namespace NOX {
  namespace Parameter {
    class List;
  }
  namespace StatusTest {
    class Generic;
  }
}

//-----------------------------------------------------------------------------
// Class         : N_NLS_NonLinearSolver
// Purpose       : Nonlinear Solver Abstract Class
// Creator       : Tammy Kolda, SNL, 8950
// Creation Date : 2/5/02
//-----------------------------------------------------------------------------

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

  class Interface : public Xyce::Nonlinear::NonLinearSolver
  {

    public:

      Interface(const Xyce::IO::CmdParse & cp);
      ~Interface();

      bool setOptions(const Xyce::Util::OptionBlock& OB);
      bool setTranOptions(const Xyce::Util::OptionBlock& OB);
      bool setHBOptions(const Xyce::Util::OptionBlock& OB);
      bool setNLPOptions(const Xyce::Util::OptionBlock& OB);

      bool setLocaOptions(const Xyce::Util::OptionBlock& OB);
      bool setICOptions(const Xyce::Util::OptionBlock& OB);
      bool setNodeSetOptions(const Xyce::Util::OptionBlock& OB);
      bool initializeAll();

      int solve (Xyce::Nonlinear::NonLinearSolver * nlsTmpPtr = NULL);

      int spiceStrategy ( ParameterSet* paramsPtr );

      int stdNewtonSolve ( ParameterSet* paramsPtr );
      int naturalParameterContinuationSolve ( ParameterSet* paramsPtr );
      int mosfetContinuationSolve ( ParameterSet* paramsPtr );
      int mosfetContinuationSolve2 ( ParameterSet* paramsPtr );
      int mosfetContinuationSolve3 ( ParameterSet* paramsPtr );
      int mosfetContinuationSolve4 ( ParameterSet* paramsPtr );
      int mosfetContinuationSolve5 ( ParameterSet* paramsPtr );
      int mosfetContinuationSolve6 ( ParameterSet* paramsPtr );
      int blockGainscaleMosfetSolve ( ParameterSet* paramsPtr );
      int gminSteppingSolve ( ParameterSet* paramsPtr );
      int pseudoTransientSolve ( ParameterSet* paramsPtr );
      int artificialParameterHomotopy ( ParameterSet* paramsPtr );
      int sourceSteppingSolve ( ParameterSet* paramsPtr );

      int takeFirstSolveStep (Xyce::Nonlinear::NonLinearSolver * nlsTmpPtr = NULL);
      int takeOneSolveStep ();

      Teuchos::RCP<N_NLS_LOCA::Group> getSolutionGroup ();

      int getNumIterations() const;
      double getMaxNormF() const;
      int getMaxNormFindex() const;

      int getDebugLevel() const;
      bool getScreenOutputFlag () const;
      double getDebugMinTime() const;
      double getDebugMaxTime() const;
      int getDebugMinTimeStep() const;
      int getDebugMaxTimeStep() const;
      bool getMMFormat () const;

      // Returns the continuation step number if available.
      int getContinuationStep() const;

      // Returns the parameter number:
      int getParameterNumber() const;

      // Returns true if this is the first continuation param
      bool isFirstContinuationParam() const;

      // Returns true if this is the first solve has been completed
      bool isFirstSolveComplete() const;
      bool getLocaFlag ();
      void setAnalysisMode(Xyce::Nonlinear::AnalysisMode mode);
      void resetAll (Xyce::Nonlinear::AnalysisMode mode);

      bool computeF();

      bool computeJacobian();

      bool applyJacobian(const Xyce::Linear::Vector& input, Xyce::Linear::Vector& result);

      bool computeNewton(Teuchos::ParameterList& p);

      bool computeDfDpMulti (const std::vector< int > & paramIDs, 
                             NOX::Abstract::MultiVector & dfdp, 
                             bool isValidF);

      Xyce::Loader::NonlinearEquationLoader& getLoader() const;

    protected:
      // Resets the stepper by destroying and reallocating.
      void resetStepper(const Teuchos::RCP<LOCA::GlobalData>& gd,
          const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
          const Teuchos::RCP<NOX::StatusTest::Generic>& test,
          const Teuchos::RCP<Teuchos::ParameterList>& p);

      bool icCont (ParameterSet* paramsPtr);
      bool icCont3 (ParameterSet* paramsPtr);

      bool nodesetCont0 (ParameterSet* paramsPtr);
      bool nodesetCont1 (ParameterSet* paramsPtr);
    private:

      // Parameters for DC_OP
      ParameterSet dcParams_;

      bool ICspecified_;
      bool NODESETspecified_;

      // Parameters for Transient
      ParameterSet transientParams_;

      // Parameters for HB: 
      ParameterSet hbParams_;

      // Parameters for NLPoisson: 
      ParameterSet nlpParams_;

      // Shared system
      SharedSystem* sharedSystemPtr_;

      // Global data for loca groups
      Teuchos::RCP<LOCA::GlobalData> globalDataPtr_;

      // LOCA Wrapper Status Tests
      Teuchos::RCP<LOCA::StatusTest::Wrapper> locaTransientStatusTestPtr_;
      Teuchos::RCP<LOCA::StatusTest::Wrapper> locaDCOpStatusTestPtr_;
      Teuchos::RCP<LOCA::StatusTest::Wrapper> locaStatusTestPtr_;

      Teuchos::RCP<LOCA::StatusTest::Wrapper> locaHBStatusTestPtr_;
      Teuchos::RCP<LOCA::StatusTest::Wrapper> locaDC_NLPStatusTestPtr_;

      // Nox group
      Teuchos::RCP<N_NLS_LOCA::Group> groupPtr_;

      // NOX Solver
      Teuchos::RCP<NOX::Solver::Generic> solverPtr_;

      // LOCA Stepper
      Teuchos::RCP<LOCA::Stepper> stepperPtr_;

      // Current analysis mode
      Xyce::Nonlinear::AnalysisMode mode_;

      // Whether or not we should use the current analysis mode
      bool usemode_;

      // save the parameters mode.
      Xyce::Nonlinear::AnalysisMode lastParametersMode_;
      Xyce::Nonlinear::AnalysisMode parametersMode_;

      bool copiedGroupFlag_;

      // Keep track of whether to set the linear solver tolerance
      // (Only do this if adaptive forcing is on)
      bool setAZ_Tol_DC;
      bool setAZ_Tol_Transient;

      //are we on the first LOCA continuation parameter?
      bool isFirstContinuationParam_;

      //is first solve completed? 
      bool firstSolveComplete_;

      //parameter index
      int iParam_;
  };

}}} // namespace N_NLS_NOX

#endif


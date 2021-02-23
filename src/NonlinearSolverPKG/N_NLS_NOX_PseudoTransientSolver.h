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

//-------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
//
//
//
//-------------------------------------------------------------------------


#ifndef NOX_SOLVER_PSEUDOTRANSIENTSOLVER_H
#define NOX_SOLVER_PSEUDOTRANSIENTSOLVER_H

#include <N_NLS_fwd.h>

#include "NOX_Solver_Generic.H"	         // base class

#include "NOX_LineSearch_Generic.H"      // class data element
#include "NOX_Direction_Generic.H"       // class data element

#include "NOX_Utils.H"		         // class data element
#include "NOX_StatusTest_FiniteValue.H"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"     

#ifdef Xyce_NOX_SOLVERSTATS
#include "NOX_SolverStats.hpp"
namespace NOX {
  class Observer;
}
#else
#include "NOX_Solver_PrePostOperator.H"  // class data element
#endif

namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

  class PseudoTransientBased : public NOX::Solver::Generic {

public:

  //! Constructor
  PseudoTransientBased(const Teuchos::RCP<AugmentLinSys>& als,
		       const Teuchos::RCP<NOX::Abstract::Group>& grp, 
		       const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
		       const Teuchos::RCP<Teuchos::ParameterList>& params,
		       double initialStepSize,
		       double minStepSize,
		       double maxStepSize);

  //! Destructor
  virtual ~PseudoTransientBased() {}

  virtual void reset(const NOX::Abstract::Vector& initial_guess);
  virtual void reset(const NOX::Abstract::Vector& initial_guess,
                     const Teuchos::RCP<NOX::StatusTest::Generic>& test);
#ifdef Xyce_NOX_SOLVERSTATS
  virtual void reset();
  virtual NOX::StatusTest::StatusType getStatus() const;
#else
  virtual NOX::StatusTest::StatusType getStatus();
#endif
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;

  //! Return the line search step size from the current iteration
  virtual double getStepSize() const;

  //! Return the pseudo transient step size.
  virtual double getPseudoTransientStepSize() const;

  //! Return a RCP to the solution group
  virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const;
  
  //! Return a RCP to the previous solution group
  virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const;
  
  //! Return a RCP to the solver parameters.
  virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const;

#ifdef Xyce_NOX_SOLVERSTATS
  //! Return a RCP to the solver statistics.
  virtual Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const { return Teuchos::null; }
#endif

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Global Data.
  Teuchos::RCP<NOX::GlobalData> globalData;

  //! RCP to the strategy for augmenting the linear system.
  Teuchos::RCP<AugmentLinSys> augmentLSStrategy;

  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnPtr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Teuchos::RCP<NOX::Abstract::Group> oldSolnPtr;	
  //! Previous solution reference.
  NOX::Abstract::Group& oldSoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Teuchos::RCP<NOX::Abstract::Vector> dirPtr;
  //! Current search direction.reference.
  NOX::Abstract::Vector& dir;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testPtr;		

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;	

  //! Utils
  NOX::Utils& utils;

  //! Linesearch. 
  Teuchos::RCP<NOX::LineSearch::Generic> lineSearch; 

  //! %Search %Direction. 
  Teuchos::RCP<NOX::Direction::Generic> direction; 

  //! Current step.
  double step_;			

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

#ifdef Xyce_NOX_SOLVERSTATS
  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> prePostOperator;
#else
  //! Pointer to a user defined NOX::Solver::PrePostOperator object.
  Teuchos::RCP<NOX::Solver::PrePostOperator> prePostOperator;
#endif

  double initialStepSize_;
  double minStepSize_;
  double maxStepSize_;
  double stepSize_;
  double previousStepSize_;
  double scaleFactor_;

  N_NLS_LOCA::Group* group_;
  N_NLS_LOCA::Group* previousGroup_;
  
  NOX::StatusTest::FiniteValue fvTest_;

  //! Type of check to use for status tests.
  NOX::StatusTest::CheckType checkType_;

};
}}} // namespace NOX

#endif


//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Particular Status Test Based on the Weighed Norm of the
//                  Update
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
//
//
//
//-------------------------------------------------------------------------
#ifndef Xyce_N_NLS_NOX_XyceTests_h
#define Xyce_N_NLS_NOX_XyceTests_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_NLS_ReturnCodes.h>
#include <N_UTL_FeatureTest.h>
#include <N_TIA_DataStore.h>

// ----------   NOX Includes   ----------

#include "NOX_StatusTest_Generic.H"	// base class
#include "NOX_StatusTest_FiniteValue.H"

// ---------- Namespace Declarations ----------

// N_NLS namespace is for the Xyce Nonlinear Solver Package
namespace Xyce {
namespace Nonlinear {
namespace N_NLS_NOX {

//-----------------------------------------------------------------------------
// Class         : N_NLS_NOX::XyceTests
//
// Purpose       :
//
//      NOX Vector Interface for Xyce vectors.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class XyceTests : public NOX::StatusTest::Generic {

public:

  //---------------------------------------------------------------------------
  // Function      : XyceTests (constructor)
  //
  // Purpose       : Constructs a NOX-compatiable status test based
  //                 the weighted norm of the update.
  //
  //---------------------------------------------------------------------------
  XyceTests(
    Parallel::Machine   comm,
    bool isTransient, 
    double normF,
    double machPrec,
    Xyce::TimeIntg::DataStore* data_store,
    double epsilon_a, 
    double epsilon_r, 
    double tol,
    int maxIters,
    double convRate,
    double relConvRate,
    double maxConvRate,
    double stagnationTol,
    int maxBadSteps,
    int checkDeviceConvergence,
    double smallUpdateTol,
    Xyce::Loader::NonlinearEquationLoader* loader, 
    Xyce::Linear::Solver *                 lsolver,
    bool maskingFlag,
    Xyce::Linear::Vector * maskVectorPtr);
  
  //---------------------------------------------------------------------------
  // Function      : Destructor
  //---------------------------------------------------------------------------
  ~XyceTests();

  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType checkStatus(
      const NOX::Solver::Generic& problem, 
      NOX::StatusTest::CheckType checkType);
  
  //---------------------------------------------------------------------------
  // Purpose       : Test stopping criterion given the current
  //                 nonlinear problem
  //---------------------------------------------------------------------------
  NOX::StatusTest::StatusType getStatus() const { return status_; };
  
  //---------------------------------------------------------------------------
  // Purpose       : Get the return code to send to the time stepper.
  //           
  //---------------------------------------------------------------------------
  int getXyceReturnCode() const;
  
  //---------------------------------------------------------------------------
  // Purpose       : Output formatted description of stopping test to
  //                 output stream.
  //---------------------------------------------------------------------------
  std::ostream& print(std::ostream& stream, int indent = 0) const;
  
  //---------------------------------------------------------------------------
  // Purpose       : Set a specific set of return codes to be used.
  //---------------------------------------------------------------------------
  void setReturnCodes (const Xyce::Nonlinear::ReturnCodes & retCodesTmp);

  double getMaxNormF() const;

  int getMaxNormFindex () const;
  
protected:
  Parallel::Machine             comm_;
  NOX::StatusTest::StatusType   status_;
  int returnTest_;
  bool isTransient_;
  int niters_;
  
  Xyce::Nonlinear::ReturnCodes retCodes_;
  
  //Test #0 
  //****************************************
  NOX::StatusTest::FiniteValue finiteTest_;
    
  //Test #1 
  //****************************************
  int maxNormFindex_;
  double maxNormF_;
  double normF_curr_;
  double requestedMaxNormF_;
  double requestedMachPrecTol_;
  Xyce::Linear::Vector* pWeightsVectorPtr_;

  // Test #2
  //****************************************
  Xyce::TimeIntg::DataStore* dsPtr_;
  Xyce::Linear::Vector* weightsVectorPtr_;
  Xyce::Linear::Vector* updateVectorPtr_;
  const double epsilon_a_;
  const double epsilon_r_;
  const double tol_;
  double weightedUpdate_;

  // Test #3
  //****************************************
  // Maximum number of nonlinear iterations allowed
  int maxIters_;

  // ||F(x_current)|| / ||F(x_previous)||
  const double requestedConvRate_;
  double currentConvRate_;

  //  ||F(x)|| / ||F(x_init)||
  const double requestedRelativeConvRate_;
  double currentRelativeConvRate_;

  // Initial norm of the RHS used to calculate ratio
  double normResidualInit_;
  
  // test #4
  double smallUpdateTol_;

  // Test #6
  const double maxConvRate_;
    
  // Test #7
  int lastIteration_;
  int badStepCount_;
  const int maxBadSteps_;
  double minConvRate_;
  const double stagnationTol_;
    
  // Xyce return Code
  int xyceReturnCode_;

  // For Device Specific Convergence cirteria
  int checkDeviceConvergence_;
  Xyce::Loader::NonlinearEquationLoader* loaderPtr_;
  Xyce::Linear::Solver * lasSolverPtr_; 

  bool maskingFlag_;
  Xyce::Linear::Vector* weightMaskVectorPtr_;

  bool allDevicesConverged_;
  bool innerDevicesConverged_;

}; // class XyceTests

inline void XyceTests::setReturnCodes 
  (const Xyce::Nonlinear::ReturnCodes & retCodesTmp)
{
  retCodes_ = retCodesTmp;
}


}}} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_XyceTests_h


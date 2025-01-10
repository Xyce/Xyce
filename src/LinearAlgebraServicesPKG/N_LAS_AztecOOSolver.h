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
//
// Purpose        : AztecOO Interface
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/18/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_AztecOOSolver_h
#define Xyce_N_LAS_AztecOOSolver_h

#include <Teuchos_RCP.hpp>
#include <N_LAS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_LAS_TransformTool.h>
#include <N_LAS_Solver.h>

class AztecOO;

class Epetra_LinearProblem;

class Ifpack_CrsRiluk;
class Ifpack_IlukGraph;

#include <N_UTL_Param.h>

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : AztecOOSolver
// Purpose       : AztecOO Interface
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/04
//-----------------------------------------------------------------------------
class AztecOOSolver : public Solver
{

public:

  //Constructors
  AztecOOSolver(
    Problem &     problem,
    Util::OptionBlock & options);

  // Destructor
  ~AztecOOSolver();

  //Control of solver options
  bool setOptions( const Util::OptionBlock & OB );
  bool setDefaultOptions();
  bool setParam( const Util::Param & param );
  bool setPreconditioner( const Teuchos::RCP<Preconditioner>& precond )
  { precond_ = precond; isPrecSet_ = false; return true; }

  //Retrieve of solver information
  bool getInfo( Util::Param & info );

  //Actual Solve Call
  int doSolve( bool reuse_factors, bool transpose = false );

  //Residual Value
  double residual() { return linearResidual_; }

  // These are all accessor functions to get at AztecOO options.  Please see
  // the Aztec User's Guide for complete definitions and values.  Each option
  // can has accessors as well as a reset to default.

  // Preconditioner Flag
  const int & getPreCond() const
  { return preCond_; }
  void setPreCond(const int & value)
  { preCond_ = value; setAztecOption_( "AZ_precond", preCond_ ); }
  void resetPreCond()
  { setPreCond( preCondDefault_ ); }

  // Solver Flag
  const int & getSolver() const
  { return solverType_; }
  void setSolver(const int & value)
  { solverType_ = value; setAztecOption_( "AZ_solver", solverType_ ); }
  void resetSolver()
  { setSolver( solverDefault_ ); }

  // Scaling Flag
  const int & getScaling() const
  { return scaling_; }
  void setScaling(const int & value)
  { scaling_ = value; setAztecOption_( "AZ_scaling", scaling_); }
  void resetScaling()
  { setScaling( scalingDefault_ ); }

  // Subdomain Solve Flag
  const int & getSubdomainSolve() const
  { return subdomainSolve_; }
  void setSubdomainSolve(const int & value)
  { subdomainSolve_ = value; setAztecOption_( "AZ_subdomain_solve", subdomainSolve_ ); }
  void resetSubdomainSolve()
  { setSubdomainSolve( subdomainSolveDefault_ ); }

  // Convergence Test Flag
  const int & getConvergence() const
  { return convergence_; }
  void setConvergence(const int & value)
  { convergence_ = value; setAztecOption_( "AZ_conv", convergence_ ); }
  void resetConvergence()
  { setConvergence( convergenceDefault_ ); }

  // Output Flag
  const int & getOutput() const
  { return output_; }
  void setOutput(const int & value)
  { output_ = value; setAztecOption_( "AZ_output", output_ ); }
  void resetOutput()
  { setOutput( outputDefault_ ); }

  // Diagnostic Flag
  const int & getDiagnostics() const
  { return diagnostics_; }
  void setDiagnostics(const int & value)
  { diagnostics_ = value; setAztecOption_( "AZ_diagnostics", diagnostics_ ); }
  void resetDiagnostics()
  { setDiagnostics( diagnosticsDefault_ ); }

  // Factorization Reuse
  const int & getPrecalc() const
  { return precalc_; }
  void setPrecalc(const int & value)
  { precalc_ = value; setAztecOption_( "AZ_pre_calc", precalc_ ); }
  void resetPrecalc()
  { setPrecalc( precalcDefault_ ); }

  // Maximum Iterations
  const int & getMaxIter() const
  { return maxIter_; }
  void setMaxIter(const int & value)
  { maxIter_ = value; setAztecOption_( "AZ_max_iter", maxIter_ ); }
  void resetMaxIter()
  { setMaxIter( maxIterDefault_ ); }

  // Domain Decomposition Overlap
  const int & getOverlap() const
  { return overlap_; }
  void setOverlap(const int & value)
  { overlap_ = value; setAztecOption_( "AZ_overlap", overlap_ ); }
  void resetOverlap()
  { setOverlap( overlapDefault_ ); }

  // Krylov-subspace size for GMRES
  const int & getKSpace() const
  { return KSpace_; }
  void setKSpace(const int & value)
  { KSpace_ = value; setAztecOption_( "AZ_kspace", KSpace_ ); }
  void resetKSpace()
  { setKSpace( KSpaceDefault_ ); }

  // RCM reorder
  const int & getReorder() const
  { return reorder_; }
  void setReorder(const int & value)
  { reorder_ = value; setAztecOption_( "AZ_reorder", reorder_ ); }
  void resetReorder()
  { setReorder( reorderDefault_ ); }

  // factorization information-retention flag
  const int & getKeepInfo() const
  { return keepInfo_; }
  void setKeepInfo(const int & value)
  { keepInfo_ = value; setAztecOption_( "AZ_keep_info", keepInfo_ ); }
  void resetKeepInfo()
  { setKeepInfo( keepInfoDefault_ ); }

  // GMRES orthogonalization scheme
  const int & getOrthog() const
  { return orthog_; }
  void setOrthog(const int & value)
  { orthog_ = value; setAztecOption_( "AZ_orthog", orthog_ ); }
  void resetOrthog()
  { setOrthog( orthogDefault_ ); }

  // Linear Convergence Tolerance
  const double & getTolerance() const
  { return tolerance_; }
  void setTolerance(const double & value)
  { tolerance_ = value; setAztecParam_( "AZ_tol", tolerance_ ); }
  void resetTolerance()
  { setTolerance( toleranceDefault_ ); }

  // ILU/ILUT drop tolerance
  const double & getDrop() const
  { return drop_; }
  void setDrop(const double & value)
  { drop_ = value; setAztecParam_( "AZ_drop", drop_ ); }
  void resetDrop()
  { setDrop( dropDefault_ ); }

  // ILUT fill-factor
  const double & getILUTFill() const
  { return ilutFill_; }
  void setILUTFill(const double & value)
  { ilutFill_ = value; setAztecParam_( "AZ_ilut_fill", ilutFill_ ); }
  void resetILUTFill()
  { setILUTFill( ilutFillDefault_ ); }

  // Diagonal-Shifting Relative Threshold
  const double & getRThresh() const
  { return rThresh_; }
  void setRThresh(const double & value)
  { rThresh_ = value; setAztecParam_( "AZ_rthresh", rThresh_ ); }
  void resetRThresh()
  { setRThresh( rThreshDefault_ ); }

  // Diagonal-Shifting Absolute Threshold
  const double & getAThresh() const
  { return aThresh_; }
  void setAThresh(const double & value)
  { aThresh_ = value; setAztecParam_( "AZ_athresh", aThresh_ ); }
  void resetAThresh()
  { setAThresh( aThreshDefault_ ); }

  // Achieved Linear Residual
  const double & getLinearResidual() const
  { return linearResidual_; }
  // Number of Linear Iterations
  const unsigned int & getNumLinearIters() const
  { return numLinearIters_; }

private:

  // For large problems the number of Krylov vectors may need to be resized for AztecOO.
  bool reduceKSpace_;
  int maxKSpace_;

  // Set Aztec options
  bool setAztecOption_(const char * paramName, const int val);
  // Set Aztec params
  bool setAztecParam_(const char * paramName, const double val);

  // Trilinos problem difficulty value - this is a temporary solution.
  int probDiff_;

  // Iterative solver parameters - these are largely taken from Aztec
  // options/parameters.
  // Preconditioner
  int preCond_;
  // Solution algorithm
  int solverType_;
  // Scaling algorithm
  int scaling_;
  // Solver for subdomain solves for domain decomposition preconditioners are
  // used
  int subdomainSolve_;
  // Residual expression used for convergence test
  int convergence_;
  // Output information flag
  int output_;
  // Diagnostic information flag
  int diagnostics_;
  // Use previous factorization info?
  int precalc_;
  // Max iterations
  int maxIter_;
  // Submatrices factored with domain decomposition
  int overlap_;
  // Krylov subspace size for restarted GMRES
  int KSpace_;
  // Reordering flag
  int reorder_;
  // Keep matrix factorization information flag
  int keepInfo_;
  // GMRES orthogonalization type flag
  int orthog_;

  // Tolerance for convergence tests
  double tolerance_;
  // Drop tolerance for LU or ILUT preconditioners
  double drop_;
  // Fill-factor for ILUT preconditioner
  double ilutFill_;
  // Diagonal shifting relative threshold
  double rThresh_;
  // Diagonal shifting absolute threshold
  double aThresh_;

  // Defaults for parameters above - see comment above for descriptions.
  // !!!NOTE: these values are based on the Aztec include files on or about
  // June, 2001.  We will eventually need a more consistent way of setting
  // these - SAH

  const int preCondDefault_;
  const int solverDefault_;
  const int scalingDefault_;
  const int subdomainSolveDefault_;
  const int convergenceDefault_;
  const int outputDefault_;
  const int diagnosticsDefault_;
  const int precalcDefault_;
  const int maxIterDefault_;
  const int KSpaceDefault_;
  const int reorderDefault_;
  const int keepInfoDefault_;
  const int orthogDefault_;
  const int overlapDefault_;

  const double toleranceDefault_;
  const double dropDefault_;
  const double ilutFillDefault_;
  const double rThreshDefault_;
  const double aThreshDefault_;

  // Linear residual returned from the linear solver
  double linearResidual_;
  // Number of linear iterations returned from the linear solver
  unsigned int numLinearIters_;

  //Output linear system every outputLS_ calls
  int outputLS_;
  int outputBaseLS_;

  // Linear Problem
  Epetra_LinearProblem * problem_;

  // Option block pointer
  Util::OptionBlock * options_;

  // Flag for ifpack use
  bool useAztecPrecond_;
  bool isPrecSet_;

  // Flag for change of params
  bool updatedParams_;

  // Pointer to the Aztec linear solver
  AztecOO * solver_;

  // Transform support
  Teuchos::RCP<Transform> transform_;
  Epetra_LinearProblem * tProblem_;

  Teuchos::RCP<Preconditioner> precond_;

  // Timing tool
  Util::Timer * timer_;

  bool setAztecCntl_( const Util::Param & param );

  static int num_AZ_options_;
  static const char * AZ_options_[];

  static int num_AZ_params_;
  static const char * AZ_params_[];

};

} // namespace Linear
} // namespace Xyce

#endif // Xyce_N_LAS_AztecOOSolver_h

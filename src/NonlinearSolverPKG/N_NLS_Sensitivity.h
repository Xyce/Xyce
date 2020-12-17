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
// Purpose        : This file contains the sensitivity class.   It mostly
//                  manages the calculations of direct (and possibly later,
//                  adjoint) sensitivities.
//
// Special Notes  : The main reason that this class is derived from
//                  N_NLS_NonLinearSolver is that this class needs to
//                  do a series of linear solves, using the jacobian
//                  matrix.  This seemed similar enough to the requirements
//                  of a nonlinear solver to have one derived off the other.
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/30/02
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_Sensitivity_h
#define Xyce_N_NLS_Sensitivity_h

#include<vector>

#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_MachDepParams.h>
#include <N_NLS_NonLinearSolver.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : objective function data
// Purpose       :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/17/2015
//-----------------------------------------------------------------------------
class objectiveFunctionData
{
public:  
  objectiveFunctionData() :
    numExpVars(0),
    expVal(0.0),
    objFuncString(""),
    expPtr(0),
    numDdt(0),
    objFuncEval(0.0),
    dOdp(0.0)
  {};

public:
  int            numExpVars;  // this is just the size of the numVarDerivs array
  std::vector<std::string> expVarNames; // need this to get the GIDs
  std::vector<int>    expVarGIDs;  // need this even with new expression lib, to populate the dOdX vector

  std::vector<double> expVarDerivs; // this is returned by expPtr->evaluate
  double         expVal; // this is returned by expPtr->evaluate

  std::vector<int> globalParamVariableStencil;

  std::string objFuncString;

  Util::Expression * expPtr;

  int numDdt;

  double objFuncEval;// value of the evaluated objective function.
  double dOdp;

  Linear::Vector* dOdXVectorPtr; // size of solution vector.
};

bool evaluateObjFuncs ( 
    std::vector<objectiveFunctionData*> & objVec, 
    Parallel::Communicator & comm,
    Loader::NonlinearEquationLoader & nlEquLoader_,
    TimeIntg::DataStore & dataStore,
    TimeIntg::StepErrorControl & sec,
    std::string & netlistFilename);

void setupObjectiveFunctions (
    Teuchos::RCP<Xyce::Util::baseExpressionGroup> & exprGroup,
    std::vector<objectiveFunctionData*> & objVec,
    IO::OutputMgr & output_manager, Linear::System & lasSys,
    const IO::CmdParse &cp,
    bool checkTimeDeriv=true);

void setupObjectiveFuncGIDs (std::vector<objectiveFunctionData*> & objVec, Parallel::Communicator& comm, 
    Topo::Topology & top, IO::OutputMgr & output_manager);

void applyHocevarDelayTerms(
    std::vector<objectiveFunctionData*> & objVec,
    std::vector<objectiveFunctionData*> & objTimeDerivVec,
    TimeIntg::DataStore & dataStore
    );

//-----------------------------------------------------------------------------
// Class         : Sensitivity
// Purpose       : This class handles many aspects of sensitivity analysis.
//
// Special Notes : It is a little strange to derive this class from the 
//                 base nonlinear solver.  However, doing so makes it easy for 
//                 sensitivity analysis to use exactly the same Jacobian, 
//                 linear solver, LU factors, etc. that were used by the Newton 
//                 solver.  So, for convenience it was derived that way.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
class Sensitivity : public NonLinearSolver
{
public:
  static void populateMetadata(IO::PkgOptionsMgr &options_manager);

  Sensitivity (
    NonLinearSolver *           nls, 
    Topo::Topology &            topTmp,
    const IO::CmdParse &        cp,
    TimeIntg::StepErrorControl & secTmp);

  ~Sensitivity ();

   bool icSensitivity (
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec);

   int solve (NonLinearSolver * nlsTmpPtr=NULL) {return -1;};
   int solve (
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec);

   int solveDirect  ();
   int solveAdjoint ();

   int solveTransientAdjoint (bool timePoint,
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec);

   std::ostream& stdOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities,
       std::ostream& os
       );

   void fileOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities
       );

   bool calcObjFuncDerivs ();
   bool calcObjFuncTimeDerivs ();

   bool setOptions(const Util::OptionBlock& OB);
   bool setSensitivityOptions(const Util::OptionBlock& OB);
   bool doAllocations();

   bool setTranOptions(const Util::OptionBlock& OB) {return true;}
   bool setHBOptions(const Util::OptionBlock& OB)   {return true;}
   bool setNLPOptions(const Util::OptionBlock& OB)  {return true;}

   void resetNLS(NonLinearSolver * nls);

   // Note, many of the following are here b/c they are purely
   // virtual functions of the nonlinear solver class.
   // They don't have much meaning here.  It may turn out that
   // having this class derive off of the NonLinearSolver
   // class doesn't make much sense.  If so, I'll change it later. ERK

   int getNumIterations() const;

   int getDebugLevel() const;
   bool getScreenOutputFlag() const;
   double getDebugMinTime() const;
   double getDebugMaxTime() const;
   int getDebugMinTimeStep() const;
   int getDebugMaxTimeStep() const;
   bool getMMFormat () const;

   double getMaxNormF() const { return nls_->getMaxNormF(); }
   int getMaxNormFindex() const { return nls_->getMaxNormFindex (); }

   int getContinuationStep() const;
   int getParameterNumber() const;
   bool isFirstContinuationParam() const;
   bool isFirstSolveComplete() const;
   void setAnalysisMode(AnalysisMode mode);

   int getNumSensParams () { return numSensParams_;}

private:
  int debugLevel_;
  int solutionSize_;
  bool solveDirectFlag_;
  bool solveAdjointFlag_;
  bool outputScaledFlag_; // include scaled sensitivities in IO 
  bool outputUnscaledFlag_; // include unscaled sensitivities in IO
  int maxParamStringSize_;

  bool stdOutputFlag_;
  bool fileOutputFlag_;
  int numSolves_;

  bool objFuncGiven_;
  bool objFuncGIDsetup_;
  std::vector<objectiveFunctionData*> objFuncDataVec_;

  bool objFuncTimeDerivGIDsetup_;
  std::vector<objectiveFunctionData*> objFuncTimeDerivDataVec_;

  // finite difference variables
  int difference_;
  double sqrtEta_;
  bool sqrtEtaGiven_;
  bool forceFD_;
  bool forceDeviceFD_;
  bool forceAnalytic_;
  bool newLowMem_;
  bool sparseAdjointStorage_;
  bool computeDelays_;
  bool timeDerivsSetup_;

  bool reuseFactors_;

  Linear::Vector * lambdaVectorPtr_;
  Linear::Vector * savedRHSVectorPtr_;
  Linear::Vector * savedNewtonVectorPtr_;

  NonLinearSolver * nls_;

  Topo::Topology & top_;

  TimeIntg::StepErrorControl & sec;

  int numSensParams_;
  int numObjectives_;
  std::vector<std::string> paramNameVec_;

  // Current analysis mode
  Xyce::Nonlinear::AnalysisMode mode_;
};

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getNumIterations
// Purpose       : doesn't do anything, is just a placeholder.
// Special Notes : This one may be needed later, I'm not sure.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
inline int Sensitivity::getNumIterations() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline int Sensitivity::getContinuationStep() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline int Sensitivity::getParameterNumber() const
{
  return 0;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline bool Sensitivity::isFirstContinuationParam() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Dummy function since homotopy doesn't work with old solver
//-----------------------------------------------------------------------------
inline bool Sensitivity::isFirstSolveComplete() const
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setAnalysisMode
// Purpose       : doesn't do anything, is just a placeholder.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
inline void Sensitivity::setAnalysisMode(AnalysisMode mode)
{
  mode_ = mode;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getDebugLevel
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
inline int Sensitivity::getDebugLevel() const
{
  return -100;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::getScreenOutputFlag
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
inline bool Sensitivity::getScreenOutputFlag () const
{
  return false;
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMinTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double Sensitivity::getDebugMinTime() const
{
  return 0.0;
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMaxTime
//
// Return Type   : double
//---------------------------------------------------------------------------
inline double Sensitivity::getDebugMaxTime() const
{
  return Util::MachineDependentParams::DoubleMax();
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMinTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int Sensitivity::getDebugMinTimeStep() const
{
  return 0;
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getDebugMaxTimeStep
//
// Return Type   : int
//---------------------------------------------------------------------------
inline int Sensitivity::getDebugMaxTimeStep() const
{
  return Util::MachineDependentParams::IntMax();
}

//---------------------------------------------------------------------------
// Function      : Sensitivity::getMMFormat
//
// Return Type   : bool
//---------------------------------------------------------------------------
inline bool Sensitivity::getMMFormat () const
{
  return false;
}

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_Sensitivity_h


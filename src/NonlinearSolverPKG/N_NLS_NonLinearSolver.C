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
// Purpose        : Body file which declares an interface common to all
//                  supported nonlinear solver algorithms.  The Manager class
//                  uses this interface to call a concrete algorithm.
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

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>

#include <N_ANP_AnalysisManager.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_PrecondFactory.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Solver.h>
#include <N_LAS_TranSolverFactory.h>
#include <N_LAS_System.h>
#include <N_LAS_SystemHelpers.h>
#include <N_LAS_Vector.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_PDS_Manager.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

#include <expressionGroup.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::NonLinearSolver
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
NonLinearSolver::NonLinearSolver(const IO::CmdParse &cp)
  : commandLine_(cp),
    netlistFilename_(""),
    rhsVectorPtr_(0),
    jacobianMatrixPtr_(0),
    gradVectorPtr_(0),
    NewtonVectorPtr_(0),
    solWtVectorPtr_(0),
    maskPNormWeights_(0),
    lasSysPtr_(0),
    lasSolverFactoryPtr_(0),
    lasPrecFactoryPtr_(0),
    linsolOptionBlockPtr_(0),
    nonlinearEquationLoader_(0),
    tlnPtr_(0),
    nonlinearParameterManager_(0),
    outMgrPtr_(0),
    initialConditionsManager_(0),
    pdsMgrPtr_(0),
    dsPtr_(0),
    numJacobianLoads_(0),
    numJacobianFactorizations_(0),
    numLinearSolves_(0),
    numFailedLinearSolves_(0),
    numResidualLoads_(0),
    totalNumLinearIters_(0),
    totalLinearSolveTime_(0.0),
    totalResidualLoadTime_(0.0),
    totalJacobianLoadTime_(0.0),
    debugTimeFlag_(true),
    contStep_(0),
    analysisManager_(0),
    outputStepNumber_(0),
    matrixFreeFlag_(false)
{
  NonLinearSolver::resetCountersAndTimers_();

  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netlistFilename_ = commandLine_.getArgumentValue("netlist");
  }
}


//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::~NonLinearSolver
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
NonLinearSolver::~NonLinearSolver()
{
  if (gradVectorPtr_)
  {
    delete gradVectorPtr_;
    gradVectorPtr_ = 0;
  }

  if (solWtVectorPtr_)
  {
    delete solWtVectorPtr_;
    solWtVectorPtr_ = 0;
  }

  if (maskPNormWeights_)
  {
    delete maskPNormWeights_;
    maskPNormWeights_ = 0;
  }

  if (linsolOptionBlockPtr_)
  {
    delete linsolOptionBlockPtr_;
    linsolOptionBlockPtr_ = 0;
  }

  if (matrixFreeFlag_)
  {
    delete lasProblemPtr_;
    lasProblemPtr_ = 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setLinsolOptions
//
// Purpose       : Passes option block corresponding to "LINSOL" onto
//                 nonlinear solver.
//
// Special Notes:  These options are saved to be passed onto the object from
//                 registerLinearSolver() in the initializeAll() function.
//                 This is only called if "NONLIN" is present in the
//                 circuit file.
//
// Return Type   : boolean
//
// - Input Arguments -
//
//    OB         : Option block containing options corresponding to
//                 "LINSOL" in the netlist.
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/9/00
//-----------------------------------------------------------------------------
bool NonLinearSolver::setLinsolOptions(const Util::OptionBlock & OB)
{
  linsolOptionBlockPtr_ = new Util::OptionBlock(OB);
  return (linsolOptionBlockPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setICOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/17/2007
//-----------------------------------------------------------------------------
bool NonLinearSolver::setICOptions(const Util::OptionBlock& OB)
{
  Report::DevelFatal0() <<  ".IC options not supported for this nonlinear solver.  Use nox instead. ";
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setNodeSetOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 09/25/2007
//-----------------------------------------------------------------------------
bool NonLinearSolver::setNodeSetOptions(const Util::OptionBlock& OB)
{
  Report::DevelFatal0() <<  ".NODESET options not supported for this nonlinear solver.  Use nox instead. ";
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool NonLinearSolver::setLocaOptions (const Util::OptionBlock& OB)
{
  Report::DevelFatal0() <<  "NonLinearSolver::setLocaOptions - not implemented for this solver. Use nox instead. ";
  return true;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::setTwoLevelLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool NonLinearSolver::setTwoLevelLocaOptions (const Util::OptionBlock& OB)
{
  Report::DevelFatal0() <<  "NonLinearSolver::setTwoLevelLocaOptions - not implemented for this solver.  Use nox instead.";
  return true;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::setTwoLevelOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//---------------------------------------------------------------------------
bool NonLinearSolver::setTwoLevelOptions (const Util::OptionBlock& OB)
{
  return true;
}

//---------------------------------------------------------------------------
// Function      : NonLinearSolver::setTwoLevelTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//---------------------------------------------------------------------------
bool NonLinearSolver::setTwoLevelTranOptions (const Util::OptionBlock& OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerNonlinearEquationLoader(Loader::NonlinearEquationLoader* nonlinear_equation_loader)
{
  nonlinearEquationLoader_ = nonlinear_equation_loader;
  return (nonlinearEquationLoader_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerLinearSystem
// Purpose       :
// Special Notes :
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerLinearSystem(Linear::System* tmp_LasSysPtr)
{
  lasSysPtr_ = tmp_LasSysPtr;
  return (lasSysPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerSolverFactory
// Purpose       :
// Special Notes :
// Return Type   : boolean
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerSolverFactory(const Linear::SolverFactory *tmp_LasSolverPtr)
{
  lasSolverFactoryPtr_ = tmp_LasSolverPtr;
  return lasSolverFactoryPtr_;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerPrecondFactory
// Purpose       :
// Special Notes :
// Return Type   : boolean
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerPrecondFactory(const Linear::PrecondFactory *tmp_LasPrecPtr)
{
  lasPrecFactoryPtr_ = tmp_LasPrecPtr;
  return lasPrecFactoryPtr_;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerParallelMgr
// Purpose       :
// Special Notes :
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/8/2013
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerParallelMgr(Parallel::Manager * pdsMgrPtr)
{
  pdsMgrPtr_ = pdsMgrPtr;
  return (pdsMgrPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerAnalysisManager
// Purpose       :
// Special Notes :
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerAnalysisManager(Analysis::AnalysisManager* tmp_anaIntPtr)
{
  analysisManager_ = tmp_anaIntPtr;
  return (analysisManager_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/23/03
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerOutputMgr (IO::OutputMgr * outPtr)
{
  outMgrPtr_ = outPtr;
  return (outMgrPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerInitialConditionsManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/23/03
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerInitialConditionsManager(IO::InitialConditionsManager * ptr)
{
  initialConditionsManager_ = ptr;
  return (initialConditionsManager_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerExpressionGroup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/2/2020
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerExpressionGroup(Teuchos::RCP<Xyce::Util::baseExpressionGroup> & group)
{
  expressionGroup_ = group;
  return ( !(Teuchos::is_null(expressionGroup_)) );
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerTwoLevelSolver
// Purpose       : This function is called in the event that the two-level
//                 Newton method has been invoked.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/24/02
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerTwoLevelSolver
  (TwoLevelNewton * tmp_tlnPtr)
{
  tlnPtr_ = tmp_tlnPtr;
  return (tlnPtr_ != 0);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::registerTIADataStore
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool NonLinearSolver::registerTIADataStore(TimeIntg::DataStore * ds_tmp)
{
  dsPtr_ = ds_tmp;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::initializeAll
//
// Purpose       : Called after all register and set functions.
//                 Once the various registrations have taken place,
//                 this function sets the remaining pointers.
//
// Special Notes:  This function obtains the solution, temporary solution and
//                 rhs vectors from the LAS system class.
//
// Return Type   : boolean
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool NonLinearSolver::initializeAll()
{
  bool bsuccess = true;

  // Check the registerLinearSystem has been successfully called.
  if (lasSysPtr_ == 0)
    return false;

  // get rhs vector:
  rhsVectorPtr_ = lasSysPtr_->getRHSVector();
  bsuccess = bsuccess && (rhsVectorPtr_ != 0);

  // get newton vector:
  NewtonVectorPtr_ = lasSysPtr_->getNewtonVector();
  bsuccess = bsuccess && (NewtonVectorPtr_ != 0);

  // get jacobian:
  if (!matrixFreeFlag_)
  {
    jacobianMatrixPtr_ = lasSysPtr_->getJacobianMatrix();
    bsuccess = bsuccess && (jacobianMatrixPtr_ != 0);
  }

  Linear::Builder & bld_ = lasSysPtr_->builder();

  // create gradient vector:
  gradVectorPtr_ = bld_.createVector();
  bsuccess = bsuccess && (gradVectorPtr_ != 0);

  // create solution weighting vector:
  solWtVectorPtr_ = bld_.createVector();
  bsuccess = bsuccess && (solWtVectorPtr_ != 0);

  if( !linsolOptionBlockPtr_ )
  {
    linsolOptionBlockPtr_ = new Util::OptionBlock();
    linsolOptionBlockPtr_->addParam(Util::Param("TYPE", "DEFAULT"));
  }

  if (!matrixFreeFlag_)
  {
    // Normal full matrix linear solver options
    lasProblemPtr_ = lasSysPtr_->getLinearProblem();
  }
  else
  {
    // Matrix free harmonic balance linear problem 
    lasProblemPtr_ = Xyce::Linear::createProblem( this,
                                                  NewtonVectorPtr_,
                                                  rhsVectorPtr_ );
  }

  if ( Teuchos::is_null(lasSolverRCPtr_) )
  {
    if (lasSolverFactoryPtr_)
    {
      lasSolverRCPtr_ = Teuchos::rcp( lasSolverFactoryPtr_->create( *linsolOptionBlockPtr_, *lasProblemPtr_ , commandLine_) );
    }
    else
    {
      Linear::TranSolverFactory lasFactory;
      lasSolverRCPtr_ = Teuchos::rcp( lasFactory.create( *linsolOptionBlockPtr_, *lasProblemPtr_ , commandLine_) );
    }
  }

  // If a preconditioner factory has been provided by the analysis package,
  // use it to generate a preconditioner for the linear solver.
  if (lasPrecFactoryPtr_) 
  {
    Teuchos::RCP<Linear::Preconditioner> precond = lasPrecFactoryPtr_->create( Teuchos::rcp( lasSysPtr_, false ) );
    lasSolverRCPtr_->setPreconditioner( precond );
  }

  if (DEBUG_NONLINEAR)
    Xyce::dout() << "size of solution vector: " << NewtonVectorPtr_->globalLength() << std::endl
                 << "End of NonLinearSolver::initializeAll\n";

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::debugOutput1_
// Purpose       : Write Jacobian to the file matrix.(n).txt, and the residual
//                 to the file rhs.(n).txt
//
// Special Notes : These are objects that will be availabled prior to the
//                 linear solve performed at each Newton step.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/00
//-----------------------------------------------------------------------------
void NonLinearSolver::debugOutput1(
  Linear::Matrix &        jacobian,
  Linear::Vector &        rhs)
{
  setNonlinearDumpDebugLevel(getDebugLevel());
  int newtStep = getNumIterations();
  int screenOutput = getScreenOutputFlag ();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  if (!debugTimeFlag_ || !isActive(Diag::NONLINEAR_DUMP_MASK)) return;

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;

  if (isActive(Diag::NONLINEAR_DUMP_PARAM_NUMBER))
  {
    sprintf(filename1, "matrix_%03d_%03d_%03d_%03d.txt", outputStepNumber_, paramNumber, contStep, newtStep);
  }
  if (isActive(Diag::NONLINEAR_DUMP_STEP))
  {
    sprintf(filename1, "matrix_%03d_%03d.txt", outputStepNumber_, newtStep);
  }
  if (isActive(Diag::NONLINEAR_DUMP))
  {
    sprintf(filename1, "matrix_%03d.txt", newtStep);
  }

  jacobian.writeToFile(filename1, false, getMMFormat () );

  if (screenOutput == 1)
  {
    Xyce::dout() << "\n\t***** Jacobian matrix:" << std::endl;
    jacobian.print(Xyce::dout());
  }

  if (isActive(Diag::NONLINEAR_DUMP_PARAM_NUMBER))
  {
    sprintf(filename2, "rhs_%03d_%03d_%03d_%03d.txt", outputStepNumber_, paramNumber, contStep, newtStep);
  }
  if (isActive(Diag::NONLINEAR_DUMP_STEP))
  {
    sprintf(filename2, "rhs_%03d_%03d.txt", outputStepNumber_, newtStep);
  }
  else
  {
    sprintf(filename2, "rhs_%03d.txt", newtStep);
  }

  if (screenOutput == 1)
  {
    Xyce::dout() << "\n\t***** RHS vector:" << std::endl;

    rhs.print(Xyce::dout());
  }

  rhs.writeToFile(filename2);

  debugOutputDAE ();

}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::debugOutputDAE
// Purpose       : Write DAE vectors and matrices to output files.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/27/03
//-----------------------------------------------------------------------------
void NonLinearSolver::debugOutputDAE()
{
  setNonlinearDumpDebugLevel(getDebugLevel());

  int newtStep = getNumIterations();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;

  char filename4[256]; for (int ich = 0; ich < 256; ++ich) filename4[ich] = 0;
  char filename6[256]; for (int ich = 0; ich < 256; ++ich) filename6[ich] = 0;
  char filename6b[256];for (int ich = 0; ich < 256; ++ich) filename6b[ich] = 0;
  char filename7[256]; for (int ich = 0; ich < 256; ++ich) filename7[ich] = 0;
  char filename8[256]; for (int ich = 0; ich < 256; ++ich) filename8[ich] = 0;
  char filename9[256]; for (int ich = 0; ich < 256; ++ich) filename9[ich] = 0;

  Linear::Matrix *dQdx    = dsPtr_->dQdxMatrixPtr;
  Linear::Matrix *dFdx    = dsPtr_->dFdxMatrixPtr;

  Linear::Vector *daeQ    = dsPtr_->daeQVectorPtr;
  Linear::Vector *daeF    = dsPtr_->daeFVectorPtr;
  Linear::Vector *daeB    = dsPtr_->daeBVectorPtr;

  Linear::Vector *daeFlim = lasSysPtr_->getdFdxdVpVector ();
  Linear::Vector *daeQlim = lasSysPtr_->getdQdxdVpVector ();

  if (isActive(Diag::NONLINEAR_DUMP_PARAM_NUMBER))
  {
    sprintf(filename1, "dQdx_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);
    sprintf(filename2, "dFdx_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);

    sprintf(filename4, "daeQ_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);
    sprintf(filename6, "daeF_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);
    sprintf(filename6b,"daeB_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);

    sprintf(filename8, "daeQlim_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);
    sprintf(filename9, "daeFlim_%03d_%03d_%03d_%03d.txt"    , outputStepNumber_, paramNumber, contStep, newtStep);
  }
  else if (isActive(Diag::NONLINEAR_DUMP_STEP))
  {
    sprintf(filename1, "dQdx_%03d_%03d.txt"    , outputStepNumber_, newtStep);
    sprintf(filename2, "dFdx_%03d_%03d.txt"    , outputStepNumber_, newtStep);

    sprintf(filename4, "daeQ_%03d_%03d.txt"    , outputStepNumber_, newtStep);
    sprintf(filename6, "daeF_%03d_%03d.txt"    , outputStepNumber_, newtStep);
    sprintf(filename6b,"daeB_%03d_%03d.txt"    , outputStepNumber_, newtStep);

    sprintf(filename8, "daeQlim_%03d_%03d.txt"    , outputStepNumber_, newtStep);
    sprintf(filename9, "daeFlim_%03d_%03d.txt"    , outputStepNumber_, newtStep);
  }
  else
  {
    sprintf(filename1, "dQdx_%03d.txt"    , newtStep);
    sprintf(filename2, "dFdx_%03d.txt"    , newtStep);

    sprintf(filename4, "daeQ_%03d.txt"    , newtStep);
    sprintf(filename6, "daeF_%03d.txt"    , newtStep);
    sprintf(filename6b,"daeB_%03d.txt"    , newtStep);

    sprintf(filename8, "daeQlim_%03d.txt"    , newtStep);
    sprintf(filename9, "daeFlim_%03d.txt"    , newtStep);
  }

  // write the matrices:
  dQdx->writeToFile (filename1, false, getMMFormat () );
  dFdx->writeToFile (filename2, false, getMMFormat () );

  // write the vectors:
  daeQ->writeToFile(filename4);
  daeF->writeToFile(filename6);
  daeB->writeToFile(filename6b);
  daeQlim->writeToFile(filename8);
  daeFlim->writeToFile(filename9);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::outputDAEvectors
// Purpose       : Write DAE vectors to output files.
//
// Special Notes : Similar to debugOutputDAE but doesn't check the debuglevel, 
//                 and only outputs the F, Q and B vectors.
//
//                 Also is only called after time step is 100% converged.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/2022
//-----------------------------------------------------------------------------
void NonLinearSolver::outputDAEvectors()
{
  int outputStepNum = 0;
  if (analysisManager_) outputStepNum = analysisManager_->getStepNumber() + 1;

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename4[256]; for (int ich = 0; ich < 256; ++ich) filename4[ich] = 0;
  char filename6[256]; for (int ich = 0; ich < 256; ++ich) filename6[ich] = 0;
  char filename6b[256];for (int ich = 0; ich < 256; ++ich) filename6[ich] = 0;

  Linear::Vector *daeQ    = dsPtr_->daeQVectorPtr;
  Linear::Vector *daeF    = dsPtr_->daeFVectorPtr;
  Linear::Vector *daeB    = dsPtr_->daeBVectorPtr;
  Linear::Vector *xVec    = dsPtr_->nextSolutionPtr;

  Linear::Vector *daeFlim = lasSysPtr_->getdFdxdVpVector ();
  Linear::Vector *daeQlim = lasSysPtr_->getdQdxdVpVector ();

  sprintf(filename1, "%s_soln_%03d.txt", netlistFilename_.c_str(), outputStepNum);
  sprintf(filename4, "%s_daeQ_%03d.txt", netlistFilename_.c_str(), outputStepNum);
  sprintf(filename6, "%s_daeF_%03d.txt", netlistFilename_.c_str(), outputStepNum);
  sprintf(filename6b,"%s_daeB_%03d.txt", netlistFilename_.c_str(), outputStepNum);

  // write the vectors:
  xVec->writeToFile(filename1);
  daeQ->writeToFile(filename4);
  daeF->writeToFile(filename6);
  daeB->writeToFile(filename6b);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::outputDAEmatrices
// Purpose       : Write DAE matrices to output files.
//
// Special Notes : Similar to debugOutputDAE but doesn't check the debuglevel, 
//                 and only outputs the F and Q matrices.
//
//                 Only called after time step is 100% converged.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/2022
//-----------------------------------------------------------------------------
void NonLinearSolver::outputDAEmatrices()
{
  int outputStepNum = 0;
  if (analysisManager_) outputStepNum = analysisManager_->getStepNumber() + 1;

  char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
  char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;

  Linear::Matrix *dQdx    = dsPtr_->dQdxMatrixPtr;
  Linear::Matrix *dFdx    = dsPtr_->dFdxMatrixPtr;

  sprintf(filename1, "%s_dQdx_%03d.txt", netlistFilename_.c_str(), outputStepNum);
  sprintf(filename2, "%s_dFdx_%03d.txt", netlistFilename_.c_str(), outputStepNum);

  // write the matrices:
  dQdx->writeToFile (filename1, false, getMMFormat () );
  dFdx->writeToFile (filename2, false, getMMFormat () );
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::debugOutput3_
// Purpose       : Write out the update vector and the new solution.
//
// Special Notes : These are objects that will be availabled *after* the
//                 linear solve performed at each Newton step.  That
//                 differentiates this function from debugOutput1_.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/00
//-----------------------------------------------------------------------------
void NonLinearSolver::debugOutput3(
  Linear::Vector &        dxVector,
  Linear::Vector &        xVector)
{
  setNonlinearDumpDebugLevel(getDebugLevel());

  int nlStep = getNumIterations();
  int contStep = getContinuationStep();
  int paramNumber = getParameterNumber ();

  if (!debugTimeFlag_ || !isActive(Diag::NONLINEAR_DUMP_MASK)) return;

  char filename[256];  for (int ich = 0; ich < 256; ++ich) filename[ich] = 0;

  if (isActive(Diag::NONLINEAR_DUMP_PARAM_NUMBER))
  {
    sprintf(filename, "update_%03d_%03d_%03d_%03d.txt", outputStepNumber_, paramNumber, contStep, nlStep);
  }
  else if (isActive(Diag::NONLINEAR_DUMP_STEP))
  {
    sprintf(filename, "update_%03d_%03d.txt", outputStepNumber_, nlStep);
  }
  else
  {
    sprintf(filename, "update_%03d.txt", nlStep);
  }
  dxVector.writeToFile(filename);

  if (isActive(Diag::NONLINEAR_DUMP_PARAM_NUMBER))
  {
    sprintf(filename, "solution_%03d_%03d_%03d_%03d.txt", outputStepNumber_, paramNumber, contStep, nlStep);
  }
  if (isActive(Diag::NONLINEAR_DUMP_STEP))
  {
    sprintf(filename, "solution_%03d_%03d.txt", outputStepNumber_, nlStep);
  }
  if (isActive(Diag::NONLINEAR_DUMP))
  {
    sprintf(filename, "solution_%03d.txt", nlStep);
  }
  xVector.writeToFile(filename);
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::resetCountersAndTimers_
//
// Purpose       : Resets all the counters and timers in this object.
//
// Scope         : protected
// Creator       : Tamara G. Kolda, SNL, CSMR (8950)
//                 Eric Keiter, SNL, Parallel Computational Sciences. (9233)
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
void NonLinearSolver::resetCountersAndTimers_()
{
  numJacobianLoads_ = 0;
  numJacobianFactorizations_ = 0;
  numLinearSolves_ = 0;
  numFailedLinearSolves_ = 0;
  numResidualLoads_ = 0;
  totalNumLinearIters_ = 0;
  totalLinearSolveTime_ = 0.0;
  totalResidualLoadTime_ = 0.0;
  totalJacobianLoadTime_ = 0.0;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setX0_()
//
// Purpose       : This should be called at the beginning of each nonlinear
//                 iteration. Copies information from nextSolVector (and
//                 related vectors that are important but hidden from the
//                 nonlinear solver) into tmpSolVector.
//
// Return Type   : boolean
// Scope         : protected
// Creator       : Tamara G. Kolda, SNL, CSMR (8950)
//                 Eric Keiter, SNL, Parallel Computational Sciences. (9233)
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool NonLinearSolver::setX0_()
{

  dsPtr_->equateTmpVectors();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::rhs_()
//
// Purpose       : Calculates the RHS corresponding to the current solution
//                 vector. More specifically, it fills in the content of
//                 RHSVectorPtr_ based on the contents of nextSolVectorPtr_.
//
// Special Notes : The rhsVectorPtr_ is really the NEGATIVE of F(x).
//
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------

bool NonLinearSolver::rhs_()
{
  Stats::StatTop _residualStat("Residual");
  Stats::TimeBlock _residualTimer(_residualStat);

  nonlinearEquationLoader_->loadRHS();
  ++numResidualLoads_;
  totalResidualLoadTime_ += nonlinearEquationLoader_->getResidualTime();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::jacobian_()
//
// Purpose       : Calculates the Jacobian corresponding to the current
//                 solution vector. More specifically, it fills in the
//                 content of jacobianMatrixPtr_ based on the contents of
//                 nextSolVectorPtr_.
//
// Return Type   : boolean
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool NonLinearSolver::jacobian_()
{
  Stats::StatTop _jacobianStat("Jacobian");
  Stats::TimeBlock _jacobianTimer(_jacobianStat);

  nonlinearEquationLoader_->loadJacobian();
  ++numJacobianLoads_;
  totalJacobianLoadTime_ += nonlinearEquationLoader_->getJacobianTime();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::applyJacobian()
//
// Purpose       : Applies the Jacobian corresponding to the current
//                 solution vector.
//
// Return Type   : boolean
// Scope         : private
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool NonLinearSolver::applyJacobian(const Linear::Vector& input, Linear::Vector& result) 
{
  Stats::StatTop _jacobianStat("Apply Jacobian");
  Stats::TimeBlock _jacobianTime(_jacobianStat);

  nonlinearEquationLoader_->applyJacobian(input, result);
  ++numJacobianLoads_;
  totalJacobianLoadTime_ += nonlinearEquationLoader_->getJacobianTime();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::newton_
//
// Purpose       : Calculates the Newton direction corresponding to the
//                 current RHS and Jacobian matrix.
//
// Return Type   : boolean
//
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool NonLinearSolver::newton_()
{
  int solutionStatus = lasSolverRCPtr_->solve( false );

  totalLinearSolveTime_ += lasSolverRCPtr_->solutionTime();
  ++numLinearSolves_;

  if( lasSolverRCPtr_->isIterative() )
  {
    Util::Param param( "Iterations", 0 );
    lasSolverRCPtr_->getInfo( param );
    totalNumLinearIters_ += param.getImmutableValue<int>();

    if( solutionStatus ) ++numFailedLinearSolves_;
  }
  else
  {
    Util::Param param( "Refactored", 0 );
    lasSolverRCPtr_->getInfo( param );
    if( param.getImmutableValue<int>() ) ++numJacobianFactorizations_;
    if( solutionStatus ) ++numFailedLinearSolves_;
  }

  if( solutionStatus ) return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::gradient_()
// Purpose       : Calculates the Gradient direction corresponding to the
//                 current RHS and Jacobian matrix.
//                 Computes gradient using jacobianMatrixPtr_ and
//                 the rhsVectorPtr_. On output, gradVectorPtr_ contains
//                 the gradient of 0.5 * ||F(x)||^2.
//
// Special Notes : The rhsVectorPtr_ is really the NEGATIVE of F(x).
//
// Return Type   : boolean
// Scope         : private
// Creator       : Tamara G. Kolda, SNL, Compuational Sciences and
//                 Mathematics Research Department
// Creation Date : 06/19/01
//-----------------------------------------------------------------------------
bool NonLinearSolver::gradient_()
{
  // Compute gradient = Jacobian' * RHS
  bool transpose = true;
  jacobianMatrixPtr_->matvec(transpose, *rhsVectorPtr_, *gradVectorPtr_);

  // We have to scale by -1 because the rhsVectorPtr_ is really the
  // NEGATIVE of F(X). This gives us gradient = Jacobian' * F(X).
  gradVectorPtr_->scale(-1.0);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::getCouplingMode ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Compuational Sciences and
// Creation Date : 12/05/02
//-----------------------------------------------------------------------------
TwoLevelNewtonMode NonLinearSolver::getCouplingMode ()
{
  return FULL_PROBLEM;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::setDebugFlags
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Compuational Sciences
// Creation Date : 09/23/01
//-----------------------------------------------------------------------------
void
NonLinearSolver::setDebugFlags(
  int           output_step_number,
  double        time)
{
  outputStepNumber_ = output_step_number;

  debugTimeFlag_ =
    (time >= getDebugMinTime() && time <= getDebugMaxTime())
    && (output_step_number >= getDebugMinTimeStep() && output_step_number <= getDebugMaxTimeStep());

  if (tlnPtr_ != 0)
    contStep_ = tlnPtr_->getContStepNumber();
  else
    contStep_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::apply
// Purpose       : Apply matrix free operator with Linear::MultiVectors
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 9/4/08
//-----------------------------------------------------------------------------
int NonLinearSolver::apply(
  const Linear::MultiVector& X, 
  Linear::MultiVector& Y,
  Teuchos::ETransp mode ) const
{
  if (mode == Teuchos::TRANS || mode == Teuchos::CONJ_TRANS)
    return -1;

  // Have to const cast here, fix const correctness in NL/Loader layer later.
  NonLinearSolver* nc_this = const_cast<NonLinearSolver* >( this );

  bool status = true;
  for (int i=0 ; i<X.numVectors() ; ++i)
  {
    Teuchos::RCP<const Linear::Vector> X_i = Teuchos::rcp( X.getVectorView(i) );
    Teuchos::RCP<Linear::Vector> Y_i = Teuchos::rcp( Y.getNonConstVectorView(i) );
    bool localStatus = nc_this->applyJacobian(*X_i, *Y_i);
    status = status && localStatus;
  }

  if (status)
    return 0;
  else
    return -1;
}

//-----------------------------------------------------------------------------
// Function      : NonLinearSolver::getPNormWeights
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi K. Thornquist, SNL
// Creation Date : 10/27/21
//-----------------------------------------------------------------------------
Linear::Vector* NonLinearSolver::getPNormWeights()
{
  // Create the mask if it hasn't been made before.
  if (!maskPNormWeights_)
  {
    int length = dsPtr_->nextSolutionPtr->localLength();
    Linear::Vector& mask = *(lasSysPtr_->getDeviceMaskVector());

    maskPNormWeights_ = lasSysPtr_->builder().createVector();
    maskPNormWeights_->putScalar( 1.0 );

    for (int i = 0; i < length; ++i)
    {
      if (mask[i] == 0.0)
        (*(maskPNormWeights_))[i] = Util::MachineDependentParams::MachineBig();
    }
  }

  return maskPNormWeights_;
}


} // namespace Nonlinear
} // namespace Xyce

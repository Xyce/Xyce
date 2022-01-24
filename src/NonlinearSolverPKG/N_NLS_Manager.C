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

//-------------------------------------------------------------------------
//
// Purpose        : Body for the nonlinear solver manager class implementation
//                  which will allow for the creation and selection of NLS
//                  algorithms.
//
// Special Notes  : GOF Strategy Pattern
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_NLS_ConductanceExtractor.h>
#include <N_NLS_DampedNewton.h>
#include <N_NLS_Manager.h>
#include <N_NLS_NOX_Interface.h>
#include <N_NLS_NonLinInfo.h>
#include <N_NLS_Sensitivity.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_TIA_DataStore.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>

// ----------  Trilinos Includes   ----------

namespace Xyce {
namespace Nonlinear {

enum {
  OPTION_BLOCK_DCOP,
  OPTION_BLOCK_TRANSIENT,
  OPTION_BLOCK_NLP,
  OPTION_BLOCK_HB,
  OPTION_BLOCK_LINSOL,
  OPTION_BLOCK_LOCA,
  OPTION_BLOCK_TWO_LEVEL_LOCA,
  OPTION_BLOCK_TWO_LEVEL,
  OPTION_BLOCK_TWO_LEVEL_TRAN,
  OPTION_BLOCK_IC,
  OPTION_BLOCK_NODESET,
  OPTION_BLOCK_SENS,
  OPTION_BLOCK_SENSITIVITY
};


//-----------------------------------------------------------------------------
// Function      : Manager::Manager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
Manager::Manager(
  const IO::CmdParse & command_line)
  : commandLine_(command_line),
    nonlinearSolver_(0),
    conductanceExtractorPtr_(0),
    nlsSensitivityPtr_(0),
    lasSolverPtr_(0),
    lasPrecPtr_(0),
    matrixFreeFlag_(false),
    noxFlag_(false),
    noxFlagInner_(false),
    noxFlagTransient_(false),
    optionBlockMap_(),
    initializeAllFlag_(false),
    retCodes_()
    //,
    //exprPtr_(0) // ERK. why?
{}


//-----------------------------------------------------------------------------
// Function      : Manager::~Manager
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

Manager::~Manager()
{
  delete nonlinearSolver_;
  delete nlsSensitivityPtr_;
  delete conductanceExtractorPtr_;
  //delete exprPtr_;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/29/00
//-----------------------------------------------------------------------------
bool Manager::setOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_DCOP] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setTranOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/05/01
//-----------------------------------------------------------------------------
bool Manager::setTranOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_TRANSIENT] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setNLPOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/07/2015
//-----------------------------------------------------------------------------
bool Manager::setNLPOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_NLP] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setHBOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/03/2009
//-----------------------------------------------------------------------------
bool Manager::setHBOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_HB] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::getHBOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/03/2009
//-----------------------------------------------------------------------------
Util::OptionBlock &
Manager::getHBOptions()
{
  return optionBlockMap_[OPTION_BLOCK_HB];
}

//-----------------------------------------------------------------------------
// Function      : Manager::setLinSolOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/9/00
//-----------------------------------------------------------------------------
bool Manager::setLinSolOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_LINSOL] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool Manager::setLocaOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_LOCA] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setTwoLevelLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool Manager::setTwoLevelLocaOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_TWO_LEVEL_LOCA] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setTwoLevelOptions
// Purpose       : This option setter will create a new class structure.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool Manager::setTwoLevelOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_TWO_LEVEL] = option_block;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setTwoLevelTranOptions
// Purpose       : This option setter will create a new class structure.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool Manager::setTwoLevelTranOptions(const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_TWO_LEVEL_TRAN] = option_block;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : setReturnCodeOption
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/18/02
//-----------------------------------------------------------------------------
bool
Manager::setReturnCodeOption(
  const Util::Param &           param)
{
  if (param.uTag() == "NLNEARCONV" )
  {
    // nlNearConvFlag = static_cast<bool> (param.getImmutableValue<int>());
    Nonlinear::ReturnCodes codes = getReturnCodes();
    codes.nearConvergence = param.getImmutableValue<int>() ? 3 : -3;
    setReturnCodes(codes);
  }
  else if (param.uTag() == "NLSMALLUPDATE" )
  {
    // nlSmallUpdateFlag = static_cast<bool> (param.getImmutableValue<int>());
    Nonlinear::ReturnCodes codes = getReturnCodes();
    codes.smallUpdate = param.getImmutableValue<int>() ? 4 : -4;
    setReturnCodes(codes);
  }
  else
  {
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::using Nox_
// Purpose       : This function determines if we are using NOX or not.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
void Manager::usingNox()
{
  noxFlag_ = true;
  noxFlagInner_ = true;
  noxFlagTransient_ = false;

  const Util::OptionBlock &option_blockdcop = optionBlockMap_[OPTION_BLOCK_DCOP];
  const Util::OptionBlock &option_blocktran = optionBlockMap_[OPTION_BLOCK_TRANSIENT];

  // scan for changes to the dcop values of noxFlag_
  for (Util::ParamList::const_iterator it = option_blockdcop.begin(); it != option_blockdcop.end(); ++ it)
  {
    if ((*it).uTag() == "NOX")
    {
      noxFlag_ = (*it).getImmutableValue<int>();
      noxFlagInner_ = noxFlag_;
    }
  }

  // now check for a nox flag on the .options twolevel line, if it exists.
  OptionBlockMap::const_iterator obmIter = optionBlockMap_.find(OPTION_BLOCK_TWO_LEVEL) ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_blocktwoLevel = obmIter->second;
    for (Util::ParamList::const_iterator it = option_blocktwoLevel.begin(); it != option_blocktwoLevel.end(); ++ it)
    {
      if ((*it).uTag() == "NOX")
      {
        noxFlagInner_ = (*it).getImmutableValue<int>();
      }
    }
  }


  // scan for changes to the transient value for noxFlagTransient_
  for (Util::ParamList::const_iterator it = option_blocktran.begin(); it != option_blocktran.end(); ++ it)
  {
    if ((*it).uTag() == "NOX")
    {
      noxFlagTransient_ = (*it).getImmutableValue<int>();
    }
  }

  obmIter = optionBlockMap_.find(OPTION_BLOCK_HB) ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_blockhb = obmIter->second;
    for (Util::ParamList::const_iterator it = option_blockhb.begin(); it != option_blockhb.end(); ++ it)
    {
      if ((*it).uTag() == "NOX")
      {
        noxFlag_ = (*it).getImmutableValue<int>();
        noxFlagInner_ = noxFlag_;
      }
    }
  }

  // now check if the command line has specified nox.  The command line
  // overrides the netlist.
  if( commandLine_.getArgumentValue( "-nox" ) == "off" )
  {
    noxFlag_ = true;  // Need to use NOX for DC solve
    noxFlagInner_ = false;
    noxFlagTransient_ = false;
  }
  else if( commandLine_.getArgumentValue( "-nox" ) == "on" )
  {
    noxFlag_ = true;
    noxFlagInner_ = true;
    noxFlagTransient_ = true;
  }

  if (DEBUG_NONLINEAR)
    Xyce::dout() << "noxFlag is: " << (noxFlag_ ? "true" : "false") << std::endl
                 << "noxFlagTransient is: " <<  (noxFlagTransient_ ? "true" : "false") << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : Manager::allocateSolver
// Purpose       : This function determines which solver to allocate, and
//                 allocates it.
//
//                 Right now the possibilities are:
//
//                    DampedNewton
//                    NOXInterface
//                    TwoLevelNewton
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/28/02
//-----------------------------------------------------------------------------
bool Manager::allocateSolver(
  Analysis::AnalysisManager &           analysis_manager,
  Loader::NonlinearEquationLoader &     nonlinear_equation_loader,
  Linear::System &                      linear_system,
  TimeIntg::DataStore &                 data_store,
  Parallel::Manager &                   parallel_manager,
  IO::InitialConditionsManager &        initial_conditions_manager,
  IO::OutputMgr &                       output_manager)
{
  bool bsuccess = true;
  bool bs1 = true;

  // determine if we are using NOX or not.
  usingNox();

  // if ".options nonlin two-level" appeared in the netlist, then
  // allocate the two-level solver.  Otherwise allocate one of the single
  // level solvers.

  delete nonlinearSolver_;
  if (optionBlockMap_.find(OPTION_BLOCK_TWO_LEVEL) != optionBlockMap_.end()
      || optionBlockMap_.find(OPTION_BLOCK_TWO_LEVEL_TRAN) != optionBlockMap_.end())
  {
    twoLevelNewtonFlag_ = true;
    nonlinearSolver_ = new TwoLevelNewton(noxFlag_, noxFlagInner_, commandLine_);
  }
  else
  {
    twoLevelNewtonFlag_ = false;
    if (noxFlag_)
    {
      nonlinearSolver_ = new N_NLS_NOX::Interface(commandLine_);
    }
    else
    {
      nonlinearSolver_ = new DampedNewton(commandLine_);
    }
  }

  // now register everything, now that the solver class is set up.
  bs1 = nonlinearSolver_->registerLinearSystem(&linear_system);    bsuccess = bsuccess && bs1;
  bs1 = nonlinearSolver_->registerAnalysisManager(&analysis_manager); bsuccess = bsuccess && bs1;
  bs1 = nonlinearSolver_->registerNonlinearEquationLoader(&nonlinear_equation_loader); bsuccess = bsuccess && bs1;
  bs1 = nonlinearSolver_->registerTIADataStore(&data_store);        bsuccess = bsuccess && bs1;

  if (lasSolverPtr_)
    bs1 = nonlinearSolver_->registerSolverFactory(lasSolverPtr_); bsuccess = bsuccess && bs1;

  if (lasPrecPtr_)
    bs1 = nonlinearSolver_->registerPrecondFactory(lasPrecPtr_); bsuccess = bsuccess && bs1;
  
  OptionBlockMap::const_iterator it = optionBlockMap_.find(OPTION_BLOCK_DCOP) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_TRANSIENT) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setTranOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_NLP) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setNLPOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_HB) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setHBOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_LOCA) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setLocaOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_TWO_LEVEL) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setTwoLevelOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_TWO_LEVEL_TRAN) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setTwoLevelTranOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_TWO_LEVEL_LOCA) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setTwoLevelLocaOptions (option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_LINSOL) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setLinsolOptions(option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_IC) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setICOptions(option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_NODESET) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = (*it).second;
    bs1 = nonlinearSolver_->setNodeSetOptions(option_block);
    bsuccess = bsuccess && bs1;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Manager::initializeAll
// Purpose       : This can only be called after the linear system
//                 (Linear::System) has been registered, and after all the
//                 options have been set.
//
// Special Notes : This function obtains the solution, temporary solution and
//                 rhs vectors from the LAS system class.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool Manager::initializeAll(
  Analysis::AnalysisManager &           analysis_manager,
  Loader::NonlinearEquationLoader &     nonlinear_equation_loader,
  Linear::System &                      linear_system,
  TimeIntg::DataStore &                 data_store,
  Parallel::Manager &                   parallel_manager,
  IO::InitialConditionsManager &        initial_conditions_manager,
  IO::OutputMgr &                       output_manager,
  Topo::Topology &                      topology)
{
  bool bsuccess = true;
  bool bs1 = true;
  bool bs2 = true;

  bs1 = allocateSolver(analysis_manager, nonlinear_equation_loader, linear_system, data_store, parallel_manager, initial_conditions_manager, output_manager);
  bsuccess = bsuccess && bs1;

  nonlinearSolver_->setMatrixFreeFlag(matrixFreeFlag_);
  nonlinearSolver_->registerParallelMgr(&parallel_manager);
  nonlinearSolver_->registerInitialConditionsManager(&initial_conditions_manager);
  nonlinearSolver_->registerOutputMgr(&output_manager);
  bs2 = nonlinearSolver_->initializeAll();  bsuccess = bsuccess && bs2;

  nonlinearSolver_->setReturnCodes(retCodes_);

  initializeAllFlag_ = true;

  if (!conductanceExtractorPtr_)
  {
    conductanceExtractorPtr_ = new ConductanceExtractor(*nonlinearSolver_, topology);
  }

  // the sensitivity pointer isn't necessarily allocated yet, even if this is a .SENS
  // run.  However, if running .SENS with .STEP, this step is necessary to avoid
  // seg faults.
  if (nlsSensitivityPtr_)
  {
    nlsSensitivityPtr_->resetNLS( nonlinearSolver_ );
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Manager::allocateTranSolver
// Purpose       : Allocate a different solver for transient problems if
//                 the user has requested a different one.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
void Manager::allocateTranSolver(
  Analysis::AnalysisManager &           analysis_manager,
  Loader::NonlinearEquationLoader &     nonlinear_equation_loader,
  Linear::System &                      linear_system,
  TimeIntg::DataStore &                 data_store,
  Parallel::Manager &                   parallel_manager,
  IO::OutputMgr &                       output_manager,
  Topo::Topology &                      topology)
{
  bool bsuccess = true;
  bool bs1 = true;

  // Always perform the reallocation of the nonlinear solver.
  // if ((noxFlag_ && !noxFlagTransient_) || (!noxFlag_&& noxFlagTransient_))
  {
    // Save the linear solver from the DCOP solve and pass it to the transient NLS
    Teuchos::RCP<Linear::Solver> lasSolverRCPtr = nonlinearSolver_->getLinearSolver();
 
    delete nonlinearSolver_;
    delete conductanceExtractorPtr_;

    if (noxFlagTransient_)
    {
      nonlinearSolver_ = new N_NLS_NOX::Interface(commandLine_);
    }
    else
    {
      nonlinearSolver_ = new DampedNewton(commandLine_);
    }

    // Pass the old linear solver to the new nonlinear solver.
    nonlinearSolver_->setLinearSolver( lasSolverRCPtr );

    OptionBlockMap::const_iterator it = optionBlockMap_.find(OPTION_BLOCK_LINSOL) ;
    if ( it != optionBlockMap_.end() )
    {
      const Util::OptionBlock &option_block = (*it).second;
      bs1 = nonlinearSolver_->setLinsolOptions(option_block);
      bsuccess = bsuccess && bs1;
    }

    it = optionBlockMap_.find(OPTION_BLOCK_TRANSIENT) ;
    if (it != optionBlockMap_.end())
    {
      const Util::OptionBlock &option_block = (*it).second;
      bs1 = nonlinearSolver_->setTranOptions (option_block);
      bsuccess = bsuccess && bs1;
    }

    // now register everything, now that the solver class is set up.
    nonlinearSolver_->registerLinearSystem(&linear_system);
    nonlinearSolver_->registerAnalysisManager(&analysis_manager);
    nonlinearSolver_->registerNonlinearEquationLoader(&nonlinear_equation_loader);
    nonlinearSolver_->registerTIADataStore(&data_store);
    nonlinearSolver_->registerOutputMgr(&output_manager);
    nonlinearSolver_->registerParallelMgr(&parallel_manager);

    nonlinearSolver_->setMatrixFreeFlag(matrixFreeFlag_);
    nonlinearSolver_->initializeAll();
    nonlinearSolver_->setReturnCodes(retCodes_);

    initializeAllFlag_ = true;

    if (nlsSensitivityPtr_)
    {
      nlsSensitivityPtr_->resetNLS( nonlinearSolver_ );
    }
    conductanceExtractorPtr_ = new ConductanceExtractor (*nonlinearSolver_, topology);
  }
}

//-----------------------------------------------------------------------------
// Function      : Manager::getNonLinInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Parallel Computational Sciences
// Creation Date : 12/05/02
//-----------------------------------------------------------------------------
NonLinInfo Manager::getNonLinInfo() const
{
  NonLinInfo nlInfo;

  nlInfo.newtonIter    = nonlinearSolver_->getNumIterations();
  nlInfo.twoLevelNewtonCouplingMode  = nonlinearSolver_->getCouplingMode ();
  nlInfo.locaFlag      = nonlinearSolver_->getLocaFlag ();

  if (nlInfo.locaFlag)
  {
    nlInfo.continuationStep       = nonlinearSolver_->getContinuationStep ();
    nlInfo.firstContinuationParam = nonlinearSolver_->isFirstContinuationParam ();
    nlInfo.firstSolveComplete     = nonlinearSolver_->isFirstSolveComplete ();
  }

  return nlInfo;
}

//-----------------------------------------------------------------------------
// Function      : solve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int Manager::solve()
{
  return nonlinearSolver_->solve();
}

//-----------------------------------------------------------------------------
// Function      : setAnalysisMode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
void Manager::setAnalysisMode(AnalysisMode mode)
{
  nonlinearSolver_->setAnalysisMode(mode);
  if (nlsSensitivityPtr_)
  {
    nlsSensitivityPtr_->setAnalysisMode(mode);
  }
}

//-----------------------------------------------------------------------------
// Function      : resetAll
// Purpose       : like setAnalysisMode, but will also result in NOX
//                 resetting a few extra things.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
void Manager::resetAll(AnalysisMode mode)
{
  nonlinearSolver_->resetAll (mode);
}

//-----------------------------------------------------------------------------
// Function      : Manager::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/18/02
//-----------------------------------------------------------------------------
bool Manager::setSensOptions (const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_SENS] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setSensitivityOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/18/02
//-----------------------------------------------------------------------------
bool Manager::setSensitivityOptions (const Util::OptionBlock & option_block)
{
  optionBlockMap_[OPTION_BLOCK_SENSITIVITY] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::enableSensitivity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/30/03
//-----------------------------------------------------------------------------
bool Manager::enableSensitivity(
  TimeIntg::DataStore & data_store,  
  TimeIntg::StepErrorControl & sec,
  Parallel::Manager &   parallel_manager,
  Topo::Topology &      topology,
  IO::OutputMgr & output_manager,
  int & numSensParamsTmp)
{
  bool bsuccess = true;

  if (!nlsSensitivityPtr_)
  {
    Stats::StatTop _sensitivityStat("Setup");

    bool b1 = setupSensitivity(
        data_store, sec,
        parallel_manager, topology, output_manager, numSensParamsTmp);
    bsuccess = bsuccess && b1;
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Manager::icSensitivity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Manager::icSensitivity (
  std::vector<double> & objectiveVec,
  std::vector<double> & dOdpVec,
  std::vector<double> & dOdpAdjVec,
  std::vector<double> & scaled_dOdpVec,
  std::vector<double> & scaled_dOdpAdjVec)
{
  if (!nlsSensitivityPtr_)
  {
    Report::DevelFatal0().in("Manager::isSensitivity") <<  "Manager::enableSensitivity must be called first";
    return false;
  }

  return nlsSensitivityPtr_->icSensitivity(objectiveVec, dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec);
}

//-----------------------------------------------------------------------------
// Function      : Manager::calcSensitivity
// Purpose       : This is the controller function for performing a direct
//                 sensitivity calculation.  It is generally called from
//                 the time integration package, as only that package really
//                 knows when to do it... (I may change this later.)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
bool Manager::calcSensitivity
(
  std::vector<double> & objectiveVec,
  std::vector<double> & dOdpVec,
  std::vector<double> & dOdpAdjVec,
  std::vector<double> & scaled_dOdpVec,
  std::vector<double> & scaled_dOdpAdjVec)
{
  bool bsuccess = true;

  if (!nlsSensitivityPtr_)
  {
    Report::DevelFatal0().in("Manager::calcSensitivity") <<  "Manager::enableSensitivity must be called first";
    return false;
  }

  bsuccess = nlsSensitivityPtr_->solve(objectiveVec, dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Manager::calcTransientAdjoint
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Manager::calcTransientAdjoint (bool timePoint,
      std::vector<double> & objectiveVec,
      std::vector<double> & dOdpVec, std::vector<double> & dOdpAdjVec,
      std::vector<double> & scaled_dOdpVec, std::vector<double> & scaled_dOdpAdjVec)
{
  bool bsuccess = true;

  if (!nlsSensitivityPtr_)
  {
    Report::DevelFatal0().in("Manager::calcTransientAdjiont") <<  "Manager::enableSensitivity must be called first";
    return false;
  }

  bsuccess = nlsSensitivityPtr_->solveTransientAdjoint(timePoint,
      objectiveVec, dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setupSensitivity
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
bool Manager::setupSensitivity(
  TimeIntg::DataStore & data_store,  
  TimeIntg::StepErrorControl & sec,
  Parallel::Manager &   parallel_manager,
  Topo::Topology &      topology,
  IO::OutputMgr & output_manager,
  int & numSensParams_tmp)
{
  if (nlsSensitivityPtr_)
  {
    Report::DevelFatal0().in("Manager::setupSensitivity") <<  "Manager::enableSensitivity may only be called once";
    return false;
  }

  bool bsuccess = true;
  bool bs1 = true;

  nlsSensitivityPtr_ = new Sensitivity(nonlinearSolver_, topology, commandLine_,sec);
  bs1 = nlsSensitivityPtr_->registerExpressionGroup(expressionGroup_); bsuccess = bsuccess && bs1;
  bs1 = nlsSensitivityPtr_->registerParallelMgr(&parallel_manager); bsuccess = bsuccess && bs1;
  bs1 = nlsSensitivityPtr_->registerTIADataStore(&data_store); bsuccess = bsuccess && bs1;
  bs1 = nlsSensitivityPtr_->registerOutputMgr(&output_manager); bsuccess = bsuccess && bs1;

  OptionBlockMap::const_iterator it = optionBlockMap_.find(OPTION_BLOCK_SENS) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = it->second;
    bs1 = nlsSensitivityPtr_->setOptions(option_block);
    bsuccess = bsuccess && bs1;
  }

  it = optionBlockMap_.find(OPTION_BLOCK_SENSITIVITY) ;
  if ( it != optionBlockMap_.end() )
  {
    const Util::OptionBlock &option_block = it->second;
    bs1 = nlsSensitivityPtr_->setSensitivityOptions(option_block);
    bsuccess = bsuccess && bs1;
  }

  bs1 = nlsSensitivityPtr_->doAllocations();
  bsuccess = bsuccess && bs1;

  numSensParams_tmp = nlsSensitivityPtr_->getNumSensParams ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setReturnCodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
void Manager::setReturnCodes(const ReturnCodes & retCodeTmp)
{
  retCodes_ = retCodeTmp;
  if (initializeAllFlag_)
  {
    nonlinearSolver_->setReturnCodes(retCodes_);
  }
}

//-----------------------------------------------------------------------------
// Function      : Manager::getReturnCodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and MEMS Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
const ReturnCodes &Manager::getReturnCodes() const
{
  // the ReturnCodes structure is very simple, so just
  // return a copy of it.
  return retCodes_;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setICOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool Manager::setICOptions (const Util::OptionBlock& option_block )
{
  optionBlockMap_[OPTION_BLOCK_IC] = option_block;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Manager::setNodeSetOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool Manager::setNodeSetOptions (const Util::OptionBlock& option_block )
{
  optionBlockMap_[OPTION_BLOCK_NODESET] = option_block;
  return true;
}

namespace {

void
populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("NONLIN-NLP");

    parameters.insert(Util::ParamMap::value_type("NLSTRATEGY", Util::Param("NLSTRATEGY", 0)));
    parameters.insert(Util::ParamMap::value_type("SEARCHMETHOD", Util::Param("SEARCHMETHOD", 0)));
    parameters.insert(Util::ParamMap::value_type("NOX", Util::Param("NOX", 1)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-3)));
    parameters.insert(Util::ParamMap::value_type("DELTAXTOL", Util::Param("DELTAXTOL", 1.0)));
    parameters.insert(Util::ParamMap::value_type("SMALLUPDATETOL", Util::Param("SMALLUPDATETOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("RHSTOL", Util::Param("RHSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("MAXSTEP", Util::Param("MAXSTEP", 200)));
    parameters.insert(Util::ParamMap::value_type("MAXSEARCHSTEP", Util::Param("MAXSEARCHSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("IN_FORCING", Util::Param("IN_FORCING", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_TOL", Util::Param("AZ_TOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("MATRIXMARKET", Util::Param("MATRIXMARKET", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 99999999)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 1.0E99)));
    parameters.insert(Util::ParamMap::value_type("SCREENOUTPUT", Util::Param("SCREENOUTPUT", 0)));
    parameters.insert(Util::ParamMap::value_type("USEMASKING", Util::Param("USEMASKING", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEPTYPE", Util::Param("RECOVERYSTEPTYPE", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEP", Util::Param("RECOVERYSTEP", 1.0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATION", Util::Param("CONTINUATION", 0)));
    parameters.insert(Util::ParamMap::value_type("ENFORCEDEVICECONV", Util::Param("ENFORCEDEVICECONV", 1)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("NONLIN-TRAN");

    parameters.insert(Util::ParamMap::value_type("NLSTRATEGY", Util::Param("NLSTRATEGY", 0)));
    parameters.insert(Util::ParamMap::value_type("SEARCHMETHOD", Util::Param("SEARCHMETHOD", 0)));
    parameters.insert(Util::ParamMap::value_type("NOX", Util::Param("NOX", 1)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("DELTAXTOL", Util::Param("DELTAXTOL", 0.33)));
    parameters.insert(Util::ParamMap::value_type("SMALLUPDATETOL", Util::Param("SMALLUPDATETOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("RHSTOL", Util::Param("RHSTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("MAXSTEP", Util::Param("MAXSTEP", 20)));
    parameters.insert(Util::ParamMap::value_type("MAXSEARCHSTEP", Util::Param("MAXSEARCHSTEP", 2)));
    parameters.insert(Util::ParamMap::value_type("IN_FORCING", Util::Param("IN_FORCING", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_TOL", Util::Param("AZ_TOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("MATRIXMARKET", Util::Param("MATRIXMARKET", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 99999999)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 1.0E99)));
    parameters.insert(Util::ParamMap::value_type("SCREENOUTPUT", Util::Param("SCREENOUTPUT", 0)));
    parameters.insert(Util::ParamMap::value_type("USEMASKING", Util::Param("USEMASKING", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEPTYPE", Util::Param("RECOVERYSTEPTYPE", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEP", Util::Param("RECOVERYSTEP", 1.0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATION", Util::Param("CONTINUATION", 0)));
    parameters.insert(Util::ParamMap::value_type("ENFORCEDEVICECONV", Util::Param("ENFORCEDEVICECONV", 0)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("NONLIN-HB");

    parameters.insert(Util::ParamMap::value_type("NLSTRATEGY", Util::Param("NLSTRATEGY", 0)));
    parameters.insert(Util::ParamMap::value_type("SEARCHMETHOD", Util::Param("SEARCHMETHOD", 0)));
    parameters.insert(Util::ParamMap::value_type("NOX", Util::Param("NOX", 0)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-9)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("DELTAXTOL", Util::Param("DELTAXTOL", 0.33)));
    parameters.insert(Util::ParamMap::value_type("SMALLUPDATETOL", Util::Param("SMALLUPDATETOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("RHSTOL", Util::Param("RHSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("MAXSTEP", Util::Param("MAXSTEP", 200)));
    parameters.insert(Util::ParamMap::value_type("MAXSEARCHSTEP", Util::Param("MAXSEARCHSTEP", 2)));
    parameters.insert(Util::ParamMap::value_type("IN_FORCING", Util::Param("IN_FORCING", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_TOL", Util::Param("AZ_TOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("MATRIXMARKET", Util::Param("MATRIXMARKET", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 99999999)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 1.0E99)));
    parameters.insert(Util::ParamMap::value_type("SCREENOUTPUT", Util::Param("SCREENOUTPUT", 0)));
    parameters.insert(Util::ParamMap::value_type("USEMASKING", Util::Param("USEMASKING", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEPTYPE", Util::Param("RECOVERYSTEPTYPE", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEP", Util::Param("RECOVERYSTEP", 1.0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATION", Util::Param("CONTINUATION", 0)));
    parameters.insert(Util::ParamMap::value_type("ENFORCEDEVICECONV", Util::Param("ENFORCEDEVICECONV", 0)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("NONLIN-TWOLEVEL");

    parameters.insert(Util::ParamMap::value_type("NLSTRATEGY", Util::Param("NLSTRATEGY", 0)));
    parameters.insert(Util::ParamMap::value_type("SEARCHMETHOD", Util::Param("SEARCHMETHOD", 0)));
    parameters.insert(Util::ParamMap::value_type("NOX", Util::Param("NOX", 1)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-3)));
    parameters.insert(Util::ParamMap::value_type("DELTAXTOL", Util::Param("DELTAXTOL", 1.0)));
    parameters.insert(Util::ParamMap::value_type("SMALLUPDATETOL", Util::Param("SMALLUPDATETOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("RHSTOL", Util::Param("RHSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("MAXSTEP", Util::Param("MAXSTEP", 200)));
    parameters.insert(Util::ParamMap::value_type("MAXSEARCHSTEP", Util::Param("MAXSEARCHSTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("IN_FORCING", Util::Param("IN_FORCING", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_TOL", Util::Param("AZ_TOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("MATRIXMARKET", Util::Param("MATRIXMARKET", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 99999999)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 1.0E99)));
    parameters.insert(Util::ParamMap::value_type("SCREENOUTPUT", Util::Param("SCREENOUTPUT", 0)));
    parameters.insert(Util::ParamMap::value_type("USEMASKING", Util::Param("USEMASKING", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEPTYPE", Util::Param("RECOVERYSTEPTYPE", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEP", Util::Param("RECOVERYSTEP", 1.0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATION", Util::Param("CONTINUATION", 0)));
    parameters.insert(Util::ParamMap::value_type("ENFORCEDEVICECONV", Util::Param("ENFORCEDEVICECONV", 0)));
    parameters.insert(Util::ParamMap::value_type("ALGORITHM", Util::Param("ALGORITHM", 0)));
    parameters.insert(Util::ParamMap::value_type("MAXCONTSTEPS", Util::Param("MAXCONTSTEPS", 0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATIONFLAG", Util::Param("CONTINUATIONFLAG", 1)));
    parameters.insert(Util::ParamMap::value_type("INNERFAIL", Util::Param("INNERFAIL", 1)));
    parameters.insert(Util::ParamMap::value_type("EXITWITHFAILURE", Util::Param("EXITWITHFAILURE", 1)));
    parameters.insert(Util::ParamMap::value_type("FULLNEWTONENFORCE", Util::Param("FULLNEWTONENFORCE", 1)));
    parameters.insert(Util::ParamMap::value_type("CONPARAM", Util::Param("CONPARAM", "VA:V0")));
    parameters.insert(Util::ParamMap::value_type("VOLTLIMTOL", Util::Param("VOLTLIMTOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("REUSEFACTORS", Util::Param("REUSEFACTORS",1)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("NONLIN-TWOLEVEL-TRAN");

    parameters.insert(Util::ParamMap::value_type("NLSTRATEGY", Util::Param("NLSTRATEGY", 0)));
    parameters.insert(Util::ParamMap::value_type("SEARCHMETHOD", Util::Param("SEARCHMETHOD", 0)));
    parameters.insert(Util::ParamMap::value_type("NOX", Util::Param("NOX", 1)));
    parameters.insert(Util::ParamMap::value_type("ABSTOL", Util::Param("ABSTOL", 1.0E-6)));
    parameters.insert(Util::ParamMap::value_type("RELTOL", Util::Param("RELTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("DELTAXTOL", Util::Param("DELTAXTOL", 0.33)));
    parameters.insert(Util::ParamMap::value_type("SMALLUPDATETOL", Util::Param("SMALLUPDATETOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("RHSTOL", Util::Param("RHSTOL", 1.0E-2)));
    parameters.insert(Util::ParamMap::value_type("MAXSTEP", Util::Param("MAXSTEP", 20)));
    parameters.insert(Util::ParamMap::value_type("MAXSEARCHSTEP", Util::Param("MAXSEARCHSTEP", 2)));
    parameters.insert(Util::ParamMap::value_type("IN_FORCING", Util::Param("IN_FORCING", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_TOL", Util::Param("AZ_TOL", 1.0E-12)));
    parameters.insert(Util::ParamMap::value_type("MATRIXMARKET", Util::Param("MATRIXMARKET", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 1)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIMESTEP", Util::Param("DEBUGMINTIMESTEP", 0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIMESTEP", Util::Param("DEBUGMAXTIMESTEP", 99999999)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMINTIME", Util::Param("DEBUGMINTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("DEBUGMAXTIME", Util::Param("DEBUGMAXTIME", 1.0E99)));
    parameters.insert(Util::ParamMap::value_type("SCREENOUTPUT", Util::Param("SCREENOUTPUT", 0)));
    parameters.insert(Util::ParamMap::value_type("USEMASKING", Util::Param("USEMASKING", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEPTYPE", Util::Param("RECOVERYSTEPTYPE", 0)));
    parameters.insert(Util::ParamMap::value_type("RECOVERYSTEP", Util::Param("RECOVERYSTEP", 1.0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATION", Util::Param("CONTINUATION", 0)));
    parameters.insert(Util::ParamMap::value_type("ENFORCEDEVICECONV", Util::Param("ENFORCEDEVICECONV", 0)));
    parameters.insert(Util::ParamMap::value_type("ALGORITHM", Util::Param("ALGORITHM", 0)));
    parameters.insert(Util::ParamMap::value_type("MAXCONTSTEPS", Util::Param("MAXCONTSTEPS", 0)));
    parameters.insert(Util::ParamMap::value_type("CONTINUATIONFLAG", Util::Param("CONTINUATIONFLAG", 1)));
    parameters.insert(Util::ParamMap::value_type("INNERFAIL", Util::Param("INNERFAIL", 1)));
    parameters.insert(Util::ParamMap::value_type("EXITWITHFAILURE", Util::Param("EXITWITHFAILURE", 1)));
    parameters.insert(Util::ParamMap::value_type("FULLNEWTONENFORCE", Util::Param("FULLNEWTONENFORCE", 1)));
    parameters.insert(Util::ParamMap::value_type("CONPARAM", Util::Param("CONPARAM", "VA:V0")));
    parameters.insert(Util::ParamMap::value_type("VOLTLIMTOL", Util::Param("VOLTLIMTOL", 1.0e-6)));
    parameters.insert(Util::ParamMap::value_type("REUSEFACTORS", Util::Param("REUSEFACTORS",1)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("LOCA");

    parameters.insert(Util::ParamMap::value_type("STEPPER", Util::Param("STEPPER", "NATURAL")));
    parameters.insert(Util::ParamMap::value_type("PREDICTOR", Util::Param("PREDICTOR", "CONSTANT")));
    parameters.insert(Util::ParamMap::value_type("STEPCONTROL", Util::Param("STEPCONTROL", "CONSTANT")));

    parameters.insert(Util::ParamMap::value_type("CONPARAM", Util::Param("CONPARAM", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("INITIALVALUE", Util::Param("INITIALVALUE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MAXVALUE", Util::Param("MAXVALUE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MINVALUE", Util::Param("MINVALUE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("INITIALSTEPSIZE", Util::Param("INITIALSTEPSIZE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MAXSTEPSIZE", Util::Param("MAXSTEPSIZE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MINSTEPSIZE", Util::Param("MINSTEPSIZE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("AGGRESSIVENESS", Util::Param("AGGRESSIVENESS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("BIFPARAM", Util::Param("BIFPARAM", "VA:V0")));
    parameters.insert(Util::ParamMap::value_type("MAXSTEPS", Util::Param("MAXSTEPS", 20)));
    parameters.insert(Util::ParamMap::value_type("MAXNLITERS", Util::Param("MAXNLITERS", 20)));
    parameters.insert(Util::ParamMap::value_type("PARAMLIST", Util::Param("PARAMLIST", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("VOLTAGELIST", Util::Param("VOLTAGELIST", "DOFS")));
    parameters.insert(Util::ParamMap::value_type("VOLTAGESCALEFACTOR", Util::Param("VOLTAGESCALEFACTOR", 1.0)));
    parameters.insert(Util::ParamMap::value_type("RESIDUALCONDUCTANCE", Util::Param("RESIDUALCONDUCTANCE", 0.0)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("TWOLEVEL-LOCA");

    parameters.insert(Util::ParamMap::value_type("STEPPER", Util::Param("STEPPER", "NATURAL")));
    parameters.insert(Util::ParamMap::value_type("PREDICTOR", Util::Param("PREDICTOR", "CONSTANT")));
    parameters.insert(Util::ParamMap::value_type("STEPCONTROL", Util::Param("STEPCONTROL", "CONSTANT")));
    parameters.insert(Util::ParamMap::value_type("CONPARAM", Util::Param("CONPARAM", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("INITIALVALUE", Util::Param("INITIALVALUE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MAXVALUE", Util::Param("MAXVALUE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MINVALUE", Util::Param("MINVALUE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("INITIALSTEPSIZE", Util::Param("INITIALSTEPSIZE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MAXSTEPSIZE", Util::Param("MAXSTEPSIZE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("MINSTEPSIZE", Util::Param("MINSTEPSIZE", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("AGGRESSIVENESS", Util::Param("AGGRESSIVENESS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("BIFPARAM", Util::Param("BIFPARAM", "VA:V0")));
    parameters.insert(Util::ParamMap::value_type("MAXSTEPS", Util::Param("MAXSTEPS", 20)));
    parameters.insert(Util::ParamMap::value_type("MAXNLITERS", Util::Param("MAXNLITERS", 20)));
    parameters.insert(Util::ParamMap::value_type("PARAMLIST", Util::Param("PARAMLIST", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("VOLTAGELIST", Util::Param("VOLTAGELIST", "DOFS")));
    parameters.insert(Util::ParamMap::value_type("VOLTAGESCALEFACTOR", Util::Param("VOLTAGESCALEFACTOR", 1.0)));
    parameters.insert(Util::ParamMap::value_type("RESIDUALCONDUCTANCE", Util::Param("RESIDUALCONDUCTANCE", 0.0)));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("LINSOL");

    parameters.insert(Util::ParamMap::value_type("AZ_max_iter", Util::Param("AZ_max_iter", 200)));
    parameters.insert(Util::ParamMap::value_type("AZ_precond", Util::Param("AZ_precond", 14)));
    parameters.insert(Util::ParamMap::value_type("AZ_solver", Util::Param("AZ_solver", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_conv", Util::Param("AZ_conv", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_pre_calc", Util::Param("AZ_pre_calc", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_keep_info", Util::Param("AZ_keep_info", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_orthog", Util::Param("AZ_orthog", 1)));
    parameters.insert(Util::ParamMap::value_type("AZ_subdomain_solve", Util::Param("AZ_subdomain_solve", 9)));
    parameters.insert(Util::ParamMap::value_type("AZ_ilut_fill", Util::Param("AZ_ilut_fill", 3.0)));
    parameters.insert(Util::ParamMap::value_type("AZ_drop", Util::Param("AZ_drop", 1.0E-3)));
    parameters.insert(Util::ParamMap::value_type("AZ_reorder", Util::Param("AZ_reorder", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_scaling", Util::Param("AZ_scaling", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_kspace", Util::Param("AZ_kspace", 50)));
    parameters.insert(Util::ParamMap::value_type("AZ_tol", Util::Param("AZ_tol", 1.0E-9)));
    parameters.insert(Util::ParamMap::value_type("AZ_output", Util::Param("AZ_output", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_diagnostics", Util::Param("AZ_diagnostics", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_overlap", Util::Param("AZ_overlap", 0)));
    parameters.insert(Util::ParamMap::value_type("AZ_rthresh", Util::Param("AZ_rthresh", 1.0001)));
    parameters.insert(Util::ParamMap::value_type("AZ_athresh", Util::Param("AZ_athresh", 1.0E-4)));
    parameters.insert(Util::ParamMap::value_type("AZ_filter", Util::Param("AZ_filter", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TR_partition", Util::Param("TR_partition", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_partition_type", Util::Param("TR_partition_type", "HYPERGRAPH")));
#ifdef Xyce_SHYLU
    parameters.insert(Util::ParamMap::value_type("ShyLU_rthresh", Util::Param("ShyLU_rthresh", 1.0E-3)));
#endif
    parameters.insert(Util::ParamMap::value_type("TR_reindex", Util::Param("TR_reindex", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_solvermap", Util::Param("TR_solvermap", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_amd", Util::Param("TR_amd", 1)));
    parameters.insert(Util::ParamMap::value_type("TR_btf", Util::Param("TR_btf", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf", Util::Param("TR_global_btf", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf_droptol", Util::Param("TR_global_btf_droptol", 1.0E-16)));
    parameters.insert(Util::ParamMap::value_type("TR_global_btf_verbose", Util::Param("TR_global_btf_verbose", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_amd", Util::Param("TR_global_amd", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_global_amd_verbose", Util::Param("TR_global_amd_verbose", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_singleton_filter", Util::Param("TR_singleton_filter", 0)));
    parameters.insert(Util::ParamMap::value_type("adaptive_solve", Util::Param("adaptive_solve", 0)));
    parameters.insert(Util::ParamMap::value_type("use_aztec_precond", Util::Param("use_aztec_precond", 1)));
    parameters.insert(Util::ParamMap::value_type("use_ifpack_factory", Util::Param("use_ifpack_factory", 0)));
    parameters.insert(Util::ParamMap::value_type("ifpack_type", Util::Param("ifpack_type", "Amesos")));
    parameters.insert(Util::ParamMap::value_type("diag_perturb", Util::Param("diag_perturb", 0.0)));
    parameters.insert(Util::ParamMap::value_type("TR_rcm", Util::Param("TR_rcm", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale", Util::Param("TR_scale", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_left", Util::Param("TR_scale_left", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_right", Util::Param("TR_scale_right", 0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_exp", Util::Param("TR_scale_exp", 1.0)));
    parameters.insert(Util::ParamMap::value_type("TR_scale_iter", Util::Param("TR_scale_iter", 0)));
    parameters.insert(Util::ParamMap::value_type("TYPE", Util::Param("TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("PREC_TYPE", Util::Param("PREC_TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("IR_SOLVER_TYPE", Util::Param("IR_SOLVER_TYPE", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("IR_SOLVER_TOL", Util::Param("IR_SOLVER_TOL", "DEFAULT")));
    parameters.insert(Util::ParamMap::value_type("BELOS_SOLVER_TYPE", Util::Param("BELOS_SOLVER_TYPE", "Block GMRES")));
    parameters.insert(Util::ParamMap::value_type("KLU_REPIVOT", Util::Param("KLU_REPIVOT", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_LS", Util::Param("OUTPUT_LS", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_BASE_LS", Util::Param("OUTPUT_BASE_LS", 1)));
    parameters.insert(Util::ParamMap::value_type("OUTPUT_FAILED_LS", Util::Param("OUTPUT_FAILED_LS", 1)));
  }
}

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool registerPkgOptionsMgr(Manager &manager, IO::PkgOptionsMgr &options_manager)
{
  populateMetadata(options_manager);
  Nonlinear::Sensitivity::populateMetadata(options_manager);
  Nonlinear::NLParams::populateMetadata(options_manager);
  
  options_manager.addCommandProcessor("SENS", IO::createRegistrationOptions(manager, &Manager::setSensOptions));
  options_manager.addCommandProcessor("IC", IO::createRegistrationOptions(manager, &Manager::setICOptions));
  options_manager.addCommandProcessor("NODESET", IO::createRegistrationOptions(manager, &Manager::setNodeSetOptions));

  options_manager.addOptionsProcessor("NONLIN", IO::createRegistrationOptions(manager, &Manager::setOptions));
  options_manager.addOptionsProcessor("NONLIN-TRAN", IO::createRegistrationOptions(manager, &Manager::setTranOptions));
  options_manager.addOptionsProcessor("NONLIN-NLP", IO::createRegistrationOptions(manager, &Manager::setNLPOptions));
  options_manager.addOptionsProcessor("NONLIN-HB", IO::createRegistrationOptions(manager, &Manager::setHBOptions));
  options_manager.addOptionsProcessor("LINSOL", IO::createRegistrationOptions(manager, &Manager::setLinSolOptions));
  options_manager.addOptionsProcessor("LOCA", IO::createRegistrationOptions(manager, &Manager::setLocaOptions));
  options_manager.addOptionsProcessor("SENSITIVITY", IO::createRegistrationOptions(manager, &Manager::setSensitivityOptions));
  options_manager.addOptionsProcessor("NONLIN-TWOLEVEL", IO::createRegistrationOptions(manager, &Manager::setTwoLevelOptions));
  options_manager.addOptionsProcessor("NONLIN-TWOLEVEL-TRAN", IO::createRegistrationOptions(manager, &Manager::setTwoLevelTranOptions));

  return true;
}

} // namespace Nonlinear
} // namespace Xyce

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

//-----------------------------------------------------------------------------
// Purpose       : MOR analysis functions.
// Special Notes :
// Creator       : Ting Mei
// Creation Date :  7/11
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iomanip>


// ----------   Xyce Includes   ----------

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_MOR.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ANP_Report.h>

#include <N_LAS_SystemHelpers.h>
#include <N_LAS_BlockSystemHelpers.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_MOROperators.h>
#include <N_LAS_EpetraMatrix.h>
#include <N_LAS_EpetraBlockMatrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_System.h>
#include <N_LAS_Graph.h>

#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_OutputROM.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_LOA_Loader.h>

#include <N_NLS_Manager.h>
#include <N_NLS_ReturnCodes.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_PDS_ParHelpers.h>
#include <N_PDS_EpetraParMap.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_NoTimeIntegration.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_TIA_StepErrorControl.h>

#include <N_TOP_Topology.h>

#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Factory.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_Timer.h>

// ----------   Other Includes   ----------

#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LinearProblem.h>
#include <Amesos.h>

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresIter.hpp>
#include <BelosDGKSOrthoManager.hpp>
#include <BelosStatusTestMaxIters.hpp>
#include <BelosOutputManager.hpp>
#include <BelosEpetraAdapter.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_BLAS.hpp>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : MOR::MOR( AnalysisManager * )
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
MOR::MOR(
  AnalysisManager &                     analysis_manager,
  Linear::System &                      linear_system,
  Nonlinear::Manager &                  nonlinear_manager,
  Loader::Loader &                      loader,
  Topo::Topology &                      topology,
  IO::InitialConditionsManager &        initial_conditions_manager)
  : AnalysisBase(analysis_manager, "MOR"),
    comm_(analysis_manager.getPDSManager()->getPDSComm()->comm()),
    analysisManager_(analysis_manager),
    loader_(loader),
    linearSystem_(linear_system),
    nonlinearManager_(nonlinear_manager),
    topology_(topology),
    initialConditionsManager_(initial_conditions_manager),
    outputManagerAdapter_(analysis_manager.getOutputManagerAdapter()),
    outputMOR_(analysisManager_.getNetlistFilename()),
    ROMsize_(-1),
    morMethod_("PRIMA"),
    morSaveRedSys_(false),
    morCompOrigTF_(false),
    morCompRedTF_(false),
    morCompType_("DEC"),
    morCompNP_(10),
    morCompFStart_(1.0),
    morAutoSize_(false),
    morMaxSize_(-1),
    morMaxFreq_(1.0e9),
    morCompFStop_(1.0),
    morExpPoint_(0.0),
    morScaleFactor_(1.0),
    morScaleType_(0),
    morScaleFactor1_(1.0),
    morSparsificationType_(0),
    isROMSparse_(false),
    morEvalSize_(0),
    numPorts_(0),
    isSingleFreq_(false),
    stepMult_(0.0),
    fStep_(0.0),
    currentFreq_(0.0),
    s0_(0.0)
{}

//-----------------------------------------------------------------------------
// Function      : MOR::~MOR()
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Ting Mei
// Creation Date : 5/11
//-----------------------------------------------------------------------------
MOR::~MOR()
{}

//-----------------------------------------------------------------------------
// Function      : MOR::setAnalysisParams
// Purpose       :
// Special Notes : These are from the .MOR statement.
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 6/11
//-----------------------------------------------------------------------------
bool MOR::setAnalysisParams(const Util::OptionBlock & paramsBlock)
{
  for (Util::ParamList::const_iterator it = paramsBlock.begin(), end = paramsBlock.end(); it != end; ++it)
  {
    if ((*it).uTag() == "PORTLIST")
    {
      portList_ = (*it).getValue<std::vector<std::string> >();
    }
  }

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    dout() << std::endl
           << section_divider << std::endl
           <<" MOR simulation parameters" << std::endl
           << " size = " << ROMsize_ << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::setMOROptions
// Purpose       :
// Special Notes : These are from '.options mor'
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 05/31/12
//-----------------------------------------------------------------------------
bool
MOR::setMOROptions(
  const Util::OptionBlock &     option_block)
{
  for (Util::ParamList::const_iterator it = option_block.begin(), end = option_block.end(); it != end; ++it)
  {
    bool value_set =
      Util::setValue(*it, "METHOD", morMethod_)
      || Util::setValue(*it, "SAVEREDSYS", morSaveRedSys_)
      || Util::setValue(*it, "COMPORIGTF", morCompOrigTF_)
      || Util::setValue(*it, "COMPREDTF", morCompRedTF_)
      || Util::setValue(*it, "COMPTYPE", morCompType_)
      || Util::setValue(*it, "COMPNP", morCompNP_)
      || Util::setValue(*it, "COMPFSTART", morCompFStart_)
      || Util::setValue(*it, "AUTOSIZE", morAutoSize_)
      || Util::setValue(*it, "MAXSIZE", morMaxSize_)
      || Util::setValue(*it, "MAXFREQ", morMaxFreq_)
      || Util::setValue(*it, "SIZE", ROMsize_)
      || Util::setValue(*it, "COMPFSTOP", morCompFStop_)
      || Util::setValue(*it, "EXPPOINT", morExpPoint_)
      || Util::setValue(*it, "SCALETYPE", morScaleType_)
      || Util::setValue(*it, "SCALEFACTOR", morScaleFactor_)
      || Util::setValue(*it, "SCALEFACTOR1", morScaleFactor1_)
      || Util::setValue(*it, "SPARSIFICATIONTYPE", morSparsificationType_)
      || Util::setValue(*it, "SUBCKTS", subcircuitNames_);

    if (value_set)
      ;
    else
    {
      Report::UserError0() << (*it).uTag() << " is not a recognized model-order reduction option.";
    }
  }

  // If we are computing the transfer function, make sure the frequency range is valid.
  if (morCompOrigTF_ || morCompRedTF_)
  {
    if (morCompFStop_ < morCompFStart_)
    {
      Report::UserError() << ".options mor COMPFSTART = " << morCompFStart_ << " > " << morCompFStop_ << " = COMPFSTOP!";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::run()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool MOR::doRun()
{
  bool bsuccess = doInit() && reduceSystem();

  if (morCompOrigTF_)
  {
    // Evaluate original system
    bsuccess = bsuccess && evalOrigTransferFunction();
  }

  // Reset the output adapter, just in case the reduced system needs to output
  // its transfer functions.
  // outputManagerAdapter_.resetOutputMORTF();
  outputMOR_.reset();

  if (morCompRedTF_)
  {
    // Evaluate reduced system
    bsuccess = bsuccess && evalRedTransferFunction();
  }

  bsuccess = bsuccess && doFinish();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::init()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool MOR::doInit()
{
  bool bsuccess = true;

  // Compute expansion point in transformed space.
  s0_ =  2.0 * M_PI * morExpPoint_;

  if (morCompOrigTF_ || morCompRedTF_)
  {
    morEvalSize_ = setupSweepParam_();
  }

  // Get set to do the operating point.
  baseIntegrationMethod_ = TimeIntg::NoTimeIntegration::type;
  analysisManager_.createTimeIntegratorMethod(tiaParams_, baseIntegrationMethod_);

  stepNumber            = 0;
  setDoubleDCOPEnabled(loader_.isPDESystem());

  // set initial guess, if there is one to be set.
  // this setInitialGuess call is to up an initial guess in the
  // devices that have them (usually PDE devices).  This is different than
  // the "intializeProblem" call, which sets IC's.  (initial conditions are
  // different than initial guesses.
  loader_.setInitialGuess (analysisManager_.getDataStore()->nextSolutionPtr);

  // If available, set initial solution
  // setInputOPFlag(
  //   outputManagerAdapter_.setupInitialConditions(
  //     *analysisManager_.getDataStore()->nextSolutionPtr,
  //     *linearSystem_.getFlagSolVector()));

  setInputOPFlag(
    initialConditionsManager_.setupInitialConditions(outputManagerAdapter_.getComm(),
                                                     topology_.getSolutionNodeNameMap(),
                                                     outputManagerAdapter_.getAliasNodeMap(),
                                                     *analysisManager_.getDataStore()->nextSolutionPtr,
                                                     *linearSystem_.getFlagSolVector()));

  // Set a constant history for operating point calculation
  analysisManager_.getDataStore()->setConstantHistory();
  analysisManager_.getWorkingIntegrationMethod().obtainCorrectorDeriv();

  // solving for DC op
  doHandlePredictor();
  loader_.updateSources();
  analysisManager_.getStepErrorControl().newtonConvergenceStatus = nonlinearManager_.solve();
  analysisManager_.getWorkingIntegrationMethod().stepLinearCombo ();
  gatherStepStatistics(stats_, nonlinearManager_.getNonlinearSolver(), analysisManager_.getStepErrorControl().newtonConvergenceStatus);
  analysisManager_.getStepErrorControl().evaluateStepError(loader_, tiaParams_);

  if ( analysisManager_.getStepErrorControl().newtonConvergenceStatus <= 0)
  {
    Report::UserError() << "Solving for DC operating point failed, cannot continue MOR analysis";
    return false;
  }

  // Create B matrix stamp
  std::vector<int> tempVec;
  loader_.getBMatrixEntries(tempVec, bMatPosEntriesVec_);

  // Determine how many global ports there are by performing a sumAll()
  int hsize = tempVec.size();
  Parallel::AllReduce(comm_, MPI_SUM, &hsize, &numPorts_, 1);

  // Check that the listed ports on the .mor line are the same number as the B matrix entries.
  if ( numPorts_ != (int)portList_.size() )
  {
    Report::UserError() << "Number of specified ports in .MOR line is inconsistent with number of voltage sources";
    return false;
  }

  // Getting the GIDs for the port list in the order specified in the .MOR line
  // This is used to permute the bMatEntriesVec_ in the order specified by the user.
  std::vector<int> gidPosEntries( hsize );
  for (int i=0; i<hsize; ++i)
  {
    std::vector<int> svGIDList1, dummyList;
    char type1;
    topology_.getNodeSVarGIDs(NodeID(portList_[i],_VNODE), svGIDList1, dummyList, type1);

    // Grab the GID for this port.
    gidPosEntries[i] = svGIDList1.front();
  }

  // Use the base map to get the global IDs
  // Get the parallel maps for the original system
  Parallel::Manager &pdsManager = *analysisManager_.getPDSManager();
  RCP<Parallel::ParMap> BaseMap = rcp(pdsManager.getParallelMap( Parallel::SOLUTION ), false);

  // Find the voltage source corresponding to each port and place the LID in the bMatEntriesVec_.
  bMatEntriesVec_.resize( hsize );
  for (int i=0; i<hsize; ++i)
  {
    int gid = gidPosEntries[i];
    bool found = false;
    for (int j=0; j<hsize; ++j)
    {
      if (gid == BaseMap->localToGlobalIndex(bMatPosEntriesVec_[j]))
      {
        bMatEntriesVec_[i] = tempVec[j];
        found = true;
        break;
      }
    }
    if (!found)
    {
      Report::UserError() << "Did not find voltage source corresponding to port";
      return false;
    }
  }

  if (morCompOrigTF_)
  {
    // Resize transfer function matrices
    origH_.shape(numPorts_, numPorts_);
  }
  if (morCompRedTF_)
  {
    // Resize transfer function matrices
    redH_.shape(numPorts_, numPorts_);
  }

  // Create C and G matrices from DCOP solution
  analysisManager_.getDataStore()->daeQVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->daeFVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->daeBVectorPtr->putScalar(0.0);

  analysisManager_.getDataStore()->dFdxdVpVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->dQdxdVpVectorPtr->putScalar(0.0);
  analysisManager_.getDataStore()->dQdxMatrixPtr->put(0.0);
  analysisManager_.getDataStore()->dFdxMatrixPtr->put(0.0);

  loader_.updateState
                ((analysisManager_.getDataStore()->nextSolutionPtr),
               (analysisManager_.getDataStore()->currSolutionPtr),
               (analysisManager_.getDataStore()->lastSolutionPtr),
               (analysisManager_.getDataStore()->nextStatePtr),
               (analysisManager_.getDataStore()->currStatePtr),
               (analysisManager_.getDataStore()->lastStatePtr),
               (analysisManager_.getDataStore()->nextStorePtr),
               (analysisManager_.getDataStore()->currStorePtr)
               );

  loader_.loadDAEVectors
              ((analysisManager_.getDataStore()->nextSolutionPtr),
               (analysisManager_.getDataStore()->currSolutionPtr),
               (analysisManager_.getDataStore()->lastSolutionPtr),
               (analysisManager_.getDataStore()->nextStatePtr),
               (analysisManager_.getDataStore()->currStatePtr),
               (analysisManager_.getDataStore()->lastStatePtr),
               (analysisManager_.getDataStore()->nextStateDerivPtr),
               (analysisManager_.getDataStore()->nextStorePtr),
               (analysisManager_.getDataStore()->currStorePtr),
               (analysisManager_.getDataStore()->nextLeadCurrentPtr),
               (analysisManager_.getDataStore()->nextLeadCurrentQPtr),
               (analysisManager_.getDataStore()->nextLeadDeltaVPtr),
               (analysisManager_.getDataStore()->daeQVectorPtr),
               (analysisManager_.getDataStore()->daeFVectorPtr),
               (analysisManager_.getDataStore()->daeBVectorPtr),
               (analysisManager_.getDataStore()->dFdxdVpVectorPtr),
               (analysisManager_.getDataStore()->dQdxdVpVectorPtr) );

  loader_.loadDAEMatrices(analysisManager_.getDataStore()->nextSolutionPtr,
      analysisManager_.getDataStore()->nextStatePtr, analysisManager_.getDataStore()->nextStateDerivPtr,
      analysisManager_.getDataStore()->nextStorePtr,
      analysisManager_.getDataStore()->dQdxMatrixPtr,  analysisManager_.getDataStore()->dFdxMatrixPtr);

  CPtr_ = rcp(dynamic_cast<Linear::EpetraMatrix *>(analysisManager_.getDataStore()->dQdxMatrixPtr), false);
  GPtr_ = rcp(dynamic_cast<Linear::EpetraMatrix *>(analysisManager_.getDataStore()->dFdxMatrixPtr), false);

  ///  Xyce::dout() << "Branch nodes: " << std::endl;
  ///  for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
  ///  {
  ///    Xyce::dout() << "Node " << i << " : " << bMatEntriesVec_[i] << std::endl;
  ///  }

  ///  Xyce::dout() << "Printing GPtr: " << std::endl;
  ///  GPtr_->print(Xyce::dout());
  ///  Xyce::dout() << "Printing CPtr: " << std::endl;
  ///  CPtr_->print(Xyce::dout());


  // Storage for row extraction
  int length=1, numEntries=0;
  std::vector<int> colIndices(length);
  std::vector<double> coeffs(length);

  for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
  {
     // Get the number of non-zero entries in this row.
     numEntries = GPtr_->getLocalRowLength( bMatEntriesVec_[i] );
     if ( numEntries != 1 )
     {
       Report::UserError0() << "Supposed voltage source row has too many entries, cannot continue MOR analysis";
       return false;
     }

     // Extract out rows of G based on indices in bMatEntriesVec_.
     GPtr_->getLocalRowCopy(bMatEntriesVec_[i], length, numEntries, &coeffs[0], &colIndices[0]);

     // If the coefficient for this voltage source is positive, make it negative.
     if ( coeffs[0] > 0.0 )
     {
       coeffs[0] *= -1.0;
       GPtr_->putLocalRow(bMatEntriesVec_[i], length, &coeffs[0], &colIndices[0]);
     }
  }


  ///  Xyce::dout() << "Printing GPtr (after scaling): " << std::endl;
  ///  GPtr_->print(Xyce::dout());

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::reduceSystem()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 6/4/2012
//-----------------------------------------------------------------------------
bool MOR::reduceSystem()
{
  // At this time, the reduceSystem() method will compute the reduced-order system using
  // multi-port PRIMA.  This method should allow the user multiple algorithms, defined by
  // the morMethod_ option, which would be best implemented using a factory design.

  bool bsuccess = true;

  // Get the parallel maps for the original system
  Parallel::Manager &pdsManager = *analysisManager_.getPDSManager();

  Parallel::ParMap &BaseMap = *pdsManager.getParallelMap( Parallel::SOLUTION );

  // Create matrix for (G + s0 * C)
  if (s0_ != 0.0)
  {
    sCpG_MatrixPtr_ = rcp( new Linear::EpetraMatrix( &CPtr_->epetraObj(), false ) );
    sCpG_MatrixPtr_->scale( s0_ );
    sCpG_MatrixPtr_->add( *GPtr_ );
  }
  else
  {
    sCpG_MatrixPtr_ = rcp( new Linear::EpetraMatrix( &GPtr_->epetraObj(), false ) );
  }

  // Create multivector for B and R
  RPtr_ = rcp( Xyce::Linear::createMultiVector( BaseMap, numPorts_ ) );
  RPtr_->putScalar( 0.0 );
  BPtr_ = rcp( Xyce::Linear::createMultiVector( BaseMap, numPorts_ ) );
  for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
  {
    BPtr_->setElementByGlobalIndex( BaseMap.localToGlobalIndex(bMatEntriesVec_[j]), -1.0, j );
  }

  ///  Xyce::dout() << "Printing out BPtr" << std::endl;
  ///  BPtr_->print(Xyce::dout());

  ///  Xyce::dout() << "Printing out sCpG" << std::endl;
  ///  sCpG_MatrixPtr_->print(Xyce::dout());

  // Create linear problem for (G + s0 * C)
  origProblem_ = rcp(new Epetra_LinearProblem(&sCpG_MatrixPtr_->epetraObj(), &RPtr_->epetraObj(), &BPtr_->epetraObj() ) );

  // Create solver object for this linear problem, which will be used to generate the projection basis
  Amesos amesosFactory;
  origSolver_ = rcp( amesosFactory.Create( "Klu", *origProblem_ ) );

  // Solve for R = inv(G + s0*C)*B
  // Perform symbolic factorization.
  int linearStatus = origSolver_->SymbolicFactorization();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  // Perform numeric factorization
  linearStatus = origSolver_->NumericFactorization();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  // Perform solve for R = inv(G + s0*C)*B
  linearStatus = origSolver_->Solve();
  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos solve exited with error: " << linearStatus << std::endl;
    bsuccess = false;
  }

  // Create an Epetra_Operator object to apply the operator inv(G + s0*C)*C
  RCP<Epetra_Operator> COpPtr_ = rcp( &CPtr_->epetraObj(), false );
  RCP<Linear::AmesosGenOp> AOp = rcp( new Linear::AmesosGenOp( origSolver_, COpPtr_ ) );

  // Check to see if the requested size of the ROM is valid.
  // An orthogonal basis cannot be generated that is larger than the dimension of the original system
  int kblock = 0;

  if (morAutoSize_)
  {
    if (morMaxSize_ != -1)
      ROMsize_ =  morMaxSize_;
    else
      ROMsize_ = RPtr_->globalLength()/4;
  }
  else
  {

    if (ROMsize_ == -1)
    {
      Report::UserError() << "Automatic Sizing is OFF. Please specify the ROM dimension";
    }
  }


  if (ROMsize_ > RPtr_->globalLength())
  {
    kblock = (int)(RPtr_->globalLength() / numPorts_);
    UserWarning(*this) << "Requested reduced-order model dimension is larger than original system dimension, resizing to original system dimension";
  }
  else
  {
    kblock = (int)(ROMsize_ / numPorts_);
  }
  int k = kblock * numPorts_;

  // Resize the projection matrices
  redG_.shape(k, k);
  redC_.shape(k, k);
  redB_.shape(k, numPorts_);
  redL_.shape(k, numPorts_);

  // ---------------------------------------------------------------------
  // Now use Belos to compute the basis vectors for K_k(inv(G + s0*C)*C, R)
  // ---------------------------------------------------------------------

  // Helpful typedefs for the templates
  typedef double                            ST;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>      MVT;

  // Output manager.
  Belos::OutputManager<ST> printer;

  // Status test.
  Belos::StatusTestMaxIters<ST, MV, OP> maxIterTest( kblock );

  // Orthogonalization manager.
  Belos::DGKSOrthoManager<ST, MV, OP> orthoMgr;

  // Linear Problem.
  // Reuse RPtr_.  We need the basis vectors for K(inv(G + s0*C)*C, R)
  Teuchos::RCP<Linear::MultiVector> temp =
    Teuchos::rcp( Xyce::Linear::createMultiVector( BaseMap, numPorts_ ) );
  Belos::LinearProblem<ST, MV, OP > problem( AOp,
                                             rcp( &temp->epetraObj(), false ),
                                             rcp( &RPtr_->epetraObj(), false ));
  problem.setProblem();

  // Create parameter list.
  Teuchos::ParameterList params;
  params.set("Num Blocks", kblock);
  params.set("Block Size", numPorts_ );

  // Create Krylov subspace iteration from Belos
  Belos::BlockGmresIter<ST, MV, OP> krylovIter( rcp( &problem, false ),
                                                rcp( &printer, false ),
                                                rcp( &maxIterTest, false ),
                                                rcp( &orthoMgr, false ),
                                                params );

  // Get a matrix to hold the orthonormalization coefficients.
  RCP<Teuchos::SerialDenseMatrix<int,ST> > z
    = rcp( new Teuchos::SerialDenseMatrix<int,ST>( numPorts_, numPorts_ ) );

  // Orthonormalize the Krylov kernel vectors V
  RCP<MV> V = rcp( new Epetra_MultiVector( RPtr_->epetraObj() ) );
  orthoMgr.normalize( *V, z );

  // Set the new state and initialize the solver to use R as the kernel vectors.
  Belos::GmresIterationState<ST,MV> initState;
  initState.V = V;
  initState.z = z;
  initState.curDim = 0;
  krylovIter.initializeGmres(initState);

  // Have the solver iterate until the basis size is computed
  try {
    krylovIter.iterate();
  }
  catch (const Belos::GmresIterationOrthoFailure &e) {
    // This might happen if the basis size is the same as the operator's dimension.
  }
  catch (const std::exception &e) {
    bsuccess = false; // Anything else is a failure.
  }

  // Get the basis vectors back from the iteration object
  Belos::GmresIterationState<ST,MV> newState = krylovIter.getState();

  // Create a view into the first k vectors of newState.V
  //  std::vector<int> indices( k );

  // search for a size for automated MOR

  RCP<MV> W = rcp( new Epetra_MultiVector( V->Map(), k ) );

  int low = 1, high = k;
  if (morAutoSize_)
  {

    double  abstol = 1e-6;
    double reltol = 1e-3;

    int mid;

    while ((high - low) > 0 )  
    {
      mid = (low + high)/2;

      std::vector<int> indices(mid);

 // Resize the projection matrices
      redG_.shape(mid, mid);
      redC_.shape(mid, mid);
      redB_.shape(mid, numPorts_);
      redL_.shape(mid, numPorts_);

      for (int i=0; i<mid; ++i) { indices[i] = i; }

      V = rcp( new Epetra_MultiVector( View, *newState.V, &indices[0], mid) );

      W = V;

// project the system
      Linear::MultiVector xyceV( &*V, false );
      Teuchos::RCP<Linear::MultiVector> temp2 =
        Teuchos::rcp( Xyce::Linear::createMultiVector( BaseMap, mid ) );
      Linear::MultiVector xyceW( &*W, false );

    // G * V
      GPtr_->matvec( false, xyceV, *temp2 );
    // V' * G * V

      MVT::MvTransMv( 1.0, xyceW.epetraObj(), temp2->epetraObj(), redG_ );

    // C * V
      CPtr_->matvec( false, xyceV, *temp2 );
    // V' * C * V
      MVT::MvTransMv( 1.0, xyceW.epetraObj(), temp2->epetraObj(), redC_ );
    //redC_.print( Xyce::dout() );

    // V' * B
      MVT::MvTransMv( 1.0, xyceW.epetraObj(), BPtr_->epetraObj(), redB_ );

      MVT::MvTransMv( 1.0, xyceV.epetraObj(), BPtr_->epetraObj(), redL_ );

      //redB_.print( Xyce::dout() );

      // calculate transfer functions:
      isSingleFreq_ = true;

      currentFreq_  = morMaxFreq_;
      origH_.shape(numPorts_, numPorts_);
      redH_.shape(numPorts_, numPorts_);

      evalOrigTransferFunction();
      evalRedTransferFunction();


      Teuchos::SerialDenseMatrix<int, double>  H_diff, totaltol, errOverTol;

      totaltol.shape(numPorts_, numPorts_);
      H_diff.shape(numPorts_, numPorts_);

 
      errOverTol.shape(numPorts_, numPorts_);
      for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
      {
        for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
        {
          totaltol(i,j) =  reltol*abs(  origH_(i,j)) +  abstol;

          H_diff(i,j) =  abs( origH_(i,j) - redH_(i, j));
      
          errOverTol(i, j) = H_diff(i,j)/totaltol(i,j);
        }
      }

      double totalErrOverTol =   errOverTol.normFrobenius()/numPorts_;


      if( totalErrOverTol < 1)  
        high = mid;
      else
        low = mid+1;

    }

    isSingleFreq_ = false;

    ROMsize_ =  high;

    if (ROMsize_ > RPtr_->globalLength())
    {
      kblock = (int)(RPtr_->globalLength() / numPorts_);
      UserWarning(*this) << "Requested reduced-order model dimension is larger than original system dimension, resizing to original system dimension";
    }
    else
    {
      kblock = (int)(ROMsize_ / numPorts_);
    }
  
    k = kblock * numPorts_;

  // Resize the projection matrices
    redG_.shape(k, k);
    redC_.shape(k, k);
    redB_.shape(k, numPorts_);
    redL_.shape(k, numPorts_);

  }
  std::vector<int> indices(k);
  for (int i=0; i<k; ++i) { indices[i] = i; }

  V = rcp( new Epetra_MultiVector( View, *newState.V, &indices[0], k) );
  W = V;

  // ---------------------------------------------------------------------
  // Sparsify the reduced system, if requested.
  // ---------------------------------------------------------------------
  if (morSparsificationType_)
  {

    Linear::MultiVector xyceV( &*V, false );
    Teuchos::RCP<Linear::MultiVector> temp2 = 
      Teuchos::rcp( Xyce::Linear::createMultiVector( BaseMap, k ) );

    // G * V
    GPtr_->matvec( false, xyceV, *temp2 );
    // V' * G * V
    MVT::MvTransMv( 1.0, xyceV.epetraObj(), temp2->epetraObj(), redG_ );
    //redG_.print( Xyce::dout() );

    // C * V
    CPtr_->matvec( false, xyceV, *temp2 );
    // V' * C * V
    MVT::MvTransMv( 1.0, xyceV.epetraObj(), temp2->epetraObj(), redC_ );
    //redC_.print( Xyce::dout() );

    // V' * B
    MVT::MvTransMv( 1.0, xyceV.epetraObj(), BPtr_->epetraObj(), redB_ );
    //redB_.print( Xyce::dout() );

    bsuccess = bsuccess & sparsifyRedSystem_();
  }

  // -----------------------------------------------------------------------------
  // Scale the basis vectors, if requested, before computing the projected system
  // -----------------------------------------------------------------------------

  // Create the scaling vector.
  Teuchos::SerialDenseMatrix<int,ST> Xhatscale( k, 1 );

  int scaleType = morScaleType_;

  if ( scaleType == 1 )
  {
    for (int i=0; i<k; ++i)
    {
      Xhatscale( i, 0 ) = morScaleFactor_;
    }
  }

//  if ( scaleType == 2 )
//  {

//  }

if ( scaleType == 2 || scaleType == 3  || scaleType ==4)
  {
    Epetra_MultiVector Xmag( V->Map(), 1 );

    if ( scaleType == 2 )
      Xmag.PutScalar( morScaleFactor_ );

    if ( scaleType == 4)
      Xmag.PutScalar(1.0);

    MVT::MvTransMv( 1.0, *V, Xmag, Xhatscale );

//    Xhatscale.print(Xyce::dout());
    for (int i=0; i<k; ++i)
    {

      if ( scaleType == 2 )
      {
        if ( fabs(Xhatscale( i, 0 )) >  morScaleFactor_)
        {
          Xhatscale( i, 0 ) = morScaleFactor_ / fabs( Xhatscale( i, 0 ) );
        }
        else
        {
          Xhatscale( i, 0 ) =  morScaleFactor1_ / fabs( Xhatscale( i, 0 ) );
        }
      }

      if ( scaleType == 3 )
        Xhatscale( i, 0 ) = morScaleFactor_ * fabs( Xhatscale( i, 0 ) );

      if ( scaleType == 4 )
      {
        if ( fabs(Xhatscale( i, 0 )) > 1.0 )
        {
          Xhatscale( i, 0 ) = morScaleFactor_ / fabs( Xhatscale( i, 0 ) );
        }
        else
        {
          Xhatscale( i, 0 ) =  morScaleFactor1_ / fabs( Xhatscale( i, 0 ) );
        }
      }


    }
  }

//  Xhatscale.print(Xyce::dout());

  // Scale the computed basis vectors before projecting the original system
  Teuchos::SerialDenseMatrix<int,ST> D(k,k);
  if ( scaleType != 0 )
  {
    RCP<MV> Vtemp = rcp( new Epetra_MultiVector( V->Map(), k ) );
//    Teuchos::SerialDenseMatrix<int,ST> D(k,k);
    for (int i=0; i<k; ++i)
    {
      if (Xhatscale( i, 0 ) != 0.0)
        D(i,i) = 1.0/Xhatscale( i, 0 );
      else
        D(i,i) = 1.0;
    }
    Xyce::dout() << " the scaling matrix "  <<  std::endl;

//    D.print(Xyce::dout());

    MVT::MvTimesMatAddMv( 1.0, *V, D, 0.0, *Vtemp );
    V = Vtemp;
  }

//  Xyce::dout() << "Printing out V" << std::endl;
//  V->Print(Xyce::dout());
//  Xyce::dout() << "Printing out W" << std::endl;
//  W->Print(Xyce::dout());

  // ---------------------------------------------------------------------
  // Now use the basis vectors for K_k(inv(G + s0*C)*C, R) to compute
  // the projected system.
  // ---------------------------------------------------------------------
  Linear::MultiVector xyceV( &*V, false );
  Linear::MultiVector xyceW( &*W, false );
  Teuchos::RCP<Linear::MultiVector> temp2 =
    Teuchos::rcp( Xyce::Linear::createMultiVector( BaseMap, k ) );

  if (!morSparsificationType_)
  {
  // G * V
    GPtr_->matvec( false, xyceV, *temp2 );
  // V' * G * V
    MVT::MvTransMv( 1.0, xyceW.epetraObj(), temp2->epetraObj(), redG_ );
    //redG_.print( Xyce::dout() );

  // C * V
    CPtr_->matvec( false, xyceV, *temp2 );
  // V' * C * V
    MVT::MvTransMv( 1.0, xyceW.epetraObj(), temp2->epetraObj(), redC_ );

  // V' * B
    MVT::MvTransMv( 1.0, xyceW.epetraObj(), BPtr_->epetraObj(), redB_ );

    MVT::MvTransMv( 1.0, xyceV.epetraObj(), BPtr_->epetraObj(), redL_ );
  }
  else
  {
    if ( scaleType <= 1 )
    {
      redCPtr_->scale(1.0/morScaleFactor_);
      redGPtr_->scale(1.0/morScaleFactor_);

      redL_.scale(1.0/morScaleFactor_);
    }
    else
    {
      Report::UserError0() << "MOR options sparsificationType=1 can only be used with scaletype=1, other scale types have not been supported for sparsification";
      return false;
    }
  }

  // ---------------------------------------------------------------------
  // Output the projected system, redG_, redC_, and redB_.
  // ---------------------------------------------------------------------

  if (morSaveRedSys_)
  {
    // If sparsification was used, write out the Linear::Matrix object instead of the dense object.
    if (morSparsificationType_)
      IO::outputROM(outputManagerAdapter_.getComm(), analysisManager_.getNetlistFilename(), *redGPtr_, *redCPtr_, redB_, redL_ );
    else
      IO::outputROM(outputManagerAdapter_.getComm(), analysisManager_.getNetlistFilename(), redG_, redC_, redB_, redL_ );  // L = B'
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::evalOrigTransferFunction()
// Purpose       : Evaluate original transfer function.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 5/29/2012
//-----------------------------------------------------------------------------
bool MOR::evalOrigTransferFunction()
{
  bool bsuccess = true;

  createOrigLinearSystem_();

  int currentStep = 0;
  int finalStep = morEvalSize_;

  if ( isSingleFreq_ )
    finalStep = 1; 

  bool stepAttemptStatus;

  while (currentStep < finalStep)
  {

    if (!isSingleFreq_)
      updateCurrentFreq_(currentStep);

    updateOrigLinearSystemFreq_();

    stepAttemptStatus = solveOrigLinearSystem_();

    currentStep++;

    if (stepAttemptStatus)
    {
      processSuccessfulStep(true);
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      doProcessFailedStep();
    }

  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::evalRedTransferFunction()
// Purpose       : Evaluate reduced transfer function.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 5/29/2012
//-----------------------------------------------------------------------------
bool MOR::evalRedTransferFunction()
{
  bool bsuccess = true;

  createRedLinearSystem_();

  int currentStep = 0;

  int finalStep = morEvalSize_;

  if ( isSingleFreq_ )
    finalStep = 1; 

  bool stepAttemptStatus;

  while (currentStep < finalStep)
  {
    if (!isSingleFreq_)
      updateCurrentFreq_(currentStep);

    updateRedLinearSystemFreq_();

    stepAttemptStatus = solveRedLinearSystem_();

    currentStep++;

    if (stepAttemptStatus)
    {
      processSuccessfulStep(false);
    }
    else // stepAttemptStatus  (ie do this if the step FAILED)
    {
      doProcessFailedStep();
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::createOrigLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool MOR::createOrigLinearSystem_()
{
  bool bsuccess = true;

  Parallel::Manager &pdsManager = *analysisManager_.getPDSManager();

  Parallel::ParMap &BaseMap = *pdsManager.getParallelMap(Parallel::SOLUTION);
  Parallel::ParMap &oBaseMap = *pdsManager.getParallelMap(Parallel::SOLUTION_OVERLAP_GND);
  const Linear::Graph* baseFullGraph = pdsManager.getMatrixGraph(Parallel::JACOBIAN);

  int numBlocks = 2;

  std::vector<RCP<Parallel::ParMap> > blockMaps = Linear::createBlockParMaps(numBlocks, BaseMap, oBaseMap);

  std::vector<std::vector<int> > blockPattern(2);
  blockPattern[0].resize(2);
  blockPattern[0][0] = 0; blockPattern[0][1] = 1;
  blockPattern[1].resize(2);
  blockPattern[1][0] = 0; blockPattern[1][1] = 1;

  int offset = Linear::generateOffset( BaseMap );
  RCP<Linear::Graph> blockGraph = Linear::createBlockGraph( offset, blockPattern, *blockMaps[0], *baseFullGraph);

  sCpG_REFMatrixPtr_ = rcp ( Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), baseFullGraph) );

  // Load diagonal blocks of real equivalent form: (G - s0*C)
  sCpG_REFMatrixPtr_->put( 0.0 ); // Zero out whole matrix

  // Add scaled C matrix first if expansion point is not zero.
  if (s0_ != 0.0)
  {
    sCpG_REFMatrixPtr_->block( 0, 0 ).add(*CPtr_);
    sCpG_REFMatrixPtr_->block( 0, 0 ).scale(-s0_);
    sCpG_REFMatrixPtr_->block( 1, 1 ).add(*CPtr_);
    sCpG_REFMatrixPtr_->block( 1, 1 ).scale(-s0_);
  }
  // Add G into diagonal block to obtain (G - s0*C)
  sCpG_REFMatrixPtr_->block( 0, 0 ).add(*GPtr_);
  sCpG_REFMatrixPtr_->block( 1, 1 ).add(*GPtr_);

  // Create solver factory
  Amesos amesosFactory;

  // Create block vectors
  REFBPtr_ = rcp ( Xyce::Linear::createBlockVector ( numBlocks, blockMaps[0], rcp(&BaseMap, false) ) );
  REFXPtr_ = rcp ( Xyce::Linear::createBlockVector ( numBlocks, blockMaps[0], rcp(&BaseMap, false) ) );
  REFXPtr_->putScalar( 0.0 );

  Teuchos::RCP<Linear::EpetraBlockMatrix> e_sCpG = Teuchos::rcp_dynamic_cast<Linear::EpetraBlockMatrix>(sCpG_REFMatrixPtr_);
  blockProblem_ = rcp(new Epetra_LinearProblem(&e_sCpG->epetraObj(), &REFXPtr_->epetraObj(), &REFBPtr_->epetraObj() ) );

  blockSolver_ = rcp( amesosFactory.Create( "Klu", *blockProblem_ ) );

  int linearStatus = blockSolver_->SymbolicFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : MOR::createRedLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------
bool MOR::createRedLinearSystem_()
{
  // Create a dense linear system if sparsification has not been, or will not be, performed.
  if ( !morSparsificationType_ || !isROMSparse_ )
  {
    // The reduced system is dense, so create a dense linear system for transfer function comparisons.

    // Get the reduced system size.
    int k = redG_.numRows();

    // Resize serial dense matrix for real-equivalent form.
    sCpG_redMatrix_.shape(2*k, 2*k);
    sCpG_tmpMatrix_.shape(2*k, 2*k);
    ref_redB_.shape(2*k, numPorts_);

    // First, load the diagonals.
    // Get a view of the primary diagonal block, insert (G - s0*C).
    Teuchos::SerialDenseMatrix<int, double> subMtx( Teuchos::View, sCpG_tmpMatrix_, k, k, 0, 0 );
    subMtx.assign( redC_ );
    subMtx.scale( -s0_ );
    subMtx += redG_;

    // Get a view of the lower k x k diagonal block, insert (G - s0*C).
    Teuchos::SerialDenseMatrix<int, double> subMtx2( Teuchos::View, sCpG_tmpMatrix_, k, k, k, k );
    subMtx2.assign( redC_ );
    subMtx2.scale( -s0_ );
    subMtx2 += redG_;
  }
  else
  {
    // The reduced system is sparse, so create a sparse linear system.

    int numBlocks = 2;

    RCP<Parallel::ParMap> redBlockMapPtr = Linear::createBlockParMap(numBlocks, *redMapPtr_);
 
    std::vector<std::vector<int> > blockPattern(2);
    blockPattern[0].resize(2);
    blockPattern[0][0] = 0; blockPattern[0][1] = 1;
    blockPattern[1].resize(2);
    blockPattern[1][0] = 0; blockPattern[1][1] = 1;
   
    int offset= Linear::generateOffset( *redMapPtr_ );
    RCP<Linear::Graph> blockGraph = Linear::createBlockGraph( offset, blockPattern, *redBlockMapPtr, *(redCPtr_->getGraph()) );

    sCpG_ref_redMatrixPtr_ = rcp( Linear::createBlockMatrix( numBlocks, offset, blockPattern, blockGraph.get(), redCPtr_->getGraph() ) );

    // Load diagonal blocks of real equivalent form: (G - s0*C)
    sCpG_ref_redMatrixPtr_->put( 0.0 ); // Zero out whole matrix

    // Add scaled C matrix first if expansion point is not zero.
    if (s0_ != 0.0)
    {
      sCpG_ref_redMatrixPtr_->block( 0, 0 ).add(*redCPtr_);
      sCpG_ref_redMatrixPtr_->block( 0, 0 ).scale(-s0_);
      sCpG_ref_redMatrixPtr_->block( 1, 1 ).add(*redCPtr_);
      sCpG_ref_redMatrixPtr_->block( 1, 1 ).scale(-s0_);
    }
    // Add G into diagonal block to obtain (G - s0*C)
    sCpG_ref_redMatrixPtr_->block( 0, 0 ).add(*redGPtr_);
    sCpG_ref_redMatrixPtr_->block( 1, 1 ).add(*redGPtr_);

    // Create solver factory
    Amesos amesosFactory;

    // Create a block vector
    ref_redBPtr_ = rcp ( Xyce::Linear::createBlockVector ( numBlocks, redBlockMapPtr, redMapPtr_ ) );
    ref_redXPtr_ = rcp ( Xyce::Linear::createBlockVector ( numBlocks, redBlockMapPtr, redMapPtr_ ) );
    ref_redXPtr_->putScalar( 0.0 );

    Teuchos::RCP<Linear::EpetraBlockMatrix> e_sCpG = Teuchos::rcp_dynamic_cast<Linear::EpetraBlockMatrix>(sCpG_ref_redMatrixPtr_);
    blockRedProblem_ = rcp(new Epetra_LinearProblem(&e_sCpG->epetraObj(), 
                                                    &ref_redXPtr_->epetraObj(), &ref_redBPtr_->epetraObj() ) );

    blockRedSolver_ = rcp( amesosFactory.Create( "Klu", *blockRedProblem_ ) );

    int linearStatus = blockRedSolver_->SymbolicFactorization();

    if (linearStatus != 0)
    {
      Xyce::dout() << "Amesos symbolic factorization exited with error: " << linearStatus;
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::updateOrigLinearSystemFreq_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------

bool MOR::updateOrigLinearSystemFreq_()
{
  double omega =  2.0 * M_PI * currentFreq_;

  sCpG_REFMatrixPtr_->block( 0, 1).put( 0.0);
  sCpG_REFMatrixPtr_->block( 0, 1).add(*CPtr_);
  sCpG_REFMatrixPtr_->block( 0, 1).scale(-omega);

  sCpG_REFMatrixPtr_->block(1, 0).put( 0.0);
  sCpG_REFMatrixPtr_->block(1, 0).add(*CPtr_);
  sCpG_REFMatrixPtr_->block(1, 0).scale(omega);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::updateRedLinearSystemFreq_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/20/2011
//-----------------------------------------------------------------------------

bool MOR::updateRedLinearSystemFreq_()
{
  double omega =  2.0 * M_PI * currentFreq_;

  if ( morSparsificationType_ && isROMSparse_ )
  {
    sCpG_ref_redMatrixPtr_->block( 0, 1).put( 0.0);
    sCpG_ref_redMatrixPtr_->block( 0, 1).add(*redCPtr_);
    sCpG_ref_redMatrixPtr_->block( 0, 1).scale(-omega);

    sCpG_ref_redMatrixPtr_->block(1, 0).put( 0.0);
    sCpG_ref_redMatrixPtr_->block(1, 0).add(*redCPtr_);
    sCpG_ref_redMatrixPtr_->block(1, 0).scale(omega);
  }
  else 
  {
    int k = redG_.numRows();

    // Reset reduced real equivalent form matrix to sCpG_tmpMatrix_
    sCpG_redMatrix_.assign( sCpG_tmpMatrix_ );

    // Insert off diagonals.
    Teuchos::SerialDenseMatrix<int, double> subMtx( Teuchos::View, sCpG_redMatrix_, k, k, 0, k );
    subMtx.assign( redC_ );
    subMtx.scale( -omega );

    Teuchos::SerialDenseMatrix<int, double> subMtx2( Teuchos::View, sCpG_redMatrix_, k, k, k, 0 );
    subMtx2.assign( redC_ );
    subMtx2.scale( omega );
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::solveLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------

bool MOR::solveOrigLinearSystem_()
{

  bool bsuccess = true;

  // Solve the block problem
  int linearStatus = blockSolver_->NumericFactorization();

  if (linearStatus != 0)
  {
    Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
    bsuccess = false;
  }

  // Loop over number of I/O ports here
  for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
  {
    REFBPtr_->putScalar( 0.0 );
    (REFBPtr_->block( 0 ))[bMatEntriesVec_[j]] = -1.0;

    linearStatus = blockSolver_->Solve();
    if (linearStatus != 0)
    {
      Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
      bsuccess = false;
    }

    // Compute transfer function entries for all I/O
    for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
    {
       // Populate H for all ports in L
       // L is the same as B and also a set of canonical basis vectors (e_i), so
       // we can pick off the appropriate entries of REFXPtr to place into H.
       origH_(i,j) = std::complex<double>(-(REFXPtr_->block( 0 ))[bMatEntriesVec_[i]],
                                          -(REFXPtr_->block( 1 ))[bMatEntriesVec_[i]]);
    }
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::solveRedLinearSystem_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :  Ting Mei, Heidi Thornquist, SNL
// Creation Date : 6/2011
//-----------------------------------------------------------------------------

bool MOR::solveRedLinearSystem_()
{
  bool bsuccess = true;

  if ( morSparsificationType_ && isROMSparse_ )
  {
    // Solve the block problem
    int linearStatus = blockRedSolver_->NumericFactorization();

    if (linearStatus != 0)
    {
      Xyce::dout() << "Amesos numeric factorization exited with error: " << linearStatus;
      bsuccess = false;
    }

    // Create a multivector with the reduced B and L vectors.
    Teuchos::RCP<Parallel::EpetraParMap> e_redMapPtr_ = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>( redMapPtr_ );
    Epetra_MultiVector tmp_redL( View, *e_redMapPtr_->petraMap(), redL_.values(), redL_.stride(), numPorts_ ); 
    Epetra_MultiVector tmp_redB( View, *e_redMapPtr_->petraMap(), redB_.values(), redB_.stride(), numPorts_ ); 
    Linear::MultiVector redBmv( &tmp_redB, tmp_redB.Map(), false ), redLmv( &tmp_redL, tmp_redL.Map(), false );

    // Loop over number of I/O ports here
    for (int j=0; j < numPorts_; ++j)
    {
      // Set the B vector for this port.
      ref_redBPtr_->putScalar( 0.0 );
      (ref_redBPtr_->block( 0 )) = Linear::Vector( (tmp_redB)(j), tmp_redB.Map(), false );

      linearStatus = blockRedSolver_->Solve();
      if (linearStatus != 0)
      {
        Xyce::dout() << "Amesos solve exited with error: " << linearStatus;
        bsuccess = false;
      }

      // Compute transfer function entries for all I/O
      for (int i=0; i < numPorts_; ++i)
      {
         RCP<const Linear::Vector> redLmv_i = rcp( redLmv.getVectorView( i ) );

         // L is probably not B in this case, so explicitly multiply to compute transfer function entries.
         double realPart = redLmv_i->dotProduct( ref_redXPtr_->block( 0 ) );
         double imagPart = redLmv_i->dotProduct( ref_redXPtr_->block( 1 ) );

         redH_(i,j) = std::complex<double>(realPart, imagPart);
      }
    }
  }
  else 
  {
    int k = redG_.numRows();
    Teuchos::LAPACK<int, double> lapack;

    ref_redB_.putScalar( 0.0 );
    Teuchos::SerialDenseMatrix<int, double> tmpRedB( Teuchos::View, ref_redB_, redB_.numRows(), redB_.numCols(), 0, 0 );
    tmpRedB.assign( redB_ );

    // First factor matrix using LU.
    int info = 0;
    std::vector<int> ipiv( sCpG_redMatrix_.numRows() );
    lapack.GETRF( sCpG_redMatrix_.numRows(), sCpG_redMatrix_.numCols(), sCpG_redMatrix_.values(),
                sCpG_redMatrix_.stride(), &ipiv[0], &info );
    if (info != 0)
    {
      Xyce::dout() << "LAPACK::GETRF: LU factorization failed with error: " << info << std::endl;
      bsuccess = false;
    }

    // Next solve for all ports at once using LU factors.
    const char trans = 'N';
    lapack.GETRS( trans, sCpG_redMatrix_.numRows(), ref_redB_.numCols(), sCpG_redMatrix_.values(),
                  sCpG_redMatrix_.stride(), &ipiv[0], ref_redB_.values(), ref_redB_.stride(), &info );
    if (info != 0)
    {
      Xyce::dout() << "LAPACK::GETRS: LU solve failed with error: " << info << std::endl;
      bsuccess = false;
    }  

    Teuchos::SerialDenseMatrix<int, double> tmpRedImag( numPorts_, numPorts_ ), tmpRedReal( numPorts_, numPorts_ );

    //  Compute real part.
    tmpRedReal.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, redB_, tmpRedB, 0.0 );

    // Compute imaginary part.
    Teuchos::SerialDenseMatrix<int, double> tmpRedB2( Teuchos::View, ref_redB_, redB_.numRows(), redB_.numCols(), k, 0 );
    tmpRedImag.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, redB_, tmpRedB2, 0.0 );

    // Compute transfer function entries for all I/O
    for (unsigned int j=0; j < bMatEntriesVec_.size(); ++j)
    {
      for (unsigned int i=0; i < bMatEntriesVec_.size(); ++i)
      {
        redH_(i,j) = std::complex<double>(tmpRedReal(i,j), tmpRedImag(i,j));
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::sparsifyRedSystem_()
// Purpose       : Use techniques to sparsify the projected system.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 8/20/2012
//-----------------------------------------------------------------------------
bool MOR::sparsifyRedSystem_()
{
  bool bsuccess = true;

  int k = redB_.numRows();

  // Storage vectors
  Teuchos::SerialDenseMatrix<int, double> R(redB_);
  Teuchos::SerialDenseMatrix<int, double> A(redC_);

  // Factor Ghat.
  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> tmp_Ghat( redG_ );

  int info = 0;
  std::vector<int> ipiv( tmp_Ghat.numRows() );
  lapack.GETRF( tmp_Ghat.numRows(), tmp_Ghat.numCols(), tmp_Ghat.values(),
                tmp_Ghat.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization of Ghat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute reciprocal condition estimate of Ghat.
  const char norm = '1';
  double condGhat = 0.0;
  std::vector<double> condWork( 4*k );
  std::vector<int> condIWork( k );
  lapack.GECON( norm, tmp_Ghat.numRows(), tmp_Ghat.values(), tmp_Ghat.stride(), tmp_Ghat.normOne(), &condGhat, &condWork[0], &condIWork[0], &info );
  // Xyce::dout() << "Condition estimate for Ghat is : " << 1.0/ condGhat << std::endl;
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GECON: Computing condition estimate of Ghat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute A = inv(Ghat)* Chat
  const char trans = 'N';
  lapack.GETRS( trans, tmp_Ghat.numRows(), A.numCols(), tmp_Ghat.values(),
                tmp_Ghat.stride(), &ipiv[0], A.values(), A.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve failed with error: " << info << std::endl;
    bsuccess = false;
  }
  //Xyce::dout() << "A = inv(Ghat)*Chat" << std::endl;
  //A.print(Xyce::dout());

  // Compute R = inv(Ghat)* Bhat
  lapack.GETRS( trans, tmp_Ghat.numRows(), R.numCols(), tmp_Ghat.values(),
                tmp_Ghat.stride(), &ipiv[0], R.values(), R.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve failed with error: " << info << std::endl;
    bsuccess = false;
  }
  //Xyce::dout() << "R = inv(Ghat)*Bhat" << std::endl;
  //R.print(Xyce::dout());

  // Reduce A to Schur form and compute the eigenvalues.
  // Allocate the work space. This space will be used below for calls to:
  // * GEES  (requires 3*k for real)
  // * TREVC (requires 3*k for real)
  // Furthermore, GEES requires a real array of length k
  //
  int lwork = 8*k;
  std::vector<double> work( lwork );
  std::vector<double> revals( k );
  std::vector<double> ievals( k );
  Teuchos::SerialDenseMatrix<int, double> Q(k, k);

  // Create diagonal tmpB for the eigenproblem (A, tmpB), beta should return as a vector of ones.
  const int ldvl = 1;
  double vl[ ldvl ];
  std::vector<double> beta( k );
  Teuchos::SerialDenseMatrix<int, double> tmpB( k, k );
  for (int i=0; i<k; i++)
    tmpB(i,i) = 1.0;
  lapack.GGEV( 'N', 'V', k, A.values(), A.stride(), tmpB.values(), tmpB.stride(),
               &revals[0], &ievals[0], &beta[0], &vl[0], ldvl, Q.values(), Q.stride(),
               &work[0], lwork, &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GGEV: Computing eigenvectors and values of A = inv(Ghat)*Chat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Scale the eigenvectors returned from GGEV to have 2-norm = 1
  // They are initially scaled to have inf-norm = 1 from LAPACK
  Teuchos::BLAS<int, double> blas;
  std::vector<int> rowIndex( k );
  int i = 0;

  while( i < k ) {
    if ( ievals[i] != 0.0 )
    {
      rowIndex[i] = 1;
      rowIndex[i+1] = -1;

      // Compute the 2-norm of the complex eigenvector
      double norm_r = blas.NRM2( k, Q[i], 1 );
      double norm_i = blas.NRM2( k, Q[i+1], 1 );
      double norm = lapack.LAPY2( norm_r, norm_i );

      // Scale the complex eigenvector
      for (int j=0; j<k; j++)
      {
        Q(j,i)   /= norm;
        Q(j,i+1) /= norm;
      }
      i += 2;
    }
    else
    {
      rowIndex[i] = 0;

      // Compute the 2-norm of the i-th column
      double norm = blas.NRM2( k, Q[i], 1 );

      // Scale the real eigenvector
      for (int j=0; j<k; j++)
      {
        Q(j,i) /= norm;
      }
      i++;
    }
  }

  // Xyce::dout() << "Eigenvalues of A: " << std::endl;
  // for (int i=0; i<k; ++i)
  //   Xyce::dout() << revals[i] << "\t" << ievals[i] << "\t" << beta[i] << std::endl;

  // Xyce::dout() << "Eigenvectors of A: " << std::endl;
  // Q.print(Xyce::dout());

  // Construct complex eigenvectors
  Teuchos::SerialDenseMatrix<int, std::complex<double> > V(k, k);
  for (int i=0; i<k; i++)
  {
    if (rowIndex[i] == 1)
    {
      for (int j=0; j<k; j++)
      {
        // Insert both conjugate pairs simultaneously.
        V(j,i) = std::complex<double>( Q(j,i), Q(j,i+1) );
        V(j,i+1) = std::complex<double>( Q(j,i), -Q(j,i+1) );
      }
      i++;  // Need to increment an extra value for the conjugate pair.
    }
    else  // The eigenvector is real, copy directly
    {
      for (int j=0; j<k; j++)
      {
        V(j,i) = std::complex<double>( Q(j,i), 0.0 );
      }
    }
  }

  // Factor V
  Teuchos::LAPACK<int, std::complex<double> > lapack_complex;
  lapack_complex.GETRF( V.numRows(), V.numCols(), V.values(), V.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization of eigenvectors failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute condition estimate of V.
  double condV = 0.0;
  std::vector<std::complex<double> > condCWork( 2*k );
  lapack_complex.GECON( norm, V.numRows(), V.values(), V.stride(), V.normOne(), &condV, &condCWork[0], &condWork[0], &info );
  //Xyce::dout() << "Condition estimate for V is : " << 1.0/condV << std::endl;
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GECON: Computing condition estimate of V failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute inv(V)
  std::vector<std::complex<double> > cwork( lwork );
  lapack_complex.GETRI( V.numRows(), V.values(), V.stride(), &ipiv[0], &cwork[0], lwork, &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRI: Computing inverse of V failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Convert R to a complex vector to use the multiply routine for inv(V).
  Teuchos::SerialDenseMatrix<int, std::complex<double> > tmpR( k, numPorts_ );
  for (int j=0; j<numPorts_; j++)
    for (int i=0; i<k; i++)
      tmpR(i,j) = std::complex<double>( R(i,j), 0.0 );

  // Compute inv(V) * R
  Teuchos::SerialDenseMatrix<int, std::complex<double> > tmp_redB(k ,numPorts_);
  tmp_redB.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, V, tmpR, 0.0);

  // Generate compressed real storage of V*R and inv(V)
  redL_ = redB_;  // Store copy of redB_ before destroying it below.
  Teuchos::SerialDenseMatrix<int, double> invV(k, k);
  for (int i=0; i<k; i++)
  {
    if (rowIndex[i] == 1)  // Complex conjugate pair.
    {
      for (int j=0; j<k; j++)
      {
        invV(i, j) = V(i,j).real();
        invV(i+1, j) = V(i,j).imag();
      }
      for (int j=0; j<numPorts_; j++)
      {
        redB_(i, j) = tmp_redB(i, j).real();
        redB_(i+1, j) = tmp_redB(i, j).imag();
      }
      i++;  // Increment by one to skip complex conjugate.
    }
    else  // Eigenvalue is real, so is eigenvector.
    {
      for (int j=0; j<k; j++)
        invV(i, j) = V(i,j).real();
      for (int j=0; j<numPorts_; j++)
        redB_(i, j) = tmp_redB(i, j).real();
    }
  }

  // Factor invV.
  lapack.GETRF( invV.numRows(), invV.numCols(), invV.values(), invV.stride(), &ipiv[0], &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRF: LU factorization of inv(V) failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Compute A = inv(invV)^T*Lhat
  const char trans2 = 'T';
  lapack.GETRS( trans2, invV.numRows(), redL_.numCols(), invV.values(),
                invV.stride(), &ipiv[0], redL_.values(), redL_.stride(), &info );
  if (info != 0)
  {
    Xyce::dout() << "LAPACK::GETRS: LU solve for Lhat failed with error: " << info << std::endl;
    bsuccess = false;
  }

  // Xyce::dout() << "Bhat : " << std::endl;
  // redB_.print(Xyce::dout());
  // Xyce::dout() << "Lhat : " << std::endl;
  // redL_.print(Xyce::dout());

  // Create redCPtr_ and redGPtr_
  // Generate Epetra_Map that puts all the values of the reduced system on one processor.
  Parallel::Manager &pdsManager = *analysisManager_.getPDSManager();
  Parallel::ParMap &BaseMap = *pdsManager.getParallelMap(Parallel::SOLUTION);
  Parallel::EpetraParMap &e_BaseMap = dynamic_cast<Parallel::EpetraParMap&>(BaseMap); 
  if (BaseMap.pdsComm().procID() == 0)
    redMapPtr_ = Teuchos::rcp( Parallel::createPDSParMap( k, k, 0, BaseMap.pdsComm() ) );
  else
    redMapPtr_ = Teuchos::rcp( Parallel::createPDSParMap( k, 0, 0, BaseMap.pdsComm() ) );

  Teuchos::RCP<Parallel::EpetraParMap> e_redMapPtr_ = Teuchos::rcp_dynamic_cast<Parallel::EpetraParMap>( redMapPtr_ );
  Epetra_CrsMatrix* tmpRedG = new Epetra_CrsMatrix( Copy, *e_redMapPtr_->petraMap(), 1 );
  Epetra_CrsMatrix* tmpRedC = new Epetra_CrsMatrix( Copy, *e_redMapPtr_->petraMap(), 2 );

  // Let processor 0 insert entries into tmpRedG and tmpRedC
  if (BaseMap.pdsComm().procID() == 0)
  {
    std::vector<int> index(2);
    std::vector<double> val(2);
    for (int i=0; i<k; i++)
    {
      // Insert G, which is a diagonal of ones.
      index[0] = i; val[0] = 1.0;
      tmpRedG->InsertGlobalValues( i, 1, &val[0], &index[0] );

      // Insert C, which is block diagonal, where the blocks represent conjugate eigenvalues.
      if (rowIndex[i] == 1)
      {
        index[0] = i; index[1] = i+1;
        val[0] = revals[i]; val[1] = -ievals[i];
        tmpRedC->InsertGlobalValues( i, 2, &val[0], &index[0] );
      }
      else if (rowIndex[i] == -1)
      {
        index[0] = i-1; index[1] = i;
        val[0] = ievals[i-1]; val[1] = revals[i-1];
        tmpRedC->InsertGlobalValues( i, 2, &val[0], &index[0] );
      }
      else // Must be real.
      {
        index[0] = i; val[0] = revals[i];
        tmpRedC->InsertGlobalValues( i, 1, &val[0], &index[0] );
      }
    }
  }
  // We are done inserting values.
  tmpRedC->FillComplete();
  tmpRedC->OptimizeStorage();
  tmpRedG->FillComplete();
  tmpRedG->OptimizeStorage();

  // Pass tmpRedG and tmpRedC into redGPtr_ and redGPtr_, respectively, which will manage the memory.
  redGPtr_ = rcp( new Linear::EpetraMatrix( tmpRedG ) );
  redCPtr_ = rcp( new Linear::EpetraMatrix( tmpRedC ) );

  // We now have a sparse ROM.
  if (bsuccess)
    isROMSparse_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::processSuccessfulStep()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::processSuccessfulStep(bool origSys)
{
  bool bsuccess = true;

  // Output H.
  if (!isSingleFreq_)
  {
    if (origSys)
    {
      // outputManagerAdapter_.outputMORTF (origSys, currentFreq_, origH_);
      outputMOR_.output(outputManagerAdapter_.getComm(), origSys, currentFreq_, origH_ );
    }
    else
    {
      // outputManagerAdapter_.outputMORTF (origSys, currentFreq_, redH_);
      outputMOR_.output(outputManagerAdapter_.getComm(), origSys, currentFreq_, redH_ );
    }
  }

  if ( !firstDoubleDCOPStep() )
  {
      stepNumber += 1;
      stats_.successStepsThisParameter_ += 1;
      stats_.successfulStepsTaken_ += 1;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::processFailedStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::doProcessFailedStep()
{
  bool bsuccess = true;

  stepNumber += 1;
  morEvalFailures_.push_back(stepNumber);
  stats_.failedStepsAttempted_  += 1;
  analysisManager_.getStepErrorControl().numberSuccessiveFailures += 1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : MOR::finish
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::doFinish()
{
  bool bsuccess = true;

  if (DEBUG_ANALYSIS)
    Xyce::dout() << "Calling MOR::doFinish() outputs!" << std::endl;

  if (!(morEvalFailures_.empty()))
  {
    bsuccess = false;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : MOR::doHandlePredictor
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 06/24/2013
//-----------------------------------------------------------------------------
bool MOR::doHandlePredictor()
{
  analysisManager_.getDataStore()->setErrorWtVector(tiaParams_, topology_.getVarTypes());
  analysisManager_.getWorkingIntegrationMethod().obtainPredictor();
  analysisManager_.getWorkingIntegrationMethod().obtainPredictorDeriv();

  // In case this is the upper level of a 2-level sim, tell the
  // inner solve to do its prediction:
  bool        beginIntegrationFlag = analysisManager_.getBeginningIntegrationFlag();           // system_state.beginIntegrationFlag;
  double      nextTimeStep = analysisManager_.getStepErrorControl().currentTimeStep;           // system_state.nextTimeStep;
  double      nextTime = analysisManager_.getStepErrorControl().nextTime;                      // system_state.nextTime;
  int         currentOrder = analysisManager_.getWorkingIntegrationMethod().getOrder();        // system_state.currentOrder;

  loader_.startTimeStep(beginIntegrationFlag, nextTimeStep, nextTime, currentOrder);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::updateCurrentFreq_(int stepNumber )
// Purpose       :
// Special Notes : Used for MOR analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
bool MOR::updateCurrentFreq_(int stepNumber)
{
  double fStart = morCompFStart_;

  if (equal_nocase(morCompType_, "LIN"))
  {
    currentFreq_  = fStart + static_cast<double>(stepNumber)*fStep_;
  }
  else if (equal_nocase(morCompType_, "DEC") || equal_nocase(morCompType_, "OCT"))
  {
    currentFreq_ = fStart*pow(stepMult_, static_cast<double>(stepNumber) );
  }
  else
  {
    Report::DevelFatal().in("MOR::updateCurrentFreq_") << "Unsupported STEP type";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : MOR::setupSweepParam_
// Purpose       : Processes sweep parameters.
// Special Notes : Used for MOR analysis classes.
// Scope         : public
// Creator       : Ting Mei, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
int MOR::setupSweepParam_()
{
  int fCount = 0;
  double fStart = morCompFStart_;
  double fStop = morCompFStop_;

  if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
  {
    Xyce::dout() << std::endl << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "MOR::setupSweepParam_" << std::endl;
  }

  if (equal_nocase(morCompType_, "LIN"))
  {
    if ( morCompNP_ == 1)
      fStep_ = 0;
    else
      fStep_  = (fStop - fStart)/(morCompNP_ - 1);

    fCount = morCompNP_;

    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      Xyce::dout() << "fStep   = " << fStep_  << std::endl;
    }
  }
  else if (equal_nocase(morCompType_, "DEC"))
  {
    stepMult_ = pow(10,(1/(double)morCompNP_));
    fCount   = (int)(floor(fabs(log10(fStart) - log10(fStop)) * morCompNP_ + 1));
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
    }
  }
  else if (equal_nocase(morCompType_, "OCT"))
  {
    stepMult_ = pow(2,1/(double)(morCompNP_));

    // changed to remove dependence on "log2" function, which apparently
    // doesn't exist in the math libraries of FreeBSD or the mingw
    // cross-compilation suite.   Log_2(x)=log_e(x)/log_e(2.0)
    double ln2=log(2.0);
    fCount   = floor(fabs(log(fStart) - log(fStop))/ln2 * morCompNP_ + 1);
    if (DEBUG_ANALYSIS && isActive(Diag::TIME_PARAMETERS))
    {
      Xyce::dout() << "stepMult_ = " << stepMult_  << std::endl;
    }
  }
  else
  {
    Report::DevelFatal().in("MOR::setupSweepParam_") << "Unsupported type";
  }

  // At this point, pinterval equals the total number of steps
  // for the step loop.
  return fCount;
}

namespace {

typedef Util::Factory<AnalysisBase, MOR>  MORFactoryBase;

//-----------------------------------------------------------------------------
// Class         : MORFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Thu Jan 29 12:53:02 2015
//-----------------------------------------------------------------------------
///
/// Factory for parsing MOR parameters from the netlist and creating MOR analysis.
///
class MORFactory : public MORFactoryBase
{
public:
  //-----------------------------------------------------------------------------
  // Function      : MORFactory
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:54:09 2015
  //-----------------------------------------------------------------------------
  ///
  /// Constructs the MOR analysis factory
  ///
  /// @invariant Stores the results of parsing, so if more than one of the analysis and
  /// associated lines are parsed, the second options simply overwrite the previously parsed
  /// values.
  ///
  /// @invariant The existence of the parameters specified in the constructor cannot
  /// change.
  ///
  /// @param analysis_manager 
  /// @param linear_system 
  /// @param nonlinear_manager 
  /// @param topology 
  ///
  MORFactory(
    Analysis::AnalysisManager &         analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager)
    : MORFactoryBase(),
      analysisManager_(analysis_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager)
  {}

  virtual ~MORFactory()
  {}

  //-----------------------------------------------------------------------------
  // Function      : create
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 12:59:00 2015
  //-----------------------------------------------------------------------------
  ///
  /// Create a new MOR analysis and applies the analysis and time integrator option blocks.
  ///
  /// @return new MOR analysis object
  ///
  MOR *create() const
  {
    analysisManager_.setAnalysisMode(ANP_MODE_MOR);

    MOR *mor = new MOR(analysisManager_, linearSystem_, nonlinearManager_, loader_, topology_, initialConditionsManager_);
    mor->setAnalysisParams(morAnalysisOptionBlock_);
    mor->setMOROptions(morOptsOptionBlock_);

    return mor;
  }

  //-----------------------------------------------------------------------------
  // Function      : setMORAnalysisOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the analysis parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified analysis option block.
  ///
  /// @param option_block parsed option block
  ///
  void setMORAnalysisOptionBlock(const Util::OptionBlock &option_block)
  {
    morAnalysisOptionBlock_ = option_block;
  }

  //-----------------------------------------------------------------------------
  // Function      : setMOROptsOptionBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Thu Jan 29 13:00:14 2015
  //-----------------------------------------------------------------------------
  ///
  /// Saves the MOR_OPTS parsed options block in the factory.
  ///
  /// @invariant Overwrites any previously specified MOR_OPTS option block.
  ///
  /// @param option_block parsed option block
  ///
  bool setMOROptsOptionBlock(const Util::OptionBlock &option_block)
  {
    morOptsOptionBlock_ = option_block;

    return true;
  }

public:
  AnalysisManager &                     analysisManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;

private:
  Util::OptionBlock     morAnalysisOptionBlock_;
  Util::OptionBlock     morOptsOptionBlock_;
};

// .MOR
struct MORAnalysisReg : public IO::PkgOptionsReg
{
  MORAnalysisReg(
    MORFactory &   factory )
    : factory_(factory)
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    factory_.setMORAnalysisOptionBlock(option_block);

    factory_.analysisManager_.addAnalysis(&factory_);

    return true;
  }

  MORFactory  &               factory_;
};

//-----------------------------------------------------------------------------
// Function      : extractMORData
// Purpose       : Extract the parameters from a netlist .MOR line held in
//                 parsed_line.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 5/25/12
//-----------------------------------------------------------------------------
bool
extractMORData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("MOR", Util::OptionBlock::NO_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  int linePosition = 0;   // Start of parameters on .param line.

  Util::Param parameter("", "");

  if (numFields > 1)
  {
    std::vector<std::string> portNames(numFields - 1);
    for (int i=0; i<(numFields-1); ++i)
    {
      ++linePosition;     // Advance to next parameter.

      // Insert next port name
      ExtendedString stringVal ( parsed_line[linePosition].string_ );
      stringVal.toUpper();
      portNames[i] = stringVal;
    }

    parameter.setTag( "PORTLIST" );
    parameter.setVal( portNames );
    option_block.addParam( parameter );
  }

  circuit_block.addOptions(option_block);

  return true;
}

void
populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("MOR_OPTS");

  parameters.insert(Util::ParamMap::value_type("SIZE", Util::Param("SIZE", -1)));
  parameters.insert(Util::ParamMap::value_type("METHOD", Util::Param("METHOD", "PRIMA")));
  parameters.insert(Util::ParamMap::value_type("AUTOSIZE", Util::Param("AUTOSIZE", false)));
  parameters.insert(Util::ParamMap::value_type("ANALYSISONLY", Util::Param("ANALYSISONLY", false)));
  parameters.insert(Util::ParamMap::value_type("MAXSIZE", Util::Param("MAXSIZE", -1)));
  parameters.insert(Util::ParamMap::value_type("MAXFREQ", Util::Param("MAXFREQ", 1.0e9)));
  parameters.insert(Util::ParamMap::value_type("SAVEREDSYS", Util::Param("SAVEREDSYS", false)));
  parameters.insert(Util::ParamMap::value_type("COMPORIGTF", Util::Param("COMPORIGTF", false)));
  parameters.insert(Util::ParamMap::value_type("COMPREDTF", Util::Param("COMPREDTF", false)));
  parameters.insert(Util::ParamMap::value_type("COMPTYPE", Util::Param("COMPTYPE", "DEC")));
  parameters.insert(Util::ParamMap::value_type("COMPNP", Util::Param("COMPNP", 10)));
  parameters.insert(Util::ParamMap::value_type("COMPFSTART", Util::Param("COMPFSTART", 1.0)));
  parameters.insert(Util::ParamMap::value_type("COMPFSTOP", Util::Param("COMPFSTOP", 1.0)));
  parameters.insert(Util::ParamMap::value_type("EXPPOINT", Util::Param("EXPPOINT", 0.0)));
  parameters.insert(Util::ParamMap::value_type("SCALETYPE", Util::Param("SCALETYPE", 0)));
  parameters.insert(Util::ParamMap::value_type("SCALEFACTOR", Util::Param("SCALEFACTOR", 1)));
  parameters.insert(Util::ParamMap::value_type("SCALEFACTOR1", Util::Param("SCALEFACTOR1", 0.01)));
  parameters.insert(Util::ParamMap::value_type("SPARSIFICATIONTYPE", Util::Param("SPARSIFICATIONTYPE", 0)));
  parameters.insert(Util::ParamMap::value_type("SUBCKTS", Util::Param("SUBCKTS", "VECTOR")));
}

} // namespace <unnamed>


bool
registerMORFactory(
  FactoryBlock &        factory_block)
{
  MORFactory *factory = new MORFactory(factory_block.analysisManager_, factory_block.linearSystem_, factory_block.nonlinearManager_, factory_block.loader_, factory_block.topology_, factory_block.initialConditionsManager_);

  addAnalysisFactory(factory_block, factory);

  populateMetadata(factory_block.optionsManager_);

  factory_block.optionsManager_.addCommandParser(".MOR", extractMORData);

  factory_block.optionsManager_.addCommandProcessor("MOR", new MORAnalysisReg(*factory));

  factory_block.optionsManager_.addOptionsProcessor("MOR_OPTS", IO::createRegistrationOptions(*factory, &MORFactory::setMOROptsOptionBlock));

  return true;
}

} // namespace Analysis
} // namespace Xyce

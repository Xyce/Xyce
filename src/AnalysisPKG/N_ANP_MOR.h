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
//
// Purpose        : MOR analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Ting Mei   
//
// Creation Date  : 01/11
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_MOR_h
#define Xyce_N_ANP_MOR_h

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;
#include <Teuchos_SerialDenseMatrix.hpp>

// ----------   Xyce Includes   ----------
#include <N_ANP_fwd.h>
#include <N_UTL_fwd.h>
#include <N_LAS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_PDS_fwd.h>

#include <N_ANP_AnalysisBase.h>
#include <N_ANP_RegisterAnalysis.h>
#include <N_IO_OutputMOR.h>
#include <N_UTL_FixedQueue.h>
#include <N_UTL_OptionBlock.h>

// ---------- Forward Declarations ----------
class Amesos_BaseSolver;
class Epetra_LinearProblem;

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : MOR 
// Purpose       : MOR analysis class
// Special Notes : 
// Creator       : Ting Mei
// Creation Date : 05/24/11
//-------------------------------------------------------------------------
class MOR: public AnalysisBase 
{
public:
  MOR(
    AnalysisManager &                   analysis_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager);

    ~MOR();
   
  const TimeIntg::TIAParams &getTIAParams() const
  {
    return tiaParams_;
  }

  TimeIntg::TIAParams &getTIAParams()
  {
    return tiaParams_;
  }

    // Method to set MOR options
    bool setMOROptions(const Util::OptionBlock & option_block);

    bool setAnalysisParams(const Util::OptionBlock & paramsBlock);

    bool reduceSystem();
    bool evalOrigTransferFunction();
    bool evalRedTransferFunction();

    bool processSuccessfulStep(bool origSys);

protected:
    bool doRun();
    bool doInit();
    virtual bool doProcessSuccessfulStep() { return true; }
    virtual bool doLoopProcess() { return true; }
    virtual bool doProcessFailedStep();
    virtual bool doFinish();
    virtual bool doHandlePredictor();

public:
  bool getDCOPFlag() const
  {
    return true;
  }

private:
  Parallel::Machine                     comm_;

  AnalysisManager &                     analysisManager_;
  Loader::Loader &                      loader_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  OutputMgrAdapter &                    outputManagerAdapter_;
  IO::OutputMOR                         outputMOR_;

  TimeIntg::TIAParams   tiaParams_;
  
  int ROMsize_;
    std::string morMethod_;
    bool morSaveRedSys_;
    bool morCompOrigTF_;
    bool morCompRedTF_;
    std::string morCompType_;
    int morCompNP_;
    double morCompFStart_;
    bool morAutoSize_;
    int morMaxSize_;
    double morMaxFreq_;

    double morCompFStop_;
    double morExpPoint_;
    double morScaleFactor_;
    int morScaleType_;
    double morScaleFactor1_;
    int morSparsificationType_;
    bool isROMSparse_;
    std::vector<std::string> subcircuitNames_; 
 
    bool isPaused; 
    int morEvalSize_;
    int numPorts_;
 
    std::list < int > morEvalFailures_;

    bool isSingleFreq_;
    std::vector<std::string> portList_;
 
    double stepMult_;
    double fStep_;
    double currentFreq_; 
    double s0_;

    int setupSweepParam_();

    bool updateCurrentFreq_(int stepNumber);
 
    bool createOrigLinearSystem_();
    bool createRedLinearSystem_();
    
    bool updateOrigLinearSystemFreq_();
    bool updateRedLinearSystemFreq_();

    bool solveOrigLinearSystem_();
    bool solveRedLinearSystem_();

    bool sparsifyRedSystem_();

    // Original system
    RCP<Linear::Matrix> CPtr_;
    RCP<Linear::Matrix> GPtr_;
    RCP<Linear::Matrix> sCpG_MatrixPtr_;
    RCP<Linear::MultiVector> RPtr_, BPtr_, VPtr_;
    std::vector<int> bMatEntriesVec_, bMatPosEntriesVec_;

    // Original system, real-equivalent form
    RCP<Linear::BlockMatrix> sCpG_REFMatrixPtr_;
    RCP<Linear::BlockVector> REFBPtr_;
    RCP<Linear::BlockVector> REFXPtr_; // Store solution from Amesos here.

    // Reduced system (dense)
    Teuchos::SerialDenseMatrix<int, double> redC_;
    Teuchos::SerialDenseMatrix<int, double> redG_;
    Teuchos::SerialDenseMatrix<int, double> redB_;
    Teuchos::SerialDenseMatrix<int, double> redL_;  // redL_ != redB_

    // Reduced system (sparse)
    RCP<Linear::Matrix> redCPtr_, redGPtr_;
    RCP<Parallel::ParMap> redMapPtr_;

    // Reduced system, real-equivalent form (dense)
    Teuchos::SerialDenseMatrix<int, double> sCpG_redMatrix_, sCpG_tmpMatrix_;
    Teuchos::SerialDenseMatrix<int, double> ref_redB_;

    // Reduced system, real-equivalent form (sparse)
    RCP<Linear::BlockMatrix> sCpG_ref_redMatrixPtr_;
    RCP<Linear::BlockVector> ref_redBPtr_;
    RCP<Linear::BlockVector> ref_redXPtr_; // Store solution from Amesos here.

    // Transfer functions
    Teuchos::SerialDenseMatrix<int, std::complex<double> > origH_;
    Teuchos::SerialDenseMatrix<int, std::complex<double> > redH_;

    // Original system solver objects
    RCP<Amesos_BaseSolver> blockSolver_, origSolver_;
    RCP<Epetra_LinearProblem> blockProblem_, origProblem_;

    // Reduced system solver objects (sparse)
    RCP<Amesos_BaseSolver> blockRedSolver_;
    RCP<Epetra_LinearProblem> blockRedProblem_;
};

bool registerMORFactory(FactoryBlock &factory_block);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_MOR_h

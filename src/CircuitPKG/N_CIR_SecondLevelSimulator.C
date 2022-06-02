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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/27/00
//
//
//
//
//-------------------------------------------------------------------------
//
#include <iostream>

#include <Xyce_config.h>

#include <N_CIR_SecondLevelSimulator.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_ExternalSimulationData.h>
#include <N_LAS_System.h>
#include <N_LAS_Matrix.h>
#include <N_LOA_CktLoader.h>
#include <N_NLS_Manager.h>
#include <N_NLS_ConductanceExtractor.h>
#include <N_NLS_TwoLevelPrintJac.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_WorkingIntegrationMethod.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Circuit {

Analysis::AnalysisManager *
SecondLevelSimulator::newAnalysisManager(
  const IO::CmdParse &                command_line,
  IO::RestartMgr &                    restart_manager,
  Analysis::OutputMgrAdapter &        output_manager_adapter,
  Stats::Stat                         analysis_stat)
{
  return secondLevelManager_ = new Analysis::SecondLevelManager(command_line, output_manager_adapter, analysis_stat);
}


// ---------------------------------------------------------------------------
// API METHODS NEEDED FOR Two-level Functions:
// ---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::doRegistrations_
// Purpose       : This is a derived function.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/2022
//---------------------------------------------------------------------------
bool SecondLevelSimulator::doRegistrations_()
{
  bool success2 = Xyce::Analysis::registerTwoLevelPkgOptionsMgr(*secondLevelManager_, *(getOptionsManager()) );
  bool success = Simulator::doRegistrations_();
  return (success && success2);
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/01/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::simulateStep(
  bool                                  external_initJctFlag,
  const std::map<std::string, double> & inputMap,
  std::vector<double> &                 outputVector,
  std::vector< std::vector<double> > &  jacobian,
  TimeIntg::TwoLevelError &             tlError)
{
  bool bsuccess = false;

  if (DEBUG_CIRCUIT)
    dout() << "\nsimulateStep: " << std::endl;

  // Apply the input voltages to their appropriate sources:
  for (std::map<std::string,double>::const_iterator it = inputMap.begin(), end = inputMap.end(); it != end; ++it)
  {
    std::string name = (*it).first;
    getDeviceManager().setParam(name, (*it).second);
  }

  secondLevelManager_->setExternalSolverState(getCircuitLoader(), external_initJctFlag);

  bsuccess = secondLevelManager_->runSecondLevelStep(tlError);

  // calculate the conductance:
  getNonlinearManager().getConductanceExtractor().extract(inputMap, outputVector, jacobian);

  // save a copy of these for use in outputs, etc.
  inputMap_ = inputMap;
  outputVector_ = outputVector;
  jacobian_ = jacobian;

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::simulateStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 07/16/2009
//---------------------------------------------------------------------------
bool SecondLevelSimulator::simulateStep(
  const Device::ExternalSimulationData &        ext_data,
  const std::map<std::string, double> &         inputMap,
  std::vector<double> &                         outputVector,
  std::vector< std::vector<double> > &          jacobian,
  TimeIntg::TwoLevelError &                     tlError)
{
  bool bsuccess = false;

  if (DEBUG_CIRCUIT)
    dout() << "\nsimulateStep: " << std::endl;

  bsuccess = simulateStep(false, inputMap, outputVector, jacobian, tlError);

  return bsuccess;
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::startupSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::startupSolvers()
{
  return secondLevelManager_->startupSecondLevelSolvers(getLinearSystem(), getNonlinearManager());
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::finishSolvers
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/10/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::finishSolvers ()
{
  return secondLevelManager_->finishSecondLevelSolvers();
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::homotopyStepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
void
SecondLevelSimulator::homotopyStepSuccess(
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals)
{
  secondLevelManager_->homotopyStepSuccess(paramNames, paramVals);
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::homotopyStepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/2006
//---------------------------------------------------------------------------
void SecondLevelSimulator::homotopyStepFailure()
{
  secondLevelManager_->homotopyStepFailure();
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::daeOutputs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2022
//---------------------------------------------------------------------------
void SecondLevelSimulator::daeOutputs()
{
  Xyce::TimeIntg::DataStore* dsPtr_ = secondLevelManager_->getDataStore();

  bool outputDAEvectors = secondLevelManager_->getOutputDAEvectors();
  bool outputDAEvectors_noport = secondLevelManager_->getOutputDAEvectors_noport();
  bool outputDAEmatrices = secondLevelManager_->getOutputDAEmatrices();

  std::string netlistFile;
  int outputStepNum=0;
  if (outputDAEvectors || outputDAEvectors_noport || outputDAEmatrices)
  {
    netlistFile = commandLine_.getArgumentValue("netlist");
    outputStepNum = secondLevelManager_->getStepNumber() + 1;
  }

  if (outputDAEvectors || outputDAEvectors_noport)
  {
    Linear::Vector *daeQ    = dsPtr_->daeQVectorPtr;
    Linear::Vector *daeF    = dsPtr_->daeFVectorPtr;
    Linear::Vector *daeB    = dsPtr_->daeBVectorPtr;
    Linear::Vector *xVec    = dsPtr_->nextSolutionPtr;

    if (outputDAEvectors)
    {
      char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
      char filename4[256]; for (int ich = 0; ich < 256; ++ich) filename4[ich] = 0;
      char filename6[256]; for (int ich = 0; ich < 256; ++ich) filename6[ich] = 0;
      char filename6b[256];for (int ich = 0; ich < 256; ++ich) filename6b[ich] = 0;
      sprintf(filename1, "%s_soln_%03d.txt", netlistFile.c_str(), outputStepNum);
      sprintf(filename4, "%s_daeQ_%03d.txt", netlistFile.c_str(), outputStepNum);
      sprintf(filename6, "%s_daeF_%03d.txt", netlistFile.c_str(), outputStepNum);
      sprintf(filename6b,"%s_daeB_%03d.txt", netlistFile.c_str(), outputStepNum);

      // write the vectors:
      xVec->writeToFile(filename1);
      daeQ->writeToFile(filename4);
      daeF->writeToFile(filename6);
      daeB->writeToFile(filename6b);
    }

    if (outputDAEvectors_noport) // output F and B vectors, with port Vsrcs subtracted out
    {
      char filename7[256]; for (int ich = 0; ich < 256; ++ich) filename7[ich] = 0;
      char filename7b[256];for (int ich = 0; ich < 256; ++ich) filename7b[ich] = 0;
      sprintf(filename7, "%s_delF_%03d.txt", netlistFile.c_str(), outputStepNum);
      sprintf(filename7b,"%s_delB_%03d.txt", netlistFile.c_str(), outputStepNum);

      Xyce::Linear::Vector *delB = getLinearSystem().getRHSVector();// at this point, this isn't being used
      Xyce::Linear::Vector *delF = getLinearSystem().getNewtonVector();// at this point, this isn't being used
      Xyce::Linear::Vector *nextSol =  dsPtr_->nextSolutionPtr;
      std::vector<std::string> names;
      int numElectrodes = jacobian_.size();
      std::map<std::string,double>::const_iterator iterM = inputMap_.begin();
      std::map<std::string,double>::const_iterator  endM = inputMap_.end  ();
      for (int iE1 = 0; iE1 < numElectrodes; ++iE1,++iterM) { names.push_back(iterM->first); }
      delF->putScalar(0.0); 
      delB->putScalar(0.0); 
      getDeviceManager().loadTwoLevelVsrcs(names,delF,delB,nextSol);

      delF->scale(-1.0);
      delF->update(1.0,*daeF);

      delB->scale(-1.0);
      delB->update(1.0,*daeB);

      // write the vectors:
      delF->writeToFile(filename7);
      delB->writeToFile(filename7b);
    }
  }

  if (outputDAEmatrices)
  {
    char filename1[256]; for (int ich = 0; ich < 256; ++ich) filename1[ich] = 0;
    char filename2[256]; for (int ich = 0; ich < 256; ++ich) filename2[ich] = 0;
    sprintf(filename1, "%s_dQdx_%03d.txt", netlistFile.c_str(), outputStepNum);
    sprintf(filename2, "%s_dFdx_%03d.txt", netlistFile.c_str(), outputStepNum);

    // write the matrices:
    bool useLIDs=false;
    bool mmFormat=false;
    dsPtr_->dQdxMatrixPtr->writeToFile (filename1, useLIDs, mmFormat );
    dsPtr_->dFdxMatrixPtr->writeToFile (filename2, useLIDs, mmFormat );
  }
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::reducedOutputs
//
// Purpose       : This function outputs the extracted conductance 
//                 and the port currents to files.  It expects to be called 
//                 only at the end of converged time steps (or converged DC 
//                 steps).  In other words, from the stepSucccess function.
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 06/01/2022
//---------------------------------------------------------------------------
void SecondLevelSimulator::reducedOutputs()
{
  // outputs
  bool condOutput = secondLevelManager_->getCondOutputFlag();
  bool currOutput = secondLevelManager_->getPortCurrentOutputFlag();

  int outputStepNum(0);
  std::ofstream out;
  std::string netlistFile;
  char tmp[8]; for (int ich = 0; ich < 8; ++ich) tmp[ich] = 0;

  std::vector<std::string> names;
  if (condOutput || currOutput)
  {
    outputStepNum = secondLevelManager_->getStepNumber();
    netlistFile = commandLine_.getArgumentValue("netlist");
    sprintf(tmp, "%03d", outputStepNum);

    int numElectrodes = jacobian_.size();
    std::map<std::string,double>::const_iterator iterM = inputMap_.begin();
    std::map<std::string,double>::const_iterator  endM = inputMap_.end  ();
    for (int iE1 = 0; iE1 < numElectrodes; ++iE1,++iterM) { names.push_back(iterM->first); }
  }

  if (condOutput)
  {
    std::string condFileName = netlistFile + "_conductance_" + std::string(tmp) + ".txt";
    out.open(condFileName);
    Xyce::Nonlinear::printJacobian(out,netlistFile,names,jacobian_);
    out.close();  
  }

  if (currOutput)
  {
    int colw=20;
    std::string currFileName = netlistFile + "_port_currents_" + std::string(tmp) + ".txt";
    out.open(currFileName);
    for (int ii=0;ii<outputVector_.size();++ii)
    {
      out << netlistFile << " port:"<< std::setw(15) << names[ii];
      out << std::scientific << std::setw(colw) << std::setprecision(8) <<  outputVector_[ii];
      out << std::endl;
    }
    out.close();  
  }
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::stepSuccess
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void SecondLevelSimulator::stepSuccess(Analysis::TwoLevelMode analysis)
{
  daeOutputs();
  secondLevelManager_->stepSecondLevelSuccess(analysis);
  reducedOutputs();
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::stepFailure
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
void SecondLevelSimulator::stepFailure(Analysis::TwoLevelMode analysis)
{
  secondLevelManager_->stepSecondLevelFailure(analysis);
}


//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::getInitialQnorm
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::getInitialQnorm (TimeIntg::TwoLevelError & tle)
{
  return secondLevelManager_->getSecondLevelInitialQnorm(tle);
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::getBreakPoints
// Purpose       :
// Special Notes : Used for 2-level Newton solves.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::getBreakPoints (
    std::vector<Util::BreakPoint> &breakPointTimes,
    std::vector<Util::BreakPoint> &pauseBreakPointTimes)
{
  return secondLevelManager_->getSecondLevelBreakPoints(getCircuitLoader(),
      breakPointTimes, pauseBreakPointTimes);
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::updateStateArrays
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::updateStateArrays()
{
  return true;
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::setInternalParam
// Purpose       :
// Special Notes : Used for 2-level Newton solves, with LOCA.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/18/2006
//---------------------------------------------------------------------------
bool SecondLevelSimulator::setInternalParam (const std::string & name, double val)
{
  getDeviceManager().setParam(name, val);
  return true;
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//---------------------------------------------------------------------------
bool
SecondLevelSimulator::startTimeStep(
  bool                          beginIntegrationFlag,
  double                        nextTimeStep,
  double                        nextTime,
  int                           currentOrder)
{
  return secondLevelManager_->startSecondLevelTimeStep(
    getAnalysisManager().getTIAParams(),
    getNonlinearManager(),
    beginIntegrationFlag,
    nextTimeStep,
    nextTime,
    currentOrder);
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::startTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
//
//                 ERK, 2018: adapted so it could be potentially used
//                 with Xyce-to-Xyce coupling. The motivation is so
//                 that there only be a single set of API functions,
//                 instead of the confusing 2.
//
//                 This adaptation has to do with two bits of information:
//
//                 (1) "beginIntegrationFlag" which is true not just at the first
//                      time step, but also at any breakpoint step.
//
//                 (2) setting the integration order.  For the Xyce-to-Xyce case,
//                     all Xyce objects use the same synchronized integration
//                     order (and integration method)
//
//                 These two issues are invoked optionally.
//
// Scope         : public
// Creator       : Roger Pawlowski, SNL
// Creation Date : 08/28/2009
//---------------------------------------------------------------------------
bool SecondLevelSimulator::startTimeStep(const Device::ExternalSimulationData & ext_data)
{
  bool beginIntegrationFlag = false;
  double nextTimeStep = 0.0;
  double nextTime = 0.0;
  int currentOrder = 1;

  if (ext_data.is_transient)
  {
    nextTimeStep = ext_data.current_time_step_size;
    nextTime = ext_data.current_time;
    beginIntegrationFlag = (ext_data.time_step_number == 0);
  }

  // forcing the order is important for Xyce-to-Xyce coupling, where all
  // Xyce objects use the same integration method and order.
  if (ext_data.forceOrder)
  {
    currentOrder = ext_data.imposedTimeIntegrationOrder;
  }

  // the begin integration flag is "true" when at a breakpoint, or at the
  // first time step out of the DCOP.
  if (ext_data.forceBeginningIntegration)
  {
    beginIntegrationFlag = ext_data.imposedBeginningIntegration;
  }

  return secondLevelManager_->startSecondLevelTimeStep(
    getAnalysisManager().getTIAParams(),
    getNonlinearManager(),
    beginIntegrationFlag,
    nextTimeStep,
    nextTime,
    currentOrder);
}

//---------------------------------------------------------------------------
// Function      : SecondLevelSimulator::endTimeStep
// Purpose       :
// Special Notes : Used for 2-level Newton solves with Charon.
// Scope         : public
// Creator       : Russell Hooper, SNL
// Creation Date : 10/16/2012
//---------------------------------------------------------------------------
bool SecondLevelSimulator::endTimeStep (Device::ExternalSimulationData & ext_data)
{
  // We could opt to obtain selected time integration data here or
  // allow it to be obtained higher up, eg in charon::sc::CircuitDriver
  bool bsuccess = true;

  ext_data.currentOrder         = getAnalysisManager().getWorkingIntegrationMethod().getOrder();                // tiInfo.currentOrder        ;
  ext_data.numberOfSteps        = getAnalysisManager().getWorkingIntegrationMethod().getNumberOfSteps();        // tiInfo.numberOfSteps       ;
  ext_data.usedOrder            = getAnalysisManager().getWorkingIntegrationMethod().getUsedOrder();            // tiInfo.usedOrder           ;
  ext_data.nscsco               = getAnalysisManager().getWorkingIntegrationMethod().getNscsco();               // tiInfo.nscsco              ;
  ext_data.pdt                  = getAnalysisManager().getWorkingIntegrationMethod().partialTimeDeriv();        // tiInfo.pdt                 ;
  ext_data.nextTimeStep         = getAnalysisManager().getStepErrorControl().currentTimeStep;                   // tiInfo.nextTimeStep        ;
  ext_data.currTimeStep         = getAnalysisManager().getStepErrorControl().lastTimeStep;                      // tiInfo.currTimeStep        ;
  ext_data.currentTime          = getAnalysisManager().getStepErrorControl().currentTime;                       // tiInfo.currentTime         ;
  ext_data.nextTime             = getAnalysisManager().getStepErrorControl().nextTime;                          // tiInfo.nextTime            ;
  ext_data.beginIntegrationFlag = getAnalysisManager().getBeginningIntegrationFlag();                           // tiInfo.beginIntegrationFlag;
  ext_data.finalTime            = getAnalysisManager().getStepErrorControl().finalTime;                         // tiInfo.finalTime           ;
  ext_data.startingTimeStep     = getAnalysisManager().getStepErrorControl().startingTimeStep;                  // tiInfo.startingTimeStep    ;
  ext_data.bpTol                = getAnalysisManager().getStepErrorControl().getBreakPointLess().tolerance_;    // tiInfo.bpTol               ;
  ext_data.dcopFlag             = getAnalysisManager().getAnalysisObject().getDCOPFlag();                       // tiInfo.dcopFlag            ;
  ext_data.acopFlag             = getAnalysisManager().getACOPFlag();                                           // tiInfo.acopFlag            ;
  ext_data.inputOPFlag          = getAnalysisManager().getAnalysisObject().getInputOPFlag();                    // tiInfo.inputOPFlag         ;
  ext_data.tranopFlag           = getAnalysisManager().getTranOPFlag();                                         // tiInfo.tranopFlag          ;
  ext_data.transientFlag        = getAnalysisManager().getTransientFlag();                                      // tiInfo.transientFlag       ;
  ext_data.dcsweepFlag          = getAnalysisManager().getDCSweepFlag();                                        // tiInfo.dcsweepFlag         ;
  ext_data.timeStepNumber       = getAnalysisManager().getStepNumber();                                         // tiInfo.timeStepNumber      ;
  ext_data.initTranFlag         = getAnalysisManager().getInitTranFlag();                                       // tiInfo.initTranFlag        ;
  ext_data.sweepSourceResetFlag = getAnalysisManager().getSweepSourceResetFlag();                               // tiInfo.sweepSourceResetFlag;
  ext_data.timeIntMode          = getAnalysisManager().getAnalysisObject().getIntegrationMethod();              // tiInfo.timeIntMode         ;
  ext_data.doubleDCOPStep       = getAnalysisManager().getDoubleDCOPStep();                                     // tiInfo.doubleDCOPStep      ;
  ext_data.doubleDCOPEnabled    = getAnalysisManager().getDoubleDCOPEnabled();                                  // tiInfo.doubleDCOPEnabled   ;

  return bsuccess;
}

} // namespace Circuit
} // namespace Xyce

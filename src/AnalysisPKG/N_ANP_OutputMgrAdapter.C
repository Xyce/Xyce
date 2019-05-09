//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_ANP_OutputMgrAdapter.h>

#include <N_ANP_NoiseData.h>

#include <N_DEV_DeviceMgr.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_InitialConditions.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputROM.h>
#include <N_IO_RestartMgr.h>
#include <N_DEV_Op.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::OutputMgrAdapter
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
OutputMgrAdapter::OutputMgrAdapter(
  Parallel::Machine                     comm,
  IO::OutputMgr &                       output_manager,
  IO::Measure::Manager &                measure_manager,
  IO::FourierMgr &                      fourier_manager,
  Device::DeviceMgr &                   device_manager)
  : comm_(comm),
    outputManager_(output_manager),
    measureManager_(measure_manager),
    fourierManager_(fourier_manager),
    deviceManager_(device_manager),
    tempOp_(new Device::ArtificialParameterOp("TEMP", deviceManager_, *(*deviceManager_.getArtificialParameterMap().find("TEMP")).second, "TEMP")),
    stepSweepVector_(),
    dcSweepVector_(),
    stepAnalysisStepNumber_(0),
    stepAnalysisMaxSteps_(0),
    dcAnalysisStepNumber_(0),
    dcAnalysisMaxSteps_(0)
{}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::~OutputMgrAdapter( )
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
OutputMgrAdapter::~OutputMgrAdapter()
{
  delete tempOp_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::notify
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgrAdapter::notify(
  const StepEvent &     event)
{
  if (event.state_ == StepEvent::STEP_STARTED)
    stepAnalysisStepNumber_ = event.count_;
  else if (event.state_ == StepEvent::INITIALIZE)
    stepAnalysisMaxSteps_ = event.count_;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::setStepSweepVector
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::setStepSweepVector(
  const Analysis::SweepVector &       sweep_vector)
{
  stepSweepVector_ = sweep_vector;

  // ERK this call is necessary b/c the homotopy outputters
  // don't get the sweep vector as a function argument.
  outputManager_.setStepSweepVector(sweep_vector);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::setDCSweepVector
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::setDCSweepVector(
  const Analysis::SweepVector &       sweep_vector)
{
  dcSweepVector_ = sweep_vector;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::dumpRestart
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgrAdapter::dumpRestart(
  Parallel::Communicator &      parallel_communicator,
  Topo::Topology &              topology,
  Analysis::AnalysisManager &   analysis_manager,
  const std::string &           job_name,
  bool                          pack,
  double                        current_time) const
{
  IO::dumpRestartData(
    parallel_communicator,
    topology,
    analysis_manager,
    deviceManager_,
    job_name,
    pack,
    current_time);
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::tranOutput
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::tranOutput(
  double                time,
  double                dt,
  Linear::Vector &      currSolutionPtr,
  Linear::Vector &      stateVecPtr,
  Linear::Vector &      storeVecPtr,
  Linear::Vector &      lead_current_vector,
  Linear::Vector &      junction_voltage_vector,
  Linear::Vector &      lead_current_dqdt_vector,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_,
  bool                  skipPrintLineOutput)
{
  fourierManager_.updateFourierData(comm_, time, &currSolutionPtr, &stateVecPtr, &storeVecPtr, &lead_current_vector, &junction_voltage_vector, &lead_current_dqdt_vector,
    objectiveVec_, dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_
      );

  measureManager_.updateTranMeasures(comm_, time, &currSolutionPtr, &stateVecPtr, &storeVecPtr, &lead_current_vector, &junction_voltage_vector, &lead_current_dqdt_vector);

  Util::Op::OpData op_data;
  double temp = (*tempOp_)(comm_, op_data).real();

  outputManager_.output(
    comm_, time, dt, temp,
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    dcAnalysisStepNumber_, dcAnalysisMaxSteps_, dcSweepVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector, 
    objectiveVec_,
    dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_,
    skipPrintLineOutput);

}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::tranSensitivityOutput
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2016
//-----------------------------------------------------------------------------
void OutputMgrAdapter::tranSensitivityOutput(
  double                time,
  double                dt,
  Linear::Vector &      currSolutionPtr,
  Linear::Vector &      stateVecPtr,
  Linear::Vector &      storeVecPtr,
  Linear::Vector &      lead_current_vector,
  Linear::Vector &      junction_voltage_vector,
  Linear::Vector &      lead_current_dqdt_vector,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_,
  bool                  skipPrintLineOutput)
{
  Util::Op::OpData op_data;
  double temp = (*tempOp_)(comm_, op_data).real();

  outputManager_.outputSensitivity(
    comm_, time, dt, temp, 
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    dcAnalysisStepNumber_, dcAnalysisMaxSteps_, dcSweepVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector, objectiveVec_,
    dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_,
    skipPrintLineOutput);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::dcOutput
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::dcOutput(
  int                   dcStepNumber,
  Linear::Vector &      currSolutionPtr,
  Linear::Vector &      stateVecPtr,
  Linear::Vector &      storeVecPtr,
  Linear::Vector &      lead_current_vector,
  Linear::Vector &      junction_voltage_vector,
  Linear::Vector &      lead_current_dqdt_vector,
  std::vector<double> & objectiveVec_,
  std::vector<double> & dOdpVec_,
  std::vector<double> & dOdpAdjVec_,
  std::vector<double> & scaled_dOdpVec_,
  std::vector<double> & scaled_dOdpAdjVec_)
{
  measureManager_.updateDCMeasures(comm_, dcSweepVector_, &currSolutionPtr, 
      &stateVecPtr, &storeVecPtr, &lead_current_vector, &junction_voltage_vector, 
      &lead_current_dqdt_vector);

  Util::Op::OpData op_data;
  double temp = (*tempOp_)(comm_, op_data).real();

  outputManager_.output(
    comm_, 0.0, 0.0, temp, 
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    dcStepNumber, dcAnalysisMaxSteps_, dcSweepVector_,
    currSolutionPtr, stateVecPtr, storeVecPtr, 
    lead_current_vector, junction_voltage_vector, lead_current_dqdt_vector, 
    objectiveVec_, dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);

  // dc values for the objective function call.
  double arg1 = 0.0;
  double arg2 = 0.0;

  if (dcSweepVector_.size() > 0)
  {
    arg1 = dcSweepVector_[0].currentVal;
  }
  if (dcSweepVector_.size() > 1)
  {
    arg2 = dcSweepVector_[1].currentVal;
  }

}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::finishOutput
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::finishOutput()
{
  outputManager_.finishOutput();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::finishSensitivityOutput
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 05/3/2016
//-----------------------------------------------------------------------------
void OutputMgrAdapter::finishSensitivityOutput()
{
  outputManager_.finishSensitivityOutput();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputMPDE
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputMPDE(
  double                        time,
  const std::vector<double> &   fast_time_points,
  const Linear::BlockVector &   solution_vector)
{
  outputManager_.outputMPDE(comm_, time, fast_time_points, solution_vector);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputHB_TD
// Purpose       : constructor
// Special Notes : Used for HB time-domain output such as .PRINT HB_TD lines.
//                 This is not used for .PRINT HB_STARTUP or .PRINT HB_IC 
//                 lines though.
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputHB_TD (
  const std::vector<double> &   timePoints,
  const Linear::BlockVector &   timeDomainSolutionVec,
  const Linear::BlockVector &   timeDomainLeadCurrentVec,
  const Linear::BlockVector &   timeDomainJunctionVoltageVec)

{
  outputManager_.outputHB_TD(
    comm_,
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    timePoints,
    timeDomainSolutionVec, 
    timeDomainLeadCurrentVec, 
    timeDomainJunctionVoltageVec);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputHB_FD
// Purpose       : constructor
// Special Notes : Used for HB frequency-domain output such as .PRINT HB_FD lines.
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputHB_FD (
  const std::vector<double> &   freqPoints,
  const Linear::BlockVector &   freqDomainSolutionVecReal,
  const Linear::BlockVector &   freqDomainSolutionVecImaginary,
  const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
  const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary )

{
  outputManager_.outputHB_FD(
    comm_,
    stepAnalysisStepNumber_, stepAnalysisMaxSteps_, stepSweepVector_,
    freqPoints,
    freqDomainSolutionVecReal, freqDomainSolutionVecImaginary, 
    freqDomainLeadCurrentVecReal, freqDomainLeadCurrentVecImaginary,
    freqDomainJunctionVoltageVecReal, freqDomainJunctionVoltageVecImaginary);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputAC
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputAC (
  double                frequency,
  double                fStart,
  double                fStop,
  const Linear::Vector &  solnVecRealPtr,
  const Linear::Vector &  solnVecImaginaryPtr,
  const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams)
{
  measureManager_.updateACMeasures(comm_, frequency, &solnVecRealPtr, &solnVecImaginaryPtr);
  
  outputManager_.outputAC(comm_, frequency, fStart, fStop, solnVecRealPtr, solnVecImaginaryPtr, Sparams);

}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputSensitivityAC
// Purpose       : constructor
// Special Notes : This is used to output the sensitivity information for
//                 .AC analyses
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 4/15/2019
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputSensitivityAC (
     double frequency,
     const Linear::Vector & solnVecRealPtr,
     const Linear::Vector & solnVecImaginaryPtr,
     const std::vector<double> & paramVals,
     const std::vector<std::string> & paramNameVec,
     const std::vector<std::string> & objFuncVars,
     const std::vector<double> & objectiveVec,
     const std::vector<double> & dOdpVec,
     const std::vector<double> & dOdpAdjVec,
     const std::vector<double> & scaled_dOdpVec,
     const std::vector<double> & scaled_dOdpAdjVec)
{
  outputManager_.outputSensitivityAC(comm_, frequency, solnVecRealPtr, solnVecImaginaryPtr,
			  paramVals, paramNameVec, objFuncVars, objectiveVec,
                          dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputSParams
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 03/25/2019
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputSParams (
  double                frequency,
  double                numFreq,
  std::vector<double> & Z0sVec,
  const Teuchos::SerialDenseMatrix<int, std::complex<double> > & Sparams)
{
  outputManager_.outputSParams(comm_, frequency, numFreq, Z0sVec, Sparams);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputNoise
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/26/2015
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputNoise (
  double                frequency,
  const Linear::Vector &  solnVecRealPtr,
  const Linear::Vector &  solnVecImaginaryPtr,
  double totalOutputNoiseDens_, 
  double totalInputNoiseDens_, 
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_)
{
  
  outputManager_.outputNoise(comm_, frequency, solnVecRealPtr, solnVecImaginaryPtr, 
      totalOutputNoiseDens_, totalInputNoiseDens_,noiseDataVec_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::getInitialOutputInterval
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
double OutputMgrAdapter::getInitialOutputInterval() const
{
  return outputManager_.getInitialOutputInterval();
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::getOutputIntervals
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
const IO::IntervalVector & OutputMgrAdapter::getOutputIntervals() const
{
  return outputManager_.getOutputIntervals();
}


//-----------------------------------------------------------------------------
// Function      : OutputMgrAdapter::outputHomotopy
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 01/24/08
//-----------------------------------------------------------------------------
void OutputMgrAdapter::outputHomotopy(
  const std::vector<std::string> &      paramNames,
  const std::vector<double> &           paramVals,
  Linear::Vector &                        solnVecPtr )
{
  outputManager_.outputHomotopy (comm_, paramNames, paramVals, solnVecPtr);
}

} // namespace Analysis
} // namespace Xyce

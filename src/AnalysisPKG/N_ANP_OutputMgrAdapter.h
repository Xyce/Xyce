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

#ifndef Xyce_N_ANP_OutputMgrAdapter_h
#define Xyce_N_ANP_OutputMgrAdapter_h

#include <Teuchos_SerialDenseMatrix.hpp>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Listener.h>
#include <N_ANP_StepEvent.h>
#include <N_ANP_UQSupport.h>
#include <N_LAS_BlockVector.h>
#include <N_IO_OutputMgr.h>
#include <N_UTL_Op.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : OutputMgrAdapter
// Purpose       : Inteface class for the output manager
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class OutputMgrAdapter : public Util::Listener<StepEvent>
{
public:
  OutputMgrAdapter(
    Parallel::Machine                   comm,
    IO::OutputMgr &                     output_manager,
    IO::Measure::Manager &              measure_manager,
    IO::FourierMgr &                    fourier_manager,
    Device::DeviceMgr &                 device_manager);

  virtual ~OutputMgrAdapter();

  Parallel::Machine getComm()
  {
    return comm_;
  }

  IO::OutputMgr &getOutputManager()
  {
    return outputManager_;
  }

  IO::Measure::Manager &getMeasureManager()
  {
    return measureManager_;
  }

  const IO::AliasNodeMap &getAliasNodeMap() const
  {
    return outputManager_.getAliasNodeMap();
  }

  void setStepSweepVector(const Analysis::SweepVector &sweep_vector);

  void setDCSweepVector(const Analysis::SweepVector &sweep_vector);

  const Analysis::SweepVector &getStepSweepVector() const
  {
    return stepSweepVector_;
  }

  const Analysis::SweepVector &getDCSweepVector() const
  {
    return dcSweepVector_;
  }

  int getStepAnalysisStepNumber()
  {
    return stepAnalysisStepNumber_; 
  }

  int getStepAnalysisMaxSteps()
  {
    return stepAnalysisMaxSteps_;
  }

  void setDCAnalysisStepNumber( int num )
  {
    dcAnalysisStepNumber_ = num;
  }

  void setDCAnalysisMaxSteps( int num )
  {
    dcAnalysisMaxSteps_ = num;
  }

  void setDotACSpecified( bool value)
  {
    outputManager_.setDotACSpecified(value);
  }

  void setDotNoiseSpecified( bool value)
  {
    outputManager_.setDotNoiseSpecified(value);
  }

  void setEnableEmbeddedSamplingFlag( bool value)
  {
    outputManager_.setEnableEmbeddedSamplingFlag(value);
  }

  void setEnablePCEFlag( bool value)
  {
    outputManager_.setEnablePCEFlag(value);
  }

  // used to determine, for -r output, whether a
  // .LIN analysis is being done.
  void setEnableSparCalcFlag( bool value)
  {
    outputManager_.setEnableSparCalcFlag(value);
  }

  bool getPhaseOutputUsesRadians() const
  {
    return outputManager_.getPhaseOutputUsesRadians();
  }

  double getInitialOutputInterval() const;

  const IO::IntervalVector &getOutputIntervals() const;

  void notify(const StepEvent &event);

  void dumpRestart(
    Parallel::Communicator &    parallel_communicator,
    Topo::Topology &            topology,
    Analysis::AnalysisManager & analysis_manager,
    const std::string &         job_name,
    bool                        pack,
    double                      current_time) const;

  void tranOutput(
    double time, 
    double dt, 
    Linear::Vector & currSolutionPtr, 
    Linear::Vector & stateVecPtr, 
    Linear::Vector & storeVecPtr, 
    Linear::Vector & lead_current_vector,
    Linear::Vector & junction_voltage_vector,
    Linear::Vector & lead_current_dqdt_vector,
    std::vector<double> & objectiveVec_, 
    std::vector<double> & dOdpVec_, 
    std::vector<double> & dOdpAdjVec_,
    std::vector<double> & scaled_dOdpVec_, 
    std::vector<double> & scaled_dOdpAdjVec_,
    bool skipPrintLineOutput = false);


  void tranSensitivityOutput(
    double time, 
    double dt, 
    Linear::Vector & currSolutionPtr, 
    Linear::Vector & stateVecPtr, 
    Linear::Vector & storeVecPtr, 
    Linear::Vector & lead_current_vector,
    Linear::Vector & junction_voltage_vector,
    Linear::Vector & lead_current_dqdt_vector,
    std::vector<double> & objectiveVec_, 
    std::vector<double> & dOdpVec_, 
    std::vector<double> & dOdpAdjVec_,
    std::vector<double> & scaled_dOdpVec_, 
    std::vector<double> & scaled_dOdpAdjVec_,
    bool skipPrintLineOutput = false);

  void dcOutput( 
    int dcStepNumber, 
    Linear::Vector & currSolutionPtr, Linear::Vector & stateVecPtr, Linear::Vector & storeVecPtr,
    Linear::Vector &lead_current_vector,
    Linear::Vector &junction_voltage_vector,
    Linear::Vector &lead_current_dqdt_vector,
    std::vector<double> & objectiveVec_, 
    std::vector<double> & dOdpVec_, 
    std::vector<double> & dOdpAdjVec_, 
    std::vector<double> & scaled_dOdpVec_, 
    std::vector<double> & scaled_dOdpAdjVec_);

  void finishOutput();

  void finishSensitivityOutput();

  void outputMPDE(double time, const std::vector<double> &fast_time_points, const Linear::BlockVector &solution_vector);

  // Used for HB time-domain output such as .PRINT HB_TD lines.  This is
  // not used for .PRINT HB_STARTUP or .PRINT HB_IC lines though.
  void outputHB_TD(
    const std::vector< double > & timePoints,
    const Linear::BlockVector & timeDomainSolnVec,
    const Linear::BlockVector & timeDomainLeadCurrentVec, 
    const Linear::BlockVector & timeDomainJunctionVoltageVec);

  // Used for HB frequency-domain output such as .PRINT HB_FD lines.
  void outputHB_FD(
    const std::vector< double > & freqPoints,
    const Linear::BlockVector & freqDomainSolnVecReal,
    const Linear::BlockVector & freqDomainSolnVecImaginary, 
    const Linear::BlockVector & freqDomainLeadCurrentVecReal, 
    const Linear::BlockVector & freqDomainLeadCurrentVecImaginary, 
    const Linear::BlockVector & freqDomainJunctionVoltageVecReal, 
    const Linear::BlockVector & freqDomainJunctionVoltageVecImaginary );

  void outputAC(
     double freq, 
     double fStart,
     double fStop,
     const Linear::Vector & solnVecRealPtr, 
     const Linear::Vector & solnVecImaginaryPtr,
     const Util::Op::RFparamsData & RFparams);

  void outputSensitivityAC(
     double freq,
     const Linear::Vector & solnVecRealPtr,
     const Linear::Vector & solnVecImaginaryPtr,
     const std::vector<double> & paramVals,
     const std::vector<std::string> & paramNameVec,
     const std::vector<std::string> & objFuncVars,
     const std::vector<double> & objectiveVec,
     const std::vector<double> & dOdpVec,
     const std::vector<double> & dOdpAdjVec,
     const std::vector<double> & scaled_dOdpVec,
     const std::vector<double> & scaled_dOdpAdjVec);

   void outputSParams(
     double freq,
     double numFreq,
     std::vector<double> & Z0sVec,
     const Util::Op::RFparamsData & RFparams);

  void outputNoise(
      double freq,
      const Linear::Vector & solnVecRealPtr, 
      const Linear::Vector & solnVecImaginaryPtr, 
      double totalOutputNoiseDens_, 
      double totalInputNoiseDens_, 
      const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_);

  void outputEmbeddedSampling(
      bool regressionPCEenable,
      bool projectionPCEenable,
      int  numSamples,
      const std::vector<std::string> & regressionPCEcoeffs_,
      const std::vector<std::string> & projectionPCEcoeffs_,
      const std::vector<UQ::outputFunctionData*> & outFuncDataVec_ );

  void outputPCE(
      int numQuadPoints,
      const std::vector<UQ::outputFunctionData*> & outFuncDataVec_ );

  void outputHomotopy( const std::vector<std::string> & paramNames, const std::vector<double> & paramVals, Linear::Vector & solnVecPtr );

private:
  Parallel::Machine                     comm_;
  IO::OutputMgr &                       outputManager_;
  IO::Measure::Manager &                measureManager_;
  IO::FourierMgr &                      fourierManager_;
  Device::DeviceMgr &                   deviceManager_;

  Util::Op::Operator *                  tempOp_;

  Analysis::SweepVector                 stepSweepVector_;
  Analysis::SweepVector                 dcSweepVector_;

  int stepAnalysisStepNumber_;
  int stepAnalysisMaxSteps_;
  int dcAnalysisStepNumber_;
  int dcAnalysisMaxSteps_;
};

class OutputAdapter
{
public:
  OutputAdapter(OutputMgrAdapter &adapter)
    : outputManagerAdapter_(adapter)
  {}

  virtual void outputMPDE(double time, const Linear::Vector *solution_vector)
  {}

  virtual bool outputFunkyMPDE()
  {
    return false;
  }

protected:
  OutputMgrAdapter &    outputManagerAdapter_;
};

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_OutputMgrAdapter_h

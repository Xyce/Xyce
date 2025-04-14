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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <list>
#include <string>
#include <vector>

#include <N_IO_Outputter.h>
#include <N_ANP_fwd.h>
#include <N_ANP_UQSupport.h>

#include <N_UTL_Demangle.h>
#include <N_UTL_NetlistLocation.h>
#include <N_UTL_Param.h>

#include <Teuchos_SerialDenseMatrix.hpp>

namespace Xyce {
namespace IO {
namespace Outputter {

const static int debug = false;

//-----------------------------------------------------------------------------
// Function      : Interface::setAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::setAnalysisMode(
  Analysis::Mode        analysis_mode)
{
  doSetAnalysisMode(analysis_mode);
}

//-----------------------------------------------------------------------------
// Function      : Interface::output
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::output(
  Parallel::Machine             comm,
  const Linear::Vector &          solution_vector,
  const Linear::Vector &          state_vector,
  const Linear::Vector &          store_vector,
  const Linear::Vector &          lead_current_vector,
  const Linear::Vector &          junction_voltage_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutput" << std::endl;

  doOutputTime(comm, solution_vector, state_vector, store_vector, lead_current_vector, junction_voltage_vector);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputAC
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputAC(
  Parallel::Machine             comm,
  double                        frequency,
  double                        fStart,
  double                        fStop,
  const Linear::Vector &          real_solution_vector,
  const Linear::Vector &          imaginary_solution_vector,
  const Util::Op::RFparamsData & RFparams)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputAC" << std::endl;

  doOutputFrequency(comm, frequency, fStart, fStop, 
                    real_solution_vector, imaginary_solution_vector, RFparams);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputSensitivityAC
// Purpose       : This is used to output sensitivity information for .AC
//                 analyses.  It outputs both direct and adjoint info.
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 4/15/2019
//-----------------------------------------------------------------------------
void Interface::outputSensitivityAC(
    Parallel::Machine                 comm,
    double                            frequency,
    const Linear::Vector &            real_solution_vector,
    const Linear::Vector &            imaginary_solution_vector,
    const std::vector<double> &       paramVals,
    const std::vector<std::string> &  paramNameVec,
    const std::vector<std::string> &  objFuncVars,
    const std::vector<double> &       objectiveVec,
    const std::vector<double> &       dOdpVec,
    const std::vector<double> &       dOdpAdjVec,
    const std::vector<double> &       scaled_dOdpVec,
    const std::vector<double> &       scaled_dOdpAdjVec)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputSensitivityAC" << std::endl;

  doOutputSensitivityAC(comm, frequency, real_solution_vector, imaginary_solution_vector,
                        paramVals, paramNameVec, objFuncVars, objectiveVec,
                        dOdpVec, dOdpAdjVec, scaled_dOdpVec, scaled_dOdpAdjVec);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputSParams
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 3/25/19
//-----------------------------------------------------------------------------
void Interface::outputSParams(
  Parallel::Machine             comm,
  double                        frequency,
  double                        numFreq,
  std::vector<double> &         Z0sVec,
  const Util::Op::RFparamsData & RFparams)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputSParams" << std::endl;

  doOutputSParams(comm, frequency, numFreq, Z0sVec, RFparams);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputNoise
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputNoise(
  Parallel::Machine             comm,
  double                        frequency,
  const Linear::Vector &        real_solution_vector,
  const Linear::Vector &        imaginary_solution_vector,
  double              totalOutputNoiseDens_, 
  double              totalInputNoiseDens_, 
  const std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec_)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputNoise" << std::endl;

  doOutputNoise (comm, frequency, real_solution_vector, imaginary_solution_vector,
        totalOutputNoiseDens_, totalInputNoiseDens_, noiseDataVec_);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputEmbeddedSampling
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 7/29/2019
//-----------------------------------------------------------------------------
void Interface::outputEmbeddedSampling(
  Parallel::Machine             comm,
  bool                          regressionPCEenable,
  bool                          projectionPCEenable,
  int                           numSamples,
  const std::vector<std::string> & regressionPCEcoeffs,
  const std::vector<std::string> & projectionPCEcoeffs,
  const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec_)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputEmbeddedSampling" << std::endl;

  doOutputEmbeddedSampling(comm, regressionPCEenable, projectionPCEenable,
      numSamples, regressionPCEcoeffs, projectionPCEcoeffs, outFuncDataVec_);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputPCE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Pete Sholander, SNL
// Creation Date : 9/3/2019
//-----------------------------------------------------------------------------
void Interface::outputPCE(
  Parallel::Machine             comm,
  int                           numQuadPoints,
  const std::vector<Xyce::Analysis::UQ::outputFunctionData*> & outFuncDataVec_)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputPCE" << std::endl;

  doOutputPCE(comm, numQuadPoints, outFuncDataVec_);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputHB_TD
// Purpose       : Used for HB time-domain output such as .PRINT HB_TD lines.  
//                 This is not used for .PRINT HB_STARTUP or .PRINT HB_IC lines 
//                 though.
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputHB_TD(
  Parallel::Machine             comm,
  const std::vector<double> &   timePoints,
  const Linear::BlockVector &   timeDomainSolutionVec,
  const Linear::BlockVector &   timeDomainLeadCurrentVec,
  const Linear::BlockVector &   timeDomainJunctionVoltageVec)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputHB_TD" << std::endl;

  doOutputHB_TD(comm, timePoints,
             timeDomainSolutionVec,
             timeDomainLeadCurrentVec,
             timeDomainJunctionVoltageVec);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputHB_FD
// Purpose       : Used for HB frequency-domain output such as .PRINT HB_FD lines.
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputHB_FD(
  Parallel::Machine             comm,
  const std::vector<double> &   freqPoints,
  const Linear::BlockVector &   freqDomainSolutionVecReal,
  const Linear::BlockVector &   freqDomainSolutionVecImaginary,
  const Linear::BlockVector &   freqDomainLeadCurrentVecReal,
  const Linear::BlockVector &   freqDomainLeadCurrentVecImaginary,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecReal,
  const Linear::BlockVector &   freqDomainJunctionVoltageVecImaginary)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputHB_FD" << std::endl;

  doOutputHB_FD(comm, freqPoints,
             freqDomainSolutionVecReal, freqDomainSolutionVecImaginary,
             freqDomainLeadCurrentVecReal, freqDomainLeadCurrentVecImaginary,
             freqDomainJunctionVoltageVecReal, freqDomainJunctionVoltageVecImaginary);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputMPDE
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputMPDE(
  Parallel::Machine             comm,
  double                        time,
  const std::vector<double> &   fast_time_points,
  const Linear::BlockVector &     solution_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputMPDE" << std::endl;

  doOutputMPDE(comm, time, fast_time_points, solution_vector);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputHomotopy
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputHomotopy(
  Parallel::Machine                     comm,
  const std::vector<std::string> &      parameter_names,
  const std::vector<double> &           parameter_values,
  const Linear::Vector &                  solution_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputHomotopy" << std::endl;

  doOutputHomotopy(comm, parameter_names, parameter_values, solution_vector);
}

//-----------------------------------------------------------------------------
// Function      : Interface::outputSensitivity
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::outputSensitivity(
  Parallel::Machine             comm,
  const std::vector<double> &   objective_values,
  const std::vector<double> &   direct_values,
  const std::vector<double> &   adjoint_values,
  const std::vector<double> &   scaled_direct_values,
  const std::vector<double> &   scaled_adjoint_values,
  const Linear::Vector &          solution_vector,
  const Linear::Vector &          state_vector,
  const Linear::Vector &          store_vector)
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doOutputSensitivity" << std::endl;

  doOutputSensitivity(
    comm, objective_values,
    direct_values, adjoint_values,
    scaled_direct_values, scaled_adjoint_values,
    solution_vector, state_vector, store_vector);
}

//-----------------------------------------------------------------------------
// Function      : Interface::finishOutput
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::finishOutput()
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doFinishOutput" << std::endl;

  doFinishOutput();
}

//-----------------------------------------------------------------------------
// Function      : Interface::startStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::startStep(int step, int max_step)
{
  doStartStep(step, max_step);
}

//-----------------------------------------------------------------------------
// Function      : Interface::resetIndex
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::resetIndex()
{
  doResetIndex();
}

//-----------------------------------------------------------------------------
// Function      : Interface::steppingComplete
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Dave Baur, Raytheon
// Creation Date : 
//-----------------------------------------------------------------------------
void Interface::steppingComplete()
{
  if (debug) Xyce::dout() << demangle(typeid(*this).name()) << " doSteppingComplete" << std::endl;

  doSteppingComplete();
}


} // namespace Outputter
} // namespace IO
} // namespace Xyce

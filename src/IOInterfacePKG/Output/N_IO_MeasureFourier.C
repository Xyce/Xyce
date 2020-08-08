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
// Purpose       : Transient analysis functions.
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/05/2013
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <iomanip>

#include <N_IO_MeasureFourier.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Math.h>

#include <Teuchos_ScalarTraits.hpp>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : Fourier::Fourier()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/05/2013
//-----------------------------------------------------------------------------
Fourier::Fourier(const Manager &measureMgr, const Util::OptionBlock & measureBlock)
  : Base(measureMgr, measureBlock),
  calculated_(false)
{
  // indicate that this measure type is supported and should be processed in simulation
  typeSupported_ = true;

  // updateTran() is likely to segfault if the .MEASURE line was incomplete
  checkMeasureLine();

  // check for illegal values that will cause core dumps later in the code.
  // at_ is the fundamental frequency
  if (at_ < 0 || numFreq_ <= 0 || gridSize_ <= 0 )
  {
     Report::UserError0() << name_ << " has illegal value for AT, NUMFREQ or GRIDSIZE";
  }
}

//-----------------------------------------------------------------------------
// Function      : Fourier::prepareOutputVariables()
// Purpose       : Validates that the number of output variables is legal for this
//                 measure type, and then makes the vector for those variables.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void Fourier::prepareOutputVariables() 
{
  // this measurement should have only one dependent variable.
  // Error for now if it doesn't
  numOutVars_ = outputVars_.size();

  if ( numOutVars_ > 1 )
  {
    std::string msg = "Too many dependent variables for FOUR measure, \"" + name_ + "\"";
    Report::UserError0() << msg;
  }
}

//-----------------------------------------------------------------------------
// Function      : Fourier::reset()
// Purpose       : Called when restarting a measure function.  Resets any state
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 8/28/2014
//-----------------------------------------------------------------------------
void Fourier::reset() 
{
  resetBase();
  time_.clear();
  mag_.clear();
  phase_.clear();
  outVarValues_.clear();
  calculated_ = false;
}


//-----------------------------------------------------------------------------
// Function      : Fourier::updateTran()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/05/2013
//-----------------------------------------------------------------------------
void Fourier::updateTran(
  Parallel::Machine comm,
  const double circuitTime,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
  if( !calculationDone_ && withinTimeWindow( circuitTime ) )
  {
     // we're in the time window, now we need to store the time and output values
    if( !initialized_  )
    {      
      initialized_ = true;
    }

    time_.push_back(circuitTime);
    outVarValues_.push_back(getOutputValue(comm, outputVars_[0],
                                           solnVec, stateVec, storeVec, 0,
                                           lead_current_vector,
                                           junction_voltage_vector,
                                           lead_current_dqdt_vector, 0, 0, 0, 0));
  }
}

//-----------------------------------------------------------------------------
// Function      : Fourier::getLastPeriod_()
// Purpose       : finds the indices to access the last period of simulation
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void Fourier::getLastPeriod_()
{
  // We want to do the analysis on only the last period of the transient window. So here we find the indices
  // to access the endpoints of that interval.
  period_ = 1.0/at_;
  int numPoints = time_.size();
  int prdEnd = numPoints-1;
  double endTime = time_[prdEnd];

  // Use this formulation to compute lastPrdStart to avoid issues with 32-bit arithmetic.
  lastPrdStart_ = (at_*endTime - 1.0)/at_;

  if ( Teuchos::ScalarTraits<double>::magnitude( lastPrdStart_ ) <
         Teuchos::ScalarTraits<double>::eps() )
  {
    lastPrdStart_ = 0.0;
    prdStart_ = 0;
  }
  else if (lastPrdStart_ > 0)
  {
    // Initialize prdStart_ to be the index of the last element in time_.
    // Then scan until time_[i] <= endTime - period_.
    prdStart_ = prdEnd;
    while (time_[prdStart_] > lastPrdStart_)
    {
      prdStart_--;
    }
  }
  else
  {
    std::string msg = "Error in measure \"" + name_ + "\": The period is greater than the length of the time simulation. Exiting.";
    Report::UserFatal() << msg;
  }

}

//-----------------------------------------------------------------------------
// Function      : Fourier::updateDC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void Fourier::updateDC(
  Parallel::Machine comm,
  const std::vector<Analysis::SweepParam> & dcParamsVec,
  const Linear::Vector *solnVec,
  const Linear::Vector *stateVec,
  const Linear::Vector *storeVec,
  const Linear::Vector *lead_current_vector,
  const Linear::Vector *junction_voltage_vector,
  const Linear::Vector *lead_current_dqdt_vector)
{
}

//-----------------------------------------------------------------------------
// Function      : Fourier::updateAC()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical Models & Simulation
// Creation Date : 8/7/2019
//-----------------------------------------------------------------------------
void Fourier::updateAC(
  Parallel::Machine comm,
  const double frequency,
  const double fStart,
  const double fStop,
  const Linear::Vector *solnVec,
  const Linear::Vector *imaginaryVec,
  const Util::Op::RFparamsData *RFparams)
{

}

//-----------------------------------------------------------------------------
// Function      : Fourier::interpolateData_()
// Purpose       : evaluates interpolating polynomial at equidistant time pts
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
bool Fourier::interpolateData_()
{
  double A, B, C;
  int nData, j;

  // nData is the number of data points in the time_ and outVarValues_ vectors over the last period.
  int numPoints = time_.size();
  nData = numPoints-prdStart_;
  std::vector<double> h(nData-1, 0.0);
  std::vector<double> b(nData-1, 0.0);
  std::vector<double> u(nData-1, 0.0);
  std::vector<double> v(nData-1, 0.0);
  std::vector<double> z(nData, 0.0);

  // Cubic spline interpolation. We first need to find the z's.
  for (int i = 0; i < nData-1; i++)
  {
    h[i] = time_[i+1+prdStart_]-time_[i+prdStart_];
    b[i] = (6/h[i])*(outVarValues_[i+1+prdStart_]-outVarValues_[i+prdStart_]);
  }

  u[1] = 2*(h[0]+h[1]);
  v[1] = b[1]-b[0];

  for (int i=2; i < nData-1; i++)
  {
    u[i] = 2*(h[i]+h[i-1])-((h[i-1])*(h[i-1]))/u[i-1];
    v[i] = b[i]-b[i-1]-(h[i-1]*v[i-1])/u[i-1];
  }

  z[nData-1] = 0;
  for (int i=nData-2; i > 0; i--)
  {
    z[i] = (v[i]-h[i]*z[i+1])/u[i];
  }
  z[0] = 0;

  // Compute new, equally spaced time points.
  newTime_.resize(gridSize_,0.0);
  newValues_.resize(gridSize_,0.0);
  double step = period_/gridSize_;

  newTime_[0] = lastPrdStart_;
  for (int i = 1; i < gridSize_; i++)
  {
    newTime_[i] = newTime_[i-1] + step;
  }

  // Calculate the new values at the new time points.
  for (int i = 0; i < gridSize_; i++)
  {
    j = nData-1;
    while ((newTime_[i]-time_[j+prdStart_]) < 0)
    {
      j--;
    }
    A = (z[j+1]-z[j])/(6*h[j]);
    B = z[j]/2;
    C = -(h[j]/6)*z[j+1]-(h[j]/3)*z[j]+(outVarValues_[j+1+prdStart_]-outVarValues_[j+prdStart_])/h[j];
    newValues_[i] = outVarValues_[j+prdStart_] + (newTime_[i]-time_[j+prdStart_])*(C+(newTime_[i]-time_[j+prdStart_])*(B+(newTime_[i]-time_[j+prdStart_])*A));
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Fourier::calculateFT_()
// Purpose       : performs fourier analysis
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
void Fourier::calculateFT_()
{
  double tmp;
  mag_.resize(numFreq_,0.0);
  phase_.resize(numFreq_,0.0);
  nmag_.resize(numFreq_,0.0);
  nphase_.resize(numFreq_,0.0);
  freq_.resize(numFreq_,0.0);

  for (int i=0; i < gridSize_; i++)
  {
    for (int j=0; j < numFreq_; j++)
    {
      mag_[j] += (newValues_[i]*sin(j*2.0*M_PI*i/((double) gridSize_)));
      phase_[j] += (newValues_[i]*cos(j*2.0*M_PI*i/((double) gridSize_)));
    }
  }

  double convRadDeg = 180.0/M_PI;

  mag_[0] = phase_[0]/gridSize_;
  phase_[0] = 0.0;
  thd_ = 0.0;
  for(int i = 1; i < numFreq_ ; i++)
  {
    tmp = mag_[i]*2.0 /gridSize_;
    phase_[i] *= 2.0/gridSize_;
    freq_[i] = i * at_;
    mag_[i] = sqrt(tmp*tmp+phase_[i]*phase_[i]);
    phase_[i] = atan2(phase_[i],tmp)*convRadDeg;
    nmag_[i] = mag_[i]/mag_[1];
    nphase_[i] = phase_[i]-phase_[1];
    if(i>1) thd_ += nmag_[i]*nmag_[i];
  }
  thd_ = 100*sqrt(thd_);

}

//-----------------------------------------------------------------------------
// Function      : Fourier::getMeasureResult()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
double Fourier::getMeasureResult()
{
  // Only output results if transient data is available to analyze.
  //if( initialized_ && !time_.empty() )
  // change was made to exclude case where (for example) TD = end of simulation time,
  // and hence the time vector only has one point.  The case will now produce a
  // failed measure
  if ( initialized_ && time_.size() > 1 )
  {
    getLastPeriod_();

    interpolateData_();

    calculateFT_();

    calculated_ = true;
  }
  // Total harmonic distortion will be the calculation result.
  // printMeasureResult will be overloaded to print out all the harmonic information
  return calculationResult_ = thd_;
}

//-----------------------------------------------------------------------------
// Function      : Fourier::printMeasureResult( std::ostream& os )
// Purpose       : used to print the measurement result to an output stream
//                 object, which is typically the mt0, ma0 or ms0 file
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
std::ostream& Fourier::printMeasureResult( std::ostream& os)
{

  // Only output results if transient data is available to analyze.
  //if ( !time_.empty() )
  // change was made to exclude case where (for example) TD = end of simulation time,
  // and hence the time vector only has one point.  The case will now produce a
  // failed measure
  if ( initialized_ && time_.size() > 1 )
  {
    // Compute measure.
    if (!calculated_)
    {
      this->getMeasureResult();
    }

    // The next statement allows maintains a reasonable column width even if
    // the user to enter a precision_ value of 0.  All floating point numbers
    // are output in scientific notation based on the precision_ variable.  
    // numFreq and ident are integers and are output as such.
    int colWidth = (precision_ < 5) ? 15 : 10 + precision_; 
    int ident = 10;
    os << name_ << ":  No. Harmonics: " << numFreq_ << ", THD: " 
       << std::scientific << std::setprecision(precision_) << thd_ 
       << ", Gridsize: " << gridSize_
       << ", Interpolation Type: Cubic Spline" << std::endl;

    os << std::setw(ident) << "Harmonic" << std::setw(colWidth) << "Frequency"
       << std::setw(colWidth) << "Magnitude" << std::setw(colWidth) << "Phase"
       << std::setw(colWidth) << "Norm. Mag" << std::setw(colWidth) << "Norm. Phase" << std::endl;
    for (int i = 0; i < numFreq_; i++)
    {
      os << std::setw(ident) << i 
	 << std::setw(colWidth) << freq_[i]
         << std::setw(colWidth) << mag_[i]
         << std::setw(colWidth) << phase_[i]
         << std::setw(colWidth) << nmag_[i]
         << std::setw(colWidth) << nphase_[i] << std::endl;
    }
  }
  else
  {
    // this covers a failed FOUR measure.
    os << name_ << ": FAILED" << std::endl;
  }
  
  return os;
}

//-----------------------------------------------------------------------------
// Function      : Fourier::printVerboseMeasureResult( std::ostream& os )
// Purpose       :
// Special Notes : This produces the same output, to stdout typically, as
//                 printMeasureResult() does.
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical Models & Simulation
// Creation Date : 6/21/13
//-----------------------------------------------------------------------------
std::ostream& Fourier::printVerboseMeasureResult( std::ostream& os)
{
  printMeasureResult(os);
  return os;
}

//-----------------------------------------------------------------------------
// Function      : Fourier::printMeasureWarnings
// Purpose       : prints error message related to invalid time windows, etc.
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2015
//-----------------------------------------------------------------------------
void Fourier::printMeasureWarnings(const double endSimTime, const double startSweepVal,
                                   const double endSweepVal)
{
  //no op
}

//-----------------------------------------------------------------------------
// Function      : Fourier::printMeasureWindow
// Purpose       : prints information related to measure window
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, Electrical and Microsystem Modeling
// Creation Date : 02/5/2015
//-----------------------------------------------------------------------------
std::ostream& Fourier::printMeasureWindow(std::ostream& os, const double endSimTime,
				          const double startSweepVal, const double endSweepVal)
{
  //no op
  return os;
}

} // namespace Measure
} // namespace IO
} // namespace Xyce

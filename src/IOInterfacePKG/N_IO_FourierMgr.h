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
// Purpose       : This file is a class to manage measure statements in a sim.
//
// Special Notes :
//
// Creator       : Heidi Thornquist, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 03/10/2009
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef  Xyce_N_IO_FourierMgr_H
#define Xyce_N_IO_FourierMgr_H

#include <list>
#include <string>

#include <N_IO_fwd.h>
#include <N_LAS_Vector.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Param.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : FourierMgr
// Purpose       : This is a manager class for handling .four and .fft statements
//                 in a simulation
// Special Notes :
// Creator       : Heidi Thornquist, Electrical and Microsystems Modeling
// Creation Date : 3/10/2009
//-----------------------------------------------------------------------------
class FourierMgr
{
public:
  FourierMgr(const std::string &netlist_filename);

  // Destructor
  ~FourierMgr();

  // Return true if Fourier analysis is being performed on any variables.
  bool isFourierActive() const { return (!freqNumOutputVarsMap_.empty() && !time_.empty()); }

  // add .four or .fft line from netlist to list of things to perform analysis on.
  bool addFourierAnalysis( const Util::OptionBlock & fourierLine );

  // get the sensitivity parameters, if necessary
  bool getSensVars (const Util::OptionBlock &option_block);
  bool registerSensOptions(const Util::OptionBlock & option_block);

  void fixupSensFourierParameters(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager);
  void fixupFourierParameters(Parallel::Machine comm, const Util::Op::BuilderManager &op_builder_manager);

  // Called during the simulation to update the fourier objects held by this class
  void updateFourierData(Parallel::Machine comm, const double circuitTime, const Linear::Vector *solnVec, 
                         const Linear::Vector *stateVec, const Linear::Vector * storeVec, 
                         const Linear::Vector *lead_current_vector, 
                         const Linear::Vector *junction_voltage_vector, const Linear::Vector *lead_current_dqdt_vector,
                         const std::vector<double> &         objectiveVec,
                         const std::vector<double> &         dOdpVec,
                         const std::vector<double> &         dOdpAdjVec,
                         const std::vector<double> &         scaled_dOdpVec,
                         const std::vector<double> &         scaled_dOdpAdjVec);

  void outputResults( std::ostream& outputStream );
  
  //added to help register lead currents with device manager
  std::set<std::string> getDevicesNeedingLeadCurrents() { return devicesNeedingLeadCurrents_; }

private:
  void getLastPeriod_();

  bool interpolateData_();

  void calculateFT_();

  std::ostream& printResult_( std::ostream& os );

private:
  std::string           netlistFilename_;

  // these are needed for sensitivities:
  int                   sensitivityOptions_;
  bool sensitivityRequested;
  Util::ParamList       sensitivityVariableList_;

  // Store frequencies and solution variables
  int numFreq_, gridSize_;
  bool calculated_;
  std::vector<int> outputVarsPtr_;
  std::vector<double> time_;
  std::vector<double> freqVector_;  
  std::vector<double> outputVarsValues_;
  std::vector<std::string> names_;
  Util::ParamList depSolVarIterVector_;
  Util::Op::OpList outputVars_;
  std::vector<int> prdStart_; 
  std::vector<double> lastPrdStart_;
  std::vector<double> newTime_, newValues_, mag_, phase_, nmag_, nphase_, freq_, thd_;

  // added to support multiple .FOUR lines, some of which have the same fundamental frequency
  // and some of which have other fundamental frequencies.
  std::map<double,int> freqNumOutputVarsMap_;
  std::multimap<double,Util::Param> freqSolVarMap_;
  std::multimap<double,std::string> freqNamesMap_;

  // At parsing, this SENS version of freqNumOutputVarsMap_ is mainly 
  // used to avoid redundant .FOUR SENS statement.
  // ie, if .FOUR freq SENS is specified multiple times, but "freq" is the same, then
  // it should only be evaluated 1x. (as it will apply to all sensitivities in any case)
  std::map<double,int> sensFreqNumOutputVarsMap_;

  //added to help register lead currents with device manager
  std::set<std::string> devicesNeedingLeadCurrents_; 
};

bool registerPkgOptionsMgr(FourierMgr &fourier_manager, PkgOptionsMgr &options_manager);

} // namespace IO
} // namespace Xyce

#endif  // Xyce_N_IO_FourierMgr_H

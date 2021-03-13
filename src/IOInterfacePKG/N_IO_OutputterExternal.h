// -*-c++-*-
//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : Provide output capability via call-back mechanism
//                  into external program
//
// Special Notes  :
//
// Creator        : Tom Russo
//
// Creation Date  : 20 Feb 18
//
//-----------------------------------------------------------------------------
#ifndef Xyce_N_IO_OutputterExternal_h
#define Xyce_N_IO_OutputterExternal_h

#include <N_IO_Outputter.h>
#include <N_IO_ExtOutWrapper.h>

namespace Xyce {
namespace IO {

namespace Outputter {

class OutputterExternal : public Interface
{
public:
  OutputterExternal(Parallel::Machine comm, OutputMgr &output_manager,
                 ExternalOutputWrapper * outputWrapper);

  virtual ~OutputterExternal();

private:
  // Disable copy constructor and assignment operator
  OutputterExternal(const OutputterExternal &);
  OutputterExternal &operator=(const OutputterExternal &);

public:
  virtual void doSetAnalysisMode(Analysis::Mode analysis_mode){};

  virtual void doOutputTime(
     Parallel::Machine           comm,
     const Linear::Vector &        solution_vector,
     const Linear::Vector &        state_vector,
     const Linear::Vector &        store_vector, 
     const Linear::Vector &        lead_current_vector,
     const Linear::Vector &        junction_voltage_vector) ;

  virtual void doOutputFrequency(
      Parallel::Machine comm, 
      double frequency,
      double fStart,
      double fStop, 
      const Linear::Vector &realSolutionVector, 
      const Linear::Vector &imaginarySolutionVector,
      const Util::Op::RFparamsData &RFparams);

  virtual void doOutputHB_FD(
     Parallel::Machine             comm,
     const std::vector<double> &   freqPoints,
     const Linear::BlockVector &     freqDomainSolutionVecReal,
     const Linear::BlockVector &     freqDomainSolutionVecImaginary,
     const Linear::BlockVector &     freqDomainLeadCurrentVecReal,
     const Linear::BlockVector &     freqDomainLeadCurrentVecImaginary,
     const Linear::BlockVector &     freqDomainJunctionVoltageVecReal,
     const Linear::BlockVector &     freqDomainJunctionVoltageVecImaginary);

   virtual void doOutputHB_TD(
     Parallel::Machine             comm,
     const std::vector<double> &   timePoints,
     const Linear::BlockVector &     timeDomainSolutionVec,
     const Linear::BlockVector &     timeDomainLeadCurrentVec,
     const Linear::BlockVector &     timeDomainJunctionVoltageVec);

  // Extend to support other output types by adding non-null implementations
  // of the functions in the Interface class.

  virtual void doFinishOutput();

  virtual void doStartStep(int current_step, int number_of_steps);

  virtual void doResetIndex() {index_=0;};

  virtual void doSteppingComplete();

private:
  OutputMgr &              outputManager_;
  ExternalOutputWrapper *  theOutputWrapper_;
  int                      currentStep_;
  int                      numberOfSteps_;

  Util::Op::OpList         opList_;
  std::vector<std::string> fieldNames_;
  int                      index_;
  bool                     initialized_;
};

}
}
}
#endif

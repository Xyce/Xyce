//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : External code output wrapper class
//
// Special Notes  : The purpose of this class is to sit between
//                  instances of the ExternalOutputInterface derived
//                  classes, and provide a set of higher level functions
//
// Creator        : Tom Russo, SNL, Electrical Models and Simulation
//
// Creation Date  : 2/8/2018
//
//-------------------------------------------------------------------------
#ifndef N_IO_ExtOutWrapper_H
#define N_IO_ExtOutWrapper_H

#include <vector>
#include <string>

#include <N_UTL_fwd.h>
#include <N_IO_ExtOutInterface.h>
#include <N_UTL_Param.h>
#include <N_IO_OutputTypes.h>

namespace Xyce {
namespace IO {

/// External code coupling output interface wrapper class
///
/// This concrete class is used by Xyce to call methods of external output
/// objects (derived from the ExternalOutputInterface) provided by coupled
/// codes.
///
/// This class exists so that the actual ExternalOutputInterface base class
/// can be as simple as possible and depend only on a few Xyce forward 
/// declaration header files.  This keeps the number of internal Xyce header
/// files that must be installed by "make install" to a minimum, and minimizes
/// the need for users of this interface to understand much about Xyce
/// internals.
///
/// Each object of this class "owns" a pointer to an object of the
/// ExternalOutputInterface type, and only calls the public interface of that
/// object.
class ExternalOutputWrapper
{

public:

  ExternalOutputWrapper(ExternalOutputInterface * extIntPtr);

  /// Destructor simply forgets the pointer --- it is the caller's responsibility to delete the object
  ~ExternalOutputWrapper() { extIntPtr_ = 0;};

  /// Return a param list corresponding to the requested output variables
  inline Util::ParamList & getParamList()
  {
    return paramList_;
  };

  /// Inform the external code of the actual names we're outputting
  ///
  /// This is basically the equivalent of outputting a header
  inline void outputFieldNames(std::vector<std::string> & outputNames)
  {
    extIntPtr_->outputFieldNames(outputNames);
  };


  /// Send vector of reals to the external code
  inline void outputReal(std::vector<double> & outputData)
  {
    extIntPtr_->outputReal(outputData);
  };


  /// Send vector of reals to the external code
  inline void outputComplex(std::vector<complex> & outputData)
  {
    extIntPtr_->outputComplex(outputData);
  };

  /// Inform external code that Xyce is done outputting simulation data
  inline void finishOutput()
  {
    extIntPtr_->finishOutput();
  };


  /// Inform external code that Xyce is done running a .step simulation
  inline void finishedStepping()
  {
    extIntPtr_->finishedStepping();
  };

  /// Inform external code that Xyce is now running a new step in a .step
  /// simulation
  inline void newStepOutput(int sN, int maxS)
  {
    extIntPtr_->newStepOutput(sN,maxS);
  };

  /// Get the unique name for this output connection so we can use it in
  /// error messages.
  inline std::string getName()
  {
    return extIntPtr_->getName();
  };

  /// Get the output type for this output connection 
  inline OutputType::OutputType getOutputType()
  {
    return extIntPtr_->getOutputType();
  };

private:
  void checkVars_();
  void normalizeVarList_();

  ExternalOutputInterface * extIntPtr_;
  Util::ParamList paramList_;
};

}
}
#endif

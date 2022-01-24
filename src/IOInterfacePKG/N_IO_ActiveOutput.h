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

#ifndef Xyce_N_IO_ActiveOutput_h
#define Xyce_N_IO_ActiveOutput_h

#include <N_IO_fwd.h>
#include <N_IO_PrintTypes.h>
#include <N_ANP_fwd.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : Xyce::IO::ActiveOutput
// Purpose       :
// Special Notes :
// Creator       : Dave Baur
// Creation Date : 24 Jul 2014
//-----------------------------------------------------------------------------
/// @brief provide a mechanism for analyses to change the set of active outputs
///
/// The output manager maintains a stack of "active" outputters.  This
/// stack is a vector of vectors of outputter pointers.  The vector at
/// the top of the stack (or rather "back" of the vector-of-vectors)
/// is the vector of currently selected ("active") outputters.
/// 
/// When constructed, objects of this class cause the
/// output manager to create an empty vector of
/// outputter pointers, and subsequent operations in the 
/// output manager to add or remove outputters will add pointers
/// to this vector at the "top" of the stack (or rather,
/// the "back" of the vector-of-vectors-of-pointers).
///
/// When destructed, objects of this class cause the
/// output manager to pop the last element off of the 
/// stack, thereby restoring the list of active outputters
/// to what it was before we created ourself.
/// 
/// This is used in the analysis package by the various
/// analysis types, which will create an object of this
/// class and call the "add" method with print type and analysis
/// mode.  This will cause the output manager to select
/// all appropriate ".print" lines and the outputters
/// associated with them, and add them to the currently
/// selected vector of "active" outputters.
/// 
/// The analysis code will then do its work and its output
/// and let the output manager decide where to send it.
///
/// When the analysis is complete, it lets the object of
/// ActiveOutput type go out of scope, thereby returning the
/// output manager's list of actives to whatever it was
/// before the analysis began.
/// @author Dave Baur
/// @date   7/24/2014
class ActiveOutput
{
public:
  ActiveOutput(OutputMgr &output_manager);

  ~ActiveOutput();

  void add(PrintType::PrintType print_type, Analysis::Mode analysis_mode);

  void add(Parallel::Machine comm, Analysis::Mode analysis_mode);

  void resetIndex();

  void reopenTmpFile();

  void copyTmpFileToOutputFile();

  void setStepSweepVector(const Analysis::SweepVector & step_sweep_parameters);

  void setDCSweepVector(const Analysis::SweepVector & dc_sweep_parameters);

private:
  OutputMgr &     outputManager_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_ActiveOutput_h

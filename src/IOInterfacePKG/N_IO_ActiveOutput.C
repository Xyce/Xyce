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
// Purpose        : Provide a means for analyses to select a current set of 
//                  "active" outputters
//
// Special Notes  : 
//
// Creator        : Dave Baur
//
// Creation Date  : 7/24/2014
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_ActiveOutput.h>
#include <N_IO_OutputMgr.h>
#include <N_ANP_SweepParam.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::ActiveOutput
// Purpose       : constructor
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief Select a new active output set
/// 
/// @param[in] output_manager   Reference to the output manager
///
/// At construction, instruct the output manager to push an empty list
/// of active outputters onto its stack.  Later operations will add outputters
/// to the list, selected from a map of outputters that were previously
/// allocated from the set of .print lines in the netlist.
///
/// Users of this class are typically analysis objects, which create
/// an object of this class, select outputters appropriate to the 
/// analysis type (via the "add" method), and keep it around until the 
/// analysis is complete, at which time they let it go out of scope to restore
/// the set of active outputters to what it was previously.
///
/// @author Dave Baur
/// @date 7/24/2014
ActiveOutput::ActiveOutput(OutputMgr &output_manager)
  : outputManager_(output_manager)
{
  outputManager_.pushActiveOutputters();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::~ActiveOutput
// Purpose       : destructor
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief destruct an ActiveOutput object
/// 
/// Destroying an ActiveOutput object causes the output manager to 
/// pop the currently active list of outputters off of its stack.
///
/// @author Dave Baur
/// @date 7/24/2014
ActiveOutput::~ActiveOutput()
{
  outputManager_.popActiveOutputters();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::add
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief Select active outputters of given type for given analysis mode
///
/// @param[in] print_type Output type to select
/// @param[in] analysis_mode Current analysis mode
///
///
/// When called, causes the output manager to add to its current vector
/// of active outputters all those outputters associated with .print
/// lines of the selected print_type.  
/// 
/// @author Dave Baur
/// @date 7/24/2014
void
ActiveOutput::add(PrintType::PrintType print_type, Analysis::Mode analysis_mode)
{
  outputManager_.addActiveOutputter(print_type, analysis_mode);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::resetIndex
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief Reset "Index" counter in all active outputters
///
/// Causes the output manager to reset to zero any "index" counters maintained
/// by all active outputters.  
/// 
/// @author Dave Baur
/// @date 7/24/2014
void
ActiveOutput::resetIndex()
{
  outputManager_.resetIndex();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::add
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief Trigger output initialization by output manager
///
/// @param[in] comm communicator
/// @param[in] analysis_mode Current analysis mode
///
///
/// When called, causes the output manager to run its "prepareOutput" method,
/// which allocates outputters associated with analysis_mode and adds them
/// to the current active outputter set.
/// 
/// @author Dave Baur
/// @date 7/24/2014
void
ActiveOutput::add(Parallel::Machine comm, Analysis::Mode analysis_mode)
{
  outputManager_.prepareOutput(comm, analysis_mode);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::setStepSweepVector
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief Interface to output manager setStepSweepVector method
///
/// @param[in] sweep_vector vector of sweep values
///
/// This method provides an interface to the output manager's 
/// setStepSweepVector method.
/// 
/// @author Dave Baur
/// @date 7/24/2014
void
ActiveOutput::setStepSweepVector(
  const Analysis::SweepVector & sweep_vector)
{
  outputManager_.setStepSweepVector(sweep_vector);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::setDCSweepVector
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Dave Baur
// Creation Date : 7/24/2014
//-----------------------------------------------------------------------------
///
/// @brief Interface to output manager setDCSweepVector method
///
/// @param[in] sweep_vector vector of sweep values
///
/// This method provides an interface to the output manager's 
/// setDCSweepVector method.
/// 
/// @author Dave Baur
/// @date 7/24/2014
void
ActiveOutput::setDCSweepVector(
  const Analysis::SweepVector & sweep_vector)
{
  outputManager_.setDCSweepVector(sweep_vector);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::reopenTmpFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
///
/// @brief re-open any tmp files in all active outputters
///
/// Causes the output manager to close and re-open any tmp files used
/// by any of the active outputters.  This is currently used by the
/// outputters for HB_IC information.
///
/// @author Pete Sholander
/// @date 11/19/2019
void
ActiveOutput::reopenTmpFile()
{
  outputManager_.reopenTmpFile();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ActiveOutput::copyTmpFileToOutputFile
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 11/19/2019
//-----------------------------------------------------------------------------
///
/// @brief copy any tmp files to the actual output files in all active outputters
///
/// Causes the output manager to copy any tmp files to the actual output files,
/// and then delete those tmp files.  This is currently used by the outputters
/// for HB_IC information.
///
/// @author Pete Sholander
/// @date 11/19/2019
void
ActiveOutput::copyTmpFileToOutputFile()
{
  outputManager_.copyTmpFileToOutputFile();
}

} // namespace IO
} // namespace Xyce

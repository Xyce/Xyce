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

#ifndef Xyce_N_UTL_PrintStats_h
#define Xyce_N_UTL_PrintStats_h

#include <iosfwd>

#include <N_UTL_Stats.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace Stats {

//-----------------------------------------------------------------------------
// Function      : printXML
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:21:57 2014
//-----------------------------------------------------------------------------
///
/// Write the metric to the os stream in XML format
///
/// @param os           output stream
/// @param metrics_mask make of metrics to write
/// @param checkpoint   true if these values are checkpointed
///
/// @return output stream
///
///
std::ostream &printXML(std::ostream& os, MetricsMask metrics_mask, bool checkpoint);

//-----------------------------------------------------------------------------
// Function      : printXML
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:21:57 2014
//-----------------------------------------------------------------------------
///
/// Write the metric to the os stream as a table, single processor.
///
/// @param os           output stream
/// @param root_timer   root timer
/// @param metrics_mask make of metrics to write
/// @param checkpoint   true if these values are checkpointed
///
/// @return output stream
///
std::ostream &printStatsTable(std::ostream& os, Stat root_timer, MetricsMask metrics_mask, bool timer_checkpoint);

//-----------------------------------------------------------------------------
// Function      : printXML
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Sep  3 10:21:57 2014
//-----------------------------------------------------------------------------
///
/// Write the metric to the os stream as a table, accumulated and reported on processor rank 0.
///
/// @param os                   output stream
/// @param root_timer           root timer
/// @param metrics_mask         make of metrics to write
/// @param checkpoint           true if these values are checkpointed
/// @param parallel_machine     communicator to accumulate stats
///
/// @return output stream
///
std::ostream &printStatsTable(std::ostream& os, Stat root_timer, MetricsMask metrics_mask, bool timer_checkpoint, Parallel::Machine parallel_machine);

} // namespace Stats
} // namespace Xyce

#endif // Xyce_N_UTL_PrintStats_h

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
// Purpose       : Handles the parameter data associated with sweeps.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 7/28/2020
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_SweepParamFreeFunctions_h
#define Xyce_N_ANP_SweepParamFreeFunctions_h

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include <random>

#include <N_ANP_fwd.h>
#include <N_LOA_fwd.h>
#include <N_PDS_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Param.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : parseSweepParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Wed Oct  1 09:19:28 2014
//-----------------------------------------------------------------------------
///
/// Populate the sweep params from the parameter list.
///
/// @invariant
///
/// @param sweep_param  sweep parameters objec t populate
/// @param first        begin iterator of the parameter list
/// @param last         end iterator of the parameter list
///
SweepParam parseSweepParams(Xyce::Util::ParamList::const_iterator first, Xyce::Util::ParamList::const_iterator last);

bool updateSweepParams(int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);
bool updateSweepParams(Loader::Loader &loader, int step_count, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end, bool overrideOriginal);
int setupSweepLoop(Parallel::Machine comm, Loader::Loader &loader, std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);
int setSweepLoopVals(std::vector<SweepParam>::iterator begin, std::vector<SweepParam>::iterator end);

bool processDataStatements(
    const Xyce::Util::OptionBlock & paramsBlock,
    std::map< std::string, std::vector<std::string> > & dataNamesMap,
    std::map< std::string, std::vector< std::vector<double> > > & dataTablesMap);

bool convertData(
    SweepVector & stepSweepVector,
    const std::map< std::string, std::vector<std::string> > & dataNamesMap,
    const std::map< std::string, std::vector< std::vector<double> > > & dataTablesMap
    );

bool isDataSpecified(const Xyce::Util::OptionBlock & paramsBlock);

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_SweepParamFreeFunctions_h

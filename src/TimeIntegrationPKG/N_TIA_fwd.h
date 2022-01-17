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
// Purpose        : Foward declarations
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 01/11
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_fwd_h
#define Xyce_N_TIA_fwd_h

namespace Xyce {
namespace TimeIntg {

class DataStore;
class TimeIntegrationMethod;
class StepErrorControl;
class MPDEInterface;
class OneStep;
class TIAParams;
class TwoLevelError;
class WorkingIntegrationMethod;

// For historical reasons, onestep is indexed to 7, and gear is indexed to 8.  
// All the other choices are for methods that no longer exist in the source code.
// For now they must be kept at 7 and 8 b/c there are still test cases that have 
// statements like ".options timeint method=7".
enum methodsEnum {
  NO_TIME_INTEGRATION, OBSOLETE1, OBSOLETE2, OBSOLETE3, OBSOLETE4, OBSOLETE5, OBSOLETE6,
  ONESTEP, GEAR};

} // namespace TimeIntg
} // namespace Xyce

#endif // Xyce_N_TIA_fwd_h

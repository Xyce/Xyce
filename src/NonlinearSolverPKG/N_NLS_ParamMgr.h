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

//-------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter
//
// Creation Date  : 10/25/02
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_ParamMgr_h
#define Xyce_N_NLS_ParamMgr_h

#include <vector>

#include <N_IO_fwd.h>

#include <N_NLS_NLParams.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : ParamMgr
//
// Purpose       : It is not unusual for a Xyce simulation to have several
//                 different parameter sets.  Each parameter set (for dcop,
//                 transient, etc.) is stored in a NLParams data
//                 structure.
//
//                 The management of which parameter set is currently being
//                 used is handled by this class.
//
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
class ParamMgr
{
public:
  ParamMgr(const IO::CmdParse & cp);
  ~ParamMgr();

  bool addParameterSet (AnalysisMode mode, NLParams & nlp);

  bool getParams (AnalysisMode mode, NLParams & nlp); 

  bool getCurrentParams (NLParams & nlp); 
  inline bool setAnalysisMode (AnalysisMode mode);
  inline void resetAnalysisMode();
  inline AnalysisMode getAnalysisMode() const;

protected:
private:

public:

protected:

private:
  std::vector<NLParams> paramVector_;
  AnalysisMode   currentMode_;

  bool modeToggled_;
  bool gcp_calledBefore_;
  bool paramsChanged_;
};

//-----------------------------------------------------------------------------
// Function      : ParamMgr::setAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
inline bool ParamMgr::setAnalysisMode (AnalysisMode mode)
{
  if (currentMode_ != mode) 
  {
    paramsChanged_ = true;
    currentMode_   = mode;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParamMgr::resetAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
inline void ParamMgr::resetAnalysisMode()
{
  if (currentMode_ != DC_OP) 
  {
    paramsChanged_ = true;
    currentMode_   = DC_OP;
  }
}

//-----------------------------------------------------------------------------
// Function      : ParamMgr::getAnalysisMode
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
inline AnalysisMode ParamMgr::getAnalysisMode() const
{
  return currentMode_;
}

} // namespace Nonlinear
} // namespace Xyce

#endif


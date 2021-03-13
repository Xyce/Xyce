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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/26/02
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_NLS_NLParams.h>
#include <N_NLS_ParamMgr.h>
#include <N_UTL_FeatureTest.h>

// ---------- Forward Declarations ----------

// ---------- Static Initializations ----------

namespace Xyce {
namespace Nonlinear  {

//-----------------------------------------------------------------------------
// Function      : ParamMgr::ParamMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
ParamMgr::ParamMgr (const IO::CmdParse & cp)
  : currentMode_(DC_OP),
    modeToggled_(true),
    gcp_calledBefore_(false),
    paramsChanged_(false)
{
  paramVector_.resize(NUM_MODES, NLParams(DC_OP, cp));
  paramVector_[TRANSIENT] = NLParams(TRANSIENT, cp);

  currentMode_   = DC_OP;
}

//-----------------------------------------------------------------------------
// Function      : ParamMgr::~ParamMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
ParamMgr::~ParamMgr ()
{

}

//-----------------------------------------------------------------------------
// Function      : ParamMgr::addParameterSet
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
bool ParamMgr::addParameterSet (AnalysisMode mode, NLParams & nlp)
{
  // add a check to make see if this parameter name has already been used.
  paramVector_[mode] = nlp;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : ParamMgr::getParams
// Purpose       : 
// Special Notes : might not need this one...
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
bool ParamMgr::getParams (AnalysisMode mode, NLParams & nlp)
{
  nlp = paramVector_[mode];
  return true;
} 

//-----------------------------------------------------------------------------
// Function      : ParamMgr::getCurrentparams
// Purpose       : copies the official set of current parameters into the
//                 passed reference.
//
// Special Notes : 
//
// If the toggle flag isn't set, then don't change the parameters.  
// Leave them as-is.
//
// The modeToggled flag is connected to the "resetTranNLS" setting
// from time integration.  If true, it means "use different params"
// for transient than for dcop.  In general, if resetTranNLS=0 in
// the netlist, then the code  will only use the first (DCOP) set of
// parameters and will never use any others.  
//
// It is kind of an obsolete option. When it was first put in the code,
// it was not possible to specify a separate set of options for transient.
// You could only specify one set, and that was it.  Initially, to get
// around this limitation, a set of options for transient was hardwired
// into the code.  This hardwired set became the default for transient.
// If you didn't want to use this default, you had to invoke the 
// resetTranNLS option, and the code would just use the DCOP options.
//
// Now, however, you can easily specify a separate set for transient,
// so this feature probably isn't necessary. ERK.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/02
//-----------------------------------------------------------------------------
bool ParamMgr::getCurrentParams (NLParams & nlp) 
{
  // copy over a new copy of the params if one of the following is true:
  //   1) parameters have changed, and resetTranNLS is true.
  //   2) this function has never been called before.
  bool updateFlag = ((paramsChanged_ && modeToggled_) || !gcp_calledBefore_);

  if ( updateFlag )
  { 
    nlp = paramVector_[currentMode_];
    paramsChanged_ = false;
    gcp_calledBefore_ = true;
  } 

  // if this has never been called before, or if the parameters have
  // been changed recently, then print them out.
  if (VERBOSE_NONLINEAR && updateFlag )
    nlp.printParams (lout());

  return true;
}

} // namespace Nonlinear
} // namespace Xyce

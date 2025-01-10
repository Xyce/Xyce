//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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

#include <Xyce_config.h>

#include <vector>
#include <string>

#include <N_IO_fwd.h>
#include <N_IO_OutputTypes.h>
#include <N_UTL_fwd.h>
#include <N_ERH_ErrorMgr.h>

#include <N_UTL_OutputAPIHelpers.h>
#include <N_IO_ExtOutInterface.h>
#include <N_IO_ExtOutWrapper.h>

namespace Xyce {
namespace IO {
  
//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ExternalOutputWrapper::ExternalOutputWrapper
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 02/08/2018
//-----------------------------------------------------------------------------
/// Construct a wrapper for an external output interface object
///
/// @param extIntPtr  Pointer to user-defined output object
///
/// This class saves a copy of the pointer and uses it to call methods
/// the user has implemented.  Upon construction, we check the
/// output variables the user has requested and signal to the user's
/// code whether we could parse them all.
ExternalOutputWrapper::ExternalOutputWrapper ( ExternalOutputInterface * extIntPtr)
  : extIntPtr_(extIntPtr)
{
  checkVars_();
  normalizeVarList_();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ExternalOutputWrapper::checkVars_
// Purpose       : Try to parse requested output vars, report status
// Special Notes : Saves parsed parameters in our ParamList
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 02/08/2018
//-----------------------------------------------------------------------------
void ExternalOutputWrapper::checkVars_()
{
  std::vector<std::string> outputVars;
  std::vector<bool> parseStatuses;
  extIntPtr_->requestedOutputs(outputVars);

  paramList_.clear();
  bool successfulParse=Util::stringsToParamList(outputVars,
                                                      paramList_,
                                                      parseStatuses);
  extIntPtr_->reportParseStatus(parseStatuses);

  if (!successfulParse)
  {
    // do some sort of error reporting.  Use one of the Error guys,
    // not a Fatal.
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::IO::ExternalOutputWrapper::normalizeVarList_
// Purpose       : Make sure that independent variable for time or frequency
//                 domain outputs is present in param list, add it
//                 at the beginning if not.
// Special Notes : This is similar to what is done in the parsePRINTBlock
//                 method of the output manager, but simpler because we're
//                 not making all manner of assumptions about fallbacks and
//                 unrequested ancillary outputs
// Scope         : private
// Creator       : Tom Russo
// Creation Date : 02/22/2018
//-----------------------------------------------------------------------------
void ExternalOutputWrapper::normalizeVarList_()
{
    OutputType::OutputType theOutputType=extIntPtr_->getOutputType();
    std::string expectedIndVar;
    switch(theOutputType)
    {
    case OutputType::TRAN:
    case OutputType::AC_IC:
    case OutputType::HB_TD:
    case OutputType::HB_IC:
    case OutputType::HB_STARTUP:
      expectedIndVar="TIME";
      break;
    case OutputType::AC:
    case OutputType::HB_FD:
      expectedIndVar="FREQ";
      break;
    default:
      expectedIndVar="";
    }

    // Now check if any param in our list is that expected parameter
    // We do NOT require that it be first, because the user may have
    // special needs.  Just require that it be present.
    if (expectedIndVar.size()!=0)
    {
      Util::ParamList::iterator iterParam = paramList_.begin();
      bool foundIt=false;
      for ( ; iterParam != paramList_.end() && !foundIt ; ++iterParam)
      {
        if (iterParam->tag() == expectedIndVar)
        {
          foundIt = true;
        }
      }
      if (!foundIt)
        paramList_.push_front(Util::Param(expectedIndVar,0.0));
    }
}

}
}

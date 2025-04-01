//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Creator        : David Baur
//
// Creation Date  : 1/7/2014
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <list>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <N_ERH_Message.h>
#include <N_ERH_Messenger.h>

namespace Xyce {
namespace Report {

namespace {

//-----------------------------------------------------------------------------
// Function      : prefix
// Purpose       : Construct message parameters based on the report_mask and the
//                 message content.
// Special Notes :
// Scope         : Private
// Creator       : David Baur
// Creation Date : 1/7/2014
//-----------------------------------------------------------------------------
std::ostream &prefix(std::ostream &os, unsigned report_mask)
{

  if (report_mask & MSG_USER)
    os << "Netlist ";

  if (report_mask & MSG_DEVEL)
    os << "Application ";

  // This annotation is inappropriate for a normal user.
  // If there is value in reporting whether a message is "symmetric" it
  // is to developers only.  This should either be wrapped in debugging
  // conditionals, or in developer-specific verbosity conditionals
  //  if ((report_mask & MSG_SYMMETRIC) && MSG_PARALLEL != 0)
  //    os << "symmetric ";

  if ((report_mask & MSG_TYPE_MASK) == MSG_FATAL)
    os << "error";

  if ((report_mask & MSG_TYPE_MASK) == MSG_ERROR)
    os << "error";

  if ((report_mask & MSG_TYPE_MASK) == MSG_WARNING)
    os << "warning";

  if ((report_mask & MSG_TYPE_MASK) == MSG_INFORMATION)
    os << "info";

  if ((report_mask & MSG_TYPE_MASK) == MSG_DEBUG)
    os << "debug";

  return os;
}

}  // namespace <unnamed>

MessageCode
MessageCode::s_defaultMessageCode(100000000);

Message::Message(
  MessageType   message_type,
  MessageCode   message_code)
  : messageType_(message_type),
    messageCode_(message_code),
    oss_(),
    netlistLocation_(),
    functionName_(0)
{}

Message::~Message()
{
  std::ostringstream os;

  if ((messageType_ & MSG_TERMINATE) == 0)
  {
    prefix(os, messageType_);

    if (netlistLocation_.getLineNumber() > 0)
      os << " in file " << netlistLocation_.getFilename() << " at or near line " << netlistLocation_.getLineNumber() << "\n";
    else
      os << ": ";
  }

  if (functionName_)
    os << "function " << functionName_ << ":\n";

  os << oss_.str();

  Xyce::Report::report_message(os.str().c_str(), messageType_, messageCode_);
}

} // namespace Report
} // namespace Xyce

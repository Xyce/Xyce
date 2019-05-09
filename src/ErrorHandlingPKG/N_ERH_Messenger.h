//-------------------------------------------------------------------------
//   Copyright 2002-2019 National Technology & Engineering Solutions of
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
//
//
//
//-------------------------------------------------------------------------


#ifndef Xyce_N_ERH_Messenger_h
#define Xyce_N_ERH_Messenger_h

#include <iosfwd>
#include <cstddef>

#include <N_ERH_fwd.h>
#include <N_PDS_fwd.h>

namespace Xyce {
namespace Report {

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @file
 */

/**
 * @brief Member function <b>get_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
unsigned get_message_count(unsigned message_type);

/**
 * @brief Member function <b>reset_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 */
void reset_message_count(unsigned message_type);

/**
 * @brief Member function <b>set_max_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @param max_count		an <b>unsigned int</b> ...
 *
 */
void set_max_message_count(unsigned message_type, unsigned max_count);

/**
 * @brief Member function <b>get_max_message_count</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
unsigned get_max_message_count(unsigned message_type);

/**
 * @brief Member function <b>get_message_name</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @return			a <b>std::string</b> ...
 */
const std::string &get_message_name(unsigned message_type);

/**
 * @brief Member function <b>set_message_name</b> ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @param max_count		an <b>unsigned int</b> ...
 *
 * @param name			a <b>std::string</b> const ...
 *
 */
void register_message_type(unsigned message_type, unsigned max_count, const char *name);

/**
 * @brief Function <b>reset_message_group</b> sets the count to zero of all messages in the
 * specified throttle group.
 *
 * @param throttle_group	an <b>int</b> value of the throttle group to reset.
 *
 */
void reset_throttle_group(int throttle_group);

/**
 * @brief Member function <b>report_message</b> ...
 *
 * @param message		an <b>char</b> const pointer ...
 *
 * @param message_type		an <b>unsigned int</b> ...
 *
 * @param message_code		a <b>MessageCode</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
void report_message(const char *message, unsigned message_type, const MessageCode &message_code);

/**
 * @brief Function <b>add_deferred_message</b> adds a message to the deferred message queue.
 *
 * @param message_type		an <b>int</b> value of the message type, usually WARNING or DOOMED
 *
 * @param message_id		a <b>MessageId</b> value of the message identifier
 *
 * @param throttle_cutoff	a <b>size_t</b> value to display before the message is no longer
 *                              displayed. 
 *
 * @param throttle_group	an <b>int</b> value to identify the throttle group that this message
 *                              belongs to. 
 *
 * @param header		a <b>char</b> const pointer to the message header string.
 *
 * @param aggegrate		a <b>char</b> const pointer to the message aggregation string.
 *
 */
void add_deferred_message(unsigned message_type, MessageId message_id, size_t throttle_cutoff, int throttle_group, const char *header, const char *aggegrate);

/**
 * @brief Function <b>report_deferred_messages</b> aggregates and reports the message on the root
 * processor. 
 *
 * @param comm			a <b>ParallelMachine</b> communicator.
 *
 */
void report_deferred_messages(Parallel::Machine comm);

/**
 * @brief Function <b>aggregate_messages</b> writes a message message to the output string by
 * joining the messages from each processor, in order.  Each message is separated by the specified
 * separation string.
 *
 * @param comm			a <b>ParallelMachine</b> communicator.
 *
 * @param os			a <b>std::ostream</b> reference to the output stream to receive the
 *                              aggregated message.
 *
 * @param separator		a <b>char</b> const pointer to the separation string.
 *
 */
void aggregate_messages(Parallel::Machine comm, std::ostringstream &os, const char *separator = ", ");

///
/// @}
///

} // namespace Report
} // namespace Xyce

#endif // Xyce_N_ERH_Messenger_h

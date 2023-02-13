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


#ifndef Xyce_N_UTL_ReportHandler_h
#define Xyce_N_UTL_ReportHandler_h

#include <iosfwd>
#include <string>
#include <sstream>

namespace Xyce {

///
/// @addtogroup runtime_message_detail
/// @{
///

/**
 * @brief Type definition REH is a pointer to a function of type void that takes a const
 * std::exception reference as a parameter.
 *
 */
typedef void (*REH)(const char *message, unsigned type);

/**
 * @brief Function <b>default_report_handler</b> is the default
 * error reporter for Xyce exceptions.  
 *
 * @param message		a <b>char</b> const pointer to the message to be
 *				displayed.
 *
 * @param type			an <b>int</b> value of the type of message from the
 *				enumeration <b>type</b>
 *
 */
void default_report_handler(const char *message, unsigned type);

/**
 * @brief Function <b>set_report_handler</b> sets the exception report function to be called when an
 * report_exception() is called.
 *
 * @param reh		a <b>REH</b> of the new exception reporter.
 *
 * @return		a <b>REH</b> to the previous exception reporter.
 */
REH set_report_handler(REH reh);

/**
 * @brief Function <b>report</b> calls the current exception reporter to report the message in
 * <b>x</b>.
 *
 * @param message	a <b>char</b> const pointer to the message to report.
 *
 * @param type		an <b>int</b> value of the type of message.
 *
 */
void report(const char *message, unsigned type);

///
/// @}
///

} // namespace Xyce

#endif // Xyce_N_UTL_ReportHandler_h

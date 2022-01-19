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
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_ERH_fwd_H
#define Xyce_ERH_fwd_H

#include <cstddef>

namespace Xyce {
namespace Report {

typedef std::ptrdiff_t MessageId;

struct MessageCode;

} // namespace Report
} // namespace Xyce

// The do-while is necessary to prevent usage of this macro from changing
// program semantics (e.g. dangling-else problem). The obvious implementation:
// if (expr) ; else throw ...
// is not adequate because it causes ambiguous else statements in this context:
// if (something)
//   ThrowRequire(foo);
// The compiler does not know whether the else statement that the macro inserts
// applies to the "if (something) " or the "if (expr)".
#define ThrowRequireMsg(expr, message)          \
  do {                                          \
    if (!(expr)) {                              \
      throw std::runtime_error(message);        \
    }                                           \
  } while (false)

#define ThrowRequire(expr) ThrowRequireMsg(expr, #expr)

namespace N_ERH_ErrorMgr = Xyce::Report;

#endif // Xyce_ERH_fwd_H

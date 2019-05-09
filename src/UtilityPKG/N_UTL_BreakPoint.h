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

//-----------------------------------------------------------------------------
//
// Purpose        : Provide a basic "breakpoint" class that can be used in
//                  STL sets, with appropriate operators so that the set will
//                  automatically eliminate duplicates on insert, with
//                  "duplicate" meaning "within a certain tolerance"
//
// Special Notes  : 
//
// Creator        : Tom Russo, SNL, Component Information and Models
//
// Creation Date  : 04/28/2004
//
//-----------------------------------------------------------------------------
#ifndef Xyce_N_UTL_BreakPoint_h
#define Xyce_N_UTL_BreakPoint_h

#include <N_UTL_Math.h>

namespace Xyce {
namespace Util {

struct BreakPointLess;

//-----------------------------------------------------------------------------
// ClassFunction : BreakPoint
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
class BreakPoint
{
  friend struct BreakPointLess;                 ///< Allow less compare to access values
  friend struct BreakPointEqual;                ///< Allow equal compare to access values

public:
  static const double defaultTolerance_;        ///< Default tolerance

  enum Type {SIMPLE, PAUSE};                    ///< Bre

  BreakPoint(double time = 0.0, Type type = SIMPLE)
    : time_(time),
      type_(type)
  {}

  BreakPoint(double time, int type)
    : time_(time),
      type_(type == 1 ? PAUSE : SIMPLE)
  {}

  BreakPoint(const BreakPoint& right)
    : time_(right.time_),
      type_(right.type_)
  {}

  BreakPoint &operator=(const BreakPoint &b)
  {
    time_ = b.time_;
    type_ = b.type_;

    return *this;
  }

  bool operator==(const BreakPoint& rhs) const
  {
     double diff = std::fabs(time_ - rhs.time_);
     return (diff <= defaultTolerance_);
  }

  double value() const { return time_; }
  Type bptype() const { return type_; }

  void setType (int type) { type_ = (type == 1 ? PAUSE : SIMPLE); }
  void setType (Type type) { type_ = type; }

private:
  double        time_;                  ///< Breakpoint time
  Type          type_;                  ///< SIMPLE or PAUSE
};

struct BreakPointLess
{
public:
  BreakPointLess(double tolerance) : tolerance_(tolerance) {}
  bool operator()(const BreakPoint &l, const BreakPoint &r) const { return l.time_ < r.time_ && std::fabs(r.time_ - l.time_) > tolerance_; }
  bool operator()(const BreakPoint &l, const double &d) const { return l.time_ < d && std::fabs(d - l.time_) > tolerance_; }
  double tolerance_;                            ///< Tolerance
};

struct BreakPointEqual
{
public:
  BreakPointEqual(double tolerance) : tolerance_(tolerance) {}
  bool operator()(const BreakPoint &l, const BreakPoint &r) const { return std::fabs(r.time_ - l.time_) <= tolerance_; }
  bool operator()(const BreakPoint &l, const double &d) const { return std::fabs(d - l.time_) >= tolerance_; }
  double tolerance_;                            ///< Tolerance
};

} // namespace Util
} // namespace Xyce

#endif

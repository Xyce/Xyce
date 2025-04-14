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

//-----------------------------------------------------------------------------
// Filename      : N_IO_CmdParse.h
//
// Purpose       : This file defines the command line parser class.
//
// Special Notes :
//
//
// Creator       : Eric Keiter
//
// Creation Date : 06/17/99
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_HangingResistor_h
#define Xyce_N_IO_HangingResistor_h

#include <iosfwd>
#include <string>
#include <map>

#include <N_IO_fwd.h>

namespace Xyce {
namespace IO {

class HangingResistor {
public:
  HangingResistor()
    : netlistCopy_(false),
      oneTerm_(false),
      noDCPath_(false),
      oneTermRes_(""),
      noDCPathRes_("")
  {}

  HangingResistor(const HangingResistor &right)
    : netlistCopy_(right.netlistCopy_),
      oneTerm_(right.oneTerm_),
      noDCPath_(right.noDCPath_),
      oneTermRes_(right.oneTermRes_),
      noDCPathRes_(right.noDCPathRes_)
  {}

  bool getNetlistCopy() const
  {
    return netlistCopy_;
  }

  bool getOneTerm() const
  {
    return oneTerm_;
  }

  bool getNoDCPath() const
  {
    return noDCPath_;
  }

  const std::string &getOneTermRes() const
  {
    return oneTermRes_;
  }

  const std::string &getNoDCPathRes() const
  {
    return noDCPathRes_;
  }

  void setNetlistCopy(bool netlistCopyArg) 
  {
    netlistCopy_ = netlistCopyArg;
  }

  void setOneTerm(bool oneTermArg) 
  {
    oneTerm_ = oneTermArg;
  }

  void setNoDCPath(bool noDCPathArg) 
  {
    noDCPath_ = noDCPathArg;
  }

  void setOneTermRes(const std::string & oneTermResArg) 
  {
    oneTermRes_ = oneTermResArg;
  }

  void setNoDCPathRes(const std::string & noDCPathResArg) 
  {
    noDCPathRes_ = noDCPathResArg;
  }

private:
  bool        netlistCopy_;
  bool        oneTerm_;
  bool        noDCPath_;
  std::string oneTermRes_;
  std::string noDCPathRes_;
};

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_HangingResistor_h

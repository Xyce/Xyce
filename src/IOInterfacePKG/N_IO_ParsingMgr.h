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
//
// Purpose        : Provide a parsing manager for Xyce that allows a limited
//                  set of HSPICE language features to be used in Xyce, on a
//                  user-configurable basis.
//
// Special Notes  :
//
// Creator        : Pete Sholander, SNL
//
// Creation Date  : 04/09/19
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_ParsingMgr_h
#define Xyce_N_IO_ParsingMgr_h

namespace Xyce {
namespace IO {

namespace RedefinedParamsSetting {
  enum RedefinedParamsSetting {ERROR, IGNORE, WARNING, USELAST, USELASTWARN, USEFIRST, USEFIRSTWARN};
}

//-----------------------------------------------------------------------------
// Class         : ParsingMgr
// Purpose       :
// Special Notes :
// Creator       : Pete Sholander, SNL
// Creation Date : 04/09/19
//-----------------------------------------------------------------------------
class ParsingMgr
{
public:
   ParsingMgr(const CmdParse &command_line);
  ~ParsingMgr();

private:
  ParsingMgr(const ParsingMgr & rhs);
  ParsingMgr &operator=(const ParsingMgr &rhs);

public:
  bool setParserOptions(const Util::OptionBlock & option_block);

  bool getHspiceExtFlag() const
  {
    return hspiceExtFlag_;
  }

  bool getModelBinningFlag() const
  {
    return modelBinningFlag_;
  }

  double getLengthScale () const
  {
    return lengthScale_;
  }

  char getSeparator() const
  {
    return ((useHspiceSeparator_)?('.'):(':'));
  }

  int getRedefinedParams() const
  {
    return redefinedParams_;
  }

  bool getImplicitSubcktMultiplier () const
  {
    return implicitSubcktMultiplier_;
  }

  bool getLocalVariationFlag () const
  {
    return localVariationFlag_;
  }

private:
  bool     hspiceExtFlag_;  // was the -hspice-ext command line option specified?
  bool     useHspiceUnits_; // was units or all specified for -hspice-ext ?
  bool     useHspiceMath_;  // was math or all specified for -hspice-ext ?
  bool     useHspiceSeparator_;  // was separator or all specified for -hspice-ext ?
  bool     modelBinningFlag_;
  double   lengthScale_;
  bool     redefinedParamsFlag_;
  int      redefinedParams_;
  bool     implicitSubcktMultiplierFlag_;
  bool     implicitSubcktMultiplier_;
  bool     localVariationFlagExists_;
  bool     localVariationFlag_;
};

bool registerPkgOptionsMgr(ParsingMgr & parsing_manager, PkgOptionsMgr &options_manager);

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_ParsingMgr_h

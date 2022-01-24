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

//-----------------------------------------------------------------------------
//
// Purpose        : Setup to register input options with pkg's that request them
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/28/03
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_PkgOptionsMgr_h
#define Xyce_N_IO_PkgOptionsMgr_h

// ---------- Standard Includes ----------

#include <string>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>
#include <N_IO_fwd.h>

#include <N_UTL_Param.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : PkgOptionsReg
// Purpose       : Abstract IF for pkg option registration
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/28/03
//-----------------------------------------------------------------------------
struct PkgOptionsReg
{
  virtual ~PkgOptionsReg()
  {}

  virtual bool operator()(const Util::OptionBlock & options) = 0;
};

template <class R>
struct RegisterOptions : public PkgOptionsReg
{
  typedef bool (R::*F)(const Util::OptionBlock &option_block);

  RegisterOptions(R &r, F f)
    : r_(r),
      f_(f)
  {}

  virtual ~RegisterOptions()
  {}

  bool operator()(const Util::OptionBlock &option_block)
  {
    (r_.*f_)(option_block);

    return true;
  }

  R &   r_;
  F     f_;
};

 template <class R>
   PkgOptionsReg *createRegistrationOptions(R &r, bool (R::*f)(const Util::OptionBlock &option_block))
 {
   return new IO::RegisterOptions<R>(r, f);
 }


//-----------------------------------------------------------------------------
// Class         : PkgOptionsMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/28/03
//-----------------------------------------------------------------------------
class PkgOptionsMgr
{
public:
  typedef std::multimap<std::string, PkgOptionsReg *, LessNoCase> ProcessorMap;
  typedef std::multimap<std::string, Util::OptionBlock, LessNoCase> OptionsMap;
  typedef std::map<std::string, ParseFunction, LessNoCase> CommandParserMap;

  PkgOptionsMgr();
  ~PkgOptionsMgr();

  bool addOptionsProcessor(const std::string &option_block_name, PkgOptionsReg *processor);

  bool addCommandProcessor(const std::string &option_block_name, PkgOptionsReg *processor) {
    return addOptionsProcessor(option_block_name, processor);
  }

  bool submitOptions(const Util::OptionBlock & options);

  void addCommandParser(const std::string &name, ParseFunction parse_function)
  {
    commandParserMap_[name] = parse_function;
  }

  ParseFunction getCommandParser(const std::string &name) const
  {
    CommandParserMap::const_iterator it = commandParserMap_.find(name);
    if (it != commandParserMap_.end())
      return (*it).second;
    return 0;
  }

  Util::OptionsMetadataMap::mapped_type &addOptionsMetadataMap(const std::string &name)
  {
    return optionsMetadata_[name];
  }

  const Util::OptionsMetadataMap::mapped_type *findOptionsMetadata(const std::string &name) const
  {
    Util::OptionsMetadataMap::const_iterator options_it = optionsMetadata_.find(name);
    if (options_it != optionsMetadata_.end())
      return &(*options_it).second;
    return 0;
  }

private:
  ProcessorMap                  processorMap_;
  OptionsMap                    optionsMap_;
  CommandParserMap              commandParserMap_;
  Util::OptionsMetadataMap      optionsMetadata_;
};

// Extract the parameters from a netlist .OPTIONS line held in parsedLine.
bool extractOptionsData     (PkgOptionsMgr &options_manager, CircuitBlock &circuit_block, const std::string &netlist_filename, const TokenVector &parsed_line);
bool extractDotDataStatement(PkgOptionsMgr &options_manager, CircuitBlock &circuit_block, const std::string &netlist_filename, const TokenVector &parsed_line);

struct ParamWarning : public Report::UserWarning
{
  ParamWarning(const Util::OptionBlock &option_block, const Util::Param &param);
};

struct ParamError : public Report::UserError
{
  ParamError(const Util::OptionBlock &option_block, const Util::Param &param);
};

} // namespace IO
} // namespace Xyce

#endif

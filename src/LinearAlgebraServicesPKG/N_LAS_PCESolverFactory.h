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
// Purpose        : Solver Factory for PCE 
//
// Special Notes  :
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 6/27/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_LAS_PCESolverFactory_h
#define Xyce_N_LAS_PCESolverFactory_h

// ---------- Standard Includes ----------

#include <string>

#include <N_LOA_fwd.h>
#include <N_LAS_fwd.h>

#include <N_LAS_SolverFactory.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_OptionBlock.h>

#include <Teuchos_RCP.hpp>

// ----------  Fwd Declares  -------------

namespace Xyce {
namespace Linear {

//-----------------------------------------------------------------------------
// Class         : PCESolverFactory
// Purpose       : 
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
class PCESolverFactory : public SolverFactory
{
public:
  // Basic Constructor, sets solver factory options.
  PCESolverFactory( Linear::Builder &builder );

  // Destructor
  virtual ~PCESolverFactory() {} 

  // Creates a new solver.
  Solver * create( Util::OptionBlock & options, Problem & problem, const IO::CmdParse & command_line) const;

  // Register the application system loader
  void registerPCELoader( const Teuchos::RCP<Loader::PCELoader>& pceLoaderPtr ) 
    { pceLoaderPtr_ = pceLoaderPtr; }

  // Register the PCE builder 
  void registerPCEBuilder( const Teuchos::RCP<PCEBuilder>& pceBuilder ) 
    { pceBuilderPtr_ = pceBuilder; }

  void setNumSamples(int numS)
  { numSamples_ = numS; }

  void setParameterOuterLoop (bool paramsOL)
    { paramsOuterLoop_ = paramsOL; }

private:
  Builder &                     builder_;
  Teuchos::RCP<Loader::PCELoader> pceLoaderPtr_;
  Teuchos::RCP<PCEBuilder> pceBuilderPtr_;
  Teuchos::RCP<Util::OptionBlock> optionBlock_;

  // Copy constructor.
  PCESolverFactory( const PCESolverFactory& pf );

  int numSamples_;
  bool paramsOuterLoop_;
};

} // namespace Linear
} // namespace Xyce

#endif

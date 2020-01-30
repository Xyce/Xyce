//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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
// Purpose       :
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date : 01/24/08
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_RegisterAnalysis_h
#define Xyce_N_ANP_RegisterAnalysis_h

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_LAS_fwd.h>
#include <N_LOA_fwd.h>
#include <N_NLS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_Factory.h>

namespace Xyce {
namespace Analysis {

class AnalysisCreatorRegistry {
private:
  typedef Util::Factory<AnalysisBase, void> Factory;
  typedef std::vector<Factory *> Registry;

public:
  ~AnalysisCreatorRegistry()
  {
    for (Registry::iterator it = registry_.begin(), end = registry_.end(); it != end; ++it)
      delete (*it);
  }
  
  void add(Factory *factory)
  {
    registry_.push_back(factory);
  }

  Registry::const_iterator begin() const
  {
    return registry_.begin();
  }

  Registry::const_iterator end() const
  {
    return registry_.end();
  }
  
private:
  Registry      registry_;
};

class ProcessorCreatorRegistry
{
private:
  typedef Util::Factory<ProcessorBase, void> Factory;
  typedef std::vector<Factory *> Registry;

public:
  ~ProcessorCreatorRegistry()
  {
    for (Registry::iterator it = registry_.begin(), end = registry_.end(); it != end; ++it)
      delete (*it);
  }

  void add(Factory *factory)
  {
    registry_.push_back(factory);
  }

  Registry::const_iterator begin() const
  {
    return registry_.begin();
  }

  Registry::const_iterator end() const
  {
    return registry_.end();
  }

private:
  Registry      registry_;
};


//-----------------------------------------------------------------------------
// Class         : FactoryBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Jan 27 14:49:12 2014
//-----------------------------------------------------------------------------
///
/// The FactoryBlock contains parameters needed by the analysis creation functions.  This allows additional parameter to
/// be added without the need to change the interface.
///
struct FactoryBlock
{
  //-----------------------------------------------------------------------------
  // Function      : FactoryBlock
  // Purpose       :
  // Special Notes :
  // Scope         : public
  // Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
  // Creation Date : Fri Mar 14 13:05:27 2014
  //-----------------------------------------------------------------------------
  ///
  /// The FactoryBlock constructs serves to pass data to the device factory functions
  ///
  /// @invariant These references must exist through the execution of Xyce
  ///
  FactoryBlock(
    AnalysisCreatorRegistry &           analysis_registry,
    ProcessorCreatorRegistry &          processor_registry,
    IO::PkgOptionsMgr &                 options_manager,
    AnalysisManager &                   analysis_manager,
    IO::OutputMgr &                     output_manager,
    Linear::System &                    linear_system,
    Nonlinear::Manager &                nonlinear_manager,
    Loader::Loader &                    loader,
    Device::DeviceMgr &                 device_manager,
    Linear::Builder &                   builder,
    Topo::Topology &                    topology,
    IO::InitialConditionsManager &      initial_conditions_manager,
    IO::RestartMgr &                    restart_manager)
    : analysisRegistry_(analysis_registry),
      processorRegistry_(processor_registry),
      optionsManager_(options_manager),
      analysisManager_(analysis_manager),
      outputManager_(output_manager),
      linearSystem_(linear_system),
      nonlinearManager_(nonlinear_manager),
      loader_(loader),
      deviceManager_(device_manager),
      builder_(builder),
      topology_(topology),
      initialConditionsManager_(initial_conditions_manager),
      restartManager_(restart_manager)
  {}

  AnalysisCreatorRegistry &             analysisRegistry_;
  ProcessorCreatorRegistry &            processorRegistry_;
  IO::PkgOptionsMgr &                   optionsManager_;
  AnalysisManager &                     analysisManager_;
  IO::OutputMgr &                       outputManager_;
  Linear::System &                      linearSystem_;
  Nonlinear::Manager &                  nonlinearManager_;
  Loader::Loader &                      loader_;
  Device::DeviceMgr &                   deviceManager_;
  Linear::Builder &                     builder_;
  Topo::Topology &                      topology_;
  IO::InitialConditionsManager &        initialConditionsManager_;
  IO::RestartMgr &                      restartManager_;
};

void registerAnalysisFactory(FactoryBlock &factory_block);

inline void addAnalysisFactory(FactoryBlock &factory_block, Util::Factory<AnalysisBase, void> *factory)
{
  factory_block.analysisRegistry_.add(factory);
}

inline void addProcessorFactory(FactoryBlock &factory_block, Util::Factory<ProcessorBase, void> *factory)
{
    factory_block.processorRegistry_.add(factory);
}

} // namespace Analysis
} // namespace Xyce

#endif // Xyce_N_ANP_RegisterAnalysis_h

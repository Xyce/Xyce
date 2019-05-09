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
// Purpose       : This file contains the functions which define the time
//                 domain & integration algorithm classes.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_ANP_RegisterAnalysis.h>

#include <N_ANP_AC.h>
#include <N_ANP_DCSweep.h>
#include <N_ANP_HB.h>
#include <N_ANP_MOR.h>
#include <N_ANP_MPDE.h>
#include <N_ANP_NOISE.h>
#include <N_ANP_ROL.h>
#include <N_ANP_Step.h>
#include <N_ANP_Sampling.h>
#include <N_ANP_EmbeddedSampling.h>
#include <N_ANP_Transient.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : registerAnalysisFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : David G. Baur  Raytheon  Sandia National Laboratories 1355
// Creation Date : Mon Mar 16 14:26:41 2015
//-----------------------------------------------------------------------------
///
/// Registers the analysis factories.
///
/// In the future this is intended to become a pluggable factory allowing new analysis type to be added with no invasion
/// into the rest of the analysis code.  Whether anyone will actually create an analysis outside the system is academic
/// and creating this interface will result in defining the roles and interactions between the various analyses.
///
/// @invariant
///
/// @param netlist_filename 
/// @param options_manager 
/// @param analysis_manager 
/// @param linear_system 
/// @param nonlinear_manager 
/// @param device_manager 
/// @param builder 
/// @param topology 
///
///
void registerAnalysisFactory(
  FactoryBlock &        factory_block)
{
  registerDCSweepFactory(factory_block);
  registerACFactory(factory_block);
  registerTransientFactory(factory_block);
  registerHBFactory(factory_block);
  registerMPDEFactory(factory_block);
  registerNOISEFactory(factory_block);
  registerMORFactory(factory_block);
  registerStepFactory(factory_block);
  registerSamplingFactory(factory_block);
  registerEmbeddedSamplingFactory(factory_block);
  registerROLFactory(factory_block);
}

} // namespace Analysis
} // namespace Xyce


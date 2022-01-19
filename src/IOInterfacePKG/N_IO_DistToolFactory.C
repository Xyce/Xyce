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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/12/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Xyce Includes   ----------

#include <N_IO_DistToolFactory.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_DistToolDefault.h>
#include <N_IO_DistToolFlatRoundRobin.h>
#include <N_IO_DistToolDevBalanced.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_fwd.h>
#include <N_IO_ParsingMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : DistToolFactory::create
// Purpose       :
// Special Notes : Static
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 05/24/04
//-----------------------------------------------------------------------------
DistributionTool * 
  DistToolFactory::create(
    Parallel::Communicator*                  pdsCommPtr,
    Util::OptionBlock &                      distOptions,
    CircuitBlock &                           circuitBlock,
    std::map<std::string,FileSSFPair>      & ssfMap,
    std::map<std::string, IncludeFileInfo> & iflMap,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
    const ParsingMgr                       & parsing_manager
    )
{
  int strategy = Xyce::IO::DistStrategy::DEFAULT;

  Util::ParamList::const_iterator itPI = distOptions.begin();
  Util::ParamList::const_iterator endPI = distOptions.end();
  for( ; itPI != endPI; ++itPI )
  {
    if( itPI->uTag() == "STRATEGY" )
    {
      strategy = itPI->getImmutableValue<int>();
    }
  }

  CircuitContext* circuitContext = circuitBlock.getCircuitContextPtr();
  const std::vector< std::string >& subcktList = circuitContext->getSubcktList();

  // If there aren't any subcircuits, then the round robin scheme can be used assign devices to processors.
  // This interleaves communication and computation to improve distribution performance.
  if (strategy == Xyce::IO::DistStrategy::FLAT_ROUND_ROBIN && (subcktList.size()!=0 || iflMap.size()!=0))
  {
    Report::UserWarning() << "Round-robin distribution strategy is not valid for this circuit!" << std::endl;
    strategy = Xyce::IO::DistStrategy::DEFAULT;
  }

  if (Parallel::is_parallel_run(pdsCommPtr->comm()))
  {
/*
    if (pdsCommPtr->procID() == 0)
    {
      std::cout << "Top level circuit has " << subcktList.size() << " subcircuits and " << circuitContext->getTotalDeviceCount() << " devices." << std::endl;
      std::cout << "Top level circuit has " << subcktList.size() << " subcircuits and " << circuitContext->getTotalLinearDeviceCount() << " linear devices." << std::endl;
      std::vector< std::string >::const_iterator ilIter = subcktList.begin();
      for ( ; ilIter != subcktList.end(); ilIter++ )
      {
        // Set the context for this subcircuit.
        circuitContext->setContext( *ilIter );

        Xyce::dout() << "For subcircuit " << *ilIter << " the total device count is : " << circuitContext->getTotalDeviceCount() << std::endl;
        Xyce::dout() << "For subcircuit " << *ilIter << " the total linear device count is : " << circuitContext->getTotalLinearDeviceCount() << std::endl;

        // Restore context before next loop.
        circuitContext->restorePreviousContext();
      }
      std::cout << "Circuit has " << iflMap.size() << " include files. " << std::endl;
    }
*/
    // Send the strategy to all processors. 
    pdsCommPtr->bcast( &strategy, 1, 0);
  }

  DistributionTool * ret = 0;

  // Now create the distribution tool based on the distribution strategy.
  switch ( strategy )
  {
    case DistStrategy::DEFAULT :
    {
      ret = new DistToolDefault(pdsCommPtr, circuitBlock, ssfMap, iflMap, parsing_manager);
      break;
    }
    case DistStrategy::FLAT_ROUND_ROBIN :
    {
      ret = new DistToolFlatRoundRobin(pdsCommPtr, circuitBlock, ssfMap, iflMap, externalNetlistParams, parsing_manager);
      break; 
    }
    case DistStrategy::DEVICE_BALANCED:
    {
      ret = new DistToolDevBalanced(pdsCommPtr, circuitBlock, ssfMap, iflMap, parsing_manager);
      break; 
    }
    default :
    {
      ret = new DistToolDefault(pdsCommPtr, circuitBlock, ssfMap, iflMap, parsing_manager);
      break;
    }
  }

  return ret;

}

} // namespace IO
} // namespace Xyce

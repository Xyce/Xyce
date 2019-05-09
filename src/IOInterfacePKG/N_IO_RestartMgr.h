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
// Purpose        : Control storage and writeback for restarts
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 7/19/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_RestartMgr_h
#define Xyce_N_IO_RestartMgr_h

#include <string>
#include <map>

#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace IO {

// The restart I/O should replaced with a restart data sink object.  Rather than having this object hold the data, a
// sink/source object should be created that handles this.

//-----------------------------------------------------------------------------
// Class         : RestartMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
class RestartMgr
{
public:
  RestartMgr();

  ~RestartMgr()
  {}

private:
  RestartMgr(const RestartMgr & right);
  RestartMgr &operator=(const RestartMgr & right);

public:
  // Restart flag accessor.
  bool isRestarting() const
  {
    return restartFlag_;
  }

  bool registerRestartOptions(const Util::OptionBlock & option_block, int proc_size, int proc_rank);
  bool registerTimeintOptions(const Util::OptionBlock & option_block);

  double getInitialRestartInterval() const
  {
    return initialRestartInterval_;
  }

  const IntervalVector &getRestartIntervals() const
  {
    return restartIntervalPairs_;
  }

  // bool registerNodePartitioning(std::map<std::string,int> & nodeMap )
  // {
  //   npMap_ = nodeMap;
  //   return true;
  // }

  const std::string &path() const
  {
    return restartFileName_;
  }

  const std::string &getJobName() const
  {
    return restartJobName_;
  }

  bool getPack() const
  {
    return pack_;
  }

  bool getPrintOptions() const
  {
    return printOptions_;
  }

  bool
  restoreRestartData(
    N_PDS_Comm &                  comm,
    Topo::Topology &              topology,
    Analysis::AnalysisManager &   analysis_manager,
    Device::DeviceMgr &           device_manager,
    const std::string &           path);

private:
  bool                          restartFlag_;
  std::string                   restartFileName_;
  std::string                   restartJobName_;
  double                        initialRestartInterval_;
  IntervalVector                restartIntervalPairs_;
  // std::map<std::string, int>    npMap_;
  bool                          pack_;
  bool                          printOptions_;

  Util::OptionBlock             savedTimeintOB_;
};

bool registerPkgOptionsMgr(RestartMgr &restart_manager, PkgOptionsMgr &options_manager, int proc_size, int proc_rank);

bool
dumpRestartData(
  N_PDS_Comm &                  comm,
  Topo::Topology &              topology,
  Analysis::AnalysisManager &   analysis_manager,
  Device::DeviceMgr &           device_manager,
  const std::string &           job_name,
  bool                          pack,
  double                        time);

} // namespace IO
} // namespace Xyce

#endif

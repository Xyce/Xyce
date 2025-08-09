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
// Purpose        : Source base  lasses.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Source_h
#define Xyce_N_DEV_Source_h

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceInstance.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : SourceInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
class SourceInstance : public DeviceInstance
{
public:
  SourceInstance(
     const InstanceBlock &       IB,
     ParametricData<void> &      parametric_data,
     const FactoryBlock &        factory_block);

  ~SourceInstance();

  bool analyticSensitivityAvailable (const std::string & paramName);

private:
  SourceInstance(const SourceInstance &);
  SourceInstance &operator=(const SourceInstance &);

public:
  void setFastSourceFlag(bool value);
  bool getFastSourceFlag() const;

  double period() const;

  virtual void setupBreakPoints();
  virtual bool getInstanceBreakPoints(std::vector<Util::BreakPoint> &breakPointTimes);
  virtual bool updateSource();

  virtual bool loadBVectorsforAC(double * bVecReal, double * bVecImag ) { return true; }


  virtual bool loadFreqBVector(double frequency,
                               std::vector<Util::FreqVecEntry>& BVecEntries) { return true; }

protected:
  int                   sourceType;             ///< type of source data
  SourceData *          tranSourceData_;
  SourceData *          acSourceData_;
  SourceData *          dcSourceData_;
};

} // namespace Device
} // namespace Xyce

#endif

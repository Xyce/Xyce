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

//-------------------------------------------------------------------------
//
// Purpose        :
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
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_Source.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Device {

// SourceInstance
//-----------------------------------------------------------------------------
// Function      : SourceInstance::SourceInstance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SourceInstance::SourceInstance(
  const InstanceBlock &         IB,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, parametric_data, factory_block),
    sourceType(_SIN_DATA),
    tranSourceData_(0),
    acSourceData_(0),
    dcSourceData_(0)
{}


//-----------------------------------------------------------------------------
// Function      : SourceInstance::~SourceInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
SourceInstance::~SourceInstance()
{}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::analyticSensitivityAvailable
//
// Purpose       : Independent sources have unorthodox parameter structures 
//                 for some source types.  As such it is necessary to overide
//                 the base class version of this function.
//
//                 Note, this function is specifically to support 
//                 analytical sensitivities in PWL sources.  The values in a 
//                 PWL source are associated with the "V" parameter, but come 
//                 into the device (post-parsing) with the index attached.  
//                 So, V0, V1, V2, etc.   That use case isn't handled by
//                 the base class function in N_DEV_DeviceEntity.
//
//                 Note, for sensitivity analysis to work on such parameters,
//                 it is also necessary to modify the DeviceMgr::getParamAndReduce 
//                 function to handle these types of exceptions.  
//                 As of this writing (August, 2025) that hasn't happened yet.
//                 Instead, a limited list of V# parameters have been added to
//                 the Vsrc and ISRC devices.  Once DeviceMgr::getParamAndReduce 
//                 is fixed, then those artificial params can go away.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 8/9/25
//-----------------------------------------------------------------------------
bool SourceInstance::analyticSensitivityAvailable (const std::string & paramName)
{
  bool returnValue = false;

  std::string paramName2 = ExtendedString( paramName ).toUpper();

  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    // check for special case of "V" parameters.  These are parameters that
    // are specified only as the single character "V" in the addPar calls, but
    // are allowed to be an indefinitely long vector for PWL sources.  The parser
    // assigns an index to each entry of the vector and appends that index to the name.
    // So, a PWL source with three entries will be given V0, V1, V2.

    int Vindex=-1;
    if ( std::string( paramName2 ,0,1) == "V" )
    {
      int size = paramName2.size();
      if (size > 1)
      {
        std::string indexNumberStr = std::string( paramName2, 1, size-1);

        std::size_t pos{};
        try
        {
          Vindex = std::stoi(indexNumberStr, &pos);
          //std::cout << "paramName = " << paramName << " Vindex = " << Vindex << std::endl;
          returnValue = true;
        }
        catch ( std::invalid_argument const& ex )
        {
          DevelFatal(*this).in("DeviceEntity::analyticSensitivityAvailable") << "Unrecognized parameter " << paramName;
          //std::cout << "paramName = " << paramName << " caught an error" << std::endl;
          //returnValue = false;
        }
      }
    }
    else
    {
      DevelFatal(*this).in("DeviceEntity::analyticSensitivityAvailable") << "Unrecognized parameter " << paramName;
    }
  }
  else
  {
    const Descriptor &param = *(*p_i).second;
    returnValue = param.getAnalyticSensitivityAvailable();
  }

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::setFastSourceFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/27/04
//-----------------------------------------------------------------------------
void SourceInstance::setFastSourceFlag (bool value)
{
  if (tranSourceData_ != 0)
  {
    tranSourceData_->setFastTimeScaleFlag(value);
  }
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::getFastSourceFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/1/04
//-----------------------------------------------------------------------------
bool SourceInstance::getFastSourceFlag() const
{
  bool flag = false;
  if (tranSourceData_ != 0)
  {
    flag = tranSourceData_->getFastTimeScaleFlag();
  }
  return flag;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::period
// Purpose       : Time period of periodic sources
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/31/04
//-----------------------------------------------------------------------------
double SourceInstance::period() const
{
  double per = 0.0;
  if (tranSourceData_ != 0)
  {
    per = tranSourceData_->period();
  }
  return per;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::setupBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/16/2019
//-----------------------------------------------------------------------------
void SourceInstance::setupBreakPoints()
{
  bool fast_source_flag = false;

  if (!getSolverState().blockAnalysisFlag_)
  {
    fast_source_flag = getFastSourceFlag();
  }

  if (!fast_source_flag && tranSourceData_ != 0)
  {
    tranSourceData_->setupBreakPoints();
  }
  
  return;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/05/06
//-----------------------------------------------------------------------------
bool SourceInstance::getInstanceBreakPoints(std::vector<Util::BreakPoint> & breakPointTimes)
{
  bool fast_source_flag = false;

  if (!getSolverState().blockAnalysisFlag_)
  {
    fast_source_flag = getFastSourceFlag();
  }

  bool tmpBool = true;
  if (!fast_source_flag && tranSourceData_ != 0)
  {
    tmpBool = tranSourceData_->getBreakPoints(breakPointTimes);
  }
  
  return tmpBool;
}

//-----------------------------------------------------------------------------
// Function      : SourceInstance::updateSource
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/05/06
//-----------------------------------------------------------------------------
bool SourceInstance::updateSource()
{
  if (tranSourceData_ != 0)
  {
    tranSourceData_->updateSource();
  }
  if (dcSourceData_ != 0)
  {
    dcSourceData_->updateSource();
  }
  if (acSourceData_ != 0)
  {
    acSourceData_->updateSource();
  }
  
  return true;
}

} // namespace Device
} // namespace Xyce

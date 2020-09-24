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

//-------------------------------------------------------------------------
//
// Purpose        : This class only contains the instance processParams
//                  function.  This was a big enough, confusing enough,
//                  function, that I decided to give it its own file.  Most
//                  of the functions that are called from here can be found
//                  in N_DEV_2DPDESetup.C.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/15/03
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------
#include <N_DEV_2DPDE.h>
#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_DeviceOptions.h>
#include <N_ERH_ErrorMgr.h>
#include <N_DEV_Message.h>

namespace Xyce {
namespace Device {
namespace TwoDPDE {

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
//
// Purpose       : This function processes the user-specified parameters.
//                 This includes setting up the mesh.
//
// Special Notes : There are three modes in which processParams can
//                 be called:
//
//                 1) Initialization of the class, for a standard simulation.
//
//                 2) perturbing away from the original as part of a sensitivity
//                     calculation, in order to obtain a numerical df/dp.
//
//                 3) restoring the device to its normal, unperturbed state.
//
//                 Call (1) is obviously the most common, and happens at
//                 the beginning of any simulation.  This call will have
//                 no parameter argument.
//
//                 Calls (2) and (3) only happen for sensitivity calculations,
//                 and they always come in pairs.  Anytime a (2) style call
//                 happens, it should be followed shortly by a (3) call.
//
//                 Call (2) will have a parameter argument, specifying which
//                 param to perturb.  Call (3) doesn't have a param argument,
//                 and for most devices is identical to a call (1).  However,
//                 for the PDE case, it doesn't makes sense to reallocate all
//                 the mesh stuff, so it just restores old values of the mesh
//                 from a copy.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  bool bsuccess = true;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::processDopingParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/04/04
//-----------------------------------------------------------------------------
bool Instance::processDopingParams (Param & ndParam, std::string param)
{
  ExtendedString tagES = ndParam.tag ();
  tagES.toLower ();

  // Start of the new doping section:
  //string tmpTag (tagES,0,3);
  if (std::string(tagES,0,3) == "reg" && ndParam.given())
      //if dope region param, continue:
  {
    // the new metadata doesn't have a limitation on the number of
    // regions.

    // first find the "." character.
    int periodLocation = 0;
    bool foundPeriod = false;
    for (int iloc=0;iloc<tagES.size();++iloc)
    {
      std::string singleChar(tagES,iloc,1);
      if (singleChar == ".")
      {
        periodLocation = iloc;
        foundPeriod = true;
        break;
      }
    }
    if (!foundPeriod)
    {
      DevelFatal(*this).in("Instance::processDopingParams")
        << "The region specification needs a period.";
    }

    std::string regionName(tagES,0,periodLocation);

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "tagES = "<< tagES;
      Xyce::dout() << "  regionName = "<< regionName;
      Xyce::dout() << "  param: "  << ndParam.stringValue() << std::endl;
    }

    if ( dopeInfoMap.find(regionName) == dopeInfoMap.end() )
    {
      dopeInfoMap[regionName] = new DopeInfo(getSolverState());
      dopeInfoMap[regionName]->name = regionName;
    }

    // The actual processing of the parameter is handled in the DopeInfo
    // class.
    dopeInfoMap[regionName]->processParam (ndParam, param, *this);

  }
  // End of the new doping section.

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Instance::processElectrodeParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/18/04
//-----------------------------------------------------------------------------
bool Instance::processElectrodeParams (Param & ndParam)
{

  ExtendedString tagES = ndParam.tag ();
  tagES.toLower ();

  if (std::string(tagES,0,4) == "node" && ndParam.given() )
  {
    int periodLocation = 5;

    // Now that we've determined that this is an electrode parameter,
    // ( and not one of the older "legacy" parameters )
    // check if an electrode with this name already exists in the
    // map.

    // If so, just modify the existing electrode.  If not, add it
    //  to the map, with a set of default params (created by the
    // default constructor).
    //
    // get the name by searching for the ".", just like in the case for
    // the region specification.
    // the new metadata doesn't have a limitation on the number of
    // nodes.

    // first find the "." character.
    periodLocation = 0;
    bool foundPeriod = false;
    for (int iloc=0;iloc<tagES.size();++iloc)
    {
      std::string singleChar(tagES,iloc,1);
      if (singleChar == ".")
      {
        periodLocation = iloc;
        foundPeriod = true;
        break;
      }
    }

    // if no period is found, then report an error.
    // At this point, node1 and node1bc style specifications have
    // already been trapped for.
    if (!foundPeriod)
    {
      DevelFatal(*this).in("Instance::processParams")
        << "The node specification needs a period";
    }

    std::string nodeName(tagES,0,periodLocation);
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "tagES      = "<< tagES      << std::endl;
      Xyce::dout() << "nodeName = "<< nodeName << std::endl;
    }

    if ( electrodeMap.find(nodeName) == electrodeMap.end() )
    {
      PDE_2DElectrode * elPtr = new PDE_2DElectrode ();
      elPtr->nodeName = nodeName;
      electrodeMap[nodeName] = elPtr;
    }

    // now do the actual processing of the param.
    int tagSize = tagES.size();

    if (tagSize > periodLocation+1)
    {
      std::string tmpParam(tagES,periodLocation+1,tagSize);

      if (tmpParam == "name")
      {
        DeviceInterfaceNode dINode;
        ExtendedString dIName = ndParam.stringValue();
        dIName.toUpper ();

        dINode.eName = dIName;
        dINode.nName = nodeName;
        dINode.given = ndParam.given();
        dINode.index = 0;

        if (dINode.given) ++numElectrodes;
        if (dINode.given) dIVec.push_back(dINode);
      }

      if (tmpParam == "bc")
      {
        ExtendedString nmName = ndParam.stringValue();
        nmName.toUpper ();
        if (nmName == "NEUMANN" ||
            nmName == "MIXED" )
        {
          tmpBCmap[nodeName] = nmName;
          Xyce::dout() << "found a neumann.  name = " << tagES << std::endl;
        }
      }

      if (tmpParam == "side")
      {
        if (ndParam.given())
        {
          ExtendedString tmpstr = ndParam.stringValue();
          tmpstr.toLower ();
          electrodeMap[nodeName]->side = tmpstr;
          electrodeMap[nodeName]->sideGiven = true;
        }
      }

      if (tmpParam == "start")
      {
        electrodeMap[nodeName]->start = ndParam.getImmutableValue<double>();
        electrodeMap[nodeName]->startGiven = ndParam.given ();
      }

      if (tmpParam == "end")
      {
        electrodeMap[nodeName]->end = ndParam.getImmutableValue<double>();
        electrodeMap[nodeName]->endGiven = ndParam.given ();
      }

      if (tmpParam == "material")
      {
        ExtendedString tmpName = ndParam.stringValue();
        tmpName.toLower();
        electrodeMap[nodeName]->material = tmpName;
        electrodeMap[nodeName]->materialGiven = ndParam.given ();
      }

      if (tmpParam == "oxidebndryflag")
      {
        if (ndParam.given() && ndParam.getImmutableValue<int>() == 1)
          electrodeMap[nodeName]->oxideBndryFlag = true;
        else
          electrodeMap[nodeName]->oxideBndryFlag = false;
      }

      if (tmpParam == "oxthick")
      {
        electrodeMap[nodeName]->oxthick = ndParam.getImmutableValue<double>();
      }

      if (tmpParam == "oxcharge")
      {
        electrodeMap[nodeName]->oxcharge = ndParam.getImmutableValue<double>();
      }
    }
  }

  return true;
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce

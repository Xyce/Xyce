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

//-------------------------------------------------------------------------
//
// Purpose        : This file mostly contains functions that are called
//                  once, during the initial setup phase of the 2D PDE
//                  device.  There are a couple of exceptions - the mesh
//                  resize functions are called during "re-set-up" phases
//                  of a sensitivity calculation.
//
//                  One very important setup function - processParams - is
//                  *not* in this file.  It gets its own file,
//                  N_DEV_2DPDEParam.C.
//
//                  All of the of the functions pertaining to the global
//                  and local ID setup are in this file.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/05/03
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
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_DEV_PDE_2DMesh.h>
#include <N_DEV_PDE_Electrode.h>

namespace Xyce {
namespace Device {
namespace TwoDPDE {

//-----------------------------------------------------------------------------
// Function      : Instance::doSensMeshResize
// Purpose       :
// Special Notes : Generally, this will be called for a mesh that was
//                 already scaled, so the resized mesh should be considered
//                 scaled as well.
//
//                 As should be obvious from the name, this function is
//                 designed for perturbing the mesh as part of a
//                 sensitivity calculation.  As such, a copy of the
//                 original mesh is saved, to be restored after the
//                 calculation is over.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
bool Instance::doSensMeshResize ()
{
  bool bsuccess = true;
  bool bs1 = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In Instance::doSensMeshResize." << std::endl;
  }

  // make a copy of the mesh.  This will need to be restored later.
  if (meshCopyContainerPtr == NULL)
  {
    meshCopyContainerPtr = new PDE_2DMesh (*meshContainerPtr);
  }
  else
  {
    *meshCopyContainerPtr = *meshContainerPtr;
  }

  // scale the new size:
  if (variablesScaled)
  {
    deviceLength /= scalingVars.x0;
    deviceWidth  /= scalingVars.x0;
  }

  // get the old size:
  double old_length = meshContainerPtr->getXMax () -
                      meshContainerPtr->getXMin ();

  double old_width  = meshContainerPtr->getYMax () -
                      meshContainerPtr->getYMin ();

  double lengthRatio  = deviceLength/old_length;
  double widthRatio   = deviceWidth/old_width;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  x0           = " << scalingVars.x0 << std::endl;
    Xyce::dout() << "  deviceWidth  = " << deviceWidth  << std::endl;
    Xyce::dout() << "  old_width    = " << old_width    << std::endl;
    Xyce::dout() << "  widthRatio   = " << widthRatio   << std::endl;
    Xyce::dout() << "  deviceLength = " << deviceLength << std::endl;
    Xyce::dout() << "  old_length   = " << old_length   << std::endl;
    Xyce::dout() << "  lengthRatio  = " << lengthRatio  << std::endl;
  }

  // first resize the mesh in the mesh class:
  meshContainerPtr->resizeMesh(deviceLength, deviceWidth);

  // Then update all the mesh sized stuff used in the Instance
  // class.
  meshContainerPtr->getXVector(xVec);
  meshContainerPtr->getYVector(yVec);
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    for (int i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << " x["<<i<<"] = " << xVec[i];
      Xyce::dout() << " y["<<i<<"] = " << yVec[i];
      Xyce::dout() << std::endl;
    }
  }

  // now update all the mesh stuff in the 2DPDE class:
  bs1 = setupBCEdgeAreas ();  bsuccess = bsuccess && bs1;
  bs1 = setupMinDXVector ();  bsuccess = bsuccess && bs1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with Instance::doSensMeshResize." << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::undoSensMeshResize
//
// Purpose       : This un-does the damage done by doSensMeshResize.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
bool Instance::undoSensMeshResize ()
{
  bool bsuccess = true;
  bool bs1 = true;
  // switch the mesh copy back into the official ptr.:
  PDE_2DMesh * tmpPtr;

  tmpPtr               = meshContainerPtr;
  meshContainerPtr     = meshCopyContainerPtr;
  meshCopyContainerPtr = tmpPtr;

  // Restore all the mesh sized stuff used in the Instance
  // class.
  meshContainerPtr->getXVector(xVec);
  meshContainerPtr->getYVector(yVec);

  // now update all the mesh stuff in the 2DPDE class:
  bs1 = setupBCEdgeAreas (); bsuccess = bsuccess && bs1;
  bs1 = setupMinDXVector (); bsuccess = bsuccess && bs1;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMesh
//
// Purpose       : This function should only be called once.  It handles
//                 most of the stuff associated with initializing the
//                 mesh class, and all the stuff in Instance
//                 that depends on the mesh.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
bool Instance::setupMesh ()
{

  bool bsuccess = true;

  ///////////////////////////////////////////////////////////////////
  // First straighten out the electrode map, if it exists.
  // The electrode map is one of the arguments that need to be passed into
  // the "internal" mesh setup, so it has to be corrected first.
  std::vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterV;

  if (!(electrodeMap.empty ()))
  {
    // First make the names in electrodeMap consistent with those in dIVec.
    for (iterV=first;iterV!=last; ++iterV)
    {
      if (!(iterV->given)) continue;

      if ( electrodeMap.find(iterV->nName) != electrodeMap.end () )
      {
        electrodeMap[iterV->nName]->name = iterV->eName;
      }
      else
      {
        DevelFatal(*this).in("Instance::doMeshBasedInitializations")
          << "can't find " << iterV->nName << " in the electrode Map";
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
       Xyce::dout() << "list of user-specified electrodes:" << std::endl;
       std::map<std::string, PDE_2DElectrode*>::iterator mapIter;
       std::map<std::string, PDE_2DElectrode*>::iterator mapStart = electrodeMap.begin();
       std::map<std::string, PDE_2DElectrode*>::iterator mapEnd = electrodeMap.end();

       // for ( mapIter = mapStart; mapIter != mapEnd; ++mapIter )
       // {
       //  Xyce::dout() << *(mapIter->second);
       // }
    }
  }

  ///////////////////////////////////////////////////////////////////
  // Allocate the mesh container.
  meshContainerPtr = new PDE_2DMesh(getDeviceOptions(), sgplotLevel);

  if (!given("MESHFILE"))
  {
    lout() << "No mesh file specified.  Setting meshfile=internal.msh\n" << std::endl;
  }

  ///////////////////////////////////////////////////////////////////
  // Now initialize the mesh either as internal or external.
  if (meshFileName != "internal" && meshFileName != "internal.msh")
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "Reading mesh file..." << std::endl;
    }
    usingInternalMesh = false;
    meshContainerPtr->initializeMesh (meshFileName);
    cylGeomFlag = meshContainerPtr->cylGeom;
  }
  else
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "Generating internal mesh..." << std::endl;
    }
    usingInternalMesh = true;

    std::string outputMeshFileName = outputName + ".msh";
    meshContainerPtr->initializeInternalMesh
      (numMeshPointsX,
       numMeshPointsY,
       deviceLength,
       deviceWidth,
       numElectrodes,
       outputMeshFileName,
       electrodeMap,
       cylGeomFlag);
  }

  ///////////////////////////////////////////////////////////////////
  numMeshPoints = meshContainerPtr->getNumNodes ();
  numMeshEdges  = meshContainerPtr->getNumEdges ();
  numMeshCells  = meshContainerPtr->getNumCells ();
  numMeshLabels = meshContainerPtr->getNumLabels ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n";
    Xyce::dout() << "Done setting up the mesh." << std::endl;
    Xyce::dout() << "  numMeshPoints      = " << numMeshPoints << "\n";
    Xyce::dout() << "  numMeshEdges       = " << numMeshEdges << "\n";
    Xyce::dout() << "  numMeshCells       = " << numMeshCells << "\n";
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupDINodes
// Purpose       : This function does some miscellaneous setup of the
//                 device interface nodes (boundary condition class).
//                 Mostly this is does some final, misc. cleanup.
//
//                 Should be called after setupMesh.  Should only be called once.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/15/03
//-----------------------------------------------------------------------------
bool Instance::setupDINodes ()
{
  bool bsuccess = true;

  std::vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterV;

  // loop through the boundary condition name vector and
  // check if it matches the boundary conditions given in the mesh file.

  for (iterV=first;iterV!=last; ++iterV)
  {
    ExtendedString tmpName = iterV->eName;
    tmpName.toUpper ();

    bool edgeLabelExist = meshContainerPtr->labelNameExist(tmpName);

    // If the device interface was given by the user(netlist)
    // (as a boundary condition), then if it doesn't exist,
    // then the user has made a mistake in setting up the input file and
    // the code should exit.
    if ((iterV->given))
    {
      if ( !(edgeLabelExist) )
      {
        meshContainerPtr->printLabels ();
        DevelFatal(*this).in("Instance::setupDINodes")
          << "The boundary condition label " << tmpName
          << " doesn't exist in the mesh file.\n";
      }
    }
    else // If this wasn't given, and doesn't exist, that just means it
         // is part of the default list of boundary condition names,
         // and should just be removed from the device interface vector.
         // This does NOT represent a mistake in the netlist.
    {
      if ( !(edgeLabelExist) )
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          Xyce::dout() << "Erasing DI: " << iterV->eName << std::endl;
        dIVec.erase (iterV);
      }
    }
  } // end of iterV loop.

  // Copy over a few of the material related variables from the electrode
  // class to the device interface node class.
  if (!(electrodeMap.empty ()))
  {
    // First make the names in electrodeMap consistent with those in dIVec.
    for (iterV=first;iterV!=last; ++iterV)
    {
      if (!(iterV->given)) continue;

      if ( electrodeMap.find(iterV->nName) != electrodeMap.end () )
      {
        // material stuff:
        iterV->material       = electrodeMap[iterV->nName]->material;
        iterV->materialGiven  = electrodeMap[iterV->nName]->materialGiven;
        iterV->oxideBndryFlag = electrodeMap[iterV->nName]->oxideBndryFlag;
        iterV->oxthick        = electrodeMap[iterV->nName]->oxthick;
        iterV->oxcharge       = electrodeMap[iterV->nName]->oxcharge;
      }
      else
      {
        DevelFatal(*this).in("Instance::doMeshBasedInitializations")
          << "can't find " << iterV->nName << " in the electrode Map";
      }
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << std::endl;
        Xyce::dout() << "name = " << iterV->eName << std::endl;
        Xyce::dout() << " material = " << iterV->material << std::endl;
        Xyce::dout() << " mat. given = " << iterV->materialGiven << std::endl;
        Xyce::dout() << " oxide boundary flag = " << iterV->oxideBndryFlag << std::endl;
        Xyce::dout() << " oxide thickness = " << iterV->oxthick << std::endl;
        Xyce::dout() << " oxide charge = " << iterV->oxcharge << std::endl;
      }
    }
  }

  // Loop over the boundaries, and check if each one is a
  // neumann, or mixed boundary condition for each variable.
  first = dIVec.begin ();
  last  = dIVec.end   ();
  for (iterV=first;iterV!=last; ++iterV)
  {
    ExtendedString tmpName = iterV->nName;
    tmpName.toLower ();

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      Xyce::dout() << "Testing the neumann stuff.  Name = " << tmpName << std::endl;

    if ( tmpBCmap.find(tmpName) != tmpBCmap.end () )
    {
      if (tmpBCmap[tmpName] == "NEUMANN")
      {
        iterV->neumannBCFlagV = true;
        iterV->neumannBCFlagN = true;
        iterV->neumannBCFlagP = true;
      }

      if (tmpBCmap[tmpName] == "MIXED")
      {
        iterV->neumannBCFlagV = false;
        iterV->neumannBCFlagN = true;
        iterV->neumannBCFlagP = true;
      }

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "Setting the neumann flags of " << tmpName << ":\n";
        Xyce::dout() << "  Vflag = " << iterV->neumannBCFlagV << std::endl;
        Xyce::dout() << "  Nflag = " << iterV->neumannBCFlagN << std::endl;
        Xyce::dout() << "  Pflag = " << iterV->neumannBCFlagP << std::endl;
      }
    }
  } // end of iterV loop.

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    first = dIVec.begin ();
    last  = dIVec.end   ();
    Xyce::dout() << "Final DI list: " << std::endl;
    for (iterV=first;iterV!=last; ++iterV)
    {
      Xyce::dout() << "DI name:" << iterV->eName;
      Xyce::dout() << "  The neumann flags are:" << std::endl;
      Xyce::dout() << "  Vflag =";
      if (iterV->neumannBCFlagV) Xyce::dout() <<" true." << std::endl;
      else                       Xyce::dout() <<" false." << std::endl;
      Xyce::dout() << "  Nflag =";
      if (iterV->neumannBCFlagN) Xyce::dout() <<" true." << std::endl;
      else                       Xyce::dout() <<" false." << std::endl;
      Xyce::dout() << "  Pflag =";
      if (iterV->neumannBCFlagP) Xyce::dout() <<" true." << std::endl;
      else                       Xyce::dout() <<" false." << std::endl;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::doAllocations
// Purpose       : A whole bunch of resizes.  Should be called after
//                 setupMesh.  Should only be called once.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/15/03
//-----------------------------------------------------------------------------
bool Instance::doAllocations ()
{
  bool bsuccess = true;

  // allocate conductance, capacitance array:
  condVec.resize(numElectrodes);
  capVec.resize(numElectrodes);
  for (int iE=0;iE<numElectrodes;++iE)
  {
    condVec[iE].resize(numElectrodes,0.0);
    capVec[iE].resize(numElectrodes,0.0);
  }

  // Set up a bunch of mesh-based arrays:
  // Local allocations:
  xVec.resize    (numMeshPoints);  meshContainerPtr->getXVector(xVec);
  yVec.resize    (numMeshPoints);  meshContainerPtr->getYVector(yVec);
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    for (int i=0;i<numMeshPoints;++i)
    {
      Xyce::dout() << " x["<<i<<"] = " << xVec[i];
      Xyce::dout() << " y["<<i<<"] = " << yVec[i];
      Xyce::dout() << std::endl;
    }
  }

  minDXVec.resize (numMeshPoints);

  areaVec.resize (numMeshPoints);

  VVec.resize      (numMeshPoints);
  RVec.resize      (numMeshPoints);
  SVec.resize      (numMeshPoints);
  totSrcVec.resize (numMeshPoints);
  nnVec.resize     (numMeshPoints);
  npVec.resize     (numMeshPoints);
  CVec.resize      (numMeshPoints);
  dRdpVec.resize   (numMeshPoints);
  dRdnVec.resize   (numMeshPoints);

  unVec.resize  (numMeshPoints, 0.0);
  upVec.resize  (numMeshPoints, 0.0);
  unE_Vec.resize (numMeshEdges, 0.0);
  upE_Vec.resize (numMeshEdges, 0.0);
  tnVec.resize  (numMeshPoints, 0.0);
  tpVec.resize  (numMeshPoints, 0.0);

  displPotential.resize(numMeshPoints);

  if (sgplotLevel > 0) outputVec.resize(numMeshPoints,0.0);

  //stateDispl.resize(numMeshPoints);
  li_stateDispl.resize(numMeshPoints,0);

  EfieldVec.resize (numMeshEdges);
  JnVec.resize (numMeshEdges);
  JpVec.resize (numMeshEdges);
  displCurrent.resize(numMeshEdges);

  dJndn1Vec.resize (numMeshEdges);
  dJndn2Vec.resize (numMeshEdges);
  dJndV1Vec.resize (numMeshEdges);
  dJndV2Vec.resize (numMeshEdges);

  dJpdn1Vec.resize (numMeshEdges);
  dJpdn2Vec.resize (numMeshEdges);
  dJpdV1Vec.resize (numMeshEdges);
  dJpdV2Vec.resize (numMeshEdges);

  //Vrowarray.resize (numMeshPoints,-1);  
  //Nrowarray.resize (numMeshPoints,-1);  
  //Prowarray.resize (numMeshPoints,-1);  

  boundarySten.resize(numMeshPoints,0);
  boundaryStenV.resize(numMeshPoints,0);
  boundaryStenN.resize(numMeshPoints,0);
  boundaryStenP.resize(numMeshPoints,0);
  boundaryTest.resize(numMeshPoints,0);

  li_Vrowarray.resize (numMeshPoints,0);
  li_Nrowarray.resize (numMeshPoints,0);
  li_Prowarray.resize (numMeshPoints,0);

  labelIndex.resize (numMeshPoints,0);
  labelNameVector.reserve(numMeshPoints);

  // allocate the boundary condition arrays:
  std::vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterV;

  for (iterV=first;iterV!=last; ++iterV)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
    int numPoints = labelPtr->mNodeVector.size();
    iterV->numBoundaryPoints = numPoints;
    iterV->VequVec.resize (numPoints,0.0);
    iterV->VbcVec.resize  (numPoints,0.0);
    iterV->nnbcVec.resize (numPoints,0.0);
    iterV->npbcVec.resize (numPoints,0.0);

    for (int i=0;i<iterV->numBoundaryPoints;++i)
    {
      mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
      int nodeIndex = labelPtr->mNodeVector[i];
      iterV->meshGlobalToLocal[nodeIndex] = i;
    }
  }

  // setup the aiEdge vector for tecplot:
  aiEdge.resize(numMeshEdges,0);
  aiEdge_nf.resize(numMeshEdges,0);

  // determine the label index for noflux.
  int nofluxIndex = 0;
  int numLabels = meshContainerPtr->getNumLabels();
  for (int i1=0;i1<numLabels;++i1)
  {
    mLabel * lPtr = meshContainerPtr->getLabel(i1);
    if (lPtr->name == "NOFLUX")
    {
      nofluxIndex = i1;
    }
  }

  // finish the aiEdge stuff.
  UINT iE;
  iNumPlotEdges = 0;
  iNumPlotEdges_nf = 0;
  for(iE = 0;iE<numMeshEdges;++iE)
  {
    mEdge * edgePtr = meshContainerPtr->getEdge(iE);
    UINT uLabel = edgePtr->uLabel;

    mLabel * labelPtr = meshContainerPtr->getLabel(uLabel);

    if(uLabel != -1u)
    {
      if(uLabel != nofluxIndex)
      {
        if(labelPtr->uType == TYPE_EDGE)
        {
          aiEdge[iE] = 1;
          ++iNumPlotEdges;
        }
      }
      else
      {
        if(labelPtr->uType == TYPE_EDGE)
        {
          aiEdge_nf[iE] = 1;
          ++iNumPlotEdges_nf;
        }
      }
    } // uLabel != -1u
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupLabelIndex
// Purpose       : This function sets up the vector, labelIndex
//
// Special Notes : labelIndex is a key, which maps a node index
//                 to a label index.  This is needed information, so that the
//                 load routines can determine what form rhs and
//                 Jacobian entries they should sum for a given node. - Should
//                 they set it up as a interior point, or a boundary, for
//                 example.
//
//                 By default, all nodes are assumed to belong to region label.
//                 However, some nodes will also belong to
//                 edges, or other regions, or whatever.  Sometimes a node
//                 will be on the node list of more than one region/edge.
//
//                 To resolve such descrepancies, this function assumes that
//                 an edge label will take priority over a region label.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/29/02
//-----------------------------------------------------------------------------
bool Instance::setupLabelIndex ()
{

  // first, loop over the various non-edge(region) labels, and loop over their
  // node lists.  Set the labelIndex for each of these nodes to correspond
  // the current region label.

  int i;
  std::vector<int>::iterator firstL;
  std::vector<int>::iterator lastL;
  std::vector<int>::iterator iterL;

  for (i=0;i<numMeshLabels;++i)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(i);
    if (labelPtr->uType == TYPE_EDGE) continue;

    firstL = labelPtr->mNodeVector.begin();
    lastL  = labelPtr->mNodeVector.end ();

    for (iterL=firstL; iterL!=lastL; ++iterL)
    {
       labelIndex[*iterL] = i;
    }
  }

  // next, loop through the edge labels, and their node lists, and update
  // the various label indices.  If there are conflicts, the edge index will
  // just write over a previously established  region index.

  for (i=0;i<numMeshLabels;++i)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(i);
    if (labelPtr->uType != TYPE_EDGE) continue;

    firstL = labelPtr->mNodeVector.begin();
    lastL  = labelPtr->mNodeVector.end ();

    for (iterL=firstL; iterL!=lastL; ++iterL)
    {
       labelIndex[*iterL] = i;
    }
  }

  for (i=0;i<numMeshPoints;++i)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(labelIndex[i]);
    labelNameVector.push_back(labelPtr->name);
  }

  int size = dIVec.size();
  for (i = 0;i < size; ++i)
  {
    labelDIMap[dIVec[i].eName] = i;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "In the Intance::setupLabelIndex ";
    Xyce::dout() << "  name = "<< getName() << std::endl;

    std::vector<std::string>::iterator firstM = labelNameVector.begin();
    std::vector<std::string>::iterator lastM  = labelNameVector.end ();
    std::vector<std::string>::iterator iterM;

    for (i=0, iterM=firstM;iterM!=lastM;++i,++iterM)
    {
      Xyce::dout() << "  i = "<<i<<" name = " << *iterM << std::endl;
    }
    std::map<std::string,int>::iterator firstDIM = labelDIMap.begin ();
    std::map<std::string,int>::iterator lastDIM  = labelDIMap.end ();
    std::map<std::string,int>::iterator iterDIM;
    for(iterDIM=firstDIM;iterDIM!=lastDIM;++iterDIM)
    {
      Xyce::dout() << " DI index = "<< iterDIM->second;
      Xyce::dout() << "  name = " << iterDIM->first << std::endl;
    }
  }

  if (isActive(Diag::DEVICE_PARAMETERS))
  {
    for (i=0;i<numMeshPoints;++i)
    {
      mLabel * labelPtr = meshContainerPtr->getLabel(labelIndex[i]);
      Xyce::dout() << "  labelIndex["<<i<<"] = " << labelIndex[i];
      Xyce::dout() << "  name = " << labelPtr->name << std::endl;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupBoundaryStencil
// Purpose       : This function sets up the stencil vector for boundary
//                 nodes.  In this case, a "boundary node" is a node from
//                 the mesh which is part of one of the user-specified
//                 boundary conditions.  So, this would include an
//                 electrode, like "COLLECTOR", but would not include
//                 "NOFLUX".
//
//                 If node i is a boundary node, then the value of
//                 boundarySten[i] = 1.  Otherwise boundarySten[i] = 0.
//
//                 7/15/03.  Revised to handle mixed boundary conditions,
//                 and the 3 new boundary stencils (one for each variable:
//                 (V,N,P)).
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/07/02
//-----------------------------------------------------------------------------
bool Instance::setupBoundaryStencil ()
{
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In Instance::setupBoundaryStencil." << std::endl;
  }

  for (iterDI = firstDI; iterDI!=lastDI; ++iterDI)
  {
    // loop over the nodes of this device interface node,
    // If it is an edge label, not a region label, and it is associated
    // with a boundary condition.

     if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

     std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     std::vector<int>::iterator iterI;

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int nodeIndex = *iterI;

       if (!(iterDI->neumannBCFlagV)) boundaryStenV[nodeIndex] = 1;
       if (!(iterDI->neumannBCFlagN)) boundaryStenN[nodeIndex] = 1;
       if (!(iterDI->neumannBCFlagP)) boundaryStenP[nodeIndex] = 1;

       // if this BC is dirichlet for all the variables, then set the
       // stencil.
       if (!(iterDI->neumannBCFlagV) &&
           !(iterDI->neumannBCFlagN) && !(iterDI->neumannBCFlagP))
       {
         boundarySten[nodeIndex] = 1;
       }
     }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::checkForElectrodeOverlap
// Purpose       : The purpose of this function is to make sure that there
//                 are not any nodes associated with multiple electrodes.
//
//                 If there were, the boundary conditions wouldn't make sense.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/2004
//-----------------------------------------------------------------------------
bool Instance::checkForElectrodeOverlap ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In Instance::checkForElectrodeOverlap." << std::endl;
  }

  for (int iDI=0;iDI<dIVec.size();++iDI)
  {
    // loop over the nodes of this device interface node,
    // If it is an edge label, not a region label, and it is associated
    // with a boundary condition.

     if ( !( meshContainerPtr->labelEdgeType (dIVec[iDI].eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[iDI].eName);

     std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     std::vector<int>::iterator iterI;

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int nodeIndex = *iterI;

       if (boundaryTest[nodeIndex] != 0)
       {
        UserFatal(*this)
          << "Electrodes " << dIVec[iDI].eName << " and " 
          << dIVec[boundaryTest[nodeIndex]-1].eName << " overlap";
       }

       boundaryTest[nodeIndex] = iDI+1;
     }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupNumVars
// Purpose       : The purpose of this function is to set up a few
//                 integer variables, such as numIntVars and numExtVars.
//                 These numbers are used by the registerGID , etc.
//                 functions.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/15/03
//-----------------------------------------------------------------------------
bool Instance::setupNumVars ()
{
  bool bsuccess = true;

  numIntVars    = 3*numMeshPoints;   // check this also.

  numExtVars    = numElectrodes; // This is the number of external nodes,
                                 // so it could be just about anything.

  numStateVars  = numElectrodes + numMeshPoints;
  //  numMeshPoints is part of numStateVars for displacement current.

  maxColsPerRow = 20;     // check this out later... depends on
                          // the max. NN count. (*3)

  int totalDirchlet = 0;

  // For the new boundary conditions, reduce the size of the problem by
  // the number of mesh points along the boundary
  // (If all BC are dirichlet, then *3, for each equation).
  numInterfaceMeshPoints = 0;

  std::vector<DeviceInterfaceNode>::iterator first = dIVec.begin ();
  std::vector<DeviceInterfaceNode>::iterator last  = dIVec.end   ();
  std::vector<DeviceInterfaceNode>::iterator iterV;
  for (iterV=first;iterV!=last; ++iterV)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
    numInterfaceMeshPoints += labelPtr->mNodeVector.size();
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << iterV->eName;
      Xyce::dout() << ":  numInterfaceMeshPoints = ";
      Xyce::dout() << labelPtr->mNodeVector.size ();
      Xyce::dout() << std::endl;
    }
  }

  // revising, because of possibility of mixed BC.
  for (iterV=first;iterV!=last; ++iterV)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterV->eName);
    int numPoints = labelPtr->mNodeVector.size();

    int mult = 0;

    if (!(iterV->neumannBCFlagV)) mult += 1;
    if (!(iterV->neumannBCFlagN)) mult += 1;
    if (!(iterV->neumannBCFlagP)) mult += 1;

    totalDirchlet += (mult * numPoints);
  }

  numIntVars -= totalDirchlet;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n";
    Xyce::dout() << " numInterfaceMeshPoints   = " << numInterfaceMeshPoints<< std::endl;
    Xyce::dout() << " numMeshPoints            = " << numMeshPoints<< std::endl;
    Xyce::dout() << " numElectrodes            = " << numElectrodes<< std::endl;
    Xyce::dout() << " numIntVars               = " << numIntVars<< std::endl;
    Xyce::dout() << " 3*numMeshPoints          = " << 3*numMeshPoints<< std::endl;
    Xyce::dout() << " 3*numInterfaceMeshPoints = "<<3*numInterfaceMeshPoints<< std::endl;
    Xyce::dout() << " totalDirchlet            = " << totalDirchlet<< std::endl;
    Xyce::dout() << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::allocatePDTerms.
//
// Purpose       : This function sets up and allocates a number of arrays
//                 that are needed by the function pdTerminalCurrents.
//
// Special Notes : Ordinarily, this function would have been called
//                 earlier, as I prefer to get allocations out of the way
//                 as early as possible.  
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/10/02
//-----------------------------------------------------------------------------
bool Instance::allocatePDTerms ()
{
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;
  // now do dIdX.
  for(iterDI=firstDI;iterDI!=lastDI;++iterDI)
  {
    mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

    std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
    std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
    std::vector<int>::iterator iterI;

    std::vector<EDGEINFO>::iterator firstEI;
    std::vector<EDGEINFO>::iterator lastEI;
    std::vector<EDGEINFO>::iterator iterEI;

    int cnt2  = 0;
    int nodeIndex;
    int col1;
    bool bmatch;
    int iVcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // voltage variables:
      // do the center point first.
      col1 = iterDI->Vcol[iVcol];
      if (col1 != -1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->dIdXcols.push_back(col1);
          iterDI->dIdX.push_back(0.0);
          iterDI->dQdX.push_back(0.0);
        }
      }
      ++iVcol;
      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iVcol)
      {
        col1 = iterDI->Vcol[iVcol];
        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->dIdXcols.push_back(col1);
            iterDI->dIdX.push_back(0.0);
            iterDI->dQdX.push_back(0.0);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    int iNcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // electron variables:
      // do the center point first.
      col1 = iterDI->Ncol[iNcol];
      if (col1 != -1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->dIdXcols.push_back(col1);
          iterDI->dIdX.push_back(0.0);
          iterDI->dQdX.push_back(0.0);
        }
      }
      ++iNcol;
      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iNcol)
      {
        col1 = iterDI->Ncol[iNcol];
        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->dIdXcols.push_back(col1);
            iterDI->dIdX.push_back(0.0);
            iterDI->dQdX.push_back(0.0);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    int iPcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // hole variables:
      // do the center point first.
      col1 = iterDI->Pcol[iPcol];
      if (col1 != -1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
        {
          if (iterDI->dIdXcols[cnt2] == col1)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->dIdXcols.push_back(col1);
          iterDI->dIdX.push_back(0.0);
          iterDI->dQdX.push_back(0.0);
        }
      }
      ++iPcol;
      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iPcol)
      {
        col1 = iterDI->Pcol[iPcol];
        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->dIdXcols.size();++cnt2)
          {
            if (iterDI->dIdXcols[cnt2] == col1)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->dIdXcols.push_back(col1);
            iterDI->dIdX.push_back(0.0);
            iterDI->dQdX.push_back(0.0);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    // Now add to the neighbor node array.  Assuming that any of the 3
    // columns iPcol, iNcol, iVcol will do for a validity check vs. -1,
    // so just checking the pcol.
    iPcol = 0;
    for(nodeIndex=0,iterI=firstI;iterI!=lastI;++iterI,++nodeIndex)
    {
      mNode * nodePtr = meshContainerPtr->getNode(*iterI);
      firstEI = nodePtr->edgeInfoVector.begin();
      lastEI  = nodePtr->edgeInfoVector.end ();

      // center point first.
      col1 = iterDI->Pcol[iPcol];
      int meshNode = *iterI;

      if (col1 !=-1)
      {
        // check if node == any previous nodes in the cols array
        bmatch = false;
        for (cnt2=0;cnt2<iterDI->neighborNodes.size();++cnt2)
        {
          if (iterDI->neighborNodes[cnt2] == meshNode)
          { bmatch = true; break; }
        }
        if (!bmatch)
        {
          iterDI->neighborNodes.push_back(meshNode);
        }
      }

      ++iPcol;

      // loop over the edges connected to the current node,
      // and do the neighbor point dependencies.
      for (iterEI=firstEI; iterEI!=lastEI; ++iterEI,++iPcol)
      {
        col1 = iterDI->Pcol[iPcol];
        int meshNode = iterEI->inode;

        if (col1 !=-1)
        {
          // check if node == any previous nodes in the cols array
          bmatch = false;
          for (cnt2=0;cnt2<iterDI->neighborNodes.size();++cnt2)
          {
            if (iterDI->neighborNodes[cnt2] == meshNode)
            { bmatch = true; break; }
          }
          if (!bmatch)
          {
            iterDI->neighborNodes.push_back(meshNode);
          }
        }
      } // end of nn edge loop
    } // end of node loop

    int size1 = iterDI->neighborNodes.size();
    int size3 = iterDI->dIdX.size();
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << std::endl;
      Xyce::dout() << "number of neighbor nodes for " << iterDI->eName;
      Xyce::dout() << " is " << size1 << std::endl;
      int  i;
      for (i=0;i<size1;++i)
      {
        Xyce::dout() << "neighborNodes["<<i<<"] = " << iterDI->neighborNodes[i] << std::endl;
      }
      Xyce::dout() << std::endl;
      Xyce::dout() << "dIdX size for " << iterDI->eName << "  is " << size3 << std::endl;
      for (i=0;i<size3;++i)
      {
        Xyce::dout() << "dIdX["<<i<<"] = " << iterDI->dIdXcols[i] << std::endl;
      }
    }

    // to set up dFdVckt, need to take  it through the same set of loops it
    // will be subject to in function pdTerminalCurrents.
    int numNeighbor = iterDI->neighborNodes.size();
    int iNeighbor;
    int dFdVindex = 0;
    for (iNeighbor=0;iNeighbor<numNeighbor;++iNeighbor)
    {
      int inode = iterDI->neighborNodes[iNeighbor];
      mNode * nodePtr = meshContainerPtr->getNode(inode);

      for (int iNN=0;iNN<nodePtr->cnode;++iNN)
      {
        int inodeB = nodePtr->edgeInfoVector[iNN].inode;

        // if nodeB is not a boundary node, never mind.
        if (boundaryStenV[inodeB]!=1) continue;

        // if it is a boundary node, but not part of the
        // current boundary, also never mind.
        if (labelNameVector[inodeB]!= iterDI->eName) continue;

        iterDI->dFdVckt.push_back(0.0); ++dFdVindex;
        iterDI->dFdVckt.push_back(0.0); ++dFdVindex;
        iterDI->dFdVckt.push_back(0.0); ++dFdVindex;
      }
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << iterDI->eName << ": size of dFdVckt = " << dFdVindex << std::endl;
    }
  }  // end of DI loop

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupBCEdgeAreas
//
// Purpose       : This function sets up the areas associated with each
//                 mesh node along a boundary condition edge.
//
// Special Notes : This is potentially tricky.  The mesh may be 2D
//                 cartesian or 2D cylindrical.
//
//                 Fortunately, the function "lengthAdjust", from the mesh
//                 class, is designed for getting areas - either
//                 cylindrical or cartesian.
//
//                 This whole should be moved to the mesh class later.
//
//                 This function is called before scaling is turned on, but
//                 should work either way.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/09/02
//-----------------------------------------------------------------------------
bool Instance::setupBCEdgeAreas ()
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "setupBCEdgeAreas.  name = " << getName() << std::endl;
    Xyce::dout().setf(std::ios::scientific);
  }

  // now set up the local areas for needed to interface device boundary
  // conditions to the circuit

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  for (iterDI = firstDI; iterDI!=lastDI; ++iterDI)
  {
    // loop over the nodes of this device interface node:

     if ( !( meshContainerPtr->labelEdgeType (iterDI->eName) ) ) continue;

     mLabel * labelPtr = meshContainerPtr->getLabel(iterDI->eName);

     std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     std::vector<int>::iterator iterI;

     iterDI->area       = 0.0;  // total area for the edge.

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // loop over neighbor nodes/edges to get area sum for this node.
       mNode * nodePtr = meshContainerPtr->getNode(*iterI);

       std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
       std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end  ();
       std::vector<EDGEINFO>::iterator iterEI;

       double areaLocal = 0.0; // total "area" for the this node
       double areaTmp   = 0.0; // partial area for the this node(from one edge)

       if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
       {
         Xyce::dout() << " --------------- " << std::endl;
         Xyce::dout() << "name = " << iterDI->eName;
         Xyce::dout() << "  node      = " << *iterI << std::endl;
       }

       for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
       {
         int neighbor = iterEI->inode;

         // if this edge is actually on the boundary, then sum the
         // edge length into the "area" associated with this
         // boundary node.  Check this by checking the label index
         // of the neighbor.

         int ilabel = labelIndex[neighbor];
         mLabel * labelPtr = meshContainerPtr->getLabel(ilabel);

         areaTmp = 0.0;
         if (labelPtr->name == iterDI->eName)
         // if this is along the bounadry
         {
           if (meshContainerPtr->cylGeom)
           {
             // get the midpoint location, x2:
             double x1 = xVec[*iterI];
             double y1 = yVec[*iterI];

                   double x2 = xVec[neighbor];
                   double y2 = yVec[neighbor];

             double dx = x2-x1;   dx *= 0.5;
             double dy = y2-y1;   dy *= 0.5;

             x2 = x1 + dx;
             y2 = y1 + dy;

             areaTmp = meshContainerPtr->lengthAdjust(x1,y1,x2,y2);
             areaLocal += areaTmp;
           }
           else // cartesian
           {
             areaTmp    = 0.5 * iterEI->elen;
             areaLocal += 0.5 * iterEI->elen;
           }
         }
         if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
         {
           Xyce::dout() << "  neighbor node   = " << neighbor << std::endl;
           Xyce::dout() << "  areaTmp         = " << areaTmp << std::endl;
           Xyce::dout() << "  areaLocal       = " << areaLocal << std::endl;
           Xyce::dout() << "  elen            = " << iterEI->elen << std::endl;
           Xyce::dout() << "  label name      = " << labelPtr->name << std::endl;
           Xyce::dout() << "  DI eName        = " << iterDI->eName << std::endl;
           Xyce::dout() << "---" << std::endl;
         }
       }

       if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
       {
         Xyce::dout() << " --------------- " << std::endl;
       }
       iterDI->area += areaLocal;
       iterDI->areaVector.push_back(areaLocal);

     } // iterI loop.


     if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout().setf(std::ios::scientific);
      Xyce::dout() << "  Total area for edge: " << iterDI->area << std::endl;
    }
  }  // iterDI loop.

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << section_divider << "\n";
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMinDXVector
// Purpose       : This finds the minimum edge length.  I'm not sure if
//                 there's any reason to do this anymore.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/26/02
//-----------------------------------------------------------------------------
bool Instance::setupMinDXVector ()
{

  int i;
  for (i=0;i<numMeshPoints;++i)
  {
    // loop over neighbor nodes/edges to get area sum for this node.
    mNode * nodePtr = meshContainerPtr->getNode(i);

    std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
    std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end  ();
    std::vector<EDGEINFO>::iterator iterEI;

    double tmpDX = +1.0e99;
    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
      if (tmpDX > iterEI->elen) tmpDX = iterEI->elen;
    }
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "i = " << i << "   minDX = " << tmpDX; Xyce::dout() << std::endl;
    }
    minDXVec[i] = tmpDX;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp
//
// Purpose       : This function sets up the jacobian stamp data structure,
//                 which is a symbolic, local-ID (LID) representation of the
//                 jacobian matrix structure of the device.
//
//                 The jacStamp is the structure that gets passed up to
//                 topology.
//
//                 It is similar to some of what has to happen in the
//                 registerGID's function, but it is not dependent upon
//                 knowing the global ID's for the solution variables.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/23/03
//-----------------------------------------------------------------------------
bool Instance::setupJacStamp ()
{
  bool bsuccess = true;
  int i;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In Instance::setupJacStamp" << std::endl;
    Xyce::dout() << "numIntVars = " << numIntVars << std::endl;
    Xyce::dout() << "numMeshPoints = " << numMeshPoints << std::endl;
  }

  // The number of rows in the jacStamp has to include all the *possible*
  // electrodes (external variables). This is a little confusing.
  // The PDE devices are set up so that there are:
  //   2 required nodes
  //   2 not-required, but not-optional nodes (these are fill nodes)
  //   100 optional nodes.
  //   For some reason (that I don't understand) this means that the
  //   jacStamp always has to have at least 4 external nodes.  Even
  //   if the device only has 2 nodes attached to the circuit, the
  //   jacStamp is required to have 4.
  int numPossibleElectrodes = (numExtVars>4)?numExtVars:4;
  jacStamp.resize(numIntVars + numPossibleElectrodes);

  MESHtoLID_V.resize(numMeshPoints,-1);
  MESHtoLID_N.resize(numMeshPoints,-1);
  MESHtoLID_P.resize(numMeshPoints,-1);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "size of jacStamp = " << jacStamp.size() << std::endl;
  }

  // Set up the MESHtoLID converter for internal vars.
  int LIDcounter = numPossibleElectrodes;

  for (i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    if (!(boundaryStenV[i]))
    {
      MESHtoLID_V[i] = LIDcounter; ++LIDcounter;
    }

    if (!(boundaryStenN[i]))
    {
      MESHtoLID_N[i] = LIDcounter; ++LIDcounter;
    }

    if (!(boundaryStenP[i]))
    {
      MESHtoLID_P[i] = LIDcounter; ++LIDcounter;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // Do the external variables first:
  // Loop over the boundary condition "device interface" labels.
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;
  int DIsize = dIVec.size ();
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  // note that the index, in addition to being the index into the
  // array of "DeviceInterfaceNodes" will also be the LID.
  int index;
  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    // check that this label exists, and is an edge label
    // (everything in the dIVec container should pass these tests
    //   by this point - these two if-statements are just a failsafe.)
     if ( !( meshContainerPtr->labelNameExist(dIVec[index].eName) ) ) continue;
     if ( !( meshContainerPtr->labelEdgeType (dIVec[index].eName) ) ) continue;


     // first do the external node's (row=cktnode,col=cktnode) pair:
     jacStamp[index].push_back(index);

     dIVec[index].numCrossTerms = 0;

     // next do the (row,col) pairs for other ckt nodes.  This is needed for
     // 2-level newton only.
     for (int ind2=0;ind2<DIsize;++ind2)
     {
       if (ind2==index) continue;

       jacStamp[index].push_back(ind2);
       ++(dIVec[index].numCrossTerms);
     }

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);

     std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     std::vector<int>::iterator iterI;

     std::vector<EDGEINFO>::iterator firstEI;
     std::vector<EDGEINFO>::iterator lastEI;
     std::vector<EDGEINFO>::iterator iterEI;

     int imesh;

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // now do the V edge and edge neighbor nodes.
       imesh = *iterI;
       if (!boundaryStenV[imesh])
         jacStamp[index].push_back(MESHtoLID_V[imesh]);

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         imesh = iterEI->inode;
         if (!boundaryStenV[imesh])
           jacStamp[index].push_back(MESHtoLID_V[imesh]);
       }
     }

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // now do the N nodes.
       imesh = *iterI;

       if (!boundaryStenN[imesh])
         jacStamp[index].push_back(MESHtoLID_N[imesh]);

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();
       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         imesh = iterEI->inode;
         if (!boundaryStenN[imesh])
           jacStamp[index].push_back(MESHtoLID_N[imesh]);
       }
     }

     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       // now do the P nodes.
       imesh = *iterI;
       if (!boundaryStenP[imesh])
         jacStamp[index].push_back(MESHtoLID_P[imesh]);

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       firstEI = nodePtr->edgeInfoVector.begin();
       lastEI  = nodePtr->edgeInfoVector.end ();
       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         imesh = iterEI->inode;
         if (!boundaryStenP[imesh])
           jacStamp[index].push_back(MESHtoLID_P[imesh]);
       }
     }// iterI loop
  } // index (dIVec) loop.

  /////////////////////////////////////////////////////////////////////////
  // Now do the internal variables.  (V,N,P on the mesh)
  for(i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "  mesh point i = " << i << std::endl;
    }

    mNode * nodePtr = meshContainerPtr->getNode(i);
    std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin ();
    std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end   ();
    std::vector<EDGEINFO>::iterator iterEI;

    // get the temporary LID row indices:
    int Vrow = MESHtoLID_V[i];
    int Nrow = MESHtoLID_N[i];
    int Prow = MESHtoLID_P[i];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "  Vrow = " << Vrow << std::endl;
      Xyce::dout() << "  Nrow = " << Nrow << std::endl;
      Xyce::dout() << "  Prow = " << Prow << std::endl;
    }

    // voltage row:
    if (Vrow != -1)
    {
      // First do the (row, row) pair.
      jacStamp[Vrow].push_back(Vrow);

      // loop  over the neighbor nodes of node i.
      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        // if the neighbor node is on a boundary, need to use the
        // GID of the connected ckt node, rather than the (nonexistant)
        // GID of the boundary mesh node.
        if (boundaryStenV[iterEI->inode])
        {
          // get the id:
          int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
          jacStamp[Vrow].push_back (DIindex);
        }
        else
        {
          int imesh = iterEI->inode;
          int lid = MESHtoLID_V[imesh];
          jacStamp[Vrow].push_back(lid);
        }
      }

      jacStamp[Vrow].push_back(Nrow);
      jacStamp[Vrow].push_back(Prow);
    }

    // electron row:
    if (Nrow != -1)
    {
      jacStamp[Nrow].push_back(Nrow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        int imesh = iterEI->inode;
        if (!boundaryStenN[imesh])
        {
          int lid = MESHtoLID_N[imesh];
          jacStamp[Nrow].push_back(lid);
        }
      }

      jacStamp[Nrow].push_back(Vrow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        // if the neighbor node is on a boundary, need to use the
        // GID of the connected ckt node, rather than the (nonexistant)
        // GID of the boundary mesh node.
        if (boundaryStenV[iterEI->inode])
        {
          // get the id:
          int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
          jacStamp[Nrow].push_back (DIindex);
        }
        else
        {
          int imesh = iterEI->inode;
          int lid = MESHtoLID_V[imesh];
          jacStamp[Nrow].push_back(lid);
        }
      }
      jacStamp[Nrow].push_back(Prow);
    }

    // hole row:
    if (Prow != -1)
    {
      jacStamp[Prow].push_back(Prow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        int imesh = iterEI->inode;
        if(!boundaryStenP[imesh])
        {
          int lid = MESHtoLID_P[imesh];
          jacStamp[Prow].push_back(lid);
        }
      }
      jacStamp[Prow].push_back(Vrow);

      for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
      {
        // if the neighbor node is on a boundary, need to use the
        // GID of the connected ckt node, rather than the (nonexistant)
        // GID of the boundary mesh node.
        if (boundaryStenV[iterEI->inode])
        {
          // get the id:
          int DIindex = labelDIMap[labelNameVector[iterEI->inode]];
          jacStamp[Prow].push_back (DIindex);
        }
        else
        {
          int imesh = iterEI->inode;
          int lid = MESHtoLID_V[imesh];
          jacStamp[Prow].push_back(lid);
        }
      }
      jacStamp[Prow].push_back(Nrow);
    }
  } // mesh point loop.

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    int irow, icol;
    for(irow=0;irow<jacStamp.size();++irow)
    {
      Xyce::dout() << "irow = " << irow;
      if (irow < dIVec.size())
        Xyce::dout() << "  " << dIVec[irow].eName << "  KCL";
      Xyce::dout() << std::endl;
      for (icol=0;icol<jacStamp[irow].size();++icol)
      {
        int jsTmp = jacStamp[irow][icol];
        Xyce::dout() << "   jacStamp["<<irow<<"]["<<icol<<"] = "<<jsTmp;
        Xyce::dout() << std::endl;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  // set up the intNameMap, if necessary.
  for (int i = 0; i < numMeshPoints; ++i)
  {
    if (boundarySten[i]) continue;

    int Vrow, Nrow, Prow;
    Vrow = li_Vrowarray[i];
    Nrow = li_Nrowarray[i];
    Prow = li_Prowarray[i];

    if (Vrow != -1)
    {
      std::ostringstream oss;
      oss << "_V_" << i << "_" << labelNameVector[i];
      addInternalNode(symbol_table, Vrow, getName(), oss.str());
    }

    if (Nrow != -1)
    {
      std::ostringstream oss;
      oss << "_N_" << i << "_" << labelNameVector[i];
      addInternalNode(symbol_table, Nrow, getName(), oss.str());
    }

    if (Prow != -1)
    {
      std::ostringstream oss;
      oss << "_P_" << i << "_" << labelNameVector[i];
      addInternalNode(symbol_table, Prow, getName(), oss.str());
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupRowColPairsLID
//
// Purpose       : This function performs part of the "registerLIDs"
//                 functionality, in that it sets up the Jacobian
//                 matrix (row,col) pairs.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 
//-----------------------------------------------------------------------------
void Instance::setupRowColPairsLID ()
{
  int i, j;
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "LID: Instance::setupRowColPairsLID ()" << std::endl;
    Xyce::dout() << "LID: doing some boundary condition vars:"<< std::endl;
  }

  // go back to the nodes which are boundary condition nodes - ones
  // that are connected to the external nodes, and add a few
  // (row, col) pairs .  These nodes will need some extra ones, at
  // least to  handle boundary conditions on the voltage.

  // loop over the boundary condition "device interface" labels.
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  int index;

  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    // check that this label exists, and is an edge label
    // (everything in the dIVec container should pass these tests
    //   by this point)
     if ( !( meshContainerPtr->labelNameExist(dIVec[index].eName) ) ) continue;
     if ( !( meshContainerPtr->labelEdgeType (dIVec[index].eName) ) ) continue;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "LID: Device Interface: " << dIVec[index].eName << std::endl;
    }

     // These are for satisfying the KCL of the node - it will depend on
     // all the currents of all the edge nodes.  Unfortunately, this
     // results in a potentially(probably) dense matrix row.
     //
     // Note that if this is a dielectric boundary, then no current is
     // going to come out of device at this boundary (no conduction
     // current, anyway), and that a lot of these (row, col) pairs will
     // wind up being loaded with 0's in that case.  However, it doesn't
     // hurt to have them, so I'm leaving them in.

     // first push back the (lid,lid) entry.  Then do the entries related
     // to the  nearest neighbors of edge nodes.
     //dIVec[index].col.push_back(dIVec[index].lid);

     if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
     {
       Xyce::dout() << "LID: V edge and edge neighbor ids:" << std::endl;
     }

     mLabel * labelPtr = meshContainerPtr->getLabel(dIVec[index].eName);
     std::vector<int>::iterator firstI = labelPtr->mNodeVector.begin();
     std::vector<int>::iterator lastI  = labelPtr->mNodeVector.end  ();
     std::vector<int>::iterator iterI;
     int itmp;

     // now do the V edge and edge neighbor nodes.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int inode=*iterI;
       if (boundarySten[inode]) { itmp=-1; } // ERK.this is concerning. Does this ever happen?
       else { itmp = li_Vrowarray[*iterI]; }
       dIVec[index].Vcol.push_back(itmp);

       if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
       {
         int ind1 = dIVec[index].Vcol.size()-1;
         Xyce::dout() << "LID: node="<<inode<<"  1Vcol["<<ind1<<"] = " << itmp << std::endl;
       }

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
       std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end ();
       std::vector<EDGEINFO>::iterator iterEI;

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         int inodeB=iterEI->inode;
         if (boundarySten[inodeB]) { itmp=-1; } // ERK.this is concerning. Does this ever happen?
         else { itmp = li_Vrowarray[iterEI->inode]; }
         dIVec[index].Vcol.push_back(itmp);

         if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
         {
           int ind1 = dIVec[index].Vcol.size()-1;
           Xyce::dout() << "LID: node="<<inodeB<<"  2Vcol["<<ind1<<"] = " << itmp << std::endl;
         }
       }
     } // iterI loop

     if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
     {
       Xyce::dout() << "LID: N edge and edge neighbor ids:" << std::endl;
     }

     // now do the N nodes.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int inode=*iterI;
       if (boundarySten[inode]) { itmp=-1; }
       else { itmp = li_Nrowarray[*iterI]; }
       dIVec[index].Ncol.push_back(itmp);

       if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
       {
         int ind1 = dIVec[index].Ncol.size()-1;
         Xyce::dout() << "LID: node="<<inode<<" 1Ncol["<<ind1<<"] = " << itmp << std::endl;
       }

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
       std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end ();
       std::vector<EDGEINFO>::iterator iterEI;

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         int inodeB=iterEI->inode;
         if (boundarySten[inodeB]) { itmp=-1; }
         else { itmp = li_Nrowarray[iterEI->inode]; }
         dIVec[index].Ncol.push_back(itmp);

         if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
         {
           int ind1 = dIVec[index].Ncol.size()-1;
           Xyce::dout() << "LID: node="<<inodeB<<" 2Ncol["<<ind1<<"] = " << itmp << std::endl;
         }
       }
     } // iterI loop

     if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
     {
       Xyce::dout() << "LID: P edge and edge neighbor ids:" << std::endl;
     }
     // now do the P nodes.
     for(iterI=firstI;iterI!=lastI;++iterI)
     {
       int inode=*iterI;
       if (boundarySten[inode]) { itmp=-1; }
       else { itmp = li_Prowarray[*iterI]; }
       dIVec[index].Pcol.push_back(itmp);

       if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
         {
           int ind1 = dIVec[index].Pcol.size()-1;
           Xyce::dout() << "LID: node="<<inode<<" 1Pcol["<<ind1<<"] = " << itmp << std::endl;
         }

       mNode * nodePtr = meshContainerPtr->getNode(*iterI);
       std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin();
       std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end ();
       std::vector<EDGEINFO>::iterator iterEI;

       // loop over the edges connected to the current node:
       for (iterEI=firstEI; iterEI!=lastEI; ++iterEI)
       {
         int inodeB=iterEI->inode;
         if (boundarySten[inodeB]) { itmp=-1; }
         else { itmp = li_Prowarray[iterEI->inode]; }
         dIVec[index].Pcol.push_back(itmp);

         if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
         {
           int ind1 = dIVec[index].Pcol.size()-1;
           Xyce::dout() << "LID: node="<<inodeB<<" 2Pcol["<<ind1<<"] = " << itmp << std::endl;
         }
       }
     }// iterI loop

     if (maxColsPerRow <
                (dIVec[index].Vcol.size() +
                 dIVec[index].Ncol.size() +
                 dIVec[index].Pcol.size() + 10)
        )
     {
         maxColsPerRow =
                (dIVec[index].Vcol.size() +
                 dIVec[index].Ncol.size() +
                 dIVec[index].Pcol.size() + 10); // extra 10 just in case.
     }
  }  // end of index (dIVec) loop
}

//-----------------------------------------------------------------------------
// Function      : 2DDPEInstance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter,  SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                        const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  // AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "In the Intance::registerLIDs function.  ";
    Xyce::dout() << "  name = "<< getName() << std::endl;
    Xyce::dout() << "numInt = " << numIntVars << std::endl;
    Xyce::dout() << "numExt = " << numExtVars << std::endl;
    Xyce::dout() << "numMeshPoints = " << numMeshPoints << std::endl;
  }

  // copy over the global ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // First do the external variables:
  // These will all be voltage nodes connected to the devices.
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Setting up the indices for the external circuit nodes:" << std::endl;
    Xyce::dout() << "External node list:" << std::endl;
  }

  int isizeDI = dIVec.size();

  int index;
  for(index=0; index < isizeDI; ++index)
  {
     dIVec[index].lid = extLIDVec[index];

     if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
     {
       Xyce::dout() << "   name = "<<dIVec[index].eName<<" lid = ";
       Xyce::dout() << dIVec[index].lid;
       Xyce::dout() << std::endl;
     }
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Doing internal vars, row arrays:"<< std::endl;
  }

  // Do the internal variables.  There should be a lot of these.
  int i=0;
  index = 0;

  // The interior points will be blocked (V,N,P) together.
  while (i < numMeshPoints)
  {
    if (boundarySten[i]) { ++i; continue; }

    if (!(boundaryStenV[i]))
    {
      li_Vrowarray[i] = intLIDVec[index];   // should be Vrowarray[i]
      ++index;
    }

    if (!(boundaryStenN[i]))
    {
      li_Nrowarray[i] = intLIDVec[index];
      ++index;
    }

    if (!(boundaryStenP[i]))
    {
      li_Prowarray[i] = intLIDVec[index];
      ++index;
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "doing lid row arrays for mesh point " << i << std::endl;
      if (!(boundaryStenV[i])) Xyce::dout() << "  li_Vrow = " << li_Vrowarray[i] << std::endl;
      if (!(boundaryStenN[i])) Xyce::dout() << "  li_Nrow = " << li_Nrowarray[i] << std::endl;
      if (!(boundaryStenP[i])) Xyce::dout() << "  li_Prow = " << li_Prowarray[i] << std::endl;
    }

    ++i;
  }


  setupRowColPairsLID ();

  // Make sure the cols and vals arrays are big enough,
  // based on maxColsPerRow. (the variable maxColsPerRow was
  // setup in function setupRowColPairs).
  getMatrixLoadData().initializeAll(maxColsPerRow);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\n";
    Xyce::dout() << section_divider << "\n";
    Xyce::dout() << "  In Instance::registerStateLIDs\n\n";
    Xyce::dout() << "  name             = " << getName() << "\n";
    Xyce::dout() << "  Number of State LIDs: " << numStateVars << "\n";
  }

  // Copy over the local ID lists:
  staLIDVec = staLIDVecRef;

  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end ();
  std::vector<DeviceInterfaceNode>::iterator iterDI  = firstDI;

  int i=0,j=0;
  for (; iterDI!=lastDI;++iterDI,++i)
  {
    iterDI->li_stateC = staLIDVec[i];
  }

  for (j=0;j<numMeshPoints;++j,++i)
  {
     li_stateDispl[j] = staLIDVec[i];
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "  State indices:" << "\n";
    Xyce::dout() << "\n";
    for (iterDI=firstDI; iterDI!=lastDI;++iterDI)
    {
      Xyce::dout() << "  ";
      Xyce::dout().width(12);
      Xyce::dout().setf(std::ios::right);
      Xyce::dout() << iterDI->eName;
      Xyce::dout().setf(std::ios::left);
      Xyce::dout() << "  li_stateC = " << iterDI->li_stateC;
      Xyce::dout() << std::endl;
    }

    Xyce::dout() << "  Displacement state indices:\n";
    for (j=0;j<numMeshPoints;++j,++i)
    {
       Xyce::dout() << "  edge: " << j << "  li_stateDispl = " << li_stateDispl[j];
       Xyce::dout() << std::endl;
    }

    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/23/03
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
//
// Purpose       : This function sets up the "local-ID" (LID) jacobian
//                 stamp.  This is neccessary for the DMA=on capability,
//                 which is the default.  This represents the stuff that
//                 comes in from topology.
//
// Special Notes : This needs to be consistent with the function,
//                 Instance::setupJacStamp.
//
//
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233
// Creation Date : 02/23/03
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs
    ( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs ( jacLIDVec );

  int i;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In Instance::registerJacLIDs" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  // Do the external variables first:
  // Loop over the boundary condition "device interface" labels.
  std::vector<DeviceInterfaceNode>::iterator firstDI = dIVec.begin();
  std::vector<DeviceInterfaceNode>::iterator lastDI  = dIVec.end  ();
  std::vector<DeviceInterfaceNode>::iterator iterDI;
  int DIsize = dIVec.size ();
  firstDI = dIVec.begin();
  lastDI  = dIVec.end  ();

  int index;
  for(index=0,iterDI=firstDI;iterDI!=lastDI;++iterDI,++index)
  {
    int jacRowSize = jacLIDVec[index].size();
    dIVec[index].dIdXoffset.resize(jacRowSize-dIVec[index].numCrossTerms-1);
    dIVec[index].lidOffset = jacLIDVec[index][0];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "index = " << index;
      Xyce::dout() << "  jacRowSize = " << jacRowSize;
      Xyce::dout() << "  name = " << dIVec[index].eName << std::endl;
      Xyce::dout() << " lidOffset = ";
      Xyce::dout() << dIVec[index].lidOffset << std::endl;
    }

    int nCT = dIVec[index].numCrossTerms;

    dIVec[index].crossOffsets.resize(nCT);

    for(int itmp=0;itmp<nCT;++itmp)
    {
      dIVec[index].crossOffsets[itmp] = jacLIDVec[index][itmp+1];
      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "  crossOffsets["<<itmp<<"] = ";
        Xyce::dout() << dIVec[index].crossOffsets[itmp] << std::endl;
      }
    }

    for (int ioff=nCT+1;ioff<jacRowSize;++ioff)
    {
      int tmpIndex = ioff - (nCT+1);
      dIVec[index].dIdXoffset[tmpIndex] = jacLIDVec[index][ioff];

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << " dIdXoffset["<<tmpIndex<<"] = ";
        Xyce::dout() << dIVec[index].dIdXoffset[tmpIndex] << std::endl;
      }
    }
  } // index (dIVec) loop.

  /////////////////////////////////////////////////////////////////////////
  // Now do the internal variables.  (V,N,P on the mesh)
  li_VoffsetArray.resize(numMeshPoints);
  li_NoffsetArray.resize(numMeshPoints);
  li_PoffsetArray.resize(numMeshPoints);
  for(i=0;i<numMeshPoints;++i)
  {
    if (boundarySten[i]) continue;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << " mesh point i = " << i << std::endl;
    }

    mNode * nodePtr = meshContainerPtr->getNode(i);
    std::vector<EDGEINFO>::iterator firstEI = nodePtr->edgeInfoVector.begin ();
    std::vector<EDGEINFO>::iterator lastEI  = nodePtr->edgeInfoVector.end   ();

    // get the temporary LID row indices:
    int Vrow = MESHtoLID_V[i];
    int Nrow = MESHtoLID_N[i];
    int Prow = MESHtoLID_P[i];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "   Vrow = " << Vrow << std::endl;
      Xyce::dout() << "   Nrow = " << Nrow << std::endl;
      Xyce::dout() << "   Prow = " << Prow << std::endl;
    }

    int ioff;
    // voltage row:
    if (Vrow != -1)
    {
      int VrowSize = jacLIDVec[Vrow].size();
      li_VoffsetArray[i].resize(VrowSize);
      for (ioff=0;ioff<VrowSize;++ioff)
      {
        li_VoffsetArray[i][ioff] = jacLIDVec[Vrow][ioff];
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "   li_Voffset["<<i<<"]["<<ioff<<"] = ";
          Xyce::dout() << li_VoffsetArray[i][ioff] << std::endl;
        }
      }
    }

    // electron row:
    if (Nrow != -1)
    {
      int NrowSize = jacLIDVec[Nrow].size();
      li_NoffsetArray[i].resize(NrowSize);
      for (ioff=0;ioff<NrowSize;++ioff)
      {
        li_NoffsetArray[i][ioff] = jacLIDVec[Nrow][ioff];
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "   li_Noffset["<<i<<"]["<<ioff<<"] = ";
          Xyce::dout() << li_NoffsetArray[i][ioff] << std::endl;
        }
      }
    }

    // hole row:
    if (Prow != -1)
    {
      int ProwSize = jacLIDVec[Prow].size();
      li_PoffsetArray[i].resize(ProwSize);
      for (ioff=0;ioff<ProwSize;++ioff)
      {
        li_PoffsetArray[i][ioff] = jacLIDVec[Prow][ioff];
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << "   li_Poffset["<<i<<"]["<<ioff<<"] = ";
          Xyce::dout() << li_PoffsetArray[i][ioff] << std::endl;
        }
      }
    }

  } // mesh point loop
}

} // namespace TwoDPDE
} // namespace Device
} // namespace Xyce

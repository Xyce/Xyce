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
// Purpose        : Contains the 2D mesh information.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/21/02
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <cstring>
#include <cstdio>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_PDE_2DMesh.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>

namespace {
typedef unsigned int    UINT;
}

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::PDE_2DMesh
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
PDE_2DMesh::PDE_2DMesh (const DeviceOptions & do1, int sgplotLevel1):
    cylGeom   (false),
    meshFileName ("internal"),
    externalMeshFlag (false),
    xMax      (0.0),
    yMax      (0.0),
    xMin      (0.0),
    yMin      (0.0),
    dx        (0.0),
    dy        (0.0),
    xRatio    (1.0),
    yRatio    (1.0),
    x0        (1.0),
    meshScaledFlag (false),
    vol       (0.0),
    invVol    (0.0),
    surfArea  (0.0),
    circum    (0.0),
    invCircum (0.0),
    depth     (1.0),
    numAdj    (0),
    numNodes  (0),
    numEdges  (0),
    numCells  (0),
    numLabels (0),
    numRegLabels (1), // default is 1 - has to be at least one region.
    numBndryNodes(0),
    maxNodeNN (0),
    iRecentCellLookup(0),
    dopingSet (false),
    ixMax(0),
    iyMax(0),
    nodeIndices(NULL),
    edgeIndices(NULL),
    cellIndices(NULL),
    aiBegin(NULL),
    aiEnd(NULL),
    adjInfoAllocFlag(false),
    devOptions_(&do1),
    sgplotLevel(sgplotLevel1),
    useDefaultLabels(true)
{


}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::PDE_2DMesh
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
PDE_2DMesh::PDE_2DMesh (const PDE_2DMesh & right) :
    cylGeom          (right.cylGeom),
    meshFileName     (right.meshFileName),
    externalMeshFlag (right.externalMeshFlag),
    xMax             (right.xMax),
    yMax             (right.yMax),
    xMin             (right.xMin),
    yMin             (right.yMin),
    dx               (right.dx),
    dy               (right.dy),
    xRatio           (right.xRatio),
    yRatio           (right.yRatio),
    x0               (right.x0),
    meshScaledFlag   (right.meshScaledFlag),
    vol              (right.vol),
    invVol           (right.invVol),
    surfArea         (right.surfArea),
    circum           (right.circum),
    invCircum        (right.invCircum),
    depth            (right.depth),
    numNodes         (right.numNodes),
    numEdges         (right.numEdges),
    numCells         (right.numCells),
    numLabels        (right.numLabels),
    numRegLabels     (right.numRegLabels),
    numBndryNodes    (right.numBndryNodes),
    maxNodeNN        (right.maxNodeNN),
    iRecentCellLookup(right.iRecentCellLookup),
    dopingSet        (right.dopingSet),
    mNodeVector      (right.mNodeVector),
    mEdgeVector      (right.mEdgeVector),
    mCellVector      (right.mCellVector),
    mLabelVector     (right.mLabelVector),
    dopingVector     (right.dopingVector),
    xVector          (right.xVector),
    yVector          (right.yVector),
    visitCellFlagVec (right.visitCellFlagVec),
    mLabelMap        (right.mLabelMap),
    ixMax            (right.ixMax),
    iyMax            (right.iyMax),
    afVisitedVec     (right.afVisitedVec),
    aiBegin(NULL),
    aiEnd(NULL),
    adjInfoAllocFlag (right.adjInfoAllocFlag),
    devOptions_      (right.devOptions_),
    sgplotLevel      (right.sgplotLevel)
{
  int i,j;

  // if we're copying over a mesh that was generated via the "internal"
  // capability, we need to copy over these 3 structures:
  if (!externalMeshFlag)
  {
    if (right.nodeIndices != NULL)
    {
      nodeIndices = new int*[ixMax+10];
      for (i=0;i<ixMax+10;++i)
      {
        nodeIndices[i] = new int[iyMax+10];

        for (j=0;j<iyMax+10;++j)
        {
          nodeIndices[i][j] = right.nodeIndices[i][j];
        }
      }
    }

    if (right.edgeIndices != NULL)
    {
      edgeIndices = new int*[numNodes+10];
      for (i=0;i<numNodes+10;++i)
      {
        edgeIndices[i] = new int[numNodes+10];

        for (j=0;j<numNodes+10;++j)
        {
          edgeIndices[i][j] = right.edgeIndices[i][j];
        }
      }
    }

    if (right.cellIndices != NULL)
    {
      cellIndices = new int *[ixMax+10];
      for (i=0;i<ixMax+10;++i)
      {
        cellIndices[i] = new int[iyMax+10];
        for (j=0;j<iyMax+10;++j)
          cellIndices[i][j] = right.cellIndices[i][j];
      }
    }
  } // end of externalMeshFlag if statement.

  // these arrays aren't worth copying,
  // but if they were allocated in the original class,
  // they need allocating here too.
  if (adjInfoAllocFlag)
  {
    aiBegin = new int [numRegLabels];
    aiEnd   = new int [numRegLabels];
  }

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::operator=
//
// Purpose       : assignment operator
//
// Special Notes : This is probably more thorough than it absolutely needs
//                 to be.
// Scope         :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/02
//-----------------------------------------------------------------------------
PDE_2DMesh & PDE_2DMesh::operator=(PDE_2DMesh const & rhsMesh)
{
  meshFileName     = rhsMesh.meshFileName;
  externalMeshFlag = rhsMesh.externalMeshFlag;
  xMax             = rhsMesh.xMax;
  yMax             = rhsMesh.yMax;
  xMin             = rhsMesh.xMin;
  yMin             = rhsMesh.yMin;
  dx               = rhsMesh.dx;
  dy               = rhsMesh.dy;
  xRatio           = rhsMesh.xRatio;
  yRatio           = rhsMesh.yRatio;
  x0               = rhsMesh.x0;
  meshScaledFlag   = rhsMesh.meshScaledFlag;
  vol              = rhsMesh.vol;
  invVol           = rhsMesh.invVol;
  surfArea         = rhsMesh.surfArea;
  circum           = rhsMesh.circum;
  invCircum        = rhsMesh.invCircum;
  depth            = rhsMesh.depth;
  numNodes         = rhsMesh.numNodes;
  numEdges         = rhsMesh.numEdges;
  numCells         = rhsMesh.numCells;
  numLabels        = rhsMesh.numLabels;
  numRegLabels     = rhsMesh.numRegLabels;
  numBndryNodes    = rhsMesh.numBndryNodes;
  maxNodeNN        = rhsMesh.maxNodeNN;
  iRecentCellLookup= rhsMesh.iRecentCellLookup;
  cylGeom          = rhsMesh.cylGeom;
  dopingSet        = rhsMesh.dopingSet;
  mNodeVector      = rhsMesh.mNodeVector;
  mEdgeVector      = rhsMesh.mEdgeVector;
  mCellVector      = rhsMesh.mCellVector;
  mLabelVector     = rhsMesh.mLabelVector;
  dopingVector     = rhsMesh.dopingVector;
  xVector          = rhsMesh.xVector;
  yVector          = rhsMesh.yVector;
  visitCellFlagVec = rhsMesh.visitCellFlagVec;
  mLabelMap        = rhsMesh.mLabelMap;
  ixMax            = rhsMesh.ixMax;
  iyMax            = rhsMesh.iyMax;
  devOptions_      = rhsMesh.devOptions_;
  afVisitedVec     = rhsMesh.afVisitedVec;
  adjInfoAllocFlag = rhsMesh.adjInfoAllocFlag;

  int i,j;

  // if we're copying over a mesh that was generated via the "internal"
  // capability, we need to copy over these 3 structures:
  if (!externalMeshFlag)
  {
    if (rhsMesh.nodeIndices != NULL)
    {
      nodeIndices = new int*[ixMax+10];
      for (i=0;i<ixMax+10;++i)
      {
        nodeIndices[i] = new int[iyMax+10];

        for (j=0;j<iyMax+10;++j)
        {
          nodeIndices[i][j] = rhsMesh.nodeIndices[i][j];
        }
      }
    }

    if (rhsMesh.edgeIndices != NULL)
    {
      edgeIndices = new int*[numNodes+10];
      for (i=0;i<numNodes+10;++i)
      {
        edgeIndices[i] = new int[numNodes+10];

        for (j=0;j<numNodes+10;++j)
        {
          edgeIndices[i][j] = rhsMesh.edgeIndices[i][j];
        }
      }
    }

    if (rhsMesh.cellIndices != NULL)
    {
      cellIndices = new int *[ixMax+10];
      for (i=0;i<ixMax+10;++i)
      {
        cellIndices[i] = new int[iyMax+10];
        for (j=0;j<iyMax+10;++j)
          cellIndices[i][j] = rhsMesh.cellIndices[i][j];
      }
    }
  } // end of externalMeshFlag if statement.

  // these arrays aren't worth copying,
  // but if they were allocated in the original class,
  // they need allocating here too.
  if (adjInfoAllocFlag)
  {
    aiBegin = new int [numRegLabels];
    aiEnd   = new int [numRegLabels];
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::~PDE_2DMesh
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
PDE_2DMesh::~PDE_2DMesh ()
{
  if (!externalMeshFlag)
  {
    int i;
    if( nodeIndices!=NULL)
    {
      for (i=0;i<ixMax+10;++i)
      {
        delete [] nodeIndices[i];
      }
      delete [] nodeIndices;
    }

    if( edgeIndices!=NULL)
    {
      for (i=0;i<numNodes+10;++i)
      {
        delete [] edgeIndices[i];
      }
      delete [] edgeIndices;
    }

    if( cellIndices!=NULL)
    {
      for (i=0;i<ixMax+10;++i)
      {
        delete [] cellIndices[i];
      }
      delete [] cellIndices;
    }
  }// end of externalMeshFlag


  if (aiBegin != NULL) delete [] aiBegin;
  if (aiEnd   != NULL) delete [] aiEnd;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::initializeMesh
// Purpose       : This function initializes the mesh, as well as many of
//                 the geometry-related variables.
//
// Special Notes : If a mesh file is specified , this function will attempt
//                 to read it in.  If none is specified, it will set up
//                 a very simple default cartesian mesh.
//
//                 In any case, it is expected that the file, or the
//                 default mesh setup, will contain a list of nodes,
//                 edges and cells.  The rest will be calculated later.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::initializeMesh (const std::string & meshFileName_tmp)
{
  bool bsuccess = true;
  bool tmpBool = true;

  if (true) // if mesh file is specified:  (for now no real test here...)
  {
    externalMeshFlag = true;
    meshFileName = meshFileName_tmp;
    bsuccess = readSGFMeshFile  (meshFileName_tmp);
  }

  // set up extra data structures:

  // set up cell Node lists:
  tmpBool = cellNodes();
  bsuccess = bsuccess && tmpBool;

  // set up the label node information.

  // set up node nearest neighbor lists:

  // calculate geometry information:
  tmpBool = setupGeometry ();
  bsuccess = bsuccess && tmpBool;

  visitCellFlagVec.resize(numCells,0);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
    dumpMesh();

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::readSGFMeshFile
// Purpose       : This function reads in a  *.msh file generated by the
//                 SGF meshing program, and places the information contained
//                 in the *msh file into Xyce data structures.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::readSGFMeshFile (const std::string & meshFileName_tmp)
{
  MESHHEAD meshhead;

  // open the mesh file
  FILE *nFile = fopen(meshFileName_tmp.c_str(), "r");

  if (nFile == NULL)
  {
    Report::UserFatal() << "PDE_2DMesh::readSGFMeshFile - "
                        << meshFileName_tmp 
                        << " file not found.";
  }

  // read the header
  fread(&meshhead, sizeof(MESHHEAD),1,nFile);

  // advance past the constants
  long lOffset = meshhead.cConstant * sizeof(XLATCONST) + sizeof(UINT);
  fseek(nFile, lOffset, SEEK_CUR);

  // read the label lists
  UINT cLabel;
  fread(&cLabel, sizeof(UINT),1,nFile);

  UINT i,j;

  numLabels = cLabel;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "About to read label lists:  cLabel = " << cLabel << std::endl;

  for(i = 0; i < cLabel; ++i)
  {
    XLATLABEL xlatlabel;
    mLabel xLabel;

    fread(&xlatlabel, sizeof(XLATLABEL),1,nFile);
    xLabel.name = std::string(xlatlabel.szName);
    xLabel.iIndex = xlatlabel.iIndex;
    xLabel.uType  = xlatlabel.uType;
    xLabel.cNode  = xlatlabel.cNode;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "i = " << i ;
      Xyce::dout() << "  iIndex = " << xLabel.iIndex;
      Xyce::dout() << "  uType = " << xLabel.uType;
      Xyce::dout() << "  cNode = " << xLabel.cNode;
      Xyce::dout() << "  label name = " << xLabel.name << std::endl;
    }

    xLabel.mNodeVector.reserve(xlatlabel.cNode+1);

    for (j=0;j<xlatlabel.cNode; ++j)
    {
      UINT au;
      fread(&au, sizeof(UINT),1,nFile);
      xLabel.mNodeVector.push_back(au);
    }

    xLabel.mNodeVector[xlatlabel.cNode] = -1;

    ExtendedString tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;
    mLabelVector.push_back(xLabel);
  }

  // Initialize the mesh arrays
  // (the first two mesh arrays are the x and y values, indexed by node)
  // The other stuff, if it exists, is probably stuff that was
  // defined to help refine the mesh, such as the doping concentration.
  //
  // Since the information contained in the x and y arrays are also
  // contained inside the node array, I'm going to skip that part.
  // The other arrays, however, might be useful, so I'm going to read
  // those in.
  //
  UINT cArray;
  fread(&cArray, sizeof(UINT),1,nFile);
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "About to read in the arrays.  cArray = " << cArray << std::endl;

  // the first two arrays are x and y...
  double xtmp;
  for(i = 0; i < 2; ++i)
  {
    char szArrayName[LEN_IDENT+1];
    fread(szArrayName, LEN_IDENT+1,1,nFile);
    ExtendedString tmpName(szArrayName);
    tmpName.toUpper ();
    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "i=" <<i<< "  name = " << tmpName << std::endl;
    }

    if (tmpName == "X")
    {
      xVector.reserve(meshhead.cNode);
      for (j=0;j<meshhead.cNode; ++j)
      {
        fread(&xtmp, sizeof(double),1,nFile);
        xVector.push_back(xtmp);
      }
    }
    else if (tmpName == "Y")
    {
      yVector.reserve(meshhead.cNode);
      for (j=0;j<meshhead.cNode; ++j)
      {
        fread(&xtmp, sizeof(double),1,nFile);
        yVector.push_back(xtmp);
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done reading x and y arrays.";
    Xyce::dout() << "  Reading the others, if they exist." << std::endl;
  }

  // Now read and store the other arrays, if any:
  if (meshhead.cArray > 2)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      Xyce::dout() << "reading extra arrays..." << std::endl;
    for(i = 2; i < meshhead.cArray; ++i)
    {
      char szArrayName[LEN_IDENT+1];
      fread(szArrayName, LEN_IDENT+1,1,nFile);
      ExtendedString tmpName(szArrayName);
      tmpName.toUpper ();
      if (tmpName == "C")
      {
        dopingVector.reserve(meshhead.cNode);
        for (j=0;j<meshhead.cNode; ++j)
        {
          fread(&xtmp, sizeof(double),1,nFile);
          dopingVector.push_back(xtmp);
        }
        dopingSet = true;
      }
      else
      {
        for (j=0;j<meshhead.cNode; ++j)
          fread(&xtmp, sizeof(double),1,nFile);
      }
    }
  }

  // Read the node array
  UINT cNodes;
  fread(&cNodes, sizeof(UINT),1,nFile);
  mNode xNode;
  numNodes = cNodes;

  numBndryNodes = meshhead.cBndryNode;
  numRegLabels  = meshhead.cRegLabel;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "About to read the nodes.  cNodes = " << cNodes << std::endl;

  for (i=0;i<cNodes;++i)
  {
    NODE sgfNode;
    fread(&sgfNode, sizeof(NODE),1,nFile);

    // copy this node information into the Xyce data structure, mNode:
    xNode.x = sgfNode.x;
    xNode.y = sgfNode.y;
    mNodeVector.push_back(xNode);
  }

  // read the edge array
  UINT cEdge;
  fread(&cEdge, sizeof(UINT),1,nFile);
  mEdge xEdge;
  numEdges = cEdge;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "About to read the edges.  cEdge = " << cEdge << std::endl;

  for (i=0; i< cEdge; ++i)
  {
    EDGE sgfEdge;
    fread(&sgfEdge, sizeof(EDGE),1,nFile);

    int inode1 = sgfEdge.inodeA;
    int inode2 = sgfEdge.inodeB;

    if (inode1 > inode2)
    {
      sgfEdge.inodeA = inode2;
      sgfEdge.inodeB = inode1;
    }

    // copy this edge information into the Xyce data structure, mEdge:
    xEdge.uLabel = sgfEdge.uLabel;
    xEdge.inodeA = sgfEdge.inodeA;
    xEdge.inodeB = sgfEdge.inodeB;
    xEdge.iedge  = i;

    mEdgeVector.push_back(xEdge);
  }

  // Read the triangle array
  UINT cTriangle;
  fread(&cTriangle, sizeof(UINT),1,nFile);
  mCell xCell;
  numCells = cTriangle;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "About to read the cells.  numCells = " << numCells  << std::endl;
  for (i=0;i<cTriangle;++i)
  {
    TRI Tri;
    fread(&Tri, sizeof(TRI),1,nFile);

    // copy this tri information into the Xyce data structure, mCell:
    xCell.uLabel = Tri.uLabel;

    xCell.iedgeAB = Tri.iedgeAB;
    xCell.iedgeBC = Tri.iedgeBC;
    xCell.iedgeCD = Tri.iedgeAC;
    xCell.iedgeDA = Tri.iedgeAD;

    xCell.icellAB = Tri.itriAB;
    xCell.icellBC = Tri.itriBC;
    xCell.icellCD = Tri.itriAC;
    xCell.icellDA = Tri.itriAD;

    mCellVector.push_back(xCell);
  }

  // read the adjacency information
  UINT cAdj;
  fread(&cAdj, sizeof(UINT),1,nFile);
  numAdj = cAdj;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "About to read the adjacency info.  cAdj = " << cAdj << std::endl;
  // Set up scaling factors for geometry.
  double rDistScale = 1.0;
  double X0_1       = rDistScale;
  double X0_2       = rDistScale * X0_1;
  double X0_3       = rDistScale * X0_2;

  double ScaleILEN, ScaleAREA, ScaleELEN = X0_1;

  cylGeom = meshhead.fCylGeom;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    if (meshhead.fCylGeom) Xyce::dout() << "Cylindrical Geometry\n\n";
    else                   Xyce::dout() << "Cartesian Geometry\n\n";
  }

  if (meshhead.fCylGeom)
  {
    ScaleILEN = X0_2;
    ScaleAREA = X0_3;
  }
  else
  {
    ScaleILEN = X0_1;
    ScaleAREA = X0_2;
  }

  // Note:  cAdj (numAdj) may be larger than numNodes.
  //  The reason for this is that nodes along boundaries between
  //  regions have more than one adjacency info structure.
  // ( Fix this later....)
  //
  // The formula used by simgen to get cAdj is:
  //   cAdj = cNode + cBndryNode * cRegLabel;
  //
  //   cNode      = number of nodes
  //   cBndryNode = number of nodes along boundaries between regions.
  //   cRegLabel  = number of region labels, such as SI and SIO2.
  //
  for(i = 0; i < cAdj; ++i)
  {
    int inode;
    fread(&inode, sizeof(int),1,nFile);

    NODEINFO nodeInfo;
    fread( &nodeInfo, sizeof(NODEINFO) - sizeof(EDGEINFO *),1,nFile);
    nodeInfo.Area /= ScaleAREA;
    UINT cNeighbor = nodeInfo.cNeighbor;

    mNodeVector[inode].area  = nodeInfo.Area;
    mNodeVector[inode].cnode = cNeighbor;
    mNodeVector[inode].inode = inode;
    mNodeVector[inode].numCells = nodeInfo.cTriangle;

    EDGEINFO EdgeInfo;
    for (j=0;j<cNeighbor;++j)
    {
      fread( &EdgeInfo, sizeof(EDGEINFO),1,nFile);
      mNodeVector[inode].edgeInfoVector.push_back(EdgeInfo);
    }

    for(j = 0; j < cNeighbor; ++j)
    {
      // Since I set rDistScale to 1.0, these scalings are kind of
      // irrelevant.  The distances scaling can technically be set from
      // the mesh file, which is why I left these lines in here.
      mNodeVector[inode].edgeInfoVector[j].ilen  /= ScaleILEN;
      mNodeVector[inode].edgeInfoVector[j].elen  /= ScaleELEN;
      mNodeVector[inode].edgeInfoVector[j].Area1 /= ScaleAREA;
      mNodeVector[inode].edgeInfoVector[j].Area2 /= ScaleAREA;

      int iedge = mNodeVector[inode].edgeInfoVector[j].iedge;
      mEdgeVector[iedge].elen  = mNodeVector[inode].edgeInfoVector[j].elen;
      mEdgeVector[iedge].ilen  = mNodeVector[inode].edgeInfoVector[j].ilen;
      mEdgeVector[iedge].Area1 = mNodeVector[inode].edgeInfoVector[j].Area1;
      mEdgeVector[iedge].Area2 = mNodeVector[inode].edgeInfoVector[j].Area2;
    }
  }

  // close the file
  fclose(nFile);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::writeSGFMeshFile
// Purpose       : This function writes a  *.msh file , which is the format of
//                 the SGF meshing program.  You would call this function if
//                 the mesh used was generated internally, but you wanted to
//                 use sgplot to view your results.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/04/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::writeSGFMeshFile (const std::string & meshFileName_tmp)
{
  MESHHEAD meshhead;
  int i,j;

  // open the mesh file
  FILE *nFile = fopen(meshFileName_tmp.c_str(), "w");

  strcpy(meshhead.szLogo, "SGFramework Mesh File  Version 1.0\n\x1a");
  strcpy(meshhead.szSign, "@~!__MESH__!~@");

  meshhead.cConstant  = 0; // no constants.
  meshhead.cLabel     = mLabelVector.size();
  meshhead.cArray     = 2; // arrays:  x and y  (maybe C also?)
  meshhead.cRegLabel  = numRegLabels; // usually just 1 of these...
  meshhead.cBndryNode = numBndryNodes;
  // fix this later!?  The sgplot program doesn't care about the boundary
  // node issue it does the same plotting whether it is a boundary node
  // or not.  Other SGF related programs do care, however, and for those
  // programs, the boundary nodes need to be listed first.

  meshhead.cNode      = numNodes;
  meshhead.cEdge      = numEdges;
  meshhead.cTriangle  = numCells;
  meshhead.fCylGeom   = cylGeom;
  fwrite( &meshhead, sizeof(MESHHEAD), 1, nFile);

  // write the constants (should be zero of these...)
  // sgplot will advance past them anyway, if they were there.
  fwrite( &meshhead.cConstant, sizeof(UINT), 1, nFile);

  // write the labels:
  fwrite( &meshhead.cLabel, sizeof(UINT), 1, nFile);
  for (i=0;i<meshhead.cLabel;++i)
  {
    // output label itself
    XLATLABEL xlabel;
    strcpy(xlabel.szName, mLabelVector[i].name.c_str());
    xlabel.iIndex = mLabelVector[i].iIndex;
    xlabel.uType  = mLabelVector[i].uType;
    xlabel.cNode  = mLabelVector[i].cNode;
    fwrite(&xlabel,sizeof(XLATLABEL) ,1,nFile);

    // output node list for this label
    // these are skipped by sgplot.
    for (j=0;j<xlabel.cNode;++j)
    {
      UINT istuf = mLabelVector[i].mNodeVector[j];
      fwrite(&istuf, sizeof(UINT), 1, nFile);
    }
  }

  // write the 1D arrays:
  // write the x and y coordinate arrays
  // cArray, for some reason, is actually the true number of arrays -2.
  // sgplot, etc., assumes there will always be x and y arrays, so the
  // the variable cArray is only for "extra" arrays.  Seems needlessly
  // confusing to me...
  UINT cArray = 0;
  fwrite( &cArray, sizeof(UINT),1,nFile);
  char szCoord1[LEN_IDENT+1];
  strcpy(szCoord1, "X");
  fwrite( szCoord1, (LEN_IDENT+1), 1, nFile);

  for(i = 0; i < numNodes; ++i)
  {
    fwrite( &(xVector[i]), sizeof(double), 1, nFile);
  }

  char szCoord2[LEN_IDENT+1];
  strcpy(szCoord2, "Y");
  fwrite( szCoord2, (LEN_IDENT+1), 1, nFile);
  for(i = 0; i < numNodes; ++i)
  {
    fwrite( &(yVector[i]), sizeof(double), 1, nFile);
  }

  // write the nodes:
  // As the nodes are mostly redundant with the x and y arrays, which
  // were just output, sgplot will skip over the nodes.
  fwrite( &meshhead.cNode, sizeof(UINT), 1, nFile);
  for (i=0;i<meshhead.cNode;++i)
  {
    NODE n1;
    n1.x = xVector[i];  n1.y = yVector[i];
    fwrite(&n1, sizeof(NODE), 1, nFile);
  }

  // write the edges:
  fwrite( &meshhead.cEdge, sizeof(UINT), 1, nFile);
  for (i=0;i<meshhead.cEdge;++i)
  {
    EDGE  e1;
    e1.uLabel = mEdgeVector[i].uLabel;
    e1.inodeA = mEdgeVector[i].inodeA;
    e1.inodeB = mEdgeVector[i].inodeB;
    fwrite( &e1, sizeof(EDGE), 1, nFile);
  }

  // write the triangles (cells):
  fwrite(&meshhead.cTriangle, sizeof(UINT),1,nFile);

  // valgrind reports an error in this loop for some reason. I can't
  // figure out why.
  for(i = 0; i < meshhead.cTriangle; ++i)
  {
    TRI t1;

    int tmp = mCellVector[i].uLabel;

    if (tmp < 0) t1.uLabel = -1u;
    else         t1.uLabel = static_cast<unsigned int>(tmp);

    t1.iedgeAB = mCellVector[i].iedgeAB;
    t1.iedgeBC = mCellVector[i].iedgeBC;
    t1.iedgeAC = mCellVector[i].iedgeCD;
    t1.iedgeAD = mCellVector[i].iedgeDA;

    t1.itriAB = mCellVector[i].icellAB;
    t1.itriBC = mCellVector[i].icellBC;
    t1.itriAC = mCellVector[i].icellCD;
    t1.itriAD = mCellVector[i].icellDA;

    fwrite(&t1, sizeof(TRI), 1, nFile);
  }

  // write the node adjacency table:
  // This part isn't read in by sgplot, so it probably isn't neccessary.
  // But, what the heck, maybe it will be useful at some point.
  fwrite(&numNodes, sizeof(UINT), 1, nFile);
  for (i=0;i<numNodes;++i)
  {
    UINT inode = i;
    fwrite(&inode, sizeof(UINT),1,nFile);
    NODEINFO nodeinfo;
    // fill in nodeinfo w/information from mNode class.
    nodeinfo.Area = mNodeVector[i].area;
    nodeinfo.cNeighbor = mNodeVector[i].cnode;
    nodeinfo.cTriangle = mNodeVector[i].numCells;

    fwrite(&nodeinfo,sizeof(NODEINFO)-sizeof(EDGEINFO *),1,nFile);

    for (j=0;j<mNodeVector[i].cnode;++j)
    {
      EDGEINFO edgeinfo = mNodeVector[i].edgeInfoVector[j];
      fwrite(&edgeinfo,sizeof(EDGEINFO),1,nFile);
    }
  }

  // close the file
  fclose(nFile);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::initializeInternalMesh
//
// Purpose       : This function initializes a simple cartesian mesh,
//                 as well as many of the geometry-related variables.
//
// Special Notes : This function is called if the user has specified
//                 that the mesh file name is "internal".  In other words,
//                 there is no external mesh file.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::initializeInternalMesh
    (int nx,
     int ny,
     double xlength,
     double ylength,
     int numElectrodes,
     std::string & outputMeshFileName,
     std::map<std::string,PDE_2DElectrode*> & elMap,
     bool cylFlag)
{
  bool bsuccess = true;
  bool tmpBool = true;

  cylGeom = cylFlag;

  externalMeshFlag = false;

  meshFileName = outputMeshFileName; // this variable is used for output.

  bsuccess = setupInternalMesh (nx,ny,xlength,ylength);

  // If the user has not specified electrodes, then assume a default set,
  // and call setupDefaultLabels.  Need to do a little error checking
  // first.
  if ( !(elMap.empty()) )
  {
    tmpBool = errorCheckElectrodes (numElectrodes, elMap);
    bsuccess = bsuccess && tmpBool;
  }

  if ( useDefaultLabels )
  {
    tmpBool = setupDefaultLabels(numElectrodes);
    bsuccess = bsuccess && tmpBool;
  }
  else
  {
    tmpBool = setupInternalLabels(numElectrodes, elMap);
    bsuccess = bsuccess && tmpBool;
  }

  // set up extra data structures:

  // set up cell Node lists:
  tmpBool = cellNodes();
  bsuccess = bsuccess && tmpBool;

  // set up node nearest neighbor lists:

  // calculate geometry information:
  tmpBool = setupGeometry ();
  bsuccess = bsuccess && tmpBool;

  visitCellFlagVec.resize(numCells,0);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_DUMP_VECTORS))
  {
    dumpMesh ();
  }

  if (sgplotLevel > 0)
  {
    tmpBool = writeSGFMeshFile (outputMeshFileName);
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::resizeMesh
// Purpose       : This function resizes the mesh to a new length and
//                 width.
//
// Special Notes : This function has a similar function to that of
//                 scaleMesh.  However, unlike scaleMesh, the x-scaling and
//                 y-scaling can be different.  Unfortunately, that makes
//                 the work of this function more complicated. Instead of
//                 simply scaling areas, it is neccessary to recompute
//                 them.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::resizeMesh(double xlength, double ylength)
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In PDE_2DMesh:resizeMesh." << std::endl;
  }

  double old_xlength = xMax-xMin;
  double old_ylength = yMax-yMin;

  xRatio = xlength/(old_xlength);
  yRatio = ylength/(old_ylength);

  // xMin and yMin will stay in the same place, while xMax and yMax move.
  xMax = xMin + xlength;
  yMax = yMin + ylength;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << " xMax = " << xMax << std::endl;
    Xyce::dout() << " yMax = " << yMax << std::endl;
  }

  // Adjust the x and y arrays:
  int i;
  for (i=0;i<numNodes;++i)
  {
    xVector[i] -= xMin;
    xVector[i] *= xRatio;
    xVector[i] += xMin;

    yVector[i] -= yMin;
    yVector[i] *= yRatio;
    yVector[i] += yMin;

    mNodeVector[i].x = xVector[i];
    mNodeVector[i].y = yVector[i];

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      Xyce::dout() << "\txVec["<<i<<"] = " << xVector[i];
      Xyce::dout() << "\tyVec["<<i<<"] = " << yVector[i] << std::endl;
    }
  }

  // resize everything else:
  calcAdjacencyInfo ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with PDE_2DMesh:resizeMesh." << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::setupInternalMesh
// Purpose       : This function sets up a simple default mesh.  It is only
//                 called in the event that the user has not specified a file.
//
// Special Notes : The default mesh is a nx x ny cartesian mesh.
//
//                 First, the code sets up a simple grid of nodes, edges
//                 and cells.
//
//                 After getting the grid in place, it then calculates the
//                 geometrical information, such as edge lengths,
//                 integration box areas, etc.  For a cartesian grid like
//                 this, this sort of information is pretty obvious.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::setupInternalMesh (
    int nx, int ny, double xlength, double ylength)
{
  int i,j;
  ixMax = nx;
  iyMax = ny;

  numNodes     = ixMax*iyMax;

  numRegLabels = 1;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    Xyce::dout() << "In setupInternalMesh"<<std::endl;

  xMax = xlength; xMin = 0.0;
  yMax = ylength; yMin = 0.0;

  dx   = xMax/(static_cast<double>(ixMax-1));
  dy   = yMax/(static_cast<double>(iyMax-1));

  // the "+10" is just for safety purposes.
  nodeIndices = new int*[ixMax+10];
  for (i=0;i<ixMax+10;++i)
  {
    nodeIndices[i] = new int[iyMax+10];

    for (j=0;j<iyMax+10;++j)
    {
      nodeIndices[i][j] = -1;
    }
  }

  // set up the node vector
  int nodeIndex = 0;
  numBndryNodes = 0;  // if only 1 region, this should stay zero.

  for (i=0;i<ixMax;++i)
  {
    for (j=0;j<iyMax;++j)
    {
      mNode n1;
      n1.x = dx*(static_cast<double>(i));
      n1.y = dy*(static_cast<double>(j));

      nodeIndices[i][j] = nodeIndex;

      ++nodeIndex;
      if (i==0 || j==0 || i==ixMax-1 || j==iyMax-1)
      {
        n1.edgeStatus = EDGESTATUS_EXTERIOR;
      }
      else
      {
        n1.edgeStatus = EDGESTATUS_INTERIOR;
      }

      mNodeVector.push_back(n1);
    }
  }
  numNodes = nodeIndex;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "\nDone with Node vector" << std::endl;
  }

  edgeIndices = new int*[numNodes+10];
  for (i=0;i<numNodes+10;++i)
  {
    edgeIndices[i] = new int[numNodes+10];

    for (j=0;j<numNodes+10;++j)
    {
      edgeIndices[i][j] = -1;
    }
  }

  // set up the edge vector
  int edgeIndex = 0;

  for (i=0;i<ixMax;++i)
  {
    for (j=0;j<iyMax;++j)
    {
      // edges parallel to the y-axis.
      if (j>0)
      {
        mEdge edge1;
        edge1.uLabel = -1;

        // note:  inodeA should always be smaller than inodeB!
        edge1.inodeA = nodeIndices[i][j  ];
        edge1.inodeB = nodeIndices[i][j-1];

        if (edge1.inodeA > edge1.inodeB)
        {
          edge1.inodeB = nodeIndices[i][j  ];
          edge1.inodeA = nodeIndices[i][j-1];
        }

        if (mNodeVector[edge1.inodeA].edgeStatus ==
            mNodeVector[edge1.inodeB].edgeStatus)
        {
          edge1.edgeStatus = mNodeVector[edge1.inodeA].edgeStatus;
        }

        if (edge1.inodeA == -1 || edge1.inodeB == -1)
        {
          if (DEBUG_DEVICE)
          {
            Xyce::dout() << " edge1.inodeA = " << edge1.inodeA;
            Xyce::dout() << "  edge1.inodeB = " << edge1.inodeB <<std::endl;
          }
          Report::DevelFatal() << "Failed on edge1";
        }

        edgeIndices[edge1.inodeA][edge1.inodeB] = edgeIndex;
        edgeIndices[edge1.inodeB][edge1.inodeA] = edgeIndex;

        mEdgeVector.push_back(edge1);
        ++edgeIndex;
      }

      // edges parallel to the x-axis.
      if (i>0)
      {
        mEdge edge2;
        edge2.uLabel = -1;

        // note:  inodeA should always be smaller than inodeB!
        edge2.inodeA = nodeIndices[i  ][j];
        edge2.inodeB = nodeIndices[i-1][j];

        if (edge2.inodeA > edge2.inodeB)
        {
          edge2.inodeB = nodeIndices[i  ][j];
          edge2.inodeA = nodeIndices[i-1][j];
        }

	      if (mNodeVector[edge2.inodeA].edgeStatus ==
	          mNodeVector[edge2.inodeB].edgeStatus)
	      {
	        edge2.edgeStatus = mNodeVector[edge2.inodeA].edgeStatus;
        }

        if (edge2.inodeA == -1 || edge2.inodeB == -1)
        {
          if (DEBUG_DEVICE)
          {
            Xyce::dout() << " edge2.inodeA = " << edge2.inodeA;
            Xyce::dout() << "  edge2.inodeB = " << edge2.inodeB <<std::endl;
          }
          Report::DevelFatal() << "Failed on edge2";
        }

        edgeIndices[edge2.inodeA][edge2.inodeB] = edgeIndex;
        edgeIndices[edge2.inodeB][edge2.inodeA] = edgeIndex;

        mEdgeVector.push_back(edge2);
        ++edgeIndex;
      }
    }
  }

  numEdges = edgeIndex;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with Edge vector" << std::endl;
  }

  // set up the cell vector:
  int cellIndex = 0;
  cellIndices = new int *[ixMax+10];
  for (i=0;i<ixMax+10;++i)
  {
    cellIndices[i] = new int[iyMax+10];
    for (j=0;j<iyMax+10;++j)
      cellIndices[i][j] = -1;
  }

  for (i=1; i<ixMax; ++i)
  {
    for (j=1;j<iyMax;++j)
    {
        // All these cells are cartesean rectangles:
        //
        //        A                     B
        //    (i-1, j-1) ---------  (i  , j-1)
        //         |                    |
        //         |                    |
        //         |                    |
        //         |                    |
        //         |                    |
        //    (i-1, j  ) ---------  (i  , j  )
        //        D                     C
        //
        int inodeA, inodeB, inodeC, inodeD;

        inodeA = nodeIndices[i-1][j-1];
        inodeB = nodeIndices[i  ][j-1];
        inodeC = nodeIndices[i  ][j  ];
        inodeD = nodeIndices[i-1][j  ];

        mCell c1;
        c1.iedgeAB = edgeIndices[inodeA][inodeB];
        c1.iedgeBC = edgeIndices[inodeB][inodeC];
        c1.iedgeCD = edgeIndices[inodeC][inodeD];
        c1.iedgeDA = edgeIndices[inodeD][inodeA];

        c1.inodeA = inodeA;
        c1.inodeB = inodeB;
        c1.inodeC = inodeC;
        c1.inodeD = inodeD;

        // Assumes always has SI label.  See setupInternalLabels.
        c1.uLabel = 0;

        mCellVector.push_back(c1);
        cellIndices[i][j] = cellIndex;

        ++cellIndex;
    }
  }
  numCells = cellIndex;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with Cell vector" << std::endl;
  }

  // set up the neighbor cells.  If we go off the mesh, the
  // value in the indices array is -1 anyway, so there will
  // be no problem.
  for (i=1; i<ixMax; ++i)
  {
    for (j=1;j<iyMax;++j)
    {
      int cellIndexTmp = cellIndices[i][j];
      mCellVector[cellIndexTmp].icellAB = cellIndices[i  ][j-1];
      mCellVector[cellIndexTmp].icellBC = cellIndices[i+1][j  ];
      mCellVector[cellIndexTmp].icellCD = cellIndices[i  ][j+1];
      mCellVector[cellIndexTmp].icellDA = cellIndices[i-1][j  ];
    }
  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with Cell neighbors." << std::endl;
  }

  // if neccessary re-order nodes, edges, etc.  SGF always orders
  // nodes in the following manner:
  //    1) all region boundary nodes
  //    2) all exterior boundary nodes, which are not region boundary nodes.
  //    3) all interior nodes.
  // NOT DONE YET...

  // copy over the x and y info into the xVector and yVector data
  // structures.
  xVector.reserve(numNodes);
  yVector.reserve(numNodes);
  for (i=0;i<numNodes;++i)
  {
    xVector.push_back(mNodeVector[i].x);
    yVector.push_back(mNodeVector[i].y);
  }

  // now that the grid relationships are established, now get the
  // edge lengths, areas, etc.
  calcAdjacencyInfo ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with setupInternalMesh\n";
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::setupInternalAdjacencyInfo
//
// Purpose       : This function is being phased out.  Eventually,
//                 calcAdjacencyInfo will be the one to call, once it has
//                 been adequately refactored.
//
//                 This function calculates adjacency information, but can
//                 only do it right if the mesh is cartesian.
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::setupInternalAdjacencyInfo ()
{
  int i,j;

  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "In PDE_2DMesh::setupInternalAdjacencyInfo." << std::endl;
  }

  // first do some of the edge information.
  for (i=0;i<numEdges;++i)
  {
    mEdge & edgeTmp = mEdgeVector[i];
    int inodeA = edgeTmp.inodeA;
    int inodeB = edgeTmp.inodeB;

    double DELx = fabs(xVector[inodeA]-xVector[inodeB]);
    double DELy = fabs(yVector[inodeA]-yVector[inodeB]);

    // assuming cartesian, as before.
    if (DELy > DELx)  // y-axis parallel
    {
      edgeTmp.elen = dy;
      edgeTmp.ilen = dx;
    }
    else  // x-axis parallel
    {
      edgeTmp.elen = dx;
      edgeTmp.ilen = dy;
    }

    if (edgeTmp.edgeStatus==EDGESTATUS_EXTERIOR)
    {
      edgeTmp.ilen *= 0.5;
    }

    edgeTmp.Area1 = (edgeTmp.elen * edgeTmp.ilen) * 0.25;
    edgeTmp.Area2 = 0.0; // need to fix Area2 later!
  }

  // Now do the "edgeinfo" information, which is owned by each node.
  // Get nearest neighbor node information, edgeinfo array stuff, etc.:
  // Set the area of the integration box as well.
  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "About to do the edgeinfo:" << std::endl;
  }

  for (i=0;i<ixMax;++i)
  {
    for (j=0;j<iyMax;++j)
    {
      int node = nodeIndices[i][j];
      mNodeVector[node].area = dx*dy;

      mNodeVector[node].fBndry = false;

      if (i==0 || i==ixMax-1)
      {
        mNodeVector[node].area *= 0.5;
      }

      if (j==0 || j==iyMax-1)
      {
        mNodeVector[node].area *= 0.5;
      }

      int nnTmp;
      int edgeTmp;
      int cnodeTmp = 0;
      EDGEINFO eiTmp;

      if (i<ixMax-1)
      {
        nnTmp = nodeIndices[i+1][j];
        edgeTmp = edgeIndices[node][nnTmp];
        eiTmp.inode = nnTmp;
        eiTmp.iedge = edgeTmp;
        eiTmp.elen  = mEdgeVector[edgeTmp].elen;
        eiTmp.ilen  = mEdgeVector[edgeTmp].ilen;
        mNodeVector[node].edgeInfoVector.push_back(eiTmp);
        ++cnodeTmp;
      }

      if (j<iyMax-1)
      {
        nnTmp = nodeIndices[i][j+1];
        edgeTmp = edgeIndices[node][nnTmp];
        eiTmp.inode = nnTmp;
        eiTmp.iedge = edgeTmp;
        eiTmp.elen  = mEdgeVector[edgeTmp].elen;
        eiTmp.ilen  = mEdgeVector[edgeTmp].ilen;
        mNodeVector[node].edgeInfoVector.push_back(eiTmp);
        ++cnodeTmp;
      }

      if (i>0)
      {
        nnTmp = nodeIndices[i-1][j];
        edgeTmp = edgeIndices[node][nnTmp];
        eiTmp.inode = nnTmp;
        eiTmp.iedge = edgeTmp;
        eiTmp.elen  = mEdgeVector[edgeTmp].elen;
        eiTmp.ilen  = mEdgeVector[edgeTmp].ilen;
        mNodeVector[node].edgeInfoVector.push_back(eiTmp);
        ++cnodeTmp;
      }

      if (j>0)
      {
        nnTmp = nodeIndices[i][j-1];
        edgeTmp = edgeIndices[node][nnTmp];
        eiTmp.inode = nnTmp;
        eiTmp.iedge = edgeTmp;
        eiTmp.elen  = mEdgeVector[edgeTmp].elen;
        eiTmp.ilen  = mEdgeVector[edgeTmp].ilen;
        mNodeVector[node].edgeInfoVector.push_back(eiTmp);
        ++cnodeTmp;
      }
      mNodeVector[node].cnode = cnodeTmp;

    } // j
  } // i
  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "Done doing the edgeinfo:" << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::errorCheckElectrodes
//
// Purpose       : This function checks that the user specification of
//                 electrodes is consistent.
//
// Special Notes :
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/18/03
//-----------------------------------------------------------------------------
bool PDE_2DMesh::errorCheckElectrodes
  (int numElectrodes, std::map<std::string,PDE_2DElectrode*> & elMap)
{
  bool bsuccess = true;

  // check that numElectrodes = elMap size.
  if (numElectrodes != elMap.size())
  {
    bsuccess = false;
    Report::UserFatal() <<  "Number of electrodes and number of nodes are not equal.";
  }


  // Check if the "start", "end" and "side" parameters
  // have been set or not.  Either they all should be set, for every
  // electrode, or none of them.  If all set, then use the
  // setupInternalLabels function (later), if none set, then use the
  // setupDefaultLabels function.

  bool allSet = true;
  bool allNotSet = true;
  useDefaultLabels = false;

  std::map<std::string, PDE_2DElectrode*>::iterator mapIter;
  std::map<std::string, PDE_2DElectrode*>::iterator mapStart = elMap.begin();
  std::map<std::string, PDE_2DElectrode*>::iterator mapEnd = elMap.end();

  for ( mapIter = mapStart; mapIter != mapEnd; ++mapIter )
  {
    PDE_2DElectrode & el = *(mapIter->second);

    allSet = allSet && (el.given("START") && el.given("END") && el.given("SIDE"));
    allNotSet = allNotSet &&
		 (!(el.given("START")) && !(el.given("END")) && !(el.given("SIDE")));
  }

  if (!allSet && !allNotSet)
  {
    bsuccess = false;
    Report::UserFatal() << "Some electrodes have start, end and side specified, some don't.  " << std::endl
                        << "Either specify start, end and side for all electrodes, or none.";
  }

  // if "allNotSet", then set the flag which will force the code to call
  // setupDefaultLabels.
  if (allNotSet)
  {
    useDefaultLabels = true;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::setupDefaultLabels
//
// Purpose       : This function sets up the label vector, and assumes a
//                 hardwired set of labels that are uniquely determined by
//                 the number of electrodes.
//
//                 If you have 2 electrodes, this is a SI diode with an
//                 ANODE and CATHODE and NOFLUX edge boundaries.
//
//                 If you have 3 electrodes, this is a silicon (SI) BJT,
//                 with an EMITTER, COLLECTOR, BASE, and NOFLUX.  (not
//                 set up yet).
//
//                 If you specify 4 electrodes, then this is a MOSFET, and
//                 the labels are:   SI, SOURCE, GATE, DRAIN, SUB, NOFLUX
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::setupDefaultLabels(int numberElectrodes)
{
  // set up the silicon region label.  This will always be here.
  int i,j;
  mLabel xLabel;
  ExtendedString tmpName = "SI";
  int nodeTmp = 0;

  xLabel.name = "SI";
  xLabel.iIndex = 0;
  xLabel.uType = TYPE_REGION;

  // This assumes every node in the mesh is SI, except for the edges.
  xLabel.mNodeVector.reserve(ixMax*iyMax); // just to be safe, reserve extra

  // fill in all the nodes:
  nodeTmp = 0;
  for (i=1;i<ixMax-1;++i)
  {
    for (j=1;j<iyMax-1;++j)
    {
      xLabel.mNodeVector.push_back(nodeIndices[i][j]);
      ++nodeTmp;
    }
  }
  xLabel.cNode = xLabel.mNodeVector.size ();

  tmpName = xLabel.name;
  tmpName.toUpper();
  mLabelMap[tmpName] = xLabel;

  mLabelVector.push_back(xLabel);



  // Now do the edge labels:
  if (numberElectrodes==2)
  {
    numLabels = 4; // SI, ANODE, CATHODE, NOFLUX.

    // do the Anode:
    xLabel.mNodeVector.clear();
    xLabel.name = "ANODE";
    xLabel.iIndex = 1;
    xLabel.uType = TYPE_EDGE;
    xLabel.cNode = ixMax;
    xLabel.mNodeVector.reserve(xLabel.cNode+1);

    // fill in all the nodes:
    nodeTmp = 0;
    j=0;
    for (i=0;i<ixMax;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;
    mLabelVector.push_back(xLabel);

    // do the Cathode:
    xLabel.mNodeVector.clear();
    xLabel.name = "CATHODE";
    xLabel.iIndex = 2;
    xLabel.uType = TYPE_EDGE;
    xLabel.cNode = ixMax;
    xLabel.mNodeVector.reserve(xLabel.cNode+1);

    // fill in all the nodes:
    nodeTmp = 0;
    j=iyMax-1;
    for (i=0;i<ixMax;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;
    mLabelVector.push_back(xLabel);

    // set up NOFLUX:
    xLabel.mNodeVector.clear();
    xLabel.name = "NOFLUX";
    xLabel.iIndex = 3;
    xLabel.uType = TYPE_EDGE;
    xLabel.cNode = 2*iyMax-2;
    xLabel.mNodeVector.reserve(xLabel.cNode+1);

    // fill in all the nodes:
    nodeTmp = 0;
    i=0;
    for (j=1;j<iyMax-1;++j)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }
    i=ixMax-1;
    for (j=1;j<iyMax-1;++j)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;
    mLabelVector.push_back(xLabel);

    // Also assign labels to the edges in mEdgeVector.
    // Note:  labels were already assigned to the cells in mCellVector, up
    //        in setupInternalMesh.
    //        Also, nodes don't get labels.

    j=0;
    for (i=1;i<ixMax;++i)
    {
      int n1 = nodeIndices[i  ][j];
      int n2 = nodeIndices[i-1][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 1; // ANODE index
    }

    j=iyMax-1;
    for (i=1;i<ixMax;++i)
    {
      int n1 = nodeIndices[i  ][j];
      int n2 = nodeIndices[i-1][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 2; // CATHODE index
    }

    i=0;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i][j];
      int n2 = nodeIndices[i][j-1];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 3; // NOFLUX index
    }

    i=ixMax-1;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i][j];
      int n2 = nodeIndices[i][j-1];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 3; // NOFLUX index
    }
  } // end of numberElectrodes==2
  else if (numberElectrodes==3)
  {
     // Assuming a BJT mesh:
     //   Ratios:  AB = 1/10
     //            BC = 2/10
     //            CD = 1/10
     //            DE = 3/10
     //            EF = 1/10
     //            FG = 2/10
     //
     //   So, hopefully, the number of mesh cells is
     //   easily divisible by 10.  The default mesh spacing will work.
     //
     //   emitter noflux  base    noflux    collector noflux
     //   A-----B--------C-----D-------------E-----F---------G
     //   |     |  /     |     |    |        |     |         |
     // n |-----|-       |     |   /         |     |         |
     // o |     |        |     |  /          |     |         |
     // f |-----|--------|-----|-            |     |         | noflux
     // l |     |        |     |             |     |         |
     // u |     |        |     |             |     |         |
     // x |     |        |     |             |     |         |
     //   H-----I--------J-----K-------------L-----M---------N
     //                    noflux
     //

     int iA = 0;
     int iB = static_cast<int>  (0.101 * static_cast<double>(ixMax-1));
     int iC = static_cast<int>  (0.301 * static_cast<double>(ixMax-1));
     int iD = static_cast<int>  (0.401 * static_cast<double>(ixMax-1));
     int iE = static_cast<int>  (0.701 * static_cast<double>(ixMax-1));
     int iF = static_cast<int>  (0.801 * static_cast<double>(ixMax-1));
     int iG = ixMax-1;

     numLabels = 5; // SI, EMITTER, BASE, COLLECTOR, NOFLUX

    // do the emitter .  This goes from A to B, inclusive of both A and B.
    xLabel.mNodeVector.clear();
    xLabel.name = "EMITTER";
    xLabel.iIndex = 1;
    xLabel.uType = TYPE_EDGE;
    xLabel.cNode = iB-iA+1;
    xLabel.mNodeVector.reserve(xLabel.cNode+1);

    // fill in all the nodes:
    nodeTmp = 0;
    j=iyMax-1;
    for (i=iA;i<iB+1;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;

    mLabelVector.push_back(xLabel);

    // do the base. This goes from C to D, inclusive of both C and D.
    xLabel.mNodeVector.clear();
    xLabel.name = "BASE";
    xLabel.iIndex = 2;
    xLabel.uType = TYPE_EDGE;
    xLabel.cNode = iD-iC+1;
    xLabel.mNodeVector.reserve(xLabel.cNode+1);

    // fill in all the nodes:
    nodeTmp = 0;
    j=iyMax-1;
    for (i=iC;i<iD+1;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;

    mLabelVector.push_back(xLabel);

    // do the collector This goes from E to F, inclusive of both E and F.
    xLabel.mNodeVector.clear();
    xLabel.name = "COLLECTOR";
    xLabel.iIndex = 3;
    xLabel.uType = TYPE_EDGE;
    xLabel.cNode = iF-iE+1;
    xLabel.mNodeVector.reserve(xLabel.cNode+1);

    // fill in all the nodes:
    nodeTmp = 0;
    j=iyMax-1;
    for (i=iE;i<iF+1;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;

    mLabelVector.push_back(xLabel);

    // Now do noflux.  This includes a bunch of stuff:
    // 		B to C  not including b or c
    // 		D to E  not including d or e
    // 		F to G  not including f or g
    // 		G to N  including both
    // 		N to H  including both
    // 		H to A  including H, not including A.

    xLabel.mNodeVector.clear();
    xLabel.name = "NOFLUX";
    xLabel.iIndex = 4;
    xLabel.uType = TYPE_EDGE;
    //xLabel.cNode = (iC-iB-1)+(iE-iD-1)+(iG-iF-1)+(iyMax)+(iyMax)+(ixMax-1);
    xLabel.mNodeVector.reserve(2*ixMax + 2*iyMax); // extra to be safe.

    // fill in all the nodes:
    nodeTmp = 0;
    j=iyMax-1;
    for (i=iB+1;i<iC;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    for (i=iD+1;i<iE;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    for (i=iF+1;i<iG;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    i=ixMax-1;
    for (j=0;j<iyMax-1;++j)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    i=0;
    for (j=0;j<iyMax-2;++j)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    j=0;
    for (i=0;i<ixMax-1;++i)
    {
      int inode = nodeIndices[i][j];
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
    }

    xLabel.cNode = xLabel.mNodeVector.size ();

    tmpName = xLabel.name;
    tmpName.toUpper();
    mLabelMap[tmpName] = xLabel;

    mLabelVector.push_back(xLabel);


    // Now assign labels to the edges in mEdgeVector.
    // Note:  labels were already assigned to the cells in mCellVector, up
    //        in setupInternalMesh.
    //        Also, nodes don't get labels.

    j=iyMax-1;
    for (i=iA+1;i<iB+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 1; // EMITTER index
    }

    j=iyMax-1;
    for (i=iC+1;i<iD+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 2; // BASE index
    }

    j=iyMax-1;
    for (i=iE+1;i<iF+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 3; // COLLECTOR index
    }


    j=iyMax-1;
    for (i=iB+1;i<iC+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // NOFLUX index
    }

    j=iyMax-1;
    for (i=iD+1;i<iE+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // NOFLUX index
    }

    j=iyMax-1;
    for (i=iF+1;i<iG+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // NOFLUX index
    }

    j=0;
    for (i=1;i<ixMax;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // NOFLUX index
    }

    i=0;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i][j-1];
      int n2 = nodeIndices[i][j  ];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // NOFLUX index
    }

    i=ixMax-1;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i][j-1];
      int n2 = nodeIndices[i][j  ];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // NOFLUX index
    }

  }
  else if (numberElectrodes==4)
  {
     // Assuming a MOSFET mesh:
     //
     //   Ratios:  AB = 2/10
     //            BC = 1/10
     //            CD = 4/10
     //            DE = 1/10
     //            EF = 2/10
     //
     //   So, hopefully, the number of mesh cells is
     //   easily divisible by 10.  The default mesh spacing will work.
     //
     //   Note:  For now, this "mosfet" doesn't include an oxide.
     //          Later, it probably will include one. ERK
     //          11/12/02.
     //
     //     Source              Gate               Drain
     //
     //    |--5--|   |------------3------------|   |--5--|
     //
     // -  A-----B---C-------------------------D---E-----F -
     // 6  |         |                         |         | |
     // -  |         I-------------------------J         | |
     // n  |                                             | |Noflux
     // o  |                                             | 2
     // f  |                                             | |
     // l  |                                             | |
     // u  |                                             | |
     // x  H---------------------------------------------G -
     //
     //    |----------------------1----------------------|
     //                          Bulk
     //
     int iA = 0;
     int iB = static_cast<int>  (0.101 * static_cast<double>(ixMax-1));
     int iC = static_cast<int>  (0.151 * static_cast<double>(ixMax-1));
     int iD = static_cast<int>  (0.851 * static_cast<double>(ixMax-1));
     int iE = static_cast<int>  (0.901 * static_cast<double>(ixMax-1));
     int iF = ixMax-1;

     numLabels = 6; // SI, SOURCE, GATE, DRAIN, SUB, NOFLUX

     if (DEBUG_DEVICE)
     {
       Xyce::dout() << "iA = "<< iA << std::endl;
     Xyce::dout() << "iB = "<< iB << std::endl;
     Xyce::dout() << "iC = "<< iC << std::endl;
     Xyce::dout() << "iD = "<< iD << std::endl;
     Xyce::dout() << "iE = "<< iE << std::endl;
     Xyce::dout() << "iF = "<< iF << std::endl;
     }

     // do the source.  This goes from A to B, inclusive of both A and B.
     xLabel.mNodeVector.clear();
     xLabel.name = "SOURCE";
     xLabel.iIndex = 1;
     xLabel.uType = TYPE_EDGE;
     xLabel.cNode = iB-iA+1;
     xLabel.mNodeVector.reserve(xLabel.cNode+1);

     // fill in all the nodes:
     nodeTmp = 0;
     j=iyMax-1;
     for (i=iA;i<iB+1;++i)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     tmpName = xLabel.name;
     tmpName.toUpper();
     mLabelMap[tmpName] = xLabel;

     mLabelVector.push_back(xLabel);

     // do the gate. This goes from C to D, inclusive of both C and D.
     xLabel.mNodeVector.clear();
     xLabel.name = "GATE";
     xLabel.iIndex = 2;
     xLabel.uType = TYPE_EDGE;
     xLabel.cNode = iD-iC+1;
     xLabel.mNodeVector.reserve(xLabel.cNode+1);

     // fill in all the nodes:
     nodeTmp = 0;
     j=iyMax-1;
     for (i=iC;i<iD+1;++i)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     tmpName = xLabel.name;
     tmpName.toUpper();
     mLabelMap[tmpName] = xLabel;

     mLabelVector.push_back(xLabel);

     // do the drain. This goes from E to F, inclusive of both E and F.
     xLabel.mNodeVector.clear();
     xLabel.name = "DRAIN";
     xLabel.iIndex = 3;
     xLabel.uType = TYPE_EDGE;
     xLabel.cNode = iF-iE+1;
     xLabel.mNodeVector.reserve(xLabel.cNode+1);

     // fill in all the nodes:
     nodeTmp = 0;
     j=iyMax-1;
     for (i=iE;i<iF+1;++i)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     tmpName = xLabel.name;
     tmpName.toUpper();
     mLabelMap[tmpName] = xLabel;

     mLabelVector.push_back(xLabel);

     // do the bulk.  This goes from G to H, inclusive of both G and H.
     // (note that this is the same x indices as A to F).
     xLabel.mNodeVector.clear();
     xLabel.name = "SUB";
     xLabel.iIndex = 4;
     xLabel.uType = TYPE_EDGE;
     xLabel.cNode = iF-iA+1;
     xLabel.mNodeVector.reserve(xLabel.cNode+1);

     // fill in all the nodes:
     nodeTmp = 0;
     j=0;
     for (i=iA;i<iF+1;++i)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     tmpName = xLabel.name;
     tmpName.toUpper();
     mLabelMap[tmpName] = xLabel;

     mLabelVector.push_back(xLabel);

     // do the noflux. Everything that has not been labeled yet,
     // gets this one.  BC, DE, FG, HA, not including
     // the endpoints.

     xLabel.mNodeVector.clear();
     xLabel.name = "NOFLUX";
     xLabel.iIndex = 5;
     xLabel.uType = TYPE_EDGE;
     xLabel.cNode = (iC-iB-1)+(iE-iD-1)+(iyMax-1)+(iyMax-1);
     xLabel.mNodeVector.reserve(xLabel.cNode+1);

     // fill in all the nodes:
     nodeTmp = 0;
     j=iyMax-1;
     for (i=iB+1;i<iC;++i)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     for (i=iD+1;i<iE;++i)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     i=ixMax-1;
     for (j=1;j<iyMax-1;++j)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     i=0;
     for (j=1;j<iyMax-1;++j)
     {
       int inode = nodeIndices[i][j];
       xLabel.mNodeVector.push_back(inode);
       ++nodeTmp;
     }

     tmpName = xLabel.name;
     tmpName.toUpper();
     mLabelMap[tmpName] = xLabel;

     mLabelVector.push_back(xLabel);

     // Now assign labels to the edges in mEdgeVector.
     // Note:  labels were already assigned to the cells in mCellVector, up
     //        in setupInternalMesh.
     //        Also, nodes don't get labels.

    j=iyMax-1;
    for (i=iA+1;i<iB+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 1; // SOURCE index
    }

    j=iyMax-1;
    for (i=iC+1;i<iD+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 2; // GATE index
    }

    j=iyMax-1;
    for (i=iE+1;i<iF+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 3; // DRAIN index
    }

    j=0;
    for (i=iA+1;i<iF+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 4; // SUB index
    }

    j=iyMax-1;
    for (i=iB+1;i<iC+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 5; // NOFLUX index
    }

    j=iyMax-1;
    for (i=iD+1;i<iE+1;++i)
    {
      int n1 = nodeIndices[i-1][j];
      int n2 = nodeIndices[i  ][j];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 5; // NOFLUX index
    }

    i=0;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i][j-1];
      int n2 = nodeIndices[i][j  ];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 5; // NOFLUX index
    }

    i=ixMax-1;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i][j-1];
      int n2 = nodeIndices[i][j  ];

      int e1 = edgeIndices[n1][n2];

      mEdgeVector[e1].uLabel = 5; // NOFLUX index
    }
  }
  else
  {
    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "Sorry, the internal mesh generator can't";
    Xyce::dout() << "  handle anything greater than 4 electrodes." << std::endl;
    }
  }

  return true;

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::setupInternalLabels
//
// Purpose       : This function sets up labels, based on user-specified
//                 information.  Unlike the labels that are set up by the
//                 function setupDefaultLabels, these label sets are not
//                 internally hardwired.
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::setupInternalLabels(int numberElectrodes,
	    std::map<std::string,PDE_2DElectrode*> & elMap)
{

  // Region labels of nodes-------------------------------------------------
  // set up the silicon region label.  This will always be here.
  int i,j;
  mLabel xLabel;
  ExtendedString tmpName = "SI";
  int nodeTmp = 0;

  xLabel.name = "SI";
  xLabel.iIndex = 0;
  xLabel.uType = TYPE_REGION;

  // This assumes every node in the mesh is SI, except for the edges.
  xLabel.mNodeVector.reserve(ixMax*iyMax); // just to be safe, reserve extra

  // fill in all the nodes:
  nodeTmp = 0;
  for (i=1;i<ixMax-1;++i)
  {
    for (j=1;j<iyMax-1;++j)
    {
      xLabel.mNodeVector.push_back(nodeIndices[i][j]);
      ++nodeTmp;
    }
  }
  xLabel.cNode = xLabel.mNodeVector.size ();

  tmpName = xLabel.name;
  tmpName.toUpper();
  mLabelMap[tmpName] = xLabel;

  mLabelVector.push_back(xLabel);

  // Edge labels of nodes----------------------------------------------------
  // Now do the edge labels, using the electrode Map.
  std::vector<int> done;
  done.resize(numNodes, 0);

  std::map<std::string, PDE_2DElectrode*>::iterator mapIter;
  std::map<std::string, PDE_2DElectrode*>::iterator mapStart =
    elMap.begin();
  std::map<std::string, PDE_2DElectrode*>::iterator mapEnd =
    elMap.end();

  numLabels = 1; // SI so far...
  int labelIndex = 0;
  ++labelIndex;

  for ( mapIter = mapStart; mapIter != mapEnd; ++mapIter )
  {
    PDE_2DElectrode & el = *(mapIter->second);

    // determine if this is an x-side or a y-side, then do
    // appropriate stuff:
    if (el.side == "top" || el.side == "bottom")
    {
      double start,end;
      if (el.start <= el.end)
      {
        start = el.start;
        end   = el.end;
      }
      else
      {
        start = el.end;
        end   = el.start;
      }

      // Note: This assumes that the lower-left-hand corner of the mesh is
      // at the (0.0, 0.0) origin!
      if (start < xMin) start = xMin;
      if (end   < xMin) end   = xMin;

      if (start > xMax) start = xMax;
      if (end   > xMax) end   = xMax;

      // Now make start and end relative quantities, rather than
      // absolute:
      start = (start-xMin)/(xMax-xMin);
      end   = (end  -xMin)/(xMax-xMin);

      el.iA = static_cast<int>  (start * static_cast<double>(ixMax-1));
      el.iB = static_cast<int>  (end   * static_cast<double>(ixMax-1));

      if (el.iA < 0)       el.iA = 0;
      if (el.iA > ixMax-1) el.iA = ixMax-1;
      if (el.iB < 0)       el.iB = 0;
      if (el.iB > ixMax-1) el.iB = ixMax-1;

      ++numLabels;

      ExtendedString tmpName = el.name;
      tmpName.toUpper ();

      // do the nodes.  This goes from A to B, inclusive of both A and B.
      xLabel.mNodeVector.clear();
      xLabel.name = tmpName;
      xLabel.iIndex = labelIndex; ++labelIndex;
      xLabel.uType = TYPE_EDGE;
      xLabel.cNode = el.iB-el.iA+1;
      xLabel.mNodeVector.reserve(xLabel.cNode+1);

      el.uLabel = xLabel.iIndex;

      // fill in all the nodes:
      nodeTmp = 0;
      if (el.side == "top")
      {
        j=iyMax-1;
      }
      else if (el.side == "bottom")
      {
        j=0;
      }

      for (i=el.iA;i<el.iB+1;++i)
      {
        int inode = nodeIndices[i][j];
        xLabel.mNodeVector.push_back(inode);
        done[inode] = 1;
        ++nodeTmp;
      }

      tmpName = xLabel.name;
      tmpName.toUpper();
      mLabelMap[tmpName] = xLabel;

      mLabelVector.push_back(xLabel);
    }
    else if (el.side == "right" || el.side == "left")
    {
      double start,end;
      if (el.start <= el.end)
      {
        start = el.start;
        end   = el.end;
      }
      else
      {
        start = el.end;
        end   = el.start;
      }

      // Note: This assumes that the lower-left-hand corner of the mesh is
      // at the (0.0, 0.0) origin!
      if (start < yMin) start = yMin;
      if (end   < yMin) end   = yMin;

      if (start > yMax) start = yMax;
      if (end   > yMax) end   = yMax;

      // Now make start and end relative quantities, rather than
      // absolute:
      start = (start-yMin)/(yMax-yMin);
      end   = (end  -yMin)/(yMax-yMin);

      el.iA = static_cast<int>  (start * static_cast<double>(iyMax-1));
      el.iB = static_cast<int>  (end   * static_cast<double>(iyMax-1));

      if (el.iA < 0)       el.iA = 0;
      if (el.iA > iyMax-1) el.iA = iyMax-1;
      if (el.iB < 0)       el.iB = 0;
      if (el.iB > iyMax-1) el.iB = iyMax-1;

      ++numLabels;

      ExtendedString tmpName = el.name;
      tmpName.toUpper ();

      // do the nodes.  This goes from A to B, inclusive of both A and B.
      xLabel.mNodeVector.clear();
      xLabel.name = tmpName;
      xLabel.iIndex = labelIndex; ++labelIndex;
      xLabel.uType = TYPE_EDGE;
      xLabel.cNode = el.iB-el.iA+1;
      xLabel.mNodeVector.reserve(xLabel.cNode+1);

      el.uLabel = xLabel.iIndex;

      // fill in all the nodes:
      nodeTmp = 0;
      if (el.side == "right")
      {
        i=ixMax-1;
      }
      else if (el.side == "left")
      {
        i=0;
      }

      for (j=el.iA;j<el.iB+1;++j)
      {
        int inode = nodeIndices[i][j];
        xLabel.mNodeVector.push_back(inode);
        done[inode] = 1;
        ++nodeTmp;
      }

      tmpName = xLabel.name;
      tmpName.toUpper();
      mLabelMap[tmpName] = xLabel;

      mLabelVector.push_back(xLabel);
    }
    else
    {
      Report::UserFatal() << "Electrode side specification not recognized.";
    }
  } // end of map loop:

  // now do the NOFLUX label.  This will include all external edges/nodes
  // which have not already been assigned a label.
  ++numLabels;

  xLabel.mNodeVector.clear();
  xLabel.name = "NOFLUX";
  xLabel.iIndex = labelIndex; ++labelIndex;
  xLabel.uType = TYPE_EDGE;
  xLabel.mNodeVector.reserve(2*ixMax + 2*iyMax); // extra to be safe.

  int NOFLUXIndex = xLabel.iIndex;

  nodeTmp = 0;
  i=0;
  for (j=0;j<iyMax;++j)
  {
    int inode = nodeIndices[i][j];
    if (done[inode]!=1)
    {
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
      done[inode] = 1;
    }
  }

  i=ixMax-1;
  for (j=0;j<iyMax;++j)
  {
    int inode = nodeIndices[i][j];
    if (done[inode]!=1)
    {
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
      done[inode] = 1;
    }
  }

  j=0;
  for (i=0;i<ixMax;++i)
  {
    int inode = nodeIndices[i][j];
    if (done[inode]!=1)
    {
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
      done[inode] = 1;
    }
  }

  j=iyMax-1;
  for (i=0;i<ixMax;++i)
  {
    int inode = nodeIndices[i][j];
    if (done[inode]!=1)
    {
      xLabel.mNodeVector.push_back(inode);
      ++nodeTmp;
      done[inode] = 1;
    }
  }

  xLabel.cNode = xLabel.mNodeVector.size ();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "NOFLUX nodes:" << std::endl;
    for (int nodeIndex=0;nodeIndex<xLabel.mNodeVector.size();++nodeIndex)
    {
      int iNode = xLabel.mNodeVector[nodeIndex];
      Xyce::dout() << "node="<<iNode<<" x=" << xVector[iNode];
      Xyce::dout() << " y="<<yVector[iNode] <<std::endl;
    }
  }

  tmpName = xLabel.name;
  tmpName.toUpper();
  mLabelMap[tmpName] = xLabel;

  mLabelVector.push_back(xLabel);

  // Labels of edges------------------------------------------------------
  // Note:  labels were already assigned to the cells in mCellVector, up
  //        in setupInternalMesh.
  //        Also, nodes don't get labels.
  std::vector<int> edgeDone;
  edgeDone.resize(numEdges, 0);

  for ( mapIter = mapStart; mapIter != mapEnd; ++mapIter )
  {
    PDE_2DElectrode & el = *(mapIter->second);

    if (el.side=="top")
    {
      j=iyMax-1;
      for (i=el.iA+1;i<el.iB+1;++i)
      {
        int n1 = nodeIndices[i-1][j];
        int n2 = nodeIndices[i  ][j];

        int e1 = edgeIndices[n1][n2];

        mEdgeVector[e1].uLabel = el.uLabel;
        edgeDone[e1] = 1;
      }
    }
    else if (el.side=="bottom")
    {
      j=0;
      for (i=el.iA+1;i<el.iB+1;++i)
      {
        int n1 = nodeIndices[i-1][j];
        int n2 = nodeIndices[i  ][j];

        int e1 = edgeIndices[n1][n2];

        mEdgeVector[e1].uLabel = el.uLabel;
        edgeDone[e1] = 1;
      }
    }
    else if (el.side=="left")
    {
      i=0;
      for (j=el.iA+1;j<el.iB+1;++j)
      {
        int n1 = nodeIndices[i][j-1];
        int n2 = nodeIndices[i][j  ];

        int e1 = edgeIndices[n1][n2];

        mEdgeVector[e1].uLabel = el.uLabel;
        edgeDone[e1] = 1;
      }
    }
    else if (el.side=="right")
    {
      i=ixMax-1;
      for (j=el.iA+1;j<el.iB+1;++j)
      {
        int n1 = nodeIndices[i][j-1];
        int n2 = nodeIndices[i][j  ];

        int e1 = edgeIndices[n1][n2];

        mEdgeVector[e1].uLabel = el.uLabel;
        edgeDone[e1] = 1;
      }
    }
    else
    {
      Report::UserFatal() << "Electrode side specification not recognized.";
    }
  }

  // Now do the NOFLUX edges:
  i=0;
  for (j=1;j<iyMax;++j)
  {
    int n1 = nodeIndices[i ][j  ];
    int n2 = nodeIndices[i ][j-1];
    int e1 = edgeIndices[n1][n2 ];

    if (edgeDone[e1]!=1)
    {
      mEdgeVector[e1].uLabel = NOFLUXIndex;
      edgeDone[e1] = 1;
    }
  }

  i=ixMax-1;
  for (j=1;j<iyMax;++j)
  {
    int n1 = nodeIndices[i ][j  ];
    int n2 = nodeIndices[i ][j-1];
    int e1 = edgeIndices[n1][n2 ];

    if (edgeDone[e1]!=1)
    {
      mEdgeVector[e1].uLabel = NOFLUXIndex;
      edgeDone[e1] = 1;
    }
  }

  j=0;
  for (i=1;i<ixMax;++i)
  {
    int n1 = nodeIndices[i  ][j ];
    int n2 = nodeIndices[i-1][j ];
    int e1 = edgeIndices[n1 ][n2];

    if (edgeDone[e1]!=1)
    {
      mEdgeVector[e1].uLabel = NOFLUXIndex;
      edgeDone[e1] = 1;
    }
  }

  j=iyMax-1;
  for (i=1;i<ixMax;++i)
  {
    int n1 = nodeIndices[i  ][j ];
    int n2 = nodeIndices[i-1][j ];
    int e1 = edgeIndices[n1 ][n2];

    if (edgeDone[e1]!=1)
    {
      mEdgeVector[e1].uLabel = NOFLUXIndex;
      edgeDone[e1] = 1;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "EDGEs:"<<std::endl;
    j=0;
    for (i=1;i<ixMax;++i)
    {
      int n1 = nodeIndices[i  ][j ];
      int n2 = nodeIndices[i-1][j ];
      int e1 = edgeIndices[n1 ][n2];

      Xyce::dout() << "edge: " << e1;
      Xyce::dout() << " x(A)="<<xVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " y(A)="<<yVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " x(B)="<<xVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " y(B)="<<yVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " label = " << mEdgeVector[e1].uLabel;
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;

    j=iyMax-1;
    for (i=1;i<ixMax;++i)
    {
      int n1 = nodeIndices[i  ][j ];
      int n2 = nodeIndices[i-1][j ];
      int e1 = edgeIndices[n1 ][n2];

      Xyce::dout() << "edge: " << e1;
      Xyce::dout() << " x(A)="<<xVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " y(A)="<<yVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " x(B)="<<xVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " y(B)="<<yVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " label = " << mEdgeVector[e1].uLabel;
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;

    i=0;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i  ][j  ];
      int n2 = nodeIndices[i  ][j-1];
      int e1 = edgeIndices[n1 ][n2 ];

      Xyce::dout() << "edge: " << e1;
      Xyce::dout() << " x(A)="<<xVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " y(A)="<<yVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " x(B)="<<xVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " y(B)="<<yVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " label = " << mEdgeVector[e1].uLabel;
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;

    i=ixMax-1;
    for (j=1;j<iyMax;++j)
    {
      int n1 = nodeIndices[i  ][j  ];
      int n2 = nodeIndices[i  ][j-1];
      int e1 = edgeIndices[n1 ][n2 ];

      Xyce::dout() << "edge: " << e1;
      Xyce::dout() << " x(A)="<<xVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " y(A)="<<yVector[mEdgeVector[e1].inodeA];
      Xyce::dout() << " x(B)="<<xVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " y(B)="<<yVector[mEdgeVector[e1].inodeB];
      Xyce::dout() << " label = " << mEdgeVector[e1].uLabel;
      Xyce::dout() << std::endl;
    }
    Xyce::dout() << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::setupGeometry
// Purpose       : This function sets up a lot of basic GLOBAL mesh geometry
//                 information.  For example it calculates total volume of
//                 the mesh, but not volumes around individual nodes.
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::setupGeometry ()
{

  // Information from a node loop:
  // Obtain:
  //   -  the maximum number of nearest neighbors for a single node.
  //   -  the total area(volume) of the mesh.  (note: this is a 2D mesh, so "volume"
  //      is clear-cut if dealing with cylindrical coord., but less so for cartesian.)
  //   -  the total surface area (surface length?) of the mesh. (not done yet)
  //
  std::vector<mNode>::iterator first = mNodeVector.begin ();
  std::vector<mNode>::iterator last  = mNodeVector.end   ();
  std::vector<mNode>::iterator iter;

  maxNodeNN = -999;
  vol = 0.0;

  for (iter=first;iter!=last;++iter)
  {
    if (maxNodeNN < iter->cnode) maxNodeNN = iter->cnode;
    vol += iter->area;

  }
  if (!cylGeom) vol *= depth;


  // loop over each region (lists of nodes associated with labels which have
  // been designated as "region" labels ) and obtain:
  //   - the total area(volume) of each region.
  //   - the surface length(area) of each region. (not done yet)

  std::map<std::string, mLabel>::iterator firstL = mLabelMap.begin ();
  std::map<std::string, mLabel>::iterator lastL  = mLabelMap.end   ();
  std::map<std::string, mLabel>::iterator iterL;

  for (iterL=firstL; iterL!=lastL; ++iterL)
  {
    if (iterL->second.uType == TYPE_EDGE) continue;

    std::vector<int>::iterator firstN = iterL->second.mNodeVector.begin ();
    std::vector<int>::iterator lastN  = iterL->second.mNodeVector.end   ();
    std::vector<int>::iterator iterN;

    iterL->second.vol = 0.0;

    for (iterN=firstN; iterN!=lastN; ++iterN)
    {
      iterL->second.vol += mNodeVector[*iterN].area;
    }
    if (!cylGeom) iterL->second.vol *= depth;
  }

  // calculate xMax and yMax (could also be called rMax and zMax)

  for (iter=first;iter!=last;++iter)
  {
    if (xMax < iter->x) xMax = iter->x;
    if (yMax < iter->y) yMax = iter->y;
    if (xMin > iter->x) xMin = iter->x;
    if (yMin > iter->y) yMin = iter->y;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::dumpMesh
// Purpose       : This function dumps the current mesh stored in memory to
//                 a text file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::dumpMesh ()
{
  // open output file
  std::string dbFileName = meshFileName.substr(0, meshFileName.size() - 4) + "mesh_debug.txt";
  FILE *file = fopen(dbFileName.c_str(), "wt");

  fprintf(file, "CELLS\n");
  fprintf(file, "      |     |          EDGES                     CELLS        \n");
  fprintf(file, "  i   |  u  |   AB    BC  AC/CD   AD     AB    BC  AC/CD   AD \n");
  fprintf(file, "----- | --- | ----- ----- ----- -----  ----- ----- ----- -----\n");

  for(int i = 0; i < numCells; ++i)
  {
    mCell c1 = mCellVector[i];

    fprintf(file, "%5u | ", i);
    if (c1.uLabel != -1 ) fprintf(file, "%3u | ", c1.uLabel);
    else fprintf(file, "    | ");

    if (c1.iedgeAB != -1 ) fprintf(file, "%5u ", c1.iedgeAB);
    else fprintf(file, "      ");

    if (c1.iedgeBC != -1 ) fprintf(file, "%5u ", c1.iedgeBC);
    else fprintf(file, "      ");

    if (c1.iedgeCD != -1 ) fprintf(file, "%5u ", c1.iedgeCD);
    else fprintf(file, "      ");

    if (c1.iedgeDA != -1 ) fprintf(file, "%5u ", c1.iedgeDA);
    else fprintf(file, "  n/a ");

    fprintf(file, " ");
    if (c1.icellAB != -1 ) fprintf(file, "%5u ", c1.icellAB);
    else fprintf(file, "      ");

    if (c1.icellBC != -1 ) fprintf(file, "%5u ", c1.icellBC);
    else fprintf(file, "      ");

    if (c1.icellCD != -1 ) fprintf(file, "%5u ", c1.icellCD);
    else fprintf(file, "      ");

    if (c1.icellDA != -1 ) fprintf(file, "%5u" , c1.icellDA);
    else fprintf(file, "  n/a");

    fprintf(file, "\n");
  }
  fprintf(file, "----- | --- | ----- ----- ----- -----  ----- ----- ----- -----\n");

  fprintf(file, "\nEDGES\n");
  fprintf(file, "      |     |    NODES    |              |              |              |              |\n");
  fprintf(file, "  i   |  u  |   A     B   |     ilen     |      elen    |   Area1      |   Area2      |\n");
  fprintf(file, "----- | --- | ----- ----- | ------------ | ------------ | ------------ | ------------ |\n");

  for(int i = 0; i < numEdges; ++i)
  {
    mEdge e1 = mEdgeVector[i];
    fprintf(file, "%5u | ", i);
    if (e1.uLabel != -1 ) fprintf(file, "%3u | ", e1.uLabel);
    else fprintf(file, "    | ");
    fprintf(file, "%5u %5u | %12.4e | %12.4e | %12.4e | %12.4e |\n",
                e1.inodeA, e1.inodeB, e1.ilen, e1.elen, e1.Area1, e1.Area2);
  }

  fprintf(file, "\nNODES\n");
  fprintf(file, "      |                           |              |\n");
  fprintf(file, "  i   |       x           y       |      Area    |\n");
  fprintf(file, "----- | ------------ ------------ | ------------ |\n");

  for(int i = 0; i < numNodes; ++i)
  {
    mNode n1 = mNodeVector[i];
    fprintf(file, "%5u | ", i);
    fprintf(file, "%12.4e %12.4e | %12.4e |\n", n1.x, n1.y, n1.area);
  }

  fprintf(file, "\nNodeEdgeInfo\n");
  for(int i = 0; i < numNodes; ++i)
  {
    fprintf(file, "------\n");
    fprintf(file, "  node index = %d\n", i);
    for (int j=0;j<mNodeVector[i].cnode;++j)
    {
      fprintf(file," local  edge index = %d\n", j);

      fprintf(file," global edge index = %d\n",
        mNodeVector[i].edgeInfoVector[j].iedge);

      fprintf(file," neighbor node     = %d\n",
        mNodeVector[i].edgeInfoVector[j].inode);

      fprintf(file," elen              = %12.4e\n",
        mNodeVector[i].edgeInfoVector[j].elen );

      fprintf(file," ilen              = %12.4e\n",
        mNodeVector[i].edgeInfoVector[j].ilen );
    }
    fprintf(file, "------\n");
  }

  fclose(file);

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::printLabels
// Purpose       : Prints the list of mesh labels (both cell and edge) to
//                 stdout.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::printLabels ()
{
  int i;

  std::map<std::string,mLabel>::iterator first = mLabelMap.begin ();
  std::map<std::string,mLabel>::iterator last  = mLabelMap.end   ();
  std::map<std::string,mLabel>::iterator iter  = mLabelMap.end   ();

  Xyce::dout() << std::endl;
  Xyce::dout() << "Mesh Labels:" <<std::endl;
  Xyce::dout() << "   Index   # nodes      Type   Label";
  Xyce::dout() << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    Xyce::dout().width(8);
    Xyce::dout() << iter->second.iIndex;
    Xyce::dout().width(10);
    Xyce::dout() << iter->second.cNode;

    if (iter->second.uType == TYPE_EDGE)  Xyce::dout() << "  Edge    ";
    else                                  Xyce::dout() << "  Region  ";

    Xyce::dout() << "   ";
    Xyce::dout().width(15);
    Xyce::dout() << iter->second.name << std::endl;
  }
  Xyce::dout() << std::endl;

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::outputMeshInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::outputMeshInfo ()
{

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::computeIntPB
// Purpose       : This function computes the intersection of the
//                 perpendicular bisection of the triangle given by nodes
//                 A, B and C.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::computeIntPB
   (double &x, double &y,int inodeA,int inodeB,int inodeC)
{

  mNode  nA = mNodeVector[inodeA];
  mNode  nB = mNodeVector[inodeB];
  mNode  nC = mNodeVector[inodeC];

  double x1 = nA.x;  double y1 = nA.y;
  double x2 = nB.x;  double y2 = nB.y;
  double x3 = nC.x;  double y3 = nC.y;

  double dx12 = x1 - x2;
  double dx23 = x2 - x3;
  double dx13 = x1 - x3;

  double m12 = (dx12) ? (y1 - y2) / dx12 : 0.0;
  double m23 = (dx23) ? (y2 - y3) / dx23 : 0.0;
  double m13 = (dx13) ? (y1 - y3) / dx13 : 0.0;
  double x12 = (x1 + x2) / 2.0;  double y12 = (y1 + y2) / 2.0;
  double x23 = (x2 + x3) / 2.0;  double y23 = (y2 + y3) / 2.0;
  double x13 = (x1 + x3) / 2.0;  double y13 = (y1 + y3) / 2.0;

  UINT iSmallest;
  if (fabs(dx12) < fabs(dx23))
  {
    iSmallest = (fabs(dx12) < fabs(dx13)) ? TAG12 : TAG13;
  }
  else
  {
    iSmallest = (fabs(dx23) < fabs(dx13)) ? TAG23 : TAG13;
  }

  switch (iSmallest)
  {
    case TAG12:
      x = (m13 * m23 * (y13 - y23) + m23 * x13 - m13 * x23) / (m23 - m13);
      break;

    case TAG23:
      x = (m13 * m12 * (y13 - y12) + m12 * x13 - m13 * x12) / (m12 - m13);
      break;

    case TAG13:
      x = (m12 * m23 * (y12 - y23) + m23 * x12 - m12 * x23) / (m23 - m12);
      break;
  }

  UINT iLargest;
  if (fabs(m12) > fabs(m23))
  {
    iLargest = (fabs(m12) > fabs(m13)) ? TAG12 : TAG13;
  }
  else
  {
    iLargest = (fabs(m23) > fabs(m13)) ? TAG23 : TAG13;
  }

  switch (iLargest)
  {
    case TAG12:
      y = (x12 - x) / m12 + y12;
      break;

    case TAG23:
      y = (x23 - x) / m23 + y23;
      break;

    case TAG13:
      y = (x13 - x) / m13 + y13;
      break;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In computeIntPB:\n";
    Xyce::dout() << " inodeA = " << inodeA;
    Xyce::dout() << " inodeB = " << inodeB;
    Xyce::dout() << " inodeC = " << inodeC << std::endl;
    Xyce::dout() << " x1 = " << x1;
    Xyce::dout() << " y1 = " << y1 <<std::endl;

    Xyce::dout() << " x2 = " << x2;
    Xyce::dout() << " y2 = " << y2 <<std::endl;

    Xyce::dout() << " x3 = " << x3;
    Xyce::dout() << " y3 = " << y3 <<std::endl;

    Xyce::dout() << " iSmallest = " << iSmallest << std::endl;
    Xyce::dout() << " iLargest  = " << iLargest  << std::endl;

    Xyce::dout() << " m12 = " << m12 << std::endl;
    Xyce::dout() << " m23 = " << m23 << std::endl;
    Xyce::dout() << " m13 = " << m13 << std::endl;

    Xyce::dout() << " x   = " << x << std::endl;
    Xyce::dout() << " y   = " << y << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function:     : PDE_2DMesh::lengthAdjust
// Purpose:      : This function returns a length between two points
//                 defined by (x1,y1) and (x2,y2).
//
//                 The length is adjusted cylinderical geometries.
//                 (otherwise the length would just be sqrt(dx*dx + dy*dy)).
//
//                 The integration segment
//                 becomes an integration surface by rotating the
//                 integration segment around the y axis.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
double PDE_2DMesh::lengthAdjust
    (double x1, double y1, double x2, double y2)
{
  double dx1 = x2 - x1;
  double dy1 = y2 - y1;
  double h  = sqrt(dx1 * dx1 + dy1 * dy1);
  double pi = M_PI;
  double A  = pi * (x1 + x2) * h;
  return A;
}


//-----------------------------------------------------------------------------
// Function:     : PDE_2DMesh::areaAdjust
//
//
// Purpose:      : This function calculates an area which has been adjusted
//                 for cylinderical geometries.  The surfaces become
//                 regions by rotating the surfaces around the y axis.
//
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
double PDE_2DMesh::areaAdjust
       (double x1, double y1, double x2, double y2, double x3, double y3)
{
  // sort the coordinates so that x1 <= x2 <= x3
  double xT, yT;
  if (x1 > x2) { xT = x1; yT = y1; x1 = x2; y1 = y2; x2 = xT; y2 = yT; }
  if (x2 > x3) { xT = x2; yT = y2; x2 = x3; y2 = y3; x3 = xT; y3 = yT; }
  if (x1 > x2) { xT = x1; yT = y1; x1 = x2; y1 = y2; x2 = xT; y2 = yT; }

  // calculate the slope of m13 which is never undefined since x1 != x3
  double m13 = (y1 - y3) / (x1 - x3);

  // calculate the square and cube of x2
  double x2_2 = x2 * x2;
  double x2_3 = x2 * x2_2;

  // calculate the volume from integrating from x1 to x2
  double V12 = 0.0;
  if (fabs(x1-x2) > 1.0e-14)  // @@@ need to make this relative
  {
    double m12 = (y1 - y2) / (x1 - x2);
    double x1_2 = x1 * x1;
    double x1_3 = x1 * x1_2;
    V12 = (m12 - m13) * ((x2_3 - x1_3) / 3.0 - x1 * (x2_2 - x1_2) / 2.0);
  }

  // calculate the volume from integrating from x2 to x3
  double V23 = 0.0;
  if (fabs(x2-x3) > 1.0e-14)  // @@@ need to make this relative
  {
    double m23 = (y2 - y3) / (x2 - x3);
    double x3_2 = x3 * x3;
    double x3_3 = x3 * x3_2;
    V23 = (m23 - m13) * ((x3_3 - x2_3) / 3.0 - x3 * (x3_2 - x2_2) / 2.0);
  }

  // calculate the total volume
  double pi = M_PI;
  double V = 2.0 * pi * (fabs(V12) + fabs(V23));
  return V;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::computeAngle
// Description   : This function computes the angle between the three nodes.
//                 The angle is assumed to be less than PI radians.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
double PDE_2DMesh::computeAngle (int inode1,int inode2,int inode3)
{
  mNode n1 = mNodeVector[inode1];
  mNode n2 = mNodeVector[inode2];
  mNode n3 = mNodeVector[inode3];

  double xl = n1.x - n2.x;
  double xr = n3.x - n2.x;
  double yl = n1.y - n2.y;
  double yr = n3.y - n2.y;

  double r = (xl*xr+yl*yr)/(sqrt(sq(xl)+sq(yl))*sqrt(sq(xr)+sq(yr)));
  if      (r >  1.0) r =  1.0;
  else if (r < -1.0) r = -1.0;

  double angle = acos(r);
  //double pi = M_PI;
  // if (xl*yr-xr*yl > 0) angle = 2*pi - angle;
  return angle;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::cellNodes
// Description   : This function determines the nodes of all the cells:
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::cellNodes ()
{

  // set up cell Node lists:
  std::vector<mCell>::iterator first = mCellVector.begin();
  std::vector<mCell>::iterator last  = mCellVector.end ();
  std::vector<mCell>::iterator iter;

  int i;
  for (i=0,iter=first; iter!=last ; ++i, ++iter)
  {
    mCell & cellObj = *iter;

    int iedgeDA   = cellObj.iedgeDA;

    mEdge edgeAB = mEdgeVector[cellObj.iedgeAB];
    mEdge edgeBC = mEdgeVector[cellObj.iedgeBC];
    mEdge edgeCD = mEdgeVector[cellObj.iedgeCD];

    // temporary vector of node indices:
    int inode[9];

    inode[0] = edgeAB.inodeA;
    inode[1] = edgeAB.inodeB;
    inode[2] = edgeBC.inodeA;
    inode[3] = edgeBC.inodeB;
    inode[4] = edgeCD.inodeA;
    inode[5] = edgeCD.inodeB;

    // triangular element
    if (iedgeDA == -1)
    {
      if (inode[0] == inode[2])
      {
        cellObj.mNodeVector[VERTEX_A] = inode[1];
        cellObj.mNodeVector[VERTEX_B] = inode[0];
        cellObj.mNodeVector[VERTEX_C] = inode[3];
        cellObj.inodeA = inode[1];
        cellObj.inodeB = inode[0];
        cellObj.inodeC = inode[3];
      }
      else if (inode[0] == inode[3])
      {
        cellObj.mNodeVector[VERTEX_A] = inode[1];
        cellObj.mNodeVector[VERTEX_B] = inode[0];
        cellObj.mNodeVector[VERTEX_C] = inode[2];
        cellObj.inodeA = inode[1];
        cellObj.inodeB = inode[0];
        cellObj.inodeC = inode[2];
      }
      else if (inode[1] == inode[2])
      {
        cellObj.mNodeVector[VERTEX_A] = inode[0];
        cellObj.mNodeVector[VERTEX_B] = inode[1];
        cellObj.mNodeVector[VERTEX_C] = inode[3];
        cellObj.inodeA = inode[0];
        cellObj.inodeB = inode[1];
        cellObj.inodeC = inode[3];
      }
      else
      {
        cellObj.mNodeVector[VERTEX_A] = inode[0];
        cellObj.mNodeVector[VERTEX_B] = inode[1];
        cellObj.mNodeVector[VERTEX_C] = inode[2];
        cellObj.inodeA = inode[0];
        cellObj.inodeB = inode[1];
        cellObj.inodeC = inode[2];
      }
      cellObj.mNodeVector[VERTEX_D] = -1;
      cellObj.inodeD = -1;
    }
    // rectangular element
    else
    {
      mEdge edgeDA = mEdgeVector[cellObj.iedgeDA];

      inode[6] = edgeDA.inodeA;
      inode[7] = edgeDA.inodeB;

      if ((inode[0] == inode[2]) || (inode[0] == inode[3]))
      {
        cellObj.mNodeVector[VERTEX_A] = inode[1];
        cellObj.mNodeVector[VERTEX_B] = inode[0];
        cellObj.inodeA = inode[1];
        cellObj.inodeB = inode[0];
      }
      else
      {
        cellObj.mNodeVector[VERTEX_A] = inode[0];
        cellObj.mNodeVector[VERTEX_B] = inode[1];
        cellObj.inodeA = inode[0];
        cellObj.inodeB = inode[1];
      }
      if ((inode[4] == inode[2]) || (inode[4] == inode[3]))
      {
        cellObj.mNodeVector[VERTEX_C] = inode[4];
        cellObj.mNodeVector[VERTEX_D] = inode[5];
        cellObj.inodeC = inode[4];
        cellObj.inodeD = inode[5];
      }
      else
      {
        cellObj.mNodeVector[VERTEX_C] = inode[5];
        cellObj.mNodeVector[VERTEX_D] = inode[4];
        cellObj.inodeC = inode[5];
        cellObj.inodeD = inode[4];
      }
    }
  } // end of cell loop.

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::fCCWorder
// Description   : This function determines if a rotation from segment 12 to
//                 segment 13 is in the counter-clockwise direction.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::fCCWorder (int inode1, int inode2, int inode3)
{
  mNode node1 = mNodeVector[inode1];
  mNode node2 = mNodeVector[inode2];
  mNode node3 = mNodeVector[inode3];

  double  x1     = node1.x;
  double  y1     = node1.y;
  double  x2     = node2.x;
  double  y2     = node2.y;
  double  x3     = node3.x;
  double  y3     = node3.y;

  double  x12    = x2 - x1;
  double  y12    = y2 - y1;
  double  x13    = x3 - x1;
  double  y13    = y3 - y1;

  double  d12    = sqrt(x12*x12 + y12*y12);
  double  d13    = sqrt(x13*x13 + y13*y13);

  double  arg12  = x12/d12;
  double  arg13  = x13/d13;

  if (arg12 < -1.0) arg12 = 1.0; else if (arg12 > 1.0) arg12 = 1.0;
  if (arg13 < -1.0) arg13 = 1.0; else if (arg13 > 1.0) arg13 = 1.0;

  double  a12    = acos(arg12);
  double  a13    = acos(arg13);

  double pi = M_PI;

  if (y12 < 0) a12 = 2*pi - a12;
  if (y13 < 0) a13 = 2*pi - a13;

  return (a13 > a12) ? true : false;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::elementNodes
// Purpose       : This function determines the three nodes of a triangle
//                  or the four nodes of a rectangle, whatever the
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::elementNodes (int itri, int *ainode)
{
  int iedgeAB = mCellVector[itri].iedgeAB;
  int inode1 = mEdgeVector[iedgeAB].inodeA;
  int inode2 = mEdgeVector[iedgeAB].inodeB;

  int iedgeBC = mCellVector[itri].iedgeBC;
  int inode3 = mEdgeVector[iedgeBC].inodeA;
  int inode4 = mEdgeVector[iedgeBC].inodeB;

  int iedgeDA = mCellVector[itri].iedgeDA;
  if (iedgeDA == -1) // triangular element
  {
    if (inode1 == inode3)
    {
      ainode[VERTEX_A] = inode2;
      ainode[VERTEX_B] = inode1;
      ainode[VERTEX_C] = inode4;
    }
    else if (inode1 == inode4)
    {
      ainode[VERTEX_A] = inode2;
      ainode[VERTEX_B] = inode1;
      ainode[VERTEX_C] = inode3;
    }
    else if (inode2 == inode3)
    {
      ainode[VERTEX_A] = inode1;
      ainode[VERTEX_B] = inode2;
      ainode[VERTEX_C] = inode4;
    }
    else
    {
      ainode[VERTEX_A] = inode1;
      ainode[VERTEX_B] = inode2;
      ainode[VERTEX_C] = inode3;
    }
    ainode[VERTEX_D] = -1;
  }
  else // rectangular element
  {
    int iedgeCD = mCellVector[itri].iedgeCD;
    int inode5 = mEdgeVector[iedgeCD].inodeA;
    int inode6 = mEdgeVector[iedgeCD].inodeB;

    if ((inode1 == inode3) || (inode1 == inode4))
    {
      ainode[VERTEX_A] = inode2;
      ainode[VERTEX_B] = inode1;
    }
    else
    {
      ainode[VERTEX_A] = inode1;
      ainode[VERTEX_B] = inode2;
    }

    if ((inode5 == inode3) || (inode5 == inode4))
    {
      ainode[VERTEX_C] = inode5;
      ainode[VERTEX_D] = inode6;
    }
    else
    {
      ainode[VERTEX_C] = inode6;
      ainode[VERTEX_D] = inode5;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::getElementInfo
//
// Purpose       : This function initialize the ainode, aiedge, aitri and
//                 auLabel arrays with the elements node indices,
//                 edge indices, adjacent triangle indices and label indices.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::getElementInfo
  (int itri, int *ainode, int *aiedge, int *aitri, int *auLabel)
{
  elementNodes(itri,ainode);

  int iedgeAB, iedgeBC, iedgeCD, iedgeDA;

  iedgeAB = aiedge[0] = mCellVector[itri].iedgeAB;
  iedgeBC = aiedge[1] = mCellVector[itri].iedgeBC;
  iedgeCD = aiedge[2] = mCellVector[itri].iedgeCD;
  iedgeDA = aiedge[3] = mCellVector[itri].iedgeDA;

  aitri[0] = mCellVector[itri].icellAB;
  aitri[1] = mCellVector[itri].icellBC;
  aitri[2] = mCellVector[itri].icellCD;
  aitri[3] = mCellVector[itri].icellDA;

  auLabel[0] = mEdgeVector[iedgeAB].uLabel;
  auLabel[1] = mEdgeVector[iedgeBC].uLabel;
  auLabel[2] = mEdgeVector[iedgeCD].uLabel;
  auLabel[3] = (iedgeDA != -1) ? mEdgeVector[iedgeDA].uLabel : -1;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::initNodeAdjStructure
// Purpose       : This function initializes the static adjacency structure
//		   with the neighbors of the specified node.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::initNodeAdjStructure
  (NADJ &nadj, int itri, int iVertex, int uIntLabel, bool fCW)
{
  int inode;
  int ainode[4];
  int aiedge[4];
  int aitri[4];
  int auLabel[4];

  // vertices to edge map
  int  aauVertices2Edge[4][4] =
  {   //              VERTEX_A  VERTEX_B  VERTEX_C  VERTEX_D
    { /* VERTEX_A */     -1   , EDGE_AB , EDGE_AC , EDGE_AD },
    { /* VERTEX_B */  EDGE_AB ,  -1     , EDGE_BC ,  -1     },
    { /* VERTEX_C */  EDGE_AC , EDGE_BC , -1	  , EDGE_CD },
    { /* VERTEX_D */  EDGE_AD ,  -1     , EDGE_CD , -1	    }
  };

  // initialize variables
  getElementInfo(itri, ainode, aiedge, aitri, auLabel);
  int cVertices = (ainode[VERTEX_D] != -1) ? 4 : 3;
  int cnode   = 0;
  nadj.cnode   = 0;
  nadj.inode   = ainode[iVertex];
  nadj.fBndry  = true;
  nadj.fGotAll = false;


  // determine the direction to rotate
  int fCCW = fCCWorder(ainode[VERTEX_A], ainode[VERTEX_B], ainode[VERTEX_C]);
  int nNextVertex = (fCCW == fCW) ? 1 : -1;
  int i	 = (iVertex + nNextVertex + cVertices) % cVertices;
  int i1 = aauVertices2Edge[iVertex][i];
  nadj.ainode[0] = inode = ainode[i];
  nadj.aiedge[0] = aiedge[i1];
  nadj.aielem[0] = itri;
  bool fBndry = (auLabel[i1] != uIntLabel);

  // loop around and find the nodes
  for(;;)
  { for(i = 0; i < cVertices; ++i)
    { if (ainode[i] == nadj.inode) iVertex = i;
      if (ainode[i] ==	    inode) i1	   = i;
    }
    UINT iA = (iVertex + 1) % cVertices;
    UINT iB = (iVertex + cVertices - 1) % cVertices;
    UINT iV = (i1 == iA) ? iB : iA;
    UINT iE = aauVertices2Edge[iVertex][iV];

    nadj.auLabel[cnode] = mCellVector[itri].uLabel;
    ++cnode;
    nadj.ainode[cnode] = inode = ainode[iV];
    nadj.aiedge[cnode] = aiedge[iE];
    nadj.aielem[cnode] = aitri[iE];
    fBndry = fBndry || (auLabel[iE] != uIntLabel);
    itri = aitri[iE];

    if (nadj.ainode[cnode] == nadj.ainode[0])
    { nadj.fGotAll = true;
      break;
    }

    if (itri == -1) break;

    getElementInfo(itri, ainode, aiedge, aitri, auLabel);
    cVertices = (ainode[VERTEX_D] != -1) ? 4 : 3;
  }

  nadj.fBndry = fBndry;
  nadj.cnode  = (nadj.fGotAll) ? cnode : cnode + 1;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::calcAdjacencyInfo
//
// Purpose       : This function calcultes the node adjacency information.
//                 This includes things such as integration box areas,
//                 edge lengths, bisector lengths, partial areas, etc.
//                 It works for structured and non-structured meshes.
//                 (non-structured = triangle based)
//
// Special Notes : This is a modified version of a similar function from
//                 SGF's refine program.  Eventually, it needs to be
//                 refactored a bit to take better advantage of Xyce's data
//                 structures.  At the moment it relies a little too much
//                 on SGF mesh primatives, and C-style implementation.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::calcAdjacencyInfo ()
{
  int i, j, k;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "In PDE_2DMesh::calcAdjacencyInfo" << std::endl;
  }

  // if neccessary, do some allocations:
  if (!adjInfoAllocFlag)
  {
    afVisitedVec.resize(numNodes,0);
    aiBegin = new int [numRegLabels];
    aiEnd   = new int [numRegLabels];
  }

  // initialize the visit flag array
  for (i=0;i<numNodes;++i) afVisitedVec[i] = 0;

  // loop through each element
  for(i = 0; i < numCells; ++i)
  {
    // get the 3 or 4 nodes of this cell:
    int ainode[4];
    elementNodes(i, ainode);
    int *pinode = ainode;

    // loop over the 3 or 4 nodes of the current cell.
    for(j = 0; j < 4; ++j, ++pinode)
    {
      // only process this node if it has not been visited yet.
      if ( (afVisitedVec[*pinode]==0) && (*pinode != -1) )
      {
        afVisitedVec[*pinode] = 1;

        // store the adjacent nodes in the nadj structure
        NADJ nadj;
        initNodeAdjStructure(nadj, i, j, -1, false);

        // hit an edge w/o an adjacent triangle ... need to find the nodes
        // we may have missed
        if (!nadj.fGotAll)
        {
          NADJ nadjT;
          initNodeAdjStructure(nadjT, i, j, -1, true);
          int c = nadjT.cnode - 2;
          for(k = nadj.cnode - 1; k != -1; --k)
          {
            nadj.ainode [k + c] = nadj.ainode[k];
            nadj.aiedge [k + c] = nadj.aiedge[k];
            nadj.auLabel[k + c] = nadj.auLabel[k];
          }
          for(k = 0; k < c; ++k)
          {
            nadj.ainode [k] = nadjT.ainode[nadjT.cnode - k - 1];
            nadj.aiedge [k] = nadjT.aiedge[nadjT.cnode - k - 1];
            nadj.auLabel[k] = nadjT.auLabel[nadjT.cnode - k - 2];
          }
          nadj.cnode += c;
        }

        // compute the corners of the integration box
        int c = nadj.cnode;
        double x[32];
        double y[32];

        double xA   = mNodeVector[nadj.inode].x;
        double yA   = mNodeVector[nadj.inode].y;
        if (!nadj.fGotAll)
        {
          x[0] = (xA + mNodeVector[nadj.ainode[0]].x) / 2.0;
          y[0] = (yA + mNodeVector[nadj.ainode[0]].y) / 2.0;
        }
        else
        {
          computeIntPB(x[0],y[0],nadj.inode,nadj.ainode[0],nadj.ainode[c-1]);
        }

        for(k = 0; k < c-1; ++k)
        {
          computeIntPB
          (x[k+1],y[k+1],nadj.inode,nadj.ainode[k],nadj.ainode[k+1]);
        }

        if (!nadj.fGotAll)
        {
          x[c] = (xA + mNodeVector[nadj.ainode[c-1]].x) / 2.0;
          y[c] = (yA + mNodeVector[nadj.ainode[c-1]].y) / 2.0;
        }
        else
        { x[c] = x[0];
          y[c] = y[0];
        }

        // compute the node and edge information
        NODEINFO nodeinfo;
        EDGEINFO aedgeinfo[32];
        UINT cTriangle = (nadj.fGotAll) ? c : c-1;
        nodeinfo.cNeighbor = c;
        nodeinfo.cTriangle = cTriangle;
        nodeinfo.Area = 0.0;
        for(k = 0; k < c; ++k)
        {
          double xB = mNodeVector[nadj.ainode[k]].x;
          double yB = mNodeVector[nadj.ainode[k]].y;
          double elen = sqrt(sq(xA-xB)+sq(yA-yB));
          double ilen, area;
          if (cylGeom)
          {
            ilen = lengthAdjust(x[k],y[k],x[k+1],y[k+1]);
            area = areaAdjust(xA,yA,x[k],y[k],x[k+1],y[k+1]);
          }
          else
          {
            ilen = sqrt(sq(x[k]-x[k+1])+sq(y[k]-y[k+1]));
            area = ilen * elen * 0.25;
          }
          aedgeinfo[k].iedge = nadj.aiedge[k];
          aedgeinfo[k].inode = nadj.ainode[k];
          aedgeinfo[k].ielem = nadj.aielem[k];
          aedgeinfo[k].ilen  = ilen;
          aedgeinfo[k].elen  = elen;
          aedgeinfo[k].Area1 = area;
          aedgeinfo[k].Area2 = 0.0;
          nodeinfo.Area += area;
        }

        for(k = 0; k < cTriangle; ++k)
        {
          double xB = mNodeVector[nadj.ainode[k]].x;
          double yB = mNodeVector[nadj.ainode[k]].y;
          double xC = mNodeVector[nadj.ainode[k+1]].x;
          double yC = mNodeVector[nadj.ainode[k+1]].y;
          double xo  = x[k+1];
          double yo  = y[k+1];
          double xAB = (xA+xB)/2.0;
          double yAB = (yA+yB)/2.0;
          double xAC = (xA+xC)/2.0;
          double yAC = (yA+yC)/2.0;
          double area;
          if (cylGeom)
          {
            area = areaAdjust(xA,yA,xo,yo,xAB,yAB)+
                   areaAdjust(xA,yA,xo,yo,xAC,yAC);
          }
          else
          {
            area = 0.5*
              (sqrt(sq(xAB-xA)+sq(yAB-yA))*sqrt(sq(xAB-xo)+sq(yAB-yo))+
               sqrt(sq(xAC-xA)+sq(yAC-yA))*sqrt(sq(xAC-xo)+sq(yAC-yo)));
          }
          aedgeinfo[k].Area2 = area;
        }

        // fill in mNode class w/information from nodeinfo.
        mNodeVector[*pinode].area  = nodeinfo.Area;
        mNodeVector[*pinode].cnode = nodeinfo.cNeighbor;
        mNodeVector[*pinode].numCells = nodeinfo.cTriangle;

        // fill in the edgeInfo...
        mNodeVector[*pinode].edgeInfoVector.resize
        (mNodeVector[*pinode].cnode);
        for (k=0;k<mNodeVector[*pinode].cnode;++k)
        {
          mNodeVector[*pinode].edgeInfoVector[k] = aedgeinfo[k];

          mEdge & edgeTmp =
            mEdgeVector[mNodeVector[*pinode].edgeInfoVector[k].iedge];

          edgeTmp.elen  = mNodeVector[*pinode].edgeInfoVector[k].elen;
          edgeTmp.ilen  = mNodeVector[*pinode].edgeInfoVector[k].ilen;
          edgeTmp.Area1 = mNodeVector[*pinode].edgeInfoVector[k].Area1;
          edgeTmp.Area2 = mNodeVector[*pinode].edgeInfoVector[k].Area2;
        }
        // Is this node a boundary node?
	// may need to change this if statement later. It dependes upon the
	// node list being ordered a certain way.  Either have to use a boundary
	// stencil instead, or make certain to re-order all the nodes.

	// For a boundary node, we do extra stuff - a different
	// set of edgeInfo's have to be generated for each region.  So
	// if this node is on the boundary between SI and SIO2, then
	// there needs to be an edgeinfo for SI and SIO2, in addition
	// to the default (which covers both).
	// For now, multiple regions are not supported, so this stuff
	// is commented out.  FIX LATER.
#if 0
        int l;
        if (nadj.inode < numBndryNodes)
        {
	  // compute the start and end of the regions
          int *piBegin = aiBegin;
          int *piEnd  = aiEnd;
          for(k = 0; k < cRegion; ++k, ++piBegin, ++piEnd)
          {
            *piBegin = *piEnd = -1;
          }

          int c = nadj.cnode;
          if (!nadj.fGotAll) --c;
          int *puLabel = nadj.auLabel;
          for(k = 0; k < c; ++k, ++puLabel)
          {
            int uLabel = *puLabel;
            if (aiBegin[uLabel] == -1)
            {
              aiBegin[uLabel] = k;
              aiEnd[uLabel] = k + 1;
            }
            else if (aiEnd[uLabel] == k) ++aiEnd[uLabel];
          }
          int uLabel = nadj.auLabel[0];
          puLabel = nadj.auLabel + c - 1;
          for(k = c; *puLabel == uLabel; --k, --puLabel);
          aiBegin[uLabel] = k % c;

          // write the node and edge information to disk
          piBegin = aiBegin;
          piEnd = aiEnd;
          for(k = 0; k < cRegion; ++k, ++piBegin, ++piEnd)
          {
            int iBegin = *piBegin;
            int iEnd  = *piEnd;
            int cEdge = 0;
            if (iBegin != -1)
            {
              cEdge  = (iEnd > iBegin) ? iEnd - iBegin : (c + iEnd) - iBegin;
              ++cEdge;
            }
            iEnd %= nadj.cnode;
            int cTEdge = (cEdge) ? cEdge - 1 : 0;
            nodeinfo.cNeighbor = cEdge;
            nodeinfo.cTriangle = cTEdge;
            nodeinfo.Area = 0.0;
            for(l = 0; l < cEdge; ++l)
            {
              int m = (iBegin + l) % nadj.cnode;
              double xB = mNodeVector[nadj.ainode[m]].x;
              double yB = mNodeVector[nadj.ainode[m]].y;
              double xmdpt = (xA + xB) / 2.0;
              double ymdpt = (yA + yB) / 2.0;
              double  elen   = sqrt(sq(xA-xB)+sq(yA-yB));
              double  ilen, area, x1, y1, x2, y2;
              if      (m == iBegin)
              { x1 = xmdpt;  x2 = x[m+1];
                y1 = ymdpt;  y2 = y[m+1];
              }
              else if (m == iEnd)
              { x1 = x[m];   x2 = xmdpt;
                y1 = y[m];   y2 = ymdpt;
              }
              else
              { x1 = x[m];   x2 = x[m+1];
                y1 = y[m];   y2 = y[m+1];
              }
              if (cylGeom)
              { ilen = lengthAdjust(x1,y1,x2,y2);
                area = areaAdjust(xA,yA,x1,y1,x2,y2);
              }
              else
              { ilen = sqrt(sq(x1-x2)+sq(y1-y2));
                area = ilen * elen / 4.0;
              }
              aedgeinfo[l].iedge = nadj.aiedge[m];
              aedgeinfo[l].inode = nadj.ainode[m];
              aedgeinfo[l].ielem = nadj.aielem[m];
              aedgeinfo[l].ilen   = ilen;
              aedgeinfo[l].elen   = elen;
              aedgeinfo[l].Area1 = area;
              aedgeinfo[l].Area2 = 0.0;
              nodeinfo.Area += area;
            }

            for(l = 0; l < cTEdge; ++l)
            { int m = (iBegin + l) % nadj.cnode;
              int n = (m + 1) % nadj.cnode;
              //PNODE pnodeB = danode.GetPointer(nadj.ainode[m]);
              double xB = mNodeVector[nadj.ainode[m]].x;
              double yB = mNodeVector[nadj.ainode[m]].y;
              //PNODE pnodeC = danode.GetPointer(nadj.ainode[n]);
              double xC = mNodeVector[nadj.ainode[n]].x;
              double yC = mNodeVector[nadj.ainode[n]].y;
              double xo  = x[m+1];
              double yo  = y[m+1];
              double xAB = (xA+xB)/2.0;
              double yAB = (yA+yB)/2.0;
              double xAC = (xA+xC)/2.0;
              double yAC = (yA+yC)/2.0;
              double area;
              if (cylGeom)
              { area = areaAdjust(xA,yA,xo,yo,xAB,yAB)+
                       areaAdjust(xA,yA,xo,yo,xAC,yAC);
              }
              else
              { area = 0.5*
                (sqrt(sq(xAB-xA)+sq(yAB-yA))*sqrt(sq(xAB-xo)+sq(yAB-yo))+
                 sqrt(sq(xAC-xA)+sq(yAC-yA))*sqrt(sq(xAC-xo)+sq(yAC-yo)));
              }
              aedgeinfo[l].Area2 = area;
            }

            // write the node and edge information to disk
            int inode = numNodes + nadj.inode * cRegion + k;
#if 0
            write(nFile, &inode, sizeof(INODE));
            write(nFile, &nodeinfo, sizeof(NODEINFO) - sizeof(EDGEINFO *));
            write(nFile, aedgeinfo, cEdge * sizeof(EDGEINFO));
#else
            // fill in mNode class w/information from nodeinfo.
            // Check:  is the mNodeVector big enough.   FIX THIS!
            mNodeVector[inode].area  = nodeinfo.Area;
            mNodeVector[inode].cnode = nodeinfo.cNeighbor;
            mNodeVector[inode].numCells = nodeinfo.cTriangle;

            // fill in the edgeInfo...
            mNodeVector[inode].edgeInfoVector.resize
              (mNodeVector[inode].cnode);
            for (k=0;k<mNodeVector[inode].cnode;++k)
            {
              mNodeVector[inode].edgeInfoVector[k] = aedgeinfo[k];
              mEdge & edgeTmp =
                mEdgeVector[mNodeVector[*pinode].edgeInfoVector[k].iedge];

              edgeTmp.elen  = mNodeVector[*pinode].edgeInfoVector[k].elen;
              edgeTmp.ilen  = mNodeVector[*pinode].edgeInfoVector[k].ilen;
              edgeTmp.Area1 = mNodeVector[*pinode].edgeInfoVector[k].Area1;
              edgeTmp.Area2 = mNodeVector[*pinode].edgeInfoVector[k].Area2;
            }
#endif
          } // loop over regions.
        } // boundary node if statement.
#endif
      } // if !visited if statement.
    } // loop over nodes of cell.
  } // loop over cells.

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Done with PDE_2DMesh::calcAdjacencyInfo" << std::endl;
  }

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::labelNameExist
// Description   : This function returns a "true" if the label name specified
//                 in the function argument exists in the mesh data structures.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::labelNameExist (std::string & labelName)
{
  bool bsuccess = false;

  ExtendedString tmpName = labelName;
  tmpName.toUpper ();

  if ( mLabelMap.find(tmpName) != mLabelMap.end() ) bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::labelEdgeType
// Description   : This function returns a "true" if the label name specified
//                 corresponds to an edge label (rather than a region label)
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::labelEdgeType  (std::string & labelName)
{
  ExtendedString tmpName = labelName;
  tmpName.toUpper ();

  if ( mLabelMap.find(tmpName) != mLabelMap.end() )
  {
    if (mLabelMap[tmpName].uType == TYPE_EDGE) return true;
  }

  return false;

}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::getDopingVector
// Description   :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::getDopingVector (std::vector<double> & cvec_tmp)
{
  cvec_tmp.resize(dopingVector.size(), 0.0);
  copy(dopingVector.begin(),dopingVector.end(),cvec_tmp.begin());
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::getXVector
// Description   :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::getXVector (std::vector<double> & xvec_tmp)
{
  xvec_tmp.resize(xVector.size(), 0.0);
  copy(xVector.begin(),xVector.end(),xvec_tmp.begin());
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::getYVector
// Description   :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::getYVector (std::vector<double> & yvec_tmp)
{
  yvec_tmp.resize(yVector.size(), 0.0);
  copy(yVector.begin(),yVector.end(),yvec_tmp.begin());
  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::getLabel
// Description   : returns a pointer to the label with the specified name.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/27/02
//-----------------------------------------------------------------------------
mLabel * PDE_2DMesh::getLabel (std::string & labelName)
{

  ExtendedString tmpName = labelName;
  tmpName.toUpper ();

  int index = 0;
  if ( mLabelMap.find(tmpName) != mLabelMap.end() )
  {
    index = mLabelMap[tmpName].iIndex;
  }

  return &(mLabelVector[index]);
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::scaleMesh
//
// Description   : This function scales mesh quantities such as area, ilen and
//                 elen based on the passed-in scalar.  The function accounts
//                 for the geomtry (cylindrical or cartesian)
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/29/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::scaleMesh (double xScale)
{
  if (meshScaledFlag == true)
  {
    meshScaledFlag = false;
  }
  else
  {
    meshScaledFlag = true;
  }
  // save the scale constant.
  x0 = xScale;

  double X0_1 = xScale;
  double X0_2 = xScale * X0_1;
  double X0_3 = xScale * X0_2;

  double scaleILEN, scaleArea;
  double scaleELEN = X0_1;

  if (cylGeom) { scaleILEN = X0_2; scaleArea = X0_3; }
  else         { scaleILEN = X0_1; scaleArea = X0_2; }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << " In PDE_2DMesh::scaleMesh"<<std::endl;
    Xyce::dout() << "  scaleELEN = " << scaleELEN <<std::endl;
    Xyce::dout() << "  scaleILEN = " << scaleILEN <<std::endl;
    Xyce::dout() << "  scaleArea = " << scaleArea <<std::endl;
    Xyce::dout() << section_divider << std::endl;
  }

  // Make everything reciprocals, so can use * instead of /.
  scaleELEN = 1.0/scaleELEN;
  scaleILEN = 1.0/scaleILEN;
  scaleArea = 1.0/scaleArea;

  int i;
  for (i=0;i<numNodes;++i)
  {
    mNodeVector[i].area *= scaleArea;

    std::vector<EDGEINFO>::iterator firstEI = mNodeVector[i].edgeInfoVector.begin ();
    std::vector<EDGEINFO>::iterator lastEI  = mNodeVector[i].edgeInfoVector.end ();
    std::vector<EDGEINFO>::iterator iterEI;
    for (iterEI=firstEI;iterEI!=lastEI;++iterEI)
    {
      iterEI->ilen  *= scaleILEN;
      iterEI->elen  *= scaleELEN;
      iterEI->Area1 *= scaleArea;
      iterEI->Area2 *= scaleArea;
    }
  }

  for (i=0;i<numEdges;++i)
  {
    mEdgeVector[i].ilen  *= scaleILEN;
    mEdgeVector[i].elen  *= scaleELEN;
    mEdgeVector[i].Area1 *= scaleArea;
    mEdgeVector[i].Area2 *= scaleArea;
  }

  for (i=0;i<numLabels;++i)
  {
    mLabelVector[i].vol *= scaleArea;
    mLabelVector[i].surfArea *= scaleILEN;
  }

  std::map<std::string, mLabel>::iterator firstL = mLabelMap.begin ();
  std::map<std::string, mLabel>::iterator lastL  = mLabelMap.end   ();
  std::map<std::string, mLabel>::iterator iterL;

  for (iterL=firstL; iterL!=lastL; ++iterL)
  {
    iterL->second.vol *= scaleArea;
    iterL->second.surfArea *= scaleILEN;
  }

  xMax *= scaleELEN;
  xMin *= scaleELEN;
  yMax *= scaleELEN;
  yMin *= scaleELEN;

  for (i=0;i<numNodes;++i)
  {
    xVector[i]       *= scaleELEN;
    yVector[i]       *= scaleELEN;
    mNodeVector[i].x *= scaleELEN;
    mNodeVector[i].y *= scaleELEN;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function    : PDE_2DMesh::interp
// Author      : Eric Keiter
// Scope       : public
// Description : This file returns an interpolated value for a function,
//               F, given the spatial location(r,z).
//
//-----------------------------------------------------------------------------
double PDE_2DMesh::interp(double *F, double r, double z)
{
  int iNodeA,iNodeB,iNodeC,iNodeD;

  mInterpAreaHelp intAHelp;

  double func;

  int iCell;

  int istatus = 0;
  int inode;
  findCell(r, z, istatus, inode, iCell, iRecentCellLookup);
  iRecentCellLookup = iCell;

  if (istatus ==  2) return F[inode];
  if (istatus == -1) return 0.0;

  iNodeA = mCellVector[iCell].inodeA;
  iNodeB = mCellVector[iCell].inodeB;
  iNodeC = mCellVector[iCell].inodeC;
  iNodeD = mCellVector[iCell].inodeD;

  if(iNodeD == -1)
  {
    intAHelp.x0 = xVector[iNodeA];
    intAHelp.x1 = xVector[iNodeB];
    intAHelp.x2 = xVector[iNodeC];
    intAHelp.y0 = yVector[iNodeA];
    intAHelp.y1 = yVector[iNodeB];
    intAHelp.y2 = yVector[iNodeC];
    intAHelp.f0 = F[iNodeA];
    intAHelp.f1 = F[iNodeB];
    intAHelp.f2 = F[iNodeC];
  }
  else    // need to figure out which half of rectangle.
  {       // Either:  ABC or ADC.
    intAHelp.x0 = xVector[iNodeA];
    intAHelp.x2 = xVector[iNodeC];
    intAHelp.y0 = yVector[iNodeA];
    intAHelp.y2 = yVector[iNodeC];
    intAHelp.f0 = F[iNodeA];
    intAHelp.f2 = F[iNodeC];

    double r1 = xVector[iNodeB];
    double r2 = xVector[iNodeD];
    double z1 = yVector[iNodeB];
    double z2 = yVector[iNodeD];
    double f1 = F[iNodeB];
    double f2 = F[iNodeD];

    if(  pow((r1-r),2.0)+pow((z1-z),2.0) <
         pow((r2-r),2.0)+pow((z2-z),2.0)
      )
    {
      intAHelp.x1 = r1;
      intAHelp.y1 = z1;
      intAHelp.f1 = f1;
    }
    else
    {
      intAHelp.x1 = r2;
      intAHelp.y1 = z2;
      intAHelp.f1 = f2;
    }
  }

  // call the FindCoef function
  func = intAHelp.interpReg(r,z);

  return(func);
}

//-----------------------------------------------------------------------------
// Function:     PDE_2DMesh::interpVector
// Author:       Eric Keiter
// Scope:        public
// Description:  This file returns an interpolated value for a function,
//               F, given the spatial location(r,z).
//
//               Note that the passed double precision array, F, is expected
//               to be a vector array, which corresponds to mesh edges.  The
//               number of elements in it should equal the size of mEdgeVector.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/24/02
//-----------------------------------------------------------------------------
bool PDE_2DMesh::interpVector (
  double *F,
  double r,
  double z,
  double & xvec,
  double & yvec)
{

  mInterpAreaHelp intAHelp;
  int iEdgeAB, iEdgeBC, iEdgeCD, iEdgeDA;

  double ABx, ABy, ABangle, ABFx, ABFy;
  double BCx, BCy, BCangle, BCFx, BCFy;
  double CDx, CDy, CDangle, CDFx, CDFy;
  double DAx, DAy, DAangle, DAFx, DAFy;
  double xA, xB, yA, yB;
  int inodeA, inodeB;

  double alpha = 0.0;
  int inodeC, inodeD;

  // Find the cell.
  int iCell;
  int istatus = 0;
  int inode;
  findCell(r, z, istatus, inode, iCell,iRecentCellLookup);
  iRecentCellLookup = iCell;

  xvec = 0.0;
  yvec = 0.0;

  if (istatus==-1)
  {
    xvec = 0.0;
    yvec = 0.0;
    return true;
  }

  // Is this cell a rectangle or a triangle?
  bool triFlag = (mCellVector[iCell].inodeD == -1);

  int iA = mCellVector[iCell].inodeA;
  int iB = mCellVector[iCell].inodeB;
  int iC = mCellVector[iCell].inodeC;
  int iD = mCellVector[iCell].inodeD;

  // Look at each edge of the cell.  Find the angle of each edge
  // with the x-axis.
  iEdgeAB = mCellVector[iCell].iedgeAB;
  iEdgeBC = mCellVector[iCell].iedgeBC;
  iEdgeCD = mCellVector[iCell].iedgeCD;
  if (triFlag) iEdgeDA = -1;
  else         iEdgeDA = mCellVector[iCell].iedgeDA;

  //get the midpoints, angles of the edges:
  // AB
  inodeA = mEdgeVector[iEdgeAB].inodeA; xA=mNodeVector[inodeA].x; yA=mNodeVector[inodeA].y;
  inodeB = mEdgeVector[iEdgeAB].inodeB; xB=mNodeVector[inodeB].x; yB=mNodeVector[inodeB].y;
  ABx = 0.5*(xA+xB);
  ABy = 0.5*(yA+yB);
  ABangle = compAngle(xB,yB, xA,yA, (xA+0.1), yA);
  mEdgeVector[iEdgeAB].midpoint_x = ABx;
  mEdgeVector[iEdgeAB].midpoint_y = ABy;
  mEdgeVector[iEdgeAB].angle      = ABangle;
  ABFx = F[iEdgeAB] * cos(ABangle);
  ABFy = F[iEdgeAB] * sin(ABangle);

  // BC
  inodeA = mEdgeVector[iEdgeBC].inodeA; xA=mNodeVector[inodeA].x; yA=mNodeVector[inodeA].y;
  inodeB = mEdgeVector[iEdgeBC].inodeB; xB=mNodeVector[inodeB].x; yB=mNodeVector[inodeB].y;
  BCx = 0.5*(xA+xB);
  BCy = 0.5*(yA+yB);
  BCangle = compAngle(xB,yB, xA,yA, (xA+0.1), yA);
  mEdgeVector[iEdgeBC].midpoint_x = BCx;
  mEdgeVector[iEdgeBC].midpoint_y = BCy;
  mEdgeVector[iEdgeBC].angle      = BCangle;
  BCFx = F[iEdgeBC] * cos(BCangle);
  BCFy = F[iEdgeBC] * sin(BCangle);

  // CD
  inodeA = mEdgeVector[iEdgeCD].inodeA; xA=mNodeVector[inodeA].x; yA=mNodeVector[inodeA].y;
  inodeB = mEdgeVector[iEdgeCD].inodeB; xB=mNodeVector[inodeB].x; yB=mNodeVector[inodeB].y;
  CDx = 0.5*(xA+xB);
  CDy = 0.5*(yA+yB);
  CDangle = compAngle(xB,yB, xA,yA, (xA+0.1), yA);
  mEdgeVector[iEdgeCD].midpoint_x = CDx;
  mEdgeVector[iEdgeCD].midpoint_y = CDy;
  mEdgeVector[iEdgeCD].angle      = CDangle;
  CDFx = F[iEdgeCD] * cos(CDangle);
  CDFy = F[iEdgeCD] * sin(CDangle);

  // DA
  if (!triFlag)
  {
    inodeA = mEdgeVector[iEdgeDA].inodeA; xA=mNodeVector[inodeA].x; yA=mNodeVector[inodeA].y;
    inodeB = mEdgeVector[iEdgeDA].inodeB; xB=mNodeVector[inodeB].x; yB=mNodeVector[inodeB].y;
    DAx = 0.5*(xA+xB);
    DAy = 0.5*(yA+yB);
    DAangle = compAngle(xB,yB, xA,yA, (xA+0.1), yA);
    mEdgeVector[iEdgeDA].midpoint_x = DAx;
    mEdgeVector[iEdgeDA].midpoint_y = DAy;
    mEdgeVector[iEdgeDA].angle      = DAangle;
    DAFx = F[iEdgeDA] * cos(DAangle);
    DAFy = F[iEdgeDA] * sin(DAangle);
  }
  else
  {
    DAx = 0.0;
    DAy = 0.0;
    DAangle = 0.0;
    DAFx = 0.0;
    DAFy = 0.0;
  }

  bool madeIt = false;

  if (triFlag)
  {
    // If this is a triangular cell, do a "interp" interpolation to
    // get Fx and Fy at (r,z).
    intAHelp.x0 = ABx;
    intAHelp.x1 = BCx;
    intAHelp.x2 = CDx;
    intAHelp.y0 = ABy;
    intAHelp.y1 = BCy;
    intAHelp.y2 = CDy;
    intAHelp.f0 = ABFx;
    intAHelp.f1 = BCFx;
    intAHelp.f2 = CDFx;

    xvec = intAHelp.interpReg(r,z);

    intAHelp.f0 = ABFy;
    intAHelp.f1 = BCFy;
    intAHelp.f2 = CDFy;

    yvec = intAHelp.interpReg(r,z);

    madeIt = true;
  }
  else  // rectangle:  Note this assumes that if not a triangle,
        //             the only other option is rectangle.
  {
    double alpha;
    int iA = mCellVector[iCell].inodeA;
    int iB = mCellVector[iCell].inodeB;
    int iC = mCellVector[iCell].inodeC;
    int iD = mCellVector[iCell].inodeD;

    if (mNodeVector[iA].y == mNodeVector[iB].y)  // if AB and CD are parallel to x-axis.
    {
      alpha = ((ABy <  CDy)?1.0:0.0) * (z-ABy)/(CDy-ABy) +
              ((ABy >= CDy)?1.0:0.0) * (ABy-z)/(ABy-CDy);
      xvec = (1.0-alpha)*ABFx + alpha*CDFx;

      alpha = ((BCx <  DAx)?1.0:0.0) * (r-BCx)/(DAx-BCx) +
              ((BCx >= DAx)?1.0:0.0) * (BCx-r)/(BCx-DAx);
      yvec = (1.0-alpha)*BCFy + alpha*DAFy;

      madeIt = true;
    }
    else if (mNodeVector[iA].x == mNodeVector[iB].x) // if AB and CD are parallel to y-axis.
    {
      alpha = ((BCy <  DAy)?1.0:0.0) * (z-BCy)/(DAy-BCy) +
              ((BCy >= DAy)?1.0:0.0) * (BCy-z)/(BCy-DAy);
      xvec = (1.0-alpha)*BCFx + alpha*DAFx;

      alpha = ((ABx <  CDx)?1.0:0.0) * (r-ABx)/(CDx-ABx) +
              ((ABx >= CDx)?1.0:0.0) * (ABx-r)/(ABx-CDx);
      yvec = (1.0-alpha)*ABFy + alpha*CDFy;

      madeIt = true;
    }
    else
    {
      madeIt = false;
    }
  }

  if (DEBUG_DEVICE)
  {
    // debug output:
    double rtest = ((2.1e-3)/209.0) * 22.0;
    double rtol  = fabs(rtest/1000.0);
    double ztest = -9.0e-5;
    double ztol  = fabs(ztest/1000.0);

    if (!madeIt ||
        (xvec != 0.0  && !(xvec > 0.0)  && !(xvec < 0.0)) ||
        (yvec != 0.0  && !(yvec > 0.0)  && !(yvec < 0.0))
        || (r >= (rtest-rtol) && r <= (rtest+rtol) &&
            z >= (ztest-ztol) && z <= (ztest+ztol))
        )
    {
      Xyce::dout() << Xyce::section_divider << std::endl;
      Xyce::dout() << "Vector Interpolation failed!" << std::endl;
      Xyce::dout() << std::endl;
      Xyce::dout() << "  iCell  = " << iCell << std::endl;
      Xyce::dout() << "  alpha  = " << alpha << std::endl;
      Xyce::dout() << "  number of cells on mesh: " << numCells << std::endl;
      Xyce::dout() << "  r      = " << r << std::endl;
      Xyce::dout() << "  z      = " << z << std::endl;
      Xyce::dout() << std::endl;
      Xyce::dout() << "  xvec   = " << xvec << std::endl;
      Xyce::dout() << "  yvec   = " << yvec << std::endl;

      Xyce::dout() << "  inodeA = " << iA << std::endl;
      Xyce::dout() << "  inodeB = " << iB << std::endl;
      Xyce::dout() << "  inodeC = " << iC << std::endl;
      Xyce::dout() << "  inodeD = " << iD << std::endl;

      Xyce::dout() << "  nodeA (x,y) = ("<<mNodeVector[iA].x<<", "<<mNodeVector[iA].y<<")"<<std::endl;
      Xyce::dout() << "  nodeB (x,y) = ("<<mNodeVector[iB].x<<", "<<mNodeVector[iB].y<<")"<<std::endl;
      Xyce::dout() << "  nodeC (x,y) = ("<<mNodeVector[iC].x<<", "<<mNodeVector[iC].y<<")"<<std::endl;
      if (iD != -1)
        Xyce::dout() << "  nodeD (x,y) = ("<<mNodeVector[iD].x<<", "<<mNodeVector[iD].y<<")"<<std::endl;

      double pi = M_PI;
      Xyce::dout() << "  ABangle = " << ABangle << " = " <<(ABangle/pi)<< " * PI" << std::endl;
      Xyce::dout() << "  BCangle = " << BCangle << " = " <<(BCangle/pi)<< " * PI" << std::endl;
      Xyce::dout() << "  CDangle = " << CDangle << " = " <<(CDangle/pi)<< " * PI" << std::endl;
      if (iD != -1)
        Xyce::dout() << "  DAangle = " << DAangle << " = " <<(DAangle/pi)<< " * PI" << std::endl;

      Xyce::dout() << "  ABFx = " << ABFx << std::endl;
      Xyce::dout() << "  ABFy = " << ABFy << std::endl;
      Xyce::dout() << "  BCFx = " << BCFx << std::endl;
      Xyce::dout() << "  BCFy = " << BCFy << std::endl;
      Xyce::dout() << "  CDFx = " << CDFx << std::endl;
      Xyce::dout() << "  CDFy = " << CDFy << std::endl;
      Xyce::dout() << "  DAFx = " << DAFx << std::endl;
      Xyce::dout() << "  DAFy = " << DAFy << std::endl;

      Xyce::dout() << "  F[iEdgeAB] = " << F[iEdgeAB] <<std::endl;
      Xyce::dout() << "  F[iEdgeBC] = " << F[iEdgeBC] <<std::endl;
      Xyce::dout() << "  F[iEdgeCD] = " << F[iEdgeCD] <<std::endl;
      if (iEdgeDA != -1)
        Xyce::dout() << "  F[iEdgeDA] = " << F[iEdgeDA] <<std::endl;
      Xyce::dout() << Xyce::section_divider << std::endl;
      Report::DevelFatal() << "Vector Interpolation failed";
    }

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DMesh::findCell
// Author        : Eric Keiter
// Scope         : public
//
// Description   : This function performs a search over the mesh to find the
//                 mesh cell which contains the point (r,z).  The test at each
//                 cell requires that (r,z) be compared to the 3 or 4
//                 lines(edges) which comprise the edges of the cell.
//                 If it is found to be on the correct side of all 3 or
//                 4 edges, it is declared to be the correct cell and
//                 the function exits.
//
//                 The original version of this function looped over the
//                 array of mesh cells, and stopped when it had found the
//                 correct one.  This was really slow.  This version
//                 marches accross the mesh in a linked-list fashion,
//                 going from cell to neighbor cell.  It chooses the
//                 next cell in the search based on which of the
//                 neighbors is closest to (r,z).
//
//                 This function is usually called as part of an
//                 interpolation, which is usually performed for the
//                 purpose of plotting results.
//
//                 As such, the result of the previous "findCell" search
//                 is often a very good starting cell for the current
//                 search.  The argument "iStartCell" should be the
//                 callers best guess for which cell is likely to be the
//                 correct cell.  Usually, if one uses the class
//                 variable "iRecentCellLookup" as the iStartCell argument,
//                 the search will be very quick - at most it will look at
//                 2 or 3 cells.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/27/02
//-----------------------------------------------------------------------------
void PDE_2DMesh::findCell(
    double r,
    double z,
    int & istatus,
    int & inode,
    int & iCell,
    int iStartCell)
{
  int iEdgeAB,iEdgeBC,iEdgeCD,iEdgeDA;
  int iNodeA,iNodeB,iNodeC,iNodeD;
  bool doneFlag = false;

  mInterpEdgeHelp intEHelp;

  int ixhi,ixlo,iyhi,iylo;

  istatus = 0;

  iCell = iStartCell;

  // note:  change this to a STL function call later, to speed it up.
  for (int i1=0;i1<numCells; ++i1) visitCellFlagVec[i1] = 0;

  // find the proper triangular/rectangular region
  while (!doneFlag)
  {
    iEdgeAB = mCellVector[iCell].iedgeAB;
    iEdgeBC = mCellVector[iCell].iedgeBC;
    iEdgeCD = mCellVector[iCell].iedgeCD;
    iEdgeDA = mCellVector[iCell].iedgeDA;

    iNodeA  = mCellVector[iCell].inodeA;
    iNodeB  = mCellVector[iCell].inodeB;
    iNodeC  = mCellVector[iCell].inodeC;
    iNodeD  = mCellVector[iCell].inodeD;

    // if (r,z) is exactly on one of the nodes, then the search is done.
    if(r==xVector[iNodeA] && z==yVector[iNodeA]){inode=iNodeA; istatus=2;}
    if(r==xVector[iNodeB] && z==yVector[iNodeB]){inode=iNodeB; istatus=2;}
    if(r==xVector[iNodeC] && z==yVector[iNodeC]){inode=iNodeC; istatus=2;}
    if(iNodeD != -1)
    if(r==xVector[iNodeD] && z==yVector[iNodeD]){inode=iNodeD; istatus=2;}

    if (istatus == 2) return;

    ixhi = 0;
    ixlo = 0;
    iyhi = 0;
    iylo = 0;

    intEHelp.xA = xVector[mEdgeVector[iEdgeAB].inodeA];
    intEHelp.xB = xVector[mEdgeVector[iEdgeAB].inodeB];
    intEHelp.yA = yVector[mEdgeVector[iEdgeAB].inodeA];
    intEHelp.yB = yVector[mEdgeVector[iEdgeAB].inodeB];
    intEHelp.setupEdge(r,z);

    if(intEHelp.x_hiFlag) ixhi += 1;
    if(intEHelp.x_loFlag) ixlo += 1;
    if(intEHelp.y_hiFlag) iyhi += 1;
    if(intEHelp.y_loFlag) iylo += 1;

    intEHelp.xA = xVector[mEdgeVector[iEdgeBC].inodeA];
    intEHelp.xB = xVector[mEdgeVector[iEdgeBC].inodeB];
    intEHelp.yA = yVector[mEdgeVector[iEdgeBC].inodeA];
    intEHelp.yB = yVector[mEdgeVector[iEdgeBC].inodeB];
    intEHelp.setupEdge(r,z);

    if(intEHelp.x_hiFlag) ixhi += 1;
    if(intEHelp.x_loFlag) ixlo += 1;
    if(intEHelp.y_hiFlag) iyhi += 1;
    if(intEHelp.y_loFlag) iylo += 1;

    intEHelp.xA = xVector[mEdgeVector[iEdgeCD].inodeA];
    intEHelp.xB = xVector[mEdgeVector[iEdgeCD].inodeB];
    intEHelp.yA = yVector[mEdgeVector[iEdgeCD].inodeA];
    intEHelp.yB = yVector[mEdgeVector[iEdgeCD].inodeB];
    intEHelp.setupEdge(r,z);

    if(intEHelp.x_hiFlag) ixhi += 1;
    if(intEHelp.x_loFlag) ixlo += 1;
    if(intEHelp.y_hiFlag) iyhi += 1;
    if(intEHelp.y_loFlag) iylo += 1;

    if(iNodeD != -1)
    {
      intEHelp.xA = xVector[mEdgeVector[iEdgeDA].inodeA];
      intEHelp.xB = xVector[mEdgeVector[iEdgeDA].inodeB];
      intEHelp.yA = yVector[mEdgeVector[iEdgeDA].inodeA];
      intEHelp.yB = yVector[mEdgeVector[iEdgeDA].inodeB];
      intEHelp.setupEdge(r,z);

      if(intEHelp.x_hiFlag) ixhi += 1;
      if(intEHelp.x_loFlag) ixlo += 1;
      if(intEHelp.y_hiFlag) iyhi += 1;
      if(intEHelp.y_loFlag) iylo += 1;
    }

    // Have we found the cell?
    if((ixhi >= 1) && (ixlo >= 1) && (iyhi >= 1) && (iylo >= 1))
    {
      // yes!
      doneFlag = true;
      istatus = 0;
      break;
    }
    else
    {
      // no! move to a neighbor cell.
      doneFlag = false;
      istatus = -1;

      // of the 3 or 4 neighbor cells, which one is closest to (r,z)
      double minimum = +1.0e+99;
      int iC;

      // first get distance for current cell:
      double curr_min = findMinDist(iCell,r,z);

      // neighbor AB
      int AB_iC    = mCellVector[iCell].icellAB;
      double ABmin;
      if (AB_iC != -1)
      {
        if (!(visitCellFlagVec[AB_iC]))
        {
          ABmin = findMinDist(AB_iC,r,z);
          if (ABmin <= minimum)  { minimum = ABmin; iC = AB_iC;  }
        }
      }

      // neighbor BC
      int BC_iC    = mCellVector[iCell].icellBC;
      double BCmin;
      if (BC_iC != -1)
      {
        if (!(visitCellFlagVec[BC_iC]))
        {
          BCmin = findMinDist(BC_iC,r,z);
          if (BCmin <= minimum)  { minimum = BCmin; iC = BC_iC;  }
        }
      }

      // neighbor CD
      int CD_iC    = mCellVector[iCell].icellCD;
      double CDmin;
      if (CD_iC != -1)
      {
        if (!(visitCellFlagVec[CD_iC]))
        {
          CDmin = findMinDist(CD_iC,r,z);
          if (CDmin <= minimum)  { minimum = CDmin; iC = CD_iC;  }
        }
      }

      // neighbor DA
      int DA_iC    = mCellVector[iCell].icellDA;
      double DAmin = 0.0;
      if (DA_iC != -1)
      {
        if (!(visitCellFlagVec[DA_iC]))
        {
          DAmin = findMinDist(DA_iC,r,z);
          if (DAmin < minimum)  { minimum = DAmin; iC = DA_iC;  }
        }
      }

      // now reset iCell.
      if (iC != iCell)
      {
        visitCellFlagVec[iCell] = 1;
        iCell = iC;
        istatus=0;
      }
      else // we are stuck, and can't go anywhere.
           // This cell doesn't exist, so exit.
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << std::endl;
          Xyce::dout() << "iCell = " << iCell << std::endl;
          Xyce::dout() << "r = " << r << "   z = " << z << std::endl;
          Xyce::dout() << "curr_min = " << curr_min << std::endl;
          Xyce::dout() << std::endl;

          if (AB_iC  != -1)
          {
            Xyce::dout() << "AB_iC = " << AB_iC;
            Xyce::dout() << "  ABmin = " << ABmin;
            Xyce::dout() << "  visit = " << visitCellFlagVec[AB_iC] << std::endl;
          }

          if (BC_iC  != -1)
          {
            Xyce::dout() << "BC_iC = " << BC_iC;
            Xyce::dout() << "  BCmin = " << BCmin;
            Xyce::dout() << "  visit = " << visitCellFlagVec[BC_iC] << std::endl;
          }

          if (CD_iC  != -1)
          {
            Xyce::dout() << "CD_iC = " << CD_iC;
            Xyce::dout() << "  CDmin = " << CDmin;
            Xyce::dout() << "  visit = " << visitCellFlagVec[CD_iC] << std::endl;
          }

          if (DA_iC  != -1)
          {
            Xyce::dout() << "DA_iC = " << DA_iC;
            Xyce::dout() << "  DAmin = " << DAmin;
            Xyce::dout() << "  visit = " << visitCellFlagVec[DA_iC] << std::endl;
          }
        }
        doneFlag = true;
        istatus = -1;
      }
    }
  } // end of while loop

  return;
}

//-----------------------------------------------------------------------------
// Function:     PDE_2DMesh::findMinDist
// Author:       Eric Keiter
// Scope:        public
// Description:  This function fines the minimum distance between an (r,z)
//               location and the nodes of a specified cell.  Whichever node
//               is closest, the distance between that node and (r,z)
//               is returned.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/24/02
//-----------------------------------------------------------------------------
double PDE_2DMesh::findMinDist(int iCell, double r, double z)
{
  double minDist = +1.0e99;
  double x1, y1;
  double rdist,zdist,dist;

  int inodeA = mCellVector[iCell].inodeA;
  if (inodeA != -1)
  {
    x1 = xVector[inodeA];
    y1 = yVector[inodeA];
    rdist = r-x1; zdist = z-y1;
    dist = sqrt (rdist*rdist + zdist*zdist);
    if (dist  < minDist) { minDist = dist; }
  }

  int inodeB = mCellVector[iCell].inodeB;
  if (inodeB != -1)
  {
    x1 = xVector[inodeB];
    y1 = yVector[inodeB];
    rdist = r-x1; zdist = z-y1;
    dist = sqrt (rdist*rdist + zdist*zdist);
    if (dist  < minDist) { minDist = dist; }
  }

  int inodeC = mCellVector[iCell].inodeC;
  if (inodeC != -1)
  {
    x1 = xVector[inodeC];
    y1 = yVector[inodeC];
    rdist = r-x1; zdist = z-y1;
    dist = sqrt (rdist*rdist + zdist*zdist);
    if (dist  < minDist) { minDist = dist; }
  }

  int inodeD = mCellVector[iCell].inodeD;
  if (inodeD != -1)
  {
    x1 = xVector[inodeD];
    y1 = yVector[inodeD];
    rdist = r-x1; zdist = z-y1;
    dist = sqrt (rdist*rdist + zdist*zdist);
    if (dist  < minDist) { minDist = dist; }
  }

  return minDist;
}

//-----------------------------------------------------------------------------
// Function:     PDE_2DMesh::compAngle
// Author:       Eric Keiter
// Scope:        public
// Description:  This function computes the angle between three points given
//               by coordinates (x1,y1), (x2,y2) and (x3,y3).  The angle
//               is measured from vector 12 to 32.
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/24/02
//-----------------------------------------------------------------------------
double PDE_2DMesh::compAngle(
     double x1, double y1,             // coordinate (x1,y1)
     double x2, double y2,             // coordinate (x2,y2)
     double x3, double y3)             // coordinate (x3,y3)
{
  double xl,yl;                        // coordinates of vector 12
  double xr,yr;                        // coordinates of vector 32
  double r;                            // temporary number
  double angle;                        // angle between vectors 12 & 32

  xl = x1 - x2;  yl = y1 - y2;
  xr = x3 - x2;  yr = y3 - y2;

  r = (xl*xr+yl*yr)/(sqrt(sq(xl)+sq(yl))*sqrt(sq(xr)+sq(yr)));
  if      (r >  1.0) r =  1.0;
  else if (r < -1.0) r = -1.0;
  angle = acos(r);
  if (xl*yr-xr*yl > 0) angle = 2*M_PI-angle;
  return(angle);
}

//-----------------------------------------------------------------------------
// Functions associated with mesh primitive classes
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : mNode::mNode
// Description   : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
mNode::mNode () :
  x(0.0),
  y(0.0),
  area(0.0),
  cnode(0),
  inode(-1),
  numCells(0),
  edgeStatus(EDGESTATUS_INTERIOR),
  fBndry(false),
  fGotAll(false)
{

}

//-----------------------------------------------------------------------------
// Function      : mEdge::mEdge
// Description   : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
mEdge::mEdge () :
  uLabel(-1),
  inodeA(-1),
  inodeB(-1),
  edgeStatus(EDGESTATUS_INTERIOR),
  ilen(0.0),
  elen(0.0),
  Area1(0.0),
  Area2(0.0),
  angle(0.0),
  midpoint_x(0.0),
  midpoint_y(0.0),
  iedge(-1),
  ielem(-1)
{

}

//-----------------------------------------------------------------------------
// Function      : mCell::mCell
// Description   : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
mCell::mCell () :
  uLabel (-1),
  iedgeAB (-1), iedgeBC (-1), iedgeCD (-1), iedgeDA (-1),
  icellAB (-1), icellBC (-1), icellCD (-1), icellDA (-1)
{
  mNodeVector.resize(4,-1);
}

//-----------------------------------------------------------------------------
// Function      : mLabel::mLabel
// Description   : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
mLabel::mLabel () :
    name("no name"),
    iIndex(-1),
    uType(-1),
    cNode(0),
    vol(0.0),
    surfArea(0.0)
{

}


//-----------------------------------------------------------------------------
// Function      : mInterpAreaHelp::mInterpAreaHelp
// Description   : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
mInterpAreaHelp::mInterpAreaHelp () :
    x0  (0.0),
    y0  (0.0),
    x1  (0.0),
    y1  (0.0),
    x2  (0.0),
    y2  (0.0),
    v0  (0.0),
    v1  (0.0),
    v2  (0.0),
    f0  (0.0),
    f1  (0.0),
    f2  (0.0),
    vlim(0.0),
    aa  (0.0),
    bb  (0.0),
    cc  (0.0),
    errorFlag(0),
    iend(0)
{

}

//-----------------------------------------------------------------------------
// Function      : mInterpAreaHelp::interpReg
// Description   :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
double mInterpAreaHelp::interpReg (double r, double z)
{
  findCoef();
  double f = aa*r + bb*z + cc;
  return(f);
}

//-----------------------------------------------------------------------------
// Function      : mInterpAreaHelp::findCoef
// Description   : This function returns the coefficients of the linear
//                 equation f = a*x + b*y + c.  The class needs to have
//                 been passed the values of f,x and y for three different
//                 points.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
bool mInterpAreaHelp::findCoef ()
{
  if( (y1-y0)*(x2-x1)- (y2-y1)*(x1-x0) != 0.0)
  {
    bb = ((f1-f0)*(x2-x1) - (f2-f1)*(x1-x0))/
         ((y1-y0)*(x2-x1) - (y2-y1)*(x1-x0));
  }
  else
  {
    bb = 0.0;
  }

  if(x1!=x0)
  {
    aa = (f1-f0)/(x1-x0) - bb*(y1-y0)/(x1-x0);
  }
  else
  {
    if(x2!=x1)
    {
      aa = (f2-f1)/(x2-x1) - bb*(y2-y1)/(x2-x1);
    }
    else
    {
      aa = 0.0;
    }
  }

  cc = f0 - aa*x0 - bb*y0;

  return(true);
}

//-----------------------------------------------------------------------------
// Function      : mInterpEdgeHelp::mInterpEdgeHelp
// Description   : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
mInterpEdgeHelp::mInterpEdgeHelp () :
    xA(0.0),
    yA(0.0),
    xB(0.0),
    yB(0.0),
    AA(0.0),
    BB(0.0),
    iflagx(false),
    iflagy(false),
    x_hiFlag(false),
    x_loFlag(false),
    y_hiFlag(false),
    y_loFlag(false)
{

}

//-----------------------------------------------------------------------------
// Function      : mInterpEdgeHelp::setupEdge
// Description   :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
bool mInterpEdgeHelp::setupEdge (double r, double z)
{
  double xtmp, ytmp;

  AA = 0.0;
  BB = 0.0;

  // is this a verticle line?
  if(xB == xA)
  {
    iflagx   = false;
    y_hiFlag = false;
    y_loFlag = false;

    if(  (yA <= z && yB >= z) || (yB <= z && yA >= z) )
      iflagy = true;
    else
      iflagy = false;

    if(iflagy)
    {
      xtmp = xA;
      if(xtmp ==r) {x_hiFlag = true;  x_loFlag = true; }
      if(xtmp > r) {x_hiFlag = true;  x_loFlag = false;}
      if(xtmp < r) {x_hiFlag = false; x_loFlag = true; }
    }
    else
    {
      x_hiFlag = false;
      x_loFlag = false;
    }
    return(true);
  }

  // is this a horizontal line?
  if( yB == yA ) AA = 0.0;
  else           AA = (yB-yA)/ (xB-xA);

  BB = yA - AA * xA;

  if((xA <= r && xB >= r)||(xB <= r && xA >= r))
    iflagx = true;
  else
    iflagx = false;

  if((yA <= z && yB >= z)||(yB <= z && yA >= z))
    iflagy = true;
  else
    iflagy = false;

  if(iflagx)
  {
    ytmp = AA*r + BB;
    if(ytmp ==z) {y_hiFlag = true; y_loFlag = true;}
    if(ytmp > z) {y_hiFlag = true; y_loFlag = false;}
    if(ytmp < z) {y_hiFlag = false; y_loFlag = true;}
  }
  else
  {
    y_hiFlag = false;
    y_loFlag = false;
  }

  if(iflagy)
  {
    xtmp = (z-BB)/AA;
    if(xtmp ==r) {x_hiFlag = true; x_loFlag = true;}
    if(xtmp > r) {x_hiFlag = true; x_loFlag = false;}
    if(xtmp < r) {x_hiFlag = false; x_loFlag = true;}
  }
  else
  {
    x_hiFlag = false;
    x_loFlag = false;
  }

  return(true);
}

} // namespace Device
} // namespace Xyce

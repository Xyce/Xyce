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
// Purpose        : This is the class for mesh processing/ownership.
//                  of two dimensional meshes.
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

#ifndef Xyce_N_DEV_PDE_2DMesh__h
#define Xyce_N_DEV_PDE_2DMesh__h

// ----------   Standard Includes   ----------
#include <list>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_PDEMeshContainer.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_PDE_Electrode.h>

// ----------   Preprocessor Defines ----------

// these are used by function computeIntPB:
#define TAG12 0
#define TAG23 1
#define TAG13 2

#define VERTEX_A        0               // vertex A
#define VERTEX_B        1               // vertex B
#define VERTEX_C        2               // vertex C
#define VERTEX_D        3               // vertex D

#define EDGE_AB         0               // edge AB
#define EDGE_BC         1               // edge BC
#define EDGE_AC         2               // edge AC
#define EDGE_CD         2               // edge CD
#define EDGE_AD         3               // edge AD

#define TYPE_EDGE             7         // edge label
#define TYPE_REGION           8         // region label


namespace Xyce {
namespace Device {

// ----------   Forward Declarations ----------
class mNode;
class mEdge;
class mCell;
class mLabel;
class mNodeInfo;
class mEdgeInfo;


class NADJ;
class MESHHEAD;
class XLATCONST;
class XLATLABEL;
class NODE;
class EDGE;
class TRI;
class EDGEINFO;
class NODEINFO;

// ---------- Enum Definitions ----------
// mesh types:
enum meshType {
  INTERNAL,      //
  EXTERNAL,      //
  NUMTYPE        // total number of 2D mesh types
};

//-----------------------------------------------------------------------------
// Class         : PDE_2DMesh : public PDEMeshContainer
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class PDE_2DMesh : public PDEMeshContainer
{
public:
  PDE_2DMesh (const DeviceOptions & do1, int sgplotLevel1);

  PDE_2DMesh (const PDE_2DMesh & right);

  ~PDE_2DMesh ();

  // Assignment operator.
  PDE_2DMesh & operator=(PDE_2DMesh const & rhsMesh);

  bool initializeMesh (const std::string & meshFileName_tmp);

  bool initializeInternalMesh
  (int nx, int ny,
   double xlength, double ylength,
   int numElectrodes, std::string & outputMeshFileName,
   std::map<std::string,PDE_2DElectrode*> & elMap,
   bool cylFlag);

  bool resizeMesh (double xlength, double ylength);

  // output mesh information:
  void dumpMesh       ();
  void printLabels    ();
  void outputMeshInfo ();

  int getNumNodes ();
  int getNumEdges ();
  int getNumCells ();
  int getNumLabels();

  int getMaxNodeNN ();

  double getMaxSize ();

  double getXMax ();
  double getXMin ();
  double getYMax ();
  double getYMin ();

  double interp (double *F, double r, double z);
  bool interpVector (double *F, double r, double z, double & xvec, double & yvec);

  void findCell(double r, double z, int & isuccess,
                int & inode, int & iCell, int iStartCell = 0);

  double findMinDist(int iCell, double r, double z);

  double compAngle( double x1, double y1,
                    double x2, double y2,
                    double x3, double y3);

  bool scaleMesh (double xScale);

  bool labelNameExist (std::string & labelName);
  bool labelEdgeType  (std::string & labelName);
  bool dopingVectorExist ();

  bool getDopingVector (std::vector<double> & cvec_tmp);
  bool getXVector      (std::vector<double> & xvec_tmp);
  bool getYVector      (std::vector<double> & yvec_tmp);

  double * getDopingVector ();
  double * getXVector ();
  double * getYVector ();

  mNode * getNode (int i);
  mEdge * getEdge (int i);
  mCell * getCell (int i);
  mLabel * getLabel (int i);
  mLabel * getLabel (std::string & name);

  double lengthAdjust   (double x1, double y1, double x2, double y2);

  int ** getNodeIndexVector ();

protected:
private:
  bool readSGFMeshFile (const std::string & meshFileName_tmp);

  bool writeSGFMeshFile (const std::string & meshFileName_tmp);

  bool setupInternalMesh (int nx, int ny, double xlength, double ylength);

  bool setupInternalAdjacencyInfo ();

  bool errorCheckElectrodes (int numElectrodes,
                             std::map<std::string,PDE_2DElectrode*> & elMap);

  bool setupDefaultLabels (int numberElectrodes);

  bool setupInternalLabels (int numberElectrodes,
                            std::map<std::string,PDE_2DElectrode*> & elMap);

  bool setupGeometry    ();

  bool computeIntPB     (double &x, double &y,
                         int inodeA,int inodeB,int inodeC);

  double areaAdjust     (double x1, double y1,
                         double x2, double y2,
                         double x3, double y3);

  double computeAngle   (int inode1,int inode2,int inode3);

  bool cellNodes        ();

  bool fCCWorder        (int inode1, int inode2, int inode3);

  void calcAdjacencyInfo ();

  void initNodeAdjStructure
  (NADJ &nadj, int itri, int iVertex, int uIntLabel, bool fCW);

  void getElementInfo
  (int itri, int *ainode, int *aiedge, int *aitri, int *auLabel);

  void elementNodes (int itri, int *ainode);

public:
  bool cylGeom;   // cylindrical geometry flag. refactor later.

protected:


private:
  PDE_2DMesh ();

  std::string meshFileName;

  bool externalMeshFlag;

  double xMax;
  double yMax;
  double xMin;
  double yMin;
  double dx;
  double dy;

  double xRatio, yRatio;

  double x0;
  bool   meshScaledFlag;

  double vol;       // total volume of domain
  double invVol;    // inverse of the total volume

  double surfArea;  // surface area
  double circum;    // total circumference
  double invCircum; // inverse of the total circumference.

  double depth;     // if this is a cartesian mesh, depth in z-direction.

  int numAdj;       // number of edge adjacency structures
  // (>= numNodes.  Every mNode class has at least one.)

  int numNodes;
  int numEdges;
  int numCells;
  int numLabels;
  int numRegLabels;  // Number of region labels.  (no edge labels, etc.)

  int numBndryNodes; // This is the number of nodes that sit on
  // boundaries between regions.

  int maxNodeNN;     // maximum number of nearest neighbors for a node.

  int iRecentCellLookup; // This var is the result of the most
  // recent "findCell" call (if that function
  // has indeed been called).  It can be used
  // for subsequent findCell calls, as the
  // iStartCell argument.  The findCell function
  // searches over the mesh, proceeding from
  // cell to cell via nearest neighbors.
  // If it has a good starting point, the search goes
  // much, much faster.

  bool dopingSet; // doping array was read in from mesh file.

  std::vector<mNode> mNodeVector;
  std::vector<mEdge> mEdgeVector;
  std::vector<mCell> mCellVector;
  std::vector<mLabel> mLabelVector;

  std::vector<double> dopingVector;
  std::vector<double> xVector;
  std::vector<double> yVector;

  std::vector<int> visitCellFlagVec;

  std::map<std::string, mLabel>  mLabelMap;

  // used for the internally generated mesh.
  int ixMax;
  int iyMax;
  int ** nodeIndices;
  int ** edgeIndices;
  int ** cellIndices;

  // arrays used by calcAdjacencyInfo:
  std::vector<int> afVisitedVec;
  int * aiBegin;  // change these to STL later...
  int * aiEnd;
  bool adjInfoAllocFlag;

  const DeviceOptions * devOptions_;

  int sgplotLevel;

  bool useDefaultLabels;
};

//-----------------------------------------------------------------------------
// Much of what follows are 2D mesh primitives.  If I ever get around to doing
// things in 3D I will reuse as much of this as possible, and rename things
// that can't be re-used (to have either the string "2D" or "3D" somewhere
// in the name).
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Class         : mNode
// Purpose       : mesh node class
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class mNode
{
public:
  mNode ();

protected:
private:

public:
  double  x;                 // x coordinate
  double  y;                 // y coordinate

  // node adjacency information:
  double area;               // area of integration box
  int cnode;                 // number of adjacent nodes
  int inode;                 // index of node under analysis
  int numCells;              // number of cells

  int  edgeStatus;           // if boundary between regions = 0
  // if exterior then = 1
  // if interior then = 2

  bool fBndry;               // set if node is a boundary node.
  bool fGotAll;              // set if all nodes were visited.

  std::vector<EDGEINFO> edgeInfoVector;

protected:
private:
};

//-----------------------------------------------------------------------------
// Class         : mEdge
// Purpose       : mesh edge class
//
// Special Notes : "1" and "A" are conceptually the same.
//                 "2" and "B" are conceptually the same.
//
//                 For box integration, the "local" node is always node A,
//                 while the neighbor node is always node B.
//
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class mEdge
{
public:
  mEdge ();

protected:
private:

public:
  int   uLabel;               // edge label
  int   inodeA;               // index of node A
  int   inodeB;               // index of node B

  int   edgeStatus;          // if boundary between regions = 0
  // if exterior then = 1
  // if interior then = 2

  // edge geometry information:
  double  ilen;                           // integration length
  double  elen;                           // edge length
  double  Area1;                          // partial area 1, nodeA
  double  Area2;                          // partial area 2, nodeB

  double  angle;                          // angle between edge and x-axis.
  double  midpoint_x;
  double  midpoint_y;

  int     iedge;                          // edge index
  int     ielem;                          // element index

protected:
private:
};

//-----------------------------------------------------------------------------
// Class         : mCell
// Purpose       : mesh cell class.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class mCell
{
public:
  mCell ();

protected:
private:

public:
  int   uLabel;               // region label

  int   iedgeAB;              // index of edge AB
  int   iedgeBC;              // index of edge BC
  int   iedgeCD;              // index of edge AC or CD
  int   iedgeDA;              // index of edge DA

  int   icellAB;              // index of cell adj. to edge AB
  int   icellBC;              // index of cell adj. to edge BC
  int   icellCD;              // index of cell adj. to edge AC or CD
  int   icellDA;              // index of cell adj. to edge DA

  // owned nodes:
  int   inodeA;
  int   inodeB;
  int   inodeC;
  int   inodeD;

  std::vector<int> mNodeVector;    // container of owned node indices.

protected:
private:
};


//-----------------------------------------------------------------------------
// Class         : mLabel
// Purpose       : mesh label class.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class mLabel
{
public:
  mLabel ();

protected:
private:

public:
  std::string name;      // label name
  int iIndex;       // label index
  int uType;        // label type  (region or edge)
  int cNode;        // number of nodes

  double vol;       // volume of this region. (if it is a region...)
  double surfArea;  // surface area of this region.

  std::vector<int> mNodeVector;
};


//-----------------------------------------------------------------------------
// Class         : mEdgeInfo
// Purpose       : This class contains edge information for a single edge.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class mEdgeInfo
{
public:
protected:
private:

public:
  double  ilen;              // integration length
  double  elen;              // edge length
  double  area1;             // partial area 1
  double  area2;             // partial area 2
  int   iedge;               // edge index
  int   inode;               // node index
  int   icell;               // cell index

protected:
private:
};


//-----------------------------------------------------------------------------
// Class         : mNodeInfo
// Purpose       : This class contains node information for a single node.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class mNodeInfo
{
public:
protected:
private:

public:
  double area;               // area of integration box
  int numNeighbors;          // number of neighbors
  int numCells;              // number of cells
  mEdgeInfo *mEdgeInfoPtr;   // edge information array

protected:
private:
};

//-----------------------------------------------------------------------------
// What follows are SGF mesh primitives.  Most of these are only used for
// reading in SGF style mesh files, but some  are used for calculating
// geometric mesh information.
//-----------------------------------------------------------------------------

#define LEN_IDENT           15          // identifier length

#define EDGESTATUS_BOUNDARY 0
#define EDGESTATUS_EXTERIOR 1
#define EDGESTATUS_INTERIOR 2

//-----------------------------------------------------------------------------
// Class         : MESHHEAD
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
// neighbor adjacency structure
class NADJ            // nadj
{
public:
  int   cnode;        // number of adjacent nodes
  int   inode;        // index of node under analysis
  bool  fBndry;       // set if node is a boundary node
  bool  fGotAll;      // set if all nodes were visited
  int   ainode[32];   // indices of neighboring nodes
  int   aiedge[32];   // indices of "spoke" edges
  int   aielem[32];   // triangle indices
  int   auLabel[32];  // triangle labels
};

//-----------------------------------------------------------------------------
// Class         : MESHHEAD
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class MESHHEAD
{
public:
  char szLogo[64];                      // text logo
  char szSign[16];                      // signature
  unsigned int cConstant;               // number of constants
  unsigned int cLabel;                  // number of labels
  unsigned int cArray;                  // number of 1D arrays
  unsigned int cRegLabel;               // number of region labels
  unsigned int cBndryNode;              // number of boundary nodes
  unsigned int cNode;                   // number of nodes
  unsigned int cEdge;                   // number of edges
  unsigned int cTriangle;               // number of triangles
  bool fCylGeom;                        // cylindrical geometry flag
};

//-----------------------------------------------------------------------------
// Class         : XLATCONST
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class XLATCONST
{
public:
  char szName[LEN_IDENT+1];             // constant name
  unsigned int uType;                   // constant type
  union {
    int  n;                             // value of integer constant
    double r;                           // value of real constant
  } data;
};

//-----------------------------------------------------------------------------
// Class         : XLATLABEL
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class XLATLABEL
{
public:
  char szName[LEN_IDENT+1];             // label name
  unsigned int iIndex;                  // label index
  unsigned int uType;                   // label type
  unsigned int cNode;                   // number of nodes
};

//-----------------------------------------------------------------------------
// Class         : NODE
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class NODE
{
public:
  double  x;                            // x coordinate
  double  y;                            // y coordinate
};

//-----------------------------------------------------------------------------
// Class         : EDGE
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class EDGE
{
public:
  unsigned int  uLabel;                 // edge label
  int   inodeA;                         // index of node A
  int   inodeB;                         // index of node B
};

//-----------------------------------------------------------------------------
// Class         : TRI
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class TRI
{
public:
  unsigned int  uLabel;                 // region label
  int   iedgeAB;                        // index of edge AB
  int   iedgeBC;                        // index of edge BC
  int   iedgeAC;                        // index of edge AC or CD
  int   iedgeAD;                        // index of edge AD
  int   itriAB;                         // index of tri. adj. to edge AB
  int   itriBC;                         // index of tri. adj. to edge BC
  int   itriAC;                         // index of tri. adj. to edge AC or CD
  int   itriAD;                         // index of tri. adj. to edge AD
};

//-----------------------------------------------------------------------------
// Class         : EDGEINFO
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class EDGEINFO
{
public:
  double  ilen;                         // integration length
  double  elen;                         // edge length
  double  Area1;                        // partial area 1
  double  Area2;                        // partial area 2
  int   iedge;                          // edge index
  int   inode;                          // node index
  int   ielem;                          // element index

  EDGEINFO() :
    ilen(0.0), elen(0.0), Area1(0.0), Area2(0.0), iedge(-1), inode(-1), ielem(-1)
  {};

};


//-----------------------------------------------------------------------------
// Class         : NODEINFO
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class NODEINFO
{
public:
  double Area;                          // area of integration box
  unsigned int cNeighbor;               // number of neighbors
  unsigned int cTriangle;               // number of triangles
  EDGEINFO *aedgeinfo;                  // edge information array
};


//-----------------------------------------------------------------------------
// Class         : mInterpAreaHelp
// Purpose       : This is a helper class for performing interpolations
//                 on a mesh.  It contains information about 3 points in space
//                 and performs a linear interpolation between them.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
class mInterpAreaHelp
{
public:
  double x0,y0;
  double x1,y1;
  double x2,y2;
  double v0,v1,v2;
  double f0,f1,f2;
  double vlim;
  double aa, bb, cc;
  int errorFlag;
  int iend;

public:
  mInterpAreaHelp ();
  double interpReg(double r, double z);
  bool findCoef();
};

//-----------------------------------------------------------------------------
// Class         : mInterpEdgeHelp
// Purpose       : This is a helper class for performing interpolations
//                 on a mesh.  It contains information about an edge, including
//                 the linear equation for it.  y = AA*x + BB.
//
//                 It also determines this edge's relationship to a point in
//                 space.  An interpolation routine can use this help class
//                 to determine if a point in space falls within a mesh cell
//                 or not.
//
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 09/18/02
//-----------------------------------------------------------------------------
class mInterpEdgeHelp
{
public:
  double xA;
  double yA;
  double xB;
  double yB;

  double AA;
  double BB;

  bool iflagx; // is the passed point between xA and xB of this edge?
  bool iflagy; // is the passed point between yA and yB of this edge?

  bool x_hiFlag; // is this edge above the passed point in x?
  bool x_loFlag; // is this edge below the passed point in x?
  bool y_hiFlag; // is this edge above the passed point in y?
  bool y_loFlag; // is this edge above the passed point in y?

public:
  mInterpEdgeHelp();
  bool setupEdge (double r, double z);

};


//------------------------ Inline  functions ----------------------------------


inline int PDE_2DMesh::getNumNodes () { return numNodes; }

inline int PDE_2DMesh::getNumEdges () { return numEdges; }

inline int PDE_2DMesh::getNumCells () { return numCells; }

inline int PDE_2DMesh::getNumLabels() { return numLabels; }

inline int PDE_2DMesh::getMaxNodeNN () { return maxNodeNN; }

inline double PDE_2DMesh::getMaxSize ()
{ return (((yMax-yMin)>(xMax-xMin))?yMax:xMax);  }

inline double PDE_2DMesh::getXMax () {return xMax; }
inline double PDE_2DMesh::getXMin () {return xMin; }
inline double PDE_2DMesh::getYMax () {return yMax; }
inline double PDE_2DMesh::getYMin () {return yMin; }

inline bool PDE_2DMesh::dopingVectorExist ()
{ return dopingSet; }

inline double * PDE_2DMesh::getDopingVector ()
{
  return &(dopingVector[0]);
}

inline double * PDE_2DMesh::getXVector ()
{
  return &(xVector[0]);
}

inline double * PDE_2DMesh::getYVector ()
{
  return &(yVector[0]);
}

inline mNode * PDE_2DMesh::getNode (int i)
{
  return &(mNodeVector[i]);
}

inline mEdge * PDE_2DMesh::getEdge (int i)
{
  return &(mEdgeVector[i]);
}

inline mCell * PDE_2DMesh::getCell (int i)
{
  return &(mCellVector[i]);
}

inline mLabel * PDE_2DMesh::getLabel (int i)
{
  return &(mLabelVector[i]);
}

inline int ** PDE_2DMesh::getNodeIndexVector ()
{
  return nodeIndices;
}

//-----------------------------------------------------------------------------
// Function      : sq
// Purpose       :
// Special Notes : put this somewhere else later...
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/22/02
//-----------------------------------------------------------------------------
inline double sq(double x)
{
  return(x * x);
}

//-----------------------------------------------------------------------------
// Function      : PDE_2DElectrode::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, SNL, Parallel Computational Sciences
// Creation Date : 04/18/03
//-----------------------------------------------------------------------------
inline std::ostream & operator<<(std::ostream & os, const PDE_2DElectrode & el)
{
  os << el.name << ":\n";
  os << "  node  = " << el.nodeName << "\n";
  os << "  side  = " << el.side << "\n";
  os << "  start = " << el.start << "\n";
  os << "  end   = " << el.end << "\n";
  os << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

#endif


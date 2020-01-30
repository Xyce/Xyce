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
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_TOP_fwd_h
#define Xyce_N_TOP_fwd_h

#include <list>
#include <vector>

#if defined(HAVE_UNORDERED_SET)
#include <unordered_set>
using std::unordered_set;
#elif defined(HAVE_TR1_UNORDERED_SET)
#include <tr1/unordered_set>
using std::tr1::unordered_set;
#else
#error neither unordered_set or tr1/unordered_set found
#endif

namespace Xyce {

enum NodeTYPE {_VNODE, _DNODE, _CNODE, _PNODE, _NUM_NODE_TYPES};

class NodeID;

namespace Topo {

class CktGraph;
class CktGraphBasic;
class CktGraphCreator;
class CktGraphCreatorBasic;
class CktGraphSupport;
class CktNode;
class CktNodeCreator;
class CktNode_Dev;
class CktNode_V;
class Directory;
class Indexor;
class Graph;
class Manager;
class Node;
class NodeDevBlock;
class TopoLSUtil;
class LSUtilFactory;
class SerialLSUtil;
class ParLSUtil;
class Topology;
class System;

typedef std::vector<CktNode *> CktNodeList;
typedef std::list<CktGraph *> CktGraphList;

} // namespace Topo
} // namespace Xyce

// typedef Xyce::NodeID NodeID;

#endif // Xyce_N_TOP_fwd_h

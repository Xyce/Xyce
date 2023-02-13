//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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
// Purpose        : Simple undirected graph class
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Electrical & MicroSystems
//
// Creation Date  : 8/10/06
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_UTL_Graph_H
#define Xyce_UTL_Graph_H

#include <vector>
#include <unordered_map>
using std::unordered_map;
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include <iostream>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Class         : Cmp
// Purpose       : Simple class to provide a comparison for removing keys from the
//                 Graph class
// Special Notes :
// Scope         : Public
// Creator       : Heidi K. Thornquist, SNL
//-----------------------------------------------------------------------------
template <typename Index>
class Cmp
{
public:
  Cmp( const std::vector<Index> & inputIDs )
    : numIDs_(inputIDs.size()), inputIDs_(inputIDs)
  {}

  int binary_search( const std::vector<Index> & array, Index id, int low, int high )
  {
    if (high < low)
      return -1;
    int mid = (low+high) / 2;
    if (array[mid] > id)
      return binary_search( array, id, low, mid-1 );
    else if (array[mid] < id)
      return binary_search( array, id, mid+1, high );
    else
      return mid;
  }

  // Comparison operator/functor
  bool operator()( Index id )
  {
    bool ret = false;
    int idx = binary_search( inputIDs_, id, 0, numIDs_-1 );
    if (idx > -1)
      ret = true;

    return ret;
  }

  private:

  int numIDs_;
  const std::vector<Index> inputIDs_;
};

template <typename Index>
class Cmp1
{
public:
  Cmp1( const Index inputID )
    : id_(inputID)
  {}

  bool operator()( Index id )
  {
    return ( (id == id_) ? true : false );
  }

  private:
  
  Index id_;

};


//-----------------------------------------------------------------------------
// Class         : Graph
// Purpose       : Simple undirected graph
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
class Graph
{
public:
  typedef unordered_map<Index, Key1Type> Key1Map;
  typedef unordered_map<Key1Type, Index> Index1Map;
  typedef unordered_map<Key1Type, DataType> Data1Map;
  typedef std::vector< std::vector<Index> > AdjacencyGraph;

  // Constructors
  Graph() : numRemovedNodes_(0) {}

  // Destructor
  virtual ~Graph() {}

private:
  Graph(const Graph&);
  Graph& operator=(const Graph&);
  int operator==(const Graph& right) const;
  int operator!=(const Graph& right) const;

public:
  bool insertNode(const Key1Type &key1, const std::vector<Key1Type> &adj, DataType &data);

  int numNodes() const;

  bool checkKey(const Key1Type& key) const;

  DataType& getData(const Key1Type& key);
  DataType& getData(const Index& idx);
  std::vector<DataType>& getData(const std::vector<Key1Type>& keys);

  int numAdjacent(const Key1Type& key) const;

  std::vector<Key1Type> getAdjacent(const Key1Type& key) const;

  const Key1Type& getKey(const Index& idx) const {
    typename Key1Map::const_iterator it = keys1_.find( idx );
    if (it == keys1_.end()) 
      throw std::runtime_error("Graph index not found");
    
    return (*it).second;
  }

  Index getIndex(const Key1Type & key) const {
    typename Index1Map::const_iterator it = rvsKeys1_.find(key);
    if (it == rvsKeys1_.end())
      throw std::runtime_error("Graph key not found");

    return (*it).second;
  }

  const std::vector<Index> & getAdjacentRow(const Index idx) const { return adjacencyGraph_[idx]; }

  void addToAdjacent(const Key1Type& oldkey, const Key1Type& key, std::vector<Key1Type> & newAdjVec);

  void replaceAdjacent(const Key1Type & oldKey1, const Key1Type & newKey1);

  void removeKey(const Key1Type &oldKey1);

  void removeKeys(const std::vector<Key1Type> & oldKeys1);

  std::vector<Key1Type> getSingletons();

  int checkGraphState() const;

  // Give access to the reverse index map 
  const Index1Map &getIndex1Map() { return rvsKeys1_; }

  // Give access to the data map in the event that the nodes need
  // to be traversed without an ordered list (BFT) being generated.
  const Data1Map &getData1Map() { return data1_; }

  const std::vector<Index>& getBFT();
  int generateBFT();
  int generateBFT(const Key1Type& key);

  void print(std::ostream & ostr) const;

private:
  int generateBFT_(const Index& start);

private:
  AdjacencyGraph        adjacencyGraph_;

  // if we remove any nodes from the adjacencyGraph_, then we'll have empty
  // rows.  We could erase these, but then we would have to reindx the keymaps.
  // For now we'll just count the number of removals
  int                   numRemovedNodes_;

  Key1Map               keys1_;
  Index1Map             rvsKeys1_;
  Data1Map              data1_;

  std::vector<Index>    bft_;
};

//-----------------------------------------------------------------------------
// Function      : Graph::insertNode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
bool
Graph<Key1Type, DataType, Index>::insertNode(
  const Key1Type &              key1,
  const std::vector<Key1Type> & adj,
  DataType &                    data)
{
  if (checkKey(key1))
    return false;

  std::vector<Index> ids;
  size_t size = adj.size();

  for (size_t i = 0; i < size; ++i)
    ids.push_back(rvsKeys1_[adj[i]]);

  adjacencyGraph_.push_back(ids);

  // add back edges
  Index currIndex = adjacencyGraph_.size() - 1;
  for (size_t i = 0; i < size; ++i)
  {
    std::vector<Index>& adjVec = adjacencyGraph_[ ids[i] ];
    adjVec.push_back( currIndex );
  }

  keys1_[currIndex] = key1;
  rvsKeys1_[key1] = currIndex;
  data1_[key1] = data;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Graph::numNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline int Graph<Key1Type, DataType, Index>::numNodes() const
{
  return adjacencyGraph_.size() - numRemovedNodes_;
}

//-----------------------------------------------------------------------------
// Function      : Graph::checkKey
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline bool Graph<Key1Type, DataType, Index>::checkKey
  (const Key1Type& key) const
{
  if(rvsKeys1_.count(key)) return true;
  return false;
}

//-----------------------------------------------------------------------------
// Function      : Graph::getData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
DataType &
Graph<Key1Type, DataType, Index>::getData(
  const Index&      idx)
{
  return data1_[keys1_[idx]];
}

//-----------------------------------------------------------------------------
// Function      : Graph::getData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
DataType &
Graph<Key1Type, DataType, Index>::getData(
  const Key1Type &      key)
{
  return data1_[key];
}

//-----------------------------------------------------------------------------
// Function      : Graph::getData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline std::vector<DataType>& Graph<Key1Type, DataType, Index>::getData
  (const std::vector<Key1Type>& keys)
{
  std::vector<DataType> data;
  for(int i = 0; i < keys.size(); ++i)
    data.push_back(data1_[keys[i]]);

  return data;
}

//-----------------------------------------------------------------------------
// Function      : Graph::numAdjacent
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
int
Graph<Key1Type, DataType, Index>::numAdjacent(
  const Key1Type &      key) const
{
  Index id = rvsKeys1_[key];
  return adjacencyGraph_[id].size();
}

//-----------------------------------------------------------------------------
// Function      : Graph::getAdjacent
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
std::vector<Key1Type>
Graph<Key1Type, DataType, Index>::getAdjacent(
  const Key1Type &      key) const
{
  std::vector<Key1Type> adjKeys;

  typename Index1Map::const_iterator it = rvsKeys1_.find(key);
  if (it == rvsKeys1_.end())
    return adjKeys;

  Index id = (*it).second;
  int size = adjacencyGraph_[id].size();
  for (int i = 0; i < size; ++i) {
    typename Key1Map::const_iterator it = keys1_.find(adjacencyGraph_[id][i]);
    if (it != keys1_.end())
      adjKeys.push_back((*it).second);
  }

  return adjKeys;
}

//-----------------------------------------------------------------------------
// Function      : Graph::addToAdjacent
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
void
Graph<Key1Type, DataType, Index>::addToAdjacent(
  const Key1Type &              oldkey,
  const Key1Type &              key,
  std::vector<Key1Type> &       newAdjVec)
{
  // originally I had this as setAdjacent() where it would clear and then
  // write a new adjacency info to adjacencyGraph_.
  // A flaw with using this as setAdjacent is that it can change the order of
  // edge nodes.  potentially reversing an element on another part of the circuit.
  // so rather than erasing and then setting a new adjacent list, we'll just add to what
  // is there

  // add in new values
  int extraElementsSize = newAdjVec.size();
  if (extraElementsSize)
  {
    Index id = rvsKeys1_[key];
    Index oldId = rvsKeys1_[oldkey];
    for(int i=0; i< extraElementsSize; i++)
    {
      Index edgeIndex = rvsKeys1_[newAdjVec[i]];
      adjacencyGraph_[id].push_back( edgeIndex );

      // and any new edges
      // need to do this by changing the edge's oldId to the id of the new key
      typename std::vector< Index >::iterator endLoc = adjacencyGraph_[ edgeIndex ].end();
      typename std::vector< Index >::iterator edgeLoc = find(adjacencyGraph_[ edgeIndex ].begin(), endLoc, oldId);
      if(edgeLoc == endLoc)
      {
        adjacencyGraph_[edgeIndex].push_back(id);
      }
      else
      {
        *edgeLoc = id;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Graph::replaceAdjacent
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
void
Graph<Key1Type, DataType, Index>::replaceAdjacent(
  const Key1Type &      oldKey1,
  const Key1Type &      newKey1)
{
  // get id's for these keys
  Index oldId = rvsKeys1_[ oldKey1 ];
  Index newId = rvsKeys1_[ newKey1 ];

  // loop through adjancy graph and replace any references to oldKey with
  // newKey only if newKey doesn't already exist on that row
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    typename std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
    typename std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
    // look for the old key
    typename std::vector< Index >::iterator oldKeyLoc = find(beginLoc, endLoc, oldId);
    if(oldKeyLoc != endLoc)
    {
      // found old key, so overwrite it
      *oldKeyLoc = newId;
      //
      // original logic removed the old one and added a new one
      // this could break ordering of the id's which would be bad
      // // found old key, so erase it
      // adjacencyGraph_[i].erase(oldKeyLoc);
      // // also check for new key before inserting it
      // std::vector< Index >::iterator newKeyLoc = find(beginLoc, endLoc, newId);
      // if(newKeyLoc == endLoc)
      // {
      //   // new key isn't there so add it in
      //   adjacencyGraph_[i].push_back(newId);
      // }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Graph::removeKey
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
void
Graph<Key1Type, DataType, Index>::removeKey(
  const Key1Type &      oldKey1)
{
  Index id = -1;

  typename Index1Map::const_iterator it = rvsKeys1_.find(oldKey1);
  if (it != rvsKeys1_.end())
  {
    // get keys and id
    id = it->second;

    // clear the row it uses
    // this leaves an empty row in place, but that's ok
    adjacencyGraph_[id].clear();

    // erase row that this key points to in adjacency graph
    // std::vector< std::vector< Index > >::iterator oldRowItr = adjacencyGraph_.begin() + id;
    // adjacencyGraph_.erase(oldRowItr);
    // if I erase the row it throws off the key maps.  So I need to fix up those

    // search through adjacencyGraph_ and remove all references to this key
    // then remove it from the maps as well
    int numAdjRows = adjacencyGraph_.size();
    for(int i=0; i<numAdjRows; ++i)
    {
      if(!adjacencyGraph_[i].empty())
      {
        typename std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
        typename std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
        adjacencyGraph_[i].erase( remove_if(beginLoc, endLoc, Cmp1<Index>(id)), endLoc);
      }
    }

    // clean up maps
    keys1_.erase(id);         // Key1Map keys1_;
    rvsKeys1_.erase(oldKey1); // Index1Map rvsKeys1_;
    data1_.erase(oldKey1);    // Data1Map data1_;
  
    numRemovedNodes_++;
  }
}

//-----------------------------------------------------------------------------
// Function      : Graph::removeKeys
// Purpose       :
// Special Notes : This is a collective removal from the adjacencyGraph for efficiency.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
void
Graph<Key1Type, DataType, Index>::removeKeys(
  const std::vector<Key1Type> &         oldKeys1)
{
  int numKeysToRemove = oldKeys1.size();
  std::vector<Index> ids;

  // Go through all the keys, retrieve their ids and clean up the maps
  for(int i=0; i<numKeysToRemove; ++i)
  {
    typename Index1Map::const_iterator it = rvsKeys1_.find(oldKeys1[i]);
    if (it != rvsKeys1_.end())
    {
      // get keys and id
      ids.push_back( it->second );

      // clear the row it uses
      // this leaves an empty row in place, but that's ok
      adjacencyGraph_[it->second].clear();

      // clean up maps
      keys1_.erase( it->second );         // Key1Map keys1_;
      rvsKeys1_.erase( oldKeys1[i] ); // Index1Map rvsKeys1_;
      data1_.erase( oldKeys1[i] );    // Data1Map data1_;
    }
  }

  if (ids.size())
  {
    // Sort the list of ids to expedite searching
    std::sort(ids.begin(), ids.end());

    // search through adjacencyGraph_ and remove all references to this key
    int numAdjRows = adjacencyGraph_.size();
    for(int i=0; i<numAdjRows; ++i)
    {
      if(!adjacencyGraph_[i].empty())
      {
        typename std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
        typename std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
        adjacencyGraph_[i].erase( remove_if(beginLoc, endLoc, Cmp<Index>(ids)), endLoc);
      }
    }
  }

  numRemovedNodes_ += ids.size();
}

//-----------------------------------------------------------------------------
// Function      : Graph::getSingletons
// Purpose       :
// Special Notes : This method returns the viable nodes that are graph singletons.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
std::vector<Key1Type>
Graph<Key1Type, DataType, Index>::getSingletons()
{
  std::vector<Key1Type> singletonKeys;

  int numAdjRows = adjacencyGraph_.size();
  for (int i = 0; i < numAdjRows; ++i)
  {
    if (adjacencyGraph_[i].empty() && keys1_.count(i) > 0)
      singletonKeys.push_back(keys1_[i]);
  }

  return singletonKeys;
}

//-----------------------------------------------------------------------------
// Function      : Graph::checkGraphState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
int
Graph<Key1Type, DataType, Index>::checkGraphState() const
{
  // this is for debugging
  // it checks that the maps are consisitent (i.e. they all have the same
  // keys and indices).
  // then it traverses the adjacencyGraph_ and makes sure that it only
  // points to items that are in the key maps

  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    int numNodes = adjacencyGraph_[i].size();
    for(int j=0; j< numNodes; j++)
    {
      // check index in all maps
      Index testIndex = adjacencyGraph_[i][j];
      Key1Type key1val = keys1_[ testIndex ];
      if(keys1_.count(testIndex) == 0)
      {
        return 1;
      }
      if(rvsKeys1_.count(key1val) == 0)
      {
        return 2;
      }
      if(data1_.count(key1val) == 0)
      {
        return 3;
      }
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : Graph::getBFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline 
const std::vector< Index >& 
Graph<Key1Type, DataType, Index>::getBFT()
{
  if(bft_.empty()) generateBFT();
  return bft_;
}

//-----------------------------------------------------------------------------
// Function      : Graph::generateBFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
int
Graph<Key1Type, DataType, Index>::generateBFT()
{
#ifdef Xyce_GRAPH_DEBUG
  int state = checkGraphState();
  if(state != 0)
  {
    Xyce::dout() << "Graph is in inconsistent state!  checkGraphState returned: " << state << std::endl;
  }
#endif
 
/* 
  // Select the first maximum degree node as the center of the graph
  typename AdjacencyGraph::iterator it = adjacencyGraph_.begin();
  typename AdjacencyGraph::iterator it_end = adjacencyGraph_.end();
  int currNode = 0, maxDegNode = 0; 
  int maxDeg = (*it).size();
  for( ; it != it_end; ++it, ++currNode )
  {
    if ( (*it).size() > maxDeg )
    {
      maxDegNode = currNode;
      maxDeg = (*it).size();
    }
  }
*/
 
  int level = 0; 
  if (keys1_.size())
  {
    Index firstIndex = (*(keys1_.begin())).first;
    level = generateBFT_(firstIndex);
/*
    std::cout << "The degree of the first-key root node is " << adjacencyGraph_[firstIndex].size()
              << " and the BFT returned a maximum level of " << level << std::endl;

    firstIndex = maxDegNode;
    level = generateBFT_(firstIndex);
    std::cout << "The degree of the max-degree root node is " << adjacencyGraph_[firstIndex].size()
              << " and the BFT returned a maximum level of " << level << std::endl;
*/
  }
  return level;
}

//-----------------------------------------------------------------------------
// Function      : Graph::generateBFT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
int
Graph<Key1Type, DataType, Index>::generateBFT(
  const Key1Type &      key)
{
#ifdef Xyce_GRAPH_DEBUG
  int state = checkGraphState();
  if(state != 0)
  {
    Xyce::dout() << "Graph is in inconsistent state!  checkGraphState returned: " << state << std::endl;
  }
#endif

  Index id = rvsKeys1_[key];
  return generateBFT_(id);
}

//-----------------------------------------------------------------------------
// Function      : Graph::print
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
void
Graph<Key1Type, DataType, Index>::print(std::ostream & ostr) const
{
  ostr << "-------------------- Basic Graph ----------------------------\n";
  ostr << "Adjacency Graph\n";
  ostr << "---------------\n";
  for(size_t i = 0; i < adjacencyGraph_.size(); ++i)
  {
    ostr << "Node " << i << " : ";
    for(size_t j = 0; j < adjacencyGraph_[i].size(); ++j)
      ostr << " " << adjacencyGraph_[i][j];
    ostr << std::endl;
  }
  ostr << "---------------\n";
  ostr << "Key1Map\n";
  for(typename Key1Map::const_iterator it_k1m = keys1_.begin(); it_k1m != keys1_.end(); ++it_k1m)
    ostr << it_k1m->first << ":" << it_k1m->second << std::endl;
  ostr << "-------\n";
  ostr << "Index1Map\n";
  for(typename Index1Map::const_iterator it_i1m = rvsKeys1_.begin(); it_i1m != rvsKeys1_.end(); ++it_i1m)
    ostr << it_i1m->first << ":" << it_i1m->second << std::endl;
  ostr << "-------\n";
  ostr << "Data1Map\n";
  for(typename Data1Map::const_iterator it_d1m = data1_.begin(); it_d1m != data1_.end(); ++it_d1m)
    ostr << it_d1m->first << ":" << it_d1m->second << std::endl;
  ostr << "-------\n";
  ostr << "BFT\n";
  for(size_t i = 0; i < bft_.size(); ++i)
    ostr << bft_[i] << std::endl;
  ostr << "-------\n";
  ostr << "-------------------- Basic Graph END ------------------------\n";
}

//-----------------------------------------------------------------------------
// Function      : Graph::generateBFT_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Robert J. Hoekstra, SNL
//               : Heidi Thornquist, SNL (modified search for efficiency)
// Creation Date :
//-----------------------------------------------------------------------------
template <typename Key1Type, typename DataType, typename Index>
inline
int
Graph<Key1Type, DataType, Index>::generateBFT_(
  const Index &         start)
{
  bft_.clear();

  // Work queue
  std::queue< std::pair<Index, int> > idQueue;

  // Keep track of which IDs have been found (false=not found, true=found)
  std::vector<bool> foundIds(adjacencyGraph_.size(), false);

  // Keep track of all the levels if the graph is a forest
  std::vector<int> levels;

  Index localCopyOfStart = start;

  // Before we push back "start" need to verify that it's valid
  bool startIsValid = false;
  int numAdjacentRows = adjacencyGraph_.size();
  while(!startIsValid)
  {
    if(keys1_.count(localCopyOfStart) > 0)
    {
      startIsValid=true;
    }
    else
    {
      localCopyOfStart++;
      if(localCopyOfStart > (numAdjacentRows-1))
      {
        localCopyOfStart = 0;
      }
    }
  }

  // Initialize level, root, and work queue
  int level = 0;
  unsigned int root = 0;
  idQueue.push(std::make_pair(localCopyOfStart, level));
  bft_.push_back(localCopyOfStart);
  foundIds[localCopyOfStart] = true;

  while(!idQueue.empty())
  {
    Index currId = idQueue.front().first;
    level = idQueue.front().second;
    idQueue.pop();

    typename std::vector< Index >::iterator it = adjacencyGraph_[currId].begin();
    typename std::vector< Index >::iterator it_end = adjacencyGraph_[currId].end();
    for( ; it != it_end; ++it )
    {
      Index adjId = *it; 
      if(!foundIds[adjId])
      {
        idQueue.push(std::make_pair(adjId, level+1));
        bft_.push_back(adjId);
        foundIds[adjId] = true;
      }
    }

    // If adjacencyGraph_.size() isn't the full size of the problem, then there is a forest.
    // At this point, the level should be reset to 0, the current level should be stored in "levels", 
    // and the search for the next viable root should start at "root".
    if(idQueue.empty() && ((int)bft_.size()!=numNodes()))
    {
      levels.push_back(level);    // Store level of completed tree.
      level = 0;                  // Initialize level since we are ordering a new tree.

      // Start from root, since we know everything before root has been seen
      for( ; root < adjacencyGraph_.size(); ++root)
      {
        // we only count need to find Id's for places on the graph with greater than zero size
        if(adjacencyGraph_[root].size() > 0 && !foundIds[root])
        {
          idQueue.push(std::make_pair(root, level));
          bft_.push_back(root);
          foundIds[root] = true;
          break;
        }
      }
    }
  }

  // The returned level is the maximum of all the trees, if there are multiple trees.
  if (levels.size() > 1)
  {
    typename std::vector<Index>::iterator max_level = std::max_element( levels.begin(), levels.end());
    level = *max_level;
  }

  return level;
}

} // namespace Util
} // namespace Xyce

#endif

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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra
//
// Creation Date  : 7/16/01
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef _N_UTL_Functors_h
#define _N_UTL_Functors_h 1

// ---------- Standard Includes ----------
#include <functional>
#include <map>

// ----------   Xyce Includes   ----------

//-----------------------------------------------------------------------------
// Class         : DeletePtr
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/16/01
//-----------------------------------------------------------------------------
template < typename T >
struct DeletePtr : public std::unary_function < const T *, void >
{
  void operator() (const T * ptr) const { delete ptr; }
};

//-----------------------------------------------------------------------------
// Class         : FirstOfPair
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/16/01
//-----------------------------------------------------------------------------
template < typename T, typename U >
struct FirstOfPair : public std::unary_function < const T &, const U & >
{
  const U & operator() (const T & ref) const { return ref.first; }
};

//-----------------------------------------------------------------------------
// Class         : SortContainer2
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/04
//-----------------------------------------------------------------------------
template < typename T, typename U >
void SortContainer2( T & firstContainer, U & secondContainer )
{
  typedef typename std::multimap< typename T::value_type, typename U::value_type> UTMultiMap;

  UTMultiMap SortMap;

  typename T::iterator iterT = firstContainer.begin();
  typename T::iterator endT = firstContainer.end();
  typename U::iterator iterU = secondContainer.begin();
  typename U::iterator endU = secondContainer.end();

  for( ; (iterT!=endT)||(iterU!=endU) ; ++iterT, ++iterU )
    SortMap.insert( typename UTMultiMap::value_type( *iterT, *iterU ) );

  firstContainer.clear();
  secondContainer.clear();

  typename UTMultiMap::iterator iterUTM = SortMap.begin();
  typename UTMultiMap::iterator endUTM = SortMap.end();

  for( ; iterUTM != endUTM; ++iterUTM )
  {
    firstContainer.push_back( iterUTM->first );
    secondContainer.push_back( iterUTM->second );
  }
}

//-----------------------------------------------------------------------------
// Class         : IsSorted
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/04
//-----------------------------------------------------------------------------
template < typename T >
bool IsSorted( T & container )
{
  if( container.size() < 2 ) return true;

  typename T::iterator iterT = container.begin();
  typename T::iterator endT = container.end();
  typename T::iterator iterTPlus = iterT;
  iterTPlus++;

  for( ; iterTPlus != endT; ++iterT, ++iterTPlus )
    if( !(*iterT<*iterTPlus) ) return false;

  return true;
}

//-----------------------------------------------------------------------------
// Class         : N_UTL_Less 
// Purpose       : Class implemented for use with RedStorm PGI compilers
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Org. 1437
// Creation Date : 3/27/06
//-----------------------------------------------------------------------------
template < typename T, typename ST >
class LessThan : public std::less<ST>
{
   public:
     bool operator() ( const T& val1, const ST& val2 ) const { return val1 < val2; }
};

#endif


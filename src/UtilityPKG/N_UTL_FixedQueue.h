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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek
//
// Creation Date  : 9/4/2009
//
//
//
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_FixedQueue_h
#define Xyce_N_UTL_FixedQueue_h

#include <vector>

// ----------   Xyce Includes   ----------

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------
// Class         : FixedQueue
// Purpose       : This class implements a fixed length queue.  While the
//                 STL provides a queue class, it is designed to be variable
//                 in length so that further additions even if coupled with
//                 deletions to keep the length constant will still result
//                 in memory allocations.  This class's purpose is to hold
//                 a fixed number of items from a growing source -- and only
//                 hold on the the most recient ones.  For example, we may
//                 want to trank the last 100 time steps.  Once the queue is
//                 full, any further time steps that are added simply erase
//                 the oldest one present.  Thus there isn't any memory
//                 penality when trying to trach the last n events when 
//                 there can be millions of them.
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/04/09
//-------------------------------------------------------------------------

template <class Type> class FixedQueue
{
private:
  std::vector< Type > queueData_;
  int queueLength_;
  int headIndex_;

  // technically we don't need tailIndex_ as it's just headIndex_+1 once
  // the queue is full.  I could rewrite this to avoid and extra variable.
  int tailIndex_;
  bool queueFull_;

public: 
  FixedQueue( int queueSize=0 )
    : queueLength_( queueSize ), headIndex_(0), tailIndex_(0), queueFull_(false)
  {
    queueData_.resize( queueLength_ );
  }

  int get_size() const
  { 
    return queueLength_;
  }
  
  void set_size( int queueSize )
  {
    queueData_.resize( queueSize );
    queueLength_ = queueData_.size();
    reset_queue();
  }

  Type at_from_head( int index ) const
  {
    int modIndex = (headIndex_ - index) % queueLength_;
    //Xyce::dout() << "tailIndex_ = " << tailIndex_ << " headIndex_ = " 
    //  << headIndex_ << " queueLength_ = " << queueLength_ 
    //  << " index = " << index << " modIndex = " << modIndex << std::endl;
    return queueData_[ modIndex ];
  } 

  Type at_from_tail( int index ) const
  {
    int modIndex = (index + tailIndex_) % queueLength_;
    //Xyce::dout() << "tailIndex_ = " << tailIndex_ << " headIndex_ = " 
    //  << headIndex_ << " queueLength_ = " << queueLength_ 
    //  << " index = " << index << " modIndex = " << modIndex << std::endl;
    return queueData_[ modIndex ];
  } 

  void push_back( Type value )
  {
    // increment the head index
    headIndex_++;
    // roll over value if needed
    if( headIndex_ >= queueLength_ )
    {
      headIndex_ = 0;
      queueFull_ = true;
    }
    if( queueFull_ )
    {
      // if the queue is full, then adding
      // an element means we're going to forget
      // one as well from the tail
      tailIndex_++;
      if( tailIndex_ >= queueLength_ )
      {
        tailIndex_ = 0;
      }
    }
    // sorted out the index, so store 
    // the value
    queueData_[ headIndex_ ] = value;
  }

  void reset_queue()
  {
    headIndex_ = 0;
    tailIndex_ = 0;
    queueFull_ = false;
  }
};

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_FixedQueue_h

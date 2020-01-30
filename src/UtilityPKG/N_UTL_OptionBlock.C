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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Rob Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_ERH_Message.h>
#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : compareParamLists
// Purpose       : Compares contained Param lists
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 2/6/12
//-----------------------------------------------------------------------------
bool compareParamLists(const OptionBlock &s0, const OptionBlock &s1)
{
  bool match = true;
  if( s1.size() == s0.size() )
  {
    // length of params list is the same, so they may match
    ParamList::const_iterator existingListItr = s1.begin();
    ParamList::const_iterator existingListEnd = s1.end();
    ParamList::const_iterator newListItr = s0.begin();
    ParamList::const_iterator newListEnd = s0.end();
    while( existingListItr != existingListEnd )
    {
      if (!Util::deepCompare(*existingListItr, *newListItr))
      {
        // Param objects didn't match
        // thus these two lists are not the same
        match = false;
        break;
      }
      existingListItr++;
      newListItr++;
    }
  }
  else
  {
    match = false;
  }
  return match;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
OptionBlock & OptionBlock::operator=(const OptionBlock & right)
{
  if (this != &right)
  {
    name_   = right.name_;
    expressionsAllowed_ = right.expressionsAllowed_;
    netlistLocation_ = right.netlistLocation_;
    params_ = right.params_;
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::addParam
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date :
//-----------------------------------------------------------------------------
void OptionBlock::addParam(const Util::Param &parameter)
{
  if (expressionsAllowed_ == NO_EXPRESSIONS && parameter.hasExpressionValue())
    Report::UserError0().at(netlistLocation_) << "Expressions are not supported for " << name_;

  params_.push_back(parameter);
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::removeParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void OptionBlock::removeParam(const std::string& tagName)
{
  Util::ParamList::iterator it_tpL = params_.begin();
  for ( ; it_tpL!=params_.end(); )
  {
    if (it_tpL->tag() == tagName)
    {
      it_tpL = params_.erase( it_tpL );
    }
    else
    {
      it_tpL++;
    }  
  }
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const OptionBlock & mb)
{
  os << "Option Block" << std::endl;
  os << " name:  " << mb.getName() << std::endl;
  os << " Params" << std::endl;
  os << " -------------" << std::endl;

  for (Util::ParamList::const_iterator it = mb.begin(), end = mb.end(); it != end; ++it)
  {
    os << *it;
  }
  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

} // namespace Device

//-----------------------------------------------------------------------------
// Function      : OptionBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
template<>
int
Pack<Util::OptionBlock>::packedByteCount(
  const Util::OptionBlock &    option_block)
{
  int byteCount = 0;
  int size, length, i;

  //----- count name
  length = option_block.name_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count netlistLocation
  byteCount += 2*sizeof(int);

  //----- count params
  size = option_block.params_.size();
  byteCount += sizeof(int);
  Util::ParamList::const_iterator it_tpL = option_block.params_.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
    byteCount += Xyce::packedByteCount(*it_tpL);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::pack
// Purpose       : Packs OptionBlock into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Util::OptionBlock>::pack(const Util::OptionBlock &option_block, char * buf, int bsize, int & pos, N_PDS_Comm * comm)
{
  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Packing OptionBlock: " << option_block.getName() << std::endl;

  int size, length;
  int i;

  //----- pack name
  length = option_block.name_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( option_block.name_.c_str(), length, buf, bsize, pos );

  //----- pack netlist location
  int file_number = option_block.netlistLocation_.getFileNumber();
  comm->pack( &file_number, 1, buf, bsize, pos );
  int line_number = option_block.netlistLocation_.getLineNumber();
  comm->pack( &line_number, 1, buf, bsize, pos );

  //----- pack params
  size = option_block.params_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  Util::ParamList::const_iterator it_tpL;
  for (i = 0, it_tpL = option_block.params_.begin(); i < size; ++i, ++it_tpL)
    Pack<Util::Param>::pack(*it_tpL, buf, bsize, pos, comm );

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Packed " << pos << " of " << bsize << " bytes for OptionBlock: " << option_block.getName() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : Util::OptionBlock::unpack
// Purpose       : Unpacks OptionBlock from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
template<>
void
Pack<Util::OptionBlock>::unpack(Util::OptionBlock &option_block, char * pB, int bsize,int & pos, N_PDS_Comm * comm)
{
  int size = 0;
  int length = 0;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  option_block.name_ = std::string( (pB+pos), length);
  pos += length;

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Unpacking OptionBlock: " << option_block.getName() << std::endl;

  //----- unpack netlist location
  int file_number = 0;
  comm->unpack( pB, bsize, pos, &file_number, 1 );
  option_block.netlistLocation_.setFileNumber(file_number);
  int line_number = 0;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  option_block.netlistLocation_.setLineNumber(line_number);

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );

  option_block.params_.clear();
  Util::Param param;
  for( i = 0; i < size; ++i )
  {
    Pack<Util::Param>::unpack(param, pB, bsize, pos, comm );
    option_block.params_.push_back( param );
  }

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Unpacked " << pos << " of " << bsize << " bytes for OptionBlock: " << option_block.getName() << std::endl;
}

} // namespace Xyce

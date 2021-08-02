//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <iostream>

#include <N_DEV_DeviceBlock.h>
#include <N_ERH_ErrorMgr.h>
#include <N_PDS_Comm.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : ModelBlock::ModelBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
ModelBlock::ModelBlock(const std::string & name, const std::string &type, int level)
  : name_(name),
    type_(type),
    level_(level),
    netlistLocation_()
{}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::ModelBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock::ModelBlock(const ModelBlock &right)
  : name_(right.name_),
    type_(right.type_),
    level_(right.level_),
    netlistLocation_(right.netlistLocation_),
    params(right.params)
{}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::~ModelBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock::~ModelBlock()
{}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock & ModelBlock::operator=(const ModelBlock &right)
{
  if (this != &right)
  {
    name_   = right.name_;
    type_   = right.type_;
    level_ =  right.level_;

    netlistLocation_ = right.netlistLocation_;

    params = right.params;
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 4/06/00
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream & os, const ModelBlock & mb)
{
  std::vector<Param>::const_iterator it_pL, end_pL;

  os << "Model Block" << std::endl;
  os << "Model:  " << mb.name_ << std::endl;
  os << " type:  " << mb.type_ << std::endl;
  os << " Level: " << mb.level_ << std::endl;
  os << " netlistLocation_: " << mb.netlistLocation_ << std::endl;
  os << " Tagged Params" << std::endl;
  os << " -------------" << std::endl;

  it_pL=mb.params.begin();
  end_pL=mb.params.end();
  for ( ; it_pL != end_pL; ++it_pL)
  {
    os << it_pL->tag() << "\t" << it_pL->stringValue() << std::endl;
  }

  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::clear
// Purpose       : empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
void ModelBlock::clear()
{
  name_ = "";
  type_ = "";
  level_ = 0;
  params.clear();

  netlistLocation_ = NetlistLocation();
}
//-----------------------------------------------------------------------------
// Function      : InstanceBlock::InstanceBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
InstanceBlock::InstanceBlock (const std::string &name)
  : name_(name),
    modelName_(),
    netlistLocation_(),
    iNumNodes(0),
    numIntVars(0),
    numExtVars(0),
    numStateVars(0),
    modelFlag(0),
    bsourceFlag(0),
    offFlag(0),
    off(0)
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::InstanceBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock::InstanceBlock (const InstanceBlock &right)
  : name_      (right.name_),
    modelName_(right.modelName_),
    netlistLocation_(right.netlistLocation_),
    params    (right.params),
    iNumNodes (right.iNumNodes),
    numIntVars(right.numIntVars),
    numExtVars(right.numExtVars),
    numStateVars(right.numStateVars),
    modelFlag (right.modelFlag),
    bsourceFlag(right.bsourceFlag),
    offFlag   (right.offFlag),
    off       (right.off)
{}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::~InstanceBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock::~InstanceBlock ()
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock & InstanceBlock::operator=(const InstanceBlock &right)
{
  if (this != &right)
  {
    name_      = right.name_;
    modelName_ = right.modelName_;
    iNumNodes = right.iNumNodes;
    numIntVars= right.numIntVars;
    numExtVars= right.numExtVars;
    numStateVars= right.numStateVars;
    modelFlag = right.modelFlag;
    bsourceFlag= right.bsourceFlag;
    offFlag   = right.offFlag;
    off       = right.off;
    netlistLocation_  = right.netlistLocation_;
    params    = right.params;
  }

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::clear
// Purpose       : empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
void InstanceBlock::clear ()
{
  name_ = InstanceName();
  modelName_  = "";
  iNumNodes  = 0;
  numIntVars = 0;
  numExtVars = 0;
  numStateVars = 0;
  modelFlag  = 0;
  bsourceFlag = 0;
  offFlag    = 0;
  off        = 0;
  netlistLocation_ = NetlistLocation();

  params.clear();
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 4/06/00
//-----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream & os, const InstanceBlock & ib)
{
  std::vector<Param>::const_iterator it_tpL, end_tpL;

  os << "Instance Block" << std::endl;
  os << "Name:    " << ib.name_ << std::endl;
  os << " Model:  " << ib.getModelName() << std::endl;
  os << " # Nodes: " << ib.iNumNodes << std::endl;
  os << " # Int Vars: " << ib.numIntVars << std::endl;
  os << " # Ext Vars: " << ib.numExtVars << std::endl;
  os << " # State Vars: " << ib.numStateVars << std::endl;
  os << " modelFlag: " << ib.modelFlag << std::endl;
  os << " bsourceFlag: " << ib.bsourceFlag << std::endl;
  os << " offFlag: " << ib.offFlag << std::endl;
  os << " off: " << ib.off << std::endl;
  os << " netlistLocation_: " << ib.netlistLocation_ << std::endl;
  os << " Tagged Params" << std::endl;
  os << " -------------" << std::endl;

  it_tpL=ib.params.begin();
  end_tpL=ib.params.end();
  for ( ; it_tpL != end_tpL; ++it_tpL)
  {
    os << it_tpL->tag() << "\t" << it_tpL->stringValue() << std::endl;
  }

  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

} // namespace Device

//-----------------------------------------------------------------------------
// Function      : ModelBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
template<>
int
Pack<Device::ModelBlock>::packedByteCount(const Device::ModelBlock &model_block) 
{

  int byteCount = 0;

  int size, length, i;

  //----- count name
  length = model_block.name_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count type
  length = model_block.type_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count level
  byteCount += sizeof(int);

  //----- count params
  size = model_block.params.size();
  byteCount += sizeof(int);
  std::vector<Device::Param>::const_iterator it_tpL = model_block.params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    byteCount += Pack<Device::Param>::packedByteCount(*it_tpL);
  }

  //----- count netlistFilename_
  byteCount += 2*sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::pack
// Purpose       : Packs ModelBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::ModelBlock>::pack(const Device::ModelBlock &model_block, char * buf, int bsize, int & pos, Parallel::Communicator * comm)
{

  int size, length;
  int i;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  //----- pack name
  length = model_block.name_.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( model_block.name_.c_str(), length, buf, bsize, pos );

  //----- pack type
  length = model_block.type_.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( model_block.type_.c_str(), length, buf, bsize, pos );

  //----- pack level
  comm->pack(&model_block.level_, 1, buf, bsize, pos );

  //----- pack params
  size = model_block.params.size();
  comm->pack(&size, 1, buf, bsize, pos );
  std::vector<Device::Param>::const_iterator it_tpL = model_block.params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    Pack<Device::Param>::pack(*it_tpL, buf, bsize, pos, comm );
  }

  //----- pack netlistFilename_
  int file_number = model_block.netlistLocation_.getFileNumber();
  comm->pack(&file_number, 1, buf, bsize, pos );
  int line_number = model_block.netlistLocation_.getLineNumber();
  comm->pack(&line_number, 1, buf, bsize, pos );

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Packed " << pos << " bytes for ModelBlock: " << model_block.name_ << std::endl;

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "ModelBlock::pack") << "Predicted pos does not match actual pos";
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::unpack
// Purpose       : Unpacks ModelBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::ModelBlock>::unpack(Device::ModelBlock &model_block, char * pB, int bsize,int & pos, Parallel::Communicator * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  model_block.name_ = std::string( (pB+pos), length);
  pos += length;

  //----- unpack type
  comm->unpack( pB, bsize, pos, &length, 1 );
  model_block.type_ = std::string( (pB+pos), length);
  pos += length;

  //----- unpack level
  comm->unpack( pB, bsize, pos, &model_block.level_, 1 );

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );
  model_block.params.clear();
  Device::Param dp;
  for( i = 0; i < size; ++i )
  {
    Pack<Device::Param>::unpack(dp, pB, bsize, pos, comm );
    model_block.params.push_back( dp );
  }

  //----- unpack netlistFilename_
  int file_number = 0;
  comm->unpack( pB, bsize, pos, &file_number, 1 );
  model_block.netlistLocation_.setFileNumber(file_number);
  int line_number = 0;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  model_block.netlistLocation_.setLineNumber(line_number);

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Unpacked " << pos << " bytes for ModelBlock: " << model_block.name_ << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::packedByteCount
// Purpose       : count bytes needed to pack block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
template<>
int
Pack<Device::InstanceBlock>::packedByteCount(const Device::InstanceBlock &instance_block) 
{

  int byteCount = 0;

  int size, length, i;

  //----- count name
  const std::string &name = instance_block.getInstanceName().getEncodedName();
  length = name.length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count getModelName()
  length = instance_block.getModelName().length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count params
  size = instance_block.params.size();
  byteCount += sizeof(int);
  
  std::vector<Device::Param>::const_iterator it_tpL = instance_block.params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    byteCount += Pack<Device::Param>::packedByteCount(*it_tpL);
  }

  //----- count iNumNodes
  byteCount += sizeof(int);

  //----- count numIntVars
  byteCount += sizeof(int);

  //----- countnumExtVars
  byteCount += sizeof(int);

  //----- count numStateVars
  byteCount += sizeof(int);

  //----- count modelFlag
  byteCount += sizeof(int);

  //----- count bsourceFlag
  byteCount += sizeof(int);

  //----- count offFlag
  byteCount += sizeof(int);

  //----- pack off
  byteCount += sizeof(int);

  //----- count netlistFilename_ (file number, line number)
  byteCount += 2*sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::pack
// Purpose       : Packs InstanceBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::InstanceBlock>::pack(const Device::InstanceBlock &instance_block, char * buf, int bsize, int & pos, Parallel::Communicator * comm) 
{

  int size, length;
  int i;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  //----- pack name
  const std::string &name = instance_block.getInstanceName().getEncodedName();
  length = name.length();
  comm->pack(&length, 1, buf, bsize, pos);
  comm->pack(name.c_str(), length, buf, bsize, pos);

  //----- pack getModelName()
  length = instance_block.getModelName().length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( instance_block.getModelName().c_str(), length, buf, bsize, pos );

  //----- pack params
  size = instance_block.params.size();
  comm->pack(&size, 1, buf, bsize, pos );
  std::vector<Device::Param>::const_iterator it_tpL = instance_block.params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    Pack<Device::Param>::pack(*it_tpL, buf, bsize, pos, comm );
  }

  //----- pack iNumNodes
  comm->pack(&instance_block.iNumNodes, 1, buf, bsize, pos );

  //----- pack numIntVars
  comm->pack(&instance_block.numIntVars, 1, buf, bsize, pos );

  //----- pack numExtVars
  comm->pack(&instance_block.numExtVars, 1, buf, bsize, pos );

  //----- pack numStateVars
  comm->pack(&instance_block.numStateVars, 1, buf, bsize, pos );

  //----- pack modelFlag
  i = instance_block.modelFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack bsourceFlag
  i = instance_block.bsourceFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack offFlag
  i = instance_block.offFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack off
  i = instance_block.off;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack fileNumber_
  int file_number = instance_block.netlistLocation_.getFileNumber();
  comm->pack(&file_number, 1, buf, bsize, pos );

  //----- pack lineNumber_
  int line_number = instance_block.netlistLocation_.getLineNumber();
  comm->pack(&line_number, 1, buf, bsize, pos );

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Packed " << pos << " bytes for InstanceBlock: " << instance_block.getInstanceName() << std::endl;

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    DevelFatal(*this, "InstanceBlock::pack") << "Predicted pos does not match actual pos";
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::unpack
// Purpose       : Unpacks InstanceBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
template<>
void
Pack<Device::InstanceBlock>::unpack(Device::InstanceBlock &instance_block, char * pB, int bsize, int & pos, Parallel::Communicator * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  instance_block.name_ = Device::InstanceName(std::string( (pB+pos), length));
  pos += length;

  //----- unpack getModelName()
  comm->unpack( pB, bsize, pos, &length, 1 );
  instance_block.modelName_ = std::string( (pB+pos), length);
  pos += length;

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );
  instance_block.params.clear();
  Device::Param dp;
  for( i = 0; i < size; ++i )
  {
    Pack<Device::Param>::unpack(dp, pB, bsize, pos, comm );
    instance_block.params.push_back( dp );
  }

  //----- unpack iNumNodes
  comm->unpack( pB, bsize, pos, &instance_block.iNumNodes, 1 );

  //----- unpack numIntVars
  comm->unpack( pB, bsize, pos, &instance_block.numIntVars, 1 );

  //----- unpack numExtVars
  comm->unpack( pB, bsize, pos, &instance_block.numExtVars, 1 );

  //----- unpack numStateVars
  comm->unpack( pB, bsize, pos, &instance_block.numStateVars, 1 );

  //----- unpack modelFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  instance_block.modelFlag = ( i != 0 );

  //----- unpack bsourceFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  instance_block.bsourceFlag = ( i != 0 );

  //----- unpack offFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  instance_block.offFlag = ( i != 0 );

  //----- unpack off
  comm->unpack( pB, bsize, pos, &i, 1 );
  instance_block.off = ( i != 0 );

  //----- unpack netlistFilename_
  int file_number = 0;
  comm->unpack( pB, bsize, pos, &file_number, 1 );
  instance_block.netlistLocation_.setFileNumber(file_number);

  //----- unpack lineNumber_
  int line_number = 0;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  instance_block.netlistLocation_.setLineNumber(line_number);

  if (DEBUG_TOPOLOGY)
    Xyce::dout() << "Unpacked " << pos << " bytes for InstanceBlock: " << instance_block.getInstanceName() << std::endl;
}

} // namespace Xyce

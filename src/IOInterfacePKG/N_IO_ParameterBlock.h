//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : Declare the N_IO_ParameterBlock class instantiations of
//                  which are associated with netlist .MODEL lines.
//
// Special Notes  : ERK.  It seems that the name "N_IO_ModelBlock" would have been
//                  more appropriate and less confusing, or possibly
//                  N_IO_ModelParametersBlock.  Calling it "parameters block"
//                  makes it sound much more general than it really is.
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/10/2001
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef N_IO_PARAMETERBLOCK_H
#define N_IO_PARAMETERBLOCK_H

// ---------- Standard Includes ----------

#include <string>

#include <vector>


// ----------   Xyce Includes   ----------

#include <N_IO_fwd.h>
#include <N_IO_SpiceSeparatedFieldTool.h>

#include <N_DEV_Param.h>

#include <N_UTL_Pack.h>

#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : ParameterBlock
// Purpose       :
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 09/10/2001
//-----------------------------------------------------------------------------

class ParameterBlock
{
  friend class Pack<ParameterBlock>;

public:
  // Constructor.
  ParameterBlock()
  {}

  // Constructor.
  ParameterBlock(
     std::string const& fileName,
     TokenVector
     const& parsedInputLine);

  // Copy Constructor.
  ParameterBlock(ParameterBlock const& rhsPB)
    : modelData(rhsPB.modelData),
      expressionValuedParams_(rhsPB.expressionValuedParams_)
  {}

  // Destructor.
  virtual ~ParameterBlock()
  {}

  // Public methods.

  // Prints the details of an ParameterBlock to standard out.
  void print();

  // Add the default parameters for a model from metadata.
  void addDefaultModelParameters( CircuitMetadata & metadata );

  bool hasExpressionValuedParams(){return !expressionValuedParams_.empty();};

  void setParameterValues(CircuitContext* contextPtr);
  // Check each parameter to see if it is an expression, if
  // so, evaluate the parameter.

  void setName( std::string const& name );
  void setType( std::string const& type );
  void setLevel( std::string const& level );
  void addParameter( Device::Param const& parameter );
  void addParameters( std::vector<Device::Param> const& parameters );
  void setParameter( int const& i, Device::Param const& parameter );

  NetlistLocation netlistLocation() const;
    
  const std::string& getName() const;
  const std::string& getType() const;
  int getLevel() const;
  std::vector<Device::Param> & getParams();

  Device::Param* findParameter( Device::Param const& parameter );
  int getNumberOfParameters() const;
  Device::Param getParameter( int const& i ) const;

  ParameterBlock & operator=(ParameterBlock const& rhsPB)
  {
    expressionValuedParams_ = rhsPB.expressionValuedParams_;
    modelData = rhsPB.modelData;
    return *this;
  };

  bool operator==(std::string const& name);
  bool operator!=(std::string const& name);

public:
  Device::ModelBlock modelData;

  std::map< std::string, std::vector<std::vector<Device::Param> > > inputCompositeParamVecMap;

private:
  bool extractModelData(TokenVector const& parsedInputLine);

  bool defaultApplied_;
  std::vector<Device::Param> expressionValuedParams_;

  void addDefaultCompositeModelParameters_(
    CircuitMetadata &                   metadata,
    const Device::Param &               baseParam ,
    std::map<std::string,bool> &        paramMetadataExistMap);

  // used to prevent segfaults during processing of parsedLine vector
  bool validLinePosition_(
    int              pos,
    int              lineLength,
    ExtendedString   paramName);
};

inline bool ParameterBlock::operator==(std::string const& rhsName)
{
  return (getName() == rhsName);
}

inline bool ParameterBlock::operator!=(std::string const& rhsName)
{
  return (getName() != rhsName);
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setName( std::string const& nameIn )
{
  modelData.setName(nameIn);
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setType( std::string const& typeIn )
{
  modelData.setType(typeIn);
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setLevel( std::string const& levelIn )
{
  Device::Param levelParam( "LEVEL", levelIn);
  modelData.setLevel(levelParam.getImmutableValue<int>());
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline const std::string& ParameterBlock::getName() const
{
  return modelData.getName();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline const std::string& ParameterBlock::getType() const
{
  return modelData.getType();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline int ParameterBlock::getLevel() const
{
  Device::Param levelParam( "LEVEL", modelData.getLevel() );
  return levelParam.getImmutableValue<int>();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/26/05
//-----------------------------------------------------------------------------
inline std::vector<Device::Param> & ParameterBlock::getParams()
{
  return modelData.params;
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addParameter
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::addParameter( Device::Param const& parameter )
{
  modelData.params.push_back( parameter );
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::addParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::addParameters(
   std::vector<Device::Param> const& parametersIn )
{
  modelData.params.insert( modelData.params.end(),
                           parametersIn.begin(), parametersIn.end() );
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::setParamter
// Purpose       :
// Special Notes : It is assumed that getNumberOfParemeters was called to
//                 ensure that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline void ParameterBlock::setParameter( int const& i,
                                          Device::Param const& parameter )
{
  modelData.params[i].setVal(parameter.stringValue());
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getNumberOfParameters
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline int ParameterBlock::getNumberOfParameters() const
{
  return modelData.params.size();
}

//-----------------------------------------------------------------------------
// Function      : ParameterBlock::getParameter
// Purpose       :
// Special Notes : It is assumed getNumberOfParameters was called to ensure
//                 that i is in range.
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 01/08/2001
//-----------------------------------------------------------------------------
inline Device::Param ParameterBlock::getParameter( int const& i ) const
{
  return modelData.params[i];
}

} // namespace IO
} // namespace Xyce

#endif

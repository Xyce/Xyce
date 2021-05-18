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
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef deviceExpressionGroup_H
#define deviceExpressionGroup_H
#include <Xyce_config.h>

#include<string>
#include<complex>
#include<unordered_map>

#include <N_UTL_fwd.h>
#include <N_DEV_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_IO_fwd.h>
#include <N_TIA_fwd.h>
#include <N_ANP_fwd.h>

#include <N_UTL_Op.h>

#include<newExpression.h>
#include <ExpressionType.h>
#include <expressionGroup.h>
#include <mainXyceExpressionGroup.h>
#include <deviceExpressionGroup.h>
#include <N_UTL_ExtendedString.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace Util {

#define CONSTCtoK    (273.15)  


//-----------------------------------------------------------------------------
// Class         : deviceExpressionGroup
//
// Purpose       : This group class is for expressions owned by devices 
//                   (usually via parameters)
//
// Special Notes : Originally, devices used the "mainXyceExprssionGroup" class,
//                 which was originally designed to be a single monolithic class
//                 that handled all data retrieval nececssary for Xyce 
//                 expressions.
//
//                 Using that class in the device package proved to be 
//                 problematic for several reasons.  The biggest problem
//                 was that the mainXyceExpressionGroup relied on a lot of 
//                 "AllReduce" operations in parallel, and this cannot 
//                 work for devices, since in Xyce's parallel design, 
//                 devices are only owned by a single processor.  So, 
//                 especially for solution variables, this needed to be 
//                 different for device expressions.
//
//                 This class is derived from the mainXyceExpressionGroup, 
//                 and it takes an RCP of that group as an argument to its 
//                 constructor.  That way it can re-use most of its structure, 
//                 and just replace the stuff that didn't work for devices.
//
// Creator       : Eric Keiter
// Creation Date : 8/22/2020
//-----------------------------------------------------------------------------
class deviceExpressionGroup : public mainXyceExpressionGroup
{
friend class outputsXyceExpressionGroup;
friend class ExpressionData;

public:
  deviceExpressionGroup ( const Teuchos::RCP<Xyce::Util::mainXyceExpressionGroup> & mainGroup);
  ~deviceExpressionGroup ();

  void setSolutionLIDs( 
    const std::vector<std::string> & expVarNames, 
    const std::vector<int> & expVarLIDs, int lo, int hi);

  virtual bool getSolutionVal(const std::string & nodeName, double & retval );
  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval );

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, double & retval ) 
  { return getSolutionVal(deviceName,retval); }

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, std::complex<double> & retval )
  { return getSolutionVal(deviceName,retval); }

  virtual bool getGlobalParameterVal (const std::string & paramName, double & retval );
  virtual bool getGlobalParameterVal (const std::string & paramName, std::complex<double> & retval );

private:
  std::unordered_map<std::string,int> lidMap_;
};

}
}

#endif


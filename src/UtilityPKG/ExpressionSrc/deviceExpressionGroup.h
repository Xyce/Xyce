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
// Purpose       : This is the "main" group class for connecting the new 
//                 expression library to Xyce
//
// Special Notes : Unlike the "xyceExpressionGroup" this class is not 
//                 intended to be lightweight.   My intention is for there to 
//                 be a single instance of this (or possible a small 
//                 number of copies).
//
//                 This class will contain all the machinery necessary to provide 
//                 newExpression with the information it needs to evaluate *any* 
//                 expression for Xyce.  ie, any external information, such as 
//                 solution values, global parameter values, etc.
//
//                 The most significant difference is (intended to be) the handling 
//                 of user defined functions.  i.e., .funcs.
//
//                 The old expression library handled these in a very inefficient 
//                 way, via string substitutions.
//
//                 The new expression library handles them by attaching the node 
//                 of the .func to the expression that is calling it.
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


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

#ifndef outputsXyceExpressionGroup_H
#define outputsXyceExpressionGroup_H
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

#include <newExpression.h>
#include <ExpressionType.h>
#include <expressionGroup.h>
#include <N_UTL_ExtendedString.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace Util {

#define CONSTCtoK    (273.15)  

//-----------------------------------------------------------------------------
// Class         : outputsXyceExpressionGroup
//
// Purpose       : This is the "outputs" group class for connecting the new 
//                 expression library to Xyce
//
// Special Notes : This group class is intended to be used with expressions 
//                 on the .PRINT line and similar outputs.
//
//                 As such, it has some capabilities that were originally 
//                 in the ExpressionData class and/or the various "Op" classes.
//
//                 .PRINT line outputs have the following requirements:
//                 - they don't need derivatives.
//
//                 - they aren't needed during initial setup.  So, certain 
//                    convenient data structures that come from the output 
//                    manager can actually be epxected to work.
//
//                 - they output lots of things that could (probably) never 
//                    be used in a device parameter or Bsrc.  
//
//                 - UQ operators like Gauss and Agauss probably 
//                   aren't needed.  Those operators are typically
//                   applied to inputs, not ouputs.  As such, this probably 
//                   doesn't need to be a singleton.  Instead can be 
//                   a group allocated by each ExpressionData object, etc.
//                   That will make some of the issues for it a bit simpler,
//                   probably.
//
//
// Creator       : Eric Keiter
// Creation Date : 4/29/2020
//-----------------------------------------------------------------------------
class outputsXyceExpressionGroup : public baseExpressionGroup
{
friend class mainXyceExpressionGroup;

public:

  outputsXyceExpressionGroup ( 
      Parallel::Communicator & comm, Topo::Topology & top,
      Analysis::AnalysisManager &analysis_manager,
      Device::DeviceMgr & device_manager,
      IO::OutputMgr &output_manager
      ) ;

  ~outputsXyceExpressionGroup ();

  virtual bool getSolutionVal(const std::string & nodeName, double & retval );
  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval );

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, double & retval );
  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, std::complex<double> & retval );

  virtual bool getParameterVal (const std::string & paramName, double & retval );
  virtual bool getParameterVal (const std::string & paramName, std::complex<double> & retval );

  virtual bool getInternalDeviceVar (const std::string & deviceName, double & retval );
  virtual bool getInternalDeviceVar (const std::string & deviceName, std::complex<double> & retval );

  virtual bool getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, double & retval);
  virtual bool getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & retval);

  virtual bool getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, double & retval);
  virtual bool getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & retval);

  virtual bool getONoise(double & retval);
  virtual bool getONoise(std::complex<double> & retval);

  virtual bool getINoise(double & retval);
  virtual bool getINoise(std::complex<double> & retval);

  virtual bool getPower(const std::string & tag, const std::string & deviceName, double & retval);
  virtual bool getPower(const std::string & tag, const std::string & deviceName, std::complex<double> & retval);

  virtual bool getSparam (const std::vector<int> & args, double & retval );
  virtual bool getSparam (const std::vector<int> & args, std::complex<double> & retval );

  virtual bool getYparam (const std::vector<int> & args, double & retval );
  virtual bool getYparam (const std::vector<int> & args, std::complex<double> & retval );

  virtual bool getZparam (const std::vector<int> & args, double & retval );
  virtual bool getZparam (const std::vector<int> & args, std::complex<double> & retval );

  virtual double getTimeStep ();
  virtual double getTimeStepAlpha () { return alpha_; }
  virtual double getTimeStepPrefac () { return (getTimeStepAlpha() / getTimeStep ()) ; } // FIX

  virtual double getTime();
  virtual double getTemp();
  virtual double getVT  ();
  virtual double getFreq();
  virtual double getGmin();

  virtual double getBpTol();
  virtual double getStartingTimeStep();
  virtual double getFinalTime();

  virtual unsigned int getStepNumber ();

  virtual bool getPhaseOutputUsesRadians();

  virtual void setRFParamsRequested(std::string type);

  void setAliasNodeMap( const IO::AliasNodeMap & anm ) { aliasNodeMap_ = anm; }

  void setOpData (const Util::Op::OpData & opData) { opData_ = opData; }


private:
  Parallel::Communicator & comm_;
  Topo::Topology & top_;

  Analysis::AnalysisManager & analysisManager_;
  Device::DeviceMgr & deviceManager_;

  IO::AliasNodeMap aliasNodeMap_;

  IO::OutputMgr &outputManager_;

  double time_, temp_, VT_, freq_, gmin_;
  double dt_, alpha_;

  Util::Op::OpData opData_;
};

}
}

#endif


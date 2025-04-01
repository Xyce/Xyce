//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
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
// Purpose        : Base expression group class.  
//
// Special Notes  : This is used during netlist parsing.  
//                  During parsing, the rest of the problem isn't set up yet, 
//                  so there isn't anything meaningful for a group to do.  
//                  As such, this group is nearly empty.
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//-----------------------------------------------------------------------------

#ifndef expressionGroup_H
#define expressionGroup_H

#include<vector>
#include<string>
#include<complex>

#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Util {

class newExpression;

#define CONSTCtoK    (273.15)  

//-----------------------------------------------------------------------------
// Class         : baseExpressionGroup
// Purpose       : 
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 10/28/2019
//-----------------------------------------------------------------------------
class baseExpressionGroup
{
public:

  baseExpressionGroup () {};
  virtual ~baseExpressionGroup () {};

  // do everything
  virtual bool setupGroup(newExpression &expr) {return true;};
  virtual bool putValues(newExpression & expr);

  // these do not set a value at all, they just pass thru and return true.
  virtual bool getSolutionVal(const std::string & nodeName, double & retval ) { return true; }
  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval ) { return true; }
  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, double & retval ) { return true; }
  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, std::complex<double> & retval ) { return true; }
  virtual bool getParameterVal (const std::string & paramName, double & retval ) {return true;}
  virtual bool getParameterVal (const std::string & paramName, std::complex<double> & retval ) {return true;}

  // these do set a value, and return false.  Rethink?  
  virtual bool getInternalDeviceVar (const std::string & deviceName, double & retval ) { retval=0.0; return false; }
  virtual bool getInternalDeviceVar (const std::string & deviceName, std::complex<double> & retval ) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, double & retval) { retval=0.0; return false; }
  virtual bool getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, double & retval) { retval=0.0; return false; }
  virtual bool getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getONoise(double & retval) { retval=0.0; return false; }
  virtual bool getONoise(std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getINoise(double & retval) { retval=0.0; return false; }
  virtual bool getINoise(std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getPower(const std::string & tag, const std::string & deviceName, double & retval) { retval=0.0; return false; }
  virtual bool getPower(const std::string & tag, const std::string & deviceName, std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getSparam (const std::vector<int> & args, double & retval ) { retval=0.0; return false; }
  virtual bool getSparam (const std::vector<int> & args, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getYparam (const std::vector<int> & args, double & retval ) { retval=0.0; return false; }
  virtual bool getYparam (const std::vector<int> & args, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getZparam (const std::vector<int> & args, double & retval ) { retval=0.0; return false; }
  virtual bool getZparam (const std::vector<int> & args, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }

  virtual double getTimeStep () { return 0.0; }

  virtual double getTimeStepAlpha () { return 0.0; }
  virtual double getTimeStepPrefac () { return 0.0; }

  virtual double getTime() { return 0.0; }
  virtual double getTemp() { return 0.0; }
  virtual double getVT() { return 0.0; }
  virtual double getFreq() { return 0.0; }
  virtual double getGmin() { return 0.0; }

  virtual double getBpTol() { return 0.0; }
  virtual double getStartingTimeStep() { return 0.0; }
  virtual double getFinalTime() { return 0.0; }

  virtual unsigned int getStepNumber () { return 0; }

  virtual bool getPhaseOutputUsesRadians() { return true; }

  virtual void setRFParamsRequested(std::string type) {}

private:

};

}
}

#endif


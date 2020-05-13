

#ifndef expressionGroup_H
#define expressionGroup_H

#include<string>
#include<complex>

#include <Teuchos_RCP.hpp>

namespace Xyce {
namespace Util {

class newExpression;

#define CONSTCtoK    (273.15)  

//-----------------------------------------------------------------------------
// Class         : baseExpressionGroup
// Purpose       : this is purely virtual class
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 10/28/2019
//-----------------------------------------------------------------------------
class baseExpressionGroup
{
public:

  baseExpressionGroup () {};
  virtual ~baseExpressionGroup () {};

  // ERK NOTE:  Need to have a "notify" for .STEP loops.  
  // Important for time-dependent expressions.

  // might adopt "next" "curr" "prev"? or some other (better) nomenclature?
  virtual bool getSolutionDdt (const std::string & nodeName, double & retval ) { retval=0.0; return false; }
  virtual bool getSolutionSdt (const std::string & nodeName, double & retval ) { retval=0.0; return false; }
  virtual bool getSolutionVal (const std::string & nodeName, double & retval ) { retval=0.0; return false; }

  virtual bool getCurrentVal  ( const std::string & deviceName, const std::string & designator, double & retval ) { retval=0.0; return false; }

  virtual bool getInternalDeviceVar (const std::string & deviceName, double & retval ) { retval=0.0; return false; }
  virtual bool getInternalDeviceVar (const std::string & deviceName, std::complex<double> & retval ) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getDnoNoiseDeviceVar(const std::string & deviceName, double & retval) { retval=0.0; return false; }
  virtual bool getDnoNoiseDeviceVar(const std::string & deviceName, std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getDniNoiseDeviceVar(const std::string & deviceName, double & retval) { retval=0.0; return false; }
  virtual bool getDniNoiseDeviceVar(const std::string & deviceName, std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getONoise(double & retval) { retval=0.0; return false; }
  virtual bool getONoise(std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getINoise(double & retval) { retval=0.0; return false; }
  virtual bool getINoise(std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getPower(const std::string & deviceName, double & retval) { retval=0.0; return false; }
  virtual bool getPower(const std::string & deviceName, std::complex<double> & retval) {retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getSolutionDdt (const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }
  virtual bool getSolutionSdt (const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }
  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }

  virtual bool getGlobalParameterVal (const std::string & nodeName, double & retval ) { retval=0.0; return false; } 
  virtual bool getGlobalParameterVal (const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; } 

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

private:

};

}
}

#endif


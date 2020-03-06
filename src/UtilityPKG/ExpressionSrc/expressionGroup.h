

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

  virtual bool resolveExpression (Xyce::Util::newExpression & exp){return true;};

  virtual bool isOption (const std::string & optionStr) { return false; }

  // ERK NOTE:  Need to have a "notify" for .STEP loops.  
  // Important for time-dependent expressions.

  // might adopt "next" "curr" "prev"? or some other (better) nomenclature?
  virtual bool getSolutionDdt (const std::string & nodeName, double & retval ) { retval=0.0; return false; }
  virtual bool getSolutionSdt (const std::string & nodeName, double & retval ) { retval=0.0; return false; }
  virtual bool getSolutionVal (const std::string & nodeName, double & retval ) { retval=0.0; return false; }

  // should this become obsolete eventually? 
  virtual bool getGlobalParameterVal (const std::string & nodeName, double & retval ) { retval=0.0; return false; } 

  virtual bool getSolutionDdt (const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }
  virtual bool getSolutionSdt (const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }
  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; }

  // should this become obsolete eventually? 
  virtual bool getGlobalParameterVal (const std::string & nodeName, std::complex<double> & retval ) { retval=std::complex<double>(0.0,0.0); return false; } 

  virtual double getTimeStep () { return 0.0; }

  virtual double getTimeStepAlpha () { return 0.0; }
  virtual double getTimeStepPrefac () { return 0.0; }

  virtual double getTime() { return 0.0; }
  virtual double getTemp() { return 0.0; }
  virtual double getVT() { return 0.0; }
  virtual double getFreq() { return 0.0; }

  virtual double getBpTol() { return 0.0; }

  virtual bool getFunction    (const std::string & name, Xyce::Util::newExpression & exp) {return false;};
  virtual bool getParam       (const std::string & name, Xyce::Util::newExpression & exp) {return false;};
  virtual bool getGlobalParam (const std::string & name, Xyce::Util::newExpression & exp) {return false;};

private:

};

}
}

#endif


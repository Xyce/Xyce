

#ifndef xyceExpressionGroup_H
#define xyceExpressionGroup_H

#include<string>
#include<complex>

#include <ExpressionType.h>
#include <expressionGroup.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Util {

#define CONSTCtoK    (273.15)  

//-----------------------------------------------------------------------------
// Class         : xyceExpressionGroup
//
// Purpose       : This is a prototype, developed primarily for this test program.
//
// Special Notes : When migrated to Xyce, a new "Xyce" group will need to be 
//                 developed.
//
// Creator       : Eric Keiter
// Creation Date : 2/12/2020
//-----------------------------------------------------------------------------
class xyceExpressionGroup : public baseExpressionGroup
{
public:

  xyceExpressionGroup ();
  ~xyceExpressionGroup () {};
  virtual bool resolveExpression (Xyce::Util::newExpression & exp);

  virtual bool isOption (const std::string & optionStr)
  {
    std::string tmp = optionStr;
    Xyce::Util::toLower(tmp);
    return false; //  FIX THIS
  }

  virtual bool getSolutionSdt(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toLower(tmp);
    retval = 0.0;
    return success; // FIX THIS
  }

  virtual bool getSolutionDdt(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toLower(tmp);
    retval = 0.0;
    return success; // FIX THIS
  }

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toLower(tmp);
    retval = 0.0;
    return success; // FIX THIS
  }

  // ERK NOTE:  Need to have a "notify" (or something) for .STEP loops.  
  // Important for time-dependent expressions.
  //
  // Also, for "time" related quantities, should these have names 
  // like "getCurrTimeStep", rather than "getTimeStep"?  Do we ever need much else?
  virtual double getTimeStep () { return 1.0e-8; } // WAG
  virtual double getTimeStepAlpha () { return 1.0; }
  virtual double getTimeStepPrefac () { return (getTimeStepAlpha() / getTimeStep ()) ; }

  virtual double getTime() { return 0.0;} // nd_.getTime(); }
  virtual double getTemp() { return 0.0;} // nd_.getTemp() - CONSTCtoK; }
  virtual double getVT  () { return 0.0;} // nd_.getTemp() * CONSTCtoK; }
  virtual double getFreq() { return 0.0;} // nd_.getFreq(); }

  // in a real Xyce group, need something like this:
  //solver_state.bpTol_ = analysis_manager.getStepErrorControl().getBreakPointLess().tolerance_;
  virtual double getBpTol() { return 0.0; }

#if 0
  virtual bool getFunction    (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression > & exp);
  virtual bool getParam       (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression > & exp);
  virtual bool getGlobalParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression > & exp);
#else
  virtual bool getFunction    (const std::string & name, Xyce::Util::newExpression & exp);
  virtual bool getParam       (const std::string & name, Xyce::Util::newExpression & exp);
  virtual bool getGlobalParam (const std::string & name, Xyce::Util::newExpression & exp);
#endif

private:

};

}
}

#endif


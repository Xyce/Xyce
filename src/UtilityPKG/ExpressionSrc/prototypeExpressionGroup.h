

#ifndef prototypeExpressionGroup_H
#define prototypeExpressionGroup_H

#include<string>
#include<complex>
//#include <Sacado_No_Kokkos.hpp>

#include <ExpressionType.h>
#include <expressionGroup.h>
#include <netlistData.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Util {

#define CONSTCtoK    (273.15)  

//-----------------------------------------------------------------------------
// Class         : prototypeExpressionGroup
//
// Purpose       : This is a prototype, developed primarily for this test program.
//
// Special Notes : When migrated to Xyce, a new "Xyce" group will need to be 
//                 developed.
//
// Creator       : Eric Keiter
// Creation Date : 10/28/2019
//-----------------------------------------------------------------------------
class prototypeExpressionGroup : public baseExpressionGroup
{
public:

  prototypeExpressionGroup (Xyce::Util::netlistData & nd);

  ~prototypeExpressionGroup () {};

  virtual bool resolveExpression (Teuchos::RCP<Xyce::Util::newExpression> exp);

  virtual bool isOption (const std::string & optionStr)
  {
    std::string tmp = optionStr;
    Xyce::Util::toUpper(tmp);
    return (nd_.getOptions().find(tmp) != nd_.getOptions().end());
  }

  virtual bool getSolutionSdt(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toUpper(tmp);

    retval = 0.0;
    std::unordered_map<std::string,Xyce::Util::nodeData> & nodes = nd_.getNodes ();
    if (nodes.find(tmp) != nodes.end())
    {
      int solIndex = nodes[tmp].index_;
      retval = Xyce::Util::Value(nodes[tmp].value_);
    }
    else
    {
      success=false;
    }
    return success;
  }

  virtual bool getSolutionDdt(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toUpper(tmp);

    retval = 0.0;
    std::unordered_map<std::string,Xyce::Util::nodeData> & nodes = nd_.getNodes ();
    if (nodes.find(tmp) != nodes.end())
    {
      int solIndex = nodes[tmp].index_;
      retval = Xyce::Util::Value(nodes[tmp].value_);
    }
    else
    {
      success=false;
    }
    return success;
  }

  virtual bool getSolutionVal(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toUpper(tmp);

    retval = 0.0;
    std::unordered_map<std::string,Xyce::Util::nodeData> & nodes = nd_.getNodes ();
    if (nodes.find(tmp) != nodes.end())
    {
      int solIndex = nodes[tmp].index_;
      retval = Xyce::Util::Value(nodes[tmp].value_);
    }
    else
    {
      success=false;
    }
    return success;
  }

  // ERK NOTE:  Need to have a "notify" (or something) for .STEP loops.  
  // Important for time-dependent expressions.
  //
  // Also, for "time" related quantities, should these have names 
  // like "getCurrTimeStep", rather than "getTimeStep"?  Do we ever need much else?
  virtual double getTimeStep () { return 1.0e-8; } // WAG
  virtual double getTimeStepAlpha () { return 1.0; }
  virtual double getTimeStepPrefac () { return (getTimeStepAlpha() / getTimeStep ()) ; }

  virtual double getTime() { return nd_.getTime(); }
  virtual double getTemp() { return nd_.getTemp() - CONSTCtoK; }
  virtual double getVT  () { return nd_.getTemp() * CONSTCtoK; }
  virtual double getFreq() { return nd_.getFreq(); }

  // in a real Xyce group, need something like this:
  //solver_state.bpTol_ = analysis_manager.getStepErrorControl().getBreakPointLess().tolerance_;
  virtual double getBpTol() { return 0.0; }

#if 0
  virtual bool getFunction    (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression > & exp);
  virtual bool getParam       (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression > & exp);
  virtual bool getGlobalParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression > & exp);
#else
  virtual bool getFunction    (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> exp);
  virtual bool getParam       (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> exp);
  virtual bool getGlobalParam (const std::string & name, Teuchos::RCP<Xyce::Util::newExpression> exp);
#endif

private:
  Xyce::Util::netlistData & nd_;

};

}
}

#endif


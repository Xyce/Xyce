

#ifndef mainXyceExpressionGroup_H
#define mainXyceExpressionGroup_H
#include <Xyce_config.h>

#include<string>
#include<complex>
#include<unordered_map>

#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_IO_fwd.h>
#include <N_TIA_fwd.h>

#include<newExpression.h>
#include <ExpressionType.h>
#include <expressionGroup.h>
#include <N_UTL_ExtendedString.h>

namespace Xyce {
namespace Util {

#define CONSTCtoK    (273.15)  

//-----------------------------------------------------------------------------
// Class         : mainXyceExpressionGroup
//
// Purpose       : This is the "main" group class for connecting the new 
//                 expression library to Xyce
//
// Special Notes : Unlike the "xyceExpressionGroup" this class is not 
//                 intended to be lightweight.   My intention is for there to 
//                 be a single instance of this (or possible a small 
//                 number of copies).
//
//                 This class will contain all the machinery necessary to provide newExpression with the information it needs to evaluate *any* expression for Xyce.  ie, any external information, such as solution values, global parameter values, etc.
//
//                 The most significant difference is (intended to be) the handling of user defined functions.  i.e., .funcs.
//
//                 The old expression library handled these in a very inefficient way, via string substitutions.
//
//                 The new expression library handles them by attaching the node of the .func to the expression that is calling it.
//
// Creator       : Eric Keiter
// Creation Date : 2/12/2020
//-----------------------------------------------------------------------------
class mainXyceExpressionGroup : public baseExpressionGroup
{
public:

  mainXyceExpressionGroup ( 
      N_PDS_Comm & comm, Topo::Topology & top,
      TimeIntg::DataStore & dataStore,
      const IO::AliasNodeMap & aliasNM) :
    comm_(comm),
    top_(top),
    dataStore_(dataStore),
    aliasNodeMap_(aliasNM),
  time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), dt_(0.0), alpha_(0.0)
  {};

  ~mainXyceExpressionGroup () {};

  virtual bool isOption (const std::string & optionStr)
  {
    std::string tmp = optionStr;
    Xyce::Util::toUpper(tmp);
    return false; //  FIX THIS
  }

  virtual bool getSolutionSdt(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toUpper(tmp);
    retval = 0.0;
    return success; // FIX THIS
  }

  virtual bool getSolutionDdt(const std::string & nodeName, double & retval )
  {
    bool success=true;
    std::string tmp = nodeName;
    Xyce::Util::toUpper(tmp);
    retval = 0.0;
    return success; // FIX THIS
  }


  virtual bool getSolutionVal(const std::string & nodeName, double & retval );

  virtual bool setTimeStep (double dt) { dt_ = dt; return true; } // WAG
  virtual bool setTimeStepAlpha (double a) { alpha_ = a; return true; }

  virtual double getTimeStep () { return dt_; } // WAG
  virtual double getTimeStepAlpha () { return alpha_; }
  virtual double getTimeStepPrefac () { return (getTimeStepAlpha() / getTimeStep ()) ; } // FIX

  virtual bool setTime(double t) { time_ = t; return true; } 
  virtual bool setTemp(double t) { temp_ = t; return true; } 
  virtual bool setVT  (double v) { VT_ = v; return true; } 
  virtual bool setFreq(double f) { freq_ = f; return true; } 

  virtual double getTime() { return time_;} 
  virtual double getTemp() { return temp_;} 
  virtual double getVT  () { return VT_;} 
  virtual double getFreq() { return freq_;} 

  // in a real Xyce group, need something like this:
  //solver_state.bpTol_ = analysis_manager.getStepErrorControl().getBreakPointLess().tolerance_;
  virtual double getBpTol() { return 0.0; }

private:

  N_PDS_Comm & comm_;
  Topo::Topology & top_;
  TimeIntg::DataStore & dataStore_;
  const IO::AliasNodeMap & aliasNodeMap_; // = output_manager.getAliasNodeMap().find(objVec[iobj]->expVarNames[i]);

  double time_, temp_, VT_, freq_;
  double dt_, alpha_;
};

}
}

#endif


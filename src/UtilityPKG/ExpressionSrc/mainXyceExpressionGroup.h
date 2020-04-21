

#ifndef mainXyceExpressionGroup_H
#define mainXyceExpressionGroup_H
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
      Analysis::AnalysisManager &analysis_manager,
      Device::DeviceMgr & device_manager
      ) ;

  ~mainXyceExpressionGroup ();

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

  virtual double getTimeStep ();
  virtual double getTimeStepAlpha () { return alpha_; }
  virtual double getTimeStepPrefac () { return (getTimeStepAlpha() / getTimeStep ()) ; } // FIX

  virtual double getTime();
  virtual double getTemp();
  virtual double getVT  ();
  virtual double getFreq();

  // in a real Xyce group, need something like this:
  //solver_state.bpTol_ = analysis_manager.getStepErrorControl().getBreakPointLess().tolerance_;
  virtual double getBpTol() { return 0.0; }

  void setAliasNodeMap( const IO::AliasNodeMap & anm ) { aliasNodeMapPtr_ = &anm; }

private:

  N_PDS_Comm & comm_;
  Topo::Topology & top_;

  Analysis::AnalysisManager & analysisManager_;
  Device::DeviceMgr & deviceManager_;

  const IO::AliasNodeMap * aliasNodeMapPtr_; // = output_manager.getAliasNodeMap().find(objVec[iobj]->expVarNames[i]);

  Util::Op::Operator * tempOp_;

  double time_, temp_, VT_, freq_;
  double dt_, alpha_;
};

}
}

#endif


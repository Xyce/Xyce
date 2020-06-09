

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
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace Util {

#define CONSTCtoK    (273.15)  

int findNodeIndex(
  const std::string &       name,
  const NodeNameMap &       node_map,
  const IO::AliasNodeMap &      alias_map);

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
// Creation Date : 2/12/2020
//-----------------------------------------------------------------------------
class mainXyceExpressionGroup : public baseExpressionGroup
{
friend class outputsXyceExpressionGroup;
friend class ExpressionData;

public:

  mainXyceExpressionGroup ( 
      N_PDS_Comm & comm, Topo::Topology & top,
      Analysis::AnalysisManager &analysis_manager,
      Device::DeviceMgr & device_manager,
      IO::OutputMgr &output_manager
      ) ;

  ~mainXyceExpressionGroup ();

  virtual bool getSolutionDdt (const std::string & nodeName, std::complex<double> & retval ) 
  { 
    retval=std::complex<double>(0.0,0.0);
    return false; 
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
  virtual bool getSolutionVal(const std::string & nodeName, std::complex<double> & retval );

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, double & retval ) 
  { return getSolutionVal(deviceName,retval); }

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, std::complex<double> & retval )
  { return getSolutionVal(deviceName,retval); }

  virtual bool getGlobalParameterVal (const std::string & paramName, double & retval );
  virtual bool getGlobalParameterVal (const std::string & paramName, std::complex<double> & retval );

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

  void setAliasNodeMap( const IO::AliasNodeMap & anm ) { aliasNodeMap_ = anm; }

private:

  int getSolutionGID_(const std::string & nodeName);

  N_PDS_Comm & comm_;
  Topo::Topology & top_;

  Analysis::AnalysisManager & analysisManager_;
  Device::DeviceMgr & deviceManager_;

  IO::AliasNodeMap aliasNodeMap_;

  IO::OutputMgr &outputManager_;

  double time_, temp_, VT_, freq_, gmin_;
  double dt_, alpha_;
};

}
}

#endif


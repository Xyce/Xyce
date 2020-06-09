

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
      N_PDS_Comm & comm, Topo::Topology & top,
      Analysis::AnalysisManager &analysis_manager,
      Device::DeviceMgr & device_manager,
      IO::OutputMgr &output_manager
      ) ;

  ~outputsXyceExpressionGroup ();

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

  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, double & retval );
  virtual bool getCurrentVal( const std::string & deviceName, const std::string & designator, std::complex<double> & retval );

  virtual bool getGlobalParameterVal (const std::string & paramName, double & retval );
  virtual bool getGlobalParameterVal (const std::string & paramName, std::complex<double> & retval );

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

  virtual bool getPower(const std::string & deviceName, double & retval);
  virtual bool getPower(const std::string & deviceName, std::complex<double> & retval);

  virtual double getTimeStep ();
  virtual double getTimeStepAlpha () { return alpha_; }
  virtual double getTimeStepPrefac () { return (getTimeStepAlpha() / getTimeStep ()) ; } // FIX

  virtual double getTime();
  virtual double getTemp();
  virtual double getVT  ();
  virtual double getFreq();

  virtual double getBpTol();
  virtual double getStartingTimeStep();
  virtual double getFinalTime();

  virtual unsigned int getStepNumber ();

  void setAliasNodeMap( const IO::AliasNodeMap & anm ) { aliasNodeMap_ = anm; }

  void setOpData (const Util::Op::OpData & opData) { opData_ = opData; }

private:

  int getSolutionGID_(const std::string & nodeName);

  N_PDS_Comm & comm_;
  Topo::Topology & top_;

  Analysis::AnalysisManager & analysisManager_;
  Device::DeviceMgr & deviceManager_;

  IO::AliasNodeMap aliasNodeMap_;

  IO::OutputMgr &outputManager_;

  double time_, temp_, VT_, freq_;
  double dt_, alpha_;

  Util::Op::OpData opData_;
};

}
}

#endif


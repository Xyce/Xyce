
#include <iostream>
#include <unordered_map>
#include <string>

#include <outputsXyceExpressionGroup.h>
#include <ast.h>
#include <newExpression.h>

#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Vector.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Op.h>
#include <N_DEV_Const.h>

#include <mainXyceExpressionGroup.h>

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_DeviceNameConverters.h>

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::outputsXyceExpressionGroup 
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
outputsXyceExpressionGroup::outputsXyceExpressionGroup ( 
 N_PDS_Comm & comm, Topo::Topology & top,
 Analysis::AnalysisManager &analysis_manager,
 Device::DeviceMgr & device_manager,
 IO::OutputMgr &output_manager
 ) :
 comm_(comm),
 top_(top),
 analysisManager_(analysis_manager),
 deviceManager_(device_manager),
 outputManager_(output_manager),
 time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), gmin_(0.0), dt_(0.0), alpha_(0.0)
{
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::outputsXyceExpressionGroup 
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 
//-------------------------------------------------------------------------------
outputsXyceExpressionGroup::~outputsXyceExpressionGroup ()
{
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getSolutionGID_
// Purpose       : 
// Special Notes : 
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
int outputsXyceExpressionGroup::getSolutionGID_(const std::string & nodeName)
{
  int tmpGID=-1;
  std::vector<int> svGIDList1, dummyList;
  char type1;

  std::string nodeNameUpper = nodeName;
  Xyce::Util::toUpper(nodeNameUpper);

  bool foundLocal = top_.getNodeSVarGIDs(NodeID(nodeNameUpper, Xyce::_VNODE), svGIDList1, dummyList, type1);
  bool found = static_cast<int>(foundLocal);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &found, 1);

  // if looking for this as a voltage node failed, try a "device" (i.e. current) node.  I(Vsrc)
  bool foundLocal2 = false;
  if (!found)
  {
    foundLocal2 = top_.getNodeSVarGIDs(NodeID(nodeNameUpper, Xyce::_DNODE), svGIDList1, dummyList, type1);
  }
  bool found2 = static_cast<int>(foundLocal2);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &found2, 1);

  // Check if this is a subcircuit interface node name, which would be found in the aliasNodeMap.
  // If it is then get the GID for the corresponding "real node name". See SRN Bug 1962 for 
  // more details but, as an example, these netlist lines:
  //
  //   X1 1 2 MYSUB
  //   .SUBCKT MYSUB a c 
  //   R1   a b 0.5
  //   R2   b c 0.5
  //  .ENDS
  //
  // would produce key-value pairs of <"X1:C","2"> and <"X1:A","1"> in the aliasNodeMap 
  bool foundAliasNodeLocal = false;
  if (!found && !found2)
  {
    IO::AliasNodeMap::const_iterator alias_it = aliasNodeMap_.find( nodeNameUpper );
    if (alias_it != aliasNodeMap_.end())
    {      
      foundAliasNodeLocal = top_.getNodeSVarGIDs(NodeID((*alias_it).second, Xyce::_VNODE), svGIDList1, dummyList, type1);
    }
  }
  bool foundAliasNode = static_cast<int>(foundAliasNodeLocal);
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_LOR, &foundAliasNode, 1);

  if(svGIDList1.size()==1)
  {
    tmpGID = svGIDList1.front();
  }

  return tmpGID;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getSolutionVal(const std::string & nodeName, double & retval )
{
  ParamList paramList;
  paramList.push_back(Param(std::string("V"),1  ));
  paramList.push_back(Param(      nodeName,0.0));
  Op::OpList internalDeviceVarOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getSolutionVal
// Purpose       : 
// Special Notes : std::complex<double> precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 3/20/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getSolutionVal(const std::string & nodeName, std::complex<double> & retval)
{
  ParamList paramList;
  paramList.push_back(Param(std::string("V"),1  ));
  paramList.push_back(Param(      nodeName,0.0));
  Op::OpList internalDeviceVarOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    double real = Util::Op::getValue(comm_.comm(), *(*it), opData_).real();
    double imag = Util::Op::getValue(comm_.comm(), *(*it), opData_).imag();
    std::complex<double> val = std::complex<double>(real,imag);
    variableValues.push_back( val );
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getCurrentVal
// Purpose       : retrieve the value of device current.  This can be a lead current or a Vsrc
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/26/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getCurrentVal(
    const std::string & deviceName,
    const std::string & designator,
    double & retval )
{
  ParamList paramList;
  paramList.push_back(Param(designator,1));
  paramList.push_back(Param(deviceName,0.0));
  Op::OpList internalDeviceVarOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getCurrentVal
// Purpose       : retrieve the value of device current.  This can be a lead current or a Vsrc
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/26/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getCurrentVal(
    const std::string & deviceName,
    const std::string & designator,
    std::complex<double> & retval )
{
  bool success=true;
  double tmpval;
  getCurrentVal(deviceName, designator, tmpval);
  retval = std::complex<double>(tmpval,0.0);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getGlobalParameterVal
//
// Purpose       : retrieve the value of a parameter that has been 
//                 declared to be a "var" via the make_var function.
//
// Special Notes : double precision version
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getGlobalParameterVal(const std::string &paramName, double & retval)
{
  bool success=true;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, retval);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getGlobalParameterVal
//
// Purpose       : retrieve the value of a parameter that has been 
//                 declared to be a "var" via the make_var function.
//
// Special Notes : std::complex<double> version
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getGlobalParameterVal (const std::string & paramName, std::complex<double> & retval)
{
  bool success=true;
  double tmpval;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, tmpval);
  retval = std::complex<double>(tmpval,0.0);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getInternalDeviceVar
// Purpose       : 
// Special Notes : this is largely copied from the OpBuilder file, specifically 
//                 the Util::Op::Builder::InternalVariableOpBuilder class.
//
//                 There are several things in it that don't make sense.  For example, 
//                 it goes thru a process of checking the solution vector for N() and 
//                 also supports NR(), NI(), NM(), NP(), which I don't think N would 
//                 ever be assiciated with. (but maybe I am missing something)
//
//                 I didn't keep the solution stuff. it can be restored if needed.   
//
//                 Also, it checks the state vector, which I also suspect is wrong.  
//                 But that is more plausible.
//
//                 The one that really matters is the store vector.
//
//                 Finally, for state, it does a MAX_ALL for the state index, 
//                 but it doesn't do it for the store vector index.  Is this right?
//
//                 Another comment; a danger of copying the Op classes is that they
//                 are designed to only work on proc 0 (or in serial).  
//
//                 For stuff on the .print line, the proc 0 restriction is fine, 
//                 but for stuff in Bsrc's and sensitivities it probably is not.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/24/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getInternalDeviceVar (const std::string & deviceName, double & retval )
{
  ParamList paramList;
  paramList.push_back(Param(std::string("N"),1  ));
  paramList.push_back(Param(      deviceName,0.0));
  Op::OpList internalDeviceVarOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getInternalDeviceVar
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/24/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getInternalDeviceVar (const std::string & deviceName, std::complex<double> & retval )
{
  retval=std::complex<double>(0.0,0.0);
  double tmpVal=0.0;
  bool bs1 = getInternalDeviceVar(deviceName, tmpVal);
  retval = std::complex<double>(tmpVal,0.0);
  return bs1;
}

// noise
//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getDnoNoiseDeviceVar
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, double & retval) 
{
  retval=0.0; 

  if (!analysisManager_.getNoiseFlag())
  {
    Report::UserError0() << "DNO and DNI operators only supported for .NOISE analyses";
    return false;
  }
  else
  {
    ParamList paramList;
    paramList.push_back(Param(std::string("DNO"), static_cast<int>(deviceNames.size())));
    for(int ii=0;ii<deviceNames.size();ii++) { paramList.push_back(Param(deviceNames[ii],0.0)); }
    Op::OpList dnoOps_;

    const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(dnoOps_));

    // loop over expressionOps_ to get all the values.
    std::vector<double> variableValues;
    for (Util::Op::OpList::const_iterator it = dnoOps_.begin(); it != dnoOps_.end(); ++it)
    {
      variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
    }

    retval = 0.0;
    if ( !(variableValues.empty()) )
    {
      retval = variableValues[0];
      return true;
    }
    else
    {
      return false;
    }
  }

  return true; 
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getDnoNoiseDeviceVar
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getDnoNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & retval) 
{
  double val=0.0;
  bool retBool = getDnoNoiseDeviceVar(deviceNames,val);
  retval=std::complex<double>(val,0.0); 
  return retBool; 
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getDniNoiseDeviceVar
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, double & retval) 
{
  retval=0.0; 

  if (!analysisManager_.getNoiseFlag())
  {
    Report::UserError0() << "DNO and DNI operators only supported for .NOISE analyses";
    return false;
  }
  else
  {
    ParamList paramList;
    paramList.push_back(Param(std::string("DNI"),static_cast<int>(deviceNames.size())));
    for(int ii=0;ii<deviceNames.size();ii++) { paramList.push_back(Param(deviceNames[ii],0.0)); }
    Op::OpList dniOps_;

    const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(dniOps_));

    // loop over expressionOps_ to get all the values.
    std::vector<double> variableValues;
    for (Util::Op::OpList::const_iterator it = dniOps_.begin(); it != dniOps_.end(); ++it)
    {
      variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
    }

    retval = 0.0;
    if ( !(variableValues.empty()) )
    {
      retval = variableValues[0];
      return true;
    }
    else
    {
      return false;
    }
  }

  return true; 
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getDniNoiseDeviceVar
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getDniNoiseDeviceVar(const std::vector<std::string> & deviceNames, std::complex<double> & retval) 
{
  double val=0.0;
  bool retBool = getDniNoiseDeviceVar(deviceNames,val);
  retval=std::complex<double>(val,0.0); 
  return retBool; 
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getONoise
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getONoise(double & retval) 
{ 
  retval=0.0; 
  if (!analysisManager_.getNoiseFlag())
  {
    Report::UserError0() << "ONOISE operator only supported for .NOISE analyses";
    return false;
  }
  else
  {
    ParamList paramList;
    paramList.push_back(Param(std::string("ONOISE"),0.0));
    Op::OpList onoiseVarOps_;

    const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(onoiseVarOps_));

    // loop over expressionOps_ to get all the values.
    std::vector<double> variableValues;
    for (Util::Op::OpList::const_iterator it = onoiseVarOps_.begin(); it != onoiseVarOps_.end(); ++it)
    {
      variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
    }

    retval = 0.0;
    if ( !(variableValues.empty()) )
    {
      retval = variableValues[0];
      return true;
    }
    else
    {
      return false;
    }
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getONoise
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getONoise(std::complex<double> & retval) 
{
  double val=0.0;
  bool retBool = getONoise(val);
  retval=std::complex<double>(val,0.0); 
  return retBool;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getINoise
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getINoise(double & retval) 
{ 
  retval=0.0; 
  if (!analysisManager_.getNoiseFlag())
  {
    Report::UserError0() << "INOISE operator only supported for .NOISE analyses";
    return false;
  }
  else
  {
    ParamList paramList;
    paramList.push_back(Param(std::string("INOISE"),0.0));
    Op::OpList inoiseVarOps_;

    const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(inoiseVarOps_));

    // loop over expressionOps_ to get all the values.
    std::vector<double> variableValues;
    for (Util::Op::OpList::const_iterator it = inoiseVarOps_.begin(); it != inoiseVarOps_.end(); ++it)
    {
      variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
    }

    retval = 0.0;
    if ( !(variableValues.empty()) )
    {
      retval = variableValues[0];
      return true;
    }
    else
    {
      return false;
    }
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getINoise
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/28/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getINoise(std::complex<double> & retval) 
{
  double val=0.0;
  bool retBool = getINoise(val);
  retval=std::complex<double>(val,0.0); 
  return retBool;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getPower
// Purpose       : 
// Special Notes : Does both Op setup and evaluation; should separate, so only setup 1x
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/12/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getPower(const std::string & deviceName, double & retval)
{
  ParamList paramList;
  paramList.push_back(Param(std::string("P"),1  ));
  paramList.push_back(Param(      deviceName,0.0));
  Op::OpList powerOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(powerOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = powerOps_.begin(); it != powerOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getPower
// Purpose       : 
// Special Notes : complex version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/12/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getPower(const std::string & deviceName, std::complex<double> & retval)
{
  double val=0.0;
  bool retBool = getPower(deviceName, val);
  retval=std::complex<double>(val,0.0); 
  return retBool;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getSparam
// Purpose       :
// Special Notes : double version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/14/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getSparam (const std::vector<int> & args, double & retval )
{
  return false;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getSparam
// Purpose       :
// Special Notes : complex version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/14/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getSparam (const std::vector<int> & args, std::complex<double> & retval )
{
  retval=0.0;
  ParamList paramList;

  paramList.push_back(Param(std::string("S"),static_cast<int>(args.size())));
  for(int ii=0;ii<args.size();ii++) { paramList.push_back(Param(std::to_string(args[ii]),0.0)); }
  Op::OpList sparamOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(sparamOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = sparamOps_.begin(); it != sparamOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_) );
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getYparam
// Purpose       :
// Special Notes : double version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/14/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getYparam (const std::vector<int> & args, double & retval )
{
  return false;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getYparam
// Purpose       :
// Special Notes : complex version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/14/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getYparam (const std::vector<int> & args, std::complex<double> & retval )
{
  retval=0.0;
  ParamList paramList;

  paramList.push_back(Param(std::string("Y"),static_cast<int>(args.size())));
  for(int ii=0;ii<args.size();ii++) { paramList.push_back(Param(std::to_string(args[ii]),0.0)); }
  Op::OpList sparamOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(sparamOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = sparamOps_.begin(); it != sparamOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_) );
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getZparam
// Purpose       :
// Special Notes : double version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/14/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getZparam (const std::vector<int> & args, double & retval )
{
  return false;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getZparam
// Purpose       :
// Special Notes : complex version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/14/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getZparam (const std::vector<int> & args, std::complex<double> & retval )
{
  retval=0.0;
  ParamList paramList;

  paramList.push_back(Param(std::string("Z"),static_cast<int>(args.size())));
  for(int ii=0;ii<args.size();ii++) { paramList.push_back(Param(std::to_string(args[ii]),0.0)); }
  Op::OpList sparamOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(sparamOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = sparamOps_.begin(); it != sparamOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_) );
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    return true;
  }
  else
  {
    return false;
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getTimeStep ()
{
  //dt_ = deviceManager_.getSolverState().currTimeStep_;
  dt_ = outputManager_.getCircuitTimeStep();
  return dt_;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getTime
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getTime() 
{ 
  // I would have preferred to use this but as the code is currently written it
  // is not safe.  The earliest call I would need to make to getTime happens before 
  // the stepErrorControl class has been allocated.  Unfortunately, the analysis
  // manager accessor returns an invalid reference in that case, which I can't 
  // really test for.
  //
  //const TimeIntg::StepErrorControl & secControl_ = (analysisManager_.getStepErrorControl());
  //time_ = secControl_.nextTime;
  //time_ = deviceManager_.getSolverState().currTime_;
  time_ = outputManager_.getCircuitTime();
  return time_;
} 

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getTemp
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getTemp() 
{ 
  temp_ = deviceManager_.getDeviceOptions().temp.getImmutableValue<double>() - CONSTCtoK;
  return temp_;
} 

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getVT  () 
{ 
  VT_ = (deviceManager_.getDeviceOptions().temp.getImmutableValue<double>())*CONSTKoverQ;
  return VT_;
} 

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getFreq() 
{ 
  freq_ = deviceManager_.getSolverState().currFreq_;
  return freq_;
} 

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getGmin
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getGmin() 
{ 
  gmin_ = deviceManager_.getGmin();
  return gmin_;
} 

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getBpTol()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getBpTol()
{
  return deviceManager_.getSolverState().bpTol_;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getStartingTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getStartingTimeStep()
{
  return deviceManager_.getSolverState().startingTimeStep_;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getFinalTime()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/27/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getFinalTime()
{
  return deviceManager_.getSolverState().finalTime_;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getStepNumber()
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/9/2020 
//-------------------------------------------------------------------------------
unsigned int outputsXyceExpressionGroup::getStepNumber()
{
  //return deviceManager_.getSolverState().timeStepNumber_; // either of these should work
  return analysisManager_.getStepNumber();
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getPhaseOutputUsesRadians
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/13/2020 
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getPhaseOutputUsesRadians()
{
  return outputManager_.getPhaseOutputUsesRadians();
}

} // Util
} // Xyce

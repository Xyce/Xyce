
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
 //tempOp_(new Device::ArtificialParameterOp("TEMP", deviceManager_, *(*deviceManager_.getArtificialParameterMap().find("TEMP")).second, "TEMP")),
 time_(0.0), temp_(0.0), VT_(0.0), freq_(0.0), dt_(0.0), alpha_(0.0)
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
  //delete tempOp_;
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
  retval = 0.0;

  int tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
    retval = dataStore_.nextSolutionPtr->getElementByGlobalIndex(tmpGID, 0);
  }
  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &retval, 1);
  return (tmpGID>=0);
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
  double real_val=0.0;
  double imag_val=0.0;

  int tmpGID = getSolutionGID_(nodeName);
  if (tmpGID >= 0)
  {
    const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
    real_val = dataStore_.nextSolutionPtr->getElementByGlobalIndex(tmpGID, 0);
  }

  Xyce::Parallel::AllReduce(comm_.comm(), MPI_SUM, &real_val, 1);

  // ERK.  To Do:
  // need to add logic to see if this is frequency domain situation or not.
  // If not, don't bother getting the imaginary part.
  //
  // Xyce currently uses a real equivalent form for everything, rather than std::complex.
  //
  // For now, however, just set imag to zero.

  retval = std::complex<double>(real_val,imag_val);
  return (tmpGID>=0);
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getCurrentVal
// Purpose       : retrieve the value of device current.  This can be a lead current or a Vsrc
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/26/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getCurrentVal(const std::string & deviceName, double & retval )
{
#if 0
  // try a solution variable first.
  // In the OpBuilder classes, the device name has to be modified substantially 
  // to use the various maps from the output manager.
  // But, for the function getSolutionGID_, which is used by getSolutionVal (called here),
  // this isn't necessary.
  //
  // Note;  Currently, the new Expression class will call "getSolutionVal"  first b4 calling this function.
  // So this function is really just for lead currents.
  bool solSuccess = getSolutionVal(deviceName,retval);
#endif

  // next try state and store (i.e. this is probably a lead current)
  //if(!solSuccess)
  {
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getCurrentVal
// Purpose       : retrieve the value of device current.  This can be a lead current or a Vsrc
// Special Notes : double precision version
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/26/2020
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::getCurrentVal(const std::string & deviceName, std::complex<double> & retval )
{
  bool success=false;
#if 1
  ParamList paramList;
  paramList.push_back(Param(std::string("I"),1  ));
  paramList.push_back(Param(      deviceName,0.0));
  Op::OpList internalDeviceVarOps_;

  const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  Linear::Vector *  solnVecPtr = dataStore_.nextSolutionPtr;
  Linear::Vector *  solnVecImagPtr = dataStore_.nextSolutionPtr; // this is wrong, but leaving for now
  Linear::Vector *  stateVecPtr = dataStore_.nextStatePtr;
  Linear::Vector *  stoVecPtr = dataStore_.nextStorePtr;

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    Util::Op::OpData opDataTmp(0, solnVecPtr, solnVecImagPtr, stateVecPtr, stoVecPtr, 0);
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opDataTmp).real());
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
#else
  // next try state and store (i.e. this is probably a lead current)
  //if(!solSuccess)
  {
    // Node name could be circuit_context:DeviceTypeDeviceName while internally it should be DeviceType:circuit_context:DeviceName.
    std::string tmp = deviceName;
    std::string modifiedName = Xyce::Util::xyceDeviceNameToSpiceName(tmp);
    std::string param_tag = "I"; // ERK.  Fix this.  it can have extra characters.
    static const char * const func_names[] = {"II", "IR", "IP", "IM", "IDB"};  // ERK

    // could be a device lead current "DEV_I" or a branch current.
    // so we don't have to duplicate solution vars(branch currents) in the
    // store vector, look for each type.
    bool param_func = std::find(func_names, func_names + sizeof(func_names)/sizeof(func_names[0]), param_tag) != func_names + sizeof(func_names)/sizeof(func_names[0]);
    std::string store_name = modifiedName + ":DEV_" + (param_func ? "I" : param_tag);  // if it is in the state/store vec.
    std::string leadCurrent_name = modifiedName + ":BRANCH_D";
    if (!param_func && (param_tag.length() > 1) )
    {
      leadCurrent_name += param_tag[1];
    } 

    // this if block allows for spaces in YPDE names as in I1(YPDE NAME)
    // we try to find devices based on store_name in the following blocks of code,
    // so do this modification now.
    std::string::size_type space = store_name.find_first_of(" ");
    if (space != std::string::npos)
    {
      if (space == 4 && store_name.substr(0, 4) == "YPDE")
      {
        store_name.replace(4, 1, ":");
      }
    }
    
    NodeNameMap::const_iterator it;

    // The remaining places to look for lead currents are NOT valid for either .AC
    // or .NOISE analyses.  So, return new_op=0 at this point for those two cases.
    // Note: the code can't also return a UserError0() message at this point.  Otherwise,
    // I(V1) and I(L1) will stop working in parallel for both the .AC and .NOISE cases
    // since a valid new_op is only non-NULL on one processor in parallel.
    if ( (analysisManager_.getNoiseFlag() || analysisManager_.getACFlag()) )
    {
      success=false; // reconsider this
      //return new_op;
    }
    else // These cases are valid for .DC and .TRAN
    {
    
      it = outputManager_.getBranchVarsNodeMap().find(leadCurrent_name );
      // Search lead current vector 
      if (it != outputManager_.getBranchVarsNodeMap().end())
      {
        int index = (*it).second;
        success = true;
      }
      else
      {
        // Search store
        NodeNameMap::const_iterator it = outputManager_.getStoreNodeMap().find(store_name);
        if (it != outputManager_.getStoreNodeMap().end())
        {
          int index = (*it).second;
          success = true;
        }
      }
    }
  }

  return success;
#endif
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

  const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  Linear::Vector *  solnVecPtr = dataStore_.nextSolutionPtr;
  Linear::Vector *  solnVecImagPtr = dataStore_.nextSolutionPtr; // this is wrong, but leaving for now
  Linear::Vector *  stateVecPtr = dataStore_.nextStatePtr;
  Linear::Vector *  stoVecPtr = dataStore_.nextStorePtr;

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    Util::Op::OpData opDataTmp(0, solnVecPtr, solnVecImagPtr, stateVecPtr, stoVecPtr, 0);
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opDataTmp).real());
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
bool outputsXyceExpressionGroup::getDnoNoiseDeviceVar(const std::string & deviceName, double & retval) 
{
  retval=0.0; 

  if (!analysisManager_.getNoiseFlag())
  {
    Report::UserError0() << "DNO and DNI operators only supported for .NOISE analyses";
    return false;
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
bool outputsXyceExpressionGroup::getDnoNoiseDeviceVar(const std::string & deviceName, std::complex<double> & retval) 
{
  double val=0.0;
  bool retBool = getDnoNoiseDeviceVar(deviceName,val);
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
bool outputsXyceExpressionGroup::getDniNoiseDeviceVar(const std::string & deviceName, double & retval) 
{
  retval=0.0; 

  if (!analysisManager_.getNoiseFlag())
  {
    Report::UserError0() << "DNO and DNI operators only supported for .NOISE analyses";
    return false;
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
bool outputsXyceExpressionGroup::getDniNoiseDeviceVar(const std::string & deviceName, std::complex<double> & retval) 
{
  double val=0.0;
  bool retBool = getDniNoiseDeviceVar(deviceName,val);
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

  const TimeIntg::DataStore & dataStore_ = *(analysisManager_.getDataStore());
  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(powerOps_));

  Linear::Vector *  solnVecPtr = dataStore_.nextSolutionPtr;
  Linear::Vector *  solnVecImagPtr = dataStore_.nextSolutionPtr; // this is wrong, but leaving for now
  Linear::Vector *  stateVecPtr = dataStore_.nextStatePtr;
  Linear::Vector *  stoVecPtr = dataStore_.nextStorePtr;

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = powerOps_.begin(); it != powerOps_.end(); ++it)
  {
    Util::Op::OpData opDataTmp(0, solnVecPtr, solnVecImagPtr, stateVecPtr, stoVecPtr, 0);
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opDataTmp).real());
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
// Function      : outputsXyceExpressionGroup::getTimeStep
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 4/20/2020 
//-------------------------------------------------------------------------------
double outputsXyceExpressionGroup::getTimeStep ()
{
  dt_ = deviceManager_.getSolverState().currTimeStep_;
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
  
  time_ = deviceManager_.getSolverState().currTime_;
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
#if 0
  // may not need to bother with the tempOp.  I copied it from the OutputMgrAdapter
  Util::Op::OpData op_data;
  temp_ = (*tempOp_)(comm_.comm(), op_data).real();
#else
  temp_ = deviceManager_.getDeviceOptions().temp.getImmutableValue<double>();
#endif
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

} // Util
} // Xyce

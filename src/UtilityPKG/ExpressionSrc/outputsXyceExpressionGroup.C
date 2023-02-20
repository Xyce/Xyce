//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

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
#include <N_ANP_UQSupport.h>
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
 Parallel::Communicator & comm, Topo::Topology & top,
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
  clearOps();
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void outputsXyceExpressionGroup::clearOps()
{
  for (Util::Op::OpList::const_iterator it = voltageOps_.begin(); it != voltageOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = currentOps_.begin(); it != currentOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = leadCurrentOps_.begin(); it != leadCurrentOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = internalDevVarOps_.begin(); it != internalDevVarOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = dnoNoiseDevVarOps_.begin(); it != dnoNoiseDevVarOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = dniNoiseDevVarOps_.begin(); it != dniNoiseDevVarOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = oNoiseOps_.begin(); it != oNoiseOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = iNoiseOps_.begin(); it != iNoiseOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = powerOps_.begin(); it != powerOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = sparamOps_.begin(); it != sparamOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = yparamOps_.begin(); it != yparamOps_.end(); ++it) { delete *it; }
  for (Util::Op::OpList::const_iterator it = zparamOps_.begin(); it != zparamOps_.end(); ++it) { delete *it; }

  voltageOps_.clear();
  currentOps_.clear();
  leadCurrentOps_.clear();
  internalDevVarOps_.clear();
  dnoNoiseDevVarOps_.clear();
  dniNoiseDevVarOps_.clear();
  oNoiseOps_.clear();
  iNoiseOps_.clear();
  powerOps_.clear();
  sparamOps_.clear();
  yparamOps_.clear();
  zparamOps_.clear();
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::setupGroup
//
// Purpose       : This group sets up all the output Ops.  Not to be confused with
//                 AST ops, which are a different sort of thing, but share the Op
//                 name.
//
// Special Notes : work in progress
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/23/2021
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::setupGroup(newExpression &expr)
{
  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  clearOps();

  if ( !(expr.voltOpVec_.empty()) )
  {
    ParamList paramList;

    for (int ii=0;ii<expr.voltOpVec_.size();ii++)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp
        = Teuchos::rcp_static_cast<voltageOp<usedType> > (expr.voltOpVec_[ii]);

      const std::string & node = voltOp->getVoltageNode();
      if ( !Xyce::Util::checkGroundNodeName(node) ) 
      {
        paramList.push_back(Param(std::string("V"),1  ));
        paramList.push_back(Param(node,0.0));
      }
    }

    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(voltageOps_));
  }


  if ( !(expr.currentOpVec_.empty()) )
  {
    ParamList paramList;

    for (int ii=0;ii<expr.currentOpVec_.size();ii++)
    {
      Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (expr.currentOpVec_[ii]);

      const std::string & deviceName = currOp->getCurrentDevice();
      std::string designator("I");
      paramList.push_back(Param(designator,1));
      paramList.push_back(Param(deviceName,0.0));
    }

    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(currentOps_));
  }

  if ( !(expr.leadCurrentOpVec_.empty()) )
  {
    ParamList paramList;

    for (int ii=0;ii<expr.leadCurrentOpVec_.size();ii++)
    {
      Teuchos::RCP<leadCurrentOp<usedType> > leadCurrOp = Teuchos::rcp_static_cast<leadCurrentOp<usedType> > (expr.leadCurrentOpVec_[ii]);

      const std::string & deviceName = leadCurrOp->getLeadCurrentDevice();
      const std::string & designator = leadCurrOp->getLeadCurrentDesignator();
      paramList.push_back(Param(designator,1));
      paramList.push_back(Param(deviceName,0.0));
    }

    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(leadCurrentOps_));
  }

  if ( !(expr.internalDevVarOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.internalDevVarOpVec_.size();ii++)
    {
      Teuchos::RCP<internalDevVarOp<usedType> > intVarOp = Teuchos::rcp_static_cast<internalDevVarOp<usedType> > (expr.internalDevVarOpVec_[ii]);

      const std::string & deviceName = intVarOp->getInternalVarDevice();
      paramList.push_back(Param(std::string("N"),1  ));
      paramList.push_back(Param(      deviceName,0.0));
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDevVarOps_));
  }

  if ( !(expr.dnoNoiseDevVarOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.dnoNoiseDevVarOpVec_.size();ii++)
    {
      Teuchos::RCP<dnoNoiseVarOp<usedType> > dnoOp = Teuchos::rcp_static_cast<dnoNoiseVarOp<usedType> > (expr.dnoNoiseDevVarOpVec_[ii]);

      const std::vector<std::string> & deviceNames = dnoOp->getNoiseDevices();
      paramList.push_back(Param(std::string("DNO"), static_cast<int>(deviceNames.size())));
      for(int ii=0;ii<deviceNames.size();ii++) { paramList.push_back(Param(deviceNames[ii],0.0)); }
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(dnoNoiseDevVarOps_));
  }

  if ( !(expr.dniNoiseDevVarOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.dniNoiseDevVarOpVec_.size();ii++)
    {
      Teuchos::RCP<dniNoiseVarOp<usedType> > dniOp = Teuchos::rcp_static_cast<dniNoiseVarOp<usedType> > (expr.dniNoiseDevVarOpVec_[ii]);

      const std::vector<std::string> & deviceNames = dniOp->getNoiseDevices();
      paramList.push_back(Param(std::string("DNI"), static_cast<int>(deviceNames.size())));
      for(int ii=0;ii<deviceNames.size();ii++) { paramList.push_back(Param(deviceNames[ii],0.0)); }
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(dniNoiseDevVarOps_));
  }

  if ( !(expr.oNoiseOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.oNoiseOpVec_.size();ii++)
    {
      Teuchos::RCP<oNoiseOp<usedType> > onoiseOp = Teuchos::rcp_static_cast<oNoiseOp<usedType> > (expr.oNoiseOpVec_[ii]);
      paramList.push_back(Param(std::string("ONOISE"),0.0));
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(oNoiseOps_));
  }

  if ( !(expr.iNoiseOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.iNoiseOpVec_.size();ii++)
    {
      Teuchos::RCP<iNoiseOp<usedType> > inoiseOp = Teuchos::rcp_static_cast<iNoiseOp<usedType> > (expr.iNoiseOpVec_[ii]);
      paramList.push_back(Param(std::string("INOISE"),0.0));
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(iNoiseOps_));
  }

  if ( !(expr.powerOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.powerOpVec_.size();ii++)
    {
      Teuchos::RCP<powerOp<usedType> > pwrOp = Teuchos::rcp_static_cast<powerOp<usedType> > (expr.powerOpVec_[ii]);

      const std::string & tag = pwrOp->getPowerTag(); 
      std::string tmpTag = tag;
      Xyce::Util::toUpper(tmpTag);
      const std::string & deviceName = pwrOp->getPowerDevice();

      paramList.push_back(Param(    tmpTag, 1  ));
      paramList.push_back(Param(deviceName,0.0));
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(powerOps_));
  }

  if ( !(expr.sparamOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.sparamOpVec_.size();ii++)
    {
      Teuchos::RCP<sparamOp<usedType> > sparOp = Teuchos::rcp_static_cast<sparamOp<usedType> > (expr.sparamOpVec_[ii]);
      const std::vector<int> & args = sparOp->getSparamArgs();
      paramList.push_back(Param(std::string("S"),static_cast<int>(args.size())));
      for(int ii=0;ii<args.size();ii++) { paramList.push_back(Param(std::to_string(args[ii]),0.0)); }
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(sparamOps_));
  }

  if ( !(expr.yparamOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.yparamOpVec_.size();ii++)
    {
      Teuchos::RCP<yparamOp<usedType> > yparOp = Teuchos::rcp_static_cast<yparamOp<usedType> > (expr.yparamOpVec_[ii]);
      const std::vector<int> & args = yparOp->getYparamArgs();
      paramList.push_back(Param(std::string("Y"),static_cast<int>(args.size())));
      for(int ii=0;ii<args.size();ii++) { paramList.push_back(Param(std::to_string(args[ii]),0.0)); }
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(yparamOps_));
  }

  if ( !(expr.zparamOpVec_.empty()) )
  {
    ParamList paramList;
    for (int ii=0;ii<expr.zparamOpVec_.size();ii++)
    {
      Teuchos::RCP<zparamOp<usedType> > zparOp = Teuchos::rcp_static_cast<zparamOp<usedType> > (expr.zparamOpVec_[ii]);
      const std::vector<int> & args = zparOp->getZparamArgs();
      paramList.push_back(Param(std::string("Z"),static_cast<int>(args.size())));
      for(int ii=0;ii<args.size();ii++) { paramList.push_back(Param(std::to_string(args[ii]),0.0)); }
    }
    Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(zparamOps_));
  }

  return true;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::putValues
// Purpose       : 
// Special Notes : work in progress
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/23/2021
//-------------------------------------------------------------------------------
bool outputsXyceExpressionGroup::putValues(newExpression & expr)
{
  bool noChange=true;

  if ( !(expr.voltOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = voltageOps_.begin();
    for (int ii=0;ii<expr.voltOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp
        = Teuchos::rcp_static_cast<voltageOp<usedType> > (expr.voltOpVec_[ii]);

      const std::string & node = voltOp->getVoltageNode();
      if ( !Xyce::Util::checkGroundNodeName(node) ) 
      {
        usedType & val = voltOp->getVoltageVal();
        usedType oldval = val;
        val = Util::Op::getValue(comm_.comm(), *(*it), opData_); // fix for double.  this assumes std::complex<double>
        if(val != oldval) noChange=false;
      }
    }
  }

  if ( !(expr.currentOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = currentOps_.begin();
    for (int ii=0;ii<expr.currentOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (expr.currentOpVec_[ii]);
      usedType & val = currOp->getCurrentVal();
      usedType oldval = val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_); // fix for double.  this assumes std::complex<double>
      if (val != oldval) noChange=false;
    }
  }

  if ( !(expr.leadCurrentOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = leadCurrentOps_.begin();
    for (int ii=0;ii<expr.leadCurrentOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<leadCurrentOp<usedType> > leadCurrOp = Teuchos::rcp_static_cast<leadCurrentOp<usedType> > (expr.leadCurrentOpVec_[ii]);
      usedType & val = leadCurrOp->getLeadCurrentVar();
      usedType oldval = val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_); // fix for double.  this assumes std::complex<double>
      if (val != oldval) noChange=false;
    }
  }

  if ( !(expr.internalDevVarOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = internalDevVarOps_.begin();
    for (int ii=0;ii<expr.internalDevVarOpVec_.size();ii++,it++)
    {
      Teuchos::RCP<internalDevVarOp<usedType> > intVarOp = Teuchos::rcp_static_cast<internalDevVarOp<usedType> > (expr.internalDevVarOpVec_[ii]);

      usedType & val = intVarOp->getInternalDeviceVar();
      usedType oldval = val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_); // fix for double.  this assumes std::complex<double>
      if (val != oldval) noChange=false;
    }
  }

  // This block of code was originally conceived of for retrieving global parameters.
  // However, it also is used for retrieving anything that is basically of "unknown" status.  
  // For .print line outputs, that includes things like device parameters (for example ISRC:mag)
  //
  // Originally, this block would also handle parameters that had been labled "isVar" via 
  // the "make_vars" function.  This was generally global_param of type Util::DBLE and Util::STR.
  // However,  that is not the case any more as .global_params are now handled 100% with 
  // attachments, and the make_var function is gone.
  if ( !(expr.paramOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.paramOpVec_.size();++ii)
    {
      Teuchos::RCP<paramOp<usedType> > parOp = Teuchos::rcp_static_cast<paramOp<usedType> > (expr.paramOpVec_[ii]);

      if ( !(parOp->getIsAttached()) && !(parOp->getIsConstant()) ) // if the param is a constant, or attached, then it already has its value
      {
        usedType oldval = parOp->getValue();
        usedType val = oldval;
        getParameterVal(parOp->getName(),val);
        parOp->setValue(val);

        if (val != oldval) noChange=false;
      }
    }
  }

  if ( !(expr.dnoNoiseDevVarOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = dnoNoiseDevVarOps_.begin();
    for (int ii=0;ii<expr.dnoNoiseDevVarOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<dnoNoiseVarOp<usedType> > dnoOp = Teuchos::rcp_static_cast<dnoNoiseVarOp<usedType> > (expr.dnoNoiseDevVarOpVec_[ii]);
      usedType & val=dnoOp->getNoiseVar ();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.dniNoiseDevVarOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = dniNoiseDevVarOps_.begin();
    for (int ii=0;ii<expr.dniNoiseDevVarOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<dniNoiseVarOp<usedType> > dniOp = Teuchos::rcp_static_cast<dniNoiseVarOp<usedType> > (expr.dniNoiseDevVarOpVec_[ii]);
      usedType & val=dniOp->getNoiseVar ();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.oNoiseOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = oNoiseOps_.begin();
    for (int ii=0;ii<expr.oNoiseOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<oNoiseOp<usedType> > onoiseOp = Teuchos::rcp_static_cast<oNoiseOp<usedType> > (expr.oNoiseOpVec_[ii]);
      usedType & val=onoiseOp->getNoiseVar();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.iNoiseOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = iNoiseOps_.begin();
    for (int ii=0;ii<expr.iNoiseOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<iNoiseOp<usedType> > inoiseOp = Teuchos::rcp_static_cast<iNoiseOp<usedType> > (expr.iNoiseOpVec_[ii]);
      usedType & val=inoiseOp->getNoiseVar();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.powerOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = powerOps_.begin();
    for (int ii=0;ii<expr.powerOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<powerOp<usedType> > pwrOp = Teuchos::rcp_static_cast<powerOp<usedType> > (expr.powerOpVec_[ii]);
      usedType & val=pwrOp->getPowerVal();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.sparamOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = sparamOps_.begin();
    for (int ii=0;ii<expr.sparamOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<sparamOp<usedType> > sparOp = Teuchos::rcp_static_cast<sparamOp<usedType> > (expr.sparamOpVec_[ii]);
      usedType & val=sparOp->getSparamValue();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.yparamOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = yparamOps_.begin();
    for (int ii=0;ii<expr.yparamOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<yparamOp<usedType> > yparOp = Teuchos::rcp_static_cast<yparamOp<usedType> > (expr.yparamOpVec_[ii]);
      usedType & val=yparOp->getYparamValue();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.zparamOpVec_.empty()) )
  {
    Util::Op::OpList::const_iterator it = zparamOps_.begin();
    for (int ii=0;ii<expr.zparamOpVec_.size();ii++,++it)
    {
      Teuchos::RCP<zparamOp<usedType> > zparOp = Teuchos::rcp_static_cast<zparamOp<usedType> > (expr.zparamOpVec_[ii]);
      usedType & val=zparOp->getZparamValue();
      usedType oldval=val;
      val = Util::Op::getValue(comm_.comm(), *(*it), opData_);
      if (val != oldval) noChange = false;
    }
  }

  return noChange;
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

  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    delete *it;
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
// Special Notes : std::complex<double> version
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

  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    delete *it;
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

  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    delete *it;
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
// Special Notes : std::complex<double> version
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
  ParamList paramList;
  paramList.push_back(Param(designator,1));
  paramList.push_back(Param(deviceName,0.0));
  Op::OpList internalDeviceVarOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(internalDeviceVarOps_));

  // loop over expressionOps_ to get all the values.
  //std::vector<double> variableValues;
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_));
  }

  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    delete *it;
  }

  retval = 0.0;
  if ( !(variableValues.empty()) )
  {
    retval = variableValues[0];
    success=true;
  }
  else
  {
    success=false;
  }

  return success;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getParameterVal
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
bool outputsXyceExpressionGroup::getParameterVal(const std::string &paramName, double & retval)
{
  bool success=true;
  Device::getParamAndReduce(comm_.comm(), deviceManager_, paramName, retval);
  return success;
}

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::getParameterVal
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
bool outputsXyceExpressionGroup::getParameterVal (const std::string & paramName, std::complex<double> & retval)
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

  for (Util::Op::OpList::const_iterator it = internalDeviceVarOps_.begin(); it != internalDeviceVarOps_.end(); ++it)
  {
    delete *it;
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

    for (Util::Op::OpList::const_iterator it = dnoOps_.begin(); it != dnoOps_.end(); ++it)
    {
      delete *it;
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

    for (Util::Op::OpList::const_iterator it = dniOps_.begin(); it != dniOps_.end(); ++it)
    {
      delete *it;
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

    for (Util::Op::OpList::const_iterator it = onoiseVarOps_.begin(); it != onoiseVarOps_.end(); ++it)
    {
      delete *it;
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

    for (Util::Op::OpList::const_iterator it = inoiseVarOps_.begin(); it != inoiseVarOps_.end(); ++it)
    {
      delete *it;
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
bool outputsXyceExpressionGroup::getPower(const std::string & tag, const std::string & deviceName, double & retval)
{
  std::string tmpTag = tag;
  Xyce::Util::toUpper(tmpTag);
  if (tmpTag != "P" && tmpTag != "W") { tmpTag = "P"; }

  ParamList paramList;
  paramList.push_back(Param(    tmpTag, 1  ));
  paramList.push_back(Param(deviceName,0.0));
  Op::OpList powerOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(powerOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<double> variableValues;
  for (Util::Op::OpList::const_iterator it = powerOps_.begin(); it != powerOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_).real());
  }

  for (Util::Op::OpList::const_iterator it = powerOps_.begin(); it != powerOps_.end(); ++it)
  {
    delete *it;
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
bool outputsXyceExpressionGroup::getPower(const std::string & tag, const std::string & deviceName, std::complex<double> & retval)
{
  double val=0.0;
  bool retBool = getPower(tag, deviceName, val);
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

  for (Util::Op::OpList::const_iterator it = sparamOps_.begin(); it != sparamOps_.end(); ++it)
  {
    delete *it;
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
  Op::OpList yparamOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(yparamOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = yparamOps_.begin(); it != yparamOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_) );
  }

  for (Util::Op::OpList::const_iterator it = yparamOps_.begin(); it != yparamOps_.end(); ++it)
  {
    delete *it;
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
  Op::OpList zparamOps_;

  const Util::Op::BuilderManager & op_builder_manager = outputManager_.getOpBuilderManager();
  Util::Op::makeOps(comm_.comm(), op_builder_manager, NetlistLocation(), paramList.begin(), paramList.end(), std::back_inserter(zparamOps_));

  // loop over expressionOps_ to get all the values.
  std::vector<std::complex<double> > variableValues;
  for (Util::Op::OpList::const_iterator it = zparamOps_.begin(); it != zparamOps_.end(); ++it)
  {
    variableValues.push_back( Util::Op::getValue(comm_.comm(), *(*it), opData_) );
  }

  for (Util::Op::OpList::const_iterator it = zparamOps_.begin(); it != zparamOps_.end(); ++it)
  {
    delete *it;
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
  freq_ = outputManager_.getFrequency();
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

//-------------------------------------------------------------------------------
// Function      : outputsXyceExpressionGroup::setRFParamsRequested
// Purpose       : 
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 6/15/2020 
//-------------------------------------------------------------------------------
void outputsXyceExpressionGroup::setRFParamsRequested(std::string type)
{
  analysisManager_.setRFParamsRequested(type);
}

} // Util
} // Xyce

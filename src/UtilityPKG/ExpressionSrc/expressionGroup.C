//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
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

#include "expressionGroup.h"
#include "newExpression.h"

namespace Xyce {
namespace Util {

//-------------------------------------------------------------------------------
// Function      : baseExpressionGroup::putValues
//
// Purpose       : This function puts values into different parts of the AST.
//
// Special Notes : This implementation is the "default" version, which was 
//                 originally coded in the newExpression::getValuesFromGroup 
//                 function.  However, its structure made it difficult to optimize 
//                 the different groups, which is why it got moved here.  
//
//                 Mainly, the groups didn't have enough information to plan 
//                 ahead.  So, searches had to be repeated every time.  With this 
//                 structure (in theory) the searches can be done 1x.  Also, 
//                 in the outputs group, the old structure required that "Ops" 
//                 be allocated every single time.  Now they can be allocated once.
//
//                 Now, each specific group can implement their own version of 
//                 this function, and it should be easier to make each group
//                 more efficient.
//
// Scope         :
// Creator       : Eric Keiter
// Creation Date : 5/23/2021
//-------------------------------------------------------------------------------
bool baseExpressionGroup::putValues(newExpression & expr)
{
  bool noChange=true;

  if ( !(expr.voltOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.voltOpVec_.size();ii++)
    {
      Teuchos::RCP<voltageOp<usedType> > voltOp
        = Teuchos::rcp_static_cast<voltageOp<usedType> > (expr.voltOpVec_[ii]);

      const std::string & node = voltOp->getVoltageNode();
      usedType & val = voltOp->getVoltageVal();
      usedType oldval = val;
      getSolutionVal(node, val);
      if(val != oldval) noChange=false;
    }
  }

  if ( !(expr.currentOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.currentOpVec_.size();ii++)
    {
      Teuchos::RCP<currentOp<usedType> > currOp = Teuchos::rcp_static_cast<currentOp<usedType> > (expr.currentOpVec_[ii]);
      usedType & val = currOp->getCurrentVal();
      usedType oldval = val;
      std::string simple("I");
      getCurrentVal(currOp->getCurrentDevice(),simple,val);
      if (val != oldval) noChange=false;
    }
  }

  if ( !(expr.leadCurrentOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.leadCurrentOpVec_.size();ii++)
    {
      Teuchos::RCP<leadCurrentOp<usedType> > leadCurrOp = Teuchos::rcp_static_cast<leadCurrentOp<usedType> > (expr.leadCurrentOpVec_[ii]);

      usedType val=0.0;
      usedType oldval = leadCurrOp->val();
      getCurrentVal(leadCurrOp->getLeadCurrentDevice(), leadCurrOp->getLeadCurrentDesignator() , val);
      leadCurrOp->setLeadCurrentVar ( val );

      if (val != oldval) noChange=false;
    }
  }

  if ( !(expr.internalDevVarOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.internalDevVarOpVec_.size();ii++)
    {
      Teuchos::RCP<internalDevVarOp<usedType> > intVarOp = Teuchos::rcp_static_cast<internalDevVarOp<usedType> > (expr.internalDevVarOpVec_[ii]);

      usedType val=0.0;
      usedType oldval = intVarOp->val();
      getInternalDeviceVar(intVarOp->getInternalVarDevice(),val);
      intVarOp->setInternalDeviceVar ( val );

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
    for (int ii=0;ii<expr.dnoNoiseDevVarOpVec_.size();ii++)
    {
      Teuchos::RCP<dnoNoiseVarOp<usedType> > dnoOp = Teuchos::rcp_static_cast<dnoNoiseVarOp<usedType> > (expr.dnoNoiseDevVarOpVec_[ii]);

      usedType val=0.0;
      usedType oldval=dnoOp->val();
      getDnoNoiseDeviceVar(dnoOp->getNoiseDevices(),val);
      dnoOp->setNoiseVar ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.dniNoiseDevVarOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.dniNoiseDevVarOpVec_.size();ii++)
    {
      Teuchos::RCP<dniNoiseVarOp<usedType> > dniOp = Teuchos::rcp_static_cast<dniNoiseVarOp<usedType> > (expr.dniNoiseDevVarOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=dniOp->val();
      getDniNoiseDeviceVar(dniOp->getNoiseDevices(),val);
      dniOp->setNoiseVar ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.oNoiseOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.oNoiseOpVec_.size();ii++)
    {
      Teuchos::RCP<oNoiseOp<usedType> > onoiseOp = Teuchos::rcp_static_cast<oNoiseOp<usedType> > (expr.oNoiseOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=onoiseOp->val();
      getONoise(val);
      onoiseOp->setNoiseVar ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.iNoiseOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.iNoiseOpVec_.size();ii++)
    {
      Teuchos::RCP<iNoiseOp<usedType> > inoiseOp = Teuchos::rcp_static_cast<iNoiseOp<usedType> > (expr.iNoiseOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=inoiseOp->val();
      getINoise(val);
      inoiseOp->setNoiseVar ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.powerOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.powerOpVec_.size();ii++)
    {
      Teuchos::RCP<powerOp<usedType> > pwrOp = Teuchos::rcp_static_cast<powerOp<usedType> > (expr.powerOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=pwrOp->val();
      getPower ( pwrOp->getPowerTag(), pwrOp->getPowerDevice(), val);
      pwrOp->setPowerVal ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.sparamOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.sparamOpVec_.size();ii++)
    {
      Teuchos::RCP<sparamOp<usedType> > sparOp = Teuchos::rcp_static_cast<sparamOp<usedType> > (expr.sparamOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=sparOp->val();
      getSparam (sparOp->getSparamArgs(), val);
      sparOp->setValue ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.yparamOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.yparamOpVec_.size();ii++)
    {
      Teuchos::RCP<yparamOp<usedType> > yparOp = Teuchos::rcp_static_cast<yparamOp<usedType> > (expr.yparamOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=yparOp->val();
      getYparam (yparOp->getYparamArgs(), val);
      yparOp->setValue ( val );

      if (val != oldval) noChange = false;
    }
  }

  if ( !(expr.zparamOpVec_.empty()) )
  {
    for (int ii=0;ii<expr.zparamOpVec_.size();ii++)
    {
      Teuchos::RCP<zparamOp<usedType> > zparOp = Teuchos::rcp_static_cast<zparamOp<usedType> > (expr.zparamOpVec_[ii]);
      usedType val=0.0;
      usedType oldval=zparOp->val();
      getZparam (zparOp->getZparamArgs(), val);
      zparOp->setValue ( val );

      if (val != oldval) noChange = false;
    }
  }

  return noChange;
}

}
}

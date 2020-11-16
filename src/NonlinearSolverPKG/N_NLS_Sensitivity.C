//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
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

//-------------------------------------------------------------------------
//
// Purpose        : Body for the sensitivity class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/30/02
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <algorithm>
#include <sstream>
#include <stdexcept>

// ----------   Xyce Includes   ----------

#include <N_ANP_AnalysisManager.h>
#include <N_ANP_OutputMgrAdapter.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_OptionBlock.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_SpiceSeparatedFieldTool.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_CmdParse.h>
#include <N_LAS_Builder.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Problem.h>
#include <N_LAS_Solver.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_FilteredMultiVector.h>
#include <N_LOA_NonlinearEquationLoader.h>
#include <N_NLS_Manager.h>
#include <N_NLS_Sensitivity.h>
#include <N_NLS_SensitivityResiduals.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_StepErrorControl.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_SaveIOSState.h>

#include <expressionGroup.h>
#include <newExpression.h>

// ----------   Static Declarations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// -- stand-alone functions --
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : evaluateObjFuncs
//
// Purpose       : Evaluates user-specified objective function expressions 
//                 and their derivatives with respect to the solution vector x.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
bool evaluateObjFuncs ( 
    std::vector<objectiveFunctionData*> & objVec, 
    N_PDS_Comm & comm,
    Loader::NonlinearEquationLoader & nlEquLoader,
    TimeIntg::DataStore & dataStore,
    TimeIntg::StepErrorControl & sec,
    std::string & netlistFilename
    )
{
  bool bsuccess = true;
  int i;

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->dOdXVectorPtr->putScalar(0.0);
  }

  // obtain the expression variable values.  It will only grab this value if it is owned on this processor.
  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expVarVals.resize (objVec[iobj]->numExpVars, 0.0);
    objVec[iobj]->expVarDerivs.resize (objVec[iobj]->numExpVars, 0.0);

    for (i = 0; i < objVec[iobj]->numExpVars; ++i)
    {
      double tmpVal=0.0;
      if ( objVec[iobj]->globalParamVariableStencil[i] == 1) // this is a global param
      {
        // This is a bit of overkill, as "getParamAndReduce" refers to all kinds of parameters, not just global 
        // parameters.  So this call will result in a search over a larger container than necessary.  (ON the other 
        // hand this might be useful later.  But the entity which "owns" the global parameter values is the device package.
        // The "global param" Op is a class in Xyce::Device.  I'm not using ops here, but am relying on the device package
        // being fully aware of global parmaeter values.
        nlEquLoader.getParamAndReduce(comm.comm(), objVec[iobj]->expVarNames[i], tmpVal);
      }
      else
      {
        int tmpGID=objVec[iobj]->expVarGIDs[i];
        if (tmpGID >= 0)
        {
          tmpVal = dataStore.nextSolutionPtr->getElementByGlobalIndex(tmpGID, 0);
        }
        else 
        {
          tmpVal = 0.0;
        }
        Xyce::Parallel::AllReduce(comm.comm(), MPI_SUM, &tmpVal, 1);
      }
      objVec[iobj]->expVarVals[i] = tmpVal;
    }
  }

  //get expression value and partial derivatives
  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expPtr->evaluate( 
        objVec[iobj]->expVal, 
        objVec[iobj]->expVarDerivs); 

    objVec[iobj]->objFuncEval = objVec[iobj]->expVal;
    objVec[iobj]->dOdXVectorPtr->putScalar(0.0);
    for (i=0;i<objVec[iobj]->numExpVars;++i)
    {
      int tmpGID = objVec[iobj]->expVarGIDs[i];
      double tmpDODX = objVec[iobj]->expVarDerivs[i];

      if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
      {
        Xyce::dout() 
          <<  objVec[iobj]->expVarNames[i] << "  "
          << "i="<<i<<"  gid = " << tmpGID << "  dodx = "<< tmpDODX << std::endl;
      }

      if (tmpGID >= 0)
      {
        objVec[iobj]->dOdXVectorPtr->setElementByGlobalIndex(tmpGID, tmpDODX, 0);
      }
    }

    // Assuming this is zero:
    objVec[iobj]->dOdp = 0.0;
  }

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->dOdXVectorPtr->fillComplete();

    if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
    {
      std::string filename = netlistFilename + "_dodx.txt";
      objVec[iobj]->dOdXVectorPtr->writeToFile(const_cast<char *>(filename.c_str()));
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : setupObjectiveFunctions
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/25/2019
//-----------------------------------------------------------------------------
void setupObjectiveFunctions(
  Teuchos::RCP<Xyce::Util::baseExpressionGroup> & exprGroup,
  std::vector<objectiveFunctionData*> & objVec,
  IO::OutputMgr & output_manager,
  Linear::System & lasSys,
  const IO::CmdParse &cp,
  bool checkTimeDeriv
  )
{

  Xyce::ExtendedString exprType = cp.getArgumentValue( "-expression" );
  exprType.toUpper();

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expPtr = new Util::Expression(exprGroup, objVec[iobj]->objFuncString);

    if (!(objVec[iobj]->expPtr->parsed()))
    {
      Report::UserFatal0()
        << "Objective function " << objVec[iobj]->objFuncString
        << " is not parsable.";
    }

    // resolve user-defined functions first
    {
    std::vector<std::string> global_function_names;
    objVec[iobj]->expPtr->getFuncNames(global_function_names);
 
    std::vector<std::string>::iterator it = global_function_names.begin();
    std::vector<std::string>::iterator end = global_function_names.end();
    const Util::ParamMap & context_function_map = output_manager.getMainContextFunctionMap();

    for ( ; it != end; ++it)
    {
      Util::ParamMap::const_iterator paramMapIter = context_function_map.find(*it);

      if (paramMapIter == context_function_map.end())
      {
        Report::UserError0() << "Cannot find global function definition for " << *it 
          << " in objective function expression " << objVec[iobj]->objFuncString;
        break;
      }

      const Util::Param &functionParameter = (*paramMapIter).second;

      std::string functionPrototype(functionParameter.tag());
      std::string functionBody(functionParameter.stringValue());
      Util::Expression prototypeExression(exprGroup, functionPrototype);
      std::vector<std::string> arguments = prototypeExression.getFunctionArgStringVec ();

      // in the parameter we found, pull out the RHS expression and attach
      if(functionParameter.getType() == Xyce::Util::EXPR)
      {
        Util::Expression & expToBeAttached 
          = const_cast<Util::Expression &> (functionParameter.getValue<Util::Expression>());

        // attach the node
        objVec[iobj]->expPtr->attachFunctionNode(*it, expToBeAttached);
      }
      else
      {
        std::cout << "functionParameter is not EXPR type!!!" <<std::endl;

        switch (functionParameter.getType()) 
        {
          case Xyce::Util::STR:
            std::cout <<"It is STR type: " <<  functionParameter.stringValue();
            break;
          case Xyce::Util::DBLE:
            std::cout <<"It is DBLE type: " <<  functionParameter.getImmutableValue<double>();
            break;
          case Xyce::Util::EXPR:
            std::cout <<"It is EXPR type: " << functionParameter.getValue<Util::Expression>().get_expression();
            break;
          default:
            std::cout <<"It is default type (whatever that is): " << functionParameter.stringValue();
        }
      }
    }
    }

#if 0
    objVec[iobj]->numDdt = objVec[iobj]->expPtr->getNumDdt();
    if ( checkTimeDeriv )
    {
      if (objVec[iobj]->numDdt > 0)
      {
        Report::DevelFatal() <<  "Objective function contains a time derivative, which cannot be processed.";
      }
    }
#endif

    // setup the names:
    objVec[iobj]->expVarNames.clear();

    std::vector<std::string> nodes;
    std::vector<std::string> instances;

    objVec[iobj]->expPtr->getVoltageNodes(nodes);
    objVec[iobj]->expPtr->getDeviceCurrents(instances);

    // Make the current (instance) strings all upper case.
    // The topology directory apparently requires this.
    std::vector<std::string>::iterator iter;
    for (iter=instances.begin();iter!=instances.end();++iter)
    {
      ExtendedString tmpString = *iter;
      tmpString.toUpper ();
      *iter  = tmpString;
    }

    objVec[iobj]->expVarNames.insert(objVec[iobj]->expVarNames.end(), nodes.begin(), nodes.end());
    objVec[iobj]->expVarNames.insert(objVec[iobj]->expVarNames.end(), instances.begin(), instances.end());

    // now handle params and global params
    std::vector<std::string> globalParams;
    const std::vector<std::string> & strings = objVec[iobj]->expPtr->getUnresolvedParams();

    const Util::ParamMap & context_param_map = output_manager.getMainContextParamMap();
    const Util::ParamMap & context_global_param_map = output_manager.getMainContextGlobalParamMap();
    for (int istring=0;istring<strings.size();istring++)
    {
      Util::ParamMap::const_iterator param_it = context_param_map.find(strings[istring]);

      if (param_it != context_param_map.end())
      {
        const Util::Param &replacement_param = (*param_it).second;

        if ( replacement_param.getType() == Xyce::Util::STR ||
             replacement_param.getType() == Xyce::Util::DBLE )
        {
          enumParamType paramType=DOT_PARAM;
          if (!objVec[iobj]->expPtr->make_constant(strings[istring], replacement_param.getImmutableValue<double>()),paramType)
          {
            Report::UserWarning0() << "Problem converting parameter " << strings[istring] << " to its value.";
          }
        }
        else if (replacement_param.getType() == Xyce::Util::EXPR)
        {
#if 0
          if (objVec[iobj]->expPtr->replace_var(strings[istring], replacement_param.getValue<Util::Expression>()) != 0)
          {
            std::string expressionString=objVec[iobj]->expPtr->get_expression();
            Report::UserWarning0() << "Problem inserting expression " << replacement_param.getValue<Util::Expression>().get_expression()
                                   << " as substitute for " << strings[istring] << " in expression " << expressionString;
          }
#else
          enumParamType paramType=DOT_PARAM;
          objVec[iobj]->expPtr->attachParameterNode (strings[istring], replacement_param.getValue<Util::Expression>(),paramType);
#endif
        }
      }
      else
      {
        // if this string is found in the global parameter map, then attach it to the expression
        param_it = context_global_param_map.find(strings[istring]);
        if (param_it != context_global_param_map.end())
        {
          globalParams.push_back(strings[istring]);

          if(param_it->second.getType() == Xyce::Util::EXPR)
          {
            Util::Expression & expToBeAttached = const_cast<Util::Expression &> (param_it->second.getValue<Util::Expression>());
            objVec[iobj]->expPtr->attachParameterNode(strings[istring], expToBeAttached);
          }
          else
          {
            if (!objVec[iobj]->expPtr->make_var(strings[istring]))
            {
              Report::UserWarning0() << "Problem setting global parameter " << strings[istring];
            }
          }
        }
        else
        {
          Report::UserWarning0() << "This field: " << strings[istring] 
            << " from the objective function " << objVec[iobj]->objFuncString << " is not resolvable";
        }
      }
    }

    objVec[iobj]->expVarNames.insert(objVec[iobj]->expVarNames.end(), globalParams.begin(), globalParams.end());

    objVec[iobj]->numExpVars = objVec[iobj]->expVarNames.size();

    if (objVec[iobj]->numExpVars<=0)
    {
      Report::UserFatal0()
        <<  "Objective function does not contain a resolvable solution variable.";
    }

    objVec[iobj]->dOdXVectorPtr = lasSys.builder().createVector();
  }
}

//-----------------------------------------------------------------------------
// Function      : setupObjectiveFuncGIDs
//
// Purpose       : 
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void setupObjectiveFuncGIDs (
   std::vector<objectiveFunctionData*> & objVec, 
   N_PDS_Comm & comm,
   Topo::Topology & top,
   IO::OutputMgr & output_manager)
{
  int found(0);
  int found2(0);
  int found3(0);
  int foundAliasNode(0);
  bool foundLocal(false);
  bool foundLocal2(false);
  bool foundLocal3(false);
  bool foundAliasNodeLocal(false);


  // ERK. this code is pretty silly, in at least one respect.  We already know (or could know) 
  // what type of variable each member of the expVarNames is.  This was figured out already in the 
  // setupObjectiveFunctions function, as we pulled "nodes" and "instances" and "strings" out of 
  // the expression and did specific things based on what they were.


  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    // set up the gid's:
    objVec[iobj]->expVarGIDs.resize( objVec[iobj]->numExpVars, -1);
    objVec[iobj]->globalParamVariableStencil.resize( objVec[iobj]->numExpVars, 0);

    for (int i = 0; i < objVec[iobj]->numExpVars; ++i)
    {
      std::vector<int> svGIDList1, dummyList;
      char type1;

      // look for this variable as a node first.
      foundLocal = top.getNodeSVarGIDs(NodeID(objVec[iobj]->expVarNames[i], Xyce::_VNODE), svGIDList1, dummyList, type1);
      found = static_cast<int>(foundLocal);
      Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found, 1);

      // if looking for this as a voltage node failed, try a "device" (i.e. current) node.  I(Vsrc)
      foundLocal2 = false;
      if (!found)
      {
        foundLocal2 = top.getNodeSVarGIDs(NodeID(objVec[iobj]->expVarNames[i], Xyce::_DNODE), svGIDList1, dummyList, type1);
      }
      found2 = static_cast<int>(foundLocal2);
      Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found2, 1);

      // check global param list.  If it is a global parameter, then we don't need the GID. 
      // The expression library, treats global parameters as variables.  
      // It treats non-global parameters as constants.  In the way the expression library
      // operates, if something is considered a variable, then it will compute derivatives of the
      // expression with respect to it.  
      //
      // However, as the sensitivity calculation only cares about derivatives that can be propagated thru
      // the Jacobian, these derivatives aren't needed for sensitivity calculations.
      foundLocal3 = false;
      if (!found && !found2)
      {
        const Util::ParamMap & context_global_param_map = output_manager.getMainContextGlobalParamMap();
        Util::ParamMap::const_iterator param_it = context_global_param_map.find(objVec[iobj]->expVarNames[i] ); 
        foundLocal3 = (param_it != context_global_param_map.end());
      }
      found3 = static_cast<int>(foundLocal3);
      Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &found3, 1);

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
      foundAliasNodeLocal = false;
      if (!found && !found2 && !found3)
      {
        IO::AliasNodeMap::const_iterator alias_it = output_manager.getAliasNodeMap().find(objVec[iobj]->expVarNames[i]);
        if (alias_it != output_manager.getAliasNodeMap().end())
        {      
          foundAliasNodeLocal = top.getNodeSVarGIDs(NodeID((*alias_it).second, Xyce::_VNODE), svGIDList1, dummyList, type1);
        }
      }
      foundAliasNode = static_cast<int>(foundAliasNodeLocal);
      Xyce::Parallel::AllReduce(comm.comm(), MPI_LOR, &foundAliasNode, 1);

      if (!found && !found2 && !found3 && !foundAliasNode)
      {
        Report::UserFatal() << "objective function variable not found!  Cannot find " << objVec[iobj]->expVarNames[i] ;
      }

      if (found || found2 || foundAliasNode)
      {
        int tmpGID=-1;
        if(svGIDList1.size()==1)
        {
          tmpGID = svGIDList1.front();
        }
        objVec[iobj]->expVarGIDs[i] = tmpGID;
      }

      if (found3)
      {
        objVec[iobj]->expVarGIDs[i] = -100;
        objVec[iobj]->globalParamVariableStencil[i] = 1;
      }
      else
      {
        objVec[iobj]->globalParamVariableStencil[i] = 0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : applyHocevarDelayTerms
//
// Purpose       : transform the computed solution sensitivities into 
//                 delay sensitivities
//
// Special Notes : This implements the Hocevar version of delay sensitivies
//
//                 If you want a sensitivity for when a voltage crosses a 
//                 specified value, use this approach.
//
//                 V_thresh = v(p,tau), where tau=time that v(p,tau) crosses V_thresh
//
//                 take derivative w.r.t. parameter, p.  Then this gives:
//
//                 0 = dv/dp + dv/dtau * dtau/dp
//
//                 rearrange:
//
//                 dtau/dp = -dv/dp * (v')^-1   where v' is the time derivative of v(tau)
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
void applyHocevarDelayTerms(
    std::vector<objectiveFunctionData*> & objVec,
    std::vector<objectiveFunctionData*> & objTimeDerivVec,
    TimeIntg::DataStore & ds
    )
{
  int size=objVec.size();
  int numSensParams=ds.dOdpVec_.size();

  for (int ii=0;ii<size;++ii)
  {
    double timeDeriv = objTimeDerivVec[ii]->objFuncEval;
    double tdRecip = (timeDeriv!=0.0)?(1/timeDeriv):0.0;

    for (int iparam=0; iparam< numSensParams; ++iparam)
    {
      ds.dOdpVec_[iparam] *= -tdRecip;
      ds.scaled_dOdpVec_[iparam] *= -tdRecip;
    }
  }
}

//-----------------------------------------------------------------------------
// -- Sensitivity class functions --
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Sensitivity::Sensitivity
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
Sensitivity::Sensitivity (
  NonLinearSolver *     nls,
  Topo::Topology &      topTmp,
  const IO::CmdParse &  cp,
  TimeIntg::StepErrorControl & secTmp
  )
  : NonLinearSolver(cp),
    debugLevel_(1),
    solutionSize_(0),
    solveDirectFlag_(false),
    solveAdjointFlag_(true),
    outputScaledFlag_(false),
    outputUnscaledFlag_(true),
    maxParamStringSize_(0),
    stdOutputFlag_(false),
    fileOutputFlag_(false),
    numSolves_(0),
    objFuncGiven_(false),
    objFuncGIDsetup_(false),
    objFuncTimeDerivGIDsetup_(false),
    difference_(SENS_FWD),
    sqrtEta_(1.0e-8),
    sqrtEtaGiven_(false),
    forceFD_(false),
    forceDeviceFD_(false),
    forceAnalytic_(false),
    newLowMem_(false),
    sparseAdjointStorage_(true),
    computeDelays_(false),
    timeDerivsSetup_(false),
    reuseFactors_(true),
    lambdaVectorPtr_(0),
    savedRHSVectorPtr_(0),
    savedNewtonVectorPtr_(0),
    nls_(nls),
    top_(topTmp),
    sec(secTmp),
    numSensParams_(0),
    numObjectives_(0),
    mode_(Nonlinear::DC_OP)
{
  resetNLS( nls );

  registerAnalysisManager(&nls_->getAnalysisManager());

  savedRHSVectorPtr_    = lasSysPtr_->builder().createVector();
  savedNewtonVectorPtr_ = lasSysPtr_->builder().createVector();

  lambdaVectorPtr_ = lasSysPtr_->builder().createVector();

  solutionSize_ = savedRHSVectorPtr_->localLength();
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::~Sensitivity
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
Sensitivity::~Sensitivity()
{
  delete savedRHSVectorPtr_;
  savedRHSVectorPtr_ = 0;

  delete savedNewtonVectorPtr_;
  savedNewtonVectorPtr_ = 0;

  delete lambdaVectorPtr_;
  lambdaVectorPtr_ = 0;

  for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
  {
    delete objFuncDataVec_[iobj]->dOdXVectorPtr;
    objFuncDataVec_[iobj]->dOdXVectorPtr = 0;

    delete objFuncDataVec_[iobj]->expPtr;
    objFuncDataVec_[iobj]->expPtr = 0;

    delete objFuncDataVec_[iobj];
    objFuncDataVec_[iobj] = 0;
  } 

  // For all the stuff that is to be deleted in the nonlinear solver
  // base class destructor, just set those pointers to zero because
  // they will have been deleted already.
  //
  // This is the consequence of having the sensitivity class derived
  // off of the nonlinear solver base class, but having it use all
  // the same linear solver objects, etc., as the nonlinear solver used
  // to solve the problem.
  gradVectorPtr_       = 0;
  solWtVectorPtr_      = 0;
  linsolOptionBlockPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setNonlinearSolver
// Purpose       : reset internal nonlinear solver reference
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
void Sensitivity::resetNLS( NonLinearSolver * nls )
{
  // Update reference to new nonlinear solver.
  nls_ = nls;

  // Update pointers.
  lasSysPtr_    = nls_->lasSysPtr_;
  nonlinearEquationLoader_    = nls_->nonlinearEquationLoader_;
  rhsVectorPtr_ = nls_->rhsVectorPtr_;

  NewtonVectorPtr_     = nls_->NewtonVectorPtr_;
  lasSolverRCPtr_      = nls_->lasSolverRCPtr_;
  jacobianMatrixPtr_   = nls_->jacobianMatrixPtr_;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::stdOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/18/2014
//-----------------------------------------------------------------------------
std::ostream& Sensitivity::stdOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities,
       std::ostream& os
       )
{
  Analysis::OutputMgrAdapter & outputManagerAdapter = getAnalysisManager().getOutputManagerAdapter();

  // save current stream state, and then set the stream to use scientific notation.
  // Otherwise the info for the stepped parameters may not be output in
  // scientific notation.
  basic_ios_all_saver<std::ostream::char_type> save(os);
  os.setf(std::ios::scientific);

  if ( !(outputManagerAdapter.getStepSweepVector().empty()) )
  {
    // Output for .STEP.  Step count should be output as an integer.
    // The other values should be output in scientific notation.
    os << "\nFor Step " << outputManagerAdapter.getStepAnalysisStepNumber() << ":" << std::endl;
    for (std::vector<Analysis::SweepParam>::const_iterator it = outputManagerAdapter.getStepSweepVector().begin();
         it != outputManagerAdapter.getStepSweepVector().end(); ++it)
    {
      os << it->name << " = " << it->currentVal << std::endl;
    }
  }

  for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
  {
    os << "\n"<<idString << " Sensitivities of objective function:";
    os
      << objFuncDataVec_[iobj]->objFuncString << " = " 
      << std::setw(13)<< std::scientific<< std::setprecision(4)
      << objFuncDataVec_[iobj]->objFuncEval << std::endl;


    os << std::setw(maxParamStringSize_)<<"Name";
    os << "\t"<<std::setw(13)<<"Value";
    os << "\t"<<std::setw(13)<<"Sensitivity";
    os << "\t"<<std::setw(13)<<"Normalized"<<std::endl;

    for (int iparam=0; iparam< numSensParams_; ++iparam)
    {
      os << std::setw(maxParamStringSize_)<<paramNameVec_[iparam];

      os << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << paramVals[iparam];

      int index= iobj*numSensParams_ +iparam;
      os << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << sensitivities[index];

      os << "\t" << std::setw(13)<< std::scientific<< std::setprecision(4)
        << scaled_sensitivities[index] << std::endl;
    }
  }

  if ( !(outputManagerAdapter.getStepSweepVector().empty()) )
  {
    if ( (outputManagerAdapter.getStepAnalysisStepNumber()+1) <
          outputManagerAdapter.getStepAnalysisMaxSteps() )
    {
      // add a blank line after each block of .STEP output, except for the last one
      os << std::endl;
    }
  }

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::fileOutput
// Purpose       : Dump sensitivity information out to a file.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/18/2014
//-----------------------------------------------------------------------------
void Sensitivity::fileOutput (
       std::string idString,
       std::vector<double> & paramVals,
       std::vector<double> & sensitivities,
       std::vector<double> & scaled_sensitivities)
{
  N_PDS_Comm & comm = *(pdsMgrPtr_->getPDSComm());
  int myPID = comm.procID();
  if (myPID==0)
  {
    std::ostringstream numSolvesOStr;
    numSolvesOStr << numSolves_;
    std::string dodpFileName = netlistFilename_ + numSolvesOStr.str() + "_dodp" + idString +".txt";
    FILE *fp = fopen(dodpFileName.c_str(),"w");
    for (int iparam=0;iparam< numSensParams_; ++iparam)
    {
      fprintf(fp,"\t%16.8e\n", sensitivities[iparam]);
    }
    fclose(fp);
  }
  comm.barrier();
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::icSensitivity
// Purpose       : This function is called when there is a NOOP or UIC.  The
//                 dqdp vector needs to be set up for the history to be correct.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool Sensitivity::icSensitivity ( 
     std::vector<double> & objectiveVec,
     std::vector<double> & dOdpVec, 
     std::vector<double> & dOdpAdjVec,
     std::vector<double> & scaled_dOdpVec, 
     std::vector<double> & scaled_dOdpAdjVec)
{
  if (!solveDirectFlag_ && !solveAdjointFlag_) return 1;

  TimeIntg::DataStore & ds = *dsPtr_;

  ds.dOdpVec_.clear();
  ds.dOdpAdjVec_.clear();

  ds.scaled_dOdpVec_.clear();
  ds.scaled_dOdpAdjVec_.clear();

  // first get the derivatives of the DAE RHS vector w.r.t. the
  // user-specified optimization parameters.
  loadSensitivityResiduals (difference_, forceFD_, forceDeviceFD_, forceAnalytic_, 
      sqrtEta_, netlistFilename_, 
      *dsPtr_, *nonlinearEquationLoader_, paramNameVec_, getAnalysisManager());

  calcObjFuncDerivs ();
  if (computeDelays_)
  {
    calcObjFuncTimeDerivs ();
  }

  objectiveVec.clear();
  ds.objectiveVec_.clear();
  for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
  {
    objectiveVec.push_back(objFuncDataVec_[iobj]->objFuncEval);
    ds.objectiveVec_.push_back(objFuncDataVec_[iobj]->objFuncEval);
  }

  if (solveDirectFlag_) 
  {
    ds.dOdpVec_.resize(numSensParams_*numObjectives_ ,0.0);
    ds.scaled_dOdpVec_.resize(numSensParams_*numObjectives_,0.0);

    if (outputUnscaledFlag_)
    {
      dOdpVec = ds.dOdpVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpVec = ds.scaled_dOdpVec_;
    }

    if (stdOutputFlag_)
    {
      stdOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_, Xyce::lout() );
    }
  }

  if (solveAdjointFlag_) 
  {
    ds.dOdpAdjVec_.resize(numSensParams_*numObjectives_,0.0);
    ds.scaled_dOdpAdjVec_.resize(numSensParams_*numObjectives_,0.0);

    if (outputUnscaledFlag_)
    {
      dOdpAdjVec = ds.dOdpAdjVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpAdjVec = ds.scaled_dOdpAdjVec_;
    }
    if (stdOutputFlag_  && mode_ != Nonlinear::TRANSIENT)
    {
      stdOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_, Xyce::lout() );
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Sensitivity::solve
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/21/02
//-----------------------------------------------------------------------------
int Sensitivity::solve( 
     std::vector<double> & objectiveVec,
     std::vector<double> & dOdpVec, 
     std::vector<double> & dOdpAdjVec,
     std::vector<double> & scaled_dOdpVec, 
     std::vector<double> & scaled_dOdpAdjVec)
{
  Stats::StatTop _solveStat("Sensistivity Solve");
  Stats::TimeBlock _solveTimer(_solveStat);

  if (!solveDirectFlag_ && !solveAdjointFlag_) return 1;

  TimeIntg::DataStore & ds = *dsPtr_;

  ds.dOdpVec_.clear();
  ds.dOdpAdjVec_.clear();

  ds.scaled_dOdpVec_.clear();
  ds.scaled_dOdpAdjVec_.clear();

  // It may now be neccessary to re-load the jacobian and rhs vectors.
  // It is necccessary, for example, if the two-level Newton solver is the
  // solver being used.
  nls_->enableSensitivity ();

  // first get the derivatives of the DAE RHS vector w.r.t. the
  // user-specified optimization parameters.
  loadSensitivityResiduals (difference_, 
      forceFD_, 
      forceDeviceFD_, 
      forceAnalytic_, 
      sqrtEta_, netlistFilename_, 
      *dsPtr_, *nonlinearEquationLoader_, paramNameVec_, getAnalysisManager());

  calcObjFuncDerivs ();
  if (computeDelays_)
  {
    calcObjFuncTimeDerivs ();
  }

  objectiveVec.clear();
  ds.objectiveVec_.clear();
  for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
  {
    objectiveVec.push_back(objFuncDataVec_[iobj]->objFuncEval);
    ds.objectiveVec_.push_back(objFuncDataVec_[iobj]->objFuncEval);
  }

  if (solveDirectFlag_) 
  {
    solveDirect ();

    if (computeDelays_)
    {
      applyHocevarDelayTerms( objFuncDataVec_, objFuncTimeDerivDataVec_, ds );
    }

    if (outputUnscaledFlag_)
    {
      dOdpVec = ds.dOdpVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpVec = ds.scaled_dOdpVec_;
    }
  }

  if (solveAdjointFlag_) 
  {
    solveAdjoint ();
    if (outputUnscaledFlag_)
    {
      dOdpAdjVec = ds.dOdpAdjVec_;
    }

    if (outputScaledFlag_)
    {
      scaled_dOdpAdjVec = ds.scaled_dOdpAdjVec_;
    }
  }

  numSolves_++;

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::solveDirect
//
// Purpose       : This function calculates the direct sensitivities for
//                 the user specified parameters.  It is used for both 
//                 steady-state and transient direct sensitivities.
//
//                 The ultimate goal of this function is to obtain dO/dp,
//                 where O is the objective function, and p is a
//                 user-defined optimization parameter.
//
//                 This is a bit confusing because the DAKOTA folks use
//                 a different naming convention than I have tended to
//                 use.  In Dakota's documentation, dO/dp would be referred
//                 to as df/du.  In Xyce, f is already considered to be the
//                 right-hand-side residual vector.  For clarity, this is
//                 the key between my notation and DAKOTA's:
//                        DAK     Xyce
//                         f       O
//                         y       x
//                         u       p
//                         c       f
//
//                 To obtain dOdp, the device manager is first called, and
//                 told to calculate dfdp, which is the derivative of the
//                 residual vector w.r.t. user-defined params.  It does
//                 this by numerical differentiation.  Then, after that,
//                 this function solves the equation:
//
//                 J dx/dp = df/dp    ->   dx/dp = J^-1 df/dp
//
//                 for each p.  
//
//                 J is the jacobian matrix, so J=df/dx.  (dc/dy)
//
//                 After performing these linear solves, dO/dp is to be
//                 obtained using the chain rule by:
//
//                 dO/dp = - dO/dx * dx/dp  + dO/dp
//
//                 The O in the dO/dp on the left side of this equation
//                 should have a hat over it "^", to indicate that it is
//                 different than the O on the right hand side.
//
//                 Note, this method is best if you have lots of objective
//                 functions, and a small number of parameters.  For adjoint
//                 calculations it is the other way around.
//
//                 11/19/02.
//                 It is assumed (for now) that dO/dp on the right hand side
//                 is zero, i.e., there is no direct
//                 dependence on p by O.  So, for example, if a user
//                 defined parameter to be used is the length of a MOSFET,
//                 the MOSFET length will NOT appear in the analytical
//                 expression for O.
//
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
int Sensitivity::solveDirect ()
{
  Stats::StatTop _solveDirectStat("Solve Direct");
  Stats::TimeBlock _solveDirectTimer(_solveDirectStat);

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    Xyce::dout() << std::endl
                 << "In Sensitivity::solveDirect" << std::endl;
  }

  TimeIntg::DataStore & ds = *dsPtr_;

  int iparam;
  Linear::MultiVector * dXdpPtrVector = ds.nextDXdpPtrVector;
  Linear::MultiVector * sensRHSPtrVector = ds.sensRHSPtrVector;

  // first save a copy of the rhs and newton vectors, in case we want them later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_), 0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_), 0.0);

  // via the loader, setup all the sensitivity residuals.
  nonlinearEquationLoader_->loadSensitivityResiduals();

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    std::string matrixFile = netlistFilename_ + "_directMatrix.txt";
    jacobianMatrixPtr_->writeToFile(const_cast<char *>(matrixFile.c_str()));
  }


  // Now solve the series of linear systems to get dXdp.
  for (iparam=0; iparam< numSensParams_; ++iparam)
  {
    // copy the sensitivity residual into the rhs vector location,
    // as that is the RHS vector that the linear solver expects.
    *rhsVectorPtr_ = *Teuchos::rcp(sensRHSPtrVector->getNonConstVectorView(iparam));
    Teuchos::RCP<Linear::Vector> dXdp = Teuchos::rcp( dXdpPtrVector->getNonConstVectorView(iparam) );

    lasSolverRCPtr_->solve(reuseFactors_);

    // copy the result of the linear solve into the dxdp data structure.
    *(dXdp) = *(NewtonVectorPtr_);

    // do debug output.
    if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
    {
      Xyce::dout() << "iparam="<<iparam << "\t" << paramNameVec_[iparam] <<std::endl;
      for (int k = 0; k < solutionSize_; ++k)
      {
        Xyce::dout() << "dXdp[" << std::setw(3) << k << "] = "<< std::setw(15)<< std::scientific
          << std::setprecision(8)<< (*dXdp)[k]
          <<std::endl;
      }

      std::ostringstream filename; 
      filename << netlistFilename_ << "_dxdp";
      filename << std::setw(3) << std::setfill('0') << iparam;
      filename << ".txt";
      dXdp->writeToFile(const_cast<char *>(filename.str().c_str()));
    }
  } // end of param for loop

  // Now store the DQdx*dXdp matvec and the DFdx*dXdp matvec.
  // These are not needed for steady-state sensitivities, only for transient direct.
  bool Transpose = false;
  ds.nextDQdxDXdpPtrVector->putScalar(0.0);
  ds.dQdxMatrixPtr->matvec( Transpose , *ds.nextDXdpPtrVector, *ds.nextDQdxDXdpPtrVector );

  ds.nextDFdxDXdpPtrVector->putScalar(0.0);
  ds.dFdxMatrixPtr->matvec( Transpose , *ds.nextDXdpPtrVector, *ds.nextDFdxDXdpPtrVector );

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // Now get the final dOdp's (one for each objective/param combination).
  for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
  {
    for (iparam=0; iparam< numSensParams_; ++iparam)
    {
      double tmp = objFuncDataVec_[iobj]->dOdXVectorPtr->dotProduct( *Teuchos::rcp(dXdpPtrVector->getNonConstVectorView(iparam)) );
      tmp += objFuncDataVec_[iobj]->dOdp;

      ds.dOdpVec_.push_back(tmp);

      // get scaled value.  dO/dp*(p/100)
      double normalize = ds.paramOrigVals_[iparam]/100.0;
      tmp *= normalize;
      ds.scaled_dOdpVec_.push_back(tmp);
    }
  }

  if (stdOutputFlag_)
  {
    stdOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_, Xyce::lout() );
  }

  if (fileOutputFlag_)
  {
    fileOutput(std::string("Direct"), ds.paramOrigVals_, ds.dOdpVec_, ds.scaled_dOdpVec_);
  }

  return 1;
}


//-----------------------------------------------------------------------------
// Function      : Sensitivity::calcObjFuncDerivs
//
// Purpose       : Evaluates user-specified objective function expressions 
//                 and their derivatives with respect to the solution vector x.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/02
//-----------------------------------------------------------------------------
bool Sensitivity::calcObjFuncDerivs ()
{
  if ( !objFuncGIDsetup_ )
  {
    setupObjectiveFuncGIDs ( objFuncDataVec_, *(pdsMgrPtr_->getPDSComm()), top_, *outMgrPtr_ );
    objFuncGIDsetup_ = true;
  }

  return evaluateObjFuncs (
      objFuncDataVec_, 
      *(pdsMgrPtr_->getPDSComm()), 
      *nonlinearEquationLoader_, 
      *dsPtr_, 
      sec, 
      netlistFilename_); 
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::calcObjFuncTimeDerivs
//
// Purpose       : computes time derivatives of the objective functions.
//
//                 This is done if the delay calculations are requested.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/25/2019
//-----------------------------------------------------------------------------
bool Sensitivity::calcObjFuncTimeDerivs ()
{
  if(!timeDerivsSetup_)
  {
    for(int ii=0;ii<objFuncDataVec_.size();++ii)
    {
      objectiveFunctionData * ofDataPtr = new objectiveFunctionData();

      std::string tmpString;
      std::string & ofString = objFuncDataVec_[ii]->objFuncString;

      std::string::iterator iter;

      if ( ofString[0] == '{' && ofString[ofString.size()-1] == '}')
      {
        tmpString = ofString.substr(1, ofString.size()-2);
      }
      else
      {
        tmpString = objFuncDataVec_[ii]->objFuncString;
      }

      ofDataPtr->objFuncString = "{ddt(" + tmpString + ")}";
      objFuncTimeDerivDataVec_.push_back(ofDataPtr);
    }

    bool checkTimeDeriv=false;
    setupObjectiveFunctions(expressionGroup_, objFuncTimeDerivDataVec_, *outMgrPtr_, *lasSysPtr_, commandLine_, checkTimeDeriv);
    timeDerivsSetup_ = true;
  }

  if ( !objFuncTimeDerivGIDsetup_ )
  {
    setupObjectiveFuncGIDs (
        objFuncTimeDerivDataVec_,
        *(pdsMgrPtr_->getPDSComm()),
        top_, *outMgrPtr_ );

    objFuncTimeDerivGIDsetup_ = true;
  }

  return evaluateObjFuncs (
      objFuncTimeDerivDataVec_, 
      *(pdsMgrPtr_->getPDSComm()),
      *nonlinearEquationLoader_, 
      *dsPtr_, 
      sec, 
      netlistFilename_); 
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::solveAdjoint
// Purpose       : Solves first for the vector, lambda, of adjoint variables.
//                 Afterwards, it solves for dO/dp.
//
// Special Notes : This function is only for steady-state adjoints.  For 
//                 transient adjoints other functions are involved, particularly
//                 solveTransientAdjoint, below.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/20/02
//-----------------------------------------------------------------------------
int Sensitivity::solveAdjoint ()
{
  Stats::StatTop _solveSteadyAdjointStat("Solve Steady-State Adjoint");
  Stats::TimeBlock _solveSteadyAdjointTimer(_solveSteadyAdjointStat);

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "In Sensitivity::solveAdjoint" << std::endl;
  }

  TimeIntg::DataStore & ds = *dsPtr_;

  Linear::MultiVector * dfdpPtrVector_ = ds.nextDfdpPtrVector;
  Linear::MultiVector * dbdpPtrVector_ = ds.nextDbdpPtrVector;
  Linear::MultiVector * sensRHSPtrVector = ds.sensRHSPtrVector;

  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_),0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

  bool useTran = jacobianMatrixPtr_->useTranspose ();
  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    std::string matrixFile = netlistFilename_ + "_adjointMatrix.txt";
    jacobianMatrixPtr_->writeToFile(const_cast<char *>(matrixFile.c_str()));
  }

  for (int iobj=0;iobj<objFuncDataVec_.size();++iobj)
  {
    // copy the current dOdx vector into the rhs vector data structure.
    rhsVectorPtr_->update(1.0, *(objFuncDataVec_[iobj]->dOdXVectorPtr),0.0);

    jacobianMatrixPtr_->setUseTranspose (true);
    int status = lasSolverRCPtr_->solveTranspose(reuseFactors_);
    if (status!=0)
    {
      Report::DevelFatal().in("Sensitivity::solveAdjoint")
        << "Solver failed";
    }

    // allocate the dxdp vector for this param, and
    // copy the resulting deltax vector into the dxdp data structure.
    lambdaVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

    if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
    {
      std::string filename = netlistFilename_ + "_lambda.txt";
      lambdaVectorPtr_->writeToFile(const_cast<char *>(filename.c_str()));
    }

    // Now that we have lambda, get the dOdp's by doing dot products of
    // lambda * df/dp.

    // assuming that this is steady-state, then we need to sum together dfdp and -dbdp.
    sensRHSPtrVector->linearCombo(1.0,*dfdpPtrVector_,-1.0,*dbdpPtrVector_);

    // compute final dot products for all params.
    std::vector<double> tmp( numSensParams_ );
    sensRHSPtrVector->dotProduct( *lambdaVectorPtr_, tmp );

    // do the final dot products, one for each param.
    for (int iparam=0; iparam< numSensParams_; ++iparam)
    {
      ds.dOdpAdjVec_.push_back(-1.0*tmp[iparam]);

      // get scaled value.  dO/dp*(p/100)
      double normalize = ds.paramOrigVals_[iparam]/100.0;
      ds.scaled_dOdpAdjVec_.push_back(-1.0*tmp[iparam]*normalize);
    }
  }

  // restore the useTranspose flag to the original setting. (false)
  jacobianMatrixPtr_->setUseTranspose (useTran);

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // set the sensitivity information to the screen:
  if (stdOutputFlag_  && mode_ != Nonlinear::TRANSIENT)
  {
    stdOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_, Xyce::lout() );
  }
  if (fileOutputFlag_)
  {
    fileOutput(std::string("Adjoint"), ds.paramOrigVals_, ds.dOdpAdjVec_, ds.scaled_dOdpAdjVec_);
  }

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::solveTransientAdjoint
// Purpose       : Solves first for the vector, lambda, of adjoint variables.
//                 Afterwards, it solves for dO/dp.
// Special Notes : 
//
//  This function is called directly from N_ANP_Transient, via N_NLS_Manager, etc.
//  Unlike some of the other forms of sensitivity analysis, this function
//  is NOT called from the Sensitivity::solve function.
//
//  That function is called directly from external to the NLS package.  And
//  then, after that call, the ::solve function calls either ::solveDirect or
//  ::solveAdjoint (or both).  This is how it works for DCOP direct, DCOP 
//  adjoint, and transient direct.
//
//  For transient Adjoint I've decided not to do it that way, because 
//  this transient adjoint calculation is too different from the steady-state
//  direct/adjoint calculations, and also too different from the transient
//  direct calculation.  So, it needs to be its own thing.  As such, there
//  is some setup code copied from ::solve in this function.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/8/2016
//-----------------------------------------------------------------------------
int Sensitivity::solveTransientAdjoint (bool timePoint,
       std::vector<double> & objectiveVec,
       std::vector<double> & dOdpVec, 
       std::vector<double> & dOdpAdjVec,
       std::vector<double> & scaled_dOdpVec, 
       std::vector<double> & scaled_dOdpAdjVec)
{
  Stats::StatTop _solveTransientAdjointStat("Solve Transient Adjoint Step");
  Stats::TimeBlock _solveTransientAdjointTimer(_solveTransientAdjointStat);

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "In Sensitivity::solveTransientAdjoint" << std::endl;
  }

  //---------------------------------------------------------------------------
  nls_->enableSensitivity (); 

  //----------------------------------------------------------------------------
  TimeIntg::DataStore & ds = *dsPtr_;
  Linear::Vector & lambda = *(ds.nextLambdaPtr);
  Linear::Matrix & Jac = *(ds.JMatrixPtr);

  // first save a copy of the rhs vector, in case we want it later.
  savedRHSVectorPtr_->update(1.0, *(rhsVectorPtr_),0.0);
  savedNewtonVectorPtr_->update(1.0, *(NewtonVectorPtr_),0.0);

  bool useTran = jacobianMatrixPtr_->useTranspose ();

  bool voltageLimiterStatus = nonlinearEquationLoader_->getVoltageLimiterStatus();

  // for now, can only handle single objective function.
  if (objFuncDataVec_.size() != 1)
  {
    Report::DevelFatal().in("Sensitivity::solveTransientAdjoint")
      << "Transient Adjoints can only handle a single objective function.";
  }

  {
    nonlinearEquationLoader_->setVoltageLimiterStatus(false);

    if (newLowMem_)
    {
      // re-purpose these DAE matrix structures for a 2-step Jacobian history.
      Linear::Matrix & Jac0 = *ds.dQdxMatrixPtr; // one step back
      Linear::Matrix & Jac1 = *ds.dFdxMatrixPtr; // two steps back
      Linear::Matrix & JacTmp = *ds.tmpMatrixPtr; // temp

      // Put the dFdp stuff HERE:  We need to do 2 (maxOrder 1) or 3 (maxOrder 2) time points.
      //   Do in this maxOrder:   dFdp(t_{i-2}), dFdp(t_{i-1}), dFdp(t_i) 
      //   then do Jacobian(t_i).
      //
      int index=ds.itAdjointIndex;
      int maxOrder=2; // hardwired for now.  Always assume we need 2 previous points.  Won't cost much if wrong.
      int indexMin = std::max(index-maxOrder,0);

      // load the direct sensitivity residuals (dF/dp) and save.  
      if (timePoint)  // if at the first one of these, then compute dFdp (really dqdp) for all 3(2) points
      {
        for (int it=indexMin;it<=ds.itAdjointIndex;++it)
        { 
          // set proper history
          ds.updateSolDataArraysAdjoint(it);

          nonlinearEquationLoader_->loadRHS ();
          nonlinearEquationLoader_->loadJacobian ();

          if (it==ds.itAdjointIndex-2)
          {
            Jac1.put(0.0);
            Jac1.add(Jac);
          }
          else if (it==ds.itAdjointIndex-1)
          {
            Jac0.put(0.0);
            Jac0.add(Jac);
          }

          // get the device sensitivities (df/dp, dq/dp, and db/dp)
          loadSensitivityResiduals (difference_, forceFD_, forceDeviceFD_, forceAnalytic_, 
              sqrtEta_, netlistFilename_, 
            ds, *nonlinearEquationLoader_, paramNameVec_, getAnalysisManager());

          // save it.
          if (it<ds.itAdjointIndex)
          {
            nonlinearEquationLoader_->updateSensitivityHistoryAdjoint();
          }
        }
      }
      else // otherwise just compute the "new" one, which is dFdp(t_{i-maxOrder})
      {
        if (index-maxOrder >= 0)
        {
          // restore the history for the Jacobian
          ds.updateSolDataArraysAdjoint(index-maxOrder);
          nonlinearEquationLoader_->loadRHS ();
          nonlinearEquationLoader_->loadJacobian ();

          // kinda klunky.  Rotate the matrices.  N_LAS_Matrix = operator not implemented
          JacTmp.put(0.0); JacTmp.add(Jac);   // JacTmp = Jac
          Jac.put(0.0);    Jac.add(Jac0);     // Jac = Jac0
          Jac0.put(0.0);   Jac0.add(Jac1);    // Jac0 = Jac1
          Jac1.put(0.0);   Jac1.add(JacTmp);  // Jac1 = JacTmp

          // get the device sensitivities (df/dp, dq/dp, and db/dp)
          loadSensitivityResiduals (difference_, forceFD_, forceDeviceFD_, forceAnalytic_, 
              sqrtEta_, netlistFilename_, 
            ds, *nonlinearEquationLoader_, paramNameVec_, getAnalysisManager());

          // save it, updating the history.
          nonlinearEquationLoader_->updateSensitivityHistoryAdjoint2();
        }
      }

      // assemble the device sensitivities into a DAE
      nonlinearEquationLoader_->loadFunctionDerivativesForTranAdjoint ();

      // End of the dFdp stuff 
    }
    else
    {
      // via the loader, load the Jacobian matrix.  
      nonlinearEquationLoader_->loadRHS ();
      nonlinearEquationLoader_->loadJacobian ();
    }
    
    // zero out the RHS vector, as it currently contains the contents of 
    // a conventional RHS DAE load, which we don't want for an adjoint
    // computation.
    rhsVectorPtr_->putScalar(0.0);

    // via the loader, setup adjoint sensitivity residual.
    nonlinearEquationLoader_->loadAdjointSensitivityResidual ();

    // add the current dOdx vector into the rhs vector data structure.
    // For local sensitivities, which this function assists in computing,
    // The RHS vector should only be nonzero for the time point of interest 
    // for this adjoint calculation.
    if (timePoint)
    {
      calcObjFuncDerivs ();

#if 0
      if (computeDelays_) // doesn't yet work for adjoint
      {
        calcObjFuncTimeDerivs ();
      }
#endif

      objectiveVec.clear();
      ds.objectiveVec_.clear();

      // for the other types of sensitivity, this would normally be a loop,
      // as they support multiple objective functions.  Currently, transient
      // adjoint only supports a single objective.
      int iobj=0;
      objectiveVec.push_back(objFuncDataVec_[iobj]->objFuncEval);
      ds.objectiveVec_.push_back(objFuncDataVec_[iobj]->objFuncEval);

      rhsVectorPtr_->update(1.0, *(objFuncDataVec_[iobj]->dOdXVectorPtr),1.0);
    }

    bool reuseTheFactors = false; // don't use the class variable reuseFactors_! It is wrong here.
    jacobianMatrixPtr_->setUseTranspose (true);
    int status = lasSolverRCPtr_->solveTranspose(reuseTheFactors);
    if (status!=0)
    {
      Report::DevelFatal().in("Sensitivity::solveTransientAdjoint")
        << "Solver failed";
    }

    lambda.update(1.0, *(NewtonVectorPtr_),0.0);

#if 0
    dout() << "Lambda vector: " << std::endl;
    lambda.print(dout());
#endif

    // Now that we have lambda, get the dOdp's by doing dot products of
    // lambda * df/dp.  Do dot products, one for each param.  
    int index = ds.itAdjointIndex; 

    if (!newLowMem_)
    {
      // Need the "function" sensitivity (dfdp, but without the chain rule terms present in the 
      // direct sensitivity transient residual).  For now this is saved in a history
      // vector. (but in future should be computed in real time)
      std::vector<double> tmp( numSensParams_ );
      if (sparseAdjointStorage_)
      {
        ds.sparseFunctionSensitivityHistory[index]->dotProduct( lambda, tmp );
      }
      else
      {
        ds.functionSensitivityHistory[index]->dotProduct( lambda, tmp );
      }

      // Add in contributions into dOdpAdjVec_ 
      for (int iparam=0; iparam< numSensParams_; ++iparam)
      {
        ds.dOdpAdjVec_[iparam] += tmp[iparam];
#if 0
        dout() << "solveTransientAdjoint:   tmp["<<iparam<<"] = " << tmp[iparam] << std::endl;
#endif
      }
    }
    else
    {
      for (int iparam=0; iparam< numSensParams_; ++iparam)
      {
        Teuchos::RCP<Linear::Vector> functionSens = Teuchos::rcp( ds.sensRHSPtrVector->getNonConstVectorView(iparam) );
        ds.dOdpAdjVec_[iparam] += lambda.dotProduct(*functionSens);
      }
    }
  }

  // copy it over....
  dOdpAdjVec = ds.dOdpAdjVec_;

  // restore the useTranspose flag to the original setting. (false, probably)
  jacobianMatrixPtr_->setUseTranspose (useTran);

  // Restore the RHS and Newton vectors.
  rhsVectorPtr_->update(1.0, *(savedRHSVectorPtr_),0.0);
  NewtonVectorPtr_->update(1.0, *(savedNewtonVectorPtr_),0.0);

  // Restore voltlim status
  nonlinearEquationLoader_->setVoltageLimiterStatus( voltageLimiterStatus );

  return 1;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::setOptions
//
// Purpose       : This function processes the .SENS line
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/15/02
//-----------------------------------------------------------------------------
bool Sensitivity::setOptions(const Util::OptionBlock& OB)
{
  bool bsuccess = true;
  Util::ParamList::const_iterator iter = OB.begin();
  Util::ParamList::const_iterator end   = OB.end();

  numSensParams_ = 0;

  for ( ; iter != end; ++ iter)
  {
    if ( std::string( iter->uTag() ,0,7) == "OBJFUNC") // this is a vector
    {
      objectiveFunctionData * ofDataPtr = new objectiveFunctionData();
      ofDataPtr->objFuncString = iter->stringValue();
      objFuncDataVec_.push_back(ofDataPtr);
      objFuncGiven_ = true;
    }
    else if ( std::string( iter->uTag() ,0,7) == "OBJVARS") // this is a vector
    {
       // do nothing for now.  this is for AC
    }
    else if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      ExtendedString tag = iter->stringValue();
      tag.toUpper();
      // set up the initial skeleton of the maps:
      ++numSensParams_;
      paramNameVec_.push_back(tag);
      int sz = tag.size();
      if (sz > maxParamStringSize_)
      {
        maxParamStringSize_ = sz;
      }
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() 
        << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }
  
  if (numSensParams_ == 0)
  {
   Report::UserFatal0() << "No PARAM values found on .SENS line";
  }

  numObjectives_ = objFuncDataVec_.size();

  // parse the expression now, so if there are any errors, they will come
  // up early in the simulation.
  if (objFuncGiven_)
  {
    setupObjectiveFunctions(expressionGroup_, objFuncDataVec_, *outMgrPtr_, *lasSysPtr_, commandLine_);
  }

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    std::vector<std::string>::iterator iter;
    std::vector<std::string>::iterator begin = paramNameVec_.begin();
    std::vector<std::string>::iterator end   = paramNameVec_.end  ();

    for (iter=begin;iter!=end;++iter)
    {
      Xyce::dout() << *iter<<std::endl;
    }
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Sensitivity::setSensitivityOptions
// Purpose       : This function processes the .options SENSITIVITY line
// Special Notes : this is called after setOptions
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/20/2013
//-----------------------------------------------------------------------------
bool Sensitivity::setSensitivityOptions(const Util::OptionBlock &OB)
{
  bool bsuccess = true;
  Util::ParamList::const_iterator it  = OB.begin();
  Util::ParamList::const_iterator end = OB.end();
  for ( ; it != end; ++ it)
  {
    if ((*it).uTag() == "ADJOINT")
    {
      solveAdjointFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DIRECT")
    {
      solveDirectFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTSCALED")
    {
      outputScaledFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "OUTPUTUNSCALED")
    {
      outputUnscaledFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "STDOUTPUT")
    {
      stdOutputFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DIAGNOSTICFILE")
    {
      fileOutputFlag_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEFD")
    {
      forceFD_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEDEVICEFD")
    {
      forceDeviceFD_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "FORCEANALYTIC")
    {
      forceAnalytic_= 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "NEWLOWMEM")
    {
      newLowMem_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "SPARSESTORAGE")
    {
      sparseAdjointStorage_ = 
        static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "COMPUTEDELAYS")
    {
      computeDelays_ = static_cast<bool>((*it).getImmutableValue<bool>());
    }
    else if ((*it).uTag() == "DIFFERENCE")
    {
      ExtendedString sval=(*it).stringValue();
      sval.toUpper();
      if(sval=="FORWARD")
      {
        difference_=SENS_FWD;
      }
      else if(sval=="REVERSE")
      {
        difference_=SENS_REV;
      }
      else if(sval=="CENTRAL")
      {
        difference_=SENS_CNT;
        Report::UserFatal0() << "difference=central not supported.\n";
      }
      else
      {
        Report::UserFatal0() << "difference not recognized!\n";
      }
    }
    else if ((*it).uTag() == "SQRTETA")
    {
      sqrtEta_ = (*it).getImmutableValue<double>();
      sqrtEtaGiven_ = true;
    }
    else if ((*it).uTag() == "REUSEFACTORS")
    {
      reuseFactors_ = (*it).getImmutableValue<double>();
    }
    else if (DEBUG_NONLINEAR && (*it).uTag() == "DEBUGLEVEL")
    {
      debugLevel_ = (*it).getImmutableValue<int>();
      Xyce::setSensitivityDebugLevel(debugLevel_);
    }
    else if ( std::string( (*it).uTag() ,0,17) == "ADJOINTTIMEPOINTS") // this is a vector
    {
      // do nothing
    }

    else
    {
      Xyce::Report::UserWarning() << (*it).uTag() 
        << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }

  if (!sqrtEtaGiven_)
  {
    double epsilon = fabs(Util::MachineDependentParams::MachineEpsilon());
    sqrtEta_= std::sqrt(epsilon);
  }

  // Figure out if we need to force the entire matrix and RHS to be re-evaluated 
  // every time.  (the default behavior is to only re-evaluate the nonlinear parts)
  //
  // It is safe to do this here, as this function is called from the "enableSensitivity" 
  // function, which is called before the main analysis starts, but after everything is
  // set up.  So, we know enough to figure this out now.
  //
  // Also, note that the "setOptions" function was called before this one, so the
  // parameter vector is set up. (with the parameter names)

  TimeIntg::DataStore & ds = *dsPtr_;

  // call setupOriginalParams to test if these parameters are known.(basically just a setup test)
  bool origParamsDone =  setupOriginalParams (ds, *nonlinearEquationLoader_, paramNameVec_, getAnalysisManager());

  // determine if any parameters require numerical derivatives.  If so, (or if the user has requested "forceFD"
  // then we must force full linear system evaluation.
  bool allAnalyticalAvailable = testForAnalyticDerivatives ( *nonlinearEquationLoader_, paramNameVec_, getAnalysisManager());

  if (forceFD_ || !allAnalyticalAvailable) 
  { 
    if (forceFD_)
    {
      Report::UserWarning0() << "User set FORCEFD=true, so numerical derivatives will be used for sensitivity calculations.  This also means that the full linear system must be evaluated at every step.  The default behavior is to only re-evaluate the nonlinear portions of the problem for efficiency.  The simulation may be slower as a result.";
    }
    else if (!allAnalyticalAvailable)
    {
      Report::UserWarning0() << "At least one specified sensitivity parameter lacks an analytic derivative, so numerical derivatives will be used for that parameter for sensitivity calculations.  This also means that the full linear system must be evaluated at every step.  The default behavior is to only re-evaluate the nonlinear portions of the problem for efficiency.  The simulation may be slower as a result.";
    }

    nonlinearEquationLoader_->setSeparateLoadFlag (false); 
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::doAllocations
//
// Purpose       : This is here only to call the allocateSensitivityArrays function
//                 in the data store class.  It is called from the 
//                 Xyce::Nonlinear::Manager::setupSensitivity function.
//
// Special Notes : 
//
// I wanted to just call this 1x, if possible.  It used to be called from setOptions.  
// This seemed safe b/c setOptions is always called if sensitivity analysis is used.
// But, I wanted to have the allocations be specific to direct and/or adjoint methods,
// and the setOptions function doesn't yet know what the user requested.
//
// I tried calling it from the setSensitivityOptions function, which is where the
// user sets if they want direct, adjoint or both.  This function is called
// second, after setOptions.  However, it is only called optionally.  If the user 
// just relies on defaults, then this fails.
//
// I could have called it from outside the sensitivity class, but that seemed more messy. 
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/1/2019
//-----------------------------------------------------------------------------
bool Sensitivity::doAllocations()
{
  dsPtr_->allocateSensitivityArrays(numSensParams_, solveDirectFlag_, solveAdjointFlag_);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : extractSENSData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/11/2007
//-----------------------------------------------------------------------------
bool extractSENSData(
  IO::PkgOptionsMgr &           options_manager,
  IO::CircuitBlock &            circuit_block,
  const std::string &           netlist_filename,
  const IO::TokenVector &       parsed_line)
{
  Util::OptionBlock option_block("SENS", Util::OptionBlock::ALLOW_EXPRESSIONS, netlist_filename, parsed_line[0].lineNumber_);

  int numFields = parsed_line.size();

  int parameterStartPos = 1;

  // Create an option block to temporarily store the default options.
  Util::OptionBlock defaultOptions;

  // Get the default options from metadata.
  addDefaultOptionsParameters(options_manager, defaultOptions, "SENS" );

  // Extract the parameters from parsed_line.
  int parameterEndPos = numFields - 1;
  Util::ParamList inputParameters;
  Util::Param parameter("", "");
  int intervalParameterStart = -1;
  int i = parameterStartPos;
  std::string paramBaseName;
  while (i <= parameterEndPos-1)
  {
    // Check for equal sign.
    if ( parsed_line[i+1].string_ != "=" )
    {
      // Stop after the tagged parameters have been extracted
      // from a .OPTIONS RESTART or .OPTIONS OUTPUT line, they
      // will be handled later.
      intervalParameterStart = i;
      break;
    }

    // Extract parameter name and value from parsed_line and add to
    // parameter list. Check to see if the parameter is "VECTOR"
    // valued and treat accordingly.
    parameter.set( parsed_line[i].string_, "" );
    Util::Param *parameterPtr = Util::findParameter(defaultOptions.begin(), defaultOptions.end(), parameter.tag());
    if (parameterPtr == NULL)
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "No .SENS parameter " << parameter.tag() << " found, parameter will be ignored.";
      i+= 3;
    }
    else if (parameterPtr->stringValue() != "VECTOR")
    {
      // error out if the next push_back will cause a segfault.
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".SENS line is missing one or more required parameters";
        return false;
      }
      parameter.setVal( parsed_line[i+2].string_ );
      inputParameters.push_back( parameter );
      i+= 3;

      if (i < parameterEndPos-1 && parsed_line[i].string_ == ",")
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".SENS parameter " << parameter.tag() << " is not a VECTOR, but has comma in value.";
      }
    }
    else
    {
      // We have a vector valued parameter.
      // Name the jth component of the parameter of the vector by appending
      // "j" to the parameter name.
      std::ostringstream paramName;
      std::string paramBaseName = ExtendedString(parsed_line[i].string_).toUpper();
      int j = 1;

      // used to help stop reading off the end of parsed_line
      int testSize= parsed_line.size() - 1;

      paramName << paramBaseName << j;
      i += 2;

      if ( i > testSize)
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".SENS line is missing one or more elements in a vector parameter";
        return false;
      }
      else if (parsed_line[i].string_ == ",")
      {
        //catch the invalid case of param=,value,value
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
            << ".SENS line is mis-formatted. Vector parameter begins with a comma";
        return false;
      }

      parameter.set(paramName.str(), parsed_line[i].string_);
      option_block.addParam(parameter);

      // The nominal format is param=value,value,value   
      // However, this loop should error out on invalid lines that have:
      // errors like:
      //
      //        objfunc=value,value,value,
      //        param=value,value,,value
      //
      // or basically any case where the "parameter value" is a comma,
      // which indicates that something went wrong during parsing.
      while ((i < testSize) && (parsed_line[i+1].string_ == ",") )
      {
        paramName.str("");
        ++j;
        paramName << paramBaseName << j;
        i += 2;
        if ( (i <= testSize) && (parsed_line[i].string_ != ",") )
	{ 
          parameter.set(paramName.str(), parsed_line[i].string_);
          option_block.addParam(parameter);
        }
        else
	{
          Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
            << ".SENS line is mis-formatted. Error parsing vector parameter";
          return false;
        }
      }

      ++i;

      // Also need to guard against these cases:
      //
      //     objfunc=value,value value
      //     objfunc=value,value value param=value 
      if ( (i==testSize) || (( i+1 <= testSize) && (parsed_line[i+1].string_ != "=")) )
      {
        Report::UserError0().at(netlist_filename, parsed_line[0].lineNumber_)
          << ".SENS line is mis-formatted. Possibly a missing comma in vector parameter";
        return false;
      }
    }
  }

  // For each input parameter, check that it is in the default
  // set and if so, set its value in "parameters" to the input
  // value, otherwise flag it as an unknown parameter.
  for (Util::ParamList::const_iterator it = inputParameters.begin(), end = inputParameters.end(); it != end; ++it)
  {
    Util::Param *parameterPtr = Util::findParameter(defaultOptions.begin(), defaultOptions.end(), (*it).tag());
    if ( parameterPtr != NULL )
    {
      parameterPtr->setVal(*it);
      option_block.addParam( *parameterPtr );
    }
    else
    {
      Report::UserWarning0().at(netlist_filename, parsed_line[0].lineNumber_)
        << "No .SENS parameter " << (*it).tag() << " found, parameter will be ignored.";
    }
  }

  circuit_block.addOptions(option_block);

  return true; // Only get here on success.
}

//-----------------------------------------------------------------------------
// Function      : Sensitivity::populateMetadata
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Baur
// Creation Date :
//-----------------------------------------------------------------------------
void Sensitivity::populateMetadata(
  IO::PkgOptionsMgr &   options_manager)
{
  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("SENSITIVITY");

    parameters.insert(Util::ParamMap::value_type("DEBUGLEVEL", Util::Param("DEBUGLEVEL", 0)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTLAMBDA", Util::Param("OUTPUTLAMBDA", 0)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTTRANSIENTADJOINT", Util::Param("OUTPUTTRANSIENTADJOINT", 0)));
    parameters.insert(Util::ParamMap::value_type("FULLADJOINTTIMERANGE", Util::Param("FULLADJOINTTIMERANGE", 0)));
    parameters.insert(Util::ParamMap::value_type("ADJOINT", Util::Param("ADJOINT", 0)));
    parameters.insert(Util::ParamMap::value_type("DIRECT", Util::Param("DIRECT", 0)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTSCALED", Util::Param("OUTPUTSCALED", 0)));
    parameters.insert(Util::ParamMap::value_type("OUTPUTUNSCALED", Util::Param("OUTPUTUNSCALED", 1)));
    parameters.insert(Util::ParamMap::value_type("STDOUTPUT", Util::Param("STDOUTPUT", 0)));
    parameters.insert(Util::ParamMap::value_type("DIAGNOSTICFILE", Util::Param("DIAGNOSTICFILE", 0)));
    parameters.insert(Util::ParamMap::value_type("DIFFERENCE", Util::Param("DIFFERENCE", 0)));
    parameters.insert(Util::ParamMap::value_type("SQRTETA", Util::Param("SQRTETA", 1.0e-8)));
    parameters.insert(Util::ParamMap::value_type("REUSEFACTORS", Util::Param("REUSEFACTORS",1)));
    parameters.insert(Util::ParamMap::value_type("ADJOINTBEGINTIME", Util::Param("ADJOINTBEGINTIME", 0.0)));
    parameters.insert(Util::ParamMap::value_type("ADJOINTFINALTIME", Util::Param("ADJOINTFINALTIME", 1.0e+199)));
    parameters.insert(Util::ParamMap::value_type("FORCEFD", Util::Param("FORCEFD", 0)));
    parameters.insert(Util::ParamMap::value_type("FORCEDEVICEFD", Util::Param("FORCEDEVICEFD", 0)));
    parameters.insert(Util::ParamMap::value_type("FORCEANALYTIC", Util::Param("FORCEANALYTIC", 0)));
    parameters.insert(Util::ParamMap::value_type("NEWLOWMEM", Util::Param("NEWLOWMEM", 1)));
    parameters.insert(Util::ParamMap::value_type("SPARSESTORAGE", Util::Param("SPARSESTORAGE", 1)));
    parameters.insert(Util::ParamMap::value_type("COMPUTEDELAYS", Util::Param("COMPUTEDELAYS", 1)));
    parameters.insert(Util::ParamMap::value_type("ADJOINTTIMEPOINTS", Util::Param("ADJOINTTIMEPOINTS", "VECTOR")));
  }

  {
    Util::ParamMap &parameters = options_manager.addOptionsMetadataMap("SENS");

    parameters.insert(Util::ParamMap::value_type("OBJFUNC", Util::Param("OBJFUNC", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("OBJVARS", Util::Param("OBJVARS", "VECTOR")));
    parameters.insert(Util::ParamMap::value_type("PARAM", Util::Param("PARAM", "VECTOR")));
  }

  options_manager.addCommandParser(".SENS", extractSENSData);
}

} // namespace Nonlinear
} // namespace Xyce

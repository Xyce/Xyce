//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        : This file contains a helper class for sensitivity 
//                  analysis objective functions
//
// Special Notes  : 
//
// Creator        : Eric Keiter, SNL
//
// Creation Date  : 12/16/20
//
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_ObjectiveFunctions_h
#define Xyce_N_NLS_ObjectiveFunctions_h

#include<vector>

#include <N_UTL_fwd.h>
#include <N_PDS_fwd.h>
#include <N_TOP_fwd.h>
#include <N_UTL_MachDepParams.h>
#include <N_NLS_NonLinearSolver.h>
#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_Diagnostic.h>

#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>

#include <N_LAS_System.h>
#include <N_TOP_Topology.h>

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Class         : objective function data
// Purpose       :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/17/2015
//-----------------------------------------------------------------------------
template <typename ScalarT>
  class objectiveFunctionData
{
public:  
  objectiveFunctionData() :
    numExpVars(0),
    expVal(0.0),
    objFuncString(""),
    expPtr(0),
    numDdt(0),
    objFuncEval(0.0),
    dOdp(0.0)
  {};

public:
  int            numExpVars;  // size of the numVarDerivs array
  std::vector<std::string> expVarNames; // need this to get the GIDs
  std::vector<int>    expVarGIDs;  // need this even with new expression lib, to populate the dOdX vector

  std::vector<ScalarT> expVarDerivs; // this is returned by expPtr->evaluate
  ScalarT expVal; // this is returned by expPtr->evaluate

  std::string objFuncString;

  Util::Expression * expPtr;

  int numDdt;

  ScalarT objFuncEval;// value of the evaluated objective function.
  ScalarT dOdp;

  Linear::Vector* dOdXVectorRealPtr; // size of solution vector.
  Linear::Vector* dOdXVectorImagPtr; // size of solution vector.
};

//-----------------------------------------------------------------------------
template <typename ScalarT>
inline bool evaluateObjFuncs ( 
    std::vector<objectiveFunctionData<ScalarT> *> & objVec, 
    Parallel::Machine & comm,
    Loader::NonlinearEquationLoader & nlEquLoader_,
    std::string & netlistFilename)
{
  bool bsuccess = true;
  int i;

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->dOdXVectorRealPtr->putScalar(0.0);
  }

  // obtain the expression variable values.  It will only grab this value if it is owned on this processor.
  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expVarDerivs.resize (objVec[iobj]->numExpVars, 0.0);
  }

  //get expression value and partial derivatives
  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expPtr->processSuccessfulTimeStep();
    objVec[iobj]->expPtr->evaluate( 
        objVec[iobj]->expVal, 
        objVec[iobj]->expVarDerivs); 
    objVec[iobj]->expPtr->clearOldResult();

    objVec[iobj]->objFuncEval = objVec[iobj]->expVal;
    objVec[iobj]->dOdXVectorRealPtr->putScalar(0.0);
    for (i=0;i<objVec[iobj]->numExpVars;++i)
    {
      int tmpGID = objVec[iobj]->expVarGIDs[i];
      ScalarT tmpDODX = objVec[iobj]->expVarDerivs[i];

      if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
      {
        Xyce::dout() 
          <<  objVec[iobj]->expVarNames[i] << "  "
          << "i="<<i<<"  gid = " << tmpGID << "  dodx = "<< tmpDODX << std::endl;
      }

      if (tmpGID >= 0)
      {
        objVec[iobj]->dOdXVectorRealPtr->setElementByGlobalIndex(tmpGID, std::real(tmpDODX), 0);
      }
    }

    // Assuming this is zero:
    objVec[iobj]->dOdp = 0.0;
  }

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->dOdXVectorRealPtr->fillComplete();

    if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
    {
      std::string filename = netlistFilename + "_dodx.txt";
      objVec[iobj]->dOdXVectorRealPtr->writeToFile(const_cast<char *>(filename.c_str()));
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// complex version
template <>
inline bool evaluateObjFuncs ( 
    std::vector<objectiveFunctionData<std::complex<double> > *> & objVec, 
    Parallel::Machine & comm,
    Loader::NonlinearEquationLoader & nlEquLoader_,
    std::string & netlistFilename)
{
  bool bsuccess = true;
  int i;

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->dOdXVectorRealPtr->putScalar(0.0);
    objVec[iobj]->dOdXVectorImagPtr->putScalar(0.0);
  }

  // obtain the expression variable values.  It will only grab this value if it is owned on this processor.
  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expVarDerivs.resize (objVec[iobj]->numExpVars, 0.0);
  }

  //get expression value and partial derivatives
  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->expPtr->evaluate( 
        objVec[iobj]->expVal, 
        objVec[iobj]->expVarDerivs); 

    objVec[iobj]->objFuncEval = objVec[iobj]->expVal;
    objVec[iobj]->dOdXVectorRealPtr->putScalar(0.0);
    objVec[iobj]->dOdXVectorImagPtr->putScalar(0.0);
    for (i=0;i<objVec[iobj]->numExpVars;++i)
    {
      int tmpGID = objVec[iobj]->expVarGIDs[i];
      std::complex<double> tmpDODX = objVec[iobj]->expVarDerivs[i];

      if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
      {
        Xyce::dout() 
          <<  objVec[iobj]->expVarNames[i] << "  "
          << "i="<<i<<"  gid = " << tmpGID << "  dodx = "<< tmpDODX << std::endl;
      }

      if (tmpGID >= 0)
      {
        // For the purpose of AC sensitivities (which is being developed for)
        // what is needed is the partial derivatives df/dx and df/dy, where
        // The objective function is f(z), where f and z are complex-valued.
        //        z = x + iy
        // The expression library returns df/dz.  df/dx and df/dy are given by:
        //        df/dx = df/dz * dz/dx 
        //        df/dy = df/dz * dz/dy 
        // where dz/dx and dz/dy are 1.0 and 1.0i, respectively. 
        objVec[iobj]->dOdXVectorRealPtr->setElementByGlobalIndex(tmpGID, std::real(tmpDODX), 0);
        objVec[iobj]->dOdXVectorImagPtr->setElementByGlobalIndex(tmpGID, std::real(tmpDODX), 0); 
      }
    }

    // Assuming this is zero:
    objVec[iobj]->dOdp = 0.0;
  }

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    objVec[iobj]->dOdXVectorRealPtr->fillComplete();
    objVec[iobj]->dOdXVectorImagPtr->fillComplete();

    if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
    {
      std::string filename = netlistFilename + "_dodx.txt";
      objVec[iobj]->dOdXVectorRealPtr->writeToFile(const_cast<char *>(filename.c_str()));
      objVec[iobj]->dOdXVectorImagPtr->writeToFile(const_cast<char *>(filename.c_str()));
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : setupObjectiveFunctions
// Purpose       : Allocate, parse and resolve objective func expressions
//
// Special Notes : templated so it can be used by both frequency domain 
//                 and time domain sensitivities.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/16/20
//-----------------------------------------------------------------------------
  template <typename ScalarT>
inline void setupObjectiveFunctions (
    Teuchos::RCP<Xyce::Util::baseExpressionGroup> & exprGroup,
    std::vector<objectiveFunctionData<ScalarT> *> & objVec,
    IO::OutputMgr & output_manager, Linear::System & lasSys,
    const IO::CmdParse &cp,
    bool checkTimeDeriv=true)
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
    const std::vector<std::string> & strings = objVec[iobj]->expPtr->getUnresolvedParams();

    const Util::ParamMap & context_param_map = output_manager.getMainContextParamMap();
    const Util::ParamMap & context_global_param_map = output_manager.getMainContextGlobalParamMap();
    for (int istring=0;istring<strings.size();istring++)
    {
      Util::ParamMap::const_iterator param_it = context_param_map.find(strings[istring]);

      if (param_it != context_param_map.end())
      {
        const Util::Param &replacement_param = param_it->second;

        if ( replacement_param.getType() == Xyce::Util::STR ||
             replacement_param.getType() == Xyce::Util::DBLE )
        {
          enumParamType paramType=DOT_PARAM;
          if (!objVec[iobj]->expPtr->make_constant(strings[istring], replacement_param.getImmutableValue<double>(),paramType))
          {
            Report::UserWarning0() << "Problem converting parameter " << strings[istring] << " to its value.";
          }
        }
        else if (replacement_param.getType() == Xyce::Util::EXPR)
        {
          enumParamType paramType=DOT_PARAM;
          objVec[iobj]->expPtr->attachParameterNode (strings[istring], replacement_param.getValue<Util::Expression>(),paramType);
        }
      }
      else
      {
        // if this string is found in the global parameter map, then attach it to the expression
        param_it = context_global_param_map.find(strings[istring]);
        if (param_it != context_global_param_map.end())
        {
          const Util::Param &replacement_param = param_it->second;
          if (replacement_param.getType() == Xyce::Util::EXPR)
          {
            const Util::Expression & expToBeAttached = replacement_param.getValue<Util::Expression>();
            objVec[iobj]->expPtr->attachParameterNode(strings[istring], expToBeAttached);
          }
        }
        else
        {
          Report::UserWarning0() << "This field: " << strings[istring] 
            << " from the objective function " << objVec[iobj]->objFuncString << " is not resolvable";
        }
      }
    }

    objVec[iobj]->numExpVars = objVec[iobj]->expVarNames.size();
    if (objVec[iobj]->numExpVars<=0)
    {
      Report::UserFatal0()
        <<  "Objective function does not contain a resolvable solution variable.";
    }

    objVec[iobj]->dOdXVectorRealPtr = lasSys.builder().createVector();
    objVec[iobj]->dOdXVectorImagPtr = lasSys.builder().createVector();
  }
}

//-----------------------------------------------------------------------------
// Function      : setupObjectiveFuncGIDs
// Purpose       : setup the expVarGIDs vector object.
// Special Notes : expVarNames was set up in setupObjectiveFunctions.  
//                 It should ONLY contain voltage and current variables.  
//                 It should not contain any .param or .global_params.
//                 With the old expression library, that would not have been 
//                 the case, but the newer one is better than that.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/16/20
//-----------------------------------------------------------------------------
template <typename ScalarT>
inline void setupObjectiveFuncGIDs (std::vector<objectiveFunctionData<ScalarT> *> & objVec, 
    Parallel::Machine & comm,
    Topo::Topology & top, IO::OutputMgr & output_manager)
{
  int found(0);
  int found2(0);
  int foundAliasNode(0);
  bool foundLocal(false);
  bool foundLocal2(false);
  bool foundAliasNodeLocal(false);

  for (int iobj=0;iobj<objVec.size();++iobj)
  {
    // set up the gid's:
    objVec[iobj]->expVarGIDs.resize( objVec[iobj]->numExpVars, -1);

    for (int i = 0; i < objVec[iobj]->numExpVars; ++i)
    {
      std::vector<int> svGIDList1, dummyList;
      char type1;

      // look for this variable as a node first.
      foundLocal = top.getNodeSVarGIDs(NodeID(objVec[iobj]->expVarNames[i], Xyce::_VNODE), svGIDList1, dummyList, type1);
      found = static_cast<int>(foundLocal);
      Xyce::Parallel::AllReduce(comm, MPI_LOR, &found, 1);

      // if looking for this as a voltage node failed, try a "device" (i.e. current) node.  I(Vsrc)
      foundLocal2 = false;
      if (!found)
      {
        foundLocal2 = top.getNodeSVarGIDs(NodeID(objVec[iobj]->expVarNames[i], Xyce::_DNODE), svGIDList1, dummyList, type1);
      }
      found2 = static_cast<int>(foundLocal2);
      Xyce::Parallel::AllReduce(comm, MPI_LOR, &found2, 1);

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
      if (!found && !found2)
      {
        IO::AliasNodeMap::const_iterator alias_it = output_manager.getAliasNodeMap().find(objVec[iobj]->expVarNames[i]);
        if (alias_it != output_manager.getAliasNodeMap().end())
        {      
          foundAliasNodeLocal = top.getNodeSVarGIDs(NodeID((*alias_it).second, Xyce::_VNODE), svGIDList1, dummyList, type1);
        }
      }
      foundAliasNode = static_cast<int>(foundAliasNodeLocal);
      Xyce::Parallel::AllReduce(comm, MPI_LOR, &foundAliasNode, 1);

      if (!found && !found2 && !foundAliasNode)
      {
        Report::UserFatal() << "objective function variable not found!  Cannot find " << objVec[iobj]->expVarNames[i] ;
      }

      if (found || found2 || foundAliasNode)
      {
        int tmpGID=-1;
        if(svGIDList1.size()==1) { tmpGID = svGIDList1.front(); }
        objVec[iobj]->expVarGIDs[i] = tmpGID;
      }
    }
  }
}

//-----------------------------------------------------------------------------
template <typename ScalarT>
inline void applyHocevarDelayTerms(
    std::vector<objectiveFunctionData<ScalarT> *> & objVec,
    std::vector<objectiveFunctionData<ScalarT> *> & objTimeDerivVec,
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

} // namespace Nonlinear
} // namespace Xyce

#endif // Xyce_N_NLS_Sensitivity_h


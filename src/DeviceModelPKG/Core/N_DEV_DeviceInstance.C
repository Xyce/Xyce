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
// Purpose        : Implementation of the base device instance class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/30/00
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_DeviceInstance.h>

#include <N_DEV_Const.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_NumericalJacobian.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace Device {

namespace {
static const std::vector<std::string> emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::DeviceInstance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceInstance::DeviceInstance(
  const InstanceBlock &         instance_block,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceEntity(parametric_data, factory_block.solverState_, factory_block.deviceOptions_, instance_block.getNetlistLocation()),
    name_(instance_block.getInstanceName()), 
    mlData(factory_block.matrixLoadData_),
    extData(factory_block.externData_),
    configuredForLeadCurrent(false),
    cols(factory_block.matrixLoadData_.cols),
    vals(factory_block.matrixLoadData_.vals),
    numJacPtr(NULL),
    origFlag(true),
    numIntVars(0),
    numExtVars(2),
    numStateVars(0),
    numStoreVars(0),
    loadLeadCurrent(false),
    numBranchDataVars(0),
    numBranchDataVarsIfAllocated(0),
    mergeRowColChecked(false)
{
  devConMap.resize(2);
  devConMap[0] = 1;
  devConMap[1] = 1;
  numJacPtr = new NumericalJacobian(factory_block.matrixLoadData_, factory_block.solverState_, factory_block.externData_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::~DeviceInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceInstance::~DeviceInstance()
{
  delete numJacPtr;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::enableLeadCurrentCalc
//
// Purpose       : This function configures a device to an auxiliary F & Q vector
//                 so that lead currents can be calculated for this device.c.
//
// Special Notes : This must be called soon after the constructor call
//                 before the store vector is allocated.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/22/13
//-----------------------------------------------------------------------------
void DeviceInstance::enableLeadCurrentCalc()
{
  if (!configuredForLeadCurrent)
  {
    // indicated that this device is now configured for lead current calculation
    // this avoids claiming too much space in the store vector if this function
    // is called more than once for a device.
    configuredForLeadCurrent = true;

    // set device instance flag to indicate the need to load lead current
    // data into store F & Q vectors
    loadLeadCurrent = true;

    // migrating the store vector calculation to the branch data vectors
    numBranchDataVars = numBranchDataVars + numBranchDataVarsIfAllocated;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceInstance::getDepSolnVars
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
// Special Notes : 
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & DeviceInstance::getDepSolnVars()
{
  return expVarNames;
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceInstance::getDepSolnTypes
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/15/20
//-----------------------------------------------------------------------------
const std::vector<int> & DeviceInstance::getDepSolnTypes()
{
  return expVarTypes;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepSolnLIDs
// Purpose       : Allows registration of LIDs of nodes and instances that
//                 appear in expressions that occur in the device.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/15/05
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepSolnLIDs(
  const std::vector< std::vector<int> > &       depSolnLIDVecRef)
{
  int size = expVarLIDs.size();
  if (size != depSolnLIDVecRef.size())
  {
    DevelFatal0(*this).in("DeviceInstance::registerDepSolnLIDs")
      << "Inconsistent number of LIDs returned from topology";
  }
  for (int i = 0; i < size; ++i)
  {
    if (depSolnLIDVecRef[i].size() != 1)
    {
      UserError0(*this) << "Problem with value for " << expVarNames[i]
                        << ".  This may be an incorrect usage of a lead current in place of a current through a voltage source.";
    }
    expVarLIDs[i] = depSolnLIDVecRef[i][0];
  }

  DeviceEntity::applyDepSolnLIDs(); // new for newExpression
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepSolnGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepSolnGIDs(
  const std::vector< std::vector<int> > & varList )
{
  int size = expVarGIDs.size();
  for (int i = 0; i < size; ++i)
  {
    expVarGIDs[i] = varList[i][0];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceInstance::registerJacLIDs( const JacobianStamp & jacLIDVec )
{
  if (getDeviceOptions().matrixSensitivityFlag || getDeviceOptions().testJacobianFlag)
  {
    devJacLIDs = jacLIDVec;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::testDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool DeviceInstance::testDAEMatrices(const std::vector<const std::string *> & nameVec)
{
  // if necessary, consolodate the LIDs vector.
  consolidateDevLIDs();
  return numJacPtr->testDAEMatrices(*this, nameVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::trivialStampLoader
// Purpose       : This function contains most of the original
//                 loadTrivialMatrixStamp function.  See comments above.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::trivialStampLoader (Linear::Matrix * matPtr)
{
  std::vector<int>::const_iterator firstVar;
  std::vector<int>::const_iterator lastVar;
  std::vector<int>::const_iterator iterVar;

  int localRows = matPtr->getLocalNumRows();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl
                 << "Loading trivial stamp for " << getName() << std::endl;
  }

  if (cols.size() < 1) cols.resize(1);
  if (vals.size() < 1) vals.resize(1);

  for (int i = 0; i < 2; ++i)
  {
    // do external vars first.
    if (i==0)
    {
      firstVar = extLIDVec.begin ();
      lastVar  = extLIDVec.end ();
    }
    // then do internal vars.
    else
    {
      firstVar = intLIDVec.begin ();
      lastVar  = intLIDVec.end ();
    }

    for (iterVar=firstVar; iterVar!=lastVar; ++iterVar)
    {
      int row = *iterVar;

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
      {
        Xyce::dout() << "matrix row = " << row;
      }

      if (row < 0)
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << ": NOT loading this one" << std::endl;
        }
        continue;
      }
      else
      {
        if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
        {
          Xyce::dout() << ": loading this one" << std::endl;
        }

        int count = 1;
        vals[0] = 1.0;
        cols[0] = row;

        matPtr->replaceLocalRow(row, count, &vals[0], &cols[0]);
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::loadTrivialDAE_FMatrixStamp
// Purpose       : See loadTrivialMatrixStamp - this is the same thing, except
//                 for the new-DAE F-Matrix.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::loadTrivialDAE_FMatrixStamp ()
{
  return trivialStampLoader (extData.dFdxMatrixPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::enablePDEContinuation(int &max_PDE_continuation_steps)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::disablePDEContinuation()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setPDEContinuationAlpha
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
void DeviceInstance::setPDEContinuationAlpha (double alpha)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setPDEContinuationBeta
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/04
//-----------------------------------------------------------------------------
void DeviceInstance::setPDEContinuationBeta (double beta)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::setInitialGuess ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
double DeviceInstance::getMaxTimeStepSize  ()
{
  return getDeviceOptions().defaultMaxTimeStep;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getNumericalSensitivity
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool DeviceInstance::getNumericalSensitivity ( const std::string & paramName,
                                std::vector<double> & dfdpVec,
                                std::vector<double> & dqdpVec,
                                std::vector<double> & dbdpVec,
                                std::vector<int> & FindicesVec,
                                std::vector<int> & QindicesVec,
                                std::vector<int> & BindicesVec )
{
  double origParamValue = 0.0;
  bool found = getParam(paramName, origParamValue); // const?

  if (found)
  {
    // if necessary, consolidate the LIDs vector.
    consolidateDevLIDs();

    if ( !(devLIDs.empty()))
    {
      if (!(FindicesVec.empty())) {  FindicesVec.clear(); }
      if (!(QindicesVec.empty())) {  QindicesVec.clear(); }
      if (!(BindicesVec.empty())) {  BindicesVec.clear(); }


      unordered_map<int,int> lidCheckMap;
      unordered_map<int,int> stateLidCheckMap;

      for (int i=0;i<devLIDs.size();++i)
      {
        int newID = devLIDs[i];
        if (lidCheckMap.find(newID) == lidCheckMap.end())
        {
          FindicesVec.push_back(newID);
          lidCheckMap[newID] = i;
        }
      }

      std::vector<int> stateIndicesVec;
      const IdVector & stateDevLIDs = getStaLIDVec();
      for (int i=0;i<stateDevLIDs.size();++i)
      {
        int newID = stateDevLIDs[i];
        if (stateLidCheckMap.find(newID) == stateLidCheckMap.end())
        {
          stateIndicesVec.push_back(newID);
          stateLidCheckMap[newID] = i;
        }
      } 

      QindicesVec.insert(QindicesVec.end(), FindicesVec.begin(), FindicesVec.end()); 
      BindicesVec.insert(BindicesVec.end(), FindicesVec.begin(), FindicesVec.end()); 

      Linear::Vector & Fvec          = (*extData.daeFVectorPtr);
      Linear::Vector & Qvec          = (*extData.daeQVectorPtr);
      Linear::Vector & Bvec          = (*extData.daeBVectorPtr);

      Linear::Vector & currSol      = (*extData.currSolVectorPtr);
      Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

      Linear::Vector & lastSta      = (*extData.lastStaVectorPtr);
      Linear::Vector & currSta      = (*extData.currStaVectorPtr);
      Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
      Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);

      int solSize=FindicesVec.size();
      int stateSize=stateIndicesVec.size();

      mlData.resizeSolnSizedVectors (solSize);
      mlData.resizeStateSizedVectors (stateSize);

      std::vector<double> & saveF = mlData.saveF;
      std::vector<double> & pertF = mlData.pertF;
      std::vector<double> & origF = mlData.origF;

      std::vector<double> & saveQ = mlData.saveQ;
      std::vector<double> & pertQ = mlData.pertQ;
      std::vector<double> & origQ = mlData.origQ;

      std::vector<double> & saveB = mlData.saveB;
      std::vector<double> & pertB = mlData.pertB;
      std::vector<double> & origB = mlData.origB;

      std::vector<double> & saveLastState = mlData.saveLastState;
      std::vector<double> & saveCurrState = mlData.saveCurrState;
      std::vector<double> & saveNextState = mlData.saveNextState;
      std::vector<double> & saveStateDerivs = mlData.saveStateDerivs;

      saveF.assign(saveF.size(),0.0);
      pertF.assign(pertF.size(),0.0);
      origF.assign(origF.size(),0.0);

      saveQ.assign(saveQ.size(),0.0);
      pertQ.assign(pertQ.size(),0.0);
      origQ.assign(origQ.size(),0.0);

      saveB.assign(saveB.size(),0.0);
      pertB.assign(pertB.size(),0.0);
      origB.assign(origB.size(),0.0);

      dfdpVec.resize(saveF.size(),0.0);
      dqdpVec.resize(saveQ.size(),0.0);
      dbdpVec.resize(saveB.size(),0.0);

      // save RHS information
      bool origFlag = getOrigFlag();

      for (int i=0; i<saveF.size(); ++i)
      {
        saveF[i]      = Fvec[FindicesVec[i]];
        saveQ[i]      = Qvec[QindicesVec[i]];
        saveB[i]      = Bvec[BindicesVec[i]];
      }
      const std::vector<int> & devStateLIDs = getStaLIDVec();
      for (int i=0; i<saveLastState.size(); ++i)
      {
        saveLastState[i] = lastSta[stateIndicesVec[i]];
        saveCurrState[i] = currSta[stateIndicesVec[i]];
        saveNextState[i] = nextSta[stateIndicesVec[i]];
        saveStateDerivs[i] = nextStaDeriv[stateIndicesVec[i]];
      }

      // zero out the RHS, and re-load, so that we have only the
      // elements from one device present.
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        Fvec[FindicesVec[i]] = 0.0;
        Qvec[QindicesVec[i]] = 0.0;
        Bvec[BindicesVec[i]] = 0.0;
      }
      // re-load for just this instance:
      numJacPtr->loadLocalDAEVectorsIncludingB(*this);

      // Save RHS for just this instance:
      for (int i=0 ; i<origF.size(); ++i)
      {
        origF[i] = Fvec[FindicesVec[i]];
        origQ[i] = Qvec[QindicesVec[i]];
        origB[i] = Bvec[BindicesVec[i]];
      }

      // perturb the parameter
      double epsilon = fabs(Util::MachineDependentParams::MachineEpsilon());
      double sqrtEta= std::sqrt(epsilon);
      double dP = sqrtEta * fabs( origParamValue );
      double minDouble = Util::MachineDependentParams::DoubleMin();
      if (dP < minDouble)
      {
        dP = sqrtEta;
      }

      double newParamValue = origParamValue + dP;

      if (paramName == "")
      {
        setDefaultParam(newParamValue, false);
      }
      else
      {
        setParam(paramName, newParamValue, false);
      }

      processParams (); // if this "entity" is a model, then need to
      // also do a  "processParams" on the related instances.
      processInstanceParams();

      // zero out F,Q,B
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        Fvec[FindicesVec[i]] = 0.0;
        Qvec[QindicesVec[i]] = 0.0;
        Bvec[BindicesVec[i]] = 0.0;
      }

      // load new values for F,Q,B
      numJacPtr->loadLocalDAEVectorsIncludingB(*this);

      // save the perturbed DAE vectors
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        pertF[i] = Fvec[FindicesVec[i]];
        pertQ[i] = Qvec[QindicesVec[i]];
        pertB[i] = Bvec[BindicesVec[i]];
      }

      // compute the derivatives
      for (int i=0; i<pertF.size(); ++i)
      {
        dfdpVec[i] = (pertF[i] - origF[i])/dP;
        dqdpVec[i] = (pertQ[i] - origQ[i])/dP;
        dbdpVec[i] = (pertB[i] - origB[i])/dP;
      }

      // restore everything:
      if (paramName == "")
      {
        setDefaultParam(origParamValue, false);
      }
      else
      {
        setParam(paramName, origParamValue, false);
      }

      processParams (); // if this "entity" is a model, then need to
      // also do a  "processParams" on the related instances.
      processInstanceParams();

      setOrigFlag(origFlag);
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        Fvec[FindicesVec[i]] = saveF[i];
        Qvec[QindicesVec[i]] = saveQ[i];
        Bvec[BindicesVec[i]] = saveB[i];
      }

      for (int i=0; i<saveLastState.size(); ++i)
      {
        lastSta[stateIndicesVec[i]] = saveLastState[i];
        currSta[stateIndicesVec[i]] = saveCurrState[i];
        nextSta[stateIndicesVec[i]] = saveNextState[i];
        nextStaDeriv[stateIndicesVec[i]] = saveStateDerivs[i];
      }

//
      if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
      {
        Xyce::dout() << paramName << ": ";
        Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
        Xyce::dout() << "deviceSens_dp = " << dP << std::endl;


        for (int k1=0; k1<pertF.size(); ++k1)
        {
          int ii=FindicesVec[k1];
          Xyce::dout() 
            <<"pertF["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<pertF[k1]
            <<" origF["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<origF[k1]
            <<" dfdp ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<dfdpVec[k1]
            <<std::endl;
        }

        Xyce::dout() << std::endl;
        for (int k1=0; k1<pertQ.size(); ++k1)
        {
          int ii=QindicesVec[k1];
          Xyce::dout() 
            <<"pertQ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<pertQ[k1]
            <<" origQ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<origQ[k1]
            <<" dqdp ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<dqdpVec[k1]
            <<std::endl;
        }

        Xyce::dout() << std::endl;
        for (int k1=0; k1<pertB.size(); ++k1)
        {
          int ii=BindicesVec[k1];
          Xyce::dout() 
            <<"pertB["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<pertB[k1]
            <<" origB["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<origB[k1]
            <<" dbdp ["<<std::setw(3)<<ii<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<dbdpVec[k1]
            <<std::endl;
        }
      }
//
    }
  }

  return found;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getNumericalMatrixSensitivity
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/21/2019
//-----------------------------------------------------------------------------
bool DeviceInstance::getNumericalMatrixSensitivity ( 
    const std::string & paramName,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & FindicesVec,
    std::vector<int> & QindicesVec,
    std::vector< std::vector<int> > & FjacLIDs,
    std::vector< std::vector<int> > & QjacLIDs )
{
  const std::vector< std::vector<int> > & jacStamp   = jacobianStamp();

  if(devJacLIDs.empty())
  {
    if (DEBUG_DEVICE)
    {
      Report::UserWarning() << getName() << " does not have jacLIDs available";
    }
    return true;
  }

  int numRows, numCols;
  numCols = devLIDs.size();
  numRows = jacStamp.size();

  int testSize= (numRows>numCols)?numRows:numCols;
  if (testSize > mlData.saveF.size())
  {
    mlData.resizeTestJacSolData(testSize);
    mlData.resizeTestJacQData(testSize);
  }

  int numState = getStaLIDVec().size();
  if (numState > mlData.saveCurrState.size())
  {
    mlData.resizeTestJacStateData(numState);
  }

  double origParamValue = 0.0;
  bool found = getParam(paramName, origParamValue); // const?

  if (found)
  {
    // get the vector and matrix indices -------
    // if necessary, consolidate the LIDs vector.
    consolidateDevLIDs();

    if ( !(devLIDs.empty()))
    {
      if (!(FindicesVec.empty())) {  FindicesVec.clear(); }
      if (!(QindicesVec.empty())) {  QindicesVec.clear(); }

      if (!(FjacLIDs.empty())) {  FjacLIDs.clear(); }
      if (!(QjacLIDs.empty())) {  QjacLIDs.clear(); }

      unordered_map<int,int> lidCheckMap;
      unordered_map<int,int> stateLidCheckMap;

      for (int i=0;i<devLIDs.size();++i)
      {
        int newID = devLIDs[i];

        std::vector<int> & jacRow = devJacLIDs[i]; // experiment

        if (lidCheckMap.find(newID) == lidCheckMap.end())
        {
          FindicesVec.push_back(newID);
          lidCheckMap[newID] = i;

          FjacLIDs.push_back(jacRow); // experiment
        }
      }

      std::vector<int> stateIndicesVec;
      const IdVector & stateDevLIDs = getStaLIDVec();
      for (int i=0;i<stateDevLIDs.size();++i)
      {
        int newID = stateDevLIDs[i];
        if (stateLidCheckMap.find(newID) == stateLidCheckMap.end())
        {
          stateIndicesVec.push_back(newID);
          stateLidCheckMap[newID] = i;
        }
      } 

      QindicesVec.insert(QindicesVec.end(), FindicesVec.begin(), FindicesVec.end()); 

      if (devJacLIDs.size() != devLIDs.size() ) // Houston we have a problem
      {
        Report::DevelFatal().in("DeviceInstance::getNumericalMatrixSensitivity") << "Internal Error 1";
      }

      QjacLIDs.resize(FjacLIDs.size());
      for (int i=0 ; i<numRows ; ++i)
      {
        int jCol = FjacLIDs[i].size(); QjacLIDs[i].resize(jCol);
        for (int j=0 ; j<jCol ; ++j) { QjacLIDs[i][j] = FjacLIDs[i][j]; }
      }

      // grab some references
      Linear::Vector & Fvec          = (*extData.daeFVectorPtr);
      Linear::Vector & Qvec          = (*extData.daeQVectorPtr);
      Linear::Vector & Bvec          = (*extData.daeBVectorPtr);

      Linear::Vector & currSol      = (*extData.currSolVectorPtr);
      Linear::Vector & nextSol      = (*extData.nextSolVectorPtr);

      Linear::Vector & lastSta      = (*extData.lastStaVectorPtr);
      Linear::Vector & currSta      = (*extData.currStaVectorPtr);
      Linear::Vector & nextSta      = (*extData.nextStaVectorPtr);
      Linear::Vector & nextStaDeriv = (*extData.nextStaDerivVectorPtr);

      Linear::Matrix & dQdxMat      = (*extData.dQdxMatrixPtr);
      Linear::Matrix & dFdxMat      = (*extData.dFdxMatrixPtr);

      int solSize=FindicesVec.size();
      int stateSize=stateIndicesVec.size();

      mlData.resizeSolnSizedVectors (solSize);
      mlData.resizeStateSizedVectors (stateSize);

      std::vector< std::vector<double> > & numJacF = mlData.numJac;
      std::vector< std::vector<double> > & saveJacF = mlData.saveJac;
      std::vector< std::vector<double> > & devJacF = mlData.devJac;

      std::vector< std::vector<double> > & numJacQ = mlData.numJacQ;
      std::vector< std::vector<double> > & saveJacQ = mlData.saveJacQ;
      std::vector< std::vector<double> > & devJacQ = mlData.devJacQ;

      std::vector< std::vector<int> > & stencil = mlData.stencil;

      std::vector<double> & saveF = mlData.saveF;
      std::vector<double> & saveQ = mlData.saveQ;
      std::vector<double> & saveB = mlData.saveB;

      std::vector<double> & saveLastState = mlData.saveLastState;
      std::vector<double> & saveCurrState = mlData.saveCurrState;
      std::vector<double> & saveNextState = mlData.saveNextState;
      std::vector<double> & saveStateDerivs = mlData.saveStateDerivs;

      saveF.assign(saveF.size(),0.0);
      saveQ.assign(saveQ.size(),0.0);
      saveB.assign(saveB.size(),0.0);

      // zero out the matrix stuff
      for (int i=0;i<numJacF.size();++i)
      {
        numJacF[i].assign(numJacF[i].size(),0.0);
        saveJacF[i].assign(saveJacF[i].size(),0.0);
        devJacF[i].assign(devJacF[i].size(),0.0);
        stencil[i].assign(stencil[i].size(),0);

        numJacQ[i].assign(numJacQ[i].size(),0.0);
        saveJacQ[i].assign(saveJacQ[i].size(),0.0);
        devJacQ[i].assign(devJacQ[i].size(),0.0);
      }

      d_dfdx_dp.clear();
      d_dqdx_dp.clear();
      d_dfdx_dp.resize(numJacF.size());
      d_dqdx_dp.resize(numJacF.size());
      for (int i=0;i<numJacF.size();++i)
      {
        d_dfdx_dp[i].resize(numJacF.size(),0.0);
        d_dqdx_dp[i].resize(numJacF.size(),0.0);
      }

      // save RHS information
      bool origFlag = getOrigFlag();

      for (int i=0; i<saveF.size(); ++i)
      {
        saveF[i]      = Fvec[FindicesVec[i]];
        saveQ[i]      = Qvec[FindicesVec[i]];
        saveB[i]      = Bvec[FindicesVec[i]];
      }
      const std::vector<int> & devStateLIDs = getStaLIDVec();
      for (int i=0; i<saveLastState.size(); ++i)
      {
        saveLastState[i] = lastSta[stateIndicesVec[i]];
        saveCurrState[i] = currSta[stateIndicesVec[i]];
        saveNextState[i] = nextSta[stateIndicesVec[i]];
        saveStateDerivs[i] = nextStaDeriv[stateIndicesVec[i]];
      }

      // re-load for just this instance:
      numJacPtr->loadLocalDAEVectorsIncludingB(*this);

      // save the Jacobian matrix information
      for (int i=0 ; i<numRows ; ++i)
      {
        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
          saveJacF[i][j] = valF;
          double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
          saveJacQ[i][j] = valQ;
        }
      }

      // Zeroing needs to be done after all saved values are
      // recorded because there can be multiple references
      // to the same element
      for (int i=0 ; i<numRows ; ++i)
      {
        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
          dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
        }
      }

      // Now that the original load has been zeroed out, re-load the
      // analytic contributions, to get the contributions from *just* this
      // device.
      loadDAEdQdx ();
      loadDAEdFdx ();

      for (int i=0 ; i<numRows ; ++i)
      {
        devJacF[i].assign(devJacF[i].size(),0.0);
        devJacQ[i].assign(devJacQ[i].size(),0.0);
        stencil[i].assign(stencil[i].size(),0);

        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          double valF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
          devJacF[i][jacStamp[i][j]] = valF;
          double valQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
          devJacQ[i][jacStamp[i][j]] = valQ;
          stencil[i][jacStamp[i][j]] = 1;
        }
      }

      // zero again
      for (int i=0 ; i<numRows ; ++i)
      {
        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          dFdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
          dQdxMat[devLIDs[i]][devJacLIDs[i][j]] = 0;
        }
      }

      // perturb the parameter
      double epsilon = fabs(Util::MachineDependentParams::MachineEpsilon());
      double sqrtEta= std::sqrt(epsilon);
      double dP = sqrtEta * fabs( origParamValue );
      double minDouble = Util::MachineDependentParams::DoubleMin();
      if (dP < minDouble)
      {
        dP = sqrtEta;
      }

      double newParamValue = origParamValue + dP;

      if (paramName == "")
      {
        setDefaultParam(newParamValue, false);
      }
      else
      {
        setParam(paramName, newParamValue, false);
      }

      processParams (); // if this "entity" is a model, then need to
      // also do a  "processParams" on the related instances.
      processInstanceParams();

      // load new values for F,Q,B  (this is done here b/c Jacobian elements computed here)
      numJacPtr->loadLocalDAEVectorsIncludingB(*this);

      loadDAEdQdx ();
      loadDAEdFdx ();

      // Save the new matrix, compute the derivatives
      for (int i=0 ; i<numRows ; ++i)
      {
        int jCol = devJacLIDs[i].size();
        for (int j=0 ; j<jCol ; ++j)
        {
          double pertValF = dFdxMat[devLIDs[i]][devJacLIDs[i][j]];
          numJacF[i][j] = pertValF;
          double pertValQ = dQdxMat[devLIDs[i]][devJacLIDs[i][j]];
          numJacQ[i][j] = pertValQ;

          d_dfdx_dp[i][j] = (pertValF - devJacF[i][j])/dP;
          d_dqdx_dp[i][j] = (pertValQ - devJacQ[i][j])/dP;
        }
      }

      // restore everything:
      if (paramName == "")
      {
        setDefaultParam(origParamValue, false);
      }
      else
      {
        setParam(paramName, origParamValue, false);
      }

      processParams (); // if this "entity" is a model, then need to
      // also do a  "processParams" on the related instances.
      processInstanceParams();

      setOrigFlag(origFlag);
      for (int i=0; i<FindicesVec.size(); ++i)
      {
        Fvec[FindicesVec[i]] = saveF[i];
        Qvec[FindicesVec[i]] = saveQ[i];
        Bvec[FindicesVec[i]] = saveB[i];
      }

      for (int i=0; i<saveLastState.size(); ++i)
      {
        lastSta[stateIndicesVec[i]] = saveLastState[i];
        currSta[stateIndicesVec[i]] = saveCurrState[i];
        nextSta[stateIndicesVec[i]] = saveNextState[i];
        nextStaDeriv[stateIndicesVec[i]] = saveStateDerivs[i];
      }
    } // devLIDs empty
  } // found

  return found;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getNumericalBSensVectorsforAC
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/1/2019
//-----------------------------------------------------------------------------
bool DeviceInstance::getNumericalBSensVectorsforAC ( const std::string & paramName,
          std::vector< std::complex<double> > & dbdp,
          std::vector<int> &        BindicesVec)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::jacStampMap
// Purpose       : Compute Jacobian Stamp and Map for devices that can have merged nodes
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 2/13/04
//-----------------------------------------------------------------------------
void DeviceInstance::jacStampMap(
  const JacobianStamp & stamp_parent,
  std::vector<int> &    map_parent,
  JacobianStamp &       map2_parent,
  JacobianStamp &       stamp,
  std::vector<int> &    map,
  JacobianStamp &       map2,
  int                   from,
  int                   to,
  int                   original_size)
{
  if (from <= to)
  {
    Report::DevelFatal().in("DeviceInstance::jacStampMap")
      << "From index " << from << " <= " << " to index " << to;
  }

  if (map_parent.size() == 0)
  {
    map_parent.resize(original_size);
    map2_parent.resize(original_size);
    for (int i = 0; i < original_size; ++i)
    {
      map_parent[i] = i;
      map2_parent[i].resize(stamp_parent[i].size());
      for (int j = 0; j < stamp_parent[i].size(); ++j)
      {
        map2_parent[i][j] = j;
      }
    }
  }
  map2.resize(original_size);

  // This is to merge the column that is being eliminated into the column it is merged.
  // If extra elements are present then we must increment the map2 value for later
  // elements.  There are multiple cases depending what is populated.
  for (int i = 0; i < original_size; ++i)
  {
    int f_index = -1;
    int t_index = -1;
    int p_row = map_parent[i];
    for (int j = 0; j < stamp_parent[p_row].size(); ++j)
    {
      if (stamp_parent[p_row][j] == from)
        f_index = j;
      if (stamp_parent[p_row][j] == to)
        t_index = j;
    }
    map2[i].resize(map2_parent[i].size());
    int f_index2 = -1;
    int t_index2 = -1;
    for (int j = 0; j < map2_parent[i].size(); ++j)
    {
      map2[i][j] = map2_parent[i][j];
      if (stamp_parent[p_row][map2[i][j]] == from)
        f_index2 = j;
      if (stamp_parent[p_row][map2[i][j]] == to)
        t_index2 = j;
    }
    if (f_index >= 0)
    {
      if (t_index >= 0)
      {
        for (int j = 0; j < map2[i].size(); ++j)
        {
          if (map2[i][j] > f_index)
            map2[i][j]--;
        }
        if (f_index2 >= 0 && t_index2 >= 0)
          map2[i][f_index2] = map2[i][t_index2];
      }
      else
      {
        for (int j = 0; j < map2[i].size(); ++j)
        {
          if (stamp_parent[p_row][map2[i][j]] > to && stamp_parent[p_row][map2[i][j]] < from)
            ++map2[i][j];
        }
        t_index = 0;
        for (int j = 0; j < stamp_parent[p_row].size(); ++j)
        {
          if (to > stamp_parent[p_row][j])
            ++t_index;
          else
            break;
        }
        if (f_index2 >= 0)
          map2[i][f_index2] = t_index;
      }
    }
  }
  map.resize(original_size);
  int p_size = stamp_parent.size();

  // This is to merge the row that is being eliminated into the row it is merged into
  // if extra elements are present then we must increment the map2 value for later
  // elements in the row
  JacobianStamp map2_tmp = map2;
  for (int i = 0; i < stamp_parent[from].size() && stamp_parent[from][i] < p_size-1; ++i)
  {
    bool new_col = false;
    for (int j = 0; j<stamp_parent[to].size(); ++j)
    {
      if (j == 0)
        new_col = true;
      if (stamp_parent[from][i] == stamp_parent[to][j])
      {
        new_col = false;
        break;
      }
    }
    if (new_col)
    {
      for (int j = 0; j < map2[to].size(); ++j)
      {
        if (stamp_parent[to][map2_tmp[to][j]] > stamp_parent[from][i])
        {
          ++map2[to][j];
        }
      }
    }
  }

  stamp.resize(p_size-1);

  std::vector<int> dup(original_size);

  int f_mod = from;
  for (int i = 1; i < original_size; ++i)
  {
    dup[i] = -1;
    for (int j = 0; j < i; ++j)
    {
      if (map_parent[i] == map_parent[j])
      {
        dup[i] = j;
        if (i <= f_mod)
          ++f_mod;
        break;
      }
    }
  }

  for (int i = 0; i < f_mod; ++i)
    map[i] = map_parent[i];
  map[f_mod] = map[to];
  
  for (int i = f_mod + 1; i < original_size; ++i)
    map[i] = map_parent[i]-1;

  for (int i = 1; i < original_size; ++i)
  {
    if (dup[i] >= 0)
      map[i] = map[dup[i]];
  }


  // Now that we know where the row originally came from, we can do any renumbering
  // needed in the map2 source row
  map2_tmp = map2;
  int map_2_from = from;
  if (map[from] != map[to]) {
    for (int i = from + 1; i < original_size; ++i) {
      if (map[i] == map[to]) {
        map_2_from = i;
        break;
      }
    }
    if (map_2_from == from) {
      Report::DevelFatal().in("DeviceInstance::jacStampMap") << "Internal Error 2";
    }
  }

  for (int i = 0; i < stamp_parent[to].size() && stamp_parent[to][i] < p_size-1; ++i)
  {
    bool new_col = false;
    for (int j = 0; j < stamp_parent[from].size(); ++j)
    {
      if (j == 0)
        new_col = true;
      if (stamp_parent[to][i] == stamp_parent[from][j])
      {
        new_col = false;
        break;
      }
    }
    if (new_col)
    {
      for (int j = 0; j < map2[map_2_from].size(); ++j)
      {
        if (stamp_parent[from][map2_tmp[map_2_from][j]] > stamp_parent[to][i])
        {
          ++map2[map_2_from][j];
        }
      }
    }
  }

  JacobianStamp fill(p_size);
  for (int i = 0; i < p_size; ++i)
  {
    fill[i].resize(p_size);
    for (int j = 0 ; j < p_size; ++j)
      fill[i][j] = 0;
    for (int j = 0; j < stamp_parent[i].size(); ++j)
      fill[i][stamp_parent[i][j]] = 1;
  }
  for (int i = 0; i < p_size; ++i)
  {
    fill[to][i] += fill[from][i];
  }
  for (int i = 0; i < p_size; ++i)
  {
    fill[i][to] += fill[i][from];
  }
  for (int i = from; i < p_size - 1; ++i)
  {
    for (int j = 0; j < p_size; ++j)
      fill[i][j] = fill[i+1][j];
    for (int j = 0; j < p_size; ++j)
      fill[j][i] = fill[j][i+1];
  }
  for (int i = 0; i < p_size - 1 ; ++i)
  {
    stamp[i].clear();
    for (int j = 0; j < p_size - 1; ++j)
    {
      if (fill[i][j] > 0)
        stamp[i].push_back(j);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::jacStampMap_fixOrder
//
// Purpose       : This function corrects the compressed row column array so
//                 that the column indices are in ascending order.
//
// Special Notes : The reason for this function is that the
//                 DeviceInstance::jacStampMap function
//                 implicitly requires an ordered jacStamp to work correctly.
//                 Some devices (particularly devices that have meshes),
//                 will start out with an non-ordered stamp, at least in the
//                 column entries.
//
//                 Note, this does require the "map" argument, because,
//                 unlike the function DeviceInstance::jacStampMap,
//                 because this function doesn't change the row ordering,
//                 or remove or merge any rows.
//
//                 This function only changes the column ordering in the
//                 compressed row form of the jacStamp.  It thus requires
//                 modifications to map2, which is essentially a column map.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/01/08
//-----------------------------------------------------------------------------
void
DeviceInstance::jacStampMap_fixOrder(
  const JacobianStamp & stamp_parent,
  JacobianStamp &       map2_parent,
  JacobianStamp &       stamp,
  JacobianStamp &       map2)
{
  int current_size = stamp_parent.size();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_JACSTAMP) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << Xyce::section_divider << std::endl
                 << "Begin DeviceInstance::jacStampMap_fixOrder." << std::endl
                 << Xyce::section_divider << std::endl;
  }

  // if this is the first time this function is called (for a particular device), then
  // allocate the map and set up their trivial contents.
  if (map2_parent.size() == 0)
  {
    map2_parent.resize(current_size);
    for (int i = 0; i < current_size; ++i)
    {
      map2_parent[i].resize(stamp_parent[i].size());
      for (int j = 0; j < stamp_parent[i].size(); ++j)
      {
        map2_parent[i][j] = j;
      }
    }
  }

  stamp.clear();
  map2.clear();

  // To make this simple, start out with a full, dense stamp.
  JacobianStamp denseStamp(current_size);
  for (int i = 0; i < current_size; ++i)
  {
    denseStamp[i].resize(current_size,-1);

    for (int j = 0; j < stamp_parent[i].size(); ++j)
    {
      int denseCol = stamp_parent[i][j];
      denseStamp[i][denseCol] = j;
    }
  }

  // At this point, the denseStamp has been set up.  Now use it to re-create the
  // compressed-row stamp.  By simply looping over the dense stamp, the column order
  // in the compressed row stamp will automatically be ascending.
  // map2 is set up here as well, by pulling the values out that we previously put into
  // dense stamp.
  stamp.resize(current_size);
  map2.resize(current_size);
  for (int i = 0; i < current_size; ++i)
  {
    for (int j = 0; j < current_size; ++j)
    {
      int colMapIndex = denseStamp[i][j];
      if (colMapIndex!=-1)
      {
        stamp[i].push_back(j);
      }
    }

    int stampRowLength=stamp[i].size();
    map2[i].resize(stampRowLength, 0);

    for (int j = 0, k = 0; j < current_size; ++j)
    {
      int colMapIndex = denseStamp[i][j];
      if (colMapIndex!=-1 && colMapIndex < stampRowLength)
      {
        map2[i][colMapIndex] = k;
        ++k;
      }
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_JACSTAMP) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "From inside of DeviceInstance::jacStampMap_fixOrder:" << std::endl
                 << "The original parent stamp is:" << std::endl;
    outputJacStamp(stamp_parent);
    Xyce::dout() << "The new reduced stamp is:" << std::endl;
    outputJacStamp(stamp);
    Xyce::dout() << "The dense stamp is:" << std::endl;
    outputJacStamp(denseStamp);

    Xyce::dout() << "The new map is:" << std::endl;
    outputJacStamp(map2);

    Xyce::dout() << Xyce::section_divider << std::endl
                 << "End DeviceInstance::jacStampMap_fixOrder."<<std::endl
                 << Xyce::section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::computeJacStampAndMaps
// Purpose       : Create a jacStamp and some maps for using it from a list
//                 of row/column pairs representing nonzero elements, and
//                 a list of nodes to collapse (may be zero length)
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 02/10/2016
//-----------------------------------------------------------------------------
///
/// @brief  Jacobian stamp creation function
///
/// The Jacobian is a sparse-matrix representation of the pattern of nonzero
/// elements that a device will load into the Jacobian matrix.
///
/// Given a vector of IntPairs containing the row/column of nonzero elements,
/// and a second vector of IntPairs representing what nodes are to be collapsed
/// onto other nodes, this function constructs the jacStamp in the format that
/// must be passed to topology.
///
/// The function also produces a vector (nodeMap) that tells what node a given
/// original node is mapped onto.
///
/// Also produced is a map that tells what element of the jacStamp corresponds
/// to each of the original nonzero pairs.
///
/// @param[in] jacobianElements Vector of row/column pairs of nonzeros
/// @param[in] collapsedNodes Vector of node collapse pairs (may be empty)
/// @param[out] jacStamp  The jacobian stamp after all collapses are performed
/// @param[out] nodeMap  Vector showing map of original node numbers onto final ordering
/// @param[out] pairToJacStampMap Map of original row/column pairs onto row/column of jacstamp
/// @param[in] numberOfNodes Number of nodes in original system
void DeviceInstance:: computeJacStampAndMaps(const PairVector & jacobianElements,
                                             const PairVector &collapsedNodes,
                                             JacobianStamp & jacStamp,
                                             IdVector & nodeMap,
                                             PairMap &pairToJacStampMap,
                                             const int numberOfNodes)
{
  // Set up trivial node map  (no nodes mapped onto others)
  nodeMap.resize(numberOfNodes);
  for (int i=0;i<numberOfNodes;i++)
    nodeMap[i]=i;

  int reducedNodeCount=nodeMap.size();

  for (PairVector::const_iterator iter=collapsedNodes.begin();
       iter!=collapsedNodes.end(); ++iter)
  {
    int from=(*iter).first;
    int to=(*iter).second;

    // About the only error condition we can encounter here:
    if (from < 0 || from >= numberOfNodes || to < -1 || to >= numberOfNodes)
    {
      Report::DevelFatal().in("DeviceInstance::computeJacStampAndMaps")
        << "In device " << getName() << ": 'from' or 'to' value out of bounds:"
        << "From =" << from << " to = " << to
        << ".  Range of valid values for 'from' is [0,"<<numberOfNodes-1
        << "].  Range of valid values for 'to' is [-1,"<<numberOfNodes-1
        << "]. ";
    }

    // We can only map the way we want if "from" is higher than "to."  If
    // it's not, just make "to" collapse to "from" instead.
    // Use the map for this comparison, because it is possible we have
    // already been mapped.
    if (to != -1 && nodeMap[from] < nodeMap[to])
    {
      int temp=from;
      from=to;
      to=temp;
    }

    // The from node might already have beem mapped or adjusted due to
    // prior collapsing, so get the *real* node that corresponds to the
    // from node at the current step.
    int preMapped = nodeMap[from];

    // to = '-1' means "map to ground", e.g. "remove"
    // "to" might also have been mapped to something other than its
    // original node ID thanks to prior collapsing, so make sure to
    // use its NEW node ID.
    int mappedTo=-1;
    if (to != -1)
      mappedTo= nodeMap[to];

    // Now map the node
    // If other nodes have been mapped to us, we must map all of
    // those onto mappedTo as well!
    for (int i=0;i<nodeMap.size();i++)
    {
      if (i != from && nodeMap[i]==preMapped)
        nodeMap[i]=mappedTo;
    }
    nodeMap[from] = mappedTo;
    reducedNodeCount--;

    // Now decrement all node numbers that are higher than what this node
    // used to be
    for (int i=0;i<nodeMap.size();i++)
    {
      if (nodeMap[i] > preMapped)
        --nodeMap[i];
    }
  }

  // We now need to compute the *UNIQUE* mapped pairs
  // Do this by creating a set into which we dump them.
  PairSet theMappedPairSet;

  for (PairVector::const_iterator iter=jacobianElements.begin();
       iter!=jacobianElements.end();
       ++iter)
  {
    // Sanity check the input elements:  they must always be between
    // 0 and numberOfNodes-1
    if ((*iter).first < 0  || (*iter).first >= numberOfNodes ||
        (*iter).second < 0  || (*iter).second >= numberOfNodes)
    {
      Report::DevelFatal().in("DeviceInstance::computeJacStampAndMaps")
        << "In device " << getName() << ":  Jacobian nonzero-element "
        << " (" << (*iter).first << " , " << (*iter).second << ")"
        << " is in row or column out of bounds! "
        << "Row and column must be in range from 0 to "
        << numberOfNodes-1
        << ".";
      continue;
    }

    if (nodeMap[(*iter).first] != -1 && nodeMap[(*iter).second] != -1)
      theMappedPairSet.insert( IntPair(nodeMap[(*iter).first],nodeMap[(*iter).second]) );
  }

  // Create a mapped jacstamp in the jacStamp array passed in.  We pull out
  // the members of the set (which are guaranteed unique, and will be
  // sorted thanks to our use of set instead of unordered_set.

  // First clear out any garbage that might have been there, then set size
  jacStamp.clear();
  jacStamp.resize(reducedNodeCount);
  PairMap pairToMappedJacMap;

  for (PairSet::iterator iter=theMappedPairSet.begin();
       iter!=theMappedPairSet.end();
       ++iter)
  {
    int jsRow=(*iter).first;
    jacStamp[jsRow].push_back((*iter).second);
    // Now construct a temporary map of mapped pairs to jacstamp entries
    // This is not the one we need to return, as it just tells us
    // where we are putting our unique values, not where the original
    // pairs map.
    pairToMappedJacMap[IntPair(jsRow,(*iter).second)] = IntPair(jsRow,jacStamp[jsRow].size()-1);
  }

  // Finally, create the map from original pairs to mapped jacstamp
  // this map tells us where to find the (originalNode1,originalNode2)
  // nonzero jacobian entry in the final, reduced sparse jacobian stamp.
  // Any node that has been mapped to ground will have had all
  // matrix elements in that row or column removed, and will not exist
  // in the set (and therefore not in the parToMappedJacMap map).  We
  // signal to the caller that the matrix element no longer exists by
  // setting the pair to (-1,-1).
  // It is the caller's responsibility never to attempt to access the
  // jacstamp (or the underlying data that it leads to in the various
  // matrices) when we have said that the jacstamp element is (-1,-1)
  PairMap unmappedPairToMappedJacMap;
  for (PairVector::const_iterator iter=jacobianElements.begin();
       iter!=jacobianElements.end();
       ++iter)
  {
    if ((*iter).first < 0  || (*iter).first >= numberOfNodes ||
        (*iter).second < 0  || (*iter).second >= numberOfNodes)
    {
      // We've already flagged this developer error, press on and ignore
      continue;
    }
    if (nodeMap[(*iter).first] != -1 && nodeMap[(*iter).second] != -1)
    {
      pairToJacStampMap[IntPair((*iter).first,(*iter).second)] = pairToMappedJacMap[IntPair(nodeMap[(*iter).first],nodeMap[(*iter).second])];
    }
    else
    {
      pairToJacStampMap[IntPair((*iter).first,(*iter).second)] = IntPair(-1,-1);
    }
  }

}


//-----------------------------------------------------------------------------
// Function      : DeviceInstance::outputJacStamp
// Purpose       : Output jacStamp (for debugging)
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/19/06
//-----------------------------------------------------------------------------
void
DeviceInstance::outputJacStamp(
  const JacobianStamp &         jac)
{
  for (int i = 0 ; i < jac.size(); ++i)
  {
    Xyce::dout() << "Row: " << i << " ::";
    for (int j = 0 ; j < jac[i].size(); ++j)
      Xyce::dout() << "  " << jac[i][j];

    Xyce::dout() << std::endl;
  }
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::outputJacMaps
// Purpose       : Output jacMap and jacMap2 (for debugging)
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli, Electrical & Microsystems Modeling
// Creation Date : 02/20/08
//-----------------------------------------------------------------------------
void
DeviceInstance::outputJacMaps(
  const std::vector<int>  &     jacMap,
  const JacobianStamp &         jacMap2)
{
  for (int i = 0 ; i < jacMap.size(); ++i)
  {
    Xyce::dout() << "Row " << i << ": ";
    for (int j = 0; j < jacMap2[i].size(); j++)
      Xyce::dout() << jacMap[i]<< "," << jacMap2[i][j] << " ";

    Xyce::dout() << std::endl;
  }

  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/00
//-----------------------------------------------------------------------------
bool DeviceInstance::updateTemperature(const double & temp_tmp)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool DeviceInstance::processParams ()
{
  Report::DevelFatal0().in("DeviceInstance::processParams")
    << "DeviceInstance::processParams() must be implemented for device " << getName();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setInternalState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
bool DeviceInstance::setInternalState(
  const DeviceState & state )
{
  Report::DevelFatal().in("DeviceInstance::setInternalState") << "does not exist for this device " << getName();

  return false;
}

std::ostream &
DeviceInstance::printName(std::ostream &os) const
{
  return os << "instance " << name_.getEncodedName();
}

} // namespace Device
} // namespace Xyce

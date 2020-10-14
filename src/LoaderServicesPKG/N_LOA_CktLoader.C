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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_LOA_CktLoader.C,v $
//
// Purpose        : This file contains class definitions for the loader
//                  services package.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/26/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.104 $
//
// Revision Date  : $Date: 2016/03/04 00:34:49 $
//
// Current Owner  : $Author: hkthorn $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>
#include <N_DEV_fwd.h>

#include <N_ANP_SweepParam.h>

#include <N_DEV_Algorithm.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_SolverState.h>
#include <N_LOA_CktLoader.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_FilteredMatrix.h>
#include <N_LAS_Builder.h>
#include <N_NLS_TwoLevelEnum.h>
#include <N_UTL_Stats.h>

namespace Xyce {
namespace Loader {

//-----------------------------------------------------------------------------
// Function      : CktLoader::CktLoader
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
CktLoader::CktLoader(
  Device::DeviceMgr &         device_manager,
  Linear::Builder &           builder)
  : deviceManager_(device_manager),
    builder_(builder),
    dcopState_(false),
    lindQdxMatrixPtr_(0),
    lindFdxMatrixPtr_(0),
    filtered_lindQdxMatrixPtr_(0),
    filtered_lindFdxMatrixPtr_(0)
{}

//-----------------------------------------------------------------------------
// Function      : CktLoader::~CktLoader
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/27/00
//-----------------------------------------------------------------------------
CktLoader::~CktLoader()
{
  delete lindQdxMatrixPtr_;
  delete lindFdxMatrixPtr_;
  delete filtered_lindQdxMatrixPtr_;
  delete filtered_lindFdxMatrixPtr_;
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::setParam
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::setParam(
  std::string &   name,
  double                val,
  bool overrideOriginal) 
{
  // Delete the current linear matrices, just in case the parameter affects
  // any linear devices and their Jacobian entries.
  delete lindQdxMatrixPtr_; lindQdxMatrixPtr_=0;
  delete lindFdxMatrixPtr_; lindFdxMatrixPtr_=0;
  delete filtered_lindQdxMatrixPtr_; filtered_lindQdxMatrixPtr_=0;
  delete filtered_lindFdxMatrixPtr_; filtered_lindFdxMatrixPtr_=0;

  return deviceManager_.setParam(name, val, overrideOriginal);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::setParamRandomExpressionTerms
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::setParamRandomExpressionTerms(
  std::string &   name,
  std::string &   opName,
  int             opIndex,
  //enum Util::astRandTypes astType,
  int astType,
  double                val,
  bool overrideOriginal) 
{
  // Delete the current linear matrices, just in case the parameter affects
  // any linear devices and their Jacobian entries.
  delete lindQdxMatrixPtr_; lindQdxMatrixPtr_=0;
  delete lindFdxMatrixPtr_; lindFdxMatrixPtr_=0;
  delete filtered_lindQdxMatrixPtr_; filtered_lindQdxMatrixPtr_=0;
  delete filtered_lindFdxMatrixPtr_; filtered_lindFdxMatrixPtr_=0;

  return deviceManager_.setParamRandomExpressionTerms(name, opName, opIndex, astType, val, overrideOriginal);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::analyticSensitivitiesAvailable
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::analyticSensitivitiesAvailable(
  const std::string &   name)
{
  return deviceManager_.analyticSensitivitiesAvailable(name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::numericalSensitivitiesAvailable
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::numericalSensitivitiesAvailable (const std::string & name)
{
  return deviceManager_.numericalSensitivitiesAvailable (name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getAnalyticSensitivities
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktLoader::getAnalyticSensitivities(
  std::string &   name, 
  std::vector<double> & dfdpVec, 
  std::vector<double> & dqdpVec,
  std::vector<double> & dbdpVec,
  std::vector<int> &    FindicesVec,
  std::vector<int> &    QindicesVec,
  std::vector<int> &    BindicesVec) const
{
  return deviceManager_.getAnalyticSensitivities(name, dfdpVec, dqdpVec, dbdpVec, 
                                                  FindicesVec, QindicesVec, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNumericalSensitivities
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktLoader::getNumericalSensitivities(
  std::string &         name,
  std::vector<double> & dfdpVec,
  std::vector<double> & dqdpVec,
  std::vector<double> & dbdpVec,
  std::vector<int> &    FindicesVec,
  std::vector<int> &    QindicesVec,
  std::vector<int> &    BindicesVec) const
{
  return deviceManager_.getNumericalSensitivities(name, dfdpVec, dqdpVec, dbdpVec, FindicesVec, QindicesVec, BindicesVec);
}

//
//-----------------------------------------------------------------------------
// Function      : CktLoader::analyticBVecSensAvailable
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::analyticBVecSensAvailable(const std::string & name)
{
  return deviceManager_.analyticBVecSensAvailable(name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::numericalBVecSensAvailable
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::numericalBVecSensAvailable(const std::string & name)
{
  return deviceManager_.numericalBVecSensAvailable(name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::analyticMatrixSensitivitiesAvailable
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::analyticMatrixSensitivitiesAvailable(
  const std::string &   name)
{
  return deviceManager_.analyticMatrixSensitivitiesAvailable(name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::numericalMatrixSensitivitiesAvailable
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::numericalMatrixSensitivitiesAvailable (const std::string & name)
{
  return deviceManager_.numericalMatrixSensitivitiesAvailable (name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getAnalyticBSensVectorsforAC
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/16/2019
//-----------------------------------------------------------------------------
void  CktLoader::getAnalyticBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const
{
  return deviceManager_.getAnalyticalBSensVectorsforAC (name, dbdp, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNumericalBSensVectorsforAC
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/16/2019
//-----------------------------------------------------------------------------
void  CktLoader::getNumericalBSensVectorsforAC (const std::string & name,
          std::vector< std::complex<double> > &     dbdp,
          std::vector<int> &        BindicesVec) const
{
  return deviceManager_.getNumericalBSensVectorsforAC (name, dbdp, BindicesVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getAnalyticMatrixSensitivities
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktLoader::getAnalyticMatrixSensitivities(
  const std::string &   name, 
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const
{
  return deviceManager_.getAnalyticMatrixSensitivities(name, 
                  d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNumericalMatrixSensitivities
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
void CktLoader::getNumericalMatrixSensitivities(
  const std::string &         name,
    std::vector <std::vector<double> > & d_dfdx_dp,
    std::vector <std::vector<double> > & d_dqdx_dp,
    std::vector<int> & F_lids,
    std::vector<int> & Q_lids,
    std::vector< std::vector<int> > & F_jacLIDs,
    std::vector< std::vector<int> > & Q_jacLIDs) const
{
  return deviceManager_.getNumericalMatrixSensitivities(name, 
                  d_dfdx_dp, d_dqdx_dp, F_lids, Q_lids, F_jacLIDs, Q_jacLIDs);
}
//

//-----------------------------------------------------------------------------
// Function      : CktLoader::getParamAndReduce
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
double CktLoader::getParamAndReduce(
  Parallel::Machine     comm,
  const std::string &   name) const
{
  return Device::getParamAndReduce(comm, deviceManager_, name);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getRandomParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
void CktLoader::getRandomParams(std::vector<Xyce::Analysis::SweepParam> & SamplingParams, Parallel::Communicator & parallel_comm)
{
  return deviceManager_.getRandomParams(SamplingParams,parallel_comm);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::initializeProblem
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
bool CktLoader::initializeProblem(
  Linear::Vector * nextSolVectorPtr,
  Linear::Vector * currSolVectorPtr,
  Linear::Vector * lastSolVectorPtr,
  Linear::Vector * nextStaVectorPtr,
  Linear::Vector * currStaVectorPtr,
  Linear::Vector * lastStaVectorPtr,
  Linear::Vector * StateDerivVectorPtr,
  Linear::Vector * nextStoVectorPtr,
  Linear::Vector * currStoVectorPtr,
  Linear::Vector * QVectorPtr,
  Linear::Vector * FVectorPtr,
  Linear::Vector * BVectorPtr,
  Linear::Vector * dFdxdVpVectorPtr,
  Linear::Vector * dQdxdVpVectorPtr) const
{
  return deviceManager_.setICs(
    nextSolVectorPtr,
    currSolVectorPtr,
    lastSolVectorPtr,
    nextStaVectorPtr,
    currStaVectorPtr,
    lastStaVectorPtr,
    StateDerivVectorPtr,
    nextStoVectorPtr,
    currStoVectorPtr,
    QVectorPtr,
    FVectorPtr,
    BVectorPtr,
    dFdxdVpVectorPtr,
    dQdxdVpVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadDeviceErrorWeightMask
// Purpose       : Load mask for use in calculating weighted norms.
// Special Notes : Loads the mask into the LAS_System.  Need only be called
//                 once.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 1/19/07
//-----------------------------------------------------------------------------
bool CktLoader::loadDeviceErrorWeightMask(
  Linear::Vector *      deviceMask) const
{
  return deviceManager_.loadErrorWeightMask(deviceMask);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadDAEMatrices
// Purpose       : This function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool CktLoader::loadDAEMatrices(
  Linear::Vector *        tmpSolVectorPtr,
  Linear::Vector *        tmpStateVectorPtr,
  Linear::Vector *        tmpStateDerivVectorPtr,
  Linear::Vector *        tmpStoreVectorPtr,
  Linear::Matrix *        tmpdQdxMatrixPtr,
  Linear::Matrix *        tmpdFdxMatrixPtr,
  int loadType)
{

#ifdef Xyce_PARALLEL_MPI
  tmpSolVectorPtr->importOverlap();
#endif

  bool bsuccess = true;

  if (loadType != Xyce::Device::ALL)
  {
    bsuccess &= deviceManager_.loadDAEMatrices(
      tmpSolVectorPtr,
      tmpStateVectorPtr,
      tmpStateDerivVectorPtr,
      tmpStoreVectorPtr,
      tmpdQdxMatrixPtr,
      tmpdFdxMatrixPtr,
      loadType);
  }
  else
  {
    bool isInnerLoad = deviceManager_.getSolverState().twoLevelNewtonCouplingMode == Nonlinear::INNER_PROBLEM;

    // Handle the two-level matrix load at this level to reduce code copying in
    // the device manager.
    if (isInnerLoad)
    {
      bsuccess &= deviceManager_.loadDAEMatrices(
        tmpSolVectorPtr,
        tmpStateVectorPtr,
        tmpStateDerivVectorPtr,
        tmpStoreVectorPtr,
        tmpdQdxMatrixPtr,
        tmpdFdxMatrixPtr,
        Xyce::Device::PDE);
    }
    else
    {
      if (deviceManager_.getDeviceOptions().separateLoad
          && !deviceManager_.getSolverState().dcopFlag
          && !deviceManager_.getDeviceOptions().testJacobianFlag)
      {
        // Simulator state has changed, reload linear portion of Jacobian matrices.
        // Also reload the Jacobian matrix if this is a DCOP calculation or the
        // analytic Jacobian is being tested.
        bool reloadLinDevs = false;
        if ((!lindQdxMatrixPtr_ && !lindFdxMatrixPtr_)
          || dcopState_ != deviceManager_.getSolverState().dcopFlag)
        {
          dcopState_ = deviceManager_.getSolverState().dcopFlag;
          reloadLinDevs = true;
        }
  
        if ( reloadLinDevs )
        {
          //Stats::StatTop _timeJacobianStat("Load Linear Matrix");
          //Stats::TimeBlock _timeJacobianTimer(_timeJacobianStat);

          delete lindQdxMatrixPtr_; lindQdxMatrixPtr_=0;
          delete lindFdxMatrixPtr_; lindFdxMatrixPtr_=0;
          delete filtered_lindQdxMatrixPtr_; filtered_lindQdxMatrixPtr_=0;
          delete filtered_lindFdxMatrixPtr_; filtered_lindFdxMatrixPtr_=0;

          // Load linear devices into lindQdxMatrixPtr_ and lindFdxMatrixPtr_.
          lindQdxMatrixPtr_ = builder_.createMatrix();
          lindFdxMatrixPtr_ = builder_.createMatrix();

          bsuccess &= deviceManager_.loadDAEMatrices(
            tmpSolVectorPtr,
            tmpStateVectorPtr,
            tmpStateDerivVectorPtr,
            tmpStoreVectorPtr,
            tmpdQdxMatrixPtr,
            tmpdFdxMatrixPtr,
            Xyce::Device::LINEAR);

          lindQdxMatrixPtr_->addOverlap(*tmpdQdxMatrixPtr);
          lindFdxMatrixPtr_->addOverlap(*tmpdFdxMatrixPtr);
        }

        // Load nonlinear devices.
        bsuccess &= deviceManager_.loadDAEMatrices(
          tmpSolVectorPtr,
          tmpStateVectorPtr,
          tmpStateDerivVectorPtr,
          tmpStoreVectorPtr,
          tmpdQdxMatrixPtr,
          tmpdFdxMatrixPtr,
          Xyce::Device::NONLINEAR);
  
        if ( !reloadLinDevs )
        {
          //Stats::StatTop _timeJacobianStat("Copy Linear Matrix");
          //Stats::TimeBlock _timeJacobianTimer(_timeJacobianStat);

          if (!filtered_lindQdxMatrixPtr_ && !filtered_lindFdxMatrixPtr_)
          {
            lindQdxMatrixPtr_->fillComplete();
            lindFdxMatrixPtr_->fillComplete();
            filtered_lindQdxMatrixPtr_ = new Xyce::Linear::FilteredMatrix( lindQdxMatrixPtr_, tmpSolVectorPtr->omap() );
            filtered_lindFdxMatrixPtr_ = new Xyce::Linear::FilteredMatrix( lindFdxMatrixPtr_, tmpSolVectorPtr->omap() );
          }

          // Copy lindFdxMatrixPtr_ and lindQdxMatrixPtr_ into tmp matrices.
          if ( !filtered_lindQdxMatrixPtr_->isEmpty() ) 
          {
            filtered_lindQdxMatrixPtr_->addToMatrix( *tmpdQdxMatrixPtr );
          }
          if ( !filtered_lindFdxMatrixPtr_->isEmpty() )
          { 
            filtered_lindFdxMatrixPtr_->addToMatrix( *tmpdFdxMatrixPtr );
          }
        }
      }
      else
      {
        dcopState_ = deviceManager_.getSolverState().dcopFlag;

        bsuccess &= deviceManager_.loadDAEMatrices(
          tmpSolVectorPtr,
          tmpStateVectorPtr,
          tmpStateDerivVectorPtr,
          tmpStoreVectorPtr,
          tmpdQdxMatrixPtr,
          tmpdFdxMatrixPtr,
          Xyce::Device::ALL);
      }
    }
  }

  // Tell Jacobian, fill is complete allowing accumulation if necessary
  {
    //Stats::StatTop _timeJacobianStat("Fill Complete Matrix");
    //Stats::TimeBlock _timeJacobianTimer(_timeJacobianStat);
  tmpdQdxMatrixPtr->fillComplete();
  tmpdFdxMatrixPtr->fillComplete();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool CktLoader::loadDAEVectors   (
  Linear::Vector * nextSolVectorPtr,
  Linear::Vector * currSolVectorPtr,
  Linear::Vector * lastSolVectorPtr,
  Linear::Vector * nextStaVectorPtr,
  Linear::Vector * currStaVectorPtr,
  Linear::Vector * lastStaVectorPtr,
  Linear::Vector * StateDerivVectorPtr,
  Linear::Vector * nextStoVectorPtr,
  Linear::Vector * currStoVectorPtr,
  Linear::Vector * nextLeadFVectorPtr,
  Linear::Vector * nextLeadQVectorPtr,
  Linear::Vector * nextJunctionVVectorPtr,
  Linear::Vector * QVectorPtr,
  Linear::Vector * FVectorPtr,
  Linear::Vector * BVectorPtr,
  Linear::Vector * dFdxdVpVectorPtr,
  Linear::Vector * dQdxdVpVectorPtr,
  int loadType)
{
  bool bsuccess = true;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  nextSolVectorPtr->importOverlap();
  StateDerivVectorPtr->importOverlap();
#endif

  // Allow level above CktLoader to specify types of devices to load.
  // This is used in HB.
  if (loadType != Xyce::Device::ALL)
  {
    bsuccess &= deviceManager_.loadDAEVectors(
      nextSolVectorPtr,
      currSolVectorPtr,
      lastSolVectorPtr,
      nextStaVectorPtr,
      currStaVectorPtr,
      lastStaVectorPtr,
      StateDerivVectorPtr,
      nextStoVectorPtr,
      currStoVectorPtr,
      nextLeadFVectorPtr,
      nextLeadQVectorPtr,
      nextJunctionVVectorPtr,
      QVectorPtr,
      FVectorPtr,
      BVectorPtr,
      dFdxdVpVectorPtr,
      dQdxdVpVectorPtr,
      loadType);
  }
  else
  {
    bool isInnerLoad = deviceManager_.getSolverState().twoLevelNewtonCouplingMode == Nonlinear::INNER_PROBLEM;

    // Handle the two-level matrix load at this level to reduce code copying in
    // the device manager.
    if (isInnerLoad)
    {
      bsuccess &= deviceManager_.loadDAEVectors(
        nextSolVectorPtr,
        currSolVectorPtr,
        lastSolVectorPtr,
        nextStaVectorPtr,
        currStaVectorPtr,
        lastStaVectorPtr,
        StateDerivVectorPtr,
        nextStoVectorPtr,
        currStoVectorPtr,
        nextLeadFVectorPtr,
        nextLeadQVectorPtr,
        nextJunctionVVectorPtr,
        QVectorPtr,
        FVectorPtr,
        BVectorPtr,
        dFdxdVpVectorPtr,
        dQdxdVpVectorPtr,
        Xyce::Device::PDE);
    }
    else
    {
      if (deviceManager_.getDeviceOptions().separateLoad)
      {
        // Simulator state has changed, reload linear portion of vectors.
        // Also reload the vectors if this system has any PDE devices,
        // only during the DCOP phase because of the double-DCOP calculation.
        // Furthermore, reload linear components if there are no stored matrices.
        bool reloadLinDevs = (filtered_lindQdxMatrixPtr_ && filtered_lindFdxMatrixPtr_)? false : true;
        if ( deviceManager_.getSolverState().dcopFlag || deviceManager_.getDeviceOptions().testJacobianFlag )
        {
          reloadLinDevs = true;
        }
    
        if ( reloadLinDevs )
        {
          //Stats::StatTop _timeRHSStat("Load Linear RHS");
          //Stats::TimeBlock _timeRHSTimer(_timeRHSStat);

          // Load linear devices into QVectorPtr, FVectorPtr, BVectorPtr
          bsuccess &= deviceManager_.loadDAEVectors(
            nextSolVectorPtr,
            currSolVectorPtr,
            lastSolVectorPtr,
            nextStaVectorPtr,
            currStaVectorPtr,
            lastStaVectorPtr,
            StateDerivVectorPtr,
            nextStoVectorPtr,
            currStoVectorPtr,
            nextLeadFVectorPtr,
            nextLeadQVectorPtr,
            nextJunctionVVectorPtr,
            QVectorPtr,
            FVectorPtr,
            BVectorPtr,
            dFdxdVpVectorPtr,
            dQdxdVpVectorPtr,
            Xyce::Device::ALL);
        }
        else
        {
          {
           // Stats::StatTop _timeVectorStat("Nonlinear Vector Load");
           // Stats::TimeBlock _timeVectorTimer(_timeVectorStat);
  
          // Load nonlinear devices.
          bsuccess &= deviceManager_.loadDAEVectors(
            nextSolVectorPtr,
            currSolVectorPtr,
            lastSolVectorPtr,
            nextStaVectorPtr,
            currStaVectorPtr,
            lastStaVectorPtr,
            StateDerivVectorPtr,
            nextStoVectorPtr,
            currStoVectorPtr,
            nextLeadFVectorPtr,
            nextLeadQVectorPtr,
            nextJunctionVVectorPtr,
            QVectorPtr,
            FVectorPtr,
            BVectorPtr,
            dFdxdVpVectorPtr,
            dQdxdVpVectorPtr,
            Xyce::Device::NONLINEAR);
          }

          // Acquire the same effect of the linear device load by using
          // a matrix-vector product with the solution vector and filtered
          // matrices.
          {
            //Stats::StatTop _timeVectorStat("New Linear Vector Load");
            //Stats::TimeBlock _timeVectorTimer(_timeVectorStat);
          if ( !filtered_lindFdxMatrixPtr_->isEmpty() ) 
          {
            filtered_lindFdxMatrixPtr_->axpy( *nextSolVectorPtr, *FVectorPtr );
          }

          if ( !filtered_lindQdxMatrixPtr_->isEmpty() ) 
          {
            filtered_lindQdxMatrixPtr_->axpy( *nextSolVectorPtr, *QVectorPtr );
          }
     
          // Load BVectorPtr here for linear devices.
          deviceManager_.loadBVectorsforSources();
          }
        }
      }
      else
      {
        bsuccess &= deviceManager_.loadDAEVectors(
          nextSolVectorPtr,
          currSolVectorPtr,
          lastSolVectorPtr,
          nextStaVectorPtr,
          currStaVectorPtr,
          lastStaVectorPtr,
          StateDerivVectorPtr,
          nextStoVectorPtr,
          currStoVectorPtr,
          nextLeadFVectorPtr,
          nextLeadQVectorPtr,
          nextJunctionVVectorPtr,
          QVectorPtr,
          FVectorPtr,
          BVectorPtr,
          dFdxdVpVectorPtr,
          dQdxdVpVectorPtr,
          Xyce::Device::ALL);
      } 
    }
  }

  // Update parallel if necessary
  QVectorPtr->fillComplete();
  FVectorPtr->fillComplete();
  BVectorPtr->fillComplete();
/*
  std::cout << "Q vector: " << std::endl;
  QVectorPtr->printPetraObject( std::cout );
  std::cout << "F vector: " << std::endl;
  FVectorPtr->printPetraObject( std::cout );
  std::cout << "B vector: " << std::endl;
  BVectorPtr->printPetraObject( std::cout );
*/

  dFdxdVpVectorPtr->fillComplete();
  dQdxdVpVectorPtr->fillComplete();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::updateFDIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo and Ting Mei
// Creation Date : 19 March 2018
//-----------------------------------------------------------------------------
/// Update device intermediate variables in the frequency domain
///
/// \author Tom Russo, Ting Mei
/// \date 19 March 2018
///
bool CktLoader::updateFDIntermediateVars(
    double frequency,
    std::complex<double>* freqSolVec)
{
  return deviceManager_.updateFDIntermediateVars(frequency, freqSolVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadFreqDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Jason C. Verley
// Creation Date : 07/13/17
//-----------------------------------------------------------------------------
/// Load frequency domain DAE matrices
///
/// Loads the DAE contributions to the Jacobian matrix.  To be used by
/// frequency-based analysis techniques.  Used only for devices that have a
/// frequency response.
///
/// \author Jason C. Verley
/// \date 07/13/17
///
bool CktLoader::loadFreqDAEMatrices(
    double frequency,
    std::complex<double>* freqSolVec,
    std::vector<Util::FreqMatEntry>& dFdxEntries)
{
  return deviceManager_.loadFreqDAEMatrices(frequency, freqSolVec, dFdxEntries);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadFreqDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Jason C. Verley
// Creation Date : 07/13/17
//-----------------------------------------------------------------------------
/// Load frequency domain DAE vectors
///
/// Loads the contributions to the `F` vector.  To be used by frequency-based
/// analyses techniques.  Used only for devices that have a frequency response.
///
/// \author Jason C. Verley
/// \date 07/13/17
///
bool CktLoader::loadFreqDAEVectors(
    double frequency,
    std::complex<double>* freqSolVec,
    std::vector<Util::FreqVecEntry>& FVecEntries,
    std::vector<Util::FreqVecEntry>& BVecEntries)
{
  return deviceManager_.loadFreqDAEVectors(frequency, freqSolVec,
                                           FVecEntries, BVecEntries);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool CktLoader::updateState(
  Linear::Vector * nextSolVectorPtr,
  Linear::Vector * currSolVectorPtr,
  Linear::Vector * lastSolVectorPtr,
  Linear::Vector * nextStaVectorPtr,
  Linear::Vector * currStaVectorPtr,
  Linear::Vector * lastStaVectorPtr,
  Linear::Vector * nextStoVectorPtr,
  Linear::Vector * currStoVectorPtr,
  int loadType)
{
      //Stats::StatTop _timeUpdateStat("Update State");
      //Stats::TimeBlock _timeUpdateTimer(_timeUpdateStat);
 
  if (loadType != Xyce::Device::ALL)
  {
    return deviceManager_.updateState(
      nextSolVectorPtr,
      currSolVectorPtr,
      lastSolVectorPtr,
      nextStaVectorPtr,
      currStaVectorPtr,
      lastStaVectorPtr,
      nextStoVectorPtr,
      currStoVectorPtr,
      loadType);
  }
  else
  { 
    bool isInnerLoad = deviceManager_.getSolverState().twoLevelNewtonCouplingMode == Nonlinear::INNER_PROBLEM;

    // Handle the two-level matrix load at this level to reduce code copying in
    // the device manager.
    if (isInnerLoad)
    {
      return deviceManager_.updateState(
        nextSolVectorPtr,
        currSolVectorPtr,
        lastSolVectorPtr,
        nextStaVectorPtr,
        currStaVectorPtr,
        lastStaVectorPtr,
        nextStoVectorPtr,
        currStoVectorPtr,
        Xyce::Device::PDE);
    }
    else
    {
      // Simulator state has changed, update linear devices as well as nonlinear devices.
      bool updateAllDevs = (filtered_lindQdxMatrixPtr_ && filtered_lindFdxMatrixPtr_)? false : true;
      if ( deviceManager_.getSolverState().dcopFlag || deviceManager_.getDeviceOptions().testJacobianFlag )
      {
        updateAllDevs = true;
      }
  
      if ( updateAllDevs )
      {
        return deviceManager_.updateState(
          nextSolVectorPtr,
          currSolVectorPtr,
          lastSolVectorPtr,
          nextStaVectorPtr,
          currStaVectorPtr,
          lastStaVectorPtr,
          nextStoVectorPtr,
          currStoVectorPtr,
          Xyce::Device::ALL);
      } 
      else
      {
        return deviceManager_.updateState(
          nextSolVectorPtr,
          currSolVectorPtr,
          lastSolVectorPtr,
          nextStaVectorPtr,
          currStaVectorPtr,
          lastStaVectorPtr,
          nextStoVectorPtr,
          currStoVectorPtr,
          Xyce::Device::NONLINEAR);
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadBVectorsforAC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool CktLoader::loadBVectorsforAC(
  Linear::Vector * bVecRealPtr,
  Linear::Vector * bVecImagPtr)
{
  return deviceManager_.loadBVectorsforAC (bVecRealPtr, bVecImagPtr);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadBVectorsforSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool CktLoader::loadBVectorsforSources()
{
  return deviceManager_.loadBVectorsforSources();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNumNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
int CktLoader::getNumNoiseSources ()
{
  return deviceManager_.getNumNoiseSources ();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNumNoiseDevices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
int CktLoader::getNumNoiseDevices ()
{
  return deviceManager_.getNumNoiseDevices ();
}


//-----------------------------------------------------------------------------
// Function      : CktLoader::setupNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
void CktLoader::setupNoiseSources 
  (std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec)
{
  return deviceManager_.setupNoiseSources (noiseDataVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNoiseSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/8/2014
//-----------------------------------------------------------------------------
void CktLoader::getNoiseSources 
  (std::vector<Xyce::Analysis::NoiseData*> & noiseDataVec)
{
  return deviceManager_.getNoiseSources (noiseDataVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getBMatrixStampforMOR
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool CktLoader::getBMatrixEntries(std::vector<int>& bMatEntriesVec,
                                        std::vector<int>& portVec,
                                        std::vector<double> * Z0sVec )
{
  return deviceManager_.getBMatrixEntries( bMatEntriesVec, portVec, Z0sVec );
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getInductorsEntriesforMOR
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
/*bool CktLoader::getInductorsEntriesforMOR(std::vector<int>& inductorEntriesVec)
{
  return deviceManager_.getInductorsEntriesforMOR( inductorEntriesVec );
}
*/
//-----------------------------------------------------------------------------
// Function      : CktLoader::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool CktLoader::setInitialGuess(Linear::Vector * solVectorPtr)
{
  return deviceManager_.setInitialGuess (solVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::updateSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
bool CktLoader::updateSources()
{
  return deviceManager_.updateSources();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool CktLoader::outputPlotFiles() const
{
  return deviceManager_.outputPlotFiles(false);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::finishOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/19/04
//-----------------------------------------------------------------------------
bool CktLoader::finishOutput() const
{
  return deviceManager_.finishOutput();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::isPDESystem
// Purpose       : This is an accessor to allow the time integrator to determine
//                 if the current problem includes a PDE device.  If it does,
//                 then it makes sense to use two-pass DCOP calulation.  Hence
//                 the "DoubleDCOPFlag" in the name.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
bool CktLoader::isPDESystem() const
{
  return deviceManager_.isPDESystem();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getLimiterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/11/04
//-----------------------------------------------------------------------------
bool CktLoader::getLimiterFlag ()
{
  return deviceManager_.getVoltageLimiterFlag ();
}


//-----------------------------------------------------------------------------
// Function      : CktLoader::resetBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/17/2019
//-----------------------------------------------------------------------------
void CktLoader::resetBreakPoints()// needed for 2-level
{
  return deviceManager_.resetBreakPoints();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/11/01
//-----------------------------------------------------------------------------
bool CktLoader::getBreakPoints ( 
    std::vector<Util::BreakPoint> & breakPointTimes,
    std::vector<Util::BreakPoint> & pauseBreakPointTimes
    ) const
{
  return deviceManager_.getBreakPoints(breakPointTimes, pauseBreakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getMaxTimeStepSize ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
double CktLoader::getMaxTimeStepSize ()
{
  return deviceManager_.getMaxTimeStepSize ();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int CktLoader::enablePDEContinuation()
{
  return deviceManager_.enablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool CktLoader::disablePDEContinuation ()
{
  return deviceManager_.disablePDEContinuation();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getNumInterfaceNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
void CktLoader::getNumInterfaceNodes (std::vector<int> & numINodes)
{
  deviceManager_.getNumInterfaceNodes (numINodes);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool CktLoader::loadCouplingRHS(int iSubProblem, int iCouple, Linear::Vector * dfdvPtr)
{
  return deviceManager_.loadCouplingRHS (iSubProblem, iCouple, dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool CktLoader::calcCouplingTerms (int iSubProblem, int iCouple, const Linear::Vector * dxdvPtr)
{
  return deviceManager_.calcCouplingTerms (iSubProblem, iCouple, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getHomotopyBlockSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL, Parallel Computational Sciences
// Creation Date : 01/26/2005
//-----------------------------------------------------------------------------
int CktLoader::getHomotopyBlockSize() const
{
  return deviceManager_.getHomotopyBlockSize();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::allDevsConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
bool CktLoader::allDevicesConverged(
  Parallel::Machine     comm)
{
  return deviceManager_.allDevicesConverged(comm);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/21/06
//-----------------------------------------------------------------------------
bool
CktLoader::innerDevicesConverged(
  Parallel::Machine     comm)
{
  return Device::devicesConverged(comm, deviceManager_.getDevices(Device::ExternDevice::Traits::modelType()));
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void CktLoader::stepSuccess(Xyce::Analysis::TwoLevelMode analysis)
{
  const Device::InstanceVector &extern_devices = deviceManager_.getDevices(Device::ExternDevice::Traits::modelType());

  for (Device::InstanceVector::const_iterator it = extern_devices.begin(); it != extern_devices.end(); ++it)
  {
    Device::ExternDevice::Instance &extern_device = static_cast<Device::ExternDevice::Instance &>(*(*it));

    extern_device.stepSuccess(analysis);
  }
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void CktLoader::stepFailure(Xyce::Analysis::TwoLevelMode analysis)
{
  const Device::InstanceVector &extern_devices = deviceManager_.getDevices(Device::ExternDevice::Traits::modelType());

  for (Device::InstanceVector::const_iterator it = extern_devices.begin(); it != extern_devices.end(); ++it)
  {
    Device::ExternDevice::Instance &extern_device = static_cast<Device::ExternDevice::Instance &>(*(*it));

    extern_device.stepFailure(analysis);
  }
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::acceptStep
// Purpose       : Communicate to devices that the current step has been
//                 accepted
// Special Notes : Most devices need not know.  The Transmission line does.
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void CktLoader::acceptStep ()
{
  deviceManager_.acceptStep();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
bool CktLoader::getInitialQnorm (std::vector<TimeIntg::TwoLevelError> & tleVec )
{
  return deviceManager_.getInitialQnorm(tleVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getInnerLoopErrorSums
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool CktLoader::getInnerLoopErrorSums (std::vector<TimeIntg::TwoLevelError> & tleVec) const
{
  return deviceManager_.getInnerLoopErrorSums (tleVec);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::updateStateArrays
// Purpose       : Tells the inner loop solve to update state arrays.
// Special Notes : Needed to support voltlim with loca, with 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool CktLoader::updateStateArrays ()
{
  return deviceManager_.updateStateArrays ();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::startTimeStep()
// Purpose       : Tells the inner loop solve to do its prediction.
// Special Notes : Needed to support 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool CktLoader::startTimeStep(
  bool                          beginIntegrationFlag,
  double                        nextTimeStep,
  double                        nextTime,
  int                           currentOrder)
{
  return deviceManager_.startTimeStep(beginIntegrationFlag, nextTimeStep, nextTime, currentOrder);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::setExternalSolverState
// Purpose       :
// Special Notes : Needed to support 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/21/06
//-----------------------------------------------------------------------------
void CktLoader::setExternalSolverState(bool external_initJctFlag)
{
  deviceManager_.setExternalSolverState(external_initJctFlag);
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::getVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015
//---------------------------------------------------------------------------
bool CktLoader::getVoltageLimiterStatus()
{
  return deviceManager_.getVoltageLimiterStatus();
}

//-----------------------------------------------------------------------------
// Function      : CktLoader::setVoltageLimiterStatus
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 6/22/2015 
//---------------------------------------------------------------------------
void CktLoader::setVoltageLimiterStatus(bool voltageLimterStatus)
{
  return deviceManager_.setVoltageLimiterStatus(voltageLimterStatus);
}

//---------------------------------------------------------------------------
// Function      : CktLoader::updateDependentParams () 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 8/11/2020
//---------------------------------------------------------------------------
void CktLoader::updateDependentParams () 
{
  return deviceManager_.updateDependentParams ();
}

//---------------------------------------------------------------------------
// Function      : CktLoader::resetScaledParams() 
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 8/11/2020
//---------------------------------------------------------------------------
void CktLoader::resetScaledParams() 
{
  return deviceManager_.resetScaledParams();
}

} // namespace Loader
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2024 National Technology & Engineering Solutions of
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
// Purpose        : 
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 2015-07-11
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Standard Includes   ----------

#include <algorithm>
#include <sstream>
#include <stdexcept>

// ----------   Xyce Includes   ----------

#include <N_ANP_AnalysisManager.h>
#include <N_NLS_SensitivityResiduals.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LOA_NonlinearEquationLoader.h>

#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Manager.h>
#include <N_PDS_Serial.h>

#include <N_TIA_DataStore.h>

#include <N_UTL_Algorithm.h>
#include <N_UTL_Diagnostic.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExtendedString.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_MachDepParams.h>

// ----------   Static Declarations ----------

namespace Xyce {
namespace Nonlinear {

//-----------------------------------------------------------------------------
// Function      : setupOriginalParams
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/29/2016
//-----------------------------------------------------------------------------
bool setupOriginalParams ( TimeIntg::DataStore & ds,
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_
    )
{
  ds.paramOrigVals_.clear();

  std::vector<std::string>::const_iterator firstParam = paramNameVec_.begin ();
  std::vector<std::string>::const_iterator lastParam  = paramNameVec_.end ();
  std::vector<std::string>::const_iterator iterParam;

  for (iterParam=firstParam; iterParam!=lastParam; ++iterParam)
  {
    std::string paramName(*iterParam);

    // get the original value of this parameter, to be used for scaling and/or 
    // numerical derivatives later.
    double paramOrig = 0.0;
    std::string setParamName(paramName);
    getSetParamName(paramName, setParamName); // remove curly braces;
    bool found = nonlinearEquationLoader_.getParamAndReduce(
        analysisManager_.getPDSManager()->getPDSComm()->comm(),setParamName, paramOrig);

    if (!found)
    {
      Report::DevelFatal().in("Sensitivity::setupOriginalParams")
        << "cannot find parameter "
        << paramName;
    }
    ds.paramOrigVals_.push_back(paramOrig);
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : testForAnalyticDerivatives
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 2/29/2016
//-----------------------------------------------------------------------------
bool testForAnalyticDerivatives ( 
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_
    )
{
  bool allAvailable=true;

  std::vector<std::string>::const_iterator firstParam = paramNameVec_.begin ();
  std::vector<std::string>::const_iterator lastParam  = paramNameVec_.end ();
  std::vector<std::string>::const_iterator iterParam;

  for (iterParam=firstParam; iterParam!=lastParam; ++iterParam)
  {
    std::string paramName(*iterParam);
    std::string setParamName(paramName);
    getSetParamName(paramName, setParamName); // remove curly braces;

    if ( !(nonlinearEquationLoader_.analyticSensitivitiesAvailable (setParamName)) )
    {
      allAvailable=false;
    }
  }

  return allAvailable;
}


//-----------------------------------------------------------------------------
// Function      : slowNumericalDerivatives
// Purpose       : Numerical calculation of dfdp, dqdp, and dbdp.
//
//
// Special Notes : This is a "last resort" for computing nuemrical derivatives.
//
//                 The device package can compute most numerical derivatives, 
//                 and do it relatively efficiently.  However, there are some
//                 use cases not currently covered by the device pacakge numerical
//                 differentiation.  (such as global_params).  For those use
//                 cases, this function is necessary.
//
//                 It computes dfdp, etc via vector-wide finite differences.  As 
//                 most entries in dfdp, etc, will be "zero" for a typical device
//                 parameter sensitivity, this is mostly wasted effort.  Also,
//                 if this method is used, the sparsity pattern is not known, 
//                 which means that the sparse storage won't work.
//
//                 But, it is better than nothing.
//
//                 Note, this is the way all the derivatives were computed 
//                 in my first implementation of sensitivities, circa 2002.  The
//                 code is largely the same, just moved around and now mostly 
//                 not used.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/2018
//-----------------------------------------------------------------------------
bool slowNumericalDerivatives( int iparam,
  std::vector<std::string>::const_iterator & iterParam, 
  int difference, 
  double sqrtEta_, std::string & netlistFilename_,
  TimeIntg::DataStore & ds,
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_ 
    )
{
  int maxParamStringSize_=0;

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    Xyce::dout() << std::endl << "  Calculating numerical df/dp, dq/dp and db/dp for: ";
    Xyce::dout() << *iterParam << std::endl;

    for (int ipar=0;ipar<paramNameVec_.size();ipar++)
    {
      int sz = paramNameVec_[ipar].size();
      if (sz > maxParamStringSize_) { maxParamStringSize_ = sz; }
    }
  }

  std::string paramName(*iterParam);
  std::string setParamName(paramName);
  getSetParamName(paramName, setParamName); // remove curly braces;

  // save a copy of the DAE vectors
  Linear::Vector * origFVector = ds.daeFVectorPtr->cloneCopyVector();
  Linear::Vector * origQVector = ds.daeQVectorPtr->cloneCopyVector();
  Linear::Vector * origBVector = ds.daeBVectorPtr->cloneCopyVector();

  double paramOrig = ds.paramOrigVals_[iparam];

  // now perturb the value of this parameter.
  double paramPerturbed = paramOrig;
  //double dp = sqrtEta_ * (1.0 + fabs(paramOrig));  // old way

  double dp = sqrtEta_ * fabs( paramOrig );
  double minDouble = 10.0*Util::MachineDependentParams::DoubleMin();
  if (dp < minDouble)
  {
    dp = sqrtEta_;
  }

  if (difference==SENS_FWD)
  {
    paramPerturbed += dp;
  }
  else if (difference==SENS_REV)
  {
    paramPerturbed -= dp;
  }
  else if (difference==SENS_CNT)
  {
    Report::UserFatal0() << "difference=central not supported.";
  }
  else
  {
    Report::UserFatal0() << "difference not recognized!";
  }

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    Xyce::dout() << std::setw(maxParamStringSize_)<< *iterParam
      << " dp = " << std::setw(11)<< std::scientific<< std::setprecision(4) << dp 
      << " original value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramOrig 
      << " modified value = " << std::setw(16)<< std::scientific<< std::setprecision(9) << paramPerturbed 
      <<std::endl;
  }

  if (!(nonlinearEquationLoader_.setParam (setParamName, paramPerturbed) ) )
  {
    Report::DevelFatal().in("Sensitivity::slowNumericalDerivatives")
      << "cannot find parameter " << setParamName;
  }

  // Now that the parameter has been perturbed,
  // calculate the numerical derivative.

  // Load F,Q and B.
  nonlinearEquationLoader_.loadRHS();

  // save the perturbed DAE vectors
  Linear::Vector * pertFVector = ds.daeFVectorPtr->cloneCopyVector();
  Linear::Vector * pertQVector = ds.daeQVectorPtr->cloneCopyVector();
  Linear::Vector * pertBVector = ds.daeBVectorPtr->cloneCopyVector();

  Linear::MultiVector * dfdpPtrVector = ds.nextDfdpPtrVector;
  Linear::MultiVector * dqdpPtrVector = ds.nextDqdpPtrVector;
  Linear::MultiVector * dbdpPtrVector = ds.nextDbdpPtrVector;

  // calculate the df/dp vector.  
  double rdp=1/dp;
  Teuchos::RCP<Linear::Vector> dfdpPtr = Teuchos::rcp( dfdpPtrVector->getNonConstVectorView(iparam) );
  dfdpPtr->update( 1.0, *pertFVector, -1.0, *origFVector, 0.0 );
  dfdpPtr->scale(rdp);

  // calculate the dq/dp vector.  
  Teuchos::RCP<Linear::Vector> dqdpPtr = Teuchos::rcp( dqdpPtrVector->getNonConstVectorView(iparam) );
  dqdpPtr->update( 1.0, *pertQVector, -1.0, *origQVector, 0.0 );
  dqdpPtr->scale(rdp);

  // calculate the db/dp vector.  
  Teuchos::RCP<Linear::Vector> dbdpPtr = Teuchos::rcp( dbdpPtrVector->getNonConstVectorView(iparam) );
  dbdpPtr->update( 1.0, *pertBVector, -1.0, *origBVector, 0.0 );
  dbdpPtr->scale(rdp);

  if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
  {
    Xyce::dout() << *iterParam << ": ";
    Xyce::dout().width(15); Xyce::dout().precision(7); Xyce::dout().setf(std::ios::scientific);
    Xyce::dout() << "deviceSens_dp = " << dp << std::endl;

    int solutionSize_ = pertFVector->localLength();

    for (int k1 = 0; k1 < solutionSize_; ++k1)
    {

      Xyce::dout() 
        <<"fpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*pertFVector)[k1]
        <<" forig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*origFVector)[k1]
        <<" dfdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*dfdpPtr)[k1]
        <<std::endl;
    }

    Xyce::dout() << std::endl;
    for (int k1 = 0; k1 < solutionSize_; ++k1)
    {
      Xyce::dout() 
        <<"qpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*pertQVector)[k1]
        <<" qorig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*origQVector)[k1]
        <<" dqdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*dqdpPtr)[k1]
        <<std::endl;
    }

    Xyce::dout() << std::endl ;
    for (int k1 = 0; k1 < solutionSize_; ++k1)
    {
      Xyce::dout() 
        <<"bpert["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*pertBVector)[k1]
        <<" borig["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*origBVector)[k1]
        <<" dbdp ["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*dbdpPtr)[k1]
        <<std::endl;

    }

    std::ostringstream filename; 
    filename << netlistFilename_ << "_dfdp";
    filename << std::setw(3) << std::setfill('0') << iparam;
    filename << ".txt";
    dfdpPtr->writeToFile(const_cast<char *>(filename.str().c_str()));

    filename.str("");
    filename << netlistFilename_ << "_fpert";
    filename << std::setw(3) << std::setfill('0') << iparam;
    filename << ".txt";
    pertFVector->writeToFile(const_cast<char *>(filename.str().c_str()));

    filename.str("");
    filename << netlistFilename_ << "_dqdp";
    filename << std::setw(3) << std::setfill('0') << iparam;
    filename << ".txt";
    dqdpPtr->writeToFile(const_cast<char *>(filename.str().c_str()));

    filename.str("");
    filename << netlistFilename_ << "_qpert";
    filename << std::setw(3) << std::setfill('0') << iparam;
    filename << ".txt";
    pertQVector->writeToFile(const_cast<char *>(filename.str().c_str()));

    filename.str("");
    filename << netlistFilename_ << "_dbdp";
    filename << std::setw(3) << std::setfill('0') << iparam;
    filename << ".txt";
    dbdpPtr->writeToFile(const_cast<char *>(filename.str().c_str()));

    filename.str("");
    filename << netlistFilename_ << "_bpert";
    filename << std::setw(3) << std::setfill('0') << iparam;
    filename << ".txt";
    pertBVector->writeToFile(const_cast<char *>(filename.str().c_str()));
  }

  // now reset the parameter and rhs to previous values.
  if (!( nonlinearEquationLoader_.setParam (setParamName, paramOrig) ) )
  {
    Report::DevelFatal().in("Sensitivity::loadSensitivityResiduals")
      << "cannot find parameter " << setParamName;
  }

  *(ds.daeFVectorPtr) = *origFVector;
  *(ds.daeQVectorPtr) = *origQVector;
  *(ds.daeBVectorPtr) = *origBVector;

#if 0
  // this is needed for finite differences to work with transient adjoints!
  // See bug 1080 (SON).
  std::vector<int> FindicesVec;
  std::vector<int> QindicesVec;
  std::vector<int> BindicesVec;

  if (!ds.masterIndexVectorSize[iparam])
  {
    computeSparseIndices( iparam, ds, FindicesVec, QindicesVec, BindicesVec );
  }
#endif

  delete origFVector;
  delete origQVector;
  delete origBVector;
  delete pertFVector;
  delete pertQVector;
  delete pertBVector;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : loadSensitivityResiduals
//
// Purpose       : This function is to be called after a calculation has
//                 converged.  It computes the sensitivies of different components
//                 of the residual.  Recall the DAE form of the residual:
//                
//                 F = dq/dt + f - b
//
//                 This function sets up the following derivatives:  
//
//                 df/dp, dq/dp, db/dp where p is the parameter.
//
//                 Which can ultimately be assembled to give dF/dp:
//
//                         d(dq/dp)    
//                 dF/dp = ------- +  df/dp - db/dp
//                           dt
//
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/15/02
//-----------------------------------------------------------------------------
bool loadSensitivityResiduals (
  int difference, bool forceFD_, bool forceDeviceFD_, bool forceAnalytic_, 
  double sqrtEta_, std::string & netlistFilename_,
  TimeIntg::DataStore & ds,
  Loader::NonlinearEquationLoader & nonlinearEquationLoader_,
  const std::vector<std::string> & paramNameVec_,
  const Analysis::AnalysisManager & analysisManager_
    )
{
  Stats::StatTop _solveTransientAdjointStat("loadSensitivityResiduals");
  Stats::TimeBlock _solveTransientAdjointTimer(_solveTransientAdjointStat);

  int iparam=0;
  std::string msg;

  Linear::MultiVector * dfdpPtrVector = ds.nextDfdpPtrVector;
  Linear::MultiVector * dqdpPtrVector = ds.nextDqdpPtrVector;
  Linear::MultiVector * dbdpPtrVector = ds.nextDbdpPtrVector;

  setupOriginalParams (ds, nonlinearEquationLoader_, paramNameVec_, analysisManager_);

  // it is necessary to load the Jacobian here to make sure we have the most
  // up-to-date matrix.  The Jacobian is not loaded for the final 
  // evaluation of the residual in the Newton solve.
  //
  // Also, if any of the sensitivities rely on a numerical derivative, it is necessary to force 
  // the load procedure to NOT use the separate load option (where the linear and nonlinear 
  // portions are separated for efficiency).
  bool origSeparateLoadFlag = nonlinearEquationLoader_.getSeparateLoadFlag ();
  bool allAnalyticalAvailable = testForAnalyticDerivatives ( nonlinearEquationLoader_, paramNameVec_, analysisManager_);

  if (forceFD_ || !allAnalyticalAvailable) { nonlinearEquationLoader_.setSeparateLoadFlag (false); }
  nonlinearEquationLoader_.loadJacobian ();

  std::vector<std::string>::const_iterator firstParam = paramNameVec_.begin ();
  std::vector<std::string>::const_iterator lastParam  = paramNameVec_.end ();
  std::vector<std::string>::const_iterator iterParam;

  dfdpPtrVector->putScalar(0.0);
  dqdpPtrVector->putScalar(0.0);
  dbdpPtrVector->putScalar(0.0);

  for ( iterParam=firstParam, iparam=0;
        iterParam!=lastParam; ++iterParam, ++iparam )
  {
    std::string paramName(*iterParam);
    std::string setParamName(paramName);
    getSetParamName(paramName, setParamName); // remove curly braces;

    bool analyticAvailable = nonlinearEquationLoader_.analyticSensitivitiesAvailable (setParamName);;
    bool deviceLevelNumericalAvailable = nonlinearEquationLoader_.numericalSensitivitiesAvailable (setParamName);

    std::vector<double> dfdpVec;
    std::vector<double> dqdpVec;
    std::vector<double> dbdpVec;

    std::vector<int> FindicesVec;
    std::vector<int> QindicesVec;
    std::vector<int> BindicesVec;

    if (forceAnalytic_ && !analyticAvailable )
    {
      Report::UserError0() << "loadSensitivityResiduals: Analytical derivative was requested, but not available\n";
    }

    if (forceDeviceFD_ && !deviceLevelNumericalAvailable)
    {
      Report::UserError0() << "loadSensitivityResiduals: Device level numerical derivative was requested, but not available\n";
    }

    // do old slow numerical derivative if user is forcing it, or if we have no other choice
    if (forceFD_ || (!analyticAvailable && !deviceLevelNumericalAvailable) )  
    {
      slowNumericalDerivatives( iparam, iterParam,
          difference, sqrtEta_, netlistFilename_, 
          ds, nonlinearEquationLoader_, paramNameVec_, analysisManager_);
    }
    else
    {
      if (analyticAvailable && !forceDeviceFD_)
      {
        nonlinearEquationLoader_.getAnalyticSensitivities(setParamName,
            dfdpVec,dqdpVec,dbdpVec,
            FindicesVec, QindicesVec, BindicesVec);
      }
      else if (deviceLevelNumericalAvailable && !forceAnalytic_)
      {
        nonlinearEquationLoader_.getNumericalSensitivities(setParamName,
            dfdpVec,dqdpVec,dbdpVec,
            FindicesVec, QindicesVec, BindicesVec);
      }
      else
      {
        Report::UserError0() << "loadSensitivityResiduals: ERROR, can't compute dfdp,dqdp,dbdp \n";
      }

      Teuchos::RCP<Linear::Vector> dfdpPtr = Teuchos::rcp( dfdpPtrVector->getNonConstVectorView(iparam) );
      Teuchos::RCP<Linear::Vector> dqdpPtr = Teuchos::rcp( dqdpPtrVector->getNonConstVectorView(iparam) );
      Teuchos::RCP<Linear::Vector> dbdpPtr = Teuchos::rcp( dbdpPtrVector->getNonConstVectorView(iparam) );

      // ADMS devices with nodes that collapse to ground can have "-1"
      // as their LIDs, meaning the node doesn't really exist
      {
        int Fsize=FindicesVec.size();
        for (int i=0;i<Fsize;++i)
        {
          if (FindicesVec[i] != -1)
            (*dfdpPtr)[FindicesVec[i]]  += dfdpVec[i];
        }

        int Qsize=QindicesVec.size();
        for (int i=0;i<Qsize;++i)
        {
          if (QindicesVec[i] != -1)
            (*dqdpPtr)[QindicesVec[i]]  += dqdpVec[i];
        }

        int Bsize=BindicesVec.size();
        for (int i=0;i<Bsize;++i)
        {
          if (BindicesVec[i] != -1)
            (*dbdpPtr)[BindicesVec[i]]  += dbdpVec[i];
        }
      }

      dfdpPtr->fillComplete();
      dqdpPtr->fillComplete();
      dbdpPtr->fillComplete();

      if (!(ds.masterIndexVectorSize.empty())) // if empty, then adjoints not requested
      {
        if (!ds.masterIndexVectorSize[iparam])  // check if these have been set up yet
        {
          computeSparseIndices( iparam, ds, FindicesVec, QindicesVec, BindicesVec );
        }
      }

      if (DEBUG_NONLINEAR && isActive(Diag::SENS_SOLVER))
      {
        Xyce::dout() << *iterParam << ": " <<std::endl;
        for (int k1 = 0; k1 < dfdpPtrVector->localLength(); ++k1)
        {
          Xyce::dout() 
            <<"dfdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*dfdpPtr)[k1]
            <<" dqdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*dqdpPtr)[k1]
            <<" dbdp["<<std::setw(3)<<k1<<"]= "<<std::setw(15)<<std::scientific<<std::setprecision(8)<<(*dbdpPtr)[k1]
            <<std::endl;
        }
      }
    }

  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : computeSparseIndices
//
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 2/2/2017
//-----------------------------------------------------------------------------
bool computeSparseIndices ( const int iparam, 
  TimeIntg::DataStore & ds,     
  const std::vector<int>& FindicesVec,
  const std::vector<int>& QindicesVec,
  const std::vector<int>& BindicesVec
    )
{
  std::vector<int> & masterRef = (ds.masterIndexVector [iparam]);

  int Fsize = FindicesVec.size();
  int Qsize = QindicesVec.size();
  int Bsize = BindicesVec.size();
  int size = std::max(Bsize,std::max(Fsize,Qsize));

  Linear::MultiVector * dfdpPtrVector = ds.nextDfdpPtrVector;
  int solutionSize = dfdpPtrVector->localLength();
  int gndNode = dfdpPtrVector->omap()->numLocalEntities() - 1;

  // Get parallel information.     
  int numProc =  dfdpPtrVector->pdsComm()->numProc();
  int myPID = dfdpPtrVector->pdsComm()->procID();

  if (Fsize >0 && Fsize != size)
  {
    std::cout << "THIS IS BROKEN  Fsize != size !!!!!" <<std::endl;
    exit(0);
  }

  if (Qsize >0 && Qsize != size)
  {
    std::cout << "THIS IS BROKEN  Qsize != size !!!!!" <<std::endl;
    exit(0);
  }

  if (Bsize >0 && Bsize != size)
  {
    std::cout << "THIS IS BROKEN  Bsize != size !!!!!" <<std::endl;
    exit(0);
  }

  std::vector<int> ghostNodes;
  for (int i=0;i<Fsize;++i)
  {
    int currIdx = FindicesVec[i];
    if (currIdx != -1 && currIdx != gndNode)
    {
      if (currIdx < solutionSize)
      {
        masterRef.push_back( currIdx );
      }
      else
      { 
        ghostNodes.push_back( dfdpPtrVector->omap()->localToGlobalIndex( currIdx ) );
      }
    }
  }

  for (int i=0;i<Qsize;++i)
  {
    int currIdx = QindicesVec[i];
    if (currIdx != -1 && currIdx != gndNode)
    {
      if (currIdx < solutionSize)
      {
        masterRef.push_back( currIdx );
      }
      else
      {
        ghostNodes.push_back( dfdpPtrVector->omap()->localToGlobalIndex( currIdx ) );
      }
    }
  }

  for (int i=0;i<Bsize;++i)
  {
    int currIdx = BindicesVec[i];
    if (currIdx != -1 && currIdx != gndNode)
    {
      if (currIdx < solutionSize)
      {
        masterRef.push_back( currIdx );
      }
      else
      {
        ghostNodes.push_back( dfdpPtrVector->omap()->localToGlobalIndex( currIdx ) );
      }
    }
  }

  // Create unique list of masterRef indices. 
  std::sort( masterRef.begin(), masterRef.end() ); 
  masterRef.erase( std::unique( masterRef.begin(), masterRef.end() ), masterRef.end() );

  // Process ghost nodes, if there are any.
  if (numProc > 1)
  {
    // Create unique list of ghost nodes.
    std::sort( ghostNodes.begin(), ghostNodes.end() );
    ghostNodes.erase( std::unique( ghostNodes.begin(), ghostNodes.end() ), ghostNodes.end() );

    // First communicate number of ghost nodes to all processors.
    std::vector<int> myGhostNodes( numProc, 0 ), totalGhostNodes( numProc, 0 );
    myGhostNodes[ myPID ] = ghostNodes.size();
    dfdpPtrVector->pdsComm()->sumAll( &myGhostNodes[0], &totalGhostNodes[0], numProc );

    // Communicate list of ghost nodes.
    int myStartIdx = 0, totalIdx = 0;
    for (int i=0; i<numProc; ++i)
    {
      if (i < myPID)
      {
        myStartIdx += totalGhostNodes[i];
      }
      totalIdx += totalGhostNodes[i];
    }
    std::vector<int> myGhostNodesGIDs( totalIdx, 0 ), totalGhostNodesGIDs( totalIdx, 0 );

    for (unsigned int i=0; i<ghostNodes.size(); ++i)
    {
      myGhostNodesGIDs[ myStartIdx + i ] = ghostNodes[i];
    }
    dfdpPtrVector->pdsComm()->sumAll( &myGhostNodesGIDs[0], &totalGhostNodesGIDs[0], totalIdx );
    std::sort( totalGhostNodesGIDs.begin(), totalGhostNodesGIDs.end() );
    totalGhostNodesGIDs.erase( std::unique( totalGhostNodesGIDs.begin(), totalGhostNodesGIDs.end() ), totalGhostNodesGIDs.end() );

    // Determine if masterRef should be ammended.
    for (unsigned int i=0; i< totalGhostNodesGIDs.size(); ++i)
    {
      int lid = dfdpPtrVector->pmap()->globalToLocalIndex( totalGhostNodesGIDs[i] ); 
      if (lid > -1)
      {
        std::vector<int>::iterator iter = std::find( masterRef.begin(), masterRef.end(), lid );
        if (iter == masterRef.end())
        {
          masterRef.push_back( lid );
        }
      }
    }

    // Collect the global master ref size.
    int tmpMRSize = masterRef.size(), mrSize = 0;
    dfdpPtrVector->pdsComm()->sumAll( &tmpMRSize, &mrSize, 1 );
    ds.masterIndexVectorSize[iparam] = mrSize;
  }
  else
  {
    ds.masterIndexVectorSize[iparam] = masterRef.size();
  } 

  return true;
}

} // namespace Nonlinear
} // namespace Xyce

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
// Purpose        : MPDE Specific Loader
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/11/04
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef Xyce_MPDE_Loader_H
#define Xyce_MPDE_Loader_H

// ---------- Standard Includes ----------

#include <vector>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_LAS_fwd.h>
#include <N_LOA_Loader.h>

#include <N_MPDE_State.h>
#include <N_MPDE_WarpedPhaseCondition.h>
#include <N_DEV_DeviceMgr.h>

// ---------- Forward declarations --------

class N_MPDE_Discretization;
class N_MPDE_Manager;

//-----------------------------------------------------------------------------
// Class         : N_MPDE_Loader
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_MPDE_Loader : public Xyce::Loader::Loader
{
public:
  N_MPDE_Loader(
    N_MPDE_State &                      state,
    Xyce::Device::DeviceMgr &           device_manager,
    Xyce::Loader::Loader &              loader,
    const N_MPDE_Discretization &       disc,
    const N_MPDE_Manager &              mgr,
    const N_MPDE_WarpedPhaseCondition * warpMPDEPhasePtr)
  : state_(state),
    loader_(loader),
    deviceManager_(device_manager),
    fastTimeDiscretizer_(disc),
    periodicTimesOffset_(0),
    period_(0),
    mpdeManager_(mgr),
    bOmegadQdt2Ptr_(0),
    warpMPDEPhasePtr_(warpMPDEPhasePtr),
    allDevicesAllTimePointsConverged_(true)
  {}

  // Destructor
  ~N_MPDE_Loader();

  bool allDevicesConverged(Xyce::Parallel::Machine comm);

  // Method which is called to load the new-DAE contributions 
  bool loadDAEMatrices( Xyce::Linear::Vector * X,
                        Xyce::Linear::Vector * S,
                        Xyce::Linear::Vector * dSdt,
                        Xyce::Linear::Vector * Store,
                        Xyce::Linear::Matrix * dQdx,
                        Xyce::Linear::Matrix * dFdx,
                        int loadType = Xyce::Device::ALL );

  // Method which is called to load the new-DAE vectors
  bool loadDAEVectors( Xyce::Linear::Vector * X,
                       Xyce::Linear::Vector * currX,
                       Xyce::Linear::Vector * lastX,
                       Xyce::Linear::Vector * S,
                       Xyce::Linear::Vector * currS,
                       Xyce::Linear::Vector * lastS,
                       Xyce::Linear::Vector * dSdt,
                       Xyce::Linear::Vector * Store,
                       Xyce::Linear::Vector * currStore,
                       Xyce::Linear::Vector * lastStore,
                       Xyce::Linear::Vector * nextLeadFVectorPtr,
                       Xyce::Linear::Vector * nextLeadQVectorPtr,
                       Xyce::Linear::Vector * nextJunctionVVectorPtr,
                       Xyce::Linear::Vector * Q,
                       Xyce::Linear::Vector * F,
                       Xyce::Linear::Vector * B,
                       Xyce::Linear::Vector * dFdxdVp,
                       Xyce::Linear::Vector * dQdxdVp,
                       int loadType = Xyce::Device::ALL );


  bool updateState(    Xyce::Linear::Vector * nextSolVectorPtr,
                       Xyce::Linear::Vector * currSolVectorPtr,
                       Xyce::Linear::Vector * lastSolVectorPtr,
                       Xyce::Linear::Vector * nextStaVectorPtr,
                       Xyce::Linear::Vector * currStaVectorPtr,
                       Xyce::Linear::Vector * lastStaVectorPtr,
                       Xyce::Linear::Vector * nextStoVectorPtr,
                       Xyce::Linear::Vector * currStoVectorPtr,
                       Xyce::Linear::Vector * lastStoVectorPtr,
                       int loadType = Xyce::Device::ALL );

  // Virtual method which initializes the nonlinear problem.
  virtual bool initializeProblem( Xyce::Linear::Vector * nextSolVectorPtr,
                          Xyce::Linear::Vector * currSolVectorPtr,
                          Xyce::Linear::Vector * lastSolVectorPtr,
                          Xyce::Linear::Vector * nextStaVectorPtr,
                          Xyce::Linear::Vector * currStaVectorPtr,
                          Xyce::Linear::Vector * lastStaVectorPtr,
                          Xyce::Linear::Vector * StateDerivVectorPtr,
                          Xyce::Linear::Vector * nextStoVectorPtr,
                          Xyce::Linear::Vector * currStoVectorPtr,
                          Xyce::Linear::Vector * lastStoVectorPtr,
                          Xyce::Linear::Vector * QVectorPtr,
                          Xyce::Linear::Vector * FVectorPtr,
                          Xyce::Linear::Vector * BVectorPtr,
                          Xyce::Linear::Vector * dFdxdVpVectorPtr,
                          Xyce::Linear::Vector * dQdxdVpVectorPtr) const  { return false; }

  // Assign times for fast time scale
  void setFastTimes( const std::vector<double> & times );
  
  void setPeriodFlags( const std::vector<bool> & periodicFlags );
  
  // Construct a periodic version of times_
  void constructPeriodicTimes();

  // Registration method for the next Application Vector
  void registerAppNextVec( RCP<Xyce::Linear::Vector> appNextVecPtr );

  // Registration method for the curr Application Vector
  void registerAppCurrVec( RCP<Xyce::Linear::Vector> appCurrVecPtr );

  // Registration method for the last Application Vector
  void registerAppLastVec( RCP<Xyce::Linear::Vector> appLastVecPtr );

  // Registration method for the next Application State Vector
  void registerAppNextStaVec( RCP<Xyce::Linear::Vector> appNextStaVecPtr );

  // Registration method for the curr Application State Vector
  void registerAppCurrStaVec( RCP<Xyce::Linear::Vector> appCurrStaVecPtr );

  // Registration method for the last Application State Vector
  void registerAppLastStaVec( RCP<Xyce::Linear::Vector> appLastStaVecPtr );

  // Registration method for the next Application Store Vector
  void registerAppNextStoVec( RCP<Xyce::Linear::Vector> appNextStoVecPtr );

  // Registration method for the curr Application Store Vector
  void registerAppCurrStoVec( RCP<Xyce::Linear::Vector> appCurrStoVecPtr );

  // Registration method for the curr Application Store Vector
  void registerAppLastStoVec( RCP<Xyce::Linear::Vector> appLastStoVecPtr );

  // Registration method for the Application Store Vector Lead Current Q Component 
  void registerAppNextLeadCurrentVec( RCP<Xyce::Linear::Vector> aVec );
  void registerAppLeadCurrentQVec( RCP<Xyce::Linear::Vector> aVec );
  void registerAppNextJunctionVVec( RCP<Xyce::Linear::Vector> aVec );
  
  // Registration method for the Application dQdx Matrix
  void registerAppdQdx( RCP<Xyce::Linear::Matrix> appdQdxPtr );

  // Registration method for the Application dFdx Matrix
  void registerAppdFdx( RCP<Xyce::Linear::Matrix> appdFdxPtr );

  // Registration method for the MPDE size dQdx Block Matrix
  void registerMPDEdQdx( RCP<Xyce::Linear::BlockMatrix> bmdQdxPtr );

  // Registration method for the MPDE size dFdx Block Matrix
  void registerMPDEdFdx( RCP<Xyce::Linear::BlockMatrix> bmdFdxPtr );

  // Registration method for a temp vector in constructin dq/dt2
  void registerOmegadQdt2( RCP<Xyce::Linear::BlockVector> omegadQdt2Ptr);

  void getAnalyticSensitivities(
      std::string & name, 
      std::vector<double> & dfdpVec, 
      std::vector<double> & dqdpVec,
      std::vector<double> & dbdpVec,
      std::vector<int> & FindicesVec,
      std::vector<int> & QindicesVec,
      std::vector<int> & BindicesVec) const
  {}

  bool analyticSensitivitiesAvailable (std::string & name) { return false; }

  virtual bool setParam (std::string & name, double val, bool overrideOriginal=false) { return false; }
  virtual double getParamAndReduce(Xyce::Parallel::Machine comm, const std::string & name) const { return 0.0; }

  // voltage limiter toggle functions (do nothing for now)
  bool getVoltageLimiterStatus() {return true;}
  void setVoltageLimiterStatus(bool voltageLimterStatus) {}

private :
  N_MPDE_State &                state_;                         /// MPDE State
  Xyce::Loader::Loader &        loader_;                  /// Base Application loader
  Xyce::Device::DeviceMgr &     deviceManager_;

  // discretization
  const N_MPDE_Discretization & fastTimeDiscretizer_;

  //Fast Time Scale Points
  std::vector<double> times_;
  int periodicTimesOffset_;
  std::vector<double> periodicTimes_;
  double period_;
  
  //a vector of bools to indicate if a solution variable does not
  //appear to be periodic.  Non-periodic signals will be handled
  //differently by the MPDE_Loader class
  std::vector<bool> nonPeriodic_;
  

  // Application Linear Objects
  RCP<Xyce::Linear::Vector> appNextVecPtr_;
  RCP<Xyce::Linear::Vector> appCurrVecPtr_;
  RCP<Xyce::Linear::Vector> appLastVecPtr_;

  RCP<Xyce::Linear::Vector> appNextStaVecPtr_;
  RCP<Xyce::Linear::Vector> appCurrStaVecPtr_;
  RCP<Xyce::Linear::Vector> appLastStaVecPtr_;
  RCP<Xyce::Linear::Matrix> appdQdxPtr_;
  RCP<Xyce::Linear::Matrix> appdFdxPtr_;
  RCP<Xyce::Linear::Vector> appNextStoVecPtr_;
  RCP<Xyce::Linear::Vector> appCurrStoVecPtr_;
  RCP<Xyce::Linear::Vector> appLastStoVecPtr_;
  
  Teuchos::RCP<Xyce::Linear::Vector> appNextLeadFVecPtr_;
  Teuchos::RCP<Xyce::Linear::Vector> appLeadQVecPtr_;
  Teuchos::RCP<Xyce::Linear::Vector> appNextJunctionVVecPtr_;

  // Tmp storage block matrices 
  RCP<Xyce::Linear::BlockMatrix> bmdQdxPtr_;
  RCP<Xyce::Linear::BlockMatrix> bmdFdxPtr_;

  // mpde manager:
  const N_MPDE_Manager &        mpdeManager_;

  // MPDE/WaMPDE data
  RCP<Xyce::Linear::BlockVector> bOmegadQdt2Ptr_;
  
  const N_MPDE_WarpedPhaseCondition * warpMPDEPhasePtr_;
  
  //bool N_MPDE_Loader::findNonPeriodicSignals_(const Xyce::Linear::BlockVector & solutionBlockVector );

  bool allDevicesAllTimePointsConverged_;
};


// //-----------------------------------------------------------------------------
// // Function      : N_MPDE_Loader::registerAppLoader
// // Purpose       : Registration method for the device packaage
// // Special Notes :
// // Scope         : public
// // Creator       : Robert Hoekstra, 9233, Computational Sciences
// // Creation Date : 03/12/04
// //-----------------------------------------------------------------------------
// inline void N_MPDE_Loader::registerAppLoader( Xyce::Loader::Loader *appLoaderPtr )
// { 
//   loader_ = appLoaderPtr; 
// }

// //-----------------------------------------------------------------------------
// // Function      : N_MPDE_Loader::registerMPDEDeviceInterface
// // Purpose       : Registration method for the device packaage
// // Special Notes :
// // Scope         : public
// // Creator       : Robert Hoekstra, 9233, Computational Sciences
// // Creation Date : 03/12/04
// //-----------------------------------------------------------------------------
// inline void N_MPDE_Loader::registerMPDEDeviceInterface( Xyce::Device::DeviceMgr * mpdeDevIntPtr )
// { 
//   deviceManager_ = mpdeDevIntPtr;
// }

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppNextVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader:: registerAppNextVec( RCP<Xyce::Linear::Vector> appNextVecPtr )
{ 
  appNextVecPtr_ = appNextVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 1437
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader:: registerAppCurrVec( RCP<Xyce::Linear::Vector> appCurrVecPtr )
{ 
  appCurrVecPtr_ = appCurrVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppLastVec
// Purpose       : Registration method for the Application Vector
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 1437
// Creation Date :
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader:: registerAppLastVec( RCP<Xyce::Linear::Vector> appLastVecPtr )
{ 
  appLastVecPtr_ = appLastVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppNextStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppNextStaVec( RCP<Xyce::Linear::Vector> appNextStaVecPtr )
{ 
  appNextStaVecPtr_ = appNextStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppCurrStaVec( RCP<Xyce::Linear::Vector> appCurrStaVecPtr )
{ 
  appCurrStaVecPtr_ = appCurrStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppLastStaVec
// Purpose       : Registration method for the Application State Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 01/25/07
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppLastStaVec( RCP<Xyce::Linear::Vector> appLastStaVecPtr )
{ 
  appLastStaVecPtr_ = appLastStaVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppNextStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : 
// Creation Date :
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppNextStoVec( RCP<Xyce::Linear::Vector> appNextStoVecPtr )
{ 
  appNextStoVecPtr_ = appNextStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppCurrStoVec( RCP<Xyce::Linear::Vector> appCurrStoVecPtr )
{ 
  appCurrStoVecPtr_ = appCurrStoVecPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppCurrStoVec
// Purpose       : Registration method for the Application Store Vector
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppLastStoVec( RCP<Xyce::Linear::Vector> appLastStoVecPtr )
{ 
  appLastStoVecPtr_ = appLastStoVecPtr; 
}

inline  void N_MPDE_Loader::registerAppNextLeadCurrentVec( RCP<Xyce::Linear::Vector> aVec )
{ 
  appNextLeadFVecPtr_ = aVec; 
}
                          
inline  void N_MPDE_Loader::registerAppLeadCurrentQVec( RCP<Xyce::Linear::Vector> aVec )
{ 
  appLeadQVecPtr_ = aVec; 
}
inline  void N_MPDE_Loader::registerAppNextJunctionVVec( RCP<Xyce::Linear::Vector> aVec )
{ 
  appNextJunctionVVecPtr_ = aVec; 
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerOmegadQdt2
// Purpose       : Registration method for a temp vector in constructin dq/dt2
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414, Computational Sciences
// Creation Date : 08/10/05
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerOmegadQdt2( RCP<Xyce::Linear::BlockVector> omegadQdt2Ptr)
{ 
  bOmegadQdt2Ptr_ = omegadQdt2Ptr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppdQdx
// Purpose       : Registration method for the Application dQdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppdQdx( RCP<Xyce::Linear::Matrix> appdQdxPtr )
{ 
  appdQdxPtr_ = appdQdxPtr; 
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Loader::registerAppdFdx
// Purpose       :  Registration method for the Application dFdx Matrix
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
inline void N_MPDE_Loader::registerAppdFdx( RCP<Xyce::Linear::Matrix> appdFdxPtr )
{ 
  appdFdxPtr_ = appdFdxPtr; 
}

// Registration method for the MPDE size dQdx Block Matrix
inline void N_MPDE_Loader::registerMPDEdQdx( RCP<Xyce::Linear::BlockMatrix> bmdQdxPtr )
{ 
  bmdQdxPtr_ = bmdQdxPtr; 
}

// Registration method for the MPDE size dFdx Block Matrix
inline void N_MPDE_Loader::registerMPDEdFdx( RCP<Xyce::Linear::BlockMatrix> bmdFdxPtr )
{ 
  bmdFdxPtr_ = bmdFdxPtr; 
}

#endif


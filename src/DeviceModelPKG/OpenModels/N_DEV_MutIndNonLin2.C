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

//-------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Rich Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/21/2005
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_MutIndNonLin2.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_MutIndNonLin.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_HspiceBools.h>

using Teuchos::rcp;

namespace Xyce {
namespace Device {


namespace MutIndNonLin2 {


void Traits::loadInstanceParameters(ParametricData<MutIndNonLin2::Instance> &p)
{
  p.addPar ("COUP_VAL",1.0,&MutIndNonLin2::Instance::mutualCup)
   .setGivenMember(&MutIndNonLin2::Instance::mutualCupGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Coupling coefficient");

  p.addPar ("NONLINEARCOUPLING",0.0,&MutIndNonLin2::Instance::nonlinFlag)
   .setGivenMember(&MutIndNonLin2::Instance::nonlinFlagGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Nonlinear coupling flag");

  p.addPar ("COUPLEDMutIndNonLin",std::vector<std::string>(),&MutIndNonLin2::Instance::inductorNames)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("" );

  p.addPar ("COUPLEDINDUCTANCE",std::vector<double>(),&MutIndNonLin2::Instance::inductorInductances)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("NODE1",std::vector<std::string>(),&MutIndNonLin2::Instance::inductorsNode1)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("NODE2",std::vector<std::string>(),&MutIndNonLin2::Instance::inductorsNode2)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");

  p.addPar ("COUPLING",std::vector<double>(),&MutIndNonLin2::Instance::couplingCoefficient)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Coupling coefficient");

  p.addPar ("COUPLEDINDUCTOR",std::vector<std::string>(),&MutIndNonLin2::Instance::couplingInductor)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("");
  
  p.addPar ("IC",std::vector<double>(),&MutIndNonLin2::Instance::initialCondition)
   .setUnit(U_AMP)
   .setCategory(CAT_NONE)
   .setDescription("Initial current through the inductor.");
}

void Traits::loadModelParameters(ParametricData<MutIndNonLin2::Model> &p)
{
  p.addPar ("A",1000.0,&MutIndNonLin2::Model::A)
 .setUnit(U_AMPMM1)
   .setCategory(CAT_MATERIAL)
   .setDescription("Thermal energy parameter");

  p.addPar ("AREA",0.1,&MutIndNonLin2::Model::AreaInCm2)
 .setUnit(U_CM2)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Mean magnetic cross-sectional area");

  p.addPar ("ALPHA",5.0e-5,&MutIndNonLin2::Model::Alpha)
 .setUnit(U_NONE)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Domain coupling parameter");

  p.addPar ("BETAH",0.0001,&MutIndNonLin2::Model::BetaH)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Modeling constant");

  p.addPar ("BETAM",3.125e-5,&MutIndNonLin2::Model::BetaM)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Modeling constant");

  p.addPar ("C",0.2,&MutIndNonLin2::Model::C)
 .setUnit(U_NONE)
   .setCategory(CAT_MATERIAL)
   .setDescription("Domain flesing parameter"); 
  
  /* 
  p.addPar("CLIM",0.005,&MutIndNonLin::Model::CLim)
  .setUnit(U_NONE)
  .setCategory(CAT_MATERIAL)
  .setDescription("Value below which domain flexing parameter will be treated as zero.");
  */
  
  /*
  p.addPar("DELVSCALING",1.0e3,&MutIndNonLin::Model::DeltaVScaling)
  .setUnit(U_VOLT)
  .setCategory(CAT_NONE)
  .setDescription("Smoothing coefficient for voltage difference over first inductor");
  */
  
  p.addPar ("DELV",0.1,&MutIndNonLin2::Model::DeltaV)
  .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Smoothing coefficient for voltage difference over first inductor");
   
  /*
  p.addPar("CONSTDELVSCALING",false,&MutIndNonLin::Model::UseConstantDeltaVScaling)
  .setUnit(U_VOLT)
  .setCategory(CAT_NONE)
  .setDescription("Use constant scaling factor to smooth voltage difference over first inductor");

  p.addPar("INCLUDEMEQU",true,&MutIndNonLin::Model::includeMEquation)
  .setGivenMember(&MutIndNonLin::Model::includeMEquationGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Flag to include the magnetics in the solution.");
  */
  
  p.addPar ("GAP",0.0,&MutIndNonLin2::Model::GapInCm)
 .setUnit(U_CM)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Effective air gap");

  p.addPar ("K",500.0,&MutIndNonLin2::Model::Kirr)
 .setUnit(U_AMPMM1)
   .setCategory(CAT_MATERIAL)
   .setDescription("Domain anisotropy parameter");

  p.addPar ("KIRR",500.0,&MutIndNonLin2::Model::Kirr)
 .setUnit(U_AMPMM1)
   .setCategory(CAT_MATERIAL)
   .setDescription("Domain anisotropy parameter");

  p.addPar ("MS",1.0e+6,&MutIndNonLin2::Model::Ms)
 .setUnit(U_AMPMM1)
   .setCategory(CAT_MATERIAL)
   .setDescription("Saturation magnetization");

  p.addPar ("LEVEL",0.0,&MutIndNonLin2::Model::LevelIgnored)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("for pspice compatibility -- ignored");

  p.addPar ("PACK",0.0,&MutIndNonLin2::Model::PackIgnored)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("for pspice compatibility -- ignored");

  p.addPar ("PATH",1.0,&MutIndNonLin2::Model::PathInCm)
 .setUnit(U_CM)
   .setCategory(CAT_GEOMETRY)
   .setDescription("Total mean magnetic path");

  /* not part of level 1 inductor */
  p.addPar ("VINF",1.0,&MutIndNonLin2::Model::Vinf)
 .setUnit(U_VOLT)
   .setCategory(CAT_NONE)
   .setDescription("Smoothing coefficient for voltage difference over first inductor");

  p.addPar ("TNOM",27.0,&MutIndNonLin2::Model::tnom)
 .setUnit(U_DEGC)
   .setCategory(CAT_MATERIAL)
   .setDescription("Reference temperature");

  p.addPar ("TC1",0.0,&MutIndNonLin2::Model::tempCoeff1)
   .setUnit(U_NONE)
   .setCategory(CAT_MATERIAL)
   .setDescription("First order temperature coeff.");

  p.addPar ("TC2",0.0,&MutIndNonLin2::Model::tempCoeff2)
   .setUnit(U_NONE)
   .setCategory(CAT_MATERIAL)
   .setDescription("Second order temperature coeff.");

  p.addPar ("PZEROTOL",0.1,&MutIndNonLin2::Model::pZeroTol)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Tolerance for nonlinear zero crossing");

  p.addPar ("MVARSCALING",1.0,&MutIndNonLin2::Model::mVarScaling)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("M-variable scaling.");

  p.addPar ("RVARSCALING",1.0,&MutIndNonLin2::Model::rVarScaling)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("R-variable scaling");

  p.addPar ("MEQNSCALING",1.0,&MutIndNonLin2::Model::mEqScaling)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("M-equation scaling");

  p.addPar ("REQNSCALING",1.0,&MutIndNonLin2::Model::rEqScaling)
 .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("R-equation scaling");

  p.addPar ("OUTPUTSTATEVARS",0,&MutIndNonLin2::Model::outputStateVars)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Flag to save state variables" );
  
  /* not part of level 1 inductor */
  p.addPar ("INCLUDEDELTAM",0,&MutIndNonLin2::Model::includeDeltaM)
   .setGivenMember(&MutIndNonLin2::Model::includeDeltaMGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Flag to make M calculation implicit" );

  /* not part of level 1 inductor */
  p.addPar ("USERKINTEGRATION",0,&MutIndNonLin2::Model::useRKIntegration)
   .setGivenMember(&MutIndNonLin2::Model::useRKIntegrationGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Flag to use 4th order Runge-Kutta integration for dM/dH" );

  /* not part of level 1 inductor */
  p.addPar ("USESTATEDERIV",0,&MutIndNonLin2::Model::useStateDeriv)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Flag to use state vector for derivatives" );

  /* not part of level 1 inductor */
  p.addPar ("VOLTAGELIMITERFLAG",0,&MutIndNonLin2::Model::voltageLimiterFlag)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Flag to use voltage limiting on Mag and R internal variables" );

  /* not part of level 1 inductor */
  p.addPar ("MAGLIMITTHRES",0.1,&MutIndNonLin2::Model::magLimitThres)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Threshold over which newton interation changes in Mag are limited." );

  /* not part of level 1 inductor */
  p.addPar ("RLIMITTHRES",0.1,&MutIndNonLin2::Model::rLimitThres)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Threshold over which newton interation changes in R are limited." );

  p.addPar ("FACTORMS",0,&MutIndNonLin2::Model::factorMS)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Flag to save state variables" );
   
  p.addPar("BHSIUNITS",0,&MutIndNonLin2::Model::BHSiUnits)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Flag to report B and H in SI units");
}


//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Iiter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Iiter),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    Happ(0.0),
    branchCurrentSum(0.0),
    P(0.0),
    PPreviousStep(0.0),
    dP_dM(0.0),
    dP_dBranchCurrentSum(0.0),
    dP_dV1Pos(0.0),
    dP_dV1Neg(0.0),
    mEquFval(0.0),
    MagVar(0.0),
    oldBranchCurrentSum(0.0),
    MagVarUpdate(0.0),
    lastMagUpdate(0.0),
    useRKIntegration(false),
    includeDeltaM(false),
    lastH(0.0),
    outputStateVarsFlag( false ),
    dMdtAverage_(0.0)
{
  if (DEBUG_DEVICE)
  {
    Xyce::dout() << "In Instance constructor" << std::endl;
  }

  // for a simple case of 2 leads, we have 3 internal vars (I_branch, H, M)
  // and one state var (I_branch)
  numExtVars   = 2;
  numIntVars   = 3;
  numStateVars = 2;
  setNumStoreVars(4);
  tempGiven    = false;

  const int ibev = IB.numExtVars;
  const int ibiv = IB.numIntVars;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // these two vectors are used for 4th order RK
  // the vectors are one shorted than one would think you would need because
  // current step is held by local variables (thus, this just data from steps
  // n-1, n-2 and n-3
  branchCurrentSumHistory.resize(3);
  PFunctionHistory.resize(3);

  // if the model card askes for the delta M calculation to be implicit, then
  // change the includeDeltaM flag
  if( model_.includeDeltaMGiven && (model_.includeDeltaM > 0))
  {
    includeDeltaM = true;
  }

  if( model_.useRKIntegrationGiven && (model_.useRKIntegration > 0))
  {
    useRKIntegration = true;
  }

  // now load the instance data vector
  for( int i=0; i<inductorNames.size(); ++i )
  {
    InductorInstanceData * inductorData = new InductorInstanceData();
    inductorData->name = inductorNames[i];
    inductorData->L = inductorInductances[i];
    inductorData->baseL = inductorInductances[i];
    inductorData->ICGiven = false;
    inductorData->inductorCurrentOffsets.resize( inductorNames.size() );

    instanceData.push_back( inductorData );
  }
  numInductors = instanceData.size();

  // Set-up for power calculations.  We allocate space for all of the
  // component inductors, if I(), P() or W() was requested for any of them.
  // This is somewhat inefficent if the mutual inductor has lots of component
  // inductors, but it was the minimal change to how lead current requests are
  // tracked for all of the devices.
  setNumBranchDataVars(0);    // by default don't allocate space in branch vectors   
  numBranchDataVarsIfAllocated = numInductors;  // space allocation if power is needed
  
  // set up the device connectivity map
  // each simple inductor in this mutual inductor
  // is maked as a connection (given a common, non-zero
  // value in devConMap)
  devConMap.resize(2*numInductors);
  for(int i=0, j=0; i<(2*numInductors); i+=2, j++)
  {
    devConMap[i] = devConMap[i+1] = (j+1);
  }
  
  inductorCurrents.resize( numInductors );
  inductanceVals.resize( numInductors );
  LOI.resize( numInductors );
  LO.resize( numInductors );
  for( int i=0; i<numInductors; ++i)
  {
    LO[i].resize( numInductors );
  }
  dHe_dI.resize(numInductors);
  dManp_dI.resize(numInductors);
  ddelM_dI.resize(numInductors);
  dMirrp_dI.resize(numInductors);
  dP_dI.resize( numInductors );

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // if the user has requested output of the internal vars H and M
  // then open a file for that output.
  if( model_.outputStateVars > 0 )
  {
    outputStateVarsFlag = true;
    std::string filename( "Inductor_" );
    filename.append( getName().getEncodedName() );
    filename.append( ".dat" );
    // convert any embeded ':' or '%' characters to '_'
    replace( filename.begin(), filename.end(), '%', '_' );
    replace( filename.begin(), filename.end(), ':', '_' );

    outputFileStreamPtr = rcp( new std::ofstream() );
    outputFileStreamPtr->open( filename.c_str() );
    if( !(*outputFileStreamPtr) )
    {
      UserError(*this) << "Could not open file " << filename << " for output of state variables";
    }
    else {
      (*outputFileStreamPtr).setf(std::ios::scientific, std::ios::floatfield );
      (*outputFileStreamPtr).width(20);
      (*outputFileStreamPtr).precision(12);
    }
  }

  // update internal/external/state variable counts
  numExtVars = 2*numInductors;       // 2 nodes Vin, Vout per inductor
  numIntVars = numInductors;     // branch current per inductor
  // if we're including the deltaM and deltaHapp variables, then there will be two more internal vars
  if(includeDeltaM)
  {
    numIntVars+=1;
  }

  // allocate space for InductorOffsets
  deltaMEquInductorOffsets.resize(numInductors);

  // set up the jacobian stamp
  // for an individual inductor with the two interal variables would be
  //
  //          V1   V2   Ib
  //  kcl1               1
  //  kcl2              -1
  //  branch  1    -1   L/dt
  //
  //  for a collection of these, the internal variable, branch equations,
  //  must be at the end of a given stamp row as well as the internal
  //  vars for M and R in this non-linear version.
  //
  //  So for N inductors the samp is:
  //
  //          V1  V2  V3  V4 ... V2N  I1  I2  ... IN   M
  //  kcl1                             1
  //  kcl2                            -1
  //  kcl3                                 1
  //  kcl4                                -1
  //  branch1 1   -1                 L/dt  c  ... c
  //  branch2          1  -1          c  L/dt ... c
  //  delta M                          x    x  ... x    x
  //
  //  where "c" is an induced current change and "x" are
  //  values which must be computed.

  jacStamp.resize( numExtVars + numIntVars);

  for( int i=0; i< numInductors; ++i )
  {
    //
    // allocate space
    //
    // kcl V+ node
    jacStamp[2*i].resize(1);
    // kcl V- node
    jacStamp[2*i+1].resize(1);
    if( i == 0 )
    {
      jacStamp[2*numInductors].resize(numInductors + 2);
    }
    else
    {
      jacStamp[2*numInductors + i].resize(numInductors + 2);
    }

    //
    // fill in dependency
    //
    // kcl V+ node
    jacStamp[2*i  ][0] = 2*numInductors + i;
    // kcl V- node
    jacStamp[2*i+1][0] = 2*numInductors + i;

    if( i==0 )
    {
      jacStamp[2*numInductors ][0] = 0;
      jacStamp[2*numInductors ][1] = 1;
      for( int j=0; j<numInductors; ++j )
      {
        jacStamp[2*numInductors][j+2] = 2*numInductors + j;
      }
    }
    else
    {
      jacStamp[2*numInductors + i][0] = 2*i;
      jacStamp[2*numInductors + i][1] = 2*i + 1;
      for( int j=0; j<numInductors; ++j )
      {
        jacStamp[2*numInductors + i][j+2] = 2*numInductors + j;
      }
    }
  }

  if(includeDeltaM)
  {
    // now the deltaM equation
    jacStamp[ 3*numInductors ].resize(numInductors + 1);
    // deltaM offsets to I_1 ... I_n and deltaM
    for(int i=0; i<numInductors; ++i)
    {
      jacStamp[ 3*numInductors ][i]=2*numInductors+i;
    }
    jacStamp[ 3*numInductors ][numInductors] = 3*numInductors;

  }
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::Instance----------" << std::endl;
    Xyce::dout() << "numExtVars = " << numExtVars << ", " << ibev << std::endl
      << "numIntVars = " << numIntVars << ", " << ibiv << std::endl
      << "numStateVars = " << numStateVars << std::endl
      << "numInductors = " << numInductors << std::endl
      << "jacStamp = " << std::endl;
    for( int i = 0; i<jacStamp.size(); ++i )
    {
      Xyce::dout() << "jacStamp[ " << i << " ] = { ";
      for( int j=0; j<jacStamp[i].size(); ++j)
      {
        Xyce::dout() << jacStamp[i][j];
        if( j != ( jacStamp[i].size() -1 ) )
        {
          Xyce::dout() << ", ";
        }
      }
      Xyce::dout() << " }" << std::endl;
    }
  }
  
  // set history buffer size
  const int dMdtHistoryLength=10;
  dMdtHistory_.set_size(dMdtHistoryLength);
  for( int j=0; j<dMdtHistory_.get_size(); j++)
  {
    dMdtHistory_.push_back(0.0);
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  // Close output file if we opened one
  if( outputStateVarsFlag && outputFileStreamPtr.get() && outputFileStreamPtr->is_open() )
  {
    outputFileStreamPtr->close();
  }

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  for ( ; currentInductor != endInductor ; ++currentInductor)
  {
    if (*currentInductor != NULL)
    {
      delete *currentInductor;
      *currentInductor = NULL;
    }
  }
  instanceData.clear();
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                                          const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.
  // get the current values of the inductances and currentOffsets
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  int j = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->li_Pos = extLIDVec[ i++ ];
    (*currentInductor)->li_Neg = extLIDVec[ i++ ];
    (*currentInductor)->li_Branch = intLIDVec[ j++ ];
    currentInductor++;
  }

  if(includeDeltaM)
  {
    // now get deltaHapp and deltaM
    //li_deltaHappVar = intLIDVec[ j++ ];
    li_deltaMagVar  = intLIDVec[ j++ ];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::registerLIDs------------------------" << std::endl;
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      Xyce::dout() << "Inductor [ " << i++ << " ] "
           << "   li_Pos = " << (*currentInductor)->li_Pos
           << "   li_Neg = " << (*currentInductor)->li_Neg
           << "   li_Branch = " << (*currentInductor)->li_Branch << std::endl;
      currentInductor++;
    }
    if(includeDeltaM)
    {
      Xyce::dout() << " li_deltaMagVar = " << li_deltaMagVar << std::endl;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  //for (std::vector<InductorInstanceData *>::const_iterator it = instanceData.begin(), end = instanceData.end(); it != end; ++it)
  //  addInternalNode(symbol_table, (*it)->li_Branch, getName(), (*it)->name + "_branch");
    
  int indCount=0;
  std::string baseName = getSubcircuitName(getName());
  for (std::vector<InductorInstanceData *>::const_iterator it = instanceData.begin(), end = instanceData.end(); it != end; ++it ) {
    std::string branchInductorName = baseName;
    if( branchInductorName != "" )
      branchInductorName += Xyce::Util::separator;
    branchInductorName += (*it)->name;
    InstanceName bInductorIName = InstanceName( branchInductorName );
    std::string encodedName = spiceInternalName( bInductorIName, "branch");
    addInternalNode(symbol_table, (*it)->li_Branch, getName(), (*it)->name + "_branch");
    addInternalNode(symbol_table, (*it)->li_Branch, encodedName);
    if (loadLeadCurrent)
    {
      addBranchDataNode(symbol_table, (*it)->li_branch_data, bInductorIName, "BRANCH_D");
    }
    ++indCount;
  }
  
  addStoreNode(symbol_table, li_MagVarStore, getName(), "m");
  addStoreNode(symbol_table, li_MagVarStore, getName().getEncodedName() + "_m");
  addStoreNode(symbol_table, li_RVarStore, getName().getEncodedName() + "_r");
  addStoreNode(symbol_table, li_BVarStore, getName().getEncodedName() + "_b");
  addStoreNode(symbol_table, li_HVarStore, getName().getEncodedName() + "_h");
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;
  
  li_MagVarState=staLIDVec[0];
  li_RVarState=staLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::registerStateLIDs-------------------" << std::endl;
  }
}

//----------------------------------------------------------------------------- 
// Function      : Instance::registerBranchDataLIDs 
// Purpose       : This allows P() and W() to work for the component inductors
//               : of a mutual inductor.
// Special Notes : 
// Scope         : public 
// Creator       : Pete Sholander, SNL 
// Creation Date : 3/28/17 
//----------------------------------------------------------------------------- 
void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef) 
{   
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());   
  
  if (loadLeadCurrent)
  { 
    std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();   
    std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end(); 
    int j=0;  
    for ( ; currentInductor != endInductor ; ++currentInductor)   
    {   
      (*currentInductor)->li_branch_data = branchLIDVecRef[j];
      j++;
    }    
  } 
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 8/17/2012
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // copy over the global ID lists.
  stoLIDVec = stoLIDVecRef;

  li_MagVarStore = stoLIDVec[0];
  li_RVarStore = stoLIDVec[1];
  li_BVarStore = stoLIDVec[2];
  li_HVarStore = stoLIDVec[3];
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::registerJacLIDs ----------------------------" << std::endl;

    Xyce::dout() << "jacLIDVec = " << std::endl;
    for( int i = 0; i<jacStamp.size(); ++i )
    {
      Xyce::dout() << "jacLIDVec[ " << i << " ] = { ";
      for( int j=0; j<jacLIDVec[i].size(); ++j)
      {
        Xyce::dout() << jacLIDVec[i][j];
        if( j != ( jacLIDVec[i].size() -1 ) )
        {
          Xyce::dout() << ", ";
        }
      }
      Xyce::dout() << " }" << std::endl;
    }
  }

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->APosEquBraVarOffset  = jacLIDVec[ 2*i     ][ 0 ];
    (*currentInductor)->ANegEquBraVarOffset  = jacLIDVec[ 2*i + 1 ][ 0 ];
    (*currentInductor)->vPosOffset = jacLIDVec[ 2*numInductors + i ][ 0 ];
    (*currentInductor)->vNegOffset = jacLIDVec[ 2*numInductors + i ][ 1 ];

    (*currentInductor)->ABraEquPosNodeOffset = jacLIDVec[ 2*numInductors + i ][ 0  ];
    (*currentInductor)->ABraEquNegNodeOffset = jacLIDVec[ 2*numInductors + i ][ 1  ];
    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        (*currentInductor)->ABraEquBraVarOffset  = jacLIDVec[ 2*numInductors + i ][ j + 2  ];
      }
      (*currentInductor)->inductorCurrentOffsets[ j ] = jacLIDVec[ 2*numInductors + i ][ j + 2  ];
    }

    currentInductor++;
    i++;
  }

  if(includeDeltaM)
  {
    // now get the deltaM
    for( int i=0; i<numInductors; i++)
    {
      deltaMEquInductorOffsets[i] = jacLIDVec[ 3*numInductors ][ i ];
    }
    deltaMEquDeltaMOffset = jacLIDVec[ 3*numInductors ][ numInductors ];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      Xyce::dout() << "Inductor [ " << i << " ] " << (*currentInductor)->name << std::endl
           << "   APosEquBraVarOffset = " << (*currentInductor)->APosEquBraVarOffset << std::endl
           << "   ANegEquBraVarOffset = " << (*currentInductor)->ANegEquBraVarOffset << std::endl
           << "   vPosOffset = " << (*currentInductor)->vPosOffset << std::endl
           << "   vNegOffset = " << (*currentInductor)->vNegOffset << std::endl
           << "   ABraEquPosNodeOffset = " << (*currentInductor)->ABraEquPosNodeOffset << std::endl
           << "   ABraEquNegNodeOffset = " << (*currentInductor)->ABraEquNegNodeOffset << std::endl
           << "   ABraEquBraVarOffset = " << (*currentInductor)->ABraEquBraVarOffset << std::endl
           << "   magOffset = " << (*currentInductor)->magOffset << std::endl;
      Xyce::dout() << "\tInductor branch offsets = { ";
      for( int j=0; j<numInductors ; ++j )
      {
        Xyce::dout() << (*currentInductor)->inductorCurrentOffsets[ j ] << ", ";
      }
      Xyce::dout() << "} " << std::endl;
      i++;
      currentInductor++;
    }
    Xyce::dout() << std::endl;

    if(includeDeltaM)
    {
      Xyce::dout() << "deltaMEquInductorOffsets = ";
      for( int i=0; i<numInductors; i++ )
      {
        Xyce::dout() << deltaMEquInductorOffsets[i] << " ";
      }
      Xyce::dout() //<< "deltaMEquDeltaHappOffset = " << deltaMEquDeltaHappOffset
        << " deltaMEquDeltaMOffset = " << deltaMEquDeltaMOffset << std::endl;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // Because we have saved the inductances from parameters in a local storage,
  //  we need to re-read them just in case the L values were specified as
  // dependent expressions (e.g. dependent on global params that are stepped)
  std::vector< InductorInstanceData* >::iterator
    currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator
    endInductor = instanceData.end();

  int i=0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->L = inductorInductances[i];
    (*currentInductor)->baseL = inductorInductances[i];
    ++i;
    ++currentInductor;
  }

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;

  // current temp difference from reference temp.
  double difference = temp - model_.tnom;

  std::vector< InductorInstanceData* >::iterator currentData = instanceData.begin();
  while( currentData != instanceData.end() )
  {
    double factor = 1.0 + (model_.tempCoeff1)*difference +
                          (model_.tempCoeff2)*difference*difference;
    (*currentData)->L = ((*currentData)->baseL)*factor;
    currentData++;
  }

  // now that the inductances have changed we need to update the matrix.
  updateInductanceMatrix();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : updates a set of common variables used by RHS and jacobian
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS)  && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updateIntermediateVars" << std::endl;
  }
  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);

  // some parameters in the model class that we will use often
  const double A      = model_.A;
  const double Alpha  = model_.Alpha;
  const double BetaH  = model_.BetaH;
  const double C      = model_.C;
  const double DeltaV = model_.DeltaV;
  const double Gap    = model_.Gap;
  const double Ms     = model_.Ms;
  const double Kirr   = model_.Kirr;
  const double Path   = model_.Path;
  const double Vinf   = model_.Vinf;

  double latestMag; //NoMag  = solVector[ li_MagVar ];

  //sum of currents through the inductors
  branchCurrentSum = 0.0;
  for(int i=0; i<numInductors; i++)
  {
    branchCurrentSum += solVector[instanceData[i]->li_Branch] * inductanceVals[ i ];
  }

  latestMag = MagVar  + MagVarUpdate;

  // used in voltage drop over first inductor
  double V1Pos = solVector[(instanceData[0])->li_Pos];
  double V1Neg = solVector[(instanceData[0])->li_Neg];

  double qV = (DeltaV / Vinf) * (V1Pos - V1Neg);

  double tanh_qV = 0.0;
  if (fabs(qV) < CONSTTANH_THRESH)
  {
    tanh_qV = tanh(qV);
  }
  else if (qV < 0.0)
  {
    tanh_qV = -1.0;
  }
  else
  {
    tanh_qV = 1.0;
  }

  Happ = branchCurrentSum / Path;
  double dHdt=0.0;
  double H = Happ - (Gap / Path) * latestMag;
  // not actually dHdt but we only need to test if this > 0 or < 0
  dHdt = - (Gap / Path) * MagVarUpdate;
  double He = H + Alpha * latestMag;
  /*
  if( dHdt > 0.0) 
  {
    tanh_qV = 0.0;
  }
  else if (dHdt < 0.0)
  {
    tanh_qV = 0.0;
  }
  else 
  {
    tanh_qV = 0.0;
  }
  */
  
  double Heo = BetaH*A;

  // terms that come up frequently
  double gap_path = Gap / Path;
  double He2 = He*He;
  double Heo2 = Heo*Heo;
  double sq_Heo2He2 = sqrt(Heo2 + He2);

  double delM0 = model_.BetaM * Ms;
  double Man = Ms * He / ( A + sq_Heo2He2 );
  double delM = Man - latestMag;

  double delM2 = delM*delM;
  double delM02 = delM0*delM0;
  double sq_delM02delM2 = sqrt( delM02 + delM2 );

  double delta0=1.0;
  double deltaM=1.0;

  double Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
  double Manp =  Ms*(A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
  /*
  if( ((dHdt < 0.0) && ((Manp - latestMag)>0.0)) || ((dHdt >= 0.0) && ((Manp - latestMag)<0.0) ) )
  {
    deltaM=1.0;
    Mirrp = 0.0;
  }
  */
  P = ( C * Manp + (1 - C) * Mirrp)        / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);
  //P = ( C * deltaM * (Manp-Mirrp) + Mirrp) / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);


  // at this point we have P so now we can update mag.
  //
  //  The problem is that if deltaM is too big, then we need to shrink
  //  the time step.  One way to control this is to set a max time
  //  step.  But what we really need to do is calculate deltaM and
  //  then if it's over some fraction of Ms then turn on the limiting
  //  flag (or bail on the step but I think turning on limiting is
  //  safer and if we hit maxItter with it on then we'll get that step
  //  rejected.

  if( useRKIntegration )
  {
    // use 4th order runga-kutta to estimate MagVarUpdate
    double stepLen = branchCurrentSumHistory[0] + branchCurrentSumHistory[1] + branchCurrentSumHistory[2] + (branchCurrentSum - oldBranchCurrentSum);
    MagVarUpdate = stepLen * (  PFunctionHistory[0] +
                    2*PFunctionHistory[1] +
                    2*PFunctionHistory[2] +
                      P) / 6;
  }
  else
  {
    // forward euler method
    MagVarUpdate = P * (branchCurrentSum - oldBranchCurrentSum) / model_.Path;

    // trap
    //    double MagVarUpdateWithTrap = 0.5 * (P + PPreviousStep) * (branchCurrentSum - oldBranchCurrentSum) / model_.Path;
    origFlag = true;
    if( fabs( MagVarUpdate ) > 0.25 * Ms )
    {
      // step was too big, so
      // turn on limiting
      origFlag = false;
    }
  }

  latestMag = MagVar + MagVarUpdate;

  if(includeDeltaM)
  {
    // in this case we're scaling MagVarUpdate by Ms because it's being solved
    // with the full system and big changes in
    MagVarUpdate /=model_.Ms;
  }

  double dP_Denom = 1.0 + (gap_path - Alpha)*C*Manp + gap_path * (1.0-C) * Mirrp;

  // now get dP_dI for each inductor
  for( int i=0; i<numInductors; i++)
  {

    dHe_dI[ i ] = inductanceVals[ i ] / Path;

    dManp_dI[i] = ( -Ms * He / (pow(A + sq_Heo2He2, 2.0)*sq_Heo2He2)) *
                   ( (Heo2 / (Heo2 + He2)) + (2.0*(A + Heo2 / sq_Heo2He2)/(A+sq_Heo2He2)) ) * dHe_dI[i];

    ddelM_dI[i] = (Ms / (A + sq_Heo2He2)) * (1.0 - He2/((A + sq_Heo2He2)*sq_Heo2He2)) * dHe_dI[i];
 
    dMirrp_dI[i] = (1.0/(2.0*(Kirr - Alpha*sq_delM02delM2))) *
                   (tanh_qV + delM/sq_delM02delM2 +
                     (2.0*Alpha*delM*(delM*tanh_qV +
                       sq_delM02delM2)/(2.0*(Kirr-Alpha*sq_delM02delM2)*sq_delM02delM2))) * ddelM_dI[i];

    dP_dI[i] = (1.0/dP_Denom) * (C * dManp_dI[i] + (1.0-C) * dMirrp_dI[i]) -
          ( (C*deltaM*(Manp-Mirrp) + Mirrp)/pow(dP_Denom,2.0) ) *
            ( (gap_path - Alpha)*C*dManp_dI[i] + gap_path*(1.0-C)*dMirrp_dI[i] );

  }
  
  Linear::Vector & stoVector = *(extData.nextStoVectorPtr);
  stoVector[li_MagVarStore] = latestMag; 
  stoVector[ li_HVarStore ] = model_.HCgsFactor * (Happ  - (model_.Gap / model_.Path) * latestMag);
  stoVector[ li_BVarStore ] = model_.BCgsFactor * (4.0e-7 * M_PI * (stoVector[ li_HVarStore ] + latestMag));
  //stoVector[ li_BVarStore ] = model_.BCgsFactor * (4.0e-7 * M_PI * (stoVector[ li_HVarStore ] ));
  /*
  if( (stoVector[ li_BVarStore ] > 0.0) && (stoVector[ li_HVarStore ] > lastH) )
    stoVector[ li_HVarStore ] = lastH;
  if( (stoVector[ li_BVarStore ] < 0.0) && (stoVector[ li_HVarStore ] > lastH) )
    stoVector[ li_HVarStore ] = lastH;
  lastH=stoVector[ li_HVarStore ];
  */
  // when the time derivative is calculated of this we will have R
  // R = dHapp/dt;
  
  staVector[li_RVarState] = Happ;
  staVector[li_MagVarState] = latestMag;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateInductanceMatrix()
// Purpose       : A matrix of inductances is used often enough that it
//                 calculated and stored as a member variable here
//                 If and inductance ever changes say from updating
//                 the temperature or a parameter udpate, then this
//                 routine must be called again.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::updateInductanceMatrix()
{
  std::vector< InductorInstanceData* >::iterator
    currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator
    endInductor = instanceData.end();

  // collec the inductances
  int i=0;
  while( currentInductor != endInductor )
  {
    inductanceVals[ i ] = (*currentInductor)->L;
    i++;
    currentInductor++;
  }

  double Area = model_.Area;
  double Path = model_.Path;

  // compute the inductance matrix
  for( i=0; i<numInductors; ++i)
  {
    for( int j=0; j<numInductors; ++j)
    {
      // 4.0e-7 * M_PI is a magnetic constant, the permeability of free space [Henries/m]
      LO[i][j] = mutualCup * 4.0e-7 * M_PI * (Area / Path) * inductanceVals[i] * inductanceVals[j];
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState---------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
  // udate dependent parameters
  updateIntermediateVars ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updateSecondaryState-------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector & staVector = *(extData.nextStaVectorPtr);
  Linear::Vector & staDerivVec = *(extData.nextStaDerivVectorPtr);
    
  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & stoVector = *(extData.nextStoVectorPtr);
  Linear::Vector & stoVectorCurr = *(extData.currStoVectorPtr);
  
  //double mVarScaling = model_.mVarScaling;
  // place current values of mag, H and R in state vector
  // must unscale them as the rest of the class assumes
  // that these aren't scaled yet.
  double latestMag = stoVector[li_MagVarStore];
  /*
  if( model_.includeMEquation ) 
  {
    staVector[ li_MagVarState ] = solVector[ li_MagVar ];
    latestMag = mVarScaling * solVector[ li_MagVar ];
  }
  else
  {
    latestMag = mVarScaling * staVector[ li_MagVarState ];
  }
  if( model_.factorMS )
  {
    latestMag *= model_.Ms;
  }
  */
  // place a copy of R in the store vector for lookup for output.
  stoVector[ li_RVarStore ] = staDerivVec[ li_RVarState ];
  
  //
  // need these to calculate H for B-H loops and be careful of
  // potential non-physical turning points in the B-H phase 
  // in general dB/dH should be >= 0.  If it's less than
  // zero then we have a non-physical turning point that 
  // is ok as part of the device model, but not physically
  // realistic as the change incurrent (dH) is apposing
  // the magnetic field (dB) will be an irreversible loss.
  //
  // So if dB/dH is negative, hold H constant while B 
  // changes to get the correct path along the B-H curve.
  // Unless the gap is non-zero then things are slightly 
  // more complex and we need to look at the dM/dt and dH/dt terms.
  //
  // dM/dt is equivalent to dB/dt within a scalar factor
  // dH/dt = dHapp/dt - (gap/path) dM/dt = R - (gap/path) dM/dt
  
  double lastB = stoVectorCurr[ li_BVarStore ];
  double lastH = stoVectorCurr[ li_HVarStore ];
  double calculatedH = model_.HCgsFactor * (Happ  - (model_.Gap / model_.Path) * latestMag);
  double dMdt = staDerivVec[ li_MagVarState ];
  double R = staDerivVec[ li_RVarState ];
  double dHdt = R - (model_.Gap/model_.Path)*dMdt;
  
  double Hfxn = Happ;
  
  if( model_.Gap <= 0 )
  {    
    if( ((dMdtAverage_<0.0) && (dHdt>0.0)) || ((dMdtAverage_>0.0) && (dHdt<0.0)) )
    {
      // derivatives disagree so changes in H are lost.
      Hfxn = lastH / model_.HCgsFactor;
    }
    else if( ((dMdtAverage_<0.0) && (lastH < calculatedH)) || ((dMdtAverage_>0.0) && (lastH > calculatedH)) )
    {
      // derivatives agree but H calculated may be in wrong direction 
      Hfxn = lastH / model_.HCgsFactor;
    }
  }
  else
  {
    double gapFactor =   - (model_.Gap / model_.Path) * latestMag;
    if( (fabs( gapFactor) < fabs( Hfxn )) && ( ((gapFactor < 0) && (Hfxn < 0)) || ((gapFactor > 0) && (Hfxn > 0)) ))
    {
      Hfxn += gapFactor;
    }
  }
  stoVector[ li_HVarStore ] = model_.HCgsFactor * Hfxn;
  stoVector[ li_BVarStore ] = model_.BCgsFactor * (4.0e-7 * M_PI * (stoVector[ li_HVarStore ] + latestMag));
  
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       : This function updates the value of MagVar
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 01/25/2011
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  if (!getSolverState().dcopFlag)
  {
    if(includeDeltaM)
    {
      MagVar += MagVarUpdate*model_.Ms;
    }
    else
    {
      MagVar += MagVarUpdate;
    }
    lastMagUpdate = MagVarUpdate;
    PPreviousStep = P;
    if( fabs(MagVar) > 2*model_.Ms )
    {
      MagVar = 0.0;
    }

    if( useRKIntegration )
    {
      // fill in history for RK integration of dM/dH
      for(int i=0; i<2; i++)
      {
        branchCurrentSumHistory[i] = branchCurrentSumHistory[i+1];
        PFunctionHistory[i] = PFunctionHistory[i+1];
      }
      branchCurrentSumHistory[2] = branchCurrentSum-oldBranchCurrentSum;
      PFunctionHistory[2] = PPreviousStep;
    }
    oldBranchCurrentSum = branchCurrentSum;
    Linear::Vector & staDerivVec = *(extData.nextStaDerivVectorPtr);
    dMdtHistory_.push_back( staDerivVec[ li_MagVarState ] );

    Linear::Vector & stoVector = *(extData.nextStoVectorPtr);

    // update average magnetic gradient.
    for( int j=0; j< dMdtHistory_.get_size(); j++)
    {
      dMdtAverage_ += dMdtHistory_.at_from_head(j);
    }
    dMdtAverage_ /= (dMdtHistory_.get_size());
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEQVector------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;
  Linear::Vector & solVector = *(extData.nextSolVectorPtr);

  // update LOI -- the following product
  // I = column vector of currents
  // L = row vector of inductances
  // LO = matrix = mutualCup * sqrt( L' * L )
  // LOI = column vector = mutualCup * sqrt( L' * L ) * I
  // LOI[1] = mutualCup * sqrt(L[1]*L[1])*I[1]) +
  //          mutualCup * sqrt(L[1]*L[2])*I[2]) + ...
  //          mutualCup * sqrt(L[1]*L[n])*I[n])

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    if( (getSolverState().dcopFlag) && (*currentInductor)->ICGiven == true )
    {
      inductorCurrents[ i ] = (*currentInductor)->IC;
    }
    else
    {
      inductorCurrents[ i ] = solVector[ (*currentInductor)->li_Branch ];
    }
    i++;
    currentInductor++;
  }

  for( i = 0; i < numInductors; ++i )
  {
    LOI[ i ] = 0;
    for( int j = 0; j < numInductors; ++j )
    {
      LOI[i] += LO[i][j] * inductorCurrents[j];
    }
  }

  // loop over each inductor and load it's Q vector components
  // and each inductor's contribution to the R equ.
  currentInductor = instanceData.begin();
  endInductor = instanceData.end();
  i = 0;
  while( currentInductor != endInductor )
  {

    (*daeQVecPtr)[((*currentInductor)->li_Branch)] += LOI[ i ];
    i++;
    currentInductor++;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 Same as loadRHS, but without the capacitor
//                 currents.
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEFVector------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;
  Linear::Vector & solVector = *(extData.nextSolVectorPtr);

  double Gap = model_.Gap;
  double Path = model_.Path;


  // used in scaling the branch equation;
  double mid = 1.0 + (1.0 - (Gap / Path))*P;

  // loop over each inductor and load it's F vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i=0;
  while( currentInductor != endInductor )
  {
    double current   = solVector[(*currentInductor)->li_Branch];
    double vNodePos  = solVector[(*currentInductor)->li_Pos];
    double vNodeNeg  = solVector[(*currentInductor)->li_Neg];


    (*daeFVecPtr)[((*currentInductor)->li_Pos)]    +=  current;

    (*daeFVecPtr)[((*currentInductor)->li_Neg)]    += -current;

    (*daeFVecPtr)[((*currentInductor)->li_Branch)] += -((vNodePos - vNodeNeg)/mid);
    if (loadLeadCurrent)
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;     
      double * junctionV = extData.nextJunctionVCompRawPtr;
      leadF[(*currentInductor)->li_branch_data] =  current;       
      junctionV[(*currentInductor)->li_branch_data] = (vNodePos - vNodeNeg);
    }

    currentInductor++;
    i++;
  }

  if(includeDeltaM)
  {
    // the deltaHapp equation
    //(*daeFVecPtr)[li_deltaHappVar] += solVector[li_deltaHappVar];
    //(*daeFVecPtr)[li_deltaHappVar] -= HappVarUpdate;

    // the deltaM equation
    (*daeFVecPtr)[li_deltaMagVar] += solVector[li_deltaMagVar];
    (*daeFVecPtr)[li_deltaMagVar] -= MagVarUpdate;
  }
  
   
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;
  int i;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdQdx-----------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  // loop over each inductor and load it's Q vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  i = 0;
  while( currentInductor != endInductor )
  {
    for( int j=0; j<numInductors; ++j )
    {
      (*dQdxMatPtr)[((*currentInductor)->li_Branch)]
                   [(*currentInductor)->inductorCurrentOffsets[j]] += LO[i][j];
    }
    i++;
    currentInductor++;
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdFdx----------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & lastSolVector = *(extData.lastSolVectorPtr);
  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  // pull these parameters up from the model class to make it easier
  // to view the equations.
  const double Gap = model_.Gap;
  const double Path = model_.Path;


  // loop over each inductor and load it's dFdx components
  double mid = 1.0 + (1.0 - (Gap/Path))*P;

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    // do the normal work for an inductor
    (*dFdxMatPtr)[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  +=  1.0;
    (*dFdxMatPtr)[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += -1.0;
    (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += -1.0/mid;
    (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] +=  1.0/mid;

    double delV = solVector[(*currentInductor)->li_Pos] - solVector[(*currentInductor)->li_Neg];

    for( int j = 0; j<numInductors; ++j )
    {

      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] +=
        delV * (1.0 - (Gap/Path)) * dP_dI[j]/(mid*mid);

    }

    currentInductor++;
  }

  if(includeDeltaM)
  {
    (*dFdxMatPtr)[li_deltaMagVar][deltaMEquDeltaMOffset] = 1.0;

    for( int i=0; i<numInductors; i++ )
    {
      (*dFdxMatPtr)[li_deltaMagVar][deltaMEquInductorOffsets[i]] =
        -((inductanceVals[i]*(solVector[ instanceData[i]->li_Branch ] - lastSolVector[ instanceData[i]->li_Branch ] ) * dP_dI[i])
        + P*inductanceVals[i])/(model_.Path*model_.Ms );
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       : If requested by the use in the model statement,
//                 this routine outputs values of the internal
//                 state variables M, H and R to a file
//                 named "Inductor_name.dat".  File is opened
//                 and closed in the contructor and destructor.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles(bool force_final_output)
{
  bool bsuccess = true;
  if( outputStateVarsFlag && outputFileStreamPtr.get() && (*outputFileStreamPtr) )
  {

    double latestMag = MagVar;

    if( includeDeltaM )
    {
      (*outputFileStreamPtr)
        << getSolverState().currTime_ << "  "
        << latestMag
        << std::endl;
    }
    else
    {
      (*outputFileStreamPtr)
        << getSolverState().currTime_ << "  "
        << latestMag
        << std::endl;
    }

  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(numInductors);
  for(int i=0; i<numInductors; i++)
  {
    varTypeVec[i] = 'I';
  }
}



//-----------------------------------------------------------------------------
// Function      : Instance::Pcalc
// Purpose       :
  // this is a templated function for a complicated term P(M,I_1... I_n) that relates
  // the magnetic saturation of the mutual indcutor to the individual currents
  // through the inductors.  We'll need dP_dM and this tempated function automates
  // that calculation via Sacado
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 10/13/2011
//-----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT Instance::Pcalc(
    const ScalarT & Mag, const ScalarT & CurrentSum,
    const ScalarT & Vpos, const ScalarT & Vneg)
   // independent variable M
   // independent variable Sum(n_i * I_i) windings * current through each inductor
   // independent variable Vpos, Vneg -- voltage drop over first inductor
{
     // some parameters in the model class that we will use often
     const double A      = model_.A;
     const double Alpha  = model_.Alpha;
     const double BetaH  = model_.BetaH;
     const double C      = model_.C;
     const double DeltaV = model_.DeltaV;
     const double Gap    = model_.Gap;
     const double Ms     = model_.Ms;
     const double Kirr   = model_.Kirr;
     const double Path   = model_.Path;
     const double Vinf   = model_.Vinf;

     ScalarT qV = (DeltaV / Vinf) * (Vpos - Vneg);

     ScalarT tanh_qV = 0.0;
     if (fabs(qV) < CONSTTANH_THRESH)
     {
       ScalarT tanh_qV = tanh(qV);
     }
     else if (qV < 0.0)
     {
       tanh_qV = -1.0;
     }
     else
     {
       tanh_qV = 1.0;
     }

     ScalarT Happ = CurrentSum / Path;

     ScalarT H = Happ - (Gap / Path) * Mag;
     ScalarT He = H + Alpha * Mag;

     ScalarT Heo = BetaH*A;

     // terms that come up frequently
     ScalarT gap_path = Gap / Path;
     ScalarT He2 = He*He;
     ScalarT Heo2 = Heo*Heo;
     ScalarT sq_Heo2He2 = sqrt(Heo2 + He2);

     ScalarT delM0 = model_.BetaM * Ms;
     ScalarT Man = Ms * He / ( A + sq_Heo2He2 );
     ScalarT delM = Man - Mag;

     ScalarT delM2 = delM*delM;
     ScalarT delM02 = delM0*delM0;
     ScalarT sq_delM02delM2 = sqrt( delM02 + delM2 );

     ScalarT Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
     ScalarT Manp =  Ms*(A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
     ScalarT Pval = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);

     return Pval;
}  // end of Pcalc() function


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  // scale gap, path and area from cm and cm^2 to m and m^2
  Gap = 1.0e-2 * GapInCm;
  Path = 1.0e-2 * PathInCm;
  Area = 1.0e-4 * AreaInCm2;
  
  if( BHSiUnits != 0 )
  {
    // user requested SI units over the default of CGS units.  Change
    // conversion factor to unity.
    BCgsFactor=1.0;
    HCgsFactor=1.0;
  }
  
  // Set any non-constant parameter defaults:

  if (!given("TNOM"))
  {
    tnom = getDeviceOptions().tnom;
  }


  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    A(0.0),
    Alpha(0.0),
    AreaInCm2(0.0),
    Area(0.0),
    BetaH(0.0),
    BetaM(0.0),
    C(0.0),
    DeltaV(0.0),
    GapInCm(0.0),
    Gap(0.0),
    Kirr(0.0),
    Ms(0.0),
    PathInCm(0.0),
    Path(0.0),
    Vinf(0.0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tnom(getDeviceOptions().tnom),
    HCgsFactor( 0.012566370614359 ),  // 4 pi / 1000
    BCgsFactor( 10000.0 ),
    outputStateVars(0),
    factorMS(0),
    BHSiUnits(0),
    includeDeltaM(0),
    useRKIntegration(0)
{
  setLevel(2);


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

// additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of MutIndNonLin instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("min", 2)
    .registerModelType("min", 2)
    .registerModelType("core", 2);
}

} // namespace MutIndNonLin2
} // namespace Device
} // namespace Xyce

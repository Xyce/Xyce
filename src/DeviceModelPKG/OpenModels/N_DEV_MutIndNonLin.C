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
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_MutIndNonLin.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_FeatureTest.h>
#include <N_UTL_Math.h>
#include <N_UTL_HspiceBools.h>

//This contains important constants like permitivity of free space
#include <N_DEV_Const.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

using Teuchos::rcp;

namespace Xyce {
namespace Device {

namespace MutIndNonLin {

void Traits::loadInstanceParameters(ParametricData<MutIndNonLin::Instance> &p)
{
  p.addPar("COUP_VAL",1.0,&MutIndNonLin::Instance::mutualCup)
  .setGivenMember(&MutIndNonLin::Instance::mutualCupGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Coupling coefficient");

  p.addPar("NONLINEARCOUPLING",0.0,&MutIndNonLin::Instance::nonlinFlag)
  .setGivenMember(&MutIndNonLin::Instance::nonlinFlagGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Nonlinear coupling flag");

  p.addPar("COUPLEDMutIndNonLin",std::vector<std::string>(),&MutIndNonLin::Instance::inductorNames)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("");

  p.addPar("COUPLEDINDUCTANCE",std::vector<double>(),&MutIndNonLin::Instance::inductorInductances)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("");

  p.addPar("NODE1",std::vector<std::string>(),&MutIndNonLin::Instance::inductorsNode1)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("");

  p.addPar("NODE2",std::vector<std::string>(),&MutIndNonLin::Instance::inductorsNode2)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("");

  p.addPar("COUPLING",std::vector<double>(),&MutIndNonLin::Instance::couplingCoefficient)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Coupling coefficient");

  p.addPar("COUPLEDINDUCTOR",std::vector<std::string>(),&MutIndNonLin::Instance::couplingInductor)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("");
   
  p.addPar ("IC",std::vector<double>(),&MutIndNonLin::Instance::initialCondition)
   .setUnit(U_AMP)
   .setCategory(CAT_NONE)
   .setDescription("Initial current through the inductor.");
}

void Traits::loadModelParameters(ParametricData<MutIndNonLin::Model> &p)
{
  p.addPar("A",1000.0,&MutIndNonLin::Model::A)
  .setUnit(U_AMPMM1)
  .setCategory(CAT_MATERIAL)
  .setDescription("Thermal energy parameter");

  p.addPar("AREA",0.1,&MutIndNonLin::Model::AreaInCm2)
  .setUnit(U_CM2)
  .setCategory(CAT_GEOMETRY)
  .setDescription("Mean magnetic cross-sectional area");

  p.addPar("ALPHA",5.0e-5,&MutIndNonLin::Model::Alpha)
  .setUnit(U_NONE)
  .setCategory(CAT_GEOMETRY)
  .setDescription("Domain coupling parameter");

  p.addPar("BETAH",0.0001,&MutIndNonLin::Model::BetaH)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Modeling constant");

  p.addPar("BETAM",3.125e-5,&MutIndNonLin::Model::BetaM)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Modeling constant");

  p.addPar("C",0.2,&MutIndNonLin::Model::C)
  .setUnit(U_NONE)
  .setCategory(CAT_MATERIAL)
  .setDescription("Domain flexing parameter");

  p.addPar("CLIM",0.005,&MutIndNonLin::Model::CLim)
  .setUnit(U_NONE)
  .setCategory(CAT_MATERIAL)
  .setDescription("Value below which domain flexing parameter will be treated as zero.");

  p.addPar("DELVSCALING",1.0e3,&MutIndNonLin::Model::DeltaVScaling)
  .setUnit(U_VOLT)
  .setCategory(CAT_NONE)
  .setDescription("Smoothing coefficient for voltage difference over first inductor");

  p.addPar("CONSTDELVSCALING",false,&MutIndNonLin::Model::UseConstantDeltaVScaling)
  .setUnit(U_VOLT)
  .setCategory(CAT_NONE)
  .setDescription("Use constant scaling factor to smooth voltage difference over first inductor");

  p.addPar("INCLUDEMEQU",true,&MutIndNonLin::Model::includeMEquation)
  .setGivenMember(&MutIndNonLin::Model::includeMEquationGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Flag to include the magnetics in the solution.");

  p.addPar("GAP",0.0,&MutIndNonLin::Model::GapInCm)
  .setUnit(U_CM)
  .setCategory(CAT_GEOMETRY)
  .setDescription("Effective air gap");

  p.addPar("K",500.0,&MutIndNonLin::Model::Kirr)
  .setUnit(U_AMPMM1)
  .setCategory(CAT_MATERIAL)
  .setDescription("Domain anisotropy parameter");

  p.addPar("KIRR",500.0,&MutIndNonLin::Model::Kirr)
  .setUnit(U_AMPMM1)
  .setCategory(CAT_MATERIAL)
  .setDescription("Domain anisotropy parameter");

  p.addPar("MS",1.0e+6,&MutIndNonLin::Model::Ms)
  .setUnit(U_AMPMM1)
  .setCategory(CAT_MATERIAL)
  .setDescription("Saturation magnetization");

  p.addPar("LEVEL",0.0,&MutIndNonLin::Model::LevelIgnored)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("for pspice compatibility -- ignored");

  p.addPar("PACK",0.0,&MutIndNonLin::Model::PackIgnored)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("for pspice compatibility -- ignored");

  p.addPar("PATH",1.0,&MutIndNonLin::Model::PathInCm)
  .setUnit(U_CM)
  .setCategory(CAT_GEOMETRY)
  .setDescription("Total mean magnetic path");

  p.addPar("TNOM",27.0,&MutIndNonLin::Model::tnom)
  .setUnit(U_DEGC)
  .setCategory(CAT_MATERIAL)
  .setDescription("Reference temperature");

  p.addPar("TC1",0.0,&MutIndNonLin::Model::tempCoeff1)
  .setUnit(U_NONE)
  .setCategory(CAT_MATERIAL)
  .setDescription("First order temperature coeff.");

  p.addPar("TC2",0.0,&MutIndNonLin::Model::tempCoeff2)
  .setUnit(U_NONE)
  .setCategory(CAT_MATERIAL)
  .setDescription("Second order temperature coeff.");

  p.addPar("PZEROTOL",0.1,&MutIndNonLin::Model::pZeroTol)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Tolerance for nonlinear zero crossing");

  p.addPar("MVARSCALING",1.0,&MutIndNonLin::Model::mVarScaling)
  .setGivenMember(&MutIndNonLin::Model::mVarScalingGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("M-variable scaling.");

  p.addPar("RVARSCALING",1.0,&MutIndNonLin::Model::rVarScaling)
  .setGivenMember(&MutIndNonLin::Model::rVarScalingGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("R-variable scaling");

  p.addPar("MEQNSCALING",1.0,&MutIndNonLin::Model::mEqScaling)
  .setGivenMember(&MutIndNonLin::Model::mEqScalingGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("M-equation scaling");

  p.addPar("REQNSCALING",1.0,&MutIndNonLin::Model::rEqScaling)
  .setGivenMember(&MutIndNonLin::Model::rEqScalingGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("R-equation scaling");

  p.addPar("OUTPUTSTATEVARS",0,&MutIndNonLin::Model::outputStateVars)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Flag to save state variables");

  p.addPar("FACTORMS",0,&MutIndNonLin::Model::factorMS)
  .setGivenMember(&MutIndNonLin::Model::factorMSGiven)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Flag to factor the saturation magnetization from the magnetics equation.");

  p.addPar("BHSIUNITS",0,&MutIndNonLin::Model::BHSiUnits)
  .setUnit(U_NONE)
  .setCategory(CAT_NONE)
  .setDescription("Flag to report B and H in SI units");
}


// Class Instance

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
    maxVoltageDrop(1.0e-10),
    outputStateVarsFlag( false )
{
  scalingRHS = 1.0;
  numExtVars   = 2;
  numIntVars   = 3;

  numStateVars = 2;
  setNumStoreVars(3);

  tempGiven    = false;

  const int ibev = IB.numExtVars;
  const int ibiv = IB.numIntVars;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);
  
  // look over IB params for IC data
  std::vector<Param>::const_iterator paramIt = IB.params.begin();
  for( ;paramIt != IB.params.end(); ++paramIt)
  {
    if( (paramIt->tag() == "IC") && (paramIt->getType() == Xyce::Util::STR))
    {
      // in the process of packing up the component inductors into a mutual inductor
      // whether an initial condition is given or not is lost.  So check if the
      // initial condition is nonzero and assue that zero was not given 
      initialCondition.push_back(paramIt->getImmutableValue<double>());
      if( paramIt->getImmutableValue<double>() != 0)
      {
        initialConditionGiven.push_back(true);
      }
      else
      {
        initialConditionGiven.push_back(false);
      }
    }
  }
  // now load the instance data vector
  for( int i=0; i<inductorNames.size(); ++i )
  {
    InductorInstanceData * inductorData = new InductorInstanceData();
    inductorData->name = inductorNames[i];
    inductorData->L = inductorInductances[i];
    inductorData->baseL = inductorInductances[i];
    // if this is true then the instance block had some IC data, so don't ignore it.
    if( i < initialCondition.size())
    {
      inductorData->ICGiven = initialConditionGiven[i];
      inductorData->IC=initialCondition[i];
    }
    else
    {
      inductorData->ICGiven = false;
      inductorData->IC = 0.0;
    }
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

  mEquInductorOffsets.resize( numInductors );
  rEquInductorOffsets.resize( numInductors );
  inductorCurrents.resize( numInductors );
  inductanceVals.resize( numInductors );
  LOI.resize( numInductors );
  LO.resize( numInductors );
  for( int i=0; i<numInductors; ++i)
  {
    LO[i].resize( numInductors );
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // if the user has requested output of the state variables M, H and R
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

  // size some vectors needed in loadDAEdFdx
  dHe_dI.resize( numInductors );
  dManp_dI.resize( numInductors );
  ddelM_dI.resize( numInductors );
  dMirrp_dI.resize( numInductors );
  dP_dI.resize( numInductors );

  // Check if the magnetic moment equation should be dropped form the 
  // system of equations because domain flexing is essentially zero 
  if( model_.C <= model_.CLim )
  {
    model_.includeMEquation=false;
  }

  // update internal/external/state variable counts
  numExtVars = 2*numInductors;
  int numExtraEquations = 1;  // Allways use the R equation
  if( model_.includeMEquation )
  {
    numExtraEquations++;      // Also add the M equation
  }
  numIntVars = numInductors + numExtraEquations;
  numStateVars = 5;  // extra state variables for M and dM/dt and R, H & B
  //numStateVars += 2*numInductors;  individual inductors no longer need state / store space

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
  //          V1  V2  V3  V4 ... V2N  I1  I2  ... IN  M  R
  //  kcl1                             1
  //  kcl2                            -1
  //  kcl3                                 1
  //  kcl4                                -1
  //  branch1 1   -1                 L/dt  c  ... c   x
  //  branch2 x    x   1  -1          c  L/dt ... c   x
  //  M equ   x    x                  x    x  ... x   x  x
  //  R equ                           x    x  ... x      x
  //
  //  where "c" is an induced current change and "x" are
  //  values which must be computed.

  jacStamp.resize( 2 * numInductors + numIntVars);

  for( int i=0; i< numInductors; ++i )
  {
    //
    // allocate space
    //
    // kcl V+ node
    jacStamp[2*i].resize(1);
    // kcl V- node
    jacStamp[2*i+1].resize(1);
    // branch node -- every branch needs to access the first
    // inductor's V+ and V-, so they all contribute there
    if( i == 0 )
    {
      if( model_.includeMEquation )
      { 
        jacStamp[2*numInductors].resize(numInductors + 3);
      }
      else
      {
        jacStamp[2*numInductors].resize(numInductors + 2);
      }
    }
    else
    {
      if( model_.includeMEquation )
      { 
        jacStamp[2*numInductors + i].resize(numInductors + 5);
      }
      else
      {
        jacStamp[2*numInductors + i].resize(numInductors + 4);
      }
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
      if( model_.includeMEquation )
      { 
        jacStamp[2*numInductors][numInductors+2] = 3*numInductors;
      }
    }
    else
    {
      jacStamp[2*numInductors + i][0] = 0;
      jacStamp[2*numInductors + i][1] = 1;
      jacStamp[2*numInductors + i][2] = 2*i;
      jacStamp[2*numInductors + i][3] = 2*i + 1;
      for( int j=0; j<numInductors; ++j )
      {
        jacStamp[2*numInductors + i][j+4] = 2*numInductors + j;
      }
      // only add this term if M equation is in use
      if( model_.includeMEquation )
      { 
        jacStamp[2*numInductors + i][numInductors+4] = 3*numInductors;
      }
    }
  }
 
  int rOffset=0; 
  if( model_.includeMEquation )
  {
    // now the M equation
    jacStamp[ 3*numInductors    ].resize(numInductors + 4);
    // M offsets for V+ and V- on first inductor
    jacStamp[ 3*numInductors    ][0] = 0;
    jacStamp[ 3*numInductors    ][1] = 1;
    // M offsets to each inductor's branch equ.
    for(int i=0; i<numInductors; ++i)
    {
      jacStamp[ 3*numInductors     ][i+2]=2*numInductors+i;
    }
    // M offsets to M and R
    jacStamp[ 3*numInductors    ][numInductors + 2] = 3*numInductors;
    jacStamp[ 3*numInductors    ][numInductors + 3] = 3*numInductors + 1;
    rOffset=1;
  }
 
  // and the R equations
  jacStamp[ 3*numInductors + rOffset].resize(numInductors + 1);

  // R offsets to each inductor's branch equ.
  for(int i=0; i<numInductors; ++i)
  {
    jacStamp[ 3*numInductors + rOffset ][i]=2*numInductors+i;
  }

  // R offsets to R
  jacStamp[ 3*numInductors + rOffset][numInductors ]    = 3*numInductors + rOffset;

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
  
  if( model_.UseConstantDeltaVScaling )
  {
    // set scaling to 1.0 so it can be safely factored in as a constant.
    maxVoltageDrop = 1.0;
  }

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
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

  // now get the M and R local id's
  if( model_.includeMEquation )
  {
    li_MagVar = intLIDVec[ j++ ];
  }
  li_RVar   = intLIDVec[ j++ ];

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
    Xyce::dout() << " li_MagVar = " << li_MagVar << std::endl
         << " li_RVar = " << li_RVar << std::endl;
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

  if (model_.includeMEquation)
  {
    addInternalNode(symbol_table, li_MagVar, getName(), "m");
    addInternalNode(symbol_table, li_MagVar, getName().getEncodedName() + "_m");
  }
  addInternalNode(symbol_table, li_RVar, getName(), "r");
  addInternalNode(symbol_table, li_RVar, getName().getEncodedName() + "_r");

  addStateNode(symbol_table, li_MagVarState, getName(), "ms");
  addStateNode(symbol_table, li_MagVarDerivState, getName(), "dmdt");

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
  int i = 0;

  li_MagVarState = staLIDVec[i++];
  li_MagVarDerivState = staLIDVec[i++];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Instance::registerStateLIDs-------------------" << std::endl;

    Xyce::dout() << "li_MagVarState = " << li_MagVarState << std::endl
      << "li_MagVarDerivState = " << li_MagVarDerivState << std::endl
      ;
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

  li_RVarStore = stoLIDVec[0];
  li_BVarStore = stoLIDVec[1];
  li_HVarStore = stoLIDVec[2];
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
  // int numInductors = instanceData.size();  // don't need this as it's defined at class level
  int i = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->APosEquBraVarOffset  = jacLIDVec[ 2*i     ][ 0 ];
    (*currentInductor)->ANegEquBraVarOffset  = jacLIDVec[ 2*i + 1 ][ 0 ];
    (*currentInductor)->vPosOffset = jacLIDVec[ 2*numInductors + i ][ 0 ];
    (*currentInductor)->vNegOffset = jacLIDVec[ 2*numInductors + i ][ 1 ];
    int extraOffset = 2;
    if( i == 0)
    {
      extraOffset = 0;
    }
    (*currentInductor)->ABraEquPosNodeOffset = jacLIDVec[ 2*numInductors + i ][ 0 + extraOffset ];
    (*currentInductor)->ABraEquNegNodeOffset = jacLIDVec[ 2*numInductors + i ][ 1 + extraOffset ];
    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        (*currentInductor)->ABraEquBraVarOffset  = jacLIDVec[ 2*numInductors + i ][ j + 2 + extraOffset ];
      }
      (*currentInductor)->inductorCurrentOffsets[ j ] = jacLIDVec[ 2*numInductors + i ][ j + 2 + extraOffset ];
    }
    if( model_.includeMEquation )
    {
      (*currentInductor)->magOffset = jacLIDVec[ 2*numInductors + i ][ numInductors + 2 + extraOffset ];
    }
    currentInductor++;
    i++;
  }
  
  int rOffset=0; 
  if( model_.includeMEquation )
  {
    // now get the M equation offsets
    mEquVPosOffset = jacLIDVec[ 3*numInductors ][0];
    mEquVNegOffset = jacLIDVec[ 3*numInductors ][1];
    for( i=0; i<numInductors; ++i )
    {
      mEquInductorOffsets[i] = jacLIDVec[ 3*numInductors ][ i + 2];
    }
    mEquMOffset = jacLIDVec[ 3*numInductors ][ numInductors + 2 ];
    mEquROffset = jacLIDVec[ 3*numInductors ][ numInductors + 3 ];

    rOffset=1;
  }

  // now get the R equation offsets
  for( i=0; i<numInductors; ++i )
  {
    rEquInductorOffsets[i] = jacLIDVec[ 3*numInductors + rOffset ][ i ];
  }
  rEquROffset = jacLIDVec[ 3*numInductors + rOffset ][ numInductors ];

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

    Xyce::dout() << "mEquVPosOffset = " << mEquVPosOffset << "\tmEquVNegOffset = " << mEquVNegOffset << std::endl;
    Xyce::dout() << "mEquInductorOffsets = ";
    for(i=0;i<numInductors; ++i)
    {
      Xyce::dout() << mEquInductorOffsets[i] << ", ";
    }
    Xyce::dout() << std::endl
      << "mEquMOffset = " << mEquMOffset << "\tmEquROffset = " << mEquROffset  << std::endl;

    Xyce::dout() << "rEquInductorOffsets = ";
    for(i=0;i<numInductors; ++i)
    {
      Xyce::dout() << rEquInductorOffsets[i] << ", ";
    }
    Xyce::dout() << std::endl
      << "rEquROffset = " << rEquROffset << std::endl;
  }
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
    Xyce::dout() << "Instance::updateIntermediateVars " << std::endl;
  }

  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);

  // some parameters in the model class that we will use often
  const double A      = model_.A;
  const double Alpha  = model_.Alpha;
  const double Area   = model_.Area;
  const double BetaH  = model_.BetaH;
  const double C      = model_.C;
  const double DeltaVScaling = model_.DeltaVScaling;
  const double Gap    = model_.Gap;
  const double Ms     = model_.Ms;
  const double Kirr   = model_.Kirr;
  const double Path   = model_.Path;

  const double mVarScaling = model_.mVarScaling;
  const double rVarScaling = model_.rVarScaling;

  // calculate the voltage drop over the first inductor
  // as this is needed later
  double Vpos = solVector[(instanceData[0])->li_Pos];
  double Vneg = solVector[(instanceData[0])->li_Neg];

  // voltage drop over first inductor.
  double voltageDrop= Vpos - Vneg;

  // only update maxVoltageDrop when system has converged or we may
  // get wildly wrong values.
  Linear::Vector & lastSolVector = *(extData.currSolVectorPtr);
  double lastVoltageDrop = lastSolVector[(instanceData[0])->li_Pos] - lastSolVector[(instanceData[0])->li_Neg];
  if ( (getSolverState().newtonIter == 0) && (fabs(lastVoltageDrop) > maxVoltageDrop) )
  {
    maxVoltageDrop=fabs(lastVoltageDrop);
  }

  // approximate the sgn( voltageDrop ) with
  // tanh ( scalefactor * voltageDrop / maxVoltageDrop )
  if( model_.UseConstantDeltaVScaling )
  {
    qV = DeltaVScaling * voltageDrop;
  }
  else
  {
    qV = DeltaVScaling * voltageDrop / maxVoltageDrop;
  }

  double tanh_qV = 0.0;

  if ( (fabs(qV) < CONSTTANH_THRESH) )
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

  Happ = 0.0;
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int il=0;
  while( currentInductor != endInductor )
  {
    if( (getSolverState().dcopFlag) && ((*currentInductor)->ICGiven) )
    {
      Happ += ((*currentInductor)->IC) * inductanceVals[ il ];
    }
    else
    {
      Happ += solVector[(*currentInductor)->li_Branch] * inductanceVals[ il ];
    }
    il++;
    currentInductor++;
  }
  Happ /= Path;

  double latestMag=0.0;
  if( model_.includeMEquation )
  {
    latestMag = mVarScaling * solVector[ li_MagVar ];
  }
  else
  {
    latestMag = mVarScaling * staVector[ li_MagVarState ];
  }

  if( model_.factorMS )
  {
    latestMag *= Ms;
  }

  double H = Happ - (Gap / Path) * latestMag;

  He = H + Alpha * latestMag;

  Heo = BetaH*A;

  // terms that come up frequently
  const double gap_path = Gap / Path;
  const double He2 = He*He;
  const double Heo2 = Heo*Heo;
  const double sq_Heo2He2 = sqrt(Heo2 + He2);

  delM0 = model_.BetaM * Ms;
  double Man = Ms * He / ( A + sq_Heo2He2 );
  delM = Man - latestMag;

  // terms that come up frequently
  const double delM2 = delM*delM;
  const double delM02 = delM0*delM0;
  const double sq_delM02delM2 = sqrt( delM02 + delM2 );

  if( model_.factorMS )
  {
    Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
    Manp =  Ms * (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
    
    /*
    double x = He / A;
    if( x == 0 )
    {
    Manp = 0.0;
    }
    else
    {
    Manp = Ms *
      (((1.0+std::exp(-2.0* (x+1.0e-6))) / (1.0 - std::exp(-2.0* (x+1.0e-6))) - 1.0/(x+1.0e-6))-
      ((1.0+std::exp(-2.0* (x-1.0e-6))) / (1.0 - std::exp(-2.0* (x-1.0e-6))) - 1.0/(x-1.0e-6))) / 2.0e-6;
    }
    */
      
    P = ( C * Manp + (1 - C) * Mirrp) / ((1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp)*Ms);
  }
  else
  {
    Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
    Manp =  Ms * (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
    
    /*
    double x = He / A;
    if( x== 0)
    {
    Manp = 0.0;
    }
    else
    {
    Manp = Ms *
      (((1.0+std::exp(-2.0* (x+1.0e-6))) / (1.0 - std::exp(-2.0* (x+1.0e-6))) - 1.0/(x+1.0e-6))-
       ((1.0+std::exp(-2.0* (x-1.0e-6))) / (1.0 - std::exp(-2.0* (x-1.0e-6))) - 1.0/(x-1.0e-6))) / 2.0e-6;
    }
    */
    
    P = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);
  }


  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS)  && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\tA = " << A << std::endl
         << "\tArea = " << Area << std::endl
         << "\tPath = " << Path << std::endl
         << "\tGap = " << Gap << std::endl
         << "\tC = " << C << std::endl
         << "\tVpos = " << Vpos << std::endl
         << "\tVneg = " << Vneg << std::endl
         << "\tvoltageDrop = " << voltageDrop << std::endl
         << "\tqV = " << qV << std::endl
         << "\tdelM0 = " << delM0 << std::endl
         << "\tHapp = " << Happ
         << "\tlatestMag = " << latestMag
         << "\tlatestR = " <<  rVarScaling * solVector[ li_RVar ] << std::endl
         << "\tHe = " << He << std::endl
         << "\tH = " << H << std::endl
         << "\tHeo = " << Heo << std::endl
         << "\tMan = " << Man << std::endl
         << "\tdelM = " << delM << std::endl
         << "\tMirrp = " << Mirrp << std::endl
         << "\tManp = " << Manp << std::endl
         << "\tP  = " << P << std::endl
         << "\tgetSolverState().newtonIter = " << getSolverState().newtonIter << std::endl
         << std::endl;
  }

  // now calculate important derivative quantities

  double dHe_dM =  ((Alpha - gap_path) * mVarScaling);

  double dManp_dM = ( -Ms * He / (pow(A + sq_Heo2He2, 2.0)*sq_Heo2He2)) *
                    ( (Heo2 / (Heo2 + He2)) + (2.0*(A + Heo2 / sq_Heo2He2)/(A+sq_Heo2He2)) ) * dHe_dM;

  double ddelM_dM = ( dHe_dM*Ms/(A + sq_Heo2He2) ) * (1.0 - He2 / ((A + sq_Heo2He2)*sq_Heo2He2)) - mVarScaling;

  double dMirrp_dM = (1.0/(2.0*(Kirr - Alpha*sq_delM02delM2))) *
                     (tanh_qV + delM/sq_delM02delM2 +
                       (2.0*Alpha*delM*(delM*tanh_qV + sq_delM02delM2)
                     /(2.0*(Kirr-Alpha*sq_delM02delM2)*sq_delM02delM2))) * ddelM_dM;

  double dP_Denom=0.0;
  if( model_.factorMS )
  {
    dP_Denom = 1.0 + (gap_path - Alpha)*C*Manp + gap_path * (1.0-C) * Mirrp;

    dP_dM = (1.0/dP_Denom) * (C * dManp_dM + (1.0-C) * dMirrp_dM) -
              ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
                ( (gap_path - Alpha)*C*dManp_dM + gap_path*(1.0-C)*dMirrp_dM );
    dP_dM /= Ms;
  }
  else
  {
    dP_Denom = 1.0 + (gap_path - Alpha)*C*Manp + gap_path * (1.0-C) * Mirrp;

    dP_dM = (1.0/dP_Denom) * (C * dManp_dM + (1.0-C) * dMirrp_dM) -
              ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
                ( (gap_path - Alpha)*C*dManp_dM + gap_path*(1.0-C)*dMirrp_dM );
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\tA = " << A << std::endl
      << "\tAlpha = " << Alpha << std::endl
      << "\tC = " << C << std::endl
      << "\tGap = " << Gap << std::endl
      << "\tMs = " << Ms << std::endl
      << "\tKirr = " << Kirr << std::endl
      << "\tPath = " << Path << std::endl
      << "\tHe2 = " << He2 << std::endl
      << "\tHeo2 = " << Heo2 << std::endl
      << "\tdelM2 = " << delM2 << std::endl
      << "\tdelM02 = " << delM02 << std::endl
      << "\tdHe_dM = " << dHe_dM << std::endl
      << "\tdManp_dM = " << dManp_dM << std::endl
      << "\tddelM_dM = " << ddelM_dM << std::endl
      << "\tdMirrp_dM = " << dMirrp_dM << std::endl
      << "\tdP_dM = " << dP_dM << std::endl
      << "\tdenom 1+(1-lg/lt)P = " << (1+(1-Gap/Path)*P) << std::endl;
  }

  //    % Now find (dP/dI_i): (this is nearly identical to dP/dM)
  currentInductor = instanceData.begin();

  for( int i=0; i<numInductors; ++i )
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
          ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
            ( (gap_path - Alpha)*C*dManp_dI[i] + gap_path*(1.0-C)*dMirrp_dI[i] );

  if( model_.factorMS )
  {
    dP_dI[i] /= Ms;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
       Xyce::dout() << "\tdHe_dI[ " << i << " ] =" << dHe_dI[ i ] << std::endl
            << "\tdManp_dI[ " << i << " ] = " << dManp_dI[i] << std::endl
            << "\tddelM_dI[ " << i << " ] = " << ddelM_dI[i] << std::endl
            << "\tMirrp_dI[ " << i << " ] = " << dMirrp_dI[i] << std::endl
            << "\tdP_dI[ " << i << " ] = " << dP_dI[i] << std::endl;
    }
    currentInductor++;
  }

  // Now find (dP/dV_1):
  double dMirrp_dVp = (delM * (DeltaVScaling/maxVoltageDrop) * (1.0-pow(tanh_qV,2.0))) /
                      (2.0 * (Kirr - Alpha * sq_delM02delM2));
  double dMirrp_dVn = -dMirrp_dVp;

  dP_dVp = (1.0/dP_Denom) * ((1.0-C) * dMirrp_dVp) -
            ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) * (  gap_path*(1.0-C)*dMirrp_dVp );

  if( model_.factorMS )
  {
    dP_dVp /= Ms;
  }

  dP_dVn = -dP_dVp;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\tdMirrp_dVp = " << dMirrp_dVp << std::endl
      << "\tdMirrp_dVn = " << dMirrp_dVn << std::endl
      << "\tdP_dVp = " << dP_dVp << std::endl
      << "\tdP_dVn = " << dP_dVn << std::endl;
  }

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
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState---------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  updateIntermediateVars ();

  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);
  Linear::Vector & stoVector = *(extData.nextStoVectorPtr);
  double mVarScaling = model_.mVarScaling;

  // place current values of mag, H and R in state vector
  // must unscale them as the rest of the class assumes
  // that these aren't scaled yet.
  if( model_.includeMEquation ) 
  {
    staVector[ li_MagVarState ] = solVector[ li_MagVar ];
  }
  stoVector[ li_RVarStore ] = solVector[ li_RVar ];

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

  // copy derivitive of Mag from result vector into state vector
  staVector[ li_MagVarDerivState ] = staDerivVec[ li_MagVarState ];
  
  double mVarScaling = model_.mVarScaling;
  // place current values of mag, H and R in state vector
  // must unscale them as the rest of the class assumes
  // that these aren't scaled yet.
  double latestMag = 0.0;
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
  stoVector[ li_RVarStore ] = solVector[ li_RVar ];
  
  
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
  double calculatedB = model_.BCgsFactor * (4.0e-7 * M_PI * (calculatedH + latestMag));

  double deltaH = calculatedH - lastH;
  double deltaB = calculatedB - lastB;
  double dBdH = 0;
  if( deltaH != 0.0)
  {
    dBdH = deltaB / deltaH;
  }
  
  double Hfxn = Happ  - (model_.Gap / model_.Path) * latestMag;
  if( (model_.Gap <= 0) && (dBdH < 0))
  {
    Hfxn = lastH / model_.HCgsFactor;
  }
  else
  {
    double dMdt = mVarScaling*staDerivVec[ li_MagVarState ];
    double R = solVector[ li_RVar ];
    double dHdt = R - (model_.Gap/model_.Path)*dMdt;
  
    if( (dMdt > 0) && (dHdt < 0) )
    {
      // B is increasing, so H should be increasing
      // ignore the M component. 
      Hfxn = Happ;
    }
    else if( (dMdt < 0) && (dHdt > 0) )
    {
      // B is decreasing, so H should be decreasing 
      // ignore the M component
      Hfxn = Happ;
    }
  }

  stoVector[ li_HVarStore ] = model_.HCgsFactor * Hfxn;
  stoVector[ li_BVarStore ] = model_.BCgsFactor * (4.0e-7 * M_PI * (stoVector[ li_HVarStore ] + latestMag));

  return bsuccess;
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
  double mVarScaling = model_.mVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEQVector------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector & staVector = *(extData.nextStaVectorPtr);
  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  double * qVec = extData.daeQVectorRawPtr;

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

    qVec[((*currentInductor)->li_Branch)] += LOI[ i ];

    double current = inductorCurrents[ i ];
    double windings = (*currentInductor)->L;

    qVec[ li_RVar ] += rEqScaling * current * windings;
    i++;
    currentInductor++;
  }
  
  // load M terms if needed 
  if( model_.includeMEquation )
  {
    double latestMag = mVarScaling * staVector[ li_MagVarState ];

    // M equation
    if(!getSolverState().dcopFlag)
    {
      qVec[ li_MagVar ] += mEqScaling * latestMag;
    }
  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
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
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;
  double mVarScaling = model_.mVarScaling;
  double rVarScaling = model_.rVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEFVector------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector & staVector = *(extData.nextStaVectorPtr);
  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & stoVector = *(extData.nextStoVectorPtr);

  double * fVec = extData.daeFVectorRawPtr;

  double latestR = rVarScaling * stoVector[ li_RVarStore ];
  
  if(getSolverState().dcopFlag)
  {
    //enforce R = 0 in dc op
    latestR = 0.0;
  }
  // load M terms if needed 
  if( model_.includeMEquation )
  { 
    // for the M equation
    fVec[li_MagVar] -= mEqScaling *  P * latestR / (model_.Path);

    // if |P| is near zero, then the M equation becomes dM/dt = 0, or M is
    // constant.  In this case we'll add a diagonal element for M so that
    // sole dM/dt element in dQ/dX doesn't cause a time step too small error
    // Since P is normally very large, we'll test for |P| <= 1.0.
  
    if( fabs( P ) <= model_.pZeroTol )
    {
      fVec[li_MagVar] -= mVarScaling * staVector[ li_MagVarState ];
    }
  } 

  // for the R equation
  fVec[li_RVar] -= rEqScaling * rVarScaling * stoVector[ li_RVarStore ];

  // used in scaling the branch equation;
  double mid=1.0;
  if( model_.factorMS )
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P*(model_.Ms);
  }
  else
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P;
  }

  // loop over each inductor and load it's F vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    double current   = solVector[(*currentInductor)->li_Branch];
    double ic_coef = 1.0;
    if( (getSolverState().dcopFlag) && (*currentInductor)->ICGiven == true )
    {
      current = (*currentInductor)->IC;
      ic_coef=0.0;
      solVector[(*currentInductor)->li_Branch] = current;
    }
    double vNodePos  = solVector[(*currentInductor)->li_Pos];
    double vNodeNeg  = solVector[(*currentInductor)->li_Neg];


    fVec[((*currentInductor)->li_Pos)]    +=  scalingRHS * current;

    fVec[((*currentInductor)->li_Neg)]    += -scalingRHS * current;

    fVec[((*currentInductor)->li_Branch)] += -ic_coef*((vNodePos - vNodeNeg)/mid);
    double windings = (*currentInductor)->L;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  Inductor = " << (*currentInductor)->name
           << " li_Pos = " << (*currentInductor)->li_Pos
           << " li_Neg = " << (*currentInductor)->li_Neg
           << " li_Branch = " << (*currentInductor)->li_Branch
           << "\tPos/Neg current*windings = " << scalingRHS*current*windings
           << "\tBranch = " << ((vNodePos - vNodeNeg)/mid)
           << std::endl;
    }

    if (loadLeadCurrent)
    {
      double * leadF = extData.nextLeadCurrFCompRawPtr;     
      double * junctionV = extData.nextJunctionVCompRawPtr;
      leadF[(*currentInductor)->li_branch_data] =  scalingRHS * current;       
      junctionV[(*currentInductor)->li_branch_data] = (vNodePos - vNodeNeg);
    }

    currentInductor++;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single instance.
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  double mVarScaling = model_.mVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

  int i;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdQdx-----------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;
  // update the M equation if it's needed
  if( model_.includeMEquation )
  { 
    // update M equation
    if(!getSolverState().dcopFlag)
    {
      (*dQdxMatPtr)[li_MagVar][mEquMOffset] += mEqScaling * mVarScaling;
    }
  }

  // update the R equation
  for( i = 0; i< numInductors; i++ )
  {
    (*dQdxMatPtr)[li_RVar][rEquInductorOffsets[i] ] += rEqScaling * inductanceVals[i];
  }

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
// Purpose       : Loads the F-vector contributions for a single instance.
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
  double rVarScaling = model_.rVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdFdx----------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }

  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & stoVector = *(extData.nextStoVectorPtr);
  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  // udate dependent parameters
  //updateIntermediateVars();

  // pull these parameters up from the model class to make it easier
  // to view the equations.
  const double Gap = model_.Gap;
  const double Path = model_.Path;

  // terms that come up frequently
  double latestR = rVarScaling * stoVector[ li_RVarStore ];

  // inlcude M equation if it's part of the system
  if( model_.includeMEquation )
  { 
    // terms for the M equation
    if(!getSolverState().dcopFlag)
    {
      (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ]    -= mEqScaling * dP_dM * latestR / Path;   // d/dM
      (*dFdxMatPtr)[ li_MagVar ][ mEquROffset ]    -= mEqScaling * P * rVarScaling / Path;   // d/dR
      (*dFdxMatPtr)[ li_MagVar ][ mEquVPosOffset ] -= mEqScaling * dP_dVp * latestR / Path;  // d/dV_+
      (*dFdxMatPtr)[ li_MagVar ][ mEquVNegOffset ] -= mEqScaling * dP_dVn * latestR / Path;  // d/dV_-
      for( int i = 0; i<numInductors; ++i)
      {
        (*dFdxMatPtr)[ li_MagVar ][mEquInductorOffsets[i] ] -=
           mEqScaling * dP_dI[i] * latestR / Path;                     // d/dI_i;
      }
    }
    else
    {
      // the above load for the M equation is basically zero in the dc op.  We
      // need something on the diagonal for M to make the matrix non-singular
      (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ]    += getSolverState().pdt_;
    }

    // if |P| is near zero, then the M equation becomes dM/dt = 0, or M is
    // constant.  In this case we'll add a diagonal element for M so that
    // sole dM/dt element in dQ/dX doesn't cause a time step too small error
    // Since P is normally very large, we'll test for |P| <= 1.0.
    if( fabs( P ) <= model_.pZeroTol )
    {
      (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ] += 1.0;
    }
  }

  // update the R equation

  (*dFdxMatPtr)[ li_RVar ][rEquROffset] -= rEqScaling * rVarScaling;

  // loop over each inductor and load it's dFdx components
  double mid=0.0;
  if( model_.factorMS )
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P*(model_.Ms);
  }
  else
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P;
  }

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    // do the normal work for an inductor
    if( (getSolverState().dcopFlag) && (*currentInductor)->ICGiven == true )
    {
      (*dFdxMatPtr)[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  += 0.0;
      (*dFdxMatPtr)[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += 0.0;
      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += 0.0;
      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] += 0.0;
      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquBraVarOffset)]  += 0.0;
    }
    else
    {
      (*dFdxMatPtr)[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  +=  scalingRHS;
      (*dFdxMatPtr)[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += -scalingRHS;
      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += -1.0/mid;
      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] +=  1.0/mid;
    }

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout()
       << "(*currentInductor)->li_Pos = " << (*currentInductor)->li_Pos << std::endl
       << "(*currentInductor)->li_Neg = " << (*currentInductor)->li_Neg << std::endl
       << "(*currentInductor)->li_Branch = " << (*currentInductor)->li_Branch << std::endl
       << "(*currentInductor)->APosEquBraVarOffset = " << (*currentInductor)->APosEquBraVarOffset << std::endl
       << "(*currentInductor)->ANegEquBraVarOffset = " << (*currentInductor)->ANegEquBraVarOffset << std::endl
       << "(*currentInductor)->ABraEquPosNodeOffset = " << (*currentInductor)->ABraEquPosNodeOffset << std::endl
       << "(*currentInductor)->ABraEquNegNodeOffset = " << (*currentInductor)->ABraEquNegNodeOffset << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Pos)<<"]   ["<<((*currentInductor)->APosEquBraVarOffset)<<"] =  " << scalingRHS << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Neg)<<"]   ["<<((*currentInductor)->ANegEquBraVarOffset)<<"]  =  " << -scalingRHS << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Branch)<<"]["<<((*currentInductor)->ABraEquPosNodeOffset)<<"] = " << -1/mid << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Branch)<<"]["<<((*currentInductor)->ABraEquNegNodeOffset)<<"] = " << 1/mid << std::endl;
    }

    double delV = solVector[(*currentInductor)->li_Pos] - solVector[(*currentInductor)->li_Neg];

    for( int j = 0; j<numInductors; ++j )
    {

      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] +=
        delV * (1.0 - (Gap/Path)) * dP_dI[j]/(mid*mid);

      if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "(*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] =  " << delV * (1 - (Gap/Path)) * dP_dI[j]/(mid*mid) << std::endl;
        }
    }
    if( model_.includeMEquation )
    { 
      (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->magOffset]  += delV * (1.0 - (Gap/Path)) * dP_dM/(mid*mid);
    }
    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vPosOffset] += delV * (1.0 - (Gap/Path)) * dP_dVp/(mid*mid);

    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vNegOffset] += delV * (1.0 - (Gap/Path)) * dP_dVn/(mid*mid);

    /*
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "(*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->magOffset] =  " << delV * (1 - (Gap/Path)) * dP_dM/(mid*mid) << std::endl
       << "(*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vPosOffset]  =  " << delV * (1 - (Gap/Path)) * dP_dVp/(mid*mid) << std::endl
       << "(*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vNegOffset] = " << delV * (1 - (Gap/Path)) * dP_dVn/(mid*mid) << std::endl;
    }
    */
    currentInductor++;
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
    Linear::Vector & staVector = *(extData.nextStaVectorPtr);
    Linear::Vector & stoVector = *(extData.nextStoVectorPtr);
    double mVarScaling = model_.mVarScaling;
    double rVarScaling = model_.rVarScaling;

    double latestMag = mVarScaling * staVector[ li_MagVarState ];
    if( model_.factorMS )
    {
      latestMag *= model_.Ms;
    }
    double latestR   = rVarScaling * stoVector[ li_RVarStore ];
    (*outputFileStreamPtr)
      << getSolverState().currTime_ << "  "
      << latestMag << "\t  "
      << latestR << "\t "
      << staVector[ li_BVarStore ] << "\t "
      << staVector[ li_HVarStore ]
      << std::endl;

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
  double * nextSolVector = extData.nextSolVectorRawPtr;
  double * currSolVector = extData.currSolVectorRawPtr;

  // loop over each inductor and load it's dFdx components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    if ((*currentInductor)->ICGiven)
    {
      currSolVector[(*currentInductor)->li_Branch] = (*currentInductor)->IC;
      nextSolVector[(*currentInductor)->li_Branch] = (*currentInductor)->IC;
    }
    currentInductor++;
  }
  
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
  varTypeVec.resize(numInductors+2);
  for(int i=0; i<numInductors; i++)
  {
    varTypeVec[i] = 'I';
  }
  // I don't know what should be used for non I,V vars.
  varTypeVec[numInductors] = 'I';
  varTypeVec[numInductors+1] = 'I';
}

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
  // when Ms factoring is off, scaling of M/R is still needed.
  // unless the user has specified something.  We don't want to 
  // override their settings
  if( !factorMSGiven && !mVarScalingGiven ) 
  {
    mVarScaling=1.0e3;
  }
  
  if( !factorMSGiven && !rVarScalingGiven ) 
  {
    rVarScaling=1.0e3;
  }
  
  if( !factorMSGiven && !mEqScalingGiven ) 
  {
    mEqScaling=1.0e-3;
  }
  
  if( !factorMSGiven && !mEqScalingGiven ) 
  {
    rEqScaling=1.0e-3;
  }
  
  if( factorMSGiven && !rVarScalingGiven )
  {
    rVarScaling=Ms;
  }

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
    DeltaVScaling(0.0),
    GapInCm(0.0),
    Gap(0.0),
    Kirr(0.0),
    Ms(0.0),
    PathInCm(0.0),
    Path(0.0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tnom(getDeviceOptions().tnom),
    HCgsFactor( 0.012566370614359 ),  // 4 pi / 1000
    BCgsFactor( 10000.0 ),
    outputStateVars(0),
    factorMS(0),
    UseConstantDeltaVScaling( false )
{
  setLevel(1);


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


//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->updatePrimaryState ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState (double * staDerivVec, double * stoVec)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->updateSecondaryState ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->loadDAEFVector();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEQVector();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->loadDAEdFdx ();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEdQdx ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("min", 1)
    .registerModelType("min", 1)
    .registerModelType("core", 1);
}

} // namespace MutIndNonLin
} // namespace Device
} // namespace Xyce

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
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 06/10/09
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------
// for random() function
#include <cstdlib>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_MembraneCS.h>
#include <N_DEV_MembraneHH.h>
#include <N_DEV_MembranePassive.h>
#include <N_DEV_MembraneUserDefined.h>
#include <N_DEV_NeuronPop1.h>
#include <N_DEV_Neuron_CommonEquations.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace NeuronPop1 {

void Traits::loadInstanceParameters(ParametricData<NeuronPop1::Instance> &p)
{
  p.addPar ("CTP",std::vector<std::string>(),&NeuronPop1::Instance::connectionTargetPopulation)
   .setGivenMember(&NeuronPop1::Instance::connectionTargetPopulationGiven)
   .setUnit(U_NONE)
   .setCategory(CAT_NONE)
   .setDescription("Connected Target Population list");
}

void Traits::loadModelParameters(ParametricData<NeuronPop1::Model> &p)
{
  p.addPar ("NEURONS_MAX",10,&NeuronPop1::Model::neuronsMax)
    .setGivenMember(&NeuronPop1::Model::neuronsMaxGiven)
    .setDescription("Maximum number of neurons in the device");

  p.addPar ("IC_MAX",2,&NeuronPop1::Model::internalMaxConnections)
    .setGivenMember(&NeuronPop1::Model::internalMaxConnectionsGiven)
    .setDescription("Maximum number of internal connections in the device");

  p.addPar ("EC_MAX",2,&NeuronPop1::Model::externalMaxConnections)
    .setGivenMember(&NeuronPop1::Model::externalMaxConnectionsGiven)
    .setDescription("Maximum number of external connections in the device");

  p.addPar ("NEUROGENESIS_RATE",0.0,&NeuronPop1::Model::populationNeurogenesisRate)
    .setGivenMember(&NeuronPop1::Model::populationNeurogenesisRateGiven)
    .setUnit(U_SECOND)
    .setDescription("Rate in days of GC neurogenesis in the population");

  p.addPar ("UPDATE_PERIOD",1.0,&NeuronPop1::Model::populationUpdatePeriod)
    .setGivenMember(&NeuronPop1::Model::populationUpdatePeriodGiven)
    .setUnit(U_SECOND)
    .setDescription("Time in days for population updates");

  p.addPar ("OUTPUTPOPULATIONVARS",0,&NeuronPop1::Model::outputPopulationVars)
    .setDescription("Flag to save population variables" );
}



//
// static class member inits
//

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    liNodeIn(-1),
    liNodeOut(-1),
    populationInitialized(false),
    lastPopulationUpdateTime(0.0),
    lastNeurogenesisUpdateTime(0.0),
    neuronPopSize(0),
    connectionTargetPopulationGiven(false),
    outputPopulationVarsFlag(false),
    newStateToOutput(false),
    numberOfUpdatesDone(0),
    jsOffsetNodeIn(-1),
    jsOffsetNodeOut(-1)
{
  // Set up the number of internal, external and state vars
  numExtVars   = 2;  // input and output voltage
  // numIntVars number of internal vars
  // currently we have 3 variables per neuron, Voltage, X-position, Y-Position plus internal and external connections (maximum)
  // These will all be in the state vector as: V1, V2, ..., Vn, X1, X2, ..., Xn, Y1, ..., Yn
  // We do this rather than have individual containers because it will be easier to push the
  // state vector around for working with these in parallel.

  // number of vars per neuron: 3
  numStateVars = (3 + model_.internalMaxConnections + model_.externalMaxConnections) * model_.neuronsMax;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // total up number of vars.
  int numVars = numExtVars + numIntVars;

  // need to work out what if-any real jac stap I need here.


  // set up jacStamp.  This is dependant on the membrane model.  The only part this
  // constructor really knows about is the external variables Vin and vOut

  if( jacStamp.empty() )       // redundant as jacStamp is not static for this device
  {                            // it can't be as each cable may have a different number of nodes
    jacStamp.resize(numVars);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 0;                               // NodeIn
    jacStamp[1].resize(1);
    jacStamp[1][0] = 1;                               // NodeOut
  }

  // if the user has requested output of the state variables M, H and R
  // then open a file for that output.
  if( model_.outputPopulationVars > 0 )
  {
    outputPopulationVarsFlag = true;
    std::string filename( "NeuronPop_" );
    filename.append( getName().getEncodedName() );
    filename.append( ".dat" );
    // convert any embeded ':' or '%' characters to '_'
    replace( filename.begin(), filename.end(), '%', '_' );
    replace( filename.begin(), filename.end(), ':', '_' );

    outputFileStreamPtr = rcp( new std::ofstream() );
    outputFileStreamPtr->open( filename.c_str() );
    if( !(*outputFileStreamPtr) )
    {
      UserError(*this) << "Could not open file " << filename << " for output of population variables";
    }
    else {
      (*outputFileStreamPtr).setf(std::ios::scientific, std::ios::floatfield );
      (*outputFileStreamPtr).width(20);
      (*outputFileStreamPtr).precision(12);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                            const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  liNodeIn = extLIDVec[0];
  liNodeOut = extLIDVec[1];

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
{}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;
  // this resize won't be true when we store more data in the state
  // vector. (i.e. not just voltages, but connectivities firing state)
  liNeuronPopState.resize( numStateVars );
  for( int i=0; i<numStateVars; i++ )
  {
    liNeuronPopState[i] = staLIDVec[i];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  // resize our storage location and store the results

  int numRows = jacLIDVec.size();

  jsOffsetNodeIn = jacLIDVec[0][0];
  jsOffsetNodeOut = jacLIDVec[1][0];

}

//-----------------------------------------------------------------------------
// Function      : Instance::initializePopulation
// Purpose       : Steps in initializing a population of neurons
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::initializePopulation()
{
  // population lives in the state vector, so get a reference to that.
  // note that the state vector is packed in the form of: V1, V2, ..., Vn, X1, X2, ..., Xn, Y1, ..., Yn, internal connections, external connections
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);

  // seed the random number generator
  //srandom(unsigned seed)
  // figure out how many neurons we will start with
  const double random_max = std::pow(2.0,31)-1;  // this may be the same as RAND_MAX but it's not
                                                 // clear from the man page.  Should move random
                                                 // number generation to the utl package.
  // try to get a random population from 1..(neuronsMax-1) inclusive of ends
  neuronPopSize = (rand()/random_max) * model_.neuronsMax + 1;
  if( neuronPopSize >= model_.neuronsMax)
    neuronPopSize = model_.neuronsMax;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::initializePopulation "
              << "neuronsMax = " << model_.neuronsMax
              << " neuronPopSize = " << neuronPopSize
              << std::endl;
  }

  int maxPopSize = model_.neuronsMax;
  // now initialize the locations
  for(int i=0; i<maxPopSize; i++)
  {
    // zero out voltages
    staVector[ liNeuronPopState[                            i ] ] = 0.0;
    // randomize x positions
    staVector[ liNeuronPopState[    model_.neuronsMax + i ] ] = (rand()/random_max);
    // randomize y positions
    staVector[ liNeuronPopState[2 * model_.neuronsMax + i ] ] = (rand()/random_max);
  }

  // now initialize the internal and external connections
  for(int i=0; i<maxPopSize; i++)
  {
    // randomly connect internal connections
    int numConnections = (rand()/random_max) * model_.internalMaxConnections;
    int k = 0;
    for(; k<numConnections; k++)
    {
      int postNeuron = (rand()/random_max) * maxPopSize;
      staVector[ liNeuronPopState[3 * model_.neuronsMax + i * model_.internalMaxConnections + k] ] = postNeuron;
    }
    for(int j=k; j<model_.internalMaxConnections; j++)
    {
      staVector[ liNeuronPopState[3 * model_.neuronsMax + i * model_.internalMaxConnections + j] ] = 0.0;
    }

    // zero out external connections
    for(int j=0; j<model_.externalMaxConnections; j++)
    {
      staVector[ liNeuronPopState[(3 + model_.internalMaxConnections) * model_.neuronsMax + i * model_.externalMaxConnections + j] ] = 0.0;
    }
  }

  // now rescale update period and neurogenesis rate in terms of seconds (instead of days)
  model_.populationUpdatePeriod = 24*60*60 * model_.populationUpdatePeriod;
  model_.populationNeurogenesisRate = 24*60*60 / model_.populationNeurogenesisRate;

  populationInitialized=true;
  // this flag is to signal to outputPlotFiles() function that the state of the system has
  // changed and should be output at the next call back
  newStateToOutput=true;
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePopulation
// Purpose       : Steps in updating a population of neurons
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 03/23/11
//-----------------------------------------------------------------------------
void Instance::updatePopulation()
{
  // population lives in the state vector, so get a reference to that.
  // note that the state vector is packed in the form of: V1, V2, ..., Vn, X1, X2, ..., Xn, Y1, ..., Yn
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePopulation "
              << "time = " << getSolverState().currTime_
              << std::endl;
  }

  double time = getSolverState().currTime_;

  // this check ensures that we've progressed far enough for another neurogenesis event to occur
  // and we will not exceed the maximum number of neurons allowed
  if ((model_.populationNeurogenesisRate > 0) &&
      (fabs(lastNeurogenesisUpdateTime - time) >= model_.populationNeurogenesisRate) &&
      (neuronPopSize < model_.neuronsMax))
  {
    neuronPopSize++;
    lastNeurogenesisUpdateTime=time;

    // this flag is to signal to outputPlotFiles() function that the state of the system has
    // changed and should be output at the next call back
    newStateToOutput=true;
  }

  lastPopulationUpdateTime=getSolverState().currTime_;
  numberOfUpdatesDone++;
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/22/2011
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints (std::vector<Util::BreakPoint> &breakPointTimes)
{
  // push on to the vector the next two update times
  breakPointTimes.push_back((numberOfUpdatesDone+1) * model_.populationUpdatePeriod);
  breakPointTimes.push_back((numberOfUpdatesDone+2) * model_.populationUpdatePeriod);

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;
  if( !populationInitialized )
  {
    initializePopulation();
  }


  double time = getSolverState().currTime_;

  // only do the population updates under the following conditions
  //   1. when we're at one of the population update times
  //   2. current time is after the last update time.
  //
  // I was going to do this like updateSource() in voltage/current sources, but I can't because
  // this device isn't derived from that type.  This is abit of a hack, but I'll
  // see if this works for now.

  // this check ensures that we're at least bpTol past the lastPopulationUpdateTime
  if( fabs(lastPopulationUpdateTime - time) > getSolverState().bpTol_)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  Instance::updateIntermediateVars\n";
    }



    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  Time = " << time << std::endl;
    }
    // subtract off any delay time
    // time -= TD;

    if (time > model_.populationUpdatePeriod && model_.populationUpdatePeriod != 0.0)
    {
      // repeating signal - figure out where we are in period
      double basetime = model_.populationUpdatePeriod * floor(time/model_.populationUpdatePeriod);
      time -= basetime;
    }

    // The basetime correction above could take a time right at 2 perionds and make time=0
    // Thus, we check if time is zero or if it's within bpTol of a period to find out if we
    // update the population

    if (fabs(time) < getSolverState().bpTol_ || (fabs(time - model_.populationUpdatePeriod) < getSolverState().bpTol_) )
    {
      // we're at a point to update the population
      updatePopulation();
    }
  }

  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  updateIntermediateVars();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Vector * daeQVecPtr = extData.daeQVectorPtr;

  // no Q component for the cable component of this devcie

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Vector * daeFVecPtr = extData.daeFVectorPtr;

  // take care of the input and output nodes as they are different
  (*daeFVecPtr)[liNodeIn]  += 0.0;
  (*daeFVecPtr)[liNodeOut]  += 0.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  Linear::Vector * solVectorPtr = extData.nextSolVectorPtr;
  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[liNodeIn][jsOffsetNodeIn]  +=  1.0;
  (*dFdxMatPtr)[liNodeOut][jsOffsetNodeOut]  +=  1.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       : If requested by the user output all the variables
//                 associated with the population
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/10/2011
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles(bool force_final_output)
{

  bool bsuccess = true;

  // only output when the population has changed to avoid lots of duplicate output
  // the newStateToOutput flag is set in initializePopulation() and updatePopulation()
  if( outputPopulationVarsFlag && newStateToOutput && outputFileStreamPtr.get() && (*outputFileStreamPtr) )
  {
    // population lives in the state vector, so get a reference to that.
    // note that the state vector is packed in the form of: V1, V2, ..., Vn, X1, X2, ..., Xn, Y1, ..., Yn
    Linear::Vector & solVector = *(extData.nextSolVectorPtr);
    Linear::Vector & staVector = *(extData.nextStaVectorPtr);
    Linear::Vector & staDerivVec = *(extData.nextStaDerivVectorPtr);

    // output format is
    // time, population_size, x1, ..., xn, y1, ..., yn, v1, ..., vn
    //
    (*outputFileStreamPtr)
      << getSolverState().currTime_ << ", "
      << neuronPopSize << ", ";

    // the state vector is packed in the order
    for( int i=0; i < neuronPopSize; i++)
    {
      (*outputFileStreamPtr) << staVector[ liNeuronPopState[ model_.neuronsMax + i]] << ", ";
    }

    for( int i=0; i < neuronPopSize; i++)
    {
      (*outputFileStreamPtr) << staVector[ liNeuronPopState[ 2*model_.neuronsMax + i]] << ", ";
    }

    for( int i=0; i < neuronPopSize; i++)
    {
      (*outputFileStreamPtr) << staVector[ liNeuronPopState[i] ] << ", ";
    }

    for( int i=0; i < neuronPopSize; i++)
    {
      for( int j=0; j < model_.internalMaxConnections; j++)
      {
        (*outputFileStreamPtr) << staVector[ liNeuronPopState[ 3*model_.neuronsMax + i*model_.internalMaxConnections + j]] << ", ";
      }
    }

    for( int i=0; i < neuronPopSize; i++)
    {
      for( int j=0; j < model_.externalMaxConnections; j++)
      {
        (*outputFileStreamPtr) << staVector[ liNeuronPopState[ (3 + model_.internalMaxConnections)*model_.neuronsMax + i*model_.externalMaxConnections + j]];
        if( (i != (neuronPopSize-1)) || (j != (model_.externalMaxConnections-1)) )
        {
          (*outputFileStreamPtr) << ", ";
        }
      }
    }
    (*outputFileStreamPtr) << std::endl;

    // we've output this state, so reset flag
    newStateToOutput=false;

  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}


//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    neuronsMax(0),
    neuronsMaxGiven(false),
    internalMaxConnections(0),
    internalMaxConnectionsGiven(false),
    externalMaxConnections(0),
    externalMaxConnectionsGiven(false),
    populationNeurogenesisRate(0.0),
    populationNeurogenesisRateGiven(false),
    populationUpdatePeriod(0.0),
    populationUpdatePeriodGiven(false),
    outputPopulationVars(0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Neuron instances: " << isize << std::endl;
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

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("NEURONPOP")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("neuronpop", 1)
      .registerModelType("neuronpop", 1);
  }
}

} // namespace NeuronPop1
} // namespace Device
} // namespace Xyce

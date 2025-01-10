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

//----------------------------------------------------------------------------
//
// Purpose        : This file implements the ADC digital to analog conversion
//                  device used in the integration of Xyce withe SAVANT VHDL
//                  simulator.
//
// Special Notes  : The ADC looks to the analog simulation like a resistor.
//
// Creator        : Lon Waters
//
// Creation Date  : 07/26/2002
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Standard Includes ----------
#include <algorithm>
#include <N_UTL_Math.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_ADC.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {


namespace ADC {


void Traits::loadInstanceParameters(ParametricData<ADC::Instance> &p)
{
  p.addPar ("R", 1.e+12, &ADC::Instance::R)
    .setUnit(U_OHM)
    .setDescription("internal Resistance");

  p.addPar ("WIDTH", 1, &ADC::Instance::outputBitVectorWidth_)
    .setGivenMember(&ADC::Instance::outputBitVectorWidthGiven_)
    .setUnit(U_NONE)
    .setDescription("Output bit vector width");
}

void Traits::loadModelParameters(ParametricData<ADC::Model> &p)
{
  p.addPar ("LOWERVOLTAGELIMIT", 0.0, &ADC::Model::lowerVoltageLimit_)
    .setUnit(U_VOLT)
    .setDescription("Lower limit of ADC voltage range");

  p.addPar ("UPPERVOLTAGELIMIT", 5.0, &ADC::Model::upperVoltageLimit_)
    .setUnit(U_VOLT)
    .setDescription("Upper limit of ADC voltage range");

  p.addPar ("SETTLINGTIME", 1.0e-8, &ADC::Model::settlingTime_)
    .setUnit(U_SECOND)
    .setDescription("Settling time");

  p.addPar ("WIDTH", 1, &ADC::Model::outputBitVectorWidth_)
    .setUnit(U_NONE)
    .setDescription("Output bit vector width");
}



std::vector< std::vector<int> > Instance::jacStamp;

//----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Instance::processParams ()
{
  G= (R != 0.0) ? 1.0/R : 0.0;

  // set the output bit vector width and number of quantization levels.
  if (!outputBitVectorWidthGiven_)
    outputBitVectorWidth_ = model_.outputBitVectorWidth_;
  setNumberQuantLevels();
  
  return true;
}

//----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    R(1.0e12),
    G(0.0),
    i0(0.0),
    outputBitVectorWidth_(1),
    outputBitVectorWidthGiven_(false),
    nQuantLevels_(0),
    lastOutputLevel_(0),
    li_Pos(-1),
    li_Neg(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    li_store_output_state(-1)
{
  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 0;
  setNumStoreVars(1);  // Device state will be copied to store vector for output.

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[1].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // Calculate dependent (ie computed) params and check for errors.
  // Note the quantization levels are set in processParams
  processParams ();
}

//----------------------------------------------------------------------------
// Function       : Instance::~Instance
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Instance::~Instance()
{
}

//----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  ADCInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;

  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    Xyce::dout() << section_divider << std::endl;
}

//----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 4/22/2019
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // copy over the global ID lists.
  stoLIDVec = stoLIDVecRef;

  li_store_output_state = stoLIDVec[0];
}

//----------------------------------------------------------------------------
// Function       : jacobianStamp
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 11/12/2002
//----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//----------------------------------------------------------------------------
// Function       : registerJacLIDs
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 11/12/2002
//----------------------------------------------------------------------------
void Instance::registerJacLIDs(const std::vector< std::vector<int> >& jacLIDVec)
{
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
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
   addStoreNode(symbol_table, li_store_output_state, getName().getEncodedName() + "_STATE");
}

//----------------------------------------------------------------------------
// Function       : Instance::getTVVec
// Purpose        :
// Special Notes  : WARNING -- Calling this function has the side effect of 
//                  automatically clearing the time-voltage pair vector.
//                  Use the function getAndDontClearTVVEC() if one does not
//                  need the history cleared.
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/28/2003
//----------------------------------------------------------------------------
void Instance::getTVVEC(std::vector< std::pair<double, double> > & TVVEC_Out)
{
  TVVEC_Out.clear();

  TVVEC_Out.insert(TVVEC_Out.end(), TVVEC.begin(), TVVEC.end());

  // Now that the digital simulator knows about these, let's forget them
  // and not worry about maintaining the list anymore
  TVVEC.clear();
}



//----------------------------------------------------------------------------
// Function       : Instance::getAndDontClearTVVEC
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Rich Schiek
// Creation Date  : 11/15/2021
//----------------------------------------------------------------------------
void Instance::getAndDontClearTVVEC(std::vector< std::pair<double, double> > & TVVEC_Out)
{
  TVVEC_Out.clear();

  TVVEC_Out.insert(TVVEC_Out.end(), TVVEC.begin(), TVVEC.end());
}

//----------------------------------------------------------------------------
// Function      : Instance::trimTVVEC
// Purpose       : clear out old time-voltage pairs
// Special Notes : ASSUMES the vector of t-v pairs is sorted by time!
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/10/2004
//----------------------------------------------------------------------------
void Instance::trimTVVEC(double earliestTime)
{
  std::vector< std::pair<double,double> >::iterator itVec;

  // get reference pointing to first element that exceeds earliestTime
  itVec = lower_bound(TVVEC.begin(),TVVEC.end(),std::pair<double,double>(earliestTime,0.0));
  // delete everything prior to that.
  TVVEC.erase(TVVEC.begin(),itVec);
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for this device which
//                 effectively is a single instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * solVector = extData.nextSolVectorRawPtr;
  double * daeFVec = extData.daeFVectorRawPtr;

  i0 = (solVector[li_Pos]-solVector[li_Neg])*G;

  // Load DAE F-vector
  daeFVec[li_Pos] += i0;
  daeFVec[li_Neg] += -i0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for this device which
//                 effectively is a single instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider <<std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
    Xyce::dout() << "  G = " << G << std::endl;
  }

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset] += G;

  (*dFdxMatPtr)[li_Pos][APosEquNegNodeOffset] -= G;

  (*dFdxMatPtr)[li_Neg][ANegEquPosNodeOffset] -= G;

  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += G;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/04
//-----------------------------------------------------------------------------

bool Instance::getInstanceBreakPoints ( std::vector<Util::BreakPoint> & breakPointTimes)
{
  bool bsuccess = true;

  // getBreakPoints is called after a timesteps solution is converged and
  // the step is finished.  So we can count on our voltage being the real one
  // that is valid for this time.

  // In this routine we need to check our current voltage value and see if it
  // has changed enough since our last time step to change a bit in the digital
  // output.  If so, we need to set a time/voltage pair of the current value
  // and (the current time + this device's conversion time), then set a
  // breakpoint of type PAUSE_BREAKPOINT for the that time

  double vPos(0.0);
  double vNeg(0.0);
  double deltaV(0.0), vFrac(0.0);
  double currentTime = getSolverState().currTime_;
  int newState;

  // Get the pointer to the vector of accepted solution values
  double * solVector = extData.nextSolVectorRawPtr;
  
  // reference to the stoVector.  Needed if we update the state
  Linear::Vector & stoVector = *(extData.nextStoVectorPtr);

  vPos = solVector[li_Pos];
  vNeg = solVector[li_Neg];

  deltaV = vPos-vNeg;
  //Convert the deltaV value into the integer state value.
  newState = deltaVToStateVal(deltaV);
  
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "In Instance::getInstanceBreakPoints.  deltaV = " << deltaV
	 << " vFrac = " << deltaVTovFrac(deltaV) << " nQuantLevels = " << nQuantLevels_ << std::endl
         << "  previous output level and new output state are " << lastOutputLevel_ 
         << " and " << newState << std::endl;     
  }

  if (newState != lastOutputLevel_)
  {
    // get time rounded to nearest femptosecond
    long long int timeInFS =static_cast<long long int>(
      (currentTime+model_.settlingTime_+6e-16)/1e-15);

    // since time is always advancing, we can simply push_back the pair,
    // and always be sure
    TVVEC.push_back(std::pair<double,double>(timeInFS*1e-15,deltaV));
    lastOutputLevel_ = newState;
    
    // update value in store vector so that it is available to the user for output
    stoVector[li_store_output_state]=(double) lastOutputLevel_;

    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "ADC ----------------------------------------" << std::endl;
      Xyce::dout() << "ADC Debug output.  name = " << getName() << std::endl;
      Xyce::dout() << "ADC setting pause breakpoint for " << currentTime << std::endl;
      double approxTime = timeInFS*1e-15;
      Xyce::dout() << "ADC Approximated time, to nearest femptosecond: " << approxTime << std::endl;
      Xyce::dout() << "ADC Time value pairs: " << std::endl;

      std::vector< std::pair<double, double> >::iterator beg = TVVEC.begin();
      std::vector< std::pair<double, double> >::iterator end = TVVEC.end();
      std::vector< std::pair<double, double> >::iterator itTV = beg;

      for (; itTV!=end; ++itTV)
      {
        std::pair<double, double>  & tvPair = *itTV;
        Xyce::dout() << "ADC time: " << tvPair.first << "  value: " << tvPair.second << std::endl;
      }

      Xyce::dout() << "ADC ----------------------------------------" << std::endl;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::deltaVToStateVal
// Purpose       : Convert the deltaV value into the integer state value.
// Special Notes : 
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 11/13/18
//-----------------------------------------------------------------------------
int Instance::deltaVToStateVal(double deltaV)
{
  int stateVal;
  double vFrac = deltaVTovFrac(deltaV); 

  if (vFrac < (1.0)/(nQuantLevels_) )
  {
    stateVal = 0;
  }
  else if (vFrac >= (nQuantLevels_-1.0)/(nQuantLevels_))
  {
    stateVal = nQuantLevels_ -1;
  }
  else
  {
    stateVal =  int(vFrac*nQuantLevels_);
  }

  return stateVal;
}

//-----------------------------------------------------------------------------
// Function      : Instance::deltaVTovFrac
// Purpose       : Convert the deltaV value into a fractional voltage value
// Special Notes : This function could be rolled into Instance::deltaVToStateVal(),
//                 except that the vFrac value is also used in debug statements in
//                 Instance::getInstanceBreakPoints().
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 11/13/18
//-----------------------------------------------------------------------------
double Instance::deltaVTovFrac(double deltaV)
{
  // This "upper voltage limit" and "lower voltage limit" approach is not right,
  // and needs to be replaced with a "Vref+"  node against which
  // Vin is compared, with a common negative reference (e.g. ground)
  // For now, let's always just document this failing, and tell the users
  // to wire vNeg to ground, and use 0.0 as the
  // lower limit.  Worry about doing it right later.
  double vFrac = deltaV/(model_.upperVoltageLimit_
                  - model_.lowerVoltageLimit_);
  return vFrac;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/05/08
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  double vPos(0.0), vNeg(0.0), deltaV(0.0), vFrac(0.0);
  double currentTime = getSolverState().currTime_;

  if (getSolverState().dcopFlag)
  {
    double * solVector = extData.nextSolVectorRawPtr;

    vPos = solVector[li_Pos];
    vNeg = solVector[li_Neg];
    deltaV = vPos-vNeg;

    TVVEC.push_back(std::pair<double,double>(0.0,deltaV));
  }

}

//----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Model::processParams()
{
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

//----------------------------------------------------------------------------
// Function       : Instance::getInstanceParamsMap
// Purpose        : Return a map of name,parameter for this instance
// Special Notes  : used by API for mixed signal
//                  in zero order version, these parameters happen to be
//                  model parameters, but that is not necessarily where they'll
//                  stay (the voltage limits, for example, should actually be
//                  taken from the nodal voltages of the reference nodes of the
//                  device instance)
// Scope          : public
// Creator        : Tom Russo, SNL, Component Information and Models
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Instance::getInstanceParamsMap(std::map<std::string,double>& paramsMap)
{
  paramsMap.clear();

  paramsMap["lowerVoltageLimit"] = model_.lowerVoltageLimit_;
  paramsMap["upperVoltageLimit"] = model_.upperVoltageLimit_;
  paramsMap["settlingTime"] = model_.settlingTime_;

  return true;
}

//----------------------------------------------------------------------------
// Function       : Instance::setBitVectorWidth
// Purpose        : set the number of bits and quantization levels in the 
//                  output of this ADC
// Special Notes  :
// Scope          : public
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool Instance::setBitVectorWidth(int width)
{
  outputBitVectorWidth_ = width;
  setNumberQuantLevels();

  return true;
}

//----------------------------------------------------------------------------
// Function       : Instance::setNumberQuantLevels
// Purpose        : set the number of quantization levels in the output of this ADC.
// Special Notes  : assumes that outputBitVectorWidth_ has been set to the desired
//                  value before this function is invoked.
// Scope          : public
// Creator        : Pete Sholander
// Creation Date  : 10/08/2018
//----------------------------------------------------------------------------
bool Instance::setNumberQuantLevels()
{
  nQuantLevels_ = 1;
  for (int i=0; i<outputBitVectorWidth_;++i)
    nQuantLevels_ *= 2;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::setNumberQuantLevels name = "
                 << getName() << " width = " << outputBitVectorWidth_
                 << " nQuantLevels_ = " << nQuantLevels_ << std::endl;
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : Model::Model
// Purpose        : constructor
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock&     model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
    lowerVoltageLimit_(0.0),
    upperVoltageLimit_(5.0),
    settlingTime_(1e-8),
    outputBitVectorWidth_(1)
{
  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (model_block.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//----------------------------------------------------------------------------
// Function       : Model::Model
// Purpose        : destructor
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Model::~Model()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//----------------------------------------------------------------------------
// Function       : printOutInstances
// Purpose        : debugging tool
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name\t\tmodelName\tParameters" << std::endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ {
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


// ADC Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);

    di.i0 = (solVec[di.li_Pos]-solVec[di.li_Neg])*di.G;

    fVec[di.li_Pos] += di.i0;
    fVec[di.li_Neg] += -di.i0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);

    dFdx[di.li_Pos][di.APosEquPosNodeOffset] += di.G;

    dFdx[di.li_Pos][di.APosEquNegNodeOffset] -= di.G;

    dFdx[di.li_Neg][di.ANegEquPosNodeOffset] -= di.G;

    dFdx[di.li_Neg][di.ANegEquNegNodeOffset] += di.G;
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : ::getBreakPoints
// Purpose       : getBreakPoints for all instances
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::getBreakPoints ( std::vector<Util::BreakPoint> & breakPointTimes )
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);
    bool tmpBool = di.getInstanceBreakPoints(breakPointTimes);
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;

}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}


void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("ADC")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("adc", 1)
      .registerModelType("adc", 1);
  }
}

} // namespace ADC
} // namespace Device
} // namespace Xyce

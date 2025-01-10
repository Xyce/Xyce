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
// Purpose        : This file implements the DAC digital to analog conversion
//                  device used in the integration of Xyce withe SAVANT VHDL
//                  simulator.
//
// Special Notes  :
//
// Creator        : Lon Waters
//
// Creation Date  : 07/26/2002
//
//----------------------------------------------------------------------------

#include <Xyce_config.h>

// ----------  Standard Includes ----------
#include <algorithm>

// ----------   Xyce Includes   ----------
#include <N_DEV_DAC.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {


namespace DAC {


void Traits::loadInstanceParameters(ParametricData<DAC::Instance> &p)
{
}

void Traits::loadModelParameters(ParametricData<DAC::Model> &p)
{
  // Set up double precision variables:
  p.addPar ("TR", 1.e-9, &DAC::Model::riseTime)
    .setUnit(U_SECOND)
    .setDescription("Rise Time");

  p.addPar ("TF", 1.e-9, &DAC::Model::fallTime)
    .setUnit(U_SECOND)
    .setDescription("Fall Time");

  p.addPar ("R",   0.01, &DAC::Model::R)
    .setUnit(U_OHM)
    .setDescription("Resistance");

  p.addPar ("L",  1.e-5, &DAC::Model::L)
    .setUnit(U_HENRY)
    .setDescription("Inductance");

  p.addPar ("C",    0.0, &DAC::Model::C)
    .setUnit(U_FARAD)
    .setDescription("Capacitance");
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
  const InstanceBlock & IB,
  Model & DACiter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(DACiter),
    numTVpairs_(0),
    v_pos(0),
    v_neg(0),
    i_bra(0),
    vDrop(0),
    voltage_(0),
    loc_(0),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 2;
    jacStamp[1].resize(1);
    jacStamp[1][0] = 2;
    jacStamp[2].resize(2);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

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
    Xyce::dout() << "  DACInstance::registerLIDs" << std::endl;
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

  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;

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
// Function      : Instance::loadNodeSymbols
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
void Instance::loadNodeSymbols(Util::SymbolTable &symbol_table) const
{
  addInternalNode(symbol_table, li_Bra, getName(), "branch");
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
  APosEquBraVarOffset = jacLIDVec[0][0];
  ANegEquBraVarOffset = jacLIDVec[1][0];
  ABraEquPosNodeOffset = jacLIDVec[2][0];
  ABraEquNegNodeOffset = jacLIDVec[2][1];
}

//----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;
  double * solVector = extData.nextSolVectorRawPtr;

  // Get the value for the source.
  updateVoltage(getSolverState().acceptedTime_);

  // get the value for v_pos, v_neg, i_bra
  v_pos = solVector[li_Pos];
  v_neg = solVector[li_Neg];
  i_bra  = solVector[li_Bra];
  vDrop = (v_pos-v_neg-voltage_);

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;
  bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : Instance::updateTVVEC
// Purpose        : Append the contents of newPairs to TVVEC
// Special Notes  :
// Scope          : public
// Creator        : Lisa Maynes & Lon Waters
// Creation Date  : 06/10/2003
//----------------------------------------------------------------------------
bool Instance::updateTVVEC (
  std::vector< std::pair<double, double> > const & newPairsIn )
{
  int i, last, newStart;
  double transitionTime;
  bool bsuccess = true;
  std::vector< std::pair<double,double> >::iterator itVec, itVec_end;
  std::vector< std::pair<double, double> > newPairs;
  std::vector< std::pair<double,double> >::const_iterator const_itVec;
  std::map<double, double> tv;
  std::map<double, double>::iterator tv_i, tv_end;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "In device " << getName() << std::endl
                 << "At time = " <<  getSolverState().acceptedTime_
                 << ", at beginning of Instance::updateTVVEC():\n"
                 << "   TVVEC size = " << numTVpairs_ << "\n"
                 << "   TVVEC loc  = " << loc_ << "\n"
                 << "   TVVEC contents:\n" ;
    itVec = TVVEC.begin();
    itVec_end = TVVEC.end();
    for( ; itVec != itVec_end; ++itVec )
    {
      Xyce::dout() << "   " << (*itVec).first
           << "s, " << (*itVec).second
           << "V\n";
    }
    Xyce::dout() << newPairsIn.size() << " New pairs:" << std::endl;
    for ( const_itVec = newPairsIn.begin() ; const_itVec != newPairsIn.end() ; ++const_itVec )
      Xyce::dout() << (*const_itVec).first << "  " << (*const_itVec).second << std::endl;
  }

  updateVoltage(getSolverState().acceptedTime_);

  itVec = TVVEC.begin();
  itVec_end = TVVEC.end();

  if (!TVVEC.empty())
  {
    for( ; itVec != itVec_end; ++itVec )
      tv[(*itVec).first] = (*itVec).second;
  }
  else if (!newPairsIn.empty())
  {
    // This bootstraps the processing if TVVEC is empty (e.g., on the first call to this function)
    for ( const_itVec = newPairsIn.begin() ; const_itVec != newPairsIn.end() ; ++const_itVec )
    {
      if ((*const_itVec).first >= getSolverState().acceptedTime_)
      {
        tv[(*const_itVec).first] = (*const_itVec).second;
        break;
      }
    }
  }

  if (!newPairsIn.empty())
  {
    if (getSolverState().acceptedTime_ == 0)
    {
      TVVEC.resize(0);
      if (newPairsIn.size() > 0)
        TVVEC.push_back(std::pair<double,double>(0,(*(newPairsIn.end()-1)).second));
      if (newPairsIn.size() > 1 && (*(newPairsIn.end()-1)).first > TVVEC[0].first)
        TVVEC.push_back(*(newPairsIn.end()-1));
      numTVpairs_ = TVVEC.size();
      updateVoltage(getSolverState().acceptedTime_);
      return bsuccess;
    }
    std::vector< std::pair<double, double> >::const_iterator n_i, n_end;
    std::vector< std::pair<double, double> >::iterator t_i, t_end;
    n_i = newPairsIn.begin();
    n_end = newPairsIn.end();
    for ( ; n_i != n_end ; ++n_i)
    {
      if ((*n_i).first < 0)
      {
        double d = -(*n_i).first;
        tv_i = lower_bound(tv.begin(), tv.end(), std::pair<const double, double>(d,0));
        tv.erase(tv_i,tv.end());
      }
    }
    n_i = newPairsIn.begin();
    for ( ; n_i != n_end ; ++n_i)
    {
      if ((*n_i).first >= getSolverState().acceptedTime_)
      {
        transitionTime = model_.riseTime;
        if (transitionTime > 0)
        {
          tv_i = lower_bound(tv.begin(), tv.end(), std::pair<const double, double>((*n_i).first,0));
          if (tv_i != tv.begin())
          {
            --tv_i;
            if ((*n_i).second < (*tv_i).second)
              transitionTime = model_.fallTime;
          }
          tv[(*n_i).first] = (*tv_i).second;
        }
        tv[(*n_i).first + transitionTime] = (*n_i).second;
      }
    }
    tv[getSolverState().acceptedTime_] = voltage_;
  }

  tv_i = lower_bound(tv.begin(), tv.end(), std::pair<const double, double>(getSolverState().acceptedTime_,0));
  if (tv_i != tv.begin())
    --tv_i;
  tv.erase(tv.begin(), tv_i);
  double lastTimeEntry = -1;
  TVVEC.clear();
  tv_i = tv.begin();
  tv_end = tv.end();
  for ( ; tv_i != tv_end ; ++tv_i)
  {
    if ((*tv_i).first - lastTimeEntry > 1e-15)
    {
      TVVEC.push_back(std::pair<double, double>((*tv_i).first,(*tv_i).second));
      lastTimeEntry = (*tv_i).first;
    }
  }
  numTVpairs_ = TVVEC.size();

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updateTVVEC():\n"
         << "   TVVEC size = " << numTVpairs_ << "\n"
         << "   TVVEC loc  = " << loc_ << "\n"
         << "   TVVEC contents:\n" ;
    std::vector< std::pair<double, double> >::iterator tv_i = TVVEC.begin();
    for( ; tv_i != TVVEC.end(); ++tv_i )
    {
      Xyce::dout() << "   " << (*tv_i).first
           << "s, " << (*tv_i).second
           << "V\n";
    }
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : Instance::updateVoltage
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 11/12/2002
//----------------------------------------------------------------------------
bool Instance::updateVoltage(double time)
{
  bool bsuccess = true;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << "  DACInstance::updateVoltage\n";
    Xyce::dout() << Xyce::section_divider << std::endl;
    Xyce::dout() << "  Time = " << time << std::endl;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag &&  numTVpairs_ > 0)
  {
    Xyce::dout() << "    TVVEC[numTVpairs_-1].first = "
                 << TVVEC[numTVpairs_-1].first << std::endl;
    for (int i=0 ; i<numTVpairs_ ; ++i)
    {
      Xyce::dout() << TVVEC[i].first << " :: " << TVVEC[i].second << std::endl;
    }
  }

  if( numTVpairs_ > 0 && time >= TVVEC[0].first )
  {
    if( time < TVVEC[numTVpairs_-1].first )
    {
      for( int i = 0; i < numTVpairs_ - 1; ++i )
      {
        if( time >= TVVEC[i].first && time <= TVVEC[i+1].first)
        {
          double time1 = TVVEC[i].first;
          double voltage1 = TVVEC[i].second;

          double time2 = TVVEC[i+1].first;
          double voltage2 = TVVEC[i+1].second;

          voltage_ = voltage1 + (voltage2 - voltage1) * (time - time1) / (time2 - time1);
          break;
        }
      }
    }
    else
    {
      voltage_ = TVVEC[numTVpairs_-1].second;
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "  voltage_ = " << voltage_ << std::endl;
    Xyce::dout() << Xyce::section_divider << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for this device.
//
// Special Notes : This is an algebraic constaint, and as such the resistor
//                 does make a contribution to it.
//
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Sim.
// Creation Date : 02/19/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double * daeFVec = extData.daeFVectorRawPtr;


  daeFVec[li_Pos] += i_bra;

  daeFVec[li_Neg] += -i_bra;

  daeFVec[li_Bra] += vDrop;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for this device.
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


  (*dFdxMatPtr)[li_Pos][APosEquBraVarOffset] += 1.0;

  (*dFdxMatPtr)[li_Neg][ANegEquBraVarOffset] -= 1.0;

  (*dFdxMatPtr)[li_Bra][ABraEquPosNodeOffset] += 1.0;

  (*dFdxMatPtr)[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 2/16/04
//-----------------------------------------------------------------------------

bool Instance::getInstanceBreakPoints ( std::vector<Util::BreakPoint> & breakPointTimes)
{
  bool bsuccess = true;
  double currentTime = getSolverState().currTime_;
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "In ::getInstanceBreakPoints " << std::endl;
    Xyce::dout() << " I want breakpoints.  Current time is " << currentTime << std::endl;
  }

  for (int i = 0; i < numTVpairs_ ; ++i)
  {
    // ERK:  the 1e-15 is a tolerance.  Fix for bug 1766.  Possibly use bpTol instead?
    // DNS:  No, bpTol tends to be ridiculously small and this might cause
    //       excessively small time steps.  I would agree if bpTol had a floor value,
    //       such as 1e-15.  Steps smaller than this are unjustified and numerically
    //       problematic.  For reference, light travels 0.3 micron in 1e-15 seconds.
    if (TVVEC[i].first >= currentTime - 1e-15 && model_.riseTime != 0 && model_.fallTime != 0)
    {
      breakPointTimes.push_back(Util::BreakPoint(TVVEC[i].first, Util::BreakPoint::SIMPLE));
    }
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "DAC ----------------------------------------" << std::endl;
    Xyce::dout() << "DAC getInstanceBreakPoints " << std::endl;
    Xyce::dout() << "DAC Debug output.  name = " << getName() << std::endl;
    Xyce::dout() << "DAC setting breakpoints at currentTime = " << currentTime << std::endl;
    Xyce::dout() << "DAC breakpoints: " << std::endl;

    std::vector< Util::BreakPoint  >::iterator beg = breakPointTimes.begin();
    std::vector< Util::BreakPoint  >::iterator end = breakPointTimes.end();
    std::vector< Util::BreakPoint  >::iterator itBP = beg;
    for (;itBP!=end;++itBP)
    {
      Util::BreakPoint & bp = *itBP;
      Xyce::dout() << "DAC breakpoint: " << bp.value() << std::endl;
    }

    Xyce::dout() << "DAC ----------------------------------------" << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInternalState
// Purpose       : Generates an DeviceState object and populates
//                 it with the contents of the TVVEC vector for use by
//                 restarts
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 09/03/04
//-----------------------------------------------------------------------------
DeviceState * Instance::getInternalState()
{
  int vsize,i,j;
  // allocate object to return
  DeviceState * myState = new DeviceState;

  myState->ID=getName().getEncodedName();
  vsize=TVVEC.size();
  // pack the pairs into the single vector of doubles.
  myState->data.resize(vsize*2);
  for (i=0;i<vsize;++i)
  {
    j=i*2;
    myState->data[j]=TVVEC[i].first;
    myState->data[j+1]=TVVEC[i].second;
  }

  return myState;
}
//-----------------------------------------------------------------------------
// Function      : Instance::setInternalState
// Purpose       : Reload TVVEC data from restart
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 09/03/04
//-----------------------------------------------------------------------------
bool Instance::setInternalState(const DeviceState &state)
{
  int dsize=state.data.size();
  int vsize,i,j;
  if (getName().getEncodedName() != state.ID)
  {
    Report::DevelFatal().in("DAC::Instance::setInternal") << "ID(" << state.ID << ") from restart does not match my name (" << getName() << ")";
    return false;
  }

  if (dsize%2 != 0)
  {
    UserError(*this) << "Data size from restart (" << dsize << " not a multiple of 2!";
    return false;
  }

  vsize=dsize/2;
  TVVEC.clear();
  TVVEC.resize(vsize);
  for (i=0;i<vsize;++i)
  {
    j=i*2;
    TVVEC[i].first=state.data[j];
    TVVEC[i].second=state.data[j+1];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 05/12/09
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
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
// Function       : Model::Model
// Purpose        : constructor
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/29/2002
//----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock&     MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    riseTime(1.0e-9),
    fallTime(1.0e-9),
    R(.01),
    L(0.0),
    C(0.0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

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



// DAC Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 02/25/2009
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);

    // here we update the voltage
    di.updateVoltage(getSolverState().currTime_);

    // get the value for v_pos, v_neg, i_bra
    di.v_pos = solVec[di.li_Pos];
    di.v_neg = solVec[di.li_Neg];
    di.i_bra  = solVec[di.li_Bra];
    di.vDrop = (di.v_pos-di.v_neg-di.voltage_);
    if (DEBUG_DEVICE)
    {
      Xyce::dout() << "DAC ----------------------------------------" << std::endl;
      Xyce::dout() << "DAC Debug output.  name = " << di.getName() << std::endl;
      Xyce::dout() << "DAC v_pos = " << di.v_pos <<std::endl;
      Xyce::dout() << "DAC v_neg = " << di.v_neg <<std::endl;
      Xyce::dout() << "DAC i_bra = " << di.i_bra <<std::endl;
      Xyce::dout() << "DAC vDrop = " << di.vDrop <<std::endl;
      Xyce::dout() << "DAC voltage_ = " << di.voltage_ <<std::endl;
      Xyce::dout() << "DAC ----------------------------------------" << std::endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, double * leadF, double * leadQ, double * junctionV)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);

    fVec[di.li_Pos] += di.i_bra;

    fVec[di.li_Neg] += -di.i_bra;

    fVec[di.li_Bra] += di.vDrop;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & di = *(*it);

    dFdx[di.li_Pos][di.APosEquBraVarOffset] += 1.0;

    dFdx[di.li_Neg][di.ANegEquBraVarOffset] -= 1.0;

    dFdx[di.li_Bra][di.ABraEquPosNodeOffset] += 1.0;

    dFdx[di.li_Bra][di.ABraEquNegNodeOffset] -= 1.0;
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("DAC")!=deviceMap.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("dac", 1)
      .registerModelType("dac", 1);
  }
}

} // namespace DAC
} // namespace Device
} // namespace Xyce

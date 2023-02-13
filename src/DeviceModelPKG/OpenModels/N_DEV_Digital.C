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
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 01/05/06
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Math.h>
#include <iostream>
#include <sstream>

#include <N_DEV_Const.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Digital.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Message.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {
namespace Digital {

// enables debug output for just the digital devices
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<Digital::Instance> &p)
{
// Set up double precision variables:

  // Set up non-double precision variables:
  p.addPar ("IC1", false, &Digital::Instance::ic1)
    .setUnit(U_LOGIC)
    .setDescription("Vector of initial values for output(s)");
  p.addPar ("IC2", false, &Digital::Instance::ic2);

  p.makeVector ("IC", 2);
}

void Traits::loadModelParameters(ParametricData<Digital::Model> &p)
{
  // Set up double precision variables:
  p.addPar ("VLO", 0., &Digital::Model::vlo)
    .setUnit(U_VOLT)
    .setDescription("Internal low state supply voltage");

  p.addPar ("VHI", 0., &Digital::Model::vhi)
    .setUnit(U_VOLT)
    .setDescription("Internal high state supply voltage");

  p.addPar ("VREF", 0., &Digital::Model::vref)
    .setUnit(U_VOLT)
    .setDescription("Internal reference voltage for inputs");

  p.addPar ("CLO", 1.e-6, &Digital::Model::clo)
    .setUnit(U_FARAD)
    .setDescription("Capacitance between output node and low reference");

  p.addPar ("CHI", 1.e-6, &Digital::Model::chi)
    .setUnit(U_FARAD)
    .setDescription("Capacitance between output node and high reference");

  p.addPar ("CLOAD", 1.e-6, &Digital::Model::cload)
    .setUnit(U_FARAD)
    .setDescription("Capacitance between input node and input reference");

  p.addPar ("RLOAD", 1000., &Digital::Model::rload)
    .setUnit(U_OHM)
    .setDescription("Resistance between input node and input reference");

  p.addPar ("S0RLO", 100., &Digital::Model::s0rlo)
    .setUnit(U_OHM)
    .setDescription("Low state resistance between output node and low reference");

  p.addPar ("S0RHI", 100., &Digital::Model::s0rhi)
    .setUnit(U_OHM)
    .setDescription("Low state resitance between output node and high reference");

  p.addPar ("S0TSW", 1.e-8, &Digital::Model::s0tsw)
    .setUnit(U_SECOND)
    .setDescription("Switching time transition to low state");

  p.addPar ("S0VLO", -1.5, &Digital::Model::s0vlo)
    .setUnit(U_VOLT)
    .setDescription("Minimum voltage to switch to low state");

  p.addPar ("S0VHI", 1.7, &Digital::Model::s0vhi)
    .setUnit(U_VOLT)
    .setDescription("Maximum voltage to switch to low state");

  p.addPar ("S1RLO", 100., &Digital::Model::s1rlo)
    .setUnit(U_OHM)
    .setDescription("High state resistance between output node and low reference");

  p.addPar ("S1RHI", 100., &Digital::Model::s1rhi)
    .setUnit(U_OHM)
    .setDescription("High state resistance between output node and high reference");

  p.addPar ("S1TSW", 1.e-8, &Digital::Model::s1tsw)
    .setUnit(U_SECOND)
    .setDescription("Switching time transition to high state");

  p.addPar ("S1VLO", 0.9, &Digital::Model::s1vlo)
    .setUnit(U_VOLT)
    .setDescription("Minimum voltage to switch to high state");

  p.addPar ("S1VHI", 7.0, &Digital::Model::s1vhi)
    .setUnit(U_VOLT)
    .setDescription("Maximum voltage to switch to high state");

  p.addPar ("DELAY", 1.e-8, &Digital::Model::delay)
    .setUnit(U_SECOND)
    .setDescription("Delay time of device");
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{

  // If there are any time dependent parameters, set their values for
  // the current time.

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    li_Lo(-1),
    li_Hi(-1),
    li_Ref(-1),
    breakTime(0.),
    row_Lo(-1),
    row_Hi(-1),
    row_Ref(-1),
    supportsXState_(false),
    gateInfo_(0)
{
  // added for compatibility with PSpice DIGINITSTATE
  digInitState_ = getDeviceOptions().digInitState;
  if ( digInitState_ < 0 || digInitState_ > 3 )
  {
    // If the user specifies an invalid value then use the default value.
    digInitState_ = 3;
    UserWarning(*this) << "Invalid DIGINITSTATE: " << getDeviceOptions().digInitState 
		       << ". Setting to default DIGINITSTATE: " <<  digInitState_ << std::endl;
  }  
 

  // added for debug purposes for tracking input/output state changes
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    prevInputStateChangeTime_ = 0;
    inputStateChange_ = false;
  }

  int dev_numInputs = 0;
  char dev_letter = getDeviceLetter(getName());
  numExtVars = 0;

  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated.
  //dev_numInputs = getNumInputs(getName());
  //deviceInfo_->updatePowerPinLi(li_Lo,li_Hi,li_Ref,numExtVars, 
  //	                        model_.given("VLO"), 
  //                            model_.given("VHI"), 
  //                            model_.given("VREF"));
  if (dev_letter == 'U') {
    dev_numInputs = getNumInputs(getName());

    // For U devices, DPWR and DGND are always specified on the instance line.
    // Warning message if VHI, VLO or VREF are in the model card.
    li_Lo = 0;
    li_Hi = 0;
    li_Ref = 0;
    numExtVars += 2;
    if (model_.given("VLO"))
    {
      UserWarning(*this)<< "VLO model parameter ignored for U digital device";
    }
    if (model_.given("VHI"))
    {
      UserWarning(*this)<< "VHI model parameter ignored for U digital device";
    }
    if (model_.given("VREF"))
    {
      UserWarning(*this)<< "VREF model parameter ignored for U digital device";
    }
  }
  else if (dev_letter == 'Y')
  {
    // legacy code required to support VLO, VHI and VREF variables
    // being on the instance line rather than in model card in Y devices
    UserWarning(*this)<< "Y-type digital device is deprecated. Consider using U-type digital device instead.";
    if (!model_.given("VLO"))
    {
      li_Lo = 0;
      ++numExtVars;
    }
    if (!model_.given("VHI"))
    {
      li_Hi = 0;
      ++numExtVars;
    }
    if (!model_.given("VREF"))
    {
      li_Ref = 0;
      ++numExtVars;
      } 
  }
  else
  {
    UserError(*this) << "Digital device letter must be Y or U: " << getName();
    return;
  }

  // Configure number of inputs/outputs for each device
  // Y devices are limited to 2 inputs.
  std::string dev_type = getName().getDeviceType();

  // the need for dev_letter in this function call will go away 
  // once the Y devices are deprecated
  if (dev_type == "NOT" || dev_type == "INV")
  {
    gateInfo_ = new InvData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "AND")
  {
    gateInfo_ = new AndData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "NAND")
  {
    gateInfo_ = new NandData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "OR")
  {
    gateInfo_ = new OrData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "NOR")
  {
    gateInfo_ = new NorData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "ADD")
  {
    gateInfo_ = new AddData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "XOR")
  {
    gateInfo_ = new XorData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "NXOR")
  {
    gateInfo_ = new NxorData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "DFF")
  {
    gateInfo_ = new DffData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "JKFF")
  {
    gateInfo_ = new JkffData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "TFF")
  {
    gateInfo_ = new TffData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "DLTCH")
  {
    gateInfo_ = new DltchData(dev_type,dev_letter,dev_numInputs);
  }
  else if (dev_type == "BUF")
  {
    gateInfo_ = new BufData(dev_type,dev_letter,dev_numInputs);
  }
  else
  {
    UserError(*this) << "Unknown digital device type " << dev_type;
    return;
  }
  // return the number of I/O pins to the Instance
  numInput = gateInfo_->getNumInput();
  numOutput = gateInfo_->getNumOutput();
  gate = gateInfo_->getType();
  supportsXState_ = gateInfo_->getSupportsXState();
  icGiven_.resize(numOutput);

  int iBase = numExtVars;
  int oBase = numExtVars + numInput;
  numExtVars += numInput + numOutput;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << "Digital Device " << getName() << " has iBase = " <<
      iBase << ", oBase = " << oBase << ", numExtVars = " <<
      numExtVars << std::endl;
  }

  // catch cases of AND, NAND, OR or NOR gate with only one input specified
  // or the number of nodes on the instance line does not match the
  // (N) value specified as part of the gate type (e.g, AND(4)).
  // Also catch case where gates with a fixed number of inputs have
  // wrong number of inputs
  gateInfo_->checkErrors(*this,instance_block,iBase,dev_numInputs);

  numIntVars   = 0;
  numStateVars = 4*numInput + 6*numOutput;

  li_Inp.resize(numInput);
  li_currentStateInp.resize(numInput);
  li_transitionTimeInp.resize(numInput);
  li_QinpState.resize(numInput);
  li_IinpState.resize(numInput);
  li_Out.resize(numOutput);
  li_currentStateOut.resize(numOutput);
  li_transitionTimeOut.resize(numOutput);
  li_QloState.resize(numOutput);
  li_IloState.resize(numOutput);
  li_QhiState.resize(numOutput);
  li_IhiState.resize(numOutput);

  qlo.resize(numOutput);
  ilo.resize(numOutput);
  vcaplo.resize(numOutput);
  qhi.resize(numOutput);
  ihi.resize(numOutput);
  vcaphi.resize(numOutput);
  qref.resize(numInput);
  iref.resize(numInput);
  vcapref.resize(numInput);
  rilo.resize(numOutput);
  rihi.resize(numOutput);
  riref.resize(numInput);
  currentOut.resize(numOutput);
  currentIn.resize(numInput);
  glo.resize(numOutput);
  ghi.resize(numOutput);

  qInp.resize(numInput);
  iInp.resize(numInput);
  vcapInp.resize(numInput);
  currentInp.resize(numInput);

  inpL.resize(numInput);
  iTime.resize(numInput);
  outL.resize(numOutput);
  oTime.resize(numOutput);

  // These are to store the indicies into the jacobian for the four element
  // stamps for the capacitor/resistors connected to the input/outputs. The
  // format is to have 6 int vectors for each stamp with format:
  // (row 1, col 1, col 2, row 2, col 1, col 2).  These will be replaced in
  // registerJacLIDs with the indicies into the actual jacobian matrix.

  li_jac_Ref.resize(numInput);
  li_jac_Hi.resize(numOutput);
  li_jac_Lo.resize(numOutput);

  devConMap.resize(numExtVars);
  for (int i=0 ; i<numExtVars ; ++i)
    devConMap[i] = 1;

  jacStamp.resize(numExtVars);
  int row = 0;
  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated.
  if (dev_letter == 'U')
  {
    // digital power and digital ground node are always on
    // the instance line for a U device.  The low reference
    // voltage for inputs is assumed to be the same as the
    // digital ground node.
    row_Hi = row;
    jacStamp[row].push_back(row);
    for (int i=0 ; i<numOutput ; ++i)
    {
      li_jac_Hi[i].push_back(row);
      li_jac_Hi[i].push_back(0);
      li_jac_Hi[i].push_back(jacStamp[row].size());
      li_jac_Hi[i].push_back(oBase+i);
      li_jac_Hi[i].push_back(jacStamp[oBase+i].size());
      jacStamp[row].push_back(oBase+i);
      jacStamp[oBase+i].push_back(row);
    }
    ++row;

    row_Lo = row;
    jacStamp[row].push_back(row);
    for (int i=0 ; i<numOutput ; ++i)
    {
      li_jac_Lo[i].push_back(row);
      li_jac_Lo[i].push_back(0);
      li_jac_Lo[i].push_back(jacStamp[row].size());
      li_jac_Lo[i].push_back(oBase+i);
      li_jac_Lo[i].push_back(jacStamp[oBase+i].size());
      jacStamp[row].push_back(oBase+i);
      jacStamp[oBase+i].push_back(row);
    }

    row_Ref = row;
    for (int i=0 ; i<numInput ; ++i)
    {
      li_jac_Ref[i].push_back(row);
      li_jac_Ref[i].push_back(0);
      li_jac_Ref[i].push_back(1);
      li_jac_Ref[i].push_back(iBase+i);
      li_jac_Ref[i].push_back(jacStamp[iBase+i].size());
      jacStamp[row].push_back(iBase+i);
      jacStamp[iBase+i].push_back(row);
    }
  }
  else if (dev_letter == 'Y')
  {
    // output low and output high reference voltages and input
    // low reference voltage are optionally on the instance line
    // for Y devices.
    if (li_Lo == 0)
    {
      row_Lo = row;
      jacStamp[row].push_back(row);
      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i].push_back(row);
        li_jac_Lo[i].push_back(0);
        li_jac_Lo[i].push_back(jacStamp[row].size());
        li_jac_Lo[i].push_back(oBase+i);
        li_jac_Lo[i].push_back(jacStamp[oBase+i].size());
        jacStamp[row].push_back(oBase+i);
        jacStamp[oBase+i].push_back(row);
      }
      ++row;
    }
    if (li_Hi == 0)
    {
      row_Hi = row;
      jacStamp[row].push_back(row);
      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i].push_back(row);
        li_jac_Hi[i].push_back(0);
        li_jac_Hi[i].push_back(jacStamp[row].size());
        li_jac_Hi[i].push_back(oBase+i);
        li_jac_Hi[i].push_back(jacStamp[oBase+i].size());
        jacStamp[row].push_back(oBase+i);
        jacStamp[oBase+i].push_back(row);
      }
      ++row;
    }
    if (li_Ref == 0)
    {
      row_Ref = row;
      jacStamp[row].push_back(row);
      for (int i=0 ; i<numInput ; ++i)
      {
        li_jac_Ref[i].push_back(row);
        li_jac_Ref[i].push_back(0);
        li_jac_Ref[i].push_back(jacStamp[row].size());
        li_jac_Ref[i].push_back(iBase+i);
        li_jac_Ref[i].push_back(jacStamp[iBase+i].size());
        jacStamp[row].push_back(iBase+i);
        jacStamp[iBase+i].push_back(row);
      }
      ++row;
    }
  }
  else
  {
    UserError(*this) << "Digital device letter must be Y or U: " << getName();
    return;
  }

  for (int i=0 ; i<numInput ; ++i)
  {
    if (li_Ref == 0)
    {
      li_jac_Ref[i].push_back(jacStamp[iBase+i].size());
    }
    else
    {
      li_jac_Ref[i].push_back(iBase+i);
    }
    jacStamp[iBase+i].push_back(iBase+i);
  }
  for (int i=0 ; i<numOutput ; ++i)
  {
    if (li_Lo == 0)
    {
      li_jac_Lo[i].push_back(jacStamp[oBase+i].size());
    }
    else
    {
      li_jac_Lo[i].push_back(oBase+i);
    }
    if (li_Hi == 0)
    {
      li_jac_Hi[i].push_back(jacStamp[oBase+i].size());
    }
    else
    {
      li_jac_Hi[i].push_back(oBase+i);
    }
    jacStamp[oBase+i].push_back(oBase+i);
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  //  Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  delete gateInfo_;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)

{

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      << "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs")
      << "numExt != numExtVars";
  }

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.  For the matrix indices, first do the
  // rows.

  int i=0, j;
  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated
  char dev_letter = getDeviceLetter(getName());
  if (dev_letter == 'U')
  {
    // ordering on U-device instance line is dig_power_node (DPWR) then
    // dig_ground_node (DGND).  Assume that input low-reference
    // voltage (VREF in Y devices) is equal to DGND.
    li_Hi = extLIDVec[i++];
    li_Lo = extLIDVec[i++];
    li_Ref = li_Lo;
  }
  else if (dev_letter == 'Y')
  {
    // for Y devices, ordering of output high/low reference nodes (VHI/VLO)
    // on instance line is reversed.  Input low-reference voltage (VREF),
    // VHI and VLO can either be on the instance line or in the model card.
    if (li_Lo == 0)
      li_Lo = extLIDVec[i++];
    if (li_Hi == 0)
      li_Hi = extLIDVec[i++];
    if (li_Ref == 0)
      li_Ref = extLIDVec[i++];
  }
  else
  {
    UserError(*this) << "Digital device letter must be Y or U: " << getName();
  }

  for (j=0 ; j<numInput ; ++j)
    li_Inp[j] = extLIDVec[i++];
  for (j=0 ; j<numOutput ; ++j)
    li_Out[j] = extLIDVec[i++];
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    DevelFatal(*this).in("Instance::registerStateLIDs")
      << "numSta != numStateVars";
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int j=0;

  for (int i=0 ; i<numInput ; ++i)
  {
    li_currentStateInp[i] = staLIDVec[j++];
    li_transitionTimeInp[i] = staLIDVec[j++];
    li_QinpState[i] = staLIDVec[j++];
    li_IinpState[i] = staLIDVec[j++];
  }

  for (int i=0 ; i<numOutput ; ++i)
  {
    li_currentStateOut[i] = staLIDVec[j++];
    li_transitionTimeOut[i] = staLIDVec[j++];
    li_QloState[i] = staLIDVec[j++];
    li_IloState[i] = staLIDVec[j++];
    li_QhiState[i] = staLIDVec[j++];
    li_IhiState[i] = staLIDVec[j++];
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
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  int lo_present, hi_present;

  lo_present = (row_Lo < 0) ? 0 : 1;
  hi_present = (row_Hi < 0) ? 0 : 1;

  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated.
  char dev_letter = getDeviceLetter(getName());

  if (dev_letter == 'U')
  {
    if ( (row_Lo == -1) || (row_Hi == -1) || (row_Ref == -1) )
    {
       UserError(*this) << "Internal error in Instance::registerJacLIDs() for " << getName();
    }
    else
    {
      for (int i=0 ; i<numInput ; ++i)
      {
        li_jac_Ref[i][1] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][1]];
        li_jac_Ref[i][2] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][2]];
        li_jac_Ref[i][4] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][4]];
        li_jac_Ref[i][5] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][5]];
      }

      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i][1] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][1]];
        li_jac_Lo[i][2] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][2]];
        li_jac_Lo[i][4] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][4]];
        li_jac_Lo[i][5] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][5]];
      }

      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i][1] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][1]];
        li_jac_Hi[i][2] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][2]];
        li_jac_Hi[i][4] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][4]];
        li_jac_Hi[i][5] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][5]];
      }
    }
  }
  else if (dev_letter == 'Y')
  {
    // row_Ref == -1 means that VREF is in the model card
    if (row_Ref == -1)
    {
      for (int i=0 ; i<numInput ; ++i)
      {
        li_jac_Ref[i].push_back(jacLIDVec[li_jac_Ref[i][0]][0]);
      }
    }
    else
    {
      // this for loop is identical for both U and Y devices
      for (int i=0 ; i<numInput ; ++i)
      {
        li_jac_Ref[i][1] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][1]];
        li_jac_Ref[i][2] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][2]];
        li_jac_Ref[i][4] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][4]];
        li_jac_Ref[i][5] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][5]];
      }
    }

    // row_Lo == -1 means that VLO is in the model card
    if (row_Lo == -1)
    {
      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i].push_back(jacLIDVec[li_jac_Lo[i][0]][hi_present]);
      }
    }
    else
    {
      // this for loop is identical for both U and Y devices
      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i][1] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][1]];
        li_jac_Lo[i][2] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][2]];
        li_jac_Lo[i][4] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][4]];
        li_jac_Lo[i][5] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][5]];
      }
    }

    // row_Hi == -1 means that VHI is in the model card
    if (row_Hi == -1)
    {
      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i].push_back(jacLIDVec[li_jac_Hi[i][0]][lo_present]);
      }
    }
    else
    {
      // this for loop is identical for both U and Y devices
      for (int i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i][1] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][1]];
        li_jac_Hi[i][2] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][2]];
        li_jac_Hi[i][4] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][4]];
        li_jac_Hi[i][5] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][5]];
      }
    }
  }
  else
  {
    UserError(*this) << "Digital device letter must be Y or U: " << getName();
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  //added for debug purposes, for tracking input state changes
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    inputStateChange_ = false;
  }

  double v_poslo, v_poshi, v_posref, v_neg;
  double elapsed, time, frac, lastT;
  int currentState;
  double transitionTime;
  bool changeState = false; //Genie 110812
  bool clocking = false; //Genie 111212

  Linear::Vector & solVector = *(extData.nextSolVectorPtr);
  Linear::Vector & staVector = *(extData.nextStaVectorPtr);
  Linear::Vector & oldStaVector = *(extData.currStaVectorPtr);
  Linear::Vector & oldSolVector = *(extData.currSolVectorPtr);
  
  // added so that the oldState vector can be passed into gateInfo_->evaluateTruthTable()
  std::vector<bool> oldState;
  oldState.resize(numOutput);

  // The convention in this device is to consider the 'positive' side of the
  // capacitors as the Vsrc or node that supplies the voltage, or for the
  // output, vref. The 'positive' voltages are thus the same for all
  // input/outputs
  v_poslo = (li_Lo >= 0) ? solVector[li_Lo] : model_.vlo;

  // If VHI is not specified, li_Hi >=0 and v_poshi is derived from
  // CLOAD/RLOAD by solVector.  If VHI is specified, li_Hi == -1 and
  // v_poshi is specified by the netlist model card.
  v_poshi = (li_Hi >= 0) ? solVector[li_Hi] : model_.vhi;
  v_posref = (li_Ref >= 0) ? solVector[li_Ref] : model_.vref;

  lastT = 0;

  for (int i=0 ; i<numInput ; ++i)
  {
    //initialize
    currentState = static_cast <int> (oldStaVector[li_currentStateInp[i]]);
    transitionTime = oldStaVector[li_transitionTimeInp[i]];
    changeState = false; //Genie 022013 Clear the memory of changeState.

    // obtain voltage drop accross the capacitor:
    v_neg = solVector[li_Inp[i]];

    time = getSolverState().currTime_;
    elapsed = time - transitionTime;

    vcapref[i] = v_posref-v_neg;
    riref[i] = model_.gload*(vcapref[i]);

    // Obtain the "current"  value for the charge stored in the capacitors.
    qref[i] = model_.cload*vcapref[i];

    if (getSolverState().dcopFlag)
    {
      transitionTime = 0;
      currentState = (-vcapref[i] < model_.s0vhi) ? 0 : 1;

      oldStaVector[li_currentStateInp[i]] = currentState;
      oldStaVector[li_transitionTimeInp[i]] = transitionTime;
    }

    iTime[i] = transitionTime;

    staVector[li_transitionTimeInp[i]] = transitionTime;

    if (currentState == 0)
    {
      inpL[i] = false;
      if ((-vcapref[i] > model_.s0vhi) && (-vcapref[i] > model_.s1vlo))
      {
          currentState = 1;
          changeState = true;
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
	    Xyce::dout() << "Device " << getName() << " changed state from 1 to 0 at time " << getSolverState().currTime_ << std::endl;
          }
          // This call updates the clock status if the gate has a clock pin
          // and this line is the clock pin. clocking is a local variable within
          // updatePrimaryState()
          if (gateInfo_->isClockLine(i)) { clocking = true; }
      }
    }
    else
    {
      inpL[i] = true;
      if ((-vcapref[i] < model_.s1vlo) && (-vcapref[i] < model_.s0vhi))
      {
          currentState = 0;
          changeState = true;
          if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
          {
	    Xyce::dout() << "Device " << getName() << " changed state from 1 to 0 at time " << getSolverState().currTime_ << std::endl;
          }
          // This call updates the clock status if the gate has a clock pin
          // and this line is the clock pin. clocking is a local variable within
          // updatePrimaryState()
          if (gateInfo_->isClockLine(i)) { clocking = true; }
      }
    }

    if (changeState)
    {
      double vOld, del;

      inpL[i] = (currentState == 1);
      vOld = (li_Ref >= 0) ? oldSolVector[li_Ref] : model_.vref;
      vOld = oldSolVector[li_Inp[i]] - vOld;
      if (fabs(vcapref[i]+vOld) < 1.e-12)
        del = 0;
      else
      {
        if (inpL[i])
        {
          del = getSolverState().currTimeStep_ * (-vcapref[i] - model_.s0vhi)/(-vcapref[i] - vOld);
        }
        else
        {
          del = getSolverState().currTimeStep_ * (vcapref[i] + model_.s1vlo)/(vcapref[i] + vOld);
        }
      }
      staVector[li_transitionTimeInp[i]] = time - del;
      iTime[i] = time;
    }

    staVector[li_currentStateInp[i]] = currentState;
    staVector[li_QinpState[i]] = qref[i];

    if (iTime[i] > lastT)
      lastT = iTime[i];

    // added for debug, for tracking input state changes
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
    {
      if (changeState == true && inputStateChange_ == false)
      {
        inputStateChange_ = true;
      }
    }
  } // end for loop numInput i

  // update the truth table for this iteration
  for (int i=0 ; i <numOutput ; ++i)
  {
    oldState[i] = oldStaVector[li_currentStateOut[i]];
  }
  gateInfo_->evalTruthTable(inpL, outL, oTime, lastT, model_.delay,
                              getSolverState().dcopFlag,clocking, oldState);

  bool curr;
  breakTime = 0;
  for (int i=0 ; i<numOutput ; ++i)
  {
    time = getSolverState().currTime_;
    if (getSolverState().dcopFlag)
    {
      gateInfo_->setIC(*this,i);
      
      oldStaVector[li_currentStateOut[i]] = outL[i]?1:0;
      oldStaVector[li_transitionTimeOut[i]] = time;
    }

    //current logic state of output nodes
    currentState = static_cast <int> (oldStaVector[li_currentStateOut[i]]);
    transitionTime = oldStaVector[li_transitionTimeOut[i]];

    if (currentState == 1)
      curr = true;
    else
      curr = false;

    if (curr != outL[i]) // This is executed when scopFlag is false
    {
      if (oTime[i] <= time)
      {
        currentState = 1-currentState;
        transitionTime = oTime[i];
      }
      else
      {
        if (breakTime == 0 || (breakTime > 0 && breakTime > oTime[i]))
        {
          breakTime = oTime[i];
        }
      }
    }

    staVector[li_currentStateOut[i]] = currentState;
    staVector[li_transitionTimeOut[i]] = transitionTime;

    // obtain voltage drop accross the capacitors:
    v_neg = solVector[li_Out[i]];

    elapsed = time - transitionTime;

    if (currentState == 0)
      elapsed /= model_.s0tsw;
    else if (currentState == 1)
      elapsed /= model_.s1tsw;
    else
    {
      DevelFatal(*this).in("Instance::updateSecondaryState")
        << "Unrecognized state";
    }

    frac = exp(-elapsed); //Genie 112812. This line can be omitted.
    if (transitionTime == 0)
      frac = 0;
    else
    {
      if (elapsed > 1)
        frac = 0;
      else
      {
        // This is a simple linear transition.  Since there is a
        // breakpoint at the start of the transition it is OK to
        // have a discontinuity there.
        frac = 1-elapsed;
      }
    }

    // DFF and DLTCH outputs are put into the "X" state at DCOP if the DIGINITSTATE option 
    // equals 2 and an initial condition was not specified for this output.
    if ( digInitState_ == 2 && getSolverState().dcopFlag && supportsXState_ && !icGiven_[i])
    {
      // for the "X" state have the output both pulled-up to the high-state voltage and pulled 
      // down to the low-state voltage
      glo[i] = 1/(frac*model_.s1rlo + (1-frac)*model_.s0rlo);
      ghi[i] = 1/(frac*model_.s0rhi + (1-frac)*model_.s1rhi);
    }
    else if (currentState == 0)
    {
      glo[i] = 1/(frac*model_.s1rlo + (1-frac)*model_.s0rlo);
      ghi[i] = 1/(frac*model_.s1rhi + (1-frac)*model_.s0rhi);
    }
    else
    {
      glo[i] = 1/(frac*model_.s0rlo + (1-frac)*model_.s1rlo);
      ghi[i] = 1/(frac*model_.s0rhi + (1-frac)*model_.s1rhi);
    }

    rilo[i] = glo[i]*(v_poslo-v_neg);
    rihi[i] = ghi[i]*(v_poshi-v_neg);

    vcaplo[i] = v_poslo-v_neg;
    vcaphi[i] = v_poshi-v_neg;

    // Obtain the "current"  value for the charge stored in the capacitors.
    qlo[i] = model_.clo*vcaplo[i];
    qhi[i] = model_.chi*vcaphi[i];

    staVector[li_QloState[i]] = qlo[i];
    staVector[li_QhiState[i]] = qhi[i];

  }// end for loop numOutput i

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{

  Linear::Vector * staVectorPtr = extData.nextStaVectorPtr;

  // Now that the state vector for time=0 is up-to-date, get the
  // derivative with respect to time of the charge, to obtain the
  // best estimate for the current in the capacitors.

  for (int i=0 ; i<numOutput ; ++i)
  {
    ilo[i] = (*extData.nextStaDerivVectorPtr)[li_QloState[i]];
    ihi[i] = (*extData.nextStaDerivVectorPtr)[li_QhiState[i]];

    currentOut[i] = ilo[i] + ihi[i] + rilo[i] + rihi[i];

    (*staVectorPtr)[li_IloState[i]] = ilo[i];
    (*staVectorPtr)[li_IhiState[i]] = ihi[i];
  }

  for (int i=0 ; i<numInput ; ++i)
  {
    iref[i] = (*extData.nextStaDerivVectorPtr)[li_QinpState[i]];

    currentIn[i] = iref[i] + riref[i];

    (*staVectorPtr)[li_IinpState[i]] = iref[i];
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : Add break point for anticipated digital output transition
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/07/06
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints(
  std::vector<Util::BreakPoint> &       breakPointTimes)
{
  if (breakTime > getSolverState().currTime_)
  {
    breakPointTimes.push_back(breakTime);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       : This function is to be called ONLY at the point
//                 when the time integrator has determined we've got a
//                 converged, acceptable solution and is accepting it,
//                 but before it's updated its times and rotated vectors.
//                 It is used for debugging purpose
//
// Special Notes : In SPICE this same stuff was done in the "TRAaccept" function.
//
// Scope         : public
// Creator       : Pete Sholander, SNL
// Creation Date : 06/30/14
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    // this debugging code prints out cases where the input changes occur within the 
    // the model delay, which may cause incorrect output states in some cases
    if (inputStateChange_)
    {
      double delta = getSolverState().currTime_ - prevInputStateChangeTime_;
      if (delta <= model_.delay)
      {
        Xyce::dout() << "Device " << getName() << " changed state & accepted time step at "
		     << getSolverState().currTime_ << ".  Previous state" << std::endl  
                     << "change at Device " << getName() << " at time "  
                     << prevInputStateChangeTime_  << " was shorted-spaced.  Delta= "
                     << delta << std::endl << std::endl;
      }
      prevInputStateChangeTime_ = getSolverState().currTime_;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 digital instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEQVector" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
  }

  for (int i=0 ; i<numOutput ; ++i)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  qlo[" << i << "] = " << qlo[i] << std::endl;
      Xyce::dout() << "  qhi[" << i << "] = " << qhi[i] << std::endl;
    }

    if (li_Lo >= 0)
    {

      (*extData.daeQVectorPtr)[li_Lo] += qlo[i];
    }
    if (li_Hi >= 0)
    {

      (*extData.daeQVectorPtr)[li_Hi] += qhi[i];
    }


    (*extData.daeQVectorPtr)[li_Out[i]] -= qlo[i];

    (*extData.daeQVectorPtr)[li_Out[i]] -= qhi[i];
  }

  for (int i=0 ; i<numInput ; ++i)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  qref[" << i << "] = " << qref[i] << std::endl;
    }

    if (li_Ref >= 0)
    {

      (*extData.daeQVectorPtr)[li_Ref] += qref[i];
    }


    (*extData.daeQVectorPtr)[li_Inp[i]] -= qref[i];
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 digital instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEFVector" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
  }

  for (int i=0 ; i<numOutput ; ++i)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  rilo[" << i << "] = " << rilo[i] << std::endl;
      Xyce::dout() << "  rihi[" << i << "] = " << rihi[i] << std::endl;
    }

    if (li_Lo >= 0)
    {

      (*extData.daeFVectorPtr)[li_Lo] += rilo[i];
    }
    if (li_Hi >= 0)
    {

      (*extData.daeFVectorPtr)[li_Hi] += rihi[i];
    }


    (*extData.daeFVectorPtr)[li_Out[i]] -= rilo[i];

    (*extData.daeFVectorPtr)[li_Out[i]] -= rihi[i];
  }

  for (int i=0 ; i<numInput ; ++i)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  riref[" << i << "] = " << riref[i] << std::endl;
    }

    if (li_Ref >= 0)
    {

      (*extData.daeFVectorPtr)[li_Ref] += riref[i];
    }


    (*extData.daeFVectorPtr)[li_Inp[i]] -= riref[i];
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 digital instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  Linear::Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider <<std::endl;
    Xyce::dout() << "  Instance::loadDAEdQdx" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\nLoading DIGITAL dQdx matrix\n";
    Xyce::dout() << "Capacitance lo: " << model_.clo << std::endl;
    Xyce::dout() << "Capacitance hi: " << model_.chi << std::endl;
    Xyce::dout() << "Capacitance load: " << model_.cload << std::endl;
    Xyce::dout() << "DONE DIGITAL dQdx matrix LOAD\n";
  }

  for (int i=0 ; i<numInput ; ++i)
  {
    if (row_Ref >= 0)
    {

      (*dQdxMatPtr)[li_Ref][li_jac_Ref[i][1]] += model_.cload;

      (*dQdxMatPtr)[li_Ref][li_jac_Ref[i][2]] -= model_.cload;

      (*dQdxMatPtr)[li_Inp[i]][li_jac_Ref[i][4]] -= model_.cload;

      (*dQdxMatPtr)[li_Inp[i]][li_jac_Ref[i][5]] += model_.cload;
    }
    else
    {

      (*dQdxMatPtr)[li_Inp[i]][li_jac_Ref[i][1]] += model_.cload;
    }
  }

  for (int i=0 ; i<numOutput ; ++i)
  {
    if (row_Lo >= 0)
    {

      (*dQdxMatPtr)[li_Lo][li_jac_Lo[i][1]] += model_.clo;

      (*dQdxMatPtr)[li_Lo][li_jac_Lo[i][2]] -= model_.clo;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Lo[i][4]] -= model_.clo;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Lo[i][5]] += model_.clo;
    }
    else
    {

      (*dQdxMatPtr)[li_Out[i]][li_jac_Lo[i][1]] += model_.clo;
    }
    if (row_Hi >= 0)
    {

      (*dQdxMatPtr)[li_Hi][li_jac_Hi[i][1]] += model_.chi;

      (*dQdxMatPtr)[li_Hi][li_jac_Hi[i][2]] -= model_.chi;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Hi[i][4]] -= model_.chi;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Hi[i][5]] += model_.chi;
    }
    else
    {

      (*dQdxMatPtr)[li_Out[i]][li_jac_Hi[i][1]] += model_.chi;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 digital instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//                 For digital devices this is the contribution of the resistors
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider <<std::endl;
    Xyce::dout() << "  Instance::loadDAEdFdx" << std::endl;
  }

  for (int i=0 ; i<numInput ; ++i)
  {
    if (row_Ref >= 0)
    {

      (*dFdxMatPtr)[li_Ref][li_jac_Ref[i][1]] += model_.gload;

      (*dFdxMatPtr)[li_Ref][li_jac_Ref[i][2]] -= model_.gload;

      (*dFdxMatPtr)[li_Inp[i]][li_jac_Ref[i][4]] -= model_.gload;

      (*dFdxMatPtr)[li_Inp[i]][li_jac_Ref[i][5]] += model_.gload;
    }
    else
    {

      (*dFdxMatPtr)[li_Inp[i]][li_jac_Ref[i][1]] += model_.gload;
    }
  }

  for (int i=0 ; i<numOutput ; ++i)
  {
    if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  glo[" << i << "] = " << glo[i] << std::endl;
      Xyce::dout() << "  ghi[" << i << "] = " << ghi[i] << std::endl;
    }

    if (row_Lo >= 0)
    {

      (*dFdxMatPtr)[li_Lo][li_jac_Lo[i][1]] += glo[i];

      (*dFdxMatPtr)[li_Lo][li_jac_Lo[i][2]] -= glo[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Lo[i][4]] -= glo[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Lo[i][5]] += glo[i];
    }
    else
    {

      (*dFdxMatPtr)[li_Out[i]][li_jac_Lo[i][1]] += glo[i];
    }
    if (row_Hi >= 0)
    {

      (*dFdxMatPtr)[li_Hi][li_jac_Hi[i][1]] += ghi[i];

      (*dFdxMatPtr)[li_Hi][li_jac_Hi[i][2]] -= ghi[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Hi[i][4]] -= ghi[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Hi[i][5]] += ghi[i];
    }
    else
    {

      (*dFdxMatPtr)[li_Out[i]][li_jac_Hi[i][1]] += ghi[i];
    }
  }

  return true;
}

// Class Model
//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Model::processParams ()
{

  // If there are any time dependent parameters, set their values for
  // the current time.

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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)

{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (rload == 0)
  {
    UserError(*this) << "Zero load resistance in inputs";
  }
  gload = 1/rload;

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;

  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of digital instances: " << isize << std::endl;
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


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() 
     || (deviceMap.find("U")!=deviceMap.end())
     || (deviceMap.find("INV")!=deviceMap.end())
     || (deviceMap.find("NOT")!=deviceMap.end())
     || (deviceMap.find("AND")!=deviceMap.end())
     || (deviceMap.find("NAND")!=deviceMap.end())
     || (deviceMap.find("OR")!=deviceMap.end())
     || (deviceMap.find("NOR")!=deviceMap.end())
     || (deviceMap.find("ADD")!=deviceMap.end())
     || (deviceMap.find("XOR")!=deviceMap.end())
     || (deviceMap.find("NXOR")!=deviceMap.end())
     || (deviceMap.find("DFF")!=deviceMap.end())
     || (deviceMap.find("JKFF")!=deviceMap.end())
     || (deviceMap.find("TFF")!=deviceMap.end())
     || (deviceMap.find("BUF")!=deviceMap.end())
     || (deviceMap.find("DLTCH")!=deviceMap.end())
     || (deviceMap.find("DIG")!=deviceMap.end()))
  {
    // NOT device is deprecated now
    Config<Traits>::addConfiguration()
      .registerDevice("inv", 1)
      .registerDevice("not",1)
      .registerDevice("and", 1)
      .registerDevice("nand", 1)
      .registerDevice("or", 1)
      .registerDevice("nor", 1)
      .registerDevice("add", 1)
      .registerDevice("xor", 1)
      .registerDevice("nxor", 1)
      .registerDevice("dff", 1)
      .registerDevice("jkff", 1)
      .registerDevice("tff", 1)
      .registerDevice("buf", 1)
      .registerDevice("dltch", 1)
      .registerModelType("dig", 1);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceData::DeviceData
// Purpose       : Constructor for the DeviceData Class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
DeviceData::DeviceData(const char devLetter_):
  devLetter_('U')
{
  //no op
}

//-----------------------------------------------------------------------------
// Function      : DeviceData::DeviceData
// Purpose       : destructor for the DeviceData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
DeviceData::~DeviceData()
{}

//-----------------------------------------------------------------------------
// Function      : DeviceData::updatePowerPinLi
// Purpose       : updates li_Lo, li_Hi, li_Ref and numExtVars
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
void DeviceData::updatePowerPinLi(int& li_Lo, int& li_Hi, int& li_Ref, 
	       int& numExtVars, const bool vlo_given, const bool vhi_given,
               const bool vref_given)
{
  // no op
}

//-----------------------------------------------------------------------------
// Function      : UDeviceData::UDeviceData
// Purpose       : Constructor for the UDeviceData Class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
UDeviceData::UDeviceData(const char devLetter_):
  DeviceData('U')
{
  //Xyce::dout() << "Here in UDeviceData constructor\n";
}

//-----------------------------------------------------------------------------
// Function      : UDeviceData::UDeviceData
// Purpose       : destructor for the UDeviceData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
UDeviceData::~UDeviceData()
{}

//-----------------------------------------------------------------------------
// Function      : UDeviceData::updatePowerPinLi
// Purpose       : updates li_Lo, li_Hi, li_Ref and numExtVars for the U Device
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
void UDeviceData::updatePowerPinLi(int& li_Lo, int& li_Hi, int& li_Ref, 
               int& numExtVars, const bool vlo_given, const bool vhi_given,
               const bool vref_given)
{
  // For U devices, DPWR and DGND are always specified on the instance line.
  // Warning message if VHI, VLO or VREF are in the model card.
  li_Lo = 0;
  li_Hi = 0;
  li_Ref = 0;
  numExtVars += 2;
  //if (model_.given("VLO"))
  //{
  //  UserWarning(*this)<< "VLO model parameter ignored for U digital device";
  //}
  //if (model_.given("VHI"))
  //{
  //  UserWarning(*this)<< "VHI model parameter ignored for U digital device";
  //}
  //if (model_.given("VREF"))
  //{
  //  UserWarning(*this)<< "VREF model parameter ignored for U digital device";
  //}
}

//-----------------------------------------------------------------------------
// Function      : YDeviceData::YDeviceData
// Purpose       : Constructor for the DeviceData Class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
YDeviceData::YDeviceData(const char devLetter_):
  DeviceData('Y')
{
  //no op
}

//-----------------------------------------------------------------------------
// Function      : YDeviceData::~YDeviceData
// Purpose       : destructor for the YDeviceData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
YDeviceData::~YDeviceData()
{}

//-----------------------------------------------------------------------------
// Function      : YDeviceData::updatePowerPinLi
// Purpose       : updates li_Lo, li_Hi, li_Ref and numExtVars for the Y Device
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/7/2015
//-----------------------------------------------------------------------------
void YDeviceData::updatePowerPinLi(int& li_Lo, int& li_Hi, int& li_Ref, 
          int& numExtVars, const bool vlo_given, const bool vhi_given,
               const bool vref_given)
{
  // legacy code required to support VLO, VHI and VREF variables
  // being on the instance line rather than in model card in Y devices
  if (!vlo_given)    //!model_.given("VLO"))
  {
    li_Lo = 0;
    ++numExtVars;
  }
  if (!vhi_given)     //(!model_.given("VHI"))
  {
    li_Hi = 0;
    ++numExtVars;
  }
  if (!vref_given) //(!model_.given("VREF"))
  {
    li_Ref = 0;
    ++numExtVars;
  }
}

//-----------------------------------------------------------------------------
// Function      : GateData::GateData
// Purpose       : Constructor for the GateData Class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
GateData::GateData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  gateType_(""),
  devLetter_('U'),
  ilNumInput_(0)
{
  //no op
}

//-----------------------------------------------------------------------------
// Function      : GateData::~GateData
// Purpose       : destructor for the GateData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
GateData::~GateData()
{}

//-----------------------------------------------------------------------------
// Function      : GateData::evalTruthTable
// Purpose       : eval truth table and update output time for generic gate type.
//                 This function will be virtualized once the rest are done. 
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void GateData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	    std::vector<double>& oTime, const double lastT, const double delay, 
	    const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState) 
{
  // no op
}

//-----------------------------------------------------------------------------
// Function      : GateData::checkErrors
// Purpose       : check for instance line errors
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void GateData::checkErrors(const Instance& instance, const InstanceBlock& instance_block,
                           const int& iBase, const int& dev_numInputs) 
{
  if (instance.numExtVars != instance_block.numExtVars)
  {
    UserError(instance) << "Incorrect number of nodes in digital device."
                        << " Found " << instance_block.numExtVars << ", should be " << instance.numExtVars;
  }
}

//-----------------------------------------------------------------------------
// Function      : GateData::isClockLine
// Purpose       : returns true if inputPin is the clock line
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/12/2015
//-----------------------------------------------------------------------------
bool GateData::isClockLine(const int inputPin) 
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : GateData::setIC
// Purpose       : set the initial conditions at a gate
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void GateData::setIC(Instance& instance, const int pinNum)
{
  if (pinNum == 0)
  {
    if (instance.given("IC1"))
    {
      instance.outL[pinNum] = instance.ic1;
      instance.icGiven_[0] = true;
    }
    else 
    {
      instance.icGiven_[0] = false;
    }
  }
  else if (pinNum == 1)
  {
    if (instance.given("IC2"))
    {
      instance.outL[pinNum] = instance.ic2;
      instance.icGiven_[1] = true;
    }
    else 
    {
      instance.icGiven_[1] = false;
    }
  }
  else
  {
    DevelFatal(instance).in("GateData::setIC")
      << "Insufficient initial conditions supported in digital device";
  }
}

//-----------------------------------------------------------------------------
// Function      : GateData::getNumIO
// Purpose       : access the numInput and numOutput variables
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void GateData::getNumIO(int& input, int& output)
{
  input = numInput_;
  output = numOutput_;
}

//-----------------------------------------------------------------------------
// Function      : InvData::InvData
// Purpose       : Constructor for the INV and NOT gate types
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
InvData::InvData(const std::string gateType_, const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 1;
   numOutput_ = 1;
   type_ = INV;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : InvData::~InvData
// Purpose       : destructor for the InvData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
InvData::~InvData()
{}

//-----------------------------------------------------------------------------
// Function      : InvData::evalTruthTable
// Purpose       : eval truth table and update output time for INV and NOT gate types
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void InvData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  outL[0] = !inpL[0];
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : InvData::checkErrors
// Purpose       : check for instance line warnings for INV and NOT gate types
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/10/2015
//-----------------------------------------------------------------------------
void InvData::checkErrors(const Instance& instance, const InstanceBlock& instance_block, 
                    const int& iBase, const int& dev_numInputs)
{
  if (gateType_ == "NOT")
  {
    UserWarning(instance)<< "NOT gate type (" << instance.getName() << ") is deprecated. Consider using INV instead.";
  }
  GateData::checkErrors(instance, instance_block, iBase, dev_numInputs);
}

//-----------------------------------------------------------------------------
// Function      : AndData::AndData
// Purpose       : Constructor for the AND gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
AndData::AndData(const std::string gateType_, const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = (devLetter_ == 'Y') ? 2 : ilNumInput_;
   numOutput_ = 1;
   type_ = AND;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : AndData::~AndData
// Purpose       : destructor for the AndData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
AndData::~AndData()
{}

//-----------------------------------------------------------------------------
// Function      : AndData::evalTruthTable
// Purpose       : eval truth table and update output time for AND gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void AndData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking,const std::vector<bool>& oldState) 
{
  outL[0] = !(count(inpL.begin(),inpL.end(),false) > 0);
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : AndData::checkErrors
// Purpose       : check for instance line errors for AND gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/10/2015
//-----------------------------------------------------------------------------
void AndData::checkErrors(const Instance& instance, const InstanceBlock& instance_block, 
                    const int& iBase, const int& dev_numInputs)
{
  if (instance.numInput == 1)
  {
    UserError(instance) << "this device must have more than one input.";
  }
  if ( (dev_numInputs != 0) && (instance_block.numExtVars - iBase - instance.numOutput != dev_numInputs) )
  {
    UserError(instance) << "too few I/O nodes on instance line.";
  }
  GateData::checkErrors(instance, instance_block, iBase, dev_numInputs);
}

//-----------------------------------------------------------------------------
// Function      : NandData::NandData
// Purpose       : Constructor for the NAND gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
NandData::NandData(const std::string gateType_, const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = (devLetter_ == 'Y') ? 2 : ilNumInput_;
   numOutput_ = 1;
   type_ = NAND;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : NandData::~NandData
// Purpose       : destructor for the AndData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
NandData::~NandData()
{}

//-----------------------------------------------------------------------------
// Function      : NandData::evalTruthTable
// Purpose       : eval truth table and update output time for NAND gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void NandData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  outL[0] = (count(inpL.begin(),inpL.end(),false) > 0);
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : NandData::checkErrors
// Purpose       : check for instance line errors for NAND gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/10/2015
//-----------------------------------------------------------------------------
void NandData::checkErrors(const Instance& instance, const InstanceBlock& instance_block, 
                    const int& iBase, const int& dev_numInputs)
{
  if (instance.numInput == 1)
  {
    UserError(instance) << "this device must have more than one input.";
  }
  if ( (dev_numInputs != 0) && (instance_block.numExtVars - iBase - instance.numOutput != dev_numInputs) )
  {
    UserError(instance) << "too few I/O nodes on instance line.";
  }
  GateData::checkErrors(instance, instance_block, iBase, dev_numInputs);
}

//-----------------------------------------------------------------------------
// Function      : OrData::OrData
// Purpose       : Constructor for the OR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
OrData::OrData(const std::string gateType_, const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = (devLetter_ == 'Y') ? 2 : ilNumInput_;
   numOutput_ = 1;
   type_ = OR;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : OrData::~OrData
// Purpose       : destructor for the OrData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
OrData::~OrData()
{}

//-----------------------------------------------------------------------------
// Function      : OrData::evalTruthTable
// Purpose       : eval truth table and update output time for OR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void OrData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	     std::vector<double>& oTime, const double lastT, const double delay,
	     const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState) 
{
  outL[0] = outL[0] = (count(inpL.begin(),inpL.end(),true) > 0);
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : OrData::checkErrors
// Purpose       : check for instance line errors for OR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/10/2015
//-----------------------------------------------------------------------------
void OrData::checkErrors(const Instance& instance, const InstanceBlock& instance_block, 
                    const int& iBase, const int& dev_numInputs)
{
  if (instance.numInput == 1)
  {
    UserError(instance) << "this device must have more than one input.";
  }
  if ( (dev_numInputs != 0) && (instance_block.numExtVars - iBase - instance.numOutput != dev_numInputs) )
  {
    UserError(instance) << "too few I/O nodes on instance line.";
  }
  GateData::checkErrors(instance, instance_block, iBase, dev_numInputs);
}

//-----------------------------------------------------------------------------
// Function      : NorData::NorData
// Purpose       : Constructor for the NOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
NorData::NorData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = (devLetter_ == 'Y') ? 2 : ilNumInput_;
   numOutput_ = 1;
   type_ = NOR;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : NorData::~NorData
// Purpose       : destructor for the NorData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
NorData::~NorData()
{}

//-----------------------------------------------------------------------------
// Function      : NorData::evalTruthTable
// Purpose       : eval truth table and update output time for NOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void NorData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState) 
{
  outL[0] = !(count(inpL.begin(),inpL.end(),true) > 0);
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : NorData::checkErrors
// Purpose       : check for instance line errors for NOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/10/2015
//-----------------------------------------------------------------------------
void NorData::checkErrors(const Instance& instance, const InstanceBlock& instance_block, 
                    const int& iBase, const int& dev_numInputs)
{
  if (instance.numInput == 1)
  {
    UserError(instance) << "this device must have more than one input.";
  }
  if ( (dev_numInputs != 0) && (instance_block.numExtVars - iBase - instance.numOutput != dev_numInputs) )
  {
    UserError(instance) << "too few I/O nodes on instance line.";
  }
  GateData::checkErrors(instance, instance_block, iBase, dev_numInputs);
}

//-----------------------------------------------------------------------------
// Function      : AddData::AddData
// Purpose       : Constructor for the ADD gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
AddData::AddData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 3;
   numOutput_ = 2;
   type_ = ADD;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : AddData::~AddData
// Purpose       : destructor for the AddData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
AddData::~AddData()
{}

//-----------------------------------------------------------------------------
// Function      : AddData::evalTruthTable
// Purpose       : eval truth table and update output time for ADD gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void AddData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState) 
{
  outL[0] = inpL[0] ^ inpL[1] ^ inpL[2];
  // carry-out sum
  outL[1] = (inpL[0] & inpL[1]) | (inpL[1] & inpL[2]) | (inpL[0] & inpL[2]);

  oTime[0] = lastT+delay;
  oTime[1] = lastT+delay;
}

//-----------------------------------------------------------------------------
// Function      : XorData::XorData
// Purpose       : Constructor for the XOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
XorData::XorData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 2;
   numOutput_ = 1;
   type_ = XOR;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : XorData::~XorData
// Purpose       : destructor for the XorData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
XorData::~XorData()
{}

//-----------------------------------------------------------------------------
// Function      : XorData::evalTruthTable
// Purpose       : eval truth table and update output time for XOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void XorData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  outL[0] = inpL[0] ^ inpL[1];
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : NxorData::NxorData
// Purpose       : Constructor for the NXOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
NxorData::NxorData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 2;
   numOutput_ = 1;
   type_ = NXOR;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : NxorData::~NxorData
// Purpose       : destructor for the NxorData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
NxorData::~NxorData()
{}

//-----------------------------------------------------------------------------
// Function      : NxorData::evalTruthTable
// Purpose       : eval truth table and update output time for NXOR gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void NxorData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldStaVector)
{
  outL[0] = !(inpL[0] ^ inpL[1]);
  oTime[0] = lastT + delay; // delay is model_.delay
}

//-----------------------------------------------------------------------------
// Function      : DffData::DffData
// Purpose       : Constructor for the DFF gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
DffData::DffData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 4;  //PREB, CLRB, clock, data
   numOutput_ = 2; // Q, Q_bar
   clockPin_ = 2;
   type_ = DFF;
   supportsXState_ = true;
}

//-----------------------------------------------------------------------------
// Function      : DffData::~DffData
// Purpose       : destructor for the DffData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
DffData::~DffData()
{}

//-----------------------------------------------------------------------------
// Function      : Dff::evalTruthTable
// Purpose       : eval truth table and update output time for DFF gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void DffData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	    std::vector<double>& oTime, const double lastT, const double delay,
	    const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  // DFF: in0: PREB, in1: CLRB, in2: clock, in3: data
  // DFF: out0: Q, out1: Q_bar
  // CD4013B: set = !PREB, reset = !CLRB
  // CD4013B: in0: set, in1: reset, in2: clock, in3: data
  // CD4013B: out0: Q, out1: Q_bar
  if (clocking && inpL[2] ==1) //clock rising edge 0->1
  {
    if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
    {
      outL[0] = inpL[3];  //Q = D
      outL[1] = !(inpL[3]); //Q_bar = !D
    }
  }
  else if (clocking && inpL[2] ==0) //clock falling edge 1->0
  {
    if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
    {
      //oldState0 = oldStaVector[li_currentStateOut[0]]
      //oldState1 = oldStaVector[li_currentStateOut[1]]
      outL[0] = oldState[0];
      outL[1] = oldState[1];
    }
  }
  else // no clock change
  {
    if (inpL[0] == 1 && inpL[1] == 0) //PREB=1, CLRB = 0
    {
      outL[0] = 0;
      outL[1] = 1;
    }
    else if (inpL[0] == 0 && inpL[1] == 1) //PREB = 0, CLRB = 1
    {
      outL[0] = 1;
      outL[1] = 0;
    }
    else if (inpL[0] == 0 && inpL[1] == 0) //PREB = CLRB = 0
    {
      outL[0] = 1;
      outL[1] = 1;
    }
    else if (dcopFlag) // getSolverState().dcopFlag
    { // Handle dcop calculation when enable is low.  This forces Q and Qbar to be
      // different at DCOP.  This behavior differs from PSpice, where Q and Qbar
      // would be in an indeterminate state that is about halfway between V_LO and V_HI.
      outL[0] = inpL[3];
      outL[1] = !inpL[3];
    }
    else if (outL[1] == outL[0])
    { // handles unstable state after PREB/CLRB = 0 state
      outL[1] = !outL[0];
    }
    else
    {
      // no op.  Keep outputs in current state
    }
  }

  //Xyce::dout() << "Clocking, dcop and I/O are " << clocking << "," << dcopFlag 
  //             << " and (" << inpL[0] << "," << inpL[1] 
  //             << "," << inpL[2] << "," << inpL[3] << ") and ("
  //             << outL[0] << "," << outL[1] << ")" << std::endl; 

  oTime[0] = lastT+delay; // delay = model_.delay
  oTime[1] = lastT+delay; 
}

//-----------------------------------------------------------------------------
// Function      : DffData::isClockLine
// Purpose       : returns true if inputPin is the clock line
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/12/2015
//-----------------------------------------------------------------------------
bool DffData::isClockLine(const int inputPin) 
{
  // Pin 2 is the clk pin in DFF.  This is defined in the class definition in the
  // header file.
  return inputPin == clockPin_ ? true : false;
}

//-----------------------------------------------------------------------------
// Function      : DffData::setIC
// Purpose       : set the initial conditions at a DFF gate
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/14/2015
//-----------------------------------------------------------------------------
void DffData::setIC(Instance& instance, const int pinNum)
{
  if (pinNum == 0)
  {
    if (instance.given("IC1"))
    {
      instance.outL[pinNum] = instance.ic1;
      instance.icGiven_[0] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[0] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[0] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[0] = false;
    }
  }
  else if (pinNum == 1)
  {
    if (instance.given("IC2"))
    {
      instance.outL[pinNum] = instance.ic2;
      instance.icGiven_[1] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[1] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[1] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[1] = false;
    }
  }
  else
  {
    DevelFatal(instance).in("DffData::setIC")
      << "Insufficient initial conditions supported in digital device";
  }
}

//-----------------------------------------------------------------------------
// Function      : JkffData::JkffData
// Purpose       : Constructor for the JKFF gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
JkffData::JkffData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 5;  // PREB, CLRB, CLK, J, K
   numOutput_ = 2; // Q, Q_bar
   clockPin_ = 2;
   type_ = JKFF;
   supportsXState_ = true;
}

//-----------------------------------------------------------------------------
// Function      : JkffData::~JkffData
// Purpose       : destructor for the JKFFData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
JkffData::~JkffData()
{}

//-----------------------------------------------------------------------------
// Function      : Jkff::evalTruthTable
// Purpose       : eval truth table and update output time for JKFF gate type
// Special Notes : This is still the truth table for the JKFF
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2015
//-----------------------------------------------------------------------------
void JkffData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	    std::vector<double>& oTime, const double lastT, const double delay,
	    const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  // JKFF: in0: PREB, in1: CLRB, in2: clock, in3: J, in4: K
  // JKFF: out0: Q, out1: Q_bar
 
  if (clocking && inpL[2] ==0) //clock falling edge 1->0.  PSpice JKFF uses falling-edge
  {
    if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
    {
      if (inpL[3] == 0 && inpL[4] == 0)
      {
        // no op.  hold state.
      }
      else if (inpL[3] == 0 && inpL[4] == 1)
      {
        // reset
        outL[0] = 0;
        outL[1] = 1;
      } 
      else if (inpL[3] == 1 && inpL[4] == 0)
      {
        // set
        outL[0] = 1;
        outL[1] = 0;
      }
      else if (inpL[3] == 1 && inpL[4] == 1)
      {
        // toggle
        outL[0] = !oldState[0];
        outL[1] = !outL[0];
      }
    }
  }
  else // no clock change
  {
    if (inpL[0] == 1 && inpL[1] == 0) //PREB=1, CLRB = 0
    {
      outL[0] = 0;
      outL[1] = 1;
    }
    else if (inpL[0] == 0 && inpL[1] == 1) //PREB = 0, CLRB = 1
    {
      outL[0] = 1;
      outL[1] = 0;
    }
    else if (inpL[0] == 0 && inpL[1] == 0) //PREB = CLRB = 0
    {
      outL[0] = 1;
      outL[1] = 1;
    }
    else if (dcopFlag) // getSolverState().dcopFlag
    { // Handle dcop calculation when enable is low.  This forces Q and Qbar to be
      // different at DCOP.  This behavior differs from PSpice, where Q and Qbar
      // would be in an indeterminate state that is about halfway between V_LO and V_HI.
      outL[0] = inpL[3];
      outL[1] = !inpL[3];
    }
    else if (outL[1] == outL[0])
    { // handles unstable state after PREB/CLRB = 0 state
      outL[1] = !outL[0];
    }
    else
    {
      // no op.  Keep outputs in current state
    }
  }

  //Xyce::dout() << "Clocking, dcop and I/O are " << clocking << "," << dcopFlag 
  //             << " and (" << inpL[0] << "," << inpL[1] 
  //             << "," << inpL[2] << "," << inpL[3] << ") and ("
  //             << outL[0] << "," << outL[1] << ")" << std::endl; 

  oTime[0] = lastT+delay; // delay = model_.delay
  oTime[1] = lastT+delay; 
}

//-----------------------------------------------------------------------------
// Function      : JkffData::isClockLine
// Purpose       : returns true if inputPin is the clock line
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
bool JkffData::isClockLine(const int inputPin) 
{
  // Pin 2 is the clk pin in JKFF.  This is defined in the class definition in the
  // header file.
  return inputPin == clockPin_ ? true : false;
}

//-----------------------------------------------------------------------------
// Function      : JkffData::setIC
// Purpose       : set the initial conditions at a JKFF gate
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
void JkffData::setIC(Instance& instance, const int pinNum)
{
  if (pinNum == 0)
  {
    if (instance.given("IC1"))
    {
      instance.outL[pinNum] = instance.ic1;
      instance.icGiven_[0] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[0] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[0] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[0] = false;
    }
  }
  else if (pinNum == 1)
  {
    if (instance.given("IC2"))
    {
      instance.outL[pinNum] = instance.ic2;
      instance.icGiven_[1] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[1] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[1] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[1] = false;
    }
  }
  else
  {
    DevelFatal(instance).in("Jkff::setIC")
      << "Insufficient initial conditions supported in digital device";
  }
}

//-----------------------------------------------------------------------------
// Function      : TffData::TffData
// Purpose       : Constructor for the TFF gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
TffData::TffData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 2;  //T, clock
   numOutput_ = 2; // Q, Q_bar
   clockPin_ = 1;
   type_ = TFF;
   supportsXState_ = true;
}

//-----------------------------------------------------------------------------
// Function      : TffData::~TffData
// Purpose       : destructor for the TffData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2015
//-----------------------------------------------------------------------------
TffData::~TffData()
{}

//-----------------------------------------------------------------------------
// Function      : Tff::evalTruthTable
// Purpose       : eval truth table and update output time for TFF gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void TffData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	    std::vector<double>& oTime, const double lastT, const double delay,
	    const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  // TFF: in0 T, in1: Clock
  // TFF: out0: Q, out1: Q_bar
  // Note: TFF uses rising-edge clocking
  
  if (clocking && inpL[0] && inpL[1] ==1) //Toggle bit (T) is high and clock rising edge 0->1)
  {
    // toggle outputs, if T is high and a clock pulse happened
    outL[0] = !oldState[0];
    outL[1] = !outL[0];
  }
  else if (dcopFlag) // getSolverState().dcopFlag
  { // Handle dcop calculation when enable is low.  This forces Q and Qbar to be
    // different at DCOP.  
    outL[0] = inpL[0];
    outL[1] = !inpL[0];
  }
  else
  {
     // no output changes if T is low, or if this a falling edge
  }

  // update output times
  oTime[0] = lastT+delay; // delay = model_.delay
  oTime[1] = lastT+delay; 
}

//-----------------------------------------------------------------------------
// Function      : TffData::isClockLine
// Purpose       : returns true if inputPin is the clock line
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
bool TffData::isClockLine(const int inputPin) 
{
  // Pin 1 is the clk pin in TFF.  This is defined in the class definition in the
  // header file.
  return inputPin == clockPin_ ? true : false;
}

//-----------------------------------------------------------------------------
// Function      : TffData::setIC
// Purpose       : set the initial conditions at a TFF gate
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 7/21/2016
//-----------------------------------------------------------------------------
void TffData::setIC(Instance& instance, const int pinNum)
{
  if (pinNum == 0)
  {
    if (instance.given("IC1"))
    {
      instance.outL[pinNum] = instance.ic1;
      instance.icGiven_[0] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[0] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[0] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[0] = false;
    }
  }
  else if (pinNum == 1)
  {
    if (instance.given("IC2"))
    {
      instance.outL[pinNum] = instance.ic2;
      instance.icGiven_[1] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[1] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[1] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[1] = false;
    }
  }
  else
  {
    DevelFatal(instance).in("TffData::setIC")
      << "Insufficient initial conditions supported in digital device";
  }
}

//-----------------------------------------------------------------------------
// Function      : DltchData::DltchData
// Purpose       : Constructor for the DLTCH gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
DltchData::DltchData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 4;  //PREB, CLRB, enable, data
   numOutput_ = 2; // Q, Q_bar
   type_ = DLTCH;
   supportsXState_ = true;
}

//-----------------------------------------------------------------------------
// Function      : DltchData::~DltchData
// Purpose       : destructor for the DltchData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
DltchData::~DltchData()
{}

//-----------------------------------------------------------------------------
// Function      : DltchData::evalTruthTable
// Purpose       : eval truth table and update output time for DLTCH gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void DltchData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState) 
{
  // DLTCH: in0: PREB, in1: CLRB, in2: enable, in3: data
  // DLTCH: out0: Q, out1: Q_bar
  if ((inpL[0] == 1) && (inpL[1] == 0))
  {
    outL[0] = 0;
    outL[1] = 1;
    oTime[0] = lastT+delay;  // delay = model_.delay
    oTime[1] = lastT+delay;
  }
  else if ((inpL[0] == 0) && (inpL[1] == 1))
  {
    outL[0] = 1;
    outL[1] = 0;
    oTime[0] = lastT+delay;
    oTime[1] = lastT+delay;
  }
  else if ((inpL[0] == 0) && (inpL[1] == 0))
  { // this state is unstable.  The transition out of this state is handled by the
    // else if (outL[1] == outL[0]) block below 
    outL[0] = 1;
    outL[1] = 1;
    oTime[0] = lastT+delay;
    oTime[1] = lastT+delay;
  }
  else if (inpL[2] == 1)
  { // enable line, PREB and CLRB are TRUE
    outL[0] = inpL[3];
    outL[1] = !inpL[3];
    oTime[0] = lastT+delay;
    oTime[1] = lastT+delay;
  }
  else if (dcopFlag) // getSolverState().dcopFlag
  { // Handle dcop calculation when enable is low.  This forces Q and Qbar to be
    // different at DCOP.  This behavior differs from PSpice, where Q and Qbar
    // would be in an indeterminate state that is about halfway between V_LO and V_HI.
    outL[0] = inpL[3];
    outL[1] = !inpL[3];
  }
  else if (outL[1] == outL[0])
  { // handles unstable state after PREB/CLRB = 0 state
    outL[1] = !outL[0];
    oTime[1] = lastT+delay;
  }
  else
  {
    // no op.  Keep outputs latched in current state
  }
}

//-----------------------------------------------------------------------------
// Function      : DltchData::setIC
// Purpose       : set the initial conditions at a DLTCH gate
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/14/2015
//-----------------------------------------------------------------------------
void DltchData::setIC(Instance& instance, const int pinNum)
{
  if (pinNum == 0)
  {
    if (instance.given("IC1"))
    {
      instance.outL[pinNum] = instance.ic1;
      instance.icGiven_[0] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[0] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[0] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[0] = false;
    }
  }
  else if (pinNum == 1)
  {
    if (instance.given("IC2"))
    {
      instance.outL[pinNum] = instance.ic2;
      instance.icGiven_[1] = true;
    }
    else if (instance.digInitState_ == 0)
    {
      instance.outL[pinNum] = 1;
      instance.icGiven_[1] = false;
    }
    else if (instance.digInitState_ == 1)
    {
      instance.outL[pinNum] = 0;
      instance.icGiven_[1] = false;
    }
    else  // covers digInitState_ == 2 || 3
    {
      instance.icGiven_[1] = false;
    }
  }
  else
  {
    DevelFatal(instance).in("DltchData::setIC")
      << "Insufficient initial conditions supported in digital device";
  }
}

//-----------------------------------------------------------------------------
// Function      : BufData::BufData
// Purpose       : Constructor for the Buf gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/3/2015
//-----------------------------------------------------------------------------
BufData::BufData(const std::string gateType_,const char devLetter_, const int ilNumInput_):
  GateData(gateType_,devLetter_,ilNumInput_)
{
   numInput_ = 1;  
   numOutput_ = 1; 
   type_= BUF;
   supportsXState_ = false;
}

//-----------------------------------------------------------------------------
// Function      : BufData::~BufData
// Purpose       : destructor for the BufData class
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
BufData::~BufData()
{}

//-----------------------------------------------------------------------------
// Function      : BufData::evalTruthTable
// Purpose       : eval truth table and update output time for BUF gate type
// Special Notes :
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 8/6/2015
//-----------------------------------------------------------------------------
void BufData::evalTruthTable(const std::vector<bool> inpL, std::vector<bool>& outL, 
	      std::vector<double>& oTime, const double lastT, const double delay,
	      const bool dcopFlag, const bool clocking, const std::vector<bool>& oldState)
{
  outL[0] = inpL[0];
  oTime[0] = lastT + delay; // delay is model_.delay
}

} // namespace Digital
} // namespace Device
} // namespace Xyce

//-------------------------------------------------------------------------
//   Copyright 2002-2021 National Technology & Engineering Solutions of
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
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//----------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ThermalResistor.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Resistor.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_Expression.h>
#include <N_UTL_FeatureTest.h>

namespace Xyce {
namespace Device {

namespace ThermalResistor {

void Traits::loadInstanceParameters(ParametricData<ThermalResistor::Instance> &p)
{
    // Set up double precision variables:
    p.addPar ("R",1000.0,&ThermalResistor::Instance::R)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Resistance");

    p.addPar("M", 1.0, &ThermalResistor::Instance::multiplicityFactor)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Multiplicity Factor");

    p.addPar ("L",0.0,&ThermalResistor::Instance::length)
     .setUnit(U_METER)
     .setCategory(CAT_NONE)
     .setDescription("Length of conductor");

    p.addPar ("W",0.0,&ThermalResistor::Instance::width)
     .setUnit(U_METER)
     .setCategory(CAT_NONE)
     .setDescription("Width of conductor");

    p.addPar ("A",0.0,&ThermalResistor::Instance::area)
     .setUnit(U_METER2)
     .setCategory(CAT_NONE)
     .setDescription("Area of conductor");

    p.addPar ("THERMAL_L",0.0,&ThermalResistor::Instance::thermalLength)
     .setUnit(U_METER)
     .setCategory(CAT_NONE)
     .setDescription("Length of material thermally coupled to conductor");

    p.addPar ("THERMAL_A",0.0,&ThermalResistor::Instance::thermalArea)
     .setUnit(U_METER2)
     .setCategory(CAT_NONE)
     .setDescription("Area of material thermally coupled to conductor");

     // This stuff is copied from the model:
    p.addPar ("RESISTIVITY",0.0,&ThermalResistor::Instance::resistivity)
     .setUnit(U_OHMM)
     .setCategory(CAT_NONE)
     .setDescription("Resistor material resistivity");

    p.addPar ("DENSITY",0.0,&ThermalResistor::Instance::density)
     .setUnit(U_KGMM3)
     .setCategory(CAT_NONE)
     .setDescription("Resistor material density (unused)");

    p.addPar ("HEATCAPACITY",0.0,&ThermalResistor::Instance::heatCapacity)
     .setUnit(U_JMM3KM1)
     .setCategory(CAT_NONE)
     .setDescription("Resistor material volumetric heat capacity");

    p.addPar ("THERMAL_HEATCAPACITY",0.0,&ThermalResistor::Instance::thermalHeatCapacity)
     .setUnit(U_JMM3KM1)
     .setCategory(CAT_NONE)
     .setDescription("Volumetric heat capacity of material thermally coupled to conductor");

    p.addPar ("TEMP",0.0,&ThermalResistor::Instance::temp)
     .setUnit(U_DEGC)
     .setCategory(CAT_NONE)
     .setDescription("Device temperature");

    // Set up non-double precision variables:
    p.addPar ("OUTPUTINTVARS",false, &ThermalResistor::Instance::outputInternalVarsFlag)
     .setUnit(U_NONE)
     .setCategory(CAT_CONTROL)
     .setDescription("Debug Output switch");
}

void Traits::loadModelParameters(ParametricData<ThermalResistor::Model> &p)
{
    p.addPar ("TC1",0.0,&ThermalResistor::Model::tempCoeff1)
     .setUnit(U_DEGCM1)
     .setCategory(CAT_NONE)
     .setDescription("Linear Temperature Coefficient");

    p.addPar ("TC2",0.0,&ThermalResistor::Model::tempCoeff2)
     .setUnit(U_DEGCM2)
     .setCategory(CAT_NONE)
     .setDescription("Quadratic Temperature Coefficient");

    p.addPar("TCE", 0.0, &ThermalResistor::Model::tempCoeffExp)
     .setGivenMember(&ThermalResistor::Model::tempCoeffExpModelGiven)
     .setUnit(U_PERCENTDEGCM1)
     .setCategory(CAT_NONE)
     .setDescription("Exponential Temperature Coefficient");

    p.addPar ("RSH",0.0,&ThermalResistor::Model::sheetRes)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("Sheet Resistance");

    p.addPar("R",1.0,&ThermalResistor::Model::resistanceMultiplier)
     .setUnit(U_NONE)
     .setDescription("Resistance Multiplier");

    p.addPar ("RESISTIVITY",0.0,&ThermalResistor::Model::resistivity)
     .setUnit(U_OHMM)
     .setCategory(CAT_NONE)
     .setDescription("Resistor material resistivity");

    p.addPar ("DENSITY",0.0,&ThermalResistor::Model::density)
     .setUnit(U_KGMM3)
     .setCategory(CAT_NONE)
     .setDescription("Resistor material density (unused)");

    p.addPar ("HEATCAPACITY",0.0,&ThermalResistor::Model::heatCapacity)
     .setUnit(U_JMM3KM1)
     .setCategory(CAT_NONE)
     .setDescription("Resistor material volumetric heat capacity");

    p.addPar ("THERMAL_HEATCAPACITY",0.0,&ThermalResistor::Model::thermalHeatCapacity)
     .setUnit(U_JMM3KM1)
     .setCategory(CAT_NONE)
     .setDescription("Volumetric heat capacity of material thermally coupled to conductor");

    p.addPar ("DEFW",1.e-5,&ThermalResistor::Model::defWidth)
     .setUnit(U_METER)
     .setCategory(CAT_NONE)
     .setDescription("Default Instance Width");

    p.addPar ("NARROW",0.0,&ThermalResistor::Model::narrow)
     .setUnit(U_METER)
     .setCategory(CAT_NONE)
     .setDescription("Narrowing due to side etching");

    p.addPar ("TNOM",0.0,&ThermalResistor::Model::tnom)
     .setUnit(STANDARD)
     .setCategory(CAT_NONE)
     .setDescription("Parameter Measurement Temperature");
}

std::vector< std::vector<int> > Instance::jacStamp;


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  // M must be non-negative
  if (multiplicityFactor <= 0)
  {
    UserError(*this) << "Multiplicity Factor (M) must be non-negative" << std::endl;
  }

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model &               Riter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Riter),
    R(0.0),
    multiplicityFactor(0.0),
    length(0.0),
    width(0.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    G(0.0),
    i0(0.0),
    li_Pos(-1),
    li_Neg(-1),
    tempModelEnabled(false),
    outputInternalVarsFlag(false),
    li_TempState(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0)
{
  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 0;
  setNumBranchDataVars(0);     // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1; // this is the space to allocate if lead current or power is needed.

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
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.getImmutableValue<double>();
  if (!given("W"))
    width = model_.defWidth;;

  // Nonzero value for numStateVars indicates that the self-consistent thermal
  // resistor model is being used.
  if (given("A") && given("L") &&
      ( (model_.given("HEATCAPACITY") && model_.given("RESISTIVITY"))
                || (given("HEATCAPACITY") && given("RESISTIVITY")) )
       )
  {
    numStateVars++;
    tempModelEnabled = true;
  }

  // If the instance parameters are NOT given, but the model parameters
  // ARE given, then copy the model params into the instance.  All the
  // work is actually done in the instance anyway,
  // but this stuff can be specified on the model level.
  if ( !(given("HEATCAPACITY")) && !(given("RESISTIVITY")) &&
      model_.given("HEATCAPACITY") && model_.given("RESISTIVITY")  )
  {
    resistivity = model_.resistivity;
    heatCapacity = model_.heatCapacity;
    thermalHeatCapacity = model_.thermalHeatCapacity;

    // copy over the dependent parameters.  For now, it only appears necessary
    // to copy the dependentParams vector, and not anything else like
    // the expVarLIDs vector.

    if (!(model_.getDependentParams().empty()))
    {
      const std::vector<Depend> & model_dp = model_.getDependentParams();
      int dpSize = model_dp.size();

      for (int i=0;i<dpSize;++i)
      {
        Depend dpTmp;
        dpTmp.name = model_dp[i].name;
        dpTmp.vals = model_dp[i].vals;
        dpTmp.global_params = model_dp[i].global_params;
        dpTmp.n_vars = model_dp[i].n_vars;
        dpTmp.lo_var = model_dp[i].lo_var;
        dpTmp.vectorIndex = -1;

        // dpTmp needs to point to a copy of the original expression.
        dpTmp.expr = new Util::Expression( *(model_dp[i].expr) );

        double *Dval;
        if (dpTmp.name=="RESISTIVITY")
        {
          //Dval = &Instance::resistivity;
          Dval = &resistivity;
          dpTmp.resultU.result = Dval;

        }

        if (dpTmp.name=="HEATCAPACITY")
        {
          Dval = &heatCapacity;
          dpTmp.resultU.result = Dval;
        }

        addDependentParameter(dpTmp);
      }
    }
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                           const std::vector<int> & extLIDVecRef )
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS))
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  ResistorInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
  }

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << section_divider << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Note that the resistor does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/12/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  if (numStateVars > 0)
    li_TempState = staLIDVec[0];

}

// ----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/24/2013
// ----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef)
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerBranchDataLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "branch data" vector.  As with other such vectors, the device
/// declares at construction time how many branch vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// The Resistor device uses exactly one "branch data vector" element, where
/// it keeps the "lead current" that may be used on .PRINT lines as
/// "I(R1)" for the current through resistor R1. and a junction voltage.
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012

void Instance::registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef)
{
  AssertLIDs(branchLIDVecRef.size() == getNumBranchDataVars());

  if (loadLeadCurrent)
  {
    li_branch_data= branchLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
#endif
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
  if (loadLeadCurrent)
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233.
// Creation Date : 3/05/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  if (tempModelEnabled)
  {
    Linear::Vector * staVectorPtr = extData.currStaVectorPtr;

    if (!getSolverState().dcopFlag)
    {
      if (li_TempState >= 0)
      {
        temp = (*staVectorPtr)[li_TempState];
        updateTemperature(temp);
      }
    }
  }

  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];

  // Load RHS vector element for the positive circuit node KCL equ.
  i0 = (v_pos-v_neg)*G;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = updateIntermediateVars ();

  if (tempModelEnabled)
  {
    double * staVec = extData.nextStaVectorRawPtr;
    if (li_TempState >= 0)
    {
      double dissipation = i0*i0*R;
      temp += dissipation*getSolverState().currTimeStep_/(area*length*heatCapacity +
                            thermalArea*thermalLength*thermalHeatCapacity);
      staVec[li_TempState] = temp;
    }
  }

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles(bool force_final_output)
{
  bool bsuccess = true;

  if (tempModelEnabled && outputInternalVarsFlag)
  {
    Xyce::dout().width(28); Xyce::dout().precision(20); Xyce::dout().setf(std::ios::scientific);
    Linear::Vector * sta1VectorPtr = extData.nextStaVectorPtr;
    Linear::Vector * sta2VectorPtr = extData.currStaVectorPtr;
    Xyce::dout() << "TEMP("<<getName()<<"):  " << getSolverState().currTime_ << "    "
        << ((*sta1VectorPtr)[li_TempState]-CONSTCtoK) << "    "
        << ((*sta2VectorPtr)[li_TempState]-CONSTCtoK)
        << std::endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  fVec[li_Pos] += i0;
  fVec[li_Neg] += -i0;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    double * solVec = extData.nextSolVectorRawPtr;
    leadF[li_branch_data] = i0;
    junctionV[li_branch_data] = solVec[li_Pos] - solVec[li_Neg];
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  double difference, factor, tempCorrFactor;

  if (tempModelEnabled)
  {
    // ERK; this conditional is necessary with new expression.  The thermal model doesn't
    // make sense for anything other than transient.
    if( !(getSolverState().tranopFlag) && !(getSolverState().dcopFlag) )
    {
      updateDependentParameters(temp_tmp);
    }
    R = resistivity * length / area;
    // Apply the multiplicityFactor factor M (from the instance line).
    factor = 1/multiplicityFactor;
  }
  else
  {
    // if the self-heating model isn't being used then the temperature
    // coefficients (TC1, TC2 and TCE) from the model will be used.
    if (!given("R") && numStateVars == 0)
    {
      if (model_.given("RSH") && given("L") && (model_.sheetRes!=0) &&
          (length != 0))
      {
        R = model_.sheetRes * (length - model_.narrow)
          / (width - model_.narrow);
      }
      else
      {
        R=1000;
        UserWarning(*this) << "Resistance is set to 0, setting to the default, " << R << " ohms";
      }
    }

    difference = temp_tmp - model_.tnom;

    if (model_.tempCoeffExpModelGiven)
    {
      // Exponential temperature coefficient takes precedence, if it appears either
      // in the model or on the instance line.  This formula matches what PSpice 
      // implements.
      tempCorrFactor = pow(1.01,(model_.tempCoeffExp*difference));
    }
    else
    {
      // otherwise use the linear and quadratic temperature coefficients. Their
      // default values (when not given) are zero.
      tempCorrFactor = 1.0 + (model_.tempCoeff1)*difference +
                             (model_.tempCoeff2)*difference*difference;
    }
    // Apply both the resistanceMultipler (from the model card), the multiplicity factor M 
    // (from the instance line) and the temperature correction.
    factor = model_.resistanceMultiplier*tempCorrFactor/multiplicityFactor;
  }

  if (R*factor != 0.0)
    G = 1.0/(R * factor);
  else
    G = 0.0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
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
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tempCoeffExp(0.0),
    tempCoeffExpModelGiven(false),
    sheetRes(0.0),
    resistanceMultiplier(1.0),
    defWidth(10e-6),
    narrow(0.0),
    tnom(getDeviceOptions().tnom)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  if (!given("THERMAL_HEATCAPACITY"))
    thermalHeatCapacity = heatCapacity;

  // calculate dependent (ie computed) params and check for errors:
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
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

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of Resistor Instances: " << isize << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << "\t\tR(Tnom) = " << (*iter)->R;
    os << "\tG(T) = " << (*iter)->G;
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


// ThermalResistor Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
    {
      Instance & ri = *(*it);

      bool btmp = ri.updateIntermediateVars ();
      bsuccess = bsuccess && btmp;

      if (ri.tempModelEnabled)
      {
        if (ri.li_TempState >= 0)
        {
          double dissipation = ri.i0*ri.i0*ri.R;
          double dt = ri.getSolverState().currTimeStep_;
          ri.temp += dissipation*dt/(ri.area*ri.length*ri.heatCapacity +
                                ri.thermalArea*ri.thermalLength*ri.thermalHeatCapacity);
          staVec[ri.li_TempState] = ri.temp;
        }
      }
    }
//  }

  return bsuccess;
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
    Instance & ri = *(*it);
    fVec[ri.li_Pos] += ri.i0;
    fVec[ri.li_Neg] += -ri.i0;
    if( ri.loadLeadCurrent )
    {
      leadF[ri.li_branch_data] = ri.i0;
      junctionV[ri.li_branch_data] = solVec[ri.li_Pos] - solVec[ri.li_Neg];
    }
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
    Instance & ri = *(*it);
#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *(ri.f_PosEquPosNodePtr) += ri.G;

    *(ri.f_PosEquNegNodePtr) -= ri.G;

    *(ri.f_NegEquPosNodePtr) -= ri.G;

    *(ri.f_NegEquNegNodePtr) += ri.G;
#else

    dFdx[ri.li_Pos][ri.APosEquPosNodeOffset] += ri.G;

    dFdx[ri.li_Pos][ri.APosEquNegNodeOffset] -= ri.G;

    dFdx[ri.li_Neg][ri.ANegEquPosNodeOffset] -= ri.G;

    dFdx[ri.li_Neg][ri.ANegEquNegNodeOffset] += ri.G;
#endif
  }

  return true;
}

Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void 
registerDevice(const DeviceCountMap& deviceMap, const std::set<int>& levelSet)
{
  if (deviceMap.empty() || (deviceMap.find("R")!=deviceMap.end() && levelSet.find(2)!=levelSet.end()))
  {
    Config<Traits>::addConfiguration()
      .registerDevice("r", 2)
      .registerModelType("r", 2);
  }
}

} // namespace ThermalResistor
} // namespace Device
} // namespace Xyce

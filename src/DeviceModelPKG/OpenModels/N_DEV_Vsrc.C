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
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
//
//
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceMgr.h>
#include <N_UTL_AssemblyTypes.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_FeatureTest.h>

#include <N_UTL_Math.h>

//#include <Teuchos_RCP.hpp>
//#include <N_UTL_FFTInterface.hpp>

 
#include <N_UTL_MachDepParams.h>

namespace Xyce {
namespace Device {

namespace Vsrc {

void Traits::loadInstanceParameters(ParametricData<Vsrc::Instance> &p)
{
    p.addPar ("DCV0", 0.0, &Vsrc::Instance::DCV0)
      .setOriginalValueStored(true)
      .setUnit(U_VOLT)
      .setDescription("DC Voltage")
      .setAnalyticSensitivityAvailable(true)
      .setSensitivityFunctor(&dcv0Sens);

    // Pulse parameters
    p.addPar ("V0", 0.0, &Vsrc::Instance::par0)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Offset Voltage");

    p.addPar ("V1", 0.0, &Vsrc::Instance::par0)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Initial Voltage");

    p.addPar ("V2", 0.0, &Vsrc::Instance::par1)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Pulsed Voltage");

    p.addPar ("TD", 0.0, &Vsrc::Instance::par2)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Delay");

    p.addPar ("TR", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Rise Time");

    p.addPar ("TF", 0.0, &Vsrc::Instance::par4)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Fall Time");

    p.addPar ("PW", 0.0, &Vsrc::Instance::par5)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Pulse Width");

    p.addPar ("PER", 0.0, &Vsrc::Instance::par6)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Period");

    // Sin parameters
    p.addPar ("VA", 0.0, &Vsrc::Instance::par1)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Amplitude");

    p.addPar ("FREQ", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_SECM1)
     .setCategory(CAT_NONE)
     .setDescription("Frequency");

    p.addPar ("THETA", 0.0, &Vsrc::Instance::par4)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Theta");

    p.addPar ("PHASE", 0.0, &Vsrc::Instance::par5)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Phase");

    // Exp parameters
    p.addPar ("TD1", 0.0, &Vsrc::Instance::par2)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Rise Delay Time");

    p.addPar ("TAU1", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Rise Time Constant");

    p.addPar ("TD2", 0.0, &Vsrc::Instance::par4)
      .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Fall Delay Time");

    p.addPar ("TAU2", 0.0, &Vsrc::Instance::par5)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Fall Time Constant");

    // AC parameters
    p.addPar ("ACMAG", 0.0, &Vsrc::Instance::ACMAG)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Amplitude")
     .setAnalyticACSensitivityAvailable(true)
     .setACSensitivityFunctor(&acMagSens);

    p.addPar ("ACPHASE", 0.0, &Vsrc::Instance::ACPHASE)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Phase")
     .setAnalyticACSensitivityAvailable(true)
     .setACSensitivityFunctor(&acPhaseSens);

    // SFFM parameters
    p.addPar ("FC", 0.0, &Vsrc::Instance::par2)
     .setUnit(U_SECM1)
     .setCategory(CAT_NONE)
     .setDescription("Carrier Frequency");

    p.addPar ("FS", 0.0, &Vsrc::Instance::par4)
     .setUnit(U_SECM1)
     .setCategory(CAT_NONE)
     .setDescription("Signal Frequency");

    p.addPar ("MDI", 0.0, &Vsrc::Instance::par3)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Modulation Index");

    // PWL params
    p.addPar ("R", 0.0, &Vsrc::Instance::REPEATTIME)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Repeat Time");

    p.addPar ("T", 0.0, &Vsrc::Instance::T)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Time");  // time-voltage pairs

    p.addPar ("V", 0.0, &Vsrc::Instance::V)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Voltage"); // time-voltage pairs

    // PAT params
    p.addPar ("VHI", 0.0, &Vsrc::Instance::par0)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("High Voltage Value");

    p.addPar ("VLO", 0.0, &Vsrc::Instance::par1)
     .setUnit(U_VOLT)
     .setCategory(CAT_NONE)
     .setDescription("Low Voltage Value");

    p.addPar ("TSAMPLE", 0.0, &Vsrc::Instance::par5)
     .setUnit(U_SECOND)
     .setCategory(CAT_NONE)
     .setDescription("Pulse Width");

    p.addPar ("DATA", "", &Vsrc::Instance::DATA)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Data Pattern");

    p.addPar ("RB", 1, &Vsrc::Instance::RB)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("Starting Bit When Repeating");

    // Set up exceptions (ie variables that are not doubles):
    p.addPar ("TRANSIENTSOURCETYPE", (int) _DC_DATA, &Vsrc::Instance::TRANSIENTSOURCETYPE)
     .setGivenMember(&Vsrc::Instance::TRANSIENTSOURCETYPEgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    p.addPar ("ACSOURCETYPE", (int) _AC_DATA, &Vsrc::Instance::ACSOURCETYPE)
     .setGivenMember(&Vsrc::Instance::ACSOURCETYPEgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    p.addPar ("DCSOURCETYPE", (int) _DC_DATA, &Vsrc::Instance::DCSOURCETYPE)
     .setGivenMember(&Vsrc::Instance::DCSOURCETYPEgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    p.addPar ("NUM", 0, &Vsrc::Instance::NUM)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("" );

    // port parameters
    p.addPar ("PORT", 0, &Vsrc::Instance::port)
     .setGivenMember(&Vsrc::Instance::PORTgiven)
     .setUnit(U_NONE)
     .setCategory(CAT_NONE)
     .setDescription("port number");

    p.addPar ("Z0", 50.0, &Vsrc::Instance::Z0)
         
     .setGivenMember(&Vsrc::Instance::Z0given)
     .setUnit(U_OHM)
     .setCategory(CAT_NONE)
     .setDescription("impedance");
}

void Traits::loadModelParameters(ParametricData<Vsrc::Model> &p)
{
}


std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : SourceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
    srcCurrent(0.0),
    srcDrop(0.0),
    srcBC(0.0),
    HBSpecified_(factory_block.deviceManager_.getHBSpecified()),
    ACSpecified_(factory_block.deviceManager_.getACSpecified()),
    DCV0(0.0),
    par0(0.0),
    par1(0.0),
    par2(0.0),
    par3(0.0),
    par4(0.0),
    par5(0.0),
    par6(0.0),
    REPEATTIME(),
    T(0.0),
    V(0.0),
    ACMAG(1.0),
    ACPHASE(0.0),
    DATA(""),
    RB(1),
    NUM(0),
    REPEAT(false),

    TRANSIENTSOURCETYPE(_DC_DATA),
    TRANSIENTSOURCETYPEgiven(false),
    ACSOURCETYPE(_AC_DATA),
    ACSOURCETYPEgiven(false),
    DCSOURCETYPE(_AC_DATA),
    DCSOURCETYPEgiven(false),
    gotParams(false),
    v_pos(0.0),
    v_neg(0.0),
    i_bra(0.0),

    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_branch_data(0),

    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    APosEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    ABraEquBraVarOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0),
    fPosEquPosNodePtr(0),
    fNegEquNegNodePtr(0),
    fBraEquBraVarPtr(0),
#endif
   
    firstTimeload(true),
    port(0),
    Z0(50.0),
    PORTgiven (false),
    Z0given (false),
    freqVarsLoaded(false)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;
  setNumBranchDataVars(0);             // by default don't allocate space in branch vectors
  numBranchDataVarsIfAllocated = 1;    // this is the space to allocate if lead current or power is needed.

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

    // PDE supporting stamp.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(3);
    jacStampPDE[0].resize(2);
    jacStampPDE[0][0] = 0;
    jacStampPDE[0][1] = 2;
    jacStampPDE[1].resize(2);
    jacStampPDE[1][0] = 1;
    jacStampPDE[1][1] = 2;
    jacStampPDE[2].resize(3);
    jacStampPDE[2][0] = 0;
    jacStampPDE[2][1] = 1;
    jacStampPDE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  const SolverState &solver_state = factory_block.solverState_;
  const DeviceOptions &device_options = factory_block.deviceOptions_;

  // Set any non-constant parameter defaults:
  if (ACSpecified_ && ACSOURCETYPEgiven)
  {
    acSourceData_ = new ACData(*this, IB.params, solver_state, device_options);
  }

//  if (DCSOURCETYPEgiven) // this will always be given, if the source spec was valid.
    dcSourceData_ = new ConstData(*this, IB.params, solver_state, device_options);

  if (HBSpecified_ || TRANSIENTSOURCETYPEgiven)
  {
    switch (TRANSIENTSOURCETYPE)
    {
      case _SIN_DATA:
        tranSourceData_ = new SinData(*this, IB.params, solver_state, device_options);
        break;

      case _EXP_DATA:
        tranSourceData_ = new ExpData(*this, IB.params, solver_state, device_options);
        break;

      case _PULSE_DATA:
        tranSourceData_ = new PulseData(*this, IB.params, solver_state, device_options);
        break;

      case _PWL_DATA:
        tranSourceData_ = new PWLinData(*this, IB.params, solver_state, device_options);
        break;

      case _PAT_DATA:
        tranSourceData_ = new PatData(*this, IB.params, solver_state, device_options);
        break;

      case _SFFM_DATA:
        tranSourceData_ = new SFFMData(*this, IB.params, solver_state, device_options);
        break;

      case _DC_DATA:
        tranSourceData_ = 0; // this forces  us to use the dcSourceData_ object instead
        break;

      default:
        UserFatal(*this) << "Cannot identify source data type for " << getName();
        break;
    }
  }

  processParams();


/*  if (tranSourceData_ != 0)
  {

    tranSourceData_->updateSource();
    double  val = tranSourceData_->returnSource();

    if ( !DCSOURCETYPEgiven )
    { 
//      tranSourceData_->updateSource();
//      double  val = tranSourceData_->returnSource();
      setParam("DCV0", val,  true);
    }
    else
    {
      dcSourceData_->updateSource();
      double valDC = dcSourceData_->returnSource();



      if ( valDC != val )
      {
        if (ACSpecified_)
        {
          UserWarning(*this) << "The value specified in the DC field and the value from transient specification at time 0 is not consistent. Using the value in the DC field = " << valDC << " for DCOP calculation in .AC analysis";
        }
        else
        {
          setParam("DCV0", val,  true);
          UserWarning(*this) << "The value specified in the DC field and the value from transient specification at time 0 is not consistent. Using the value from transient specification at time 0 =  " << val << " for DCOP calculation";
        }

      }

    }


  }   */

//  processParams();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();
  processParams();

  // calculate dependent (ie computed) params and check for errors:

  if( PORTgiven )
  {                  

    jacStamp[2].resize(3);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;
    jacStamp[2][2] = 2;
  }
          
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/26/03
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  if (gotParams)
  {
    if (dcSourceData_ != 0)
    {
      dcSourceData_->setParams (&DCV0);
    }
    if (acSourceData_ != 0)
    {
      acSourceData_->setParams (&ACMAG);
    }
    if (tranSourceData_ != 0)
    {
      tranSourceData_->setParams(&par0);
    }
  }
  else
  {
    if (dcSourceData_ != 0)
    {
      dcSourceData_->getParams (&DCV0);
    }
    if (acSourceData_ != 0)
    {
      acSourceData_->getParams (&ACMAG);
    }
    if (tranSourceData_ != 0)
    {
      tranSourceData_->getParams(&par0);
    }
    gotParams = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
  delete tranSourceData_;
  delete acSourceData_;
  delete dcSourceData_;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
	                                const std::vector<int> & extLIDVecRef)
{
  std::string msg;

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  VsrcInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs") << "numInt != numIntVars";
  }

  if (numExt != numExtVars)
  {
    DevelFatal(*this).in("Instance::registerLIDs") << "numExt != numExtVars";
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Bra = intLIDVec[0];

  if (DEBUG_DEVICE && isActive(Diag::DEVICE_PARAMETERS) )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
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
/// The Voltage Source device uses exactly one "branch data vector" element, where
/// it keeps the "lead current" that may be used on .PRINT lines as
/// "I(V1)" for the current through resistor V1. and a junction voltage.
///
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
  if (loadLeadCurrent)
  {
    addBranchDataNode( symbol_table, li_branch_data, getName(), "BRANCH_D");
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/21/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (getSolverState().isPDESystem_)
    return jacStampPDE;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if (getSolverState().isPDESystem_)
  {
    APosEquBraVarOffset  = jacLIDVec[0][1];
    ANegEquBraVarOffset  = jacLIDVec[1][1];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }
  else
  {
    APosEquBraVarOffset  = jacLIDVec[0][0];
    ANegEquBraVarOffset  = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];

    if( PORTgiven )
      ABraEquBraVarOffset = jacLIDVec[2][2];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  Linear::Matrix & dFdx = *(extData.dFdxMatrixPtr);

  fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);

  if( PORTgiven )
    fBraEquBraVarPtr = &(dFdx[li_Bra][ABraEquBraVarOffset]);

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       : Loads the F-vector contributions for a single vsrc instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  srcCurrent        = solVec[li_Bra];
  srcDrop           = (solVec[li_Pos]-solVec[li_Neg]);

  if( PORTgiven && !getSolverState().spAnalysisFlag_ ) 
    srcDrop  -= srcCurrent * Z0;

  fVec[li_Pos] += srcCurrent;
  fVec[li_Neg] += -srcCurrent;
  fVec[li_Bra] += srcDrop;

  if( loadLeadCurrent )
  {
    double * leadF = extData.nextLeadCurrFCompRawPtr;
    double * junctionV = extData.nextJunctionVCompRawPtr;
    leadF[li_branch_data] = srcCurrent;
    junctionV[li_branch_data] = srcDrop;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEBVector ()
{
  // Get the value for the source.
  SourceData *dataPtr  = dcSourceData_; // by default assume the DC value.
  if ((HBSpecified_ || (getSolverState().tranopFlag && (!getSolverState().locaEnabledFlag || !DCSOURCETYPEgiven ) ) || getSolverState().transientFlag || (ACSpecified_ && !DCSOURCETYPEgiven ) ) && tranSourceData_ != 0 )
  {
    dataPtr = tranSourceData_;
  }

  if (dataPtr != 0)
  {
    srcBC = dataPtr->returnSource();
  }
  else
  {
    srcBC = 0.0;
  }

  double * bVec = extData.daeBVectorRawPtr;
  bVec[li_Bra] += srcBC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadBVectorsforAC
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool Instance::loadBVectorsforAC(double * bVecReal, double * bVecImag )
{
  if (acSourceData_ != 0)
  {
    bool flag = true;
    acSourceData_->setRealFlag(flag);

    acSourceData_->updateSource ();
    srcBC = acSourceData_->returnSource();

    bVecReal[li_Bra] += srcBC;

    flag = false;
    acSourceData_->setRealFlag(flag);

    acSourceData_->updateSource ();
    srcBC = acSourceData_->returnSource();

    bVecImag[li_Bra] += srcBC;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadFreqBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadFreqBVector (double frequency,
                                std::vector<Util::FreqVecEntry>& bVec)
{
      

  if ( !freqVarsLoaded )
    calculateFDVars();

  Util::FreqVecEntry tmpEntry;

  {

    std::complex<double> tmpVal = 0.0;
    double tol = 2.0*Util::MachineDependentParams::MachinePrecision();

    switch (TRANSIENTSOURCETYPE)
    {

      case _SIN_DATA:
      {
	if (frequency == 0.0 )
	  tmpVal = std::complex<double> ( v0, 0);

	if ( fabs(frequency - freq) < (frequency * tol  + tol ) )
	  tmpVal = std::complex<double> ( 0.5*mag*sin(phase), -0.5*mag*cos(phase) );

      }
      break;

      case _PULSE_DATA:
      {

	int fIdx;

	fIdx = std::round( frequency/freq);

//	double tol = 2.0*Util::MachineDependentParams::MachinePrecision();

	if ( ( fabs(frequency - freq * fIdx) < (frequency * tol  + tol ) ) &&  ( 2* fIdx + 1 <= size_ ) )
	  tmpVal = std::complex<double> ( ftOutData_[ 2* fIdx]/size_ , ftOutData_[ 2* fIdx + 1 ]/size_);

//	std::cout << "loaded value is " << tmpVal << std::endl;
      }
      break;

      case _DC_DATA:
      {
	if (frequency == 0.0 )
	  tmpVal = std::complex<double> ( v0, 0 );

 //  	std::cout << "loaded DC value is " << tmpVal << std::endl;
      }
      break;

      default:
        UserFatal(*this) << "Cannot identify source data type for " << getName();
        break;
    }
     // Add RHS vector element for the positive circuit node KCL equ.
    tmpEntry.val = tmpVal;
    tmpEntry.lid = li_Bra;
    bVec.push_back(tmpEntry);
     
  }        
     
  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadFreqBVector
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::calculateFDVars()
{

  SourceData *dataPtr  = dcSourceData_;
  if ( HBSpecified_ && tranSourceData_ != 0 )
  {
    dataPtr =  tranSourceData_;
  }

  if ( (dataPtr != 0)  )
  {

    switch (TRANSIENTSOURCETYPE)
    {

      case _SIN_DATA:
      {
        v0 = par0;

        mag = par1;

        freq = par3;

        phase = M_PI * par5/180;
      }
      break;

      case _PULSE_DATA:
      {

	dataPtr->setUseLocalTimeFlag(true);

	double dt = std::min( {par3, par4, par5} );

	int overSampleRate = 2;

	size_ = std::round( par6/dt ) * overSampleRate;

	if ( size_ % 2 == 0)
	  size_ = size_ + 1;

	double tstep = par6/size_;
	freq = 1/par6;

	ftInData_.resize( size_ );
	ftOutData_.resize( size_ +1 );
	iftInData_.resize( size_  +1 );
	iftOutData_.resize( size_ );

	if (ftInterface_ == Teuchos::null)
	{
	  ftInterface_ = Teuchos::rcp( new N_UTL_FFTInterface<std::vector<double> >( size_ ) );
	  ftInterface_->registerVectors( ftInData_, &ftOutData_, iftInData_, &iftOutData_ );
	} 


	for ( int i=0;  i < size_; ++i )
	{

	   dataPtr->setTime(i * tstep);

	   dataPtr->updateSource();

	   ftInData_[i] = dataPtr->returnSource();

  //         std::cout <<  "ftIndata " << i << " is " << ftInData_[i] << std::endl;
	}

        ftInterface_->calculateFFT();

        dataPtr->setUseLocalTimeFlag(false);

      }
      break;

      case _DC_DATA:
      {
        v0 = DCV0;
      }
      break;

      default:
        UserFatal(*this) << "Cannot identify source data type for " << getName();
        break;
    }

  }        
     
  freqVarsLoaded =true;  

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 vsrc instance.
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

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;
  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;


  if( PORTgiven && !getSolverState().spAnalysisFlag_ ) 
    dFdx[li_Bra][ABraEquBraVarOffset] -= Z0;

  return true;
}

// end of new-DAE functions

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize  ()
{
  double maxStep = 1.0e+100;
  if (tranSourceData_ != 0)
  {
    maxStep = tranSourceData_->getMaxTimeStepSize ();
  }
  return maxStep;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    DC_TRAN (0)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
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

// Additional Declarations

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

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
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
// Vsrc Master functions:
//-----------------------------------------------------------------------------

Master::Master(
  const Configuration & configuration,
  const FactoryBlock &  factory_block,
  const SolverState &   solver_state,
  const DeviceOptions & device_options)
  : DeviceMaster<Traits>(configuration, factory_block, solver_state, device_options),
    HBSpecified_(factory_block.deviceManager_.getHBSpecified()),
    ACSpecified_(factory_block.deviceManager_.getACSpecified()),
    separateInstances_(false)
  {}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * bVec, 
                             double * leadF, double * leadQ, double * junctionV, int loadType)
{
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR ))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & vi = *(*it);

    // Get the value for the source.
    SourceData *dataPtr  = vi.dcSourceData_; // by default assume the DC value


    if ( (getSolverState().tranopFlag && !vi.DCSOURCETYPEgiven && getSolverState().locaEnabledFlag ) && vi.tranSourceData_ != 0)
    {

      if (vi.firstTimeload)
      {
        double  val = vi.tranSourceData_->returnSource();

        vi.setParam("DCV0", val,  true);

        vi.dcSourceData_->setParams (&vi.DCV0);

        vi.firstTimeload = false;
      }
    }   

    if ((HBSpecified_ || (getSolverState().tranopFlag && (!getSolverState().locaEnabledFlag ) ) || getSolverState().transientFlag || (ACSpecified_ && !vi.DCSOURCETYPEgiven ) ) && vi.tranSourceData_ != 0 )
//    if ((HBSpecified_ || (getSolverState().tranopFlag && (!getSolverState().locaEnabledFlag ||  !vi.DCSOURCETYPEgiven ) ) || getSolverState().transientFlag || (ACSpecified_ && !vi.DCSOURCETYPEgiven ) ) && vi.tranSourceData_ != 0 )
    {
      dataPtr            = vi.tranSourceData_;
    }

    if (dataPtr != 0)
    {
      vi.srcBC           = dataPtr->returnSource();
    }
    else
    {
      vi.srcBC           = 0.0;
    }

    vi.srcCurrent        = solVec[vi.li_Bra];
    vi.srcDrop           = (solVec[vi.li_Pos]-solVec[vi.li_Neg]);

    if( vi.PORTgiven && !getSolverState().spAnalysisFlag_ ) 
      vi.srcDrop  -= vi.srcCurrent * vi.Z0;

    fVec[vi.li_Pos] += vi.srcCurrent;
    fVec[vi.li_Neg] += -vi.srcCurrent;
    fVec[vi.li_Bra] += vi.srcDrop;

    bVec[vi.li_Bra] += vi.srcBC;

    if( vi.loadLeadCurrent )
    {
      leadF[vi.li_branch_data] = vi.srcCurrent;
      junctionV[vi.li_branch_data] = vi.srcDrop;
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
bool Master::loadDAEMatrices (Linear::Matrix & dFdx, Linear::Matrix & dQdx, int loadType)
{
  InstanceVector::const_iterator it, end;

  if (loadType == LINEAR_FREQ)
    loadType = LINEAR;

  if (!separateInstances_ && ( loadType == LINEAR || loadType == NONLINEAR ))
  {
    separateInstanceTypes(linearInstances_, nonlinearInstances_);
    separateInstances_ = true;
  }

  if (loadType == ALL)
  {
    it = getInstanceBegin();
    end = getInstanceEnd();
  }
  else if (loadType == LINEAR)
  {
    it = linearInstances_.begin();
    end = linearInstances_.end();
  }
  else
  {
    it = nonlinearInstances_.begin();
    end = nonlinearInstances_.end();
  }

  for ( ; it != end; ++it )
  {
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(vi.fPosEquBraVarPtr) += 1.0;
    *(vi.fNegEquBraVarPtr) -= 1.0;
    *(vi.fBraEquPosNodePtr) += 1.0;
    *(vi.fBraEquNegNodePtr) -= 1.0;

    if( vi.PORTgiven && !getSolverState().spAnalysisFlag_ )
      *(vi.fBraEquBraVarPtr) -= vi.Z0;

#else
    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;
    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;
    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;
    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;

    if( vi.PORTgiven && !getSolverState().spAnalysisFlag_ ) 
      dFdx[vi.li_Bra][vi.ABraEquBraVarOffset] -= vi.Z0;

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
  if (deviceMap.empty() || (deviceMap.find("V")!=deviceMap.end())  || (deviceMap.find("P")!=deviceMap.end() ) )
  {
    Config<Traits>::addConfiguration()
      .registerDevice("v", 1)
      .registerDevice("p", 1);

  }
}


//-----------------------------------------------------------------------------
// Function      : dcVsrcSensitivity::operator
// Purpose       : produces df/dp and dq/dp, where p=DCV0.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 7/18/2014
//-----------------------------------------------------------------------------
void dcVsrcSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector<double> & dfdp,
    std::vector<double> & dqdp,
    std::vector<double> & dbdp,
    std::vector<int> & Findices,
    std::vector<int> & Qindices,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  dbdp.resize(1);
  dbdp[0] += 1.0;
  Bindices.resize(1);
  Bindices[0] = in->li_Bra;
}

//-----------------------------------------------------------------------------
// Function      : acMagVsrcSensitivity::operator
// Purpose       : B-vector sensitivities w.r.t. Mag parameter
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/1/2019
//-----------------------------------------------------------------------------
void acMagVsrcSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector< std::complex<double> > & dbdp,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);

  double mpi = M_PI;
  double phase = in->ACPHASE;
  double realpart=std::cos(2.0*mpi*phase/360);
  double imagpart=std::sin(2.0*mpi*phase/360);

  dbdp.resize(1);
  dbdp[0] += std::complex<double> (realpart,imagpart);
  Bindices.resize(1);
  Bindices[0] = in->li_Bra;
}

//-----------------------------------------------------------------------------
// Function      : acPhaseVsrcSensitivity::operator
// Purpose       : B-vector sensitivities w.r.t. Phase parameter
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 4/1/2019
//-----------------------------------------------------------------------------
void acPhaseVsrcSensitivity::operator()(
    const ParameterBase &entity,
    const std::string & name,
    std::vector< std::complex<double> > & dbdp,
    std::vector<int> & Bindices
    ) const
{
  const ParameterBase * e1 = &entity;
  const Instance * in = dynamic_cast<const Instance *> (e1);
  double mpi = M_PI;
  double mag = in->ACMAG;
  double phase = in->ACPHASE;
  double realpart= -(2.0*mpi*mag*std::sin(2.0*mpi*phase/360))/360;
  double imagpart= +(2.0*mpi*mag*std::cos(2.0*mpi*phase/360))/360;

  dbdp.resize(1);
  dbdp[0] += std::complex<double> (realpart, imagpart);
  Bindices.resize(1);
  Bindices[0] = in->li_Bra;
}

} // namespace Vsrc
} // namespace Device
} // namespace Xyce
